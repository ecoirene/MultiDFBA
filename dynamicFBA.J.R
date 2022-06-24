################################################
# Function: dynamicFBA
#
# Performs a dynamic flux balance analysis
#
# The function dynamicFBA() is inspired by the function
# dynamicFBA() contained in the COBRA Toolbox.
# The algorithm is the same.

# this allows us to have debugging info that can be turned off in production
DEVEL <- TRUE # are we in development mode?

DBG <- function(...) {
    if (DEVEL == TRUE) {
        cat(
            deparse(sys.call(-1)), # get my own name
            "\n   ",
            ...,
            "\n\n"
        )
    }
}

info <- function(...) {
    cat(
        "INFO:",
        deparse(sys.call(-1)), # get my own name
        ...,
        "\n"
    )
}

warn <- function(...) {
    cat(
        "WARNING:",
        deparse(sys.call(-1)), # get my own name
        ...,
        "\n"
    )
}

err <- function(..., exit = FALSE) {
    cat(
        "ERROR:",
        deparse(sys.call(-1)), # get my own name
        ...,
        "\n"
    )
    # this allows us to either exit on error or return something else
    if (exit == TRUE) stop()
}


# obtain the name of the script that calls this function
this_script_name <- function() {
    return(sys.calls()[[1]][2])
}


#' Get the path of the calling script
#'
#' \code{get_scriptpath} returns the full path of the script that called this function (if any)
#' or NULL if path is not available
#'
#' @examples
#' mypath <- get_scriptpath()
#' @export
#'
# https://stackoverflow.com/questions/18000708/find-location-of-current-r-file
# DonJ, answered Sep 9, 2016 at 14:00
#
get_script_path <- function() {
    # location of script can depend on how it was invoked:
    # source() and knit() put it in sys.calls()
    path <- NULL

    if (!is.null(sys.calls())) {
        # get name of script - hope this is consisitent!
        path <- as.character(sys.call(1))[2]
        # make sure we got a file that ends in .R, .Rmd or .Rnw
        if (grepl("..+\\.[R|Rmd|Rnw]", path, perl = TRUE, ignore.case = TRUE)) {
            return(path)
        } else {
            message("Obtained value for path does not end with .R, .Rmd or .Rnw: ", path)
        }
    } else {
        # Rscript and R -f put it in commandArgs
        args <- commandArgs(trailingOnly = FALSE)
    }
    return(path)
}


#library(styler)
restyle.file <- function(file) {
    # both are the same
    # style_file(file, style = tidyverse_style, indent_by = 4)
    styler::style_file(file, transformers = tidyverse_style(indent_by = 4))
}



buildWorld <- function(model,
                       substrateRxns,
                       initConcentrations,
                       initBiomass,
                       timeStep,
                       nSteps,
                       exclUptakeRxns,
                       retOptSol = TRUE,
                       fld = FALSE,
                       verboseMode = 2,
                       ...) {

    # check if model is a list of models or a single model
    DBG("Ensuring model is a list")
    nModels <- length(model)

    if (is.list(model)) {
        for (i in 1:nModels) {
            if (is(model[[i]], "modelorg")) {
                DBG("Model number", i, "is", model[[i]]@mod_name, "\n")
            } else if (!is(model[[i]], "modelorg")) {
                warn("Model", model[[i]]@mod_name, " is not a modelorg: IGNORED")
            }
        }
    } else if (is(model, "modelorg")) {
        model <- list(model)
    } else {
        err("argument model must be of type modelorg\n")
        return(NULL)
    }

    nModels <- length(model)

    # print the models we got
    for (i in 1:nModels) {
        m <- model[[i]]
        DBG("Model", i, "is a single model with name ", m@mod_name, "\n")
    }

    DBG("Ensuring we have some excluded uptake reactions\n")
    # Uptake reactions whose substrate concentrations do not change
    if (missing(exclUptakeRxns)) {
        # this is the default value (which I wonder why it was not
        # defined in the function call definition which would have made
        # it clear to any user)
        exclUptakeRxns <- c("EX_co2(e)", "EX_o2(e)", "EX_h2o(e)", "EX_h(e)")
        if (verboseMode > 2) {
            cat("Default extra cellular uptake reactions will be used: ")
            cat(exclUptakeRxns)
        }
    }
    # else, excluded uptake reactions is a list of reaction names whose
    # concentration does not change; a model does not need to have all
    # of them, but if it has, then those reactions will not affect the
    # concentration

    DBG("Ensuring Biomass is a list\n")
    biomass_list <- list()

    if (length(initBiomass) == 1) {
        # we got only one initial biomass value, so we will use this
        # same value for all models
        cat("Setting initial biomass to", initBiomass, "for all", nModels, "models\n")
        for (i in 1:nModels) {
            biomass_list[[i]] <- initBiomass
        }
    } else if (length(initBiomass) == nModels) {
        biomass_list <- initBiomass
    } else {
        err("We need only one, or as many as models, initial biomass(es)")
        return(NULL)
    }
    for (i in 1:nModels) {
        if (!is.numeric(biomass_list[[i]])) {
            err("Initial  biomass #", i, "is not numeric")
            return(NULL)
        }
    }

    excReactInd_list <- list()
    excRxnNames_list <- list()
    concentrations_list <- list()
    originalBound_list <- list()
    uptakeBound_list <- list()
    concentrationMatrix_list <- list()
    biomassVec_list <- list()
    timeVec_list <- list()
    aboveOriginal_list <- list()
    lpmod_list <- list()


    # Find exchange reactions for each model
    for (i in 1:nModels) {
        m <- model[[i]]
        biomass <- biomass_list[[i]]

        # Find exchange reactions for this model
        DBG("Finding exchange reactions for model", i, m@mod_name)

        # get the names of each exchange reaction
        excReact <- findExchReact(m)
        # get a logical vector telling which reactions are exchange reactions
        excReactInd <- (react_id(m) %in% react_id(excReact)) # excReact$exchange

        # find uptake reactions
        uptakeRxns <- uptReact(excReact)
        uptakeRxnsInd <- is.element(react_id(m), uptakeRxns)

        # now identify which excluded uptake reactions are present in this model
        # represent extra cellular reaction with boolean vector.
        exclUptakeRxnsInd <- is.element(react_id(m), exclUptakeRxns)

        # Exclude from the exchange reactions of this model those whose
        # concentrations will not change
        excReactInd <- excReactInd & !exclUptakeRxnsInd # excInd & ~ismember(model.rxns,exclUptakeRxns);
        # Now, add them to the ByModel list
        #  this is a list that tells which exchange reactions are
        #  not excluded in this model
        excReactInd_list[[i]] <- excReactInd

        # get reaction names of the selected exchange reactions
        excRxnNames <- react_id(m)[excReactInd] # excRxnNames = model.rxns(excInd)
        # and add them to the ByModel list
        excRxnNames_list[[i]] <- excRxnNames

        DBG("Checking substrateRxns for model", i, m@mod_name)
        # Some models might lack some substrate reactions or they could
        # be named differently in each model
        #
        # FOR NOW WE WILL REQUIRE THAT ALL MODELS CONTAIN EXCHANGE
        # REACTIONS WITH THE SAME NAMES AS THE SUBSTRATERXNS.
        #
        # this finds reactions in this model that are in substrateRxns
        substrateRxnsInd <- (react_id(m) %in% substrateRxns)
        # Figure out if substrate reactions are correct: all
        # substrate reactions should be non-excluded exchange reactions.
        # This finds if a substrate reaction present in this model
        # is not an exchange reaction (which would be invalid)
        missingSub <- substrateRxnsInd & !excReactInd
        DBG("substrateRxnsInd =", substrateRxnsInd)
        DBG("missingSub (head) = (", length(missingSub), ")", head(missingSub))
        if (sum(missingSub) != 0) {
            info(sum(missingSub))
            info(react_id(m)[missingSub])
            info("Invalid substrate uptake reaction!")
        }


        ## 	***********************************************************     ##
        # Initialize concentrations
        # substrateMatchInd = intersect(excRxnNames,substrateRxns);
        concentrations <- rep(0, length(react_id(m))) # table(excRxnNames);##vector(length=length(excRxnNames),mode="numeric");

        # this was wrong in the reference implementation:
        # substrateRxnsInd may not have an index for every substrateRxns
        # (it should, but may not if a substrateRxn is not in
        # react_id(model))
        # concentrations[1:length(concentrations)]=0;
        # concentrations[substrateRxnsInd] = initConcentrations;
        for (s in 1:length(substrateRxns)) {
            if (substrateRxns[s] %in% react_id(m)) {
                concentrations[react_id(m) == substrateRxns[s]] <-
                    initConcentrations[s]
            }
        }

        # Deal with reactions for which there are no initial concentrations
        originalBound <- -lowbnd(m) # take all to be able to directly update
        noInitConcentration <- (concentrations == 0) & (lowbnd(m) < 0) # (concentrations == 0 & originalBound > 0);
        concentrations[noInitConcentration] <- 1000

        # initialize bounds
        uptakeBound <- concentrations / (biomass * timeStep)

        # Make sure bounds are not higher than what are specified in the model
        aboveOriginal <- (uptakeBound > originalBound) & (originalBound > 0)
        uptakeBound[aboveOriginal] <- originalBound[aboveOriginal]
        lowbnd(m)[excReactInd] <- -uptakeBound[excReactInd]

        # add to the per-model lists
        originalBound_list[[i]]  <- originalBound
        concentrations_list[[i]] <- concentrations
        uptakeBound_list[[i]]    <- uptakeBound
        aboveOriginal_list[[i]]  <- aboveOriginal
        DBG("conc =", dim(concentrations), "\n", head(concentrations), "\n")
        DBG("conc_list =", str(concentrations_list))

        # initialize per-model time-series
        concentrationMatrix <- matrix(concentrations[excReactInd])
        row.names(concentrationMatrix) <- react_id(m)[excReactInd]

        concentrationMatrix_list[[i]] <- concentrationMatrix
        DBG("cMatrix =", dim(concentrationMatrix), "\n", concentrationMatrix, "\n")
        # DBG("cM_list =", str(concentrationMatrix_list), '\n')

        biomassVec <- biomass
        biomassVec_list[[i]] <- biomassVec
        
        timeVec <- c(0)
        timeVec_list[[i]] <- timeVec

        # prepare calculation
        ## ------------------------------------- Prepare Problem object --------------------##
        # get OptObj instance
        # lpmod <- prepProbObj(model,
        #                         nCols      = react_num(model),
        #                         nRows      = met_num(model),
        #             #            alg        = "FBA",
        #                         solver     = solver,
        #                         method     = method,
        #                         lpdir      = lpdir
        #                         #solverParm = solverParm
        #             )
        lpmod <- sybil::sysBiolAlg(m, algorithm = "fba", ...)
        lpmod_list[[i]] <- lpmod
    }
    # at  this point we have collected
    # model_list
    # biomass_list
    # excRxnNames_list
    # excReactInd_list
    # concentrations_list
    # originalBound_list
    # uptakeBound_list
    # concentrationMatrix_list
    # biomassVec_list
    # timeVec_list
    # aboveOriginal_list
    # lpmod_list

    names(initConcentrations) <- substrateRxns

    # return a list with all of them (world)
    world <- list(
        models = model,
        substrateRxns = substrateRxns,
        concentrations = initConcentrations,
        concentrations.per.model = concentrations_list,
        biomass.per.model = biomass_list,
        timeStep = timeStep,
        nSteps = nSteps,
        excRxnNames.per.model = excRxnNames_list,
        excReactInd.per.model = excReactInd_list,
        originalBound.per.model = originalBound_list,
        aboveOriginal.per.model = aboveOriginal_list,
        uptakeBound.per.model = uptakeBound_list,
        concentrationMatrix.per.model = concentrationMatrix_list,
        biomassVec.per.model = biomassVec_list,
        lpmod.per.model = lpmod_list,
        all_fluxes.per.model = list(), # these will be filled-in
        all_stat.per.model = list(), # during the simulation
        stepNo = 0,
        noFeasibleSolution = FALSE,
        timeVec.per.model = timeVec_list
    )
    return(world)
}

printWorld <- function(world) {
    info("World is")
    # print(str(world))
    nModels <- length(world$models)

    info("SubstrateRxns")
    print(world$substrateRxns)
    info("concentrations")
    print(world$concentrations)
    for (i in 1:nModels) {
        info(paste("concentrations as seen by model", i, world$models[[i]]@mod_name))
        info(world$concentrations.per.model[[i]])

        info("concentrationMatrix for model", i, world$models[[i]]@mod_name)
        concentrationMatrix <- world$concentrationMatrix.per.model[[i]]
        print(paste("cM dim =", dim(concentrationMatrix)))
        for (t in 1:dim(concentrationMatrix)[2]) {
            info("tStep =", t, "concMat length =", length(concentrationMatrix[, t]))
            info("conMat =", concentrationMatrix[, t])
        }

        m <- world$models[[i]] # we need [[ ]] to get the actual model
        cat("model", i, m@mod_name, "\n")
        print("    Current Biomass")
        cat("    ", world$biomass.per.model[[i]], "\n") # as is this
        print("    Biomass hitory (head)")
        cat("    ", head(world$biomassVec.per.model[[i]]), "\n")
        print("    Exchange reaction names (head)")
        cat("    ", head(world$excRxnNames.per.model[[i]]), "\n")
        print("    Exchange reaction indexes (head)")
        cat("    ", head(world$excReactInd.per.model[[i]]), "\n")
        print("    Uptake bounds (head)")
        cat("    ", world$uptakeBound.per.model[[i]], "\n")
        print("    Original bounds (head)")
        cat("    ", head(world$originalBound.per.model[[i]]), "\n")
        print("    Above original bounds (head)")
        cat("    ", head(world$aboveOriginal[[i]]), "\n")
        cat("\n")
    }
}


updateWorld <- function(world, stepNo, fld, verboseMode, ...) {
    # extract common world variables
    substrateRxns <- world$substrateRxns
    common.concentrations <- world$concentrations
    timeStep <- world$timeStep
    nSteps <- nSteps
    nModels <- length(world$models)

    noFeasibleSolution <- FALSE # we do not know yet if there is a solution

    for (i in 1:nModels) {
        model <- world$models[[i]]
        DBG("Doing step", stepNo, "for model", i, model@mod_name)
        # extract per-model variables
        concentrations <- world$concentrations.per.model[[i]]
        biomass <- world$biomass.per.model[[i]]
        excRxnNames <- world$excRxnNames.per.model[[i]]
        excReactInd <- world$excReactInd.per.model[[i]]
        originalBound <- world$originalBound.per.model[[i]]
        aboveOriginal <- world$aboveOriginal.per.model[[i]] # unneeded?
        uptakeBound <- world$uptakeBound.per.model[[i]]
        concentrationMatrix <- world$concentrationMatrix.per.model[[i]]
        biomassVec <- world$biomassVec.per.model[[i]]
        lpmod <- world$lpmod.per.model[[i]]
        timeVec <- world$timeVec.per.model[[i]]
        all_stat <- c()
        all_fluxes <- matrix()
        
        #DBG(sum(is.na(concentrations)), length(concentrations), "concentrations", head(concentrations))
        #DBG(sum(is.na(substrateRxns)), length(substrateRxns), "substrateRxns", head(substrateRxns))
        #DBG(sum(is.na(common.concentrations)), length(common.concentrations), "common.concentrations", head(common.concentrations))
        # 'concentrations' needs to be updated from common.concentrations:
        # 	when it was calculated, in a previous per-model loop,
        # 	we didn't know yet the effect of posterior models, so
        # 	the values for substrateRxns are likely invalid (for
        # 	all but the last model in the list)
        for (s in 1:length(substrateRxns)) {
            if (substrateRxns[s] %in% react_id(model)) {
                concentrations[react_id(model) == substrateRxns[s]] <-
                    common.concentrations[s]
            }
        }
        #DBG(sum(is.na(concentrations)), length(concentrations), "concentrations", head(concentrations))


        # Run FBA
        sol <- sybil::optimizeProb(lpmod)
        mu <- sol$obj ## objvalue sol.f
        DBG("step", stepNo, "status", sol$stat, "obj.func", mu, "\n")
        if (length(checkSolStat(sol$stat, solver(problem(lpmod)))) != 0) { ## checkSolStat
            print("No feasible solution - nutrients exhausted\n")
            noFeasibleSolution <- TRUE
            break
        }
        all_stat <- c(all_stat, sol$stat)
        uptakeFlux <- sol$fluxes[excReactInd]
        biomass <- biomass * exp(mu * timeStep)
        # biomass = biomass*(1+mu*timeStep);
        biomassVec <- c(biomassVec, biomass)

        # waitbar(stepNo/nSteps,h);
        timeVec <- c(timeVec, stepNo * timeStep)
        if (fld) {
            if (stepNo == 1) {
                all_fluxes <- sol$fluxes
            } else {
                all_fluxes <- cbind(all_fluxes, sol$fluxes)
            }
        }
        # Update concentrations
        # This works because uptakeFlux is sol$fluxes[excReactInd], thus,
        # concentrations[excReactInd] and uptakeFlux both refer to the same
        # subet of reactions.
        deltaConc <- uptakeFlux /
            mu * biomass * (1 - exp(mu * timeStep))
        concentrations[excReactInd] <-
            concentrations[excReactInd] - deltaConc
        #DBG(sum(is.na(concentrations)), length(concentrations), "concentrations", head(concentrations))
        # concentrations = concentrations + uptakeFlux*biomass*timeStep;
        concentrations[concentrations <= 0] <- 0
        concentrationMatrix <- cbind(concentrationMatrix, concentrations[excReactInd])

        #DBG(sum(is.na(uptakeBound)), length(uptakeBound), "uptakeBound", head(uptakeBound))
        #DBG(sum(is.na(concentrations)), length(concentrations), "concentrations", head(concentrations))
        #DBG(sum(is.na(excReactInd)), length(excReactInd), "excReactInd", head(excReactInd))
        #DBG("biomass", biomass)
        #DBG("timeStep", timeStep)
        # Update bounds for uptake reactions
        uptakeBound[excReactInd] <-
            concentrations[excReactInd] /
                (biomass * timeStep)

        # This is to avoid any numerical issues
        #DBG(sum(is.na(uptakeBound)), length(uptakeBound), "uptakeBound", head(uptakeBound))
        uptakeBound[uptakeBound > 1000] <- 1000

        # Figure out if the computed bounds were above the original bounds
        #DBG("before assigning aboveOriginal")
        #DBG(sum(is.na(originalBound)), length(originalBound), "originalBound", head(originalBound))
        #DBG(sum(is.na(aboveOriginal)), length(aboveOriginal), "aboveOriginal", head(aboveOriginal))
        #DBG(sum(is.na(uptakeBound)), length(uptakeBound), "uptakeBound", head(uptakeBound))
        aboveOriginal <- (uptakeBound > originalBound) & (originalBound > 0)

        # Revert to original bounds if the newly computed rate was too high
        #DBG("assigning uptakeBound")
        #DBG(sum(is.na(aboveOriginal)), length(aboveOriginal), "aboveOriginal", head(aboveOriginal))
        #DBG(sum(is.na(originalBound)), length(originalBound), "originalBound", head(originalBound))
        #DBG(sum(is.na(uptakeBound)), length(uptakeBound), "uptakeBound", head(uptakeBound))
        uptakeBound[aboveOriginal] <- originalBound[aboveOriginal] # uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
        # this is to avoid calculation overflows
        uptakeBound <- ifelse(abs(uptakeBound) < 1e-9, 0, uptakeBound)

        ## Change lower bounds according to the result of last step
        # lowbnd(model)[excReactInd]  = -uptakeBound[excReactInd];
        uppb_tmp <- getColsUppBnds(problem(lpmod), which(excReactInd))
        changeColsBnds(problem(lpmod), which(excReactInd), lb = -uptakeBound[excReactInd], ub = uppb_tmp)

        if (verboseMode > 2) print(paste(stepNo, sep = "    ", biomass))

if (TRUE) {
        ##### UPDATE VALUES FOR NEXT ITERATION
        # update common.concentrations
        #DBG(sum(is.na(uptakeFlux)), length(uptakeFlux), "uptakeFlux", head(uptakeFlux))
        # we cannot use uptakeFlux because it is a subset of all the fluxes, and hence
        # its dimension is smaller than all the reactions, hence, we need sol$fluxes
        for (s in 1:length(substrateRxns)) {
            present.in.model <- substrateRxns[s] %in% react_id(model)
            #cat("substrateRxns[", s,"]", present.in.model, '\n')
            if (present.in.model) {
                substIndex <- react_id(model) == substrateRxns[s]
                #cat("    index in model", length(substIndex), sum(substIndex), which(TRUE == substIndex), '\n')
                #cat("    comm.conc[", s, "] =", common.concentrations[s], '\n')
                #cat("    uptakeFlux[", which(TRUE == substIndex), "] =", uptakeFlux[which(TRUE == substIndex)], uptakeFlux[substIndex], '\n')
                #cat("    mu =", mu, "biomass =", biomass, "timStep =", timeStep, '\n')
                common.concentrations[s] <-
                    common.concentrations[s] -
                    (sol$fluxes[substIndex] /
                        mu * biomass * (1 - exp(mu * timeStep)))
            }
        }
        #DBG(sum(is.na(common.concentrations)), length(common.concentrations), "common.concentrations", head(common.concentrations))
} else {
        # this is equivalent to the above but restricting the search to exchange reactions
        for (s in 1:length(substrateRxns)) {
            present.in.model <- substrateRxns[s] %in% react_id(model)[excReactInd]
            #cat("substrateRxns[", s,"]", present.in.model, '\n')
            if (present.in.model) {
                # get position of this substrate in uptakeFlux
                substIndex <- react_id(model)[excReactInd] == substrateRxns[s]
                #DBG("    index in model", length(substIndex), sum(substIndex), which(TRUE == substIndex), '\n')
                #DBG("    comm.conc[", s, "] =", common.concentrations[s], '\n')
                #DBG("    uptakeFlux[", which(TRUE == substIndex), "] =", uptakeFlux[which(TRUE == substIndex)], uptakeFlux[substIndex], '\n')
                #DBG("    mu =", mu, "biomass =", biomass, "timStep =", timeStep, '\n')
                common.concentrations[s] <-
                    common.concentrations[s] -
                    (uptakeFlux[substIndex] /
                        mu * biomass * (1 - exp(mu * timeStep)))
            }
        }
        #DBG(sum(is.na(common.concentrations)), length(common.concentrations), "common.concentrations", head(common.concentrations))
}
        #DBG("Step", stepNo, "model", i, "done\n")
        # update per-model world variables
        world$concentrations.per.model[[i]] <- concentrations
        world$biomass.per.model[[i]] <- biomass
        world$excRxnNames.per.model[[i]] <- excRxnNames # unneeded?
        world$excReactInd.per.model[[i]] <- excReactInd # unneeded?
        world$originalBound.per.model[[i]] <- originalBound # unneeded?
        world$aboveOriginal.per.model[[i]] <- aboveOriginal
        world$uptakeBound.per.model[[i]] <- uptakeBound
        world$concentrationMatrix.per.model[[i]] <- concentrationMatrix
        DBG("biomassVec = (", length(biomassVec), ")", biomassVec)
        world$biomassVec.per.model[[i]] <- biomassVec
        DBG("timeVec = (", length(timeVec), ")", timeVec)
        world$timeVec.per.model[[i]] <- timeVec
        world$lpmod.per.model[[i]] <- lpmod
        world$all_fluxes.per.model[[i]] <- all_fluxes
        world$all_stat.per.model[[i]] <- all_stat
        world$timeVec.per.model[[i]] <- timeVec
    } # end per-model loop

    #DBG("Step", stepNo, "done for all", nModels, "models\n\n")
    # update common variables
    world$concentrations <- common.concentrations
    world$substrateRxns <- substrateRxns # unneeded?
    world$timeStep <- timeStep # unneeded?
    world$nSteps <- nSteps # unneeded?
    world$stepNo <- stepNo
    world$noFeasibleSolution <- noFeasibleSolution

    return(world)
}

generateOutput <- function(world, fld, retOptSol) {
    nModels <- length(world$models)

    ###    if (nModels == 1) {
    ###        # extract variables returned
    ###        model <- world$models[[1]]
    ###        lpmod <- world$lpmod.per.model[[1]]
    ###        stepNo <- world$stepNo
    ###        all_fluxes <- world$all_fluxes.per.model[[1]]
    ###        concentrationMatrix <- world$concentrationMatrix.per.model[[1]]
    ###        excRxnNames <- world$excRxnNames.per.model[[1]]
    ###        tmVec <- world$timeVec.per.model[[1]]
    ###        biomassVec <- world$biomassVec.per.model[[1]]
    ###
    ###
    ###        ## Prepare OUTPUT
    ###        #concentrationMatrix,excRxnNames,timeVec,biomassVec
    ###        if (isTRUE(retOptSol)) {
    ###            if(is.null(all_fluxes)) all_fluxes=as.matrix(NA);
    ### 	        return (optsol_dynamicFBA(
    ###                            solver = solver(problem(lpmod)),
    ### 			    method = method(problem(lpmod)),
    ### 			    nprob  = stepNo,
    ### 			    ncols  = react_num(model),
    ### 			    nrows  = met_num(model),
    ### 			    fld    = fld,
    ### 			    all_fluxes = all_fluxes,
    ### 			    concmat=concentrationMatrix,
    ### 			    exRxn=excRxnNames,
    ### 			    tmVec=tmVec,
    ###                            bmVec=biomassVec
    ### 		      )
    ### 	      )
    ###        } else {
    ### 	    return(optsol <- list(
    ###            		    nprob  = stepNo,
    ### 			    ncols  = react_num(model),
    ### 			    nrows  = met_num(model),
    ### 			    all_fluxes = all_fluxes,
    ### 			    all_stat=all_stat,
    ### 			    concentrationMatrix=concentrationMatrix,
    ### 			    excRxnNames=excRxnNames,
    ### 			    timeVec=tmVec,
    ### 			    biomassVec=biomassVec
    ### 		   ))
    ###        }
    ###    } else { 	# length(world$models) > 1
    # we will return either a list of optsol results (one per model)
    # or a list of result lists (one per model)
    optsol_list <- list()

    for (i in 1:nModels) {
        # extract variables returned
        model <- world$models[[i]]
        lpmod <- world$lpmod.per.model[[i]]
        stepNo <- world$stepNo
        all_fluxes <- world$all_fluxes.per.model[[i]]
        all_stat <- world$all_stat.per.model[[i]]
        concentrationMatrix <- world$concentrationMatrix.per.model[[i]]
        excRxnNames <- world$excRxnNames.per.model[[i]]
        tmVec <- world$timeVec.per.model[[i]] # timeVec.per.model[[i]]
        biomassVec <- world$biomassVec[[i]]

        ## Prepare OUTPUT
        # concentrationMatrix,excRxnNames,timeVec,biomassVec
        if (isTRUE(retOptSol)) {
            if (is.null(all_fluxes)) all_fluxes <- as.matrix(NA)
            optsol <- (optsol_dynamicFBA(
                solver = solver(problem(lpmod)),
                method = method(problem(lpmod)),
                nprob = stepNo,
                ncols = react_num(model),
                nrows = met_num(model),
                fld = fld,
                all_fluxes = all_fluxes,
                concmat = concentrationMatrix,
                exRxn = excRxnNames,
                tmVec = tmVec,
                bmVec = biomassVec
            ))
        } else {
            optsol <- list(
                nprob = stepNo,
                ncols = react_num(model),
                nrows = met_num(model),
                all_fluxes = all_fluxes,
                all_stat = all_stat,
                concentrationMatrix = concentrationMatrix,
                excRxnNames = excRxnNames,
                timeVec = tmVec,
                biomassVec = biomassVec
            )
        }

        optsol_list[[i]] <- optsol
    }
    if (nModels == 1) {
        return(optsol_list[[1]])
    } # for backwards compatibility
    else {
        return(optsol_list)
    }

    ###    }
}


multiDynFBA <- function(model,
                        substrateRxns,
                        initConcentrations,
                        initBiomass,
                        timeStep,
                        nSteps,
                        exclUptakeRxns,
                        retOptSol = TRUE,
                        fld = FALSE,
                        verboseMode = 2,
                        ...) {
    # PARAMETERS:
    # ===========
    # model                 Sybil model structure (class modelorg) (it can be either a single
    #                       model of class modelorg or a list of many models
    #                       each of them of class modelorg)
    # substrateRxns         List of exchange reaction names for substrates
    #                       initially in the media that may change (e.g. not
    #                       h2o or co2)
    # initConcentrations    Initial concentrations of substrates (in the same
    #                       structure as substrateRxns)
    # initBiomass           Initial biomass (must be non zero)
    # timeStep              Time step size
    # nSteps                Maximum number of time steps
    # fld                   indicates if all fluxes at all steps will be returned.
    # retOptSol             indicates if optsol calss will be returned or simple list
    #
    # OPTIONAL PARAMETERS
    # ===================
    # exclUptakeRxns        List of uptake reactions whose substrate
    #                       concentrations do not change (Default =
    #                       {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'})
    #
    # RETURN VALUES:
    # =============
    # This function may return either a single value (if called with a single
    # model= that matches dynamicFBA return, or a list of values (if ccalled
    # with multiple models.
    #
    # In case of retOptSol == TRUE, the value or list of values returned is
    # of type optsol_dynamicFBA, and if retOptSol == FALSE, the return value
    # is a list with the elementes listed below or a list of lists with the
    # following elements:
    #
    # concentrationMatrix   Matrix of extracellular metabolite concentrations
    # excRxnNames           Names of exchange reactions for the EC metabolites
    # timeVec               Vector of time points
    # biomassVec            Vector of biomass values
    # all_fluxes            Matrix containing the fluxes of all reactions at different steps

    world <- buildWorld(model,
        substrateRxns = substrateRxns,
        initConcentrations = initConcentrations,
        initBiomass = initBiomass,
        timeStep = timeStep,
        nSteps = nSteps,
        exclUptakeRxns = exclUptakeRxns,
        retOptSol = retOptSol,
        fld = fld,
        verboseMode = verboseMode,
        ...
    )

    if (is.null(world)) {
        return(NULL)
    }

    if (verboseMode <= 10) {
        # DBG("world contains the following slots:")
        # DBG(names(world))
        # DBG(str(world))
        printWorld(world)
    }

    nModels <- length(world$models)

    ## -----------------------------------------------------------------------------##
    if (verboseMode > 2) print("Step number    Biomass\n")
    # Inititialize progress bar ...');
    # if (verboseMode == 2)  progr <- .progressBar();

    for (stepNo in 1:nSteps) {
        #    if (verboseMode == 2)  progr <- .progressBar(stepNo, nSteps, progr);
        world <- updateWorld(world, stepNo, fld, verboseMode, ...)

        if (is.null(world)) {
            return(NULL)
        } # a catastrophic error occurred

        if (world$noFeasibleSolution == TRUE) break # no need to continue
    }
    return(generateOutput(world, fld, retOptSol))
}
