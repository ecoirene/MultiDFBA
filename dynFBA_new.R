library(sybil)
library(sybilSBML)
library(sybilDynFBA)

options(error=function() {traceback(); q() })

#source("optsol_dynamicFBA-I.R", keep.source=T)
source("optsol_dynamicFBA.R", keep.source=T)
#source("dynamicFBA-I.R", keep.source=T)
source("dynamicFBA.J.R", keep.source=T)
source("dynamicFBA.orig.R", keep.source=T)

#modelName <- 'nmR_l0';
#modelName <- readline(prompt="Model name (without extension): ")
modelName <- "MODELaa";

# objective function
#obj <- 'BIOMASS';
obj <- 'BIOMASS_SCO';	# For S. coelicolor
#obj <- 'EX_AMLB';
ci <- 1;
#obj <- c('EX_BIOMASS', 'EX_AMLB');
#ci <- c(1, 1);

# time step and length of simulation
ts <- 1;	# hr
ns <- 60;	# steps

inoculum <- 0.01;

calcExtras <- TRUE

# Load model
# ----------
#
# first try SBML and if that fails, try TSV
#
xmlName <- paste(modelName, ".xml", sep="");
tsvName <- modelName;	# files must be called tsvName_{react|met|desc}.tsv


if (file.exists(xmlName)) {
	model <- readSBMLmod(xmlName);
} else if (file.exists(paste(tsvName, "_react.tsv", sep=""))) {
	# we only need the reaction table as a minimum
	model <- readTSVmod(tsvName);
} else {
 	stop('No suitable model (.xml or .tsv) found');
}


# set up output file names
# ------------------------
#
tsvName   <- paste( "mod/", modelName, sep="" );	# files will be named modelName_{desc|met|react}.tsv
sbmlName  <- paste( "mod/", modelName, "_L2v1.xml", sep="" );
logName   <- paste("out/", modelName, ".Rlog", sep="");
pngFVAplot <- paste("out/", modelName, "_FVA.png", sep="")
pngRAplot <- paste("out/", modelName, "_RA.png", sep="")
pngPPPAplot <- paste("out/", modelName, "_PPPA.png", sep="")
pngDFBAplot <- paste("out/", modelName, "_DFBA.png", sep="")

logFile <- file(logName, open="wt");
sink(logFile, type=c("output", "message"), split=TRUE);

# we read the model as SBML (xml). Save it also as TEXT
if (! file.exists(paste(tsvName, '_react.tsv', sep=''))) {
    cat('Writing file ', tsvName, '_{desc|met|react}.tsv\n')
   modelorg2tsv(model, tsvName)
}
# and in a more modern SBML format
if (! file.exists(sbmlName)) {
    cat('Writing file ', sbmlName, '\n')
    writeSBML(model, filename=sbmlName, level=2, version=1)
}

# --------------------------------------------
# This model uses as units mmol per gDW per hr
# gDW (grams Dry Weight)
# --------------------------------------------

exchngRxns <- findExchReact(model)
uptakeRxns <- uptReact(exchngRxns)

cat('\nEXCHANGE REACTIONS:\n===================\n')
print(exchngRxns)
cat('\nUPTAKE REACTIONS:\n=================\n')
print(uptakeRxns)

# Set-up variables for simulation
#================================
#
# Initial concentrations
# ----------------------
# initial concentrations are in millimoles per gram[dry weight] per hour ???
#	i.e the substrate concentration is scaled to define the amount
#	of substrate available per unit of biomass per unit of time.
#	Varma and Palsson:
#	For E.coli, the maximum O2 utilization rate was 15 mmol of O2 per
#	g (dry weight) per h.
#	For E.coli aerobic cultures, the ratio of the growth rate to the
#	biomass yield was 10.5 mmol of Glc per g (dry weight) per h
#	For E. coli anaerobic cultures the maximum glucose utilization
#	rate was determined to be 18.5 mmol per g (gry weight) per h
# 
# Our values are:
# 
nonAaMets = 	c('EX_glc(e)', 'EX_mnl(e)');
nonAaMets_mg_per_l = c( 10000,     10000 );
nonAaMets_mwt = c( 180.15588,	182.172	 );
nonAaMets_mmol_per_l = nonAaMets_mg_per_l / nonAaMets_mwt;

#
# Casamino acids using composition from D'Huys 2011
# -------------------------------------------------
#	NOTE: we'll use 1e-01 for N/A and 1e-02 for 0 values
aa =  c('ALA',	'ARG',	'ASN',	'ASP',	'CYS',	'GLU',	'GLN',	'GLY',	'HIS',	'ILE',	'LEU',	'LYS',	'MET',	'PHE',	'PRO',	'SER',	'THR',	'TRP',	'TYR',	'VAL');
exchange_aa =  c('EX_ala_L(e)',	'EX_arg_L(e)',	'EX_asn_L(e)',	'EX_asp_L(e)',	'EX_cys_L(e)',	'EX_glu_L(e)',	'EX_gln_L(e)',	'EX_gly(e)',	'EX_his_L(e)',	'EX_ile_L(e)',	'EX_leu_L(e)',	'EX_lys_L(e)',	'EX_met_L(e)',	'EX_phe_L(e)',	'EX_pro_L(e)',	'EX_ser_L(e)',	'EX_thr_L(e)',	'EX_trp_L(e)',	'EX_tyr_L(e)',	'EX_val_L(e)');
#mg_per_g = c(17.8,	-0,	0,	42.56,	-0,	102.9,	0,	12,	12.4,	14.41,	39.3,	23.36,	11.92,	19.8,	46,	20,	10.7,	0,	14.48,	31.8);
mg_per_g = c(17.8,	1e-1,	1e-2,	42.56,	1e-1,	102.9,	1e-2,	12,	12.4,	14.41,	39.3,	23.36,	11.92,	19.8,	46,	20,	10.7,	1e-2,	14.48,	31.8);
aa_mwt =   c(89.0935,	174.2017,	132.1184,	133.1032,	121.1590,	147.1299,	146.1451,	75.0669,	155.1552,	131.1736,	131.1736,	146.1882,	149.2124,	165.1900,	115.1310,	105.0930,	119.1197,	204.2262,	181.1894,	117.1469)
aa_mmol_per_g_l = mg_per_g / aa_mwt; 	# .* 1 (per 1 litre)
aa_mmol_per_5g_l = aa_mmol_per_g_l * 5;
aa_mmol_per_15g_l = aa_mmol_per_g_l * 15;
#
#substrateRxns = c(exchange_major_components_minus_aa,
#               exchange_aa);
#initConcentrations = c(concentration_of_major_components_minus_aa ,
#               aa_mmol_per_5g_l);
#exclUptakeRxns = c(major_components_in_excess,
#		minor_components_in_excess, 'EX_CO2', 'EX_O2', 'EX_H2O', 'EX_H');

#zeroRxns = c('EX_HMP', 'EX_ACAL', 'EX_FAN', 'EX_GL', 'EX_GLYN', 'EX_GN', 'EX_HYXN', 'EX_NAC', 'EX_T3', 'EX_UREA');
#zeroConc = c(1.e-5,	1.e-5,	1.e-5,	1.e-5,	1.e-5,	1.e-5,	1.e-5,	1.e-5,	1.e-5,	1.e-5);

# Set-up variables for dynamicFBA
#substrateRxns = c('EX_CO2'); 
#initConcentrations = c(10);
# 
#substrateRxns = c( nonAaMets );	# 5g GLC 5g MNT
#initConcentrations = nonAaMets_mmol_per_l;
# 0.0001 glucose and 0.0274 mannitol
#initConcentrations = c( 1.e-4, nonAaMets_mmol_per_l(2) );
#
#substrateRxns = c( nonAaMets, zeroRxns );	# 5g GLC 5g MNT
#initConcentrations = ( 1.e-4, nonAaMets_mmol_per_l(2), zeroConc );
# S. lividans grows but mannitol is exhausted very fast, hence there must be other
# carbon sources
#
#substrateRxns = c( nonAaMets, exchange_aa );
#initConcentrations = c( nonAaMets_mmol_per_l, aa_mmol_per_5g_l );
exclUptakeRxns = c()
#

#exclUptakeRxns = c(major_components_in_excess,
#		minor_components_in_excess, 'EX_PI','EX_NA','EX_NH3');


# For S. coelicolor (Kim et al.)
substrateRxns = c('EX_glc(e)', 'EX_mnl(e)');
initConcentrations = c(55.51, 54.89);		# mmol/L (10g)
#initConcentrations = c(10, 10);		# mmol/L
#substrateRxns = c('EX_glc(e)');
#initConcentrations = c(55.51);		# mmol/L (10g)
#substrateRxns = c('EX_mnl(e)');
#initConcentrations = c(54.89);		# mmol/L (10g)
substrateRxns = c('EX_glc(e)', 'EX_mnl(e)', exchange_aa);
initConcentrations = c(55.51, 54.89, aa_mmol_per_5g_l);		# mmol/L (10g)
#exclUptakeRxns = c('EX_h2o(e)', 'EX_h(e)', 'EX_o2(e)', 'EX_co2(e)');


# Objective function(s)
# ---------------------
# default is BIOMASS
#
# the objective coefficients specify the wwight given to the
# respective reaction in the rxnNameList
#
# according to the TNF XML file, the default objective is R701,
#	EX_BIOMASS (exchange of Biomass)
#
#model <- changeObjFunc(model, c('EX_BIOMASS', 'EX_AMLB'), c(0.5, 0.5))
#model <- changeObjFunc(model, 'EX_BIOMASS', c(1))
model <- changeObjFunc(model, obj, ci)



# initial Biomass
# ---------------
#
initBiomass <- inoculum;

# Simulation granularity
# ----------------------
#
# Duplication time in S.lividans should be ~4-6h
timeStep <- ts;		# set at the beginning of this file
nSteps <-   ns;		# ditto.

# graphical output
# ----------------
#plotRxns <- ('EX_GLC', 'EX_CO2', 'EX_O2', 'EX_MNT', 'EX_PI', 'EX_ALA', 'EX_AMLB');
plotRxns <- c(substrateRxns);
#plotRxns = c('EX_glc(e)', 'EX_mnl(e)');
#plotRxns <- c( 'EX_GLC',	'EX_MNT',	'EX_PI', 	'EX_AMLB');


# CALCULATIONS:
# =============          PUMP_NA_SER           PUMP_NA_THR        XCHNG_CADA_LYS 

#

if (calcExtras == TRUE) {

    # # Standard FBA 
    # ------------ 
    # 
    # Perform  ux-balance analysis (FBA) by using method optimizeProb of 
    # class modelorg . Method optimizeProb performs flux-balance analysis 
    # [Edwards et al., 2002a, Orth et al., 2010b]. 
    #
    cat('\nStandard Flux Balance Analysis(FBA)')
    cat('\n-----------------------------------\n')

    optL <- optimizeProb(model, algorithm = "fba", retOptSol = FALSE);
    #
    opt <- optimizeProb(model, algorithm = "fba", retOptSol = TRUE);
    #
    cat('\nResult\n')
    print( checkOptSol(opt) )

    cat('\nValue of the objective function:\n')
    print( lp_obj(opt) )

    # The function summaryOptsol() returns an object of class optsolSummary  
    # and needs the object of class optsol (results of simulations) and the 
    # corresponding object of class modelorg (the entire metabolic network).  
    # The generated object of class summaryOptsol contains some information 
    # about the flux distribution, substrates, products and limiting reactions. 
    #
    # The method printExchange() prints a subset of the flux distribution for 
    # the exchange reactions in the model.  Each column represents the 
    # environment of one optimization. The symbol "-" indicates
    # that the corresponding metabolite is imported (is a substrate); 
    # the symbol "+" indicates, that the corresponding metabolite is 
    # excreted (is a product).
    #
    cat('\nSummary:\n')
    sum <- summaryOptsol(opt, model)
    print(sum)
    # or only exchange reactions:
    #printExchange(sum, dense = TRUE)

    #


    # Minimize total flux
    # -------------------
    #
    # Usually, an FBA solution is not unique.  There can be many equivalent flux
    # distributions supporting the same objective value.  A method to decide for 
    # one out of these solutions is to compute the flux distribution minimizing 
    # the total absolute flux (MTF) but still supporting the objective value of 
    # the FBA solution.  At first, an objective value, for example calculated 
    # via FBA, is required
    #

    cat('\nMinimize total flux')
    cat('\n-------------------\n')

    cat('\nObtain an initial optimized value using FBA\n');
    fba <- optimizeProb(model, algorithm="fba");
    #
    cat('\nOptimized value of the objective function\n')
    print( mod_obj(fba) )
    #
    cat('\nNumber of variables (fba)\n')
    print( nvar(fluxdist(fba)) )
    #
    cat('\nUsing FBA value to obtain MTF\n')
    mtf <- optimizeProb(model, algorithm="mtf", wtobj = mod_obj(fba));
    #
    cat('\nValue of the objective function for the MTF:\n')
    print( lp_obj(mtf) )
    #
    cat('\nNumber of variables (mtf)\n')
    print( nvar(fluxdist(mtf)) )
    #
    cat('\nFlux distribution of the MTF solution\n')
    fl <- getFluxDist(mtf);

    cat('\nLength of the flux distribution:\n')
    print( length(fl) )

    cat('\nFlux distribution of ALL the reactions\n')
    print(getFluxDist(mtf,checkReactId(model, react_id(model))))
    #
    cat('\nFlux distribution of the exchange reactions\n')
    fd <- getFluxDist(mtf, exchngRxns);
    # 
    cat('\nNet flux of the exchange reactions\n');
    print( getNetFlux(fd) )
    #
    cat('Value of the objective function in the MTF model:\n')
    print( mod_obj(mtf) )
    #
    cat('\nSummary:\n')
    sum <- summaryOptsol(mtf, model)
    print(sum)
    # or only exchange reactions:
    #printExchange(sum, dense = TRUE)


    # Flux Variability Analysis (FVA)
    # -------------------------------
    # FBA only returns a single flux distribution that corresponds to maximal growth
    # under given growth conditions. However, alternate optimal solutions may exist
    # which correspond to maximal growth. FVA calculates the full range of numerical
    # values for each reaction flux within the network
    # 
    # The function fluxVar performs a flux variability analysis with a given model 
    # [Mahadevan and Schilling, 2003].  The minimum and maximum flux values for 
    # each reaction in the model are calculated, which still support a certain 
    # percentage of a given optimal functional state Z_opt
    #
    cat('\nFlux Variability Analysis (FVA)')
    cat('\n-------------------------------\n')

    opt <- fluxVar(model, percentage=80, verboseMode=0)
    #
    cat('\nSee plot to visualize minimum and maximum flux values for each reaction\n')
    png(pngFVAplot)

    plot(opt)

    dev.off()

    cat('\nSummary:\n')
    sum <- summaryOptsol(opt, model)
    print(sum)
    # or only exchange reactions:
    #printExchange(sum, dense = TRUE)


    # Robustness analysis
    # -------------------
    #
    # The function robAna performs a robustness analysis with a given model.  
    # The flux of a control reaction will be varied stepwise between the maximum 
    # and minimum value the flux of the control reaction can reach [Palsson, 2006].
    #
    #
    doRA=0
    if (doRA == 1) {
    cat('\nRobustness Analysis')
    cat('\n-------------------\n')

    cat('\nSee plot to visualize RA of EX_O2\n')
    # ideally this would be a 'for' loop
    opt <- robAna(model, ctrlreact='EX_O2', verboseMode=0)

    png(pngRAplot)

    plot(opt)

    dev.off()
    }

    # Phenotypic phase plane analysis
    # -------------------------------
    #
    # The function phpp performs a phenotypic phase plane analysis [Edwards et al.,
    # 2001, 2002b] with a given model.  The flux of two control reactions will be 
    # varied stepwise between a given maximum and minimum value.  
    #
    # We should select a meaningful pair of reactions here (e.g. BIOMASS/AMLB)
    doPPPA=0
    if (doPPPA == 1) {
    cat('\nPhenotypic Phase Plane Analysis')
    cat('\n-------------------------------\n')

    model_wo_glc <- changeUptake(model, off=FALSE)
    #
    opt <- phpp(model_wo_glc, 
	          ctrlreact = c('EX_mnl(e)', 'EX_o2(e)'),
	          redCosts = TRUE,
                  numP = 25,
                  verboseMode = 0)
    #
    png(pngPPPAplot)

    plot(opt)
    plot(opt, 'EX_glc(e)')
    plot(opt, 'EX_mnl(e)')
    plot(opt, 'EX_o2(e)')

    dev.off()
    }
    #

}

# Dynamic Flux Balance Analysis (DFBA)
# ------------------------------------
#
# Calculate concentrations of metabolites of exchange reactions at defined 
# time points given the initial concentrations. To accomplish this task 
# this function calls optimizeProb function to get the fluxes then update 
# the concentrations and the reaction boundaries ..etc
#
#
cat('\nDynamic Flux Balance Analysis (DFBA)')
cat('\n------------------------------------\n')

print(substrateRxns)
print(initConcentrations)

# Last minute changes:
# Make CO2 to be only exportable
#lowbnd(model)[react_id(model)=='EX_CO2']=0;
#uppbnd(model)[react_id(model)=='EX_CO2']=-10;

two.models <- list(model, model)
initBiomass <- initBiomass / 2
retOptSol <- TRUE

#
#df <- dynamicFBA.orig(
df <- multiDynFBA(
#df <- dynamicFBA.I(
                #model, 
                two.models,
                substrateRxns=substrateRxns, 
		initConcentrations=initConcentrations,
		initBiomass = initBiomass,
		timeStep = timeStep,
		nSteps = nSteps,
		exclUptakeRxns=exclUptakeRxns,
		retOptSol=retOptSol,
		fld=FALSE,
		verbose=3);

cat("\nWriting out final reports\n")
cat(  "-------------------------\n\n")

# esto no podemos usarlo porque df es una lista de objetos optsol_dynamicFBA, pero
# el plot definido en optsol_dynamicFBA está registrado para recibir un objeto 
# optsol_dynamicFBA, no una lista, por lo tanto "plot" cuando busque algo para
# plotear 'df' (lista) no encontrará el plot definido en optsol_dynamicFBA.R
#    x11()
#    plot(df, plotRxns=plotRxns);
#
# así que hay que volver a este enfoque.
#
# para que funcionara habría que definir una nueva clase, optsol_MDFBA, que trabajara
# con una lista de optsol_dynamicFBA, definida con un nombre propio (optsol_MDFBA) y 
# cuya función plot estuviera registrada para plotear objetos de ese tipo.
#
if (retOptSol == TRUE) {
    if (is(df, 'optsol_dynamicFBA')) {
        cat("Single model solution as optsol_dynamicFBA\n")
        # there was only one model, test old methods

    #    print(str(df))
    #    print(str(df@concentrationMatrix))
    #    print(str(df@excRxnNames))
    #    print(str(df@timeVec))
    #    print(str(df@biomassVec))
    #    print(str(df@all_fluxes))
        print('')
        print('')
        #print(str(df))

        x11()
        plot(df, plotRxns=plotRxns);

        png(pngDFBAplot, width=1000, height=1000)
        plot(df, plotRxns=plotRxns);
        dev.off()

    } else {
        if (is(df[[1]], "optsol_dynamicFBA")) {
            cat("Multiple model solution as list of optsol_dynamicFBA\n")
            # we requested retOptSol = TRUE and got a list of optsol_dynamicFBA objects

# to try later
#        nModels <- length(df)
#        # ensure all plots fit inside only one figure
#        previous.mfrow <- par("mfrow") # save so we can restore it at the end
#        plt.dim <- ceiling(sqrt(nModels))
#        par(mfrow = c(plt.dim, plt.dim))
#        png(pngDFBAplot, width=1000, height=1000)
            for (i in 1:length(df)) {
                #dfbasol <- df[i] # this is a list with one optsol element
                dfbasol <- df[[i]] # this is the optsol element in position i
    #            print(str(dfbasol))
               # there was only one model, test old methods
    #            print(str(dfbasol@concentrationMatrix))
    #            print(str(dfbasol@excRxnNames))
    #            print(str(dfbasol@timeVec))
    #            print(str(dfbasol@biomassVec))
    #            print(str(dfbasol@all_fluxes))
                print('')
                print('')
                #print(str(dfba))

                x11()
                plot(dfbasol, plotRxns=plotRxns);
                pngDFBAplot <- paste("out/", i,"_", modelName, "_DFBA.png", sep="")
                png(pngDFBAplot, width=1000, height=1000)
                plot(dfbasol, plotRxns=plotRxns);
                dev.off()
            }
# to try later
#        par(mfrow=previous.mfrow)
#	dev.off()
        } 
    }
} else {
    # if we did not get an optsol or a list of optsol objects, we got a 
    # list or a list of lists
    if ("timeVec" %in% names(df)) {
        # we had one model and got directly a list of optsol data
        print(df$nprob)
        print(df$ncols)
        print(df$nrows)
        print(str(df$all_fluxes))
        print(str(df$all_stat))
        print(str(df$concentrationMatrix))
        print(str(df$excRxnNames))
        print(str(df$timeVec))
        print(str(df$biomassVec))
        print('')
        print('')

        x11()
        plot(dfbasol$timeVec, dfbasol$biomassVec)
        png(pngDFBAplot, width=1000, height=1000)
        #my.plot(df, plotRxns=plotRxns);
        dev.off()  
    } else {
        # we got a list of optsol lists
        for (i in 1:length(df)) {
            dfbasol <- df[[i]]

            print(dfbasol$nprob)
            print(dfbasol$ncols)
            print(dfbasol$nrows)
            print(str(dfbasol$all_fluxes))
            print(str(dfbasol$all_stat))
            print(str(dfbasol$concentrationMatrix))
            print(str(dfbasol$excRxnNames))
            print(str(dfbasol$timeVec))
            print(str(dfbasol$biomassVec))
            print('')
            print('')

            x11()
            plot(dfbasol$timeVec, dfbasol$biomassVec)
            pngDFBAplot <- paste("out/", i,"_", modelName, "_DFBA.png", sep="")
            png(pngDFBAplot, width=1000, height=1000)
            #my.plot(df, plotRxns=plotRxns);
            dev.off()
        }
    }

}
cat('\nSummary:\n')
sum <- summaryOptsol(mtf, model)
print(sum)

# or only exchange reactions:

cat('\nSummary of the exchange reactions\n')
printExchange(sum, dense = TRUE)

cat('\nFlux distribution of ALL the reactions\n')
print(getFluxDist(mtf,checkReactId(model, react_id(model))))


sink()
close(logFile)
