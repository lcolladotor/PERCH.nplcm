
nplcm_fit = function (Mobs, Y, X, model_options, mcmc_options) 
{
###### define call bug function
        call.bugs <- function(data, inits, parameters, m.file, bugsmodel.dir = mcmc_options$bugsmodel.dir, 
                              winbugs.dir = mcmc_options$winbugs.dir, nitermcmc = mcmc_options$n.itermcmc, 
                              nburnin = mcmc_options$n.burnin, nthin = mcmc_options$n.thin, 
                              nchains = mcmc_options$n.chains, dic = FALSE, is.debug = mcmc_options$debugstatus, 
                              workd = result.folder, ...) 
        {
                m.file <- paste(bugsmodel.dir, m.file, sep = "")        ####### specify the model file used.
                
                f.tmp <- function()                                     ####### pass default par to bugs.
                {
                        gs <- bugs(data, inits, parameters, model.file = m.file, 
                                   working.directory = workd, n.chains = nchains, 
                                   n.iter = nitermcmc, n.burnin = nburnin, n.thin = nthin, 
                                   bugs.directory = winbugs.dir, DIC = dic, debug = is.debug, 
                                   ...)
                        gs                                              ####### MCMC par set, directory set.
                }
                bugs.try <- try(rst.bugs <- f.tmp(), silent = FALSE)    ###### error handling.
                if (class(bugs.try) == "try-error") {
                        rst.bugs <- NULL
                }
                rst.bugs
        }
####### call bug function end
        
        data_source <- c("yes", "no")[is.na(Mobs) + 1]
        names(data_source) <- names(Mobs)
        cat("Available measurements: \n")
        print(data_source)      ###### list(MBS = M_NPPCR,MSS = M_BCX,MGS = NA) => yes yes no => bronze + silver data

        if (length(model_options$M_use) != length(model_options$TPR_prior)) {
                                ###### M_use must be a subset of available measures
                stop("The number of measurement source(s) is different from\n              
                     the number of TPR prior option!\n              
                     Make them equal, and match with order!")
        }                       ###### M_use = c("BrS") => use Bronze data only


        Nd <- sum(Y == 1)       ###### number of cases
        Nu <- sum(Y == 0)       ###### number of controls

        model_data_source <- rep(NA, 3)
        names(model_data_source) <- c("MBS", "MSS", "MGS")
        model_data_source[1] <- c("no", "yes")["BrS" %in% model_options$M_use + 1]
        model_data_source[2] <- c("no", "yes")["SS" %in% model_options$M_use + 1]
        model_data_source[3] <- c("no", "yes")["GS" %in% model_options$M_use + 1]


        cat("Actual measurements used in the model: \n")
        print(model_data_source) ####### print out which set of data are used. e.g. Bronze only.


        cat("True positive rate (TPR) prior(s) for ", 
            c("MBS", "MSS", "MGS")[model_data_source == "yes"], "\n", 
            " is(are respectively) ", model_options$TPR_prior, "\n") ####### e.g. TPR_prior = c("noninformative")

        allowed_list <- model_options$allowed_list   ####### allowed list of pathogen/pathogen combinations infecting
        pathogen_list <- model_options$pathogen_list ####### complete list of pathogens
        
        Jfull <- length(pathogen_list)               
        Jallowed <- length(allowed_list)
        
        ######## ????
        template <- as.matrix(rbind(symb2I(allowed_list, pathogen_list), rep(0, Jfull)))
        ########## if allowed == pathogen, then this is a Identity Matrix plus a row of 0s.

############################################################        
####### when there is no subclass, not nested.##############
############################################################
        if (model_options$k_subclass == 1) 
        {
                cat("Number of subclasses: ", model_options$k_subclass, "\n")
                
        ##============= CASE 1.1: MBS only - NPPCR model ==================##
                if (model_data_source[1] == "yes" & 
                    model_data_source[2] == "no" & 
                    model_data_source[3] == "no") #"MBS" "MSS" "MGS"#
                {
                        MBS.case <- Mobs$MBS[Y == 1, ]
                        MBS.ctrl <- Mobs$MBS[Y == 0, ]
                        MBS <- as.matrix(rbind(MBS.case, MBS.ctrl))
                        
                ##----------------- Set non-informative priors ---------------##
                        # Etiology prior #
                        if (model_options$Eti_prior == "overall_uniform") {
                                alpha <- rep(1, Jfull) ##### Dirichlet par: (1,1,...,1)
                        }
                        # Sensitivity prior#
                        if (model_options$TPR_prior == "noninformative") {
                                alphaB <- rep(1, Jfull)
                                betaB <- rep(1, Jfull) ###### Beta par: (1,1)
                        }
                        # FPR priors are all set as dbeta(1,1)
                
                ##----------- No informative priors ?? ------------------------##
                
                        
                ##---------------- Define a BrS data specific , non-nested bug function --------------##
                        mybugs <- function(...) 
                        {
                                inits <- function() 
                                {
                                        list(thetaBS = rbeta(Jfull, 1, 1), psiBS = rbeta(Jfull, 
                                                                                         1, 1))
                                }
                                
                                data <- c("Nd", "Nu", "Jfull", "Jallowed", "alpha", 
                                          "template", "MBS", "alphaB", "betaB")
                                
                        ####### individual.pred: whether to perform individual prediction
                        ####### ppd: whether to perform posterior predictive checking
                        ####### lead to different parameter space.
                                if (mcmc_options$individual.pred == FALSE & 
                                            mcmc_options$ppd == TRUE) 
                                {
                                        parameters <- c("thetaBS", "psiBS", "pEti", "MBS.new")
                                ####### sensitivity, specificity, etio fraction, predicted Y?
                                }
                                else if (mcmc_options$individual.pred == TRUE & 
                                                 mcmc_options$ppd == TRUE) 
                                {
                                        parameters <- c("thetaBS", "psiBS", "pEti", "Icat", "MBS.new")
                                ####### sensitivity, specificity, etio fraction, predicted.etio.category?, Predicted Y?
                                }
                                else if (mcmc_options$individual.pred == TRUE & 
                                                 mcmc_options$ppd == FALSE) 
                                {
                                        parameters <- c("thetaBS", "psiBS", "pEti", "Icat")
                                }
                                else if (mcmc_options$individual.pred == FALSE & 
                                                 mcmc_options$ppd == FALSE) 
                                {
                                        parameters <- c("thetaBS", "psiBS", "pEti")
                                ####### only do parameter estimation. 
                                }
                                rst.bugs <- call.bugs(data, inits, parameters, 
                                                      ...)
                                rst.bugs
                        }
                ######################## Function definition ENDs ################################
                
                        if (mcmc_options$ppd == TRUE) 
                        {
                                gs <- mybugs("model_plcm_brsonly_ppd.bug")
#### BUGS model 1 ------------- not nested, BrS only, with posterior predictive checking
                        }
                        else 
                        {
                                gs <- mybugs("model_plcm_brsonly.bug")
#### BUGS model 2 ------------- not nested, BrS only, with out posterior predictive checking
                        }
                }
                
        ##============= CASE 1.2: MBS and MSS - NPPCR & BCX model ==================##
                if (model_data_source[1] == "yes" & 
                            model_data_source[2] == "yes" & 
                            model_data_source[3] == "no") 
                {
                ##----- BrS part : case & control
                        MBS.case <- Mobs$MBS[Y == 1, ]
                        MBS.ctrl <- Mobs$MBS[Y == 0, ]
                        MBS <- as.matrix(rbind(MBS.case, MBS.ctrl))
                        
                ##----- SS part : perfect specificity => no control needed (all assumed test negative)
                        MSS.case <- Mobs$MSS[Y == 1, ]
                        MSS.case <- as.matrix(MSS.case)
                        SS_index <- which(colMeans(is.na(MSS.case)) < 0.9) ##### missing rates < 90%
                        JSS <- length(SS_index)         ##### number of pathogens with SS data
                        MSS <- MSS.case[, SS_index]     ###### NPPCR case data
                
                ##----------------- Set non-informative Etiology prior ---------------##                
                        if (model_options$Eti_prior == "overall_uniform") {
                                alpha <- rep(1, Jfull)
                        }
                
                ##------------- Set BrS sensitivity  priors ---------------##
                        if (model_options$TPR_prior[1] == "noninformative") {
                                alphaB <- rep(1, Jfull)
                                betaB <- rep(1, Jfull)
                        }
                
                ######### Informative ############
                        else if (model_options$TPR_prior[1] == "informative") 
                        {
                                alphaB <- rep(NA, Jfull)
                                betaB <- rep(NA, Jfull)
                                temp_cat_ind <- sapply(model_options$pathogen_list, 
                                                       function(path) which(model_options$pathogen_cat$X == path))
                                temp_cat <- model_options$pathogen_cat[temp_cat_ind, ] ###### get pathogen type: V/B/F
                                
                                for (t in 1:nrow(temp_cat)) ##### apply to each pathogen 
                                {
                                ####### Bacteria: dbeta(1,1)
                                        if (temp_cat$pathogen_type[t] == "B") 
                                        {
                                                alphaB[t] <- 1
                                                betaB[t] <- 1
                                        }
                                ####### Virus or Fungus: dbeta(6,2) why?
                                        if (temp_cat$pathogen_type[t] == "V" | temp_cat$pathogen_type[t] == "F") 
                                        {
                                                alphaB[t] <- 6
                                                betaB[t] <- 2
                                        }
                                }
                        }
                
                ##------------- Set SS sensitivity priors ---------------##
                        if (model_options$TPR_prior[2] == "noninformative") {
                                alphaS <- rep(1, JSS)
                                betaS <- rep(1, JSS)
                        }
                ######## Informative #########
                        else if (model_options$TPR_prior[2] == "informative") 
                        {
                                alphaS <- rep(NA, JSS)
                                betaS <- rep(NA, JSS)
                                temp_cat_ind <- sapply(model_options$pathogen_list, 
                                                       function(path) which(model_options$pathogen_cat$X == path))
                                temp_cat <- model_options$pathogen_cat[temp_cat_ind, ] ##### should be of the same length of  pathogen list
                                
                                for (t in 1:JSS) 
                                {
                                ####### generate a beta prior matching the TPR range with 95% CI:  e.g. 0.05-0.15
                                        temp_param <- beta_parms_from_quantiles(c(0.05, 0.15), 
                                                                                p = c(0.025, 0.975), plot = FALSE)
                                        alphaS[t] <- temp_param$a
                                        betaS[t] <- temp_param$b
                                }
                        }
                ##---------------- Define a BrS + SS data specific , non-nested bug function --------------##
                        mybugs <- function(...) {
                                inits <- function() {
                                        list(thetaBS = rbeta(Jfull, 1, 1), psiBS = rbeta(Jfull, 1, 1))
                                }
                                data <- c("Nd", "Nu", "Jfull", "Jallowed", "alpha", 
                                          "template", "MBS", "JSS", "MSS", "alphaB", 
                                          "betaB", "alphaS", "betaS")
                                
                                if (mcmc_options$individual.pred == FALSE & mcmc_options$ppd == TRUE) 
                                {
                                        parameters <- c("thetaBS", "psiBS", "pEti", "thetaSS", "MBS.new")
                                }
                                else if (mcmc_options$individual.pred == TRUE & mcmc_options$ppd == TRUE) 
                                {
                                        parameters <- c("thetaBS", "psiBS", "pEti", "thetaSS", "Icat", "MBS.new")
                                }
                                else if (mcmc_options$individual.pred == TRUE & mcmc_options$ppd == FALSE) 
                                {
                                        parameters <- c("thetaBS", "psiBS", "pEti", "thetaSS", "Icat")
                                }
                                else if (mcmc_options$individual.pred == FALSE & mcmc_options$ppd == FALSE) 
                                {
                                        parameters <- c("thetaBS", "psiBS", "pEti",  "thetaSS")
                                }
                                ###### passes diffenret parameter spaces and inits.
                                rst.bugs <- call.bugs(data, inits, parameters, 
                                                      ...)
                                rst.bugs
                        }
                ############################ function ENDS ###########################################
       
#### BUGS model 3 ------------- not nested, BrS + SS, with posterior predictive checking
                        if (mcmc_options$ppd == TRUE) {
                                gs <- mybugs("model_plcm_ppd.bug")
                        }
#### BUGS model 4 ------------- not nested, BrS + SS, without posterior predictive checking
                        else {
                                gs <- mybugs("model_plcm.bug")
                        }
                }
        }

############################################################        
####### when there are subclasses, nested plcm.#############
############################################################
        else {
                cat("Number of subclasses: ", model_options$k_subclass, 
                    "\n")
                if (model_data_source[1] == "yes" & model_data_source[2] == 
                            "no" & model_data_source[3] == "no") {
                        MBS.case <- Mobs$MBS[Y == 1, ]
                        MBS.ctrl <- Mobs$MBS[Y == 0, ]
                        MBS <- as.matrix(rbind(MBS.case, MBS.ctrl))
                        K <- model_options$k_subclass
                        if (model_options$Eti_prior == "overall_uniform") {
                                alpha <- rep(1, Jfull)
                        }
                        if (model_options$TPR_prior == "noninformative") {
                                alphaB <- rep(1, Jfull)
                                betaB <- rep(1, Jfull)
                        }
                        mybugs <- function(...) {
                                inits <- function() {
                                        list(pEti = rep(1/Jallowed, Jallowed), r0 = c(rep(0.5, 
                                                                                          K - 1), NA), r1 = cbind(matrix(rep(0.5, Jfull * 
                                                                                                                                     (K - 1)), nrow = Jfull, ncol = K - 1), rep(NA, 
                                                                                                                                                                                Jfull)), alphadp0 = 1)
                                }
                                data <- c("Nd", "Nu", "Jfull", "Jallowed", "alpha", 
                                          "template", "K", "MBS", "alphaB", "betaB")
                                if (mcmc_options$individual.pred == FALSE & mcmc_options$ppd == 
                                            TRUE) {
                                        parameters <- c("pEti", "Lambda", "Eta", "alphadp0", 
                                                        "MBS.new", "ThetaBS.marg", "PsiBS.marg", 
                                                        "PsiBS.case", "ThetaBS", "PsiBS")
                                }
                                else if (mcmc_options$individual.pred == TRUE & 
                                                 mcmc_options$ppd == TRUE) {
                                        parameters <- c("pEti", "Lambda", "Eta", "alphadp0", 
                                                        "Icat", "MBS.new", "ThetaBS.marg", "PsiBS.marg", 
                                                        "PsiBS.case", "ThetaBS", "PsiBS")
                                }
                                else if (mcmc_options$individual.pred == TRUE & 
                                                 mcmc_options$ppd == FALSE) {
                                        parameters <- c("pEti", "Lambda", "Eta", "alphadp0", 
                                                        "Icat", "ThetaBS.marg", "PsiBS.marg", "PsiBS.case", 
                                                        "ThetaBS", "PsiBS")
                                }
                                else if (mcmc_options$individual.pred == FALSE & 
                                                 mcmc_options$ppd == FALSE) {
                                        parameters <- c("pEti", "Lambda", "Eta", "alphadp0", 
                                                        "ThetaBS.marg", "PsiBS.marg", "PsiBS.case", 
                                                        "ThetaBS", "PsiBS")
                                }
                                rst.bugs <- call.bugs(data, inits, parameters, 
                                                      ...)
                                rst.bugs
                        }
                        if (mcmc_options$ppd == TRUE) {
                                gs <- mybugs("model_nplcm_brsonly_ppd.bug")
                        }
                        else {
                                gs <- mybugs("model_nplcm_brsonly.bug")
                        }
                }
                if (model_data_source[1] == "yes" & model_data_source[2] == 
                            "yes" & model_data_source[3] == "no") {
                        MBS.case <- Mobs$MBS[Y == 1, ]
                        MBS.ctrl <- Mobs$MBS[Y == 0, ]
                        MBS <- as.matrix(rbind(MBS.case, MBS.ctrl))
                        MSS.case <- Mobs$MSS[Y == 1, ]
                        MSS.case <- as.matrix(MSS.case)
                        SS_index <- which(colMeans(is.na(MSS.case)) < 0.9)
                        JSS <- length(SS_index)
                        MSS <- MSS.case[, SS_index]
                        K <- model_options$k_subclass
                        if (model_options$Eti_prior == "overall_uniform") {
                                alpha <- rep(1, Jfull)
                        }
                        if (model_options$TPR_prior[1] == "noninformative") {
                                alphaB <- rep(1, Jfull)
                                betaB <- rep(1, Jfull)
                        }
                        else if (model_options$TPR_prior[1] == "informative") {
                                alphaB <- rep(NA, Jfull)
                                betaB <- rep(NA, Jfull)
                                temp_cat_ind <- sapply(model_options$pathogen_list, 
                                                       function(path) {
                                                               which(model_options$pathogen_cat$X == path)
                                                       })
                                temp_cat <- model_options$pathogen_cat[temp_cat_ind, 
                                                                       ]
                                for (t in 1:nrow(temp_cat)) {
                                        if (temp_cat$pathogen_type[t] == "B") {
                                                alphaB[t] <- 1
                                                betaB[t] <- 1
                                        }
                                        if (temp_cat$pathogen_type[t] == "V" | temp_cat$pathogen_type[t] == 
                                                    "F") {
                                                alphaB[t] <- 6
                                                betaB[t] <- 2
                                        }
                                }
                        }
                        if (model_options$TPR_prior[2] == "noninformative") {
                                alphaS <- rep(1, JSS)
                                betaS <- rep(1, JSS)
                        }
                        else if (model_options$TPR_prior[2] == "informative") {
                                alphaS <- rep(NA, JSS)
                                betaS <- rep(NA, JSS)
                                temp_cat_ind <- sapply(model_options$pathogen_list, 
                                                       function(path) {
                                                               which(model_options$pathogen_cat$X == path)
                                                       })
                                temp_cat <- model_options$pathogen_cat[temp_cat_ind, 
                                                                       ]
                                for (t in 1:JSS) {
                                        temp_param <- beta_parms_from_quantiles(c(0.05, 
                                                                                  0.15), p = c(0.025, 0.975), plot = FALSE)
                                        alphaS[t] <- temp_param$a
                                        betaS[t] <- temp_param$b
                                }
                        }
                        
                        mybugs <- function(...) {
                                inits <- function() 
                                {
                                        list(pEti = rep(1/Jallowed, Jallowed), 
                                             r0 = c(rep(0.5, K - 1), NA), 
                                             r1 = cbind(matrix(rep(0.5, Jfull * (K - 1)),
                                                               nrow = Jfull, ncol = K - 1), 
                                                        rep(NA, Jfull)), alphadp0 = 1)
                                }
                                
                                data <- c("Nd", "Nu", "Jfull", "Jallowed", "alpha", 
                                          "template", "K", "JSS", "MSS", "MBS", "alphaB", 
                                          "betaB", "alphaS", "betaS")
                                
                                if (mcmc_options$individual.pred == FALSE & mcmc_options$ppd == TRUE) 
                                {
                                        parameters <- c("pEti", "Lambda", "Eta", "alphadp0", 
                                                        "MBS.new", "ThetaBS.marg", "PsiBS.marg", 
                                                        "PsiBS.case", "ThetaBS", "PsiBS", "thetaSS")
                                }
                                else if (mcmc_options$individual.pred == TRUE &  mcmc_options$ppd == TRUE) 
                                {
                                        parameters <- c("pEti", "Lambda", "Eta", "alphadp0", 
                                                        "Icat", "MBS.new", "ThetaBS.marg", "PsiBS.marg", 
                                                        "PsiBS.case", "ThetaBS", "PsiBS", "thetaSS")
                                }
                                else if (mcmc_options$individual.pred == TRUE & mcmc_options$ppd == FALSE) 
                                {
                                        parameters <- c("pEti", "Lambda", "Eta", "alphadp0", 
                                                        "Icat", "ThetaBS.marg", "PsiBS.marg", "PsiBS.case", 
                                                        "ThetaBS", "PsiBS", "thetaSS")
                                }
                                else if (mcmc_options$individual.pred == FALSE & mcmc_options$ppd == FALSE) 
                                {
                                        parameters <- c("pEti", "Lambda", "Eta", "alphadp0", 
                                                        "ThetaBS.marg", "PsiBS.marg", "PsiBS.case", 
                                                        "ThetaBS", "PsiBS", "thetaSS")
                                }
                                rst.bugs <- call.bugs(data, inits, parameters, 
                                                      ...)
                                rst.bugs
                        }
                ########################## function ENDS ########################
                        if (mcmc_options$ppd == TRUE) 
                        {
                                gs <- mybugs("model_nplcm_ppd.bug")
                        }
                        else {
                                gs <- mybugs("model_nplcm.bug")
                        }
                }
        }
        return(gs)
}

extract_data_raw = function(Pathogen, Specimen, Test, X, Xval, MeasDir, PathCatDir, extra_covariates = NULL, silent = TRUE) 
{
        pathogen_type = read.csv(PathCatDir)
        rownames(pathogen_type) = pathogen_type[, 1]
        typeOrder = order(pathogen_type[Pathogen, 2])
        Pathogen = Pathogen[typeOrder]
        if (!silent) {
                cat("Pathogens included:")
                print(t(rbind(Pathogen, as.character(pathogen_type[Pathogen, 
                                                                   2]))))
        }
        datraw = read.csv(MeasDir)
        delete_start_with = function(s, vec) {
                ind = grep(s, substring(vec, 1, nchar(s)))
                old = vec[ind]
                vec[ind] = substring(old, nchar(s) + 1)
                return(vec)
        }
        cleanName = delete_start_with("X_", names(datraw))
        dat0 = datraw
        colnames(dat0) = cleanName
        indX = 1:nrow(dat0)
        for (j in 1:length(X)) {
                indX = indX[which(dat0[indX, X[j]] == Xval[j] & !is.na(dat0[indX, 
                                                                            X[j]]))]
        }
        dat = dat0[indX, ]
        pstGrid = apply(expand.grid(paste(Pathogen, "_", sep = ""), 
                                    Specimen, Test), 1, paste, collapse = "")
        stGrid = apply(expand.grid(Specimen, Test), 1, paste, collapse = "")
        pstTable = matrix(NA, length(Pathogen), length(stGrid))
        colnames(pstTable) = stGrid
        rownames(pstTable) = Pathogen
        for (i in 1:length(stGrid)) {
                for (j in 1:length(Pathogen)) {
                        tempName = paste(Pathogen[j], stGrid[i], sep = "_")
                        if (tempName %in% cleanName) {
                                if (sum(is.na(dat[, tempName])) != nrow(dat)) {
                                        pstTable[j, i] = TRUE
                                }
                        }
                }
        }
        notindata = which(rowSums(!is.na(pstTable)) == 0)
        if (length(notindata) > 0) {
                stop(Pathogen[notindata], "can't be found in the dataset!", 
                     "\n")
        }
        naColumns = which(colSums(is.na(pstTable)) == length(Pathogen))
        if (length(naColumns) == length(stGrid)) {
                stop("No test has available results on selected pathogens! Try other pathogens.", 
                     "\n")
        }
        else {
                if (length(naColumns) == 0) {
                        actualpstTable = pstTable
                }
                else {
                        actualpstTable = pstTable[, -naColumns, drop = FALSE]
                }
                resdat = list()
                for (j in 1:ncol(actualpstTable)) {
                        if (sum(is.na(actualpstTable[, j])) == 0) {
                                tempnm = paste(Pathogen, colnames(actualpstTable)[j], 
                                               sep = "_")
                                resdat[[j]] = dat[, tempnm]
                        }
                        else {
                                dftemp = as.data.frame(matrix(NA, nrow = nrow(dat), 
                                                              ncol = nrow(actualpstTable)))
                                colnames(dftemp) = paste(Pathogen, colnames(actualpstTable)[j], 
                                                         sep = "_")
                                dftemp[, !is.na(actualpstTable[, j])] = dat[, 
                                                                            paste(Pathogen[!is.na(actualpstTable[, j])], 
                                                                                  colnames(actualpstTable)[j], sep = "_")]
                                resdat[[j]] = dftemp
                        }
                }
                if (!is.null(extra_covariates)) {
                        for (i in seq_along(extra_covariates)) {
                                if (!extra_covariates[i] %in% colnames(dat)) {
                                        stop(extra_covariates[i], " is not in the data set!", 
                                             "\n")
                                }
                                else {
                                        resdat[[ncol(actualpstTable) + i]] = dat[, 
                                                                                 extra_covariates[i]]
                                }
                        }
                }
                curr_len = length(resdat)
                for (j in seq_along(X)) {
                        resdat[[curr_len + j]] = rep(Xval[j], nrow(dat))
                }
                names(resdat) = c(colnames(actualpstTable), extra_covariates, 
                                  X)
                return(resdat)
        }
}


################################################################################################################
# ------------- model_plcm_ssonly.bug ---------------------- #

data <- c("Nd", "Nu", "Jfull", "Jallowed", "alpha", 
          "template", "MBS", "JSS", "MSS", "alphaB", 
          "betaB", "alphaS", "betaS", "Jss.only") #################### 

parameters <- c("thetaBS", "psiBS", "pEti",  "thetaSS")

MODEL{
##=============== Likelihood: L(M;pi,theta,psi)  ==============##
        # BrS Measure data #
        for (k in 1:(Nd+Nu)) # for all subjects
        {
                for (j in (Jss.only+1):Jfull) # for all pathogens
                {
                        ind[k,j]<-equals(1,template[Icat[k],j]) ##### latent lung state indicator 1(Ilung_i == j)
                        MBS[k,j] ~ dbern(mu_bs[k,j])            ##### BS measure data likelihood
                        mu_bs[k,j]<-ind[k,j]*thetaBS[j]+(1-ind[k,j])*psiBS[j]
                        
                }
        }
        
        # SS Measure data #
        for (k in 1:Nd) ###### cases
        { 
                for (j in 1:JSS) ##### for pathogens with SS data
                {
                        MSS[k,j] ~ dbern(mu_ss[k,j])            ##### SS measure data likelihood
                        mu_ss[k,j]<-ind[k,j]*thetaSS[j]+(1-ind[k,j])*psiSS[j]
                }
        }
        
        # Latent variable: for cases
        for (k in 1:Nd) 
        {
                Icat[k] ~ dcat(pEti[1:Jallowed]) ###### latent Ilung ~ uniform multinomial( pi )
        }
        
        # Latent variable: for controls
        for (k in (Nd+1):(Nd+Nu)) 
        {
                Icat[k]<-Jallowed+1              ###### fixed category J+1 as control
        }
        
##===============  priors : Pi, theta, psi =================##
        pEti[1:Jallowed]~ddirch(alpha[])         ###### pi ~ dirichlet prior (alpha)        
        
        # bronze-standard measurement characteristics:
        for (j in (Jss.only+1):Jfull)
        {
                thetaBS[j]~dbeta(alphaB[j],betaB[j]) ###### sensitivity prior BS
                psiBS[j]~dbeta(1,1)                  ###### FPR prior   
        }
        
        
        # silver-standard measurement characteristics:
        for (j in 1:JSS)
        {
                thetaSS[j]~dbeta(alphaS[j],betaS[j]) ###### sensitivity prior SS
                psiSS[j]<-0                          ###### perfect specificity
        }
        
}#end of model


################################################################################################################
# ------------- model_nplcm.bug ---------------------- #

#######################################################
## Model for nested partial latent class model
#######################################################

MODEL{#BEGIN MODEL
        ##case measurements:
        for (i in 1:Nd){
                for (j in 1:Jfull){
                        ind[i,j] <- equals(1,template[Icat[i],j])
                        MBS[i,j]~dbern(mu_bs.bound[i,j])
                        mu_bs.bound[i,j]<-max(0.000001,min(0.999999,mu_bs[i,j]))
                        mu_bs[i,j]<-PR_BS[i,j,Z[i]]
                        
                        for (s in 1:K){
                                PR_BS[i,j,s]<-PsiBS.cut[j,s]*(1-ind[i,j])+ThetaBS[j,s]*ind[i,j]
                        }
                }
        }
        
        ## cut the feedback from case model to FPR:
        for (j in 1:Jfull){
                for (s in 1:K){
                        PsiBS.cut[j,s]<-cut(PsiBS[j,s])
                }
        }
        
        ## control measurements
        for (i in (Nd+1):(Nd+Nu)){
                for (j in 1:Jfull){
                        MBS[i,j]~dbern(mu_bs.bound[i,j])
                        mu_bs.bound[i,j] <-max(0.000001,min(0.999999,mu_bs[i,j]))
                        mu_bs[i,j]<-PsiBS[j,Z[i]]
                }
        }
        
        for (i in 1:Nd){
                for (j in 1:JSS){
                        MSS[i,j] ~ dbern(mu_ss[i,j])
                        mu_ss[i,j]<-ind[i,j]*thetaSS[j]+(1-ind[i,j])*psiSS[j]
                }
        }
        
        for (i in 1:Nd){
                Z[i]~dcat(Eta[Icat[i],1:K])
                Icat[i] ~ dcat(pEti[1:Jallowed])
        }
        
        ######################
        ##etiology prior
        ######################
        pEti[1:Jallowed]~ddirch(alpha[])
        for (i in (Nd+1):(Nd+Nu)){
                Z[i]~dcat(Lambda[1:K])
        }
        
        ####################################
        ### stick-breaking specification 2
        ####################################
        Lambda0[1]<-r0[1]
        r0[K]<-1
        for(j in 2:K) {Lambda0[j]<-r0[j]*(1-r0[j-1])*Lambda0[j-1]/r0[j-1]}
        for(k in 1:K-1){
                r0[k]~dbeta(1,alphadp0)I(0.000001,0.999999)
        }
        
        for (k in 1:K-1){Lambda[k]<-max(0.000001,min(0.999999,Lambda0[k]))}
        Lambda[K]<-1-sum(Lambda[1:(K-1)])
        
        for (s in 1:Jfull){
                Eta0[s,1]<-r1[s,1]
                r1[s,K]<-1
                for(j in 2:K) {Eta0[s,j]<-r1[s,j]*(1-r1[s,j-1])*Eta0[s,j-1]/r1[s,j-1]}
                for(k in 1:K-1){
                        r1[s,k]~dbeta(1,alphadp0)I(0.000001,0.999999)
                }
        }
        
        for (s in 1:Jfull){
                for (k in 1:K-1){Eta[s,k]<-max(0.000001,min(0.999999,Eta0[s,k]))}
                Eta[s,K]<-1-sum(Eta[s,1:(K-1)])
        }
        
        alphadp0~dgamma(.25,.25)I(0.001,20)
        
        #########################
        ## priors on TPR and FPR:
        #########################
        
        for (j in 1:Jfull){
                for (s in 1:K){
                        PsiBS[j,s]~dbeta(1,1)
                        #ThetaBS[j,s]~dbeta(1,1)
                        ThetaBS[j,s]~dbeta(alphaB[j],betaB[j])
                }
                ThetaBS.marg[j]<-inprod2(ThetaBS[j,1:K],Eta[j,1:K])
                PsiBS.marg[j]<-inprod2(PsiBS[j,1:K],Lambda[1:K])
                
                for (l in 1:Jfull){
                        PsiBS.case[j,l]<-inprod2(PsiBS[j,1:K],Eta[l,1:K])
                }
        }
        
        # silver-standard measurement characteristics:
        for (j in 1:JSS){
                thetaSS[j]~dbeta(alphaS[j],betaS[j])
                psiSS[j]<-0
        }
        
        
}#END MODEL








