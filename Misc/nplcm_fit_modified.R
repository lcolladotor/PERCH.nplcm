nplcm_fit2<-function(Mobs,Y,X,model_options,mcmc_options){#BEGIN function
        #define the generic function to call WinBUGS:
        if (.Platform$OS.type != "windows"){
                call.bugs <- function(data, inits, parameters,m.file,
                                      bugsmodel.dir = mcmc_options$bugsmodel.dir,
                                      winbugs.dir   = mcmc_options$winbugs.dir,
                                      nitermcmc     = mcmc_options$n.itermcmc,
                                      nburnin       = mcmc_options$n.burnin,
                                      nthin         = mcmc_options$n.thin,
                                      nchains       = mcmc_options$n.chains,
                                      dic = FALSE, is.debug = mcmc_options$debugstatus,
                                      workd= mcmc_options$result.folder,
                                      WINE = mcmc_options$WINE,
                                      WINEPATH = mcmc_options$WINEPATH,...) {
                        
                        m.file <- paste(bugsmodel.dir, m.file, sep="");
                        f.tmp <- function() {
                                ##winbugs
                                gs <- bugs(data, inits, parameters,
                                           model.file = m.file,
                                           working.directory=workd,
                                           n.chains = nchains,
                                           n.iter   = nitermcmc,
                                           n.burnin = nburnin,
                                           n.thin   = nthin,
                                           bugs.directory=winbugs.dir,
                                           DIC=dic,
                                           debug=is.debug,...);
                                
                                gs;
                        }
                        
                        bugs.try  <- try(rst.bugs <- f.tmp(), silent=FALSE);
                        if (class(bugs.try) == "try-error") {
                                rst.bugs <- NULL;
                        }
                        rst.bugs;
                }
        }else {
                #define the generic function to call WinBUGS:
                call.bugs <- function(data, inits, parameters,m.file,
                                      bugsmodel.dir = mcmc_options$bugsmodel.dir,
                                      winbugs.dir   = mcmc_options$winbugs.dir,
                                      nitermcmc     = mcmc_options$n.itermcmc,
                                      nburnin       = mcmc_options$n.burnin,
                                      nthin         = mcmc_options$n.thin,
                                      nchains       = mcmc_options$n.chains,
                                      dic = FALSE, is.debug = mcmc_options$debugstatus,
                                      workd= mcmc_options$result.folder,...) {
                        
                        m.file <- paste(bugsmodel.dir, m.file, sep="");
                        f.tmp <- function() {
                                ##winbugs
                                gs <- bugs(data, inits, parameters,
                                           model.file = m.file,
                                           working.directory=workd,
                                           n.chains = nchains,
                                           n.iter   = nitermcmc,
                                           n.burnin = nburnin,
                                           n.thin   = nthin,
                                           bugs.directory=winbugs.dir,
                                           DIC=dic,
                                           debug=is.debug,...);
                                
                                gs;
                        }
                        
                        bugs.try  <- try(rst.bugs <- f.tmp(), silent=FALSE);
                        if (class(bugs.try) == "try-error") {
                                rst.bugs <- NULL;
                        }
                        rst.bugs;
                }
                
        }
        
        #check sources of measurement data:
        data_source        <- c("yes","no")[is.na(Mobs)+1]
        names(data_source) <- names(Mobs)
        cat("Available measurements: \n")
        print(data_source)
        
        #compatibility checking:
        if (length(model_options$M_use)!=length(model_options$TPR_prior)){
                stop("The number of measurement source(s) is different from
               the number of TPR prior option!
               Make them equal, and match with order!")
        }
        
        #some data preparation:
        Nd <- sum(Y==1)
        Nu <- sum(Y==0)
        
        model_data_source <- rep(NA,3)
        names(model_data_source) <- c("MBS","MSS","MGS")
        model_data_source[1] <- c("no","yes")["BrS"%in%model_options$M_use+1]
        model_data_source[2] <- c("no","yes")["SS"%in%model_options$M_use+1]
        model_data_source[3] <- c("no","yes")["GS"%in%model_options$M_use+1]
        
        cat("Actual measurements used in the model: \n")
        print(model_data_source)
        
        cat("True positive rate (TPR) prior(s) for ",
            c("MBS","MSS","MGS")[model_data_source=="yes"],"\n",
            " is(are respectively) ", model_options$TPR_prior,"\n")
        
        allowed_list <- model_options$allowed_list
        pathogen_list <-model_options$pathogen_list
        
        Jfull         <- length(pathogen_list)
        Jallowed      <- length(allowed_list)
        
        Jss.only = sum(colMeans(is.na(Mobs$MBS))) # <--------------- Detian: number of SS only pathogen
        
        template <-  as.matrix(rbind(symb2I(allowed_list,
                                            pathogen_list),rep(0,Jfull)))
        
        # BEGIN fit model -------------------------------------------------------------
        if (model_options$k_subclass ==1){
                cat("Number of subclasses: ", model_options$k_subclass,"\n")
                #plcm: conditional independence model
                #plcm-BrS only:
                if (model_data_source[1]=="yes" &
                            model_data_source[2]=="no" &
                            model_data_source[3]=="no"){
                        
                        if (Jss.only>0) stop("Please remove SS only data to run BS model.") # <---- Detian: 
                        
                        MBS.case <- Mobs$MBS[Y==1,]
                        MBS.ctrl <- Mobs$MBS[Y==0,]
                        MBS      <- as.matrix(rbind(MBS.case,MBS.ctrl))
                        
                        if (model_options$Eti_prior=="overall_uniform"){
                                alpha    <-  rep(1,Jfull)
                        }
                        TPR_prior_list <- TPR_prior_set(model_options)
                        
                        alphaB <- TPR_prior_list$alphaB
                        betaB <- TPR_prior_list$betaB
                        
                        mybugs <- function(...){
                                inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                                              psiBS   = rbeta(Jfull,1,1))};
                                data       <- c("Nd","Nu","Jfull","Jallowed","alpha",
                                                "template","MBS","alphaB","betaB");
                                
                                if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
                                        parameters <- c("thetaBS","psiBS","pEti","MBS.new");
                                        
                                } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
                                        parameters <- c("thetaBS","psiBS","pEti","Icat","MBS.new");
                                        
                                } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
                                        parameters <- c("thetaBS","psiBS","pEti","Icat");
                                        
                                } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
                                        parameters <- c("thetaBS","psiBS","pEti");
                                        
                                }
                                rst.bugs   <- call.bugs(data, inits, parameters,...);
                                rst.bugs
                        }
                        
                        if (mcmc_options$ppd==TRUE){
                                gs <- mybugs("model_plcm_brsonly_ppd.bug")
                        } else {
                                gs <- mybugs("model_plcm_brsonly.bug")
                        }
                }
                
                ##plcm-BrS + SS:
                if (model_data_source[1]=="yes" &
                            model_data_source[2]=="yes" &
                            model_data_source[3]=="no"){
                        
                        MBS.case <- Mobs$MBS[Y==1,]
                        MBS.ctrl <- Mobs$MBS[Y==0,]
                        MBS      <- as.matrix(rbind(MBS.case,MBS.ctrl))
                        
                        MSS.case <- Mobs$MSS[Y==1,]
                        MSS.case <- as.matrix(MSS.case)
                        
                        SS_index <- which(colMeans(is.na(MSS.case))<0.9)#.9 is arbitrary; any number <1 will work.
                        JSS      <- length(SS_index)
                        MSS      <- MSS.case[,SS_index]
                        
                        if (model_options$Eti_prior=="overall_uniform"){
                                alpha    <-  rep(1,Jfull)
                        }
                        
                        TPR_prior_list <- TPR_prior_set(model_options)
                        
                        alphaB <- TPR_prior_list$alphaB
                        betaB <- TPR_prior_list$betaB
                        alphaS <- TPR_prior_list$alphaS
                        betaS <- TPR_prior_list$betaS
                        
                        data.in = c("Nd","Nu","Jfull","Jallowed","alpha","template",  ##########
                          "MBS","JSS","MSS",                                          # edited #
                          "alphaB","betaB","alphaS","betaS")                          ##########
                        
                        if(Jss.only>0)          
                        {
                                data.in = c(data.in,"Jss.only") # <------------- Detian: incoporate new variable          
                        }
                        
                        mybugs <- function(...){
                                inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                                              psiBS   = rbeta(Jfull,1,1))};
                                data       <- data.in
                                
                                if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
                                        parameters <- c("thetaBS","psiBS","pEti","thetaSS","MBS.new");
                                        
                                } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
                                        parameters <- c("thetaBS","psiBS","pEti","thetaSS","Icat","MBS.new")
                                        
                                } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
                                        parameters <- c("thetaBS","psiBS","pEti","thetaSS","Icat")
                                        
                                } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
                                        parameters <- c("thetaBS","psiBS","pEti","thetaSS")
                                        
                                }
                                rst.bugs   <- call.bugs(data, inits, parameters,...);
                                rst.bugs
                        }
                        
                        if (Jss.only == 0)
                        {
                                if (mcmc_options$ppd==TRUE){
                                        gs <- mybugs("model_plcm_ppd.bug")
                                } else {
                                        gs <- mybugs("model_plcm.bug")
                                }
                        }
                        else # <----------------------------------------------- Detian: if there are SS only data.
                        {
                                if (mcmc_options$ppd==TRUE)
                                {
                                        stop("PPD with SS only data to be built")
                                }
                                else
                                {
                                        print("SS only data detected:")
                                        gs <- mybugs("model_plcm_ssonly.bug")
                                }
                        }
                }
                
        }else{
                #nplcm: conditional dependence model
                cat("Number of subclasses: ", model_options$k_subclass,"\n")
                # nplcm: BrS only:
                if (model_data_source[1]=="yes" &
                            model_data_source[2]=="no" &
                            model_data_source[3]=="no"){
                        
                        MBS.case <- Mobs$MBS[Y==1,]
                        MBS.ctrl <- Mobs$MBS[Y==0,]
                        MBS      <- as.matrix(rbind(MBS.case,MBS.ctrl))
                        
                        K        <- model_options$k_subclass
                        
                        if (model_options$Eti_prior=="overall_uniform"){
                                alpha    <-  rep(1,Jfull)
                        }
                        TPR_prior_list <- TPR_prior_set(model_options)
                        
                        alphaB <- TPR_prior_list$alphaB
                        betaB <- TPR_prior_list$betaB
                        
                        mybugs <- function(...){
                                inits      <- function(){list(pEti = rep(1/Jallowed,Jallowed),
                                                              r0 = c(rep(.5,K-1),NA),
                                                              r1 = cbind(matrix(rep(.5,Jfull*(K-1)),
                                                                                nrow=Jfull,ncol=K-1),
                                                                         rep(NA,Jfull)),
                                                              alphadp0 = 1)};
                                data       <- c("Nd","Nu","Jfull","Jallowed",
                                                "alpha","template","K",
                                                "MBS","alphaB","betaB");
                                
                                if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
                                        parameters <- c("pEti","Lambda","Eta","alphadp0","MBS.new",
                                                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                                                        "ThetaBS","PsiBS")
                                        
                                } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
                                        parameters <- c("pEti","Lambda","Eta","alphadp0","Icat","MBS.new",
                                                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                                                        "ThetaBS","PsiBS")
                                        
                                } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
                                        parameters <- c("pEti","Lambda","Eta","alphadp0","Icat",
                                                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                                                        "ThetaBS","PsiBS")
                                        
                                } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
                                        parameters <- c("pEti","Lambda","Eta","alphadp0",
                                                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                                                        "ThetaBS","PsiBS")
                                        
                                }
                                rst.bugs   <- call.bugs(data, inits, parameters,...);
                                rst.bugs
                        }
                        
                        if (mcmc_options$ppd==TRUE){
                                gs <- mybugs("model_nplcm_brsonly_ppd.bug")
                        } else {
                                gs <- mybugs("model_nplcm_brsonly.bug")
                        }
                }
                
                
                ## nplcm: BrS + SS:
                if (model_data_source[1]=="yes" &
                            model_data_source[2]=="yes" &
                            model_data_source[3]=="no"){
                        
                        MBS.case <- Mobs$MBS[Y==1,]
                        MBS.ctrl <- Mobs$MBS[Y==0,]
                        MBS      <- as.matrix(rbind(MBS.case,MBS.ctrl))
                        
                        MSS.case <- Mobs$MSS[Y==1,]
                        MSS.case <- as.matrix(MSS.case)
                        
                        SS_index <- which(colMeans(is.na(MSS.case))<0.9)#.9 is arbitrary; any number <1 will work.
                        JSS      <- length(SS_index)
                        MSS      <- MSS.case[,SS_index]
                        
                        K        <- model_options$k_subclass
                        
                        if (model_options$Eti_prior=="overall_uniform"){
                                alpha    <-  rep(1,Jfull)
                        }
                        
                        TPR_prior_list <- TPR_prior_set(model_options)
                        
                        alphaB <- TPR_prior_list$alphaB
                        betaB <- TPR_prior_list$betaB
                        alphaS <- TPR_prior_list$alphaS
                        betaS <- TPR_prior_list$betaS
                        
                        mybugs <- function(...){
                                inits      <- function(){list(pEti = rep(1/Jallowed,Jallowed),
                                                              r0 = c(rep(.5,K-1),NA),
                                                              r1 = cbind(matrix(rep(.5,Jfull*(K-1)),
                                                                                nrow=Jfull,ncol=K-1),
                                                                         rep(NA,Jfull)),
                                                              alphadp0 = 1)};
                                data       <- c("Nd","Nu","Jfull","Jallowed",
                                                "alpha","template","K",
                                                "JSS","MSS",
                                                "MBS","alphaB","betaB","alphaS","betaS");
                                
                                if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
                                        parameters <- c("pEti","Lambda","Eta","alphadp0","MBS.new",
                                                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                                                        "ThetaBS","PsiBS","thetaSS")
                                        
                                } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
                                        parameters <- c("pEti","Lambda","Eta","alphadp0","Icat","MBS.new",
                                                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                                                        "ThetaBS","PsiBS","thetaSS")
                                        
                                } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
                                        parameters <- c("pEti","Lambda","Eta","alphadp0","Icat",
                                                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                                                        "ThetaBS","PsiBS","thetaSS")
                                        
                                } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
                                        parameters <- c("pEti","Lambda","Eta","alphadp0",
                                                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                                                        "ThetaBS","PsiBS","thetaSS")
                                        
                                }
                                rst.bugs   <- call.bugs(data, inits, parameters,...);
                                rst.bugs
                        }
                        
                        if (mcmc_options$ppd==TRUE){
                                gs <- mybugs("model_nplcm_ppd.bug")
                        } else {
                                gs <- mybugs("model_nplcm.bug")
                        }
                }
                
        }
        
        #END fitting model---------------------------------------------------------
        
        # return the WinBUGS fitted results:
        return(gs)
}#END function























