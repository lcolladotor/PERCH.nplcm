#' Fitting nested partially-latent class models (lower level)
#'
#'
#' Current implementation use WinBUGS in Windows system
#' (see readme file for an instruction to install WinBUGS on windows 7 or 8.)
#' Developer Note (DN):
#' 1.gold-standard measurements not implemented
#' 2.need to include Wei's introduction to installing WinBUGS for PQ analysis
#' 3.need code for borrowing information of TPRs across (pathogens,x_borrow_TPR)
#' 4.need to build in regression functionality
#' 5.allowing those pathogens that have only silver-standard measurements.
#'
#'
#'
#' @param Mobs A list of measurements as model input. The elements of the list should
#' include (according to \code{model}) \code{MBS}, \code{MSS}, and \code{MGS}.
#' Here, \code{MBS} is a data frame of bronze-standard (BrS) measurements.
#' Rows are subjects, columns are pathogens targeted in the BrS measurements
#' (e.g. nasalpharyngeal PCR). These mesaurements have imperfect sensitivity/specificty.
#' \code{MSS} is a data frame of silver-standard (SS) measurements. Rows are subjects,
#'   columns are pathogens targeted in the SS measurements (e.g. blood culture).
#'   These measurements have perfect specificity but imperfect sensitivity.
#' \code{MGS} is a data frame of gold-standard (GS) measurements. Rows are subject,
#' columns are pathogens targeted in the GS measurements. These measurements
#' have both perfect sensitivity and perfect specificity.
#'
#' @param Y A vector of binary variables specifying the disease status (1 for case;
#' 0 for control).
#' @param X A design matrix. For FPR and etiology regressions.
#' @param model_options A list of model options.
#'
#' \code{M_use} List of measurements to be used in the model;
#' \code{k_subclass}The number of nested subclasses. It equals 1 if the
#' model is plcm; its value will be useful if the model is nplcm.
#' \code{TPR_prior} the prior for the measurements used. Its length should be
#' the same with \code{M_use};
#' \code{Eti_prior} the vector of parameters in etiology prior;
#' \code{allowed_list} for the allowed list of pathogen/pathogen combinations infecting
#' the lung;
#' \code{pathogen_list} The complete list of pathogens considered in the analysis.
#' \code{pathogen_cat} The two-column dataframe that has category of pathogens: virus (V), bacteria (B)
#' and fungi (F)
#'
#'
#' @param mcmc_options A list of Markov chain Monte Carlo related options.
#'
#' \code{debugstatus} for whether to stop WinBUGS window after it finishes Gibbs sampling;
#' \code{n.chains} for the number of MCMC chains;
#' \code{n.burnin} for the number of burn-in samples;
#' \code{n.thin} keep every other n.thin samples after burn-in period;
#' \code{individual.pred} whether to perform individual prediction;
#' \code{ppd} whether to perform posterior predictive checking
#' \code{result.folder} for the path to folder storing the results;
#' \code{bugsmodel.dir} for the directory to WinBUGS model files;
#' \code{winbugs.dir} for the directory where software WinBUGS 1.4 is installed.

#' @return A WinBUGS result, fitted by function \code{bugs()} from the R2WinBUGS package.
#'
#' @export
#'
nplcm_fit<-function(Mobs,Y,X,model_options,mcmc_options){#BEGIN function
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

              mybugs <- function(...){
                    inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                                  psiBS   = rbeta(Jfull,1,1))};
                    data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template",
                                    "MBS","JSS","MSS",
                                    "alphaB","betaB","alphaS","betaS");

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

              if (mcmc_options$ppd==TRUE){
                gs <- mybugs("model_plcm_ppd.bug")
              } else {
                gs <- mybugs("model_plcm.bug")
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























