#' set true positive rate (TPR) prior ranges
#'
#' Current prior assignment let bacteria NPPCR to have uniform range, viral NPPCR
#' to have .5-1.0 range. The PCP (a fungus) NPPCR TPR is also set to be .5-1.0; PCP
#' has no blood culture measurements. Also, not all the bacteria have blood culture
#' measurments. One question is whether to use informative NPPCR range or non-informative
#' NPPCR range (0-100%).
#'
#' DN: 1.make the assignment of prior dependent on the BCX data availability
#' or species (bacteria)?
#'
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
#' @return Parameters for the TPR priors, separately for BrS
#'        and SS
#'
#' @export
TPR_prior_set <- function(model_options){

      model_data_source <- rep(NA,3)
      names(model_data_source) <- c("MBS","MSS","MGS")
      model_data_source[1] <- c("no","yes")["BrS"%in%model_options$M_use+1]
      model_data_source[2] <- c("no","yes")["SS"%in%model_options$M_use+1]
      model_data_source[3] <- c("no","yes")["GS"%in%model_options$M_use+1]

      pathogen_list <- model_options$pathogen_list
      Jfull         <- length(pathogen_list)

      # only BrS data is used:
      if (model_data_source[1]=="yes" &
            model_data_source[2]=="no" &
             model_data_source[3]=="no"){
        temp_cat_ind <- sapply(model_options$pathogen_list,
                               function(path) {which(model_options$pathogen_cat$X==path)})
        temp_cat     <- model_options$pathogen_cat[temp_cat_ind,]

        if (model_options$TPR_prior=="noninformative"){
           alphaB = rep(1,Jfull)
           betaB  = rep(1,Jfull)
        }else if (model_options$TPR_prior=="informative"){
          #virus have BrS measurement sensitivity in the range of 50-100%.
          # currently use beta (6,2). bacterial use 0-100%. currently we use beta(1,1).
          alphaB <- rep(NA,Jfull)
          betaB  <- rep(NA,Jfull)

          for (t in 1:nrow(temp_cat)){
            if (temp_cat$pathogen_type[t]=="B"){
              alphaB[t] <- 1
              betaB[t]  <- 1
            }
            if (temp_cat$pathogen_type[t]=="V" | temp_cat$pathogen_type[t]=="F"){
              alphaB[t] <- 6
              betaB[t]  <- 2
            }
          }
        }
        res <- list (alphaB = alphaB, betaB = betaB,
                     used_cat = temp_cat)
      }

      # both BrS and SS data are used:
      if (model_data_source[1]=="yes" &
            model_data_source[2]=="yes" &
              model_data_source[3]=="no"){

        temp_cat_ind <- sapply(model_options$pathogen_list,
                               function(path) {which(model_options$pathogen_cat$X==path)})
        temp_cat     <-model_options$pathogen_cat[temp_cat_ind,]

        if (model_options$TPR_prior[1]=="noninformative"){
             alphaB <- rep(1,Jfull)
             betaB  <- rep(1,Jfull)
          }else if (model_options$TPR_prior[1]=="informative"){
            #virus have BrS measurement sensitivity in the range of 50-100%.
            # currently use beta (6,2). bacterial use 0-100%. currently we use beta(1,1).
            alphaB <- rep(NA,Jfull)
            betaB  <- rep(NA,Jfull)

            for (t in 1:nrow(temp_cat)){
              if (temp_cat$pathogen_type[t]=="B"){
                alphaB[t] <- 1
                betaB[t]  <- 1
              }
              if (temp_cat$pathogen_type[t]=="V" | temp_cat$pathogen_type[t]=="F"){
                alphaB[t] <- 6
                betaB[t]  <- 2
              }
            }
          }

        if (model_options$TPR_prior[2]=="noninformative"){
            alphaS <- rep(1,JSS)
            betaS  <- rep(1,JSS)
          } else if (model_options$TPR_prior[2]=="informative"){
            #For bacteria, informative sensitivities lie in the range of
            # (.05,.15):
            alphaS <- rep(NA,JSS)
            betaS  <- rep(NA,JSS)
            temp_cat_ind <- sapply(model_options$pathogen_list,
                                   function(path) {which(model_options$pathogen_cat$X==path)})
            temp_cat     <- model_options$pathogen_cat[temp_cat_ind,]
            for (t in 1:JSS){
              temp_param <- beta_parms_from_quantiles(c(.05,.15),p=c(0.025,.975),plot=FALSE)
              alphaS[t] <- temp_param$a
              betaS[t] <- temp_param$b
            }
         }

         res <- list(alphaB = alphaB, betaB = betaB,
                    alphaS = alphaS, betaS = betaS,
                    used_cat = temp_cat)
       }

      #return results:
       res

}

