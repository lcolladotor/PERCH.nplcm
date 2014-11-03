#' Exploratory Data Analysis (EDA)
#'
#' Outputs tables and figures to visualize data. Current functionality includes:
#' 1.total no. of positives in cases and controls
#' 2.bubble plot of pairwise associations BrS measurements, separately for controls
#'   and cases.
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
#'
#' @param eda_options A list of options for EDA.
#' @return Tables and figures
#'
#' @export
#'
eda <- function(Mobs, Y, X, eda_options,pathogen_name){
  if (is.null(Mobs$MBS)){
     stop("No bronze-standard data!")
  }else{

  MBS.case <- Mobs$MBS[Y==1,]
  MBS.ctrl <- Mobs$MBS[Y==0,]
  n1       <- nrow(MBS.case)
  n0       <- nrow(MBS.ctrl)
  J        <- ncol(MBS.case)

#1. comparing the total no. of positives for cases and controls
  if (eda_options$total_positives == TRUE){
        ct1 = rowSums(MBS.case)
        ct0 = rowSums(MBS.ctrl)
        tb1 = rep(NA,J+1);names(tb1)=c(0,1:J)
        tb0 = tb1
        for (j in 1:(J+1)){
          tb1[j] = sum(ct1==j-1);tb0[j]=sum(ct0==j-1)
        }

        barplot(rbind(tb1/n1,tb0/n0),beside=TRUE,ylim=c(0,max(c(tb1/n1,tb0/n0))+.1),
                legend.text = c("case", "control"),
                xlab="no. of positives",
                ylab="sample frequency",
                main=paste(eda_options$X_names,eda_options$X_values,
                           sep="="))
  }
#2. bubble plot of pairs wise associations:
# a) unadjusted for covaraites:
  if (eda_options$bubble_plot == TRUE){
        logORmat_bubble(MBS.case,MBS.ctrl,pathogen_name)
  }


  }
}






