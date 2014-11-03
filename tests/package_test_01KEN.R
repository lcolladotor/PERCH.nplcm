rm(list=ls())
library(lubridate)
library(sets)
library(R2WinBUGS)
library(gplots)
library(RColorBrewer)
library(binom)
library(nplcm)

## data cleaning:-------------------------------------------------

# combine two sub-sites:
# 06NTH and 07STH  --> THA,
# 08MBA and 09DBA  --> BAN:
RawMeasDir <- "D:/Research_PERCH/working_version/package_test/PQ_20JUN14.csv"
PERCH_data_with_newSITE <- combine_subsites(RawMeasDir,
                             subsites_list = list(c("06NTH","07STH"),
                                                  c("08MBA","09DBA")),
                             newsites_vec  = c("THA","BAN"))

write.csv(PERCH_data_with_newSITE,"D:/Research_PERCH/working_version/package_test/PERCH_data_with_newSITE.csv")

# set the measurement directory for the data with combined site names:
MeasDir    <- "D:/Research_PERCH/working_version/package_test//PERCH_data_with_newSITE.csv"

PathCatDir <- "D:/Research_PERCH/working_version/package_test/pathogen_category.csv"

# Pathogen   <- c("HINF","MCAT","PNEU","SASP","SAUR",#with BCX measures.
#                          "BORD","C_PNEU","M_PNEU",
#                          "PCP",# up to now, all are bacteria or fungus.
#                          "ADENOVIRUS","CMV","COR_229","COR_43","COR_63","COR_HKU",
#                          "FLU_C","HBOV","HMPV_A_B",
#                          "INFLUENZA_A","INFLUENZA_B",
#                          "PARA1","PARA2","PARA3","PARA4",
#                          "PV_EV","RHINO","RSV_A_B")

pathogen10 <- c("HINF", "SASP", "SAUR",# all bacteria up to now
               "ADENOVIRUS","COR_43", "FLU_C",
               "HMPV_A_B", "PARA1","RHINO","RSV_A_B")

silver.only = c("ECOL","PAER")  ############### silver only added ################### edited by Detian

Pathogen <- c(silver.only, pathogen10) ############## silver only need to be placed in front of others 

# list the pathogen categories:
Pathogen_cat <- read.csv(PathCatDir)

Specimen  <- c("NP","B")
Test      <- c("PCR","CX")

X_strat        <- c("newSITE","CASECONT")

newALLSITES <- c("01KEN","02GAM","03MAL","04ZAM","05SAF","THA","BAN")
curr.site   <- 1
sitename    <- newALLSITES[curr.site]
Xval.not.disease.status   <- sitename #stratifying variable

extra_X = c("ENRLDATE","patid","AGECAT","HIV")

datacase <- extract_data_raw(Pathogen,Specimen,Test,
                             X_strat,c(Xval.not.disease.status,1),
                             extra_covariates = extra_X,
                             MeasDir,PathCatDir,silent=TRUE)

datactrl <- extract_data_raw(setdiff(Pathogen,silver.only),Specimen,Test, ############ edited by Detian
                             X_strat,c(Xval.not.disease.status,2),
                             extra_covariates = extra_X,
                             MeasDir,PathCatDir,silent=TRUE)

# deal with incompleteness of measurement data:

# get pathogens that have many BcX measurements:
SS_index  <- which(colMeans(is.na(datacase$BCX))<.9)
JSS       <- length(SS_index)

# get case/control indices who have complete observations on NPPCR&BCX(cases),
# or NPPCR(ctrls):
complete_case_index <- which(rowMeans(is.na(datacase$NPPCR[,-(1:length(silver.only))]))==0 & #### adjust on silver only # edited by Detian
                                        rowMeans(is.na(datacase$BCX[,SS_index]))==0)
complete_ctrl_index <- which(rowMeans(is.na(datactrl$NPPCR))==0)

# actual complete numbers of cases/controls:
Nd  <- length(complete_case_index)
Nu  <- length(complete_ctrl_index)


ctrl.SSonly.NP = as.data.frame(matrix(NA,nrow=Nu,ncol=length(silver.only))) ##### fill in space with NA ###### Detian
colnames(ctrl.SSonly.NP) = paste(silver.only,"NPPCR",sep="_")               ##################################

M_NPPCR    <- rbind(datacase$NPPCR[complete_case_index,],
                    cbind(ctrl.SSonly.NP,datactrl$NPPCR[complete_ctrl_index,])) ##############################
Y          <- c(rep(1,Nd),rep(0,Nu))
M_BCX_ctrl <- as.data.frame(matrix(NA,nrow=Nu,ncol=length(Pathogen)))
colnames(M_BCX_ctrl) <- paste(Pathogen,"BCX",sep="_")
M_BCX     <- rbind(datacase$BCX[complete_case_index,],M_BCX_ctrl)

Rdate.case <- as.Date(datacase$ENRLDATE[complete_case_index], "%d%B%Y")
Rdate.ctrl <- as.Date(datactrl$ENRLDATE[complete_ctrl_index], "%d%B%Y")

uniq.month.case <- unique(paste(month(Rdate.case),year(Rdate.case),sep="-"))
uniq.month.ctrl <- unique(paste(month(Rdate.ctrl),year(Rdate.ctrl),sep="-"))

symm.diff.dates <- as.set(uniq.month.case)%D% as.set(uniq.month.ctrl)
if (length(symm.diff.dates)!=0){
  cat("Cases and controls have different enrollment months:","\n")
  print(symm.diff.dates)
}

# observations unordered with respect to date:
datobs <- cbind(M_NPPCR,
                M_BCX,
                Y = c(rep(1,Nd),rep(0,Nu)),
                newSITE  = c(datacase$newSITE[complete_case_index],
                             datactrl$newSITE[complete_ctrl_index]),
                raw.date = c(Rdate.case,Rdate.ctrl),
                patid    = c(as.character(datacase$patid[complete_case_index]),
                             as.character(datactrl$patid[complete_ctrl_index])),
                AGECAT   = c(datacase$AGECAT[complete_case_index],
                             datactrl$AGECAT[complete_ctrl_index]),
                HIV      = c(datacase$HIV[complete_case_index],
                             datactrl$HIV[complete_ctrl_index])
                )

datobs <- as.data.frame(datobs)

# separately order cases and controls according to enrollment dates, and
# place cases' data on the top:
datobs <- datobs[c(order(Rdate.case),Nd+order(Rdate.ctrl)),]
## end of PERCH data cleaning -----------------------------------------

## BEGIN preparation of data for nplcm_fit function-------------------
# create the list of measurements
# look for pathogens that have more than one positives in BCX:
           cat("Total positives in MSS:","\n")
           print(table(rowSums(M_BCX[1:Nd,SS_index])))
           BCX_more_than_one_index <- which(rowSums(M_BCX[1:Nd,SS_index])>3)
           if (length(BCX_more_than_one_index)>0){
               cat("Removed case(s) who have more than one positive in BCX:","\n")
               print(datobs$patid[BCX_more_than_one_index])
               print(M_BCX[BCX_more_than_one_index,SS_index])

               Mobs   <- list(MBS = M_NPPCR[-BCX_more_than_one_index,],
                              MSS = M_BCX[-BCX_more_than_one_index,],
                              MGS = NA)
               Y      <- datobs$Y[-BCX_more_than_one_index]
               X      <- datobs[-BCX_more_than_one_index,c("raw.date","AGECAT","HIV")]
               Nd     <- Nd - length(BCX_more_than_one_index)
           } else {
               Mobs   <- list(MBS = M_NPPCR,
                              MSS = M_BCX,
                              MGS = NA)
               Y      <- datobs$Y
               X      <- datobs[,c("raw.date","AGECAT","HIV")]
           }

## BEGIN exploratory data analysis----------------------------------------------
eda_options <- list(X_names = list("site"),X_values = list(sitename),
                    marginal_rate_comparison = TRUE,
                    bubble_plot  = TRUE)
eda(Mobs, Y, X, eda_options,Pathogen) ########### error reported: Error in if (eda_options$total_positives == TRUE) { : argument is of length zero
## END exploratory data analysis----------------------------------------------


# to add: 1.stratifying functionality, not regression; 2. borrow TPR information
#         across the same category of pathogens
      model_options1 <- list(M_use = c("BrS"), # has to be a subset of non-NA entries in Mobs.
                            k_subclass = 1,
                            TPR_prior  = c("noninformative"),#same length as M_use.
                            Eti_prior  = "overall_uniform",
                            allowed_list  = Pathogen,
                            pathogen_list = Pathogen,
                            FPR_regression = NULL,
                            Eti_regression = NULL)

      model_options2 <- list(M_use = c("BrS","SS"), # has to be a subset of non-NA entries in Mobs.
                             k_subclass = 1,
                             TPR_prior  = c("noninformative","noninformative"),#same length as M_use.
                             Eti_prior  = "overall_uniform",
                             allowed_list  = Pathogen,
                             pathogen_list = Pathogen,
                             FPR_regression = NULL,
                             Eti_regression = NULL)

      model_options3 <- list(M_use = c("BrS"), # has to be a subset of non-NA entries in Mobs.
                             k_subclass = 2,
                             TPR_prior  = c("noninformative"),#same length as M_use.
                             Eti_prior  = "overall_uniform",
                             allowed_list  = Pathogen,
                             pathogen_list = Pathogen,
                             FPR_regression = NULL,
                             Eti_regression = NULL)

      model_options4 <- list(M_use = c("BrS","SS"), # has to be a subset of non-NA entries in Mobs.
                             k_subclass = 2,
                             TPR_prior  = c("noninformative","noninformative"),#same length as M_use.
                             Eti_prior  = "overall_uniform",
                             allowed_list  = Pathogen,
                             pathogen_list = Pathogen,
                             FPR_regression = NULL,
                             Eti_regression = NULL)

      model_options5 <- list(M_use = c("BrS","SS"), # has to be a subset of non-NA entries in Mobs.
                             k_subclass = 10,
                             TPR_prior  = c("informative","informative"),#same length as M_use.
                             Eti_prior  = "overall_uniform",
                             allowed_list  = Pathogen,
                             pathogen_list = Pathogen,
                             FPR_regression = NULL,
                             Eti_regression = NULL,
                             pathogen_cat = Pathogen_cat)

      model_options <- model_options2

      Date = gsub("-", "", Sys.Date())

      fname   <- paste0("D:\\Research_PERCH\\working_version\\",Date,"_",sitename)
      dir.create(fname)
      fullname <- fname

      ## for finer scenarios
      result.folder <- fullname
      dir.create(result.folder)

      bugs.model.dir <- "D:\\Research_PERCH\\working_version\\winbugs_model_package\\"
      winbugs.dir    <- "D:\\WinBUGS14\\"

      mcmc_options <- list(debugstatus = !TRUE,
                                n.chains   = 1,
                                n.itermcmc = 100, # 50000
                                n.burnin   = 20, # 10000
                                n.thin     = 2,
                           individual.pred = !TRUE,
                                    ppd    = !TRUE,
                             result.folder = result.folder,
                             bugsmodel.dir = bugs.model.dir,
                             winbugs.dir   = winbugs.dir)
## END preparation of data for nplcm_fit function---------------------

## WinBUGS fitting
gs <- nplcm_fit2(Mobs,Y,X,model_options,mcmc_options)

# DN: 1.incorporate results visualization
#     2.check simulation compatibility
pdf(paste0(result.folder,"\\",sitename,"_three_panel_plot.pdf"),width=11,height=10)
nplcm_plot_three_panel(DIR_NPLCM = result.folder,Mobs,Y,X,model_options)
dev.off()

pdf(paste0(result.folder,"\\",sitename,"_individual_diagnosis.pdf"),
    width=16,height=16)
par(mfrow=c(4,4))
nplcm_plot_individual_diagnosis(DIR_NPLCM = result.folder,Mobs,Y,X,model_options)
dev.off()



save.image(paste0(result.folder,"\\for_replot.RDATA"))
