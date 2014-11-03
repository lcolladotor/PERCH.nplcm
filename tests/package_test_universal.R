rm(list=ls())
library(lubridate)
library(sets)
library(R2WinBUGS)
library(gplots)
library(RColorBrewer)
library(binom)
library(coda)
library(nplcm)

newALLSITES           <- c("01KEN","02GAM","03MAL","04ZAM","05SAF","THA","BAN")
current_study_site    <- 6
sitename              <- newALLSITES[current_study_site]
RawMeasDir            <- "D:/Research_PERCH/working_version/package_test/PQ_02OCT14.csv" ################ edited by Detian

pathogen10new <- c("HINF", "SASP", "SAUR",# all bacteria up to now
                   "ECOL","PAER", # silver only pathogen  ########################## edited by Detian
                "ADENO","COR_43", "FLU_C",
                "HMPV_A_B", "PARA_1","RHINO","RSV")

pathogen27new   <- c("HINF","MCAT","PNEU_N","SASP","SAUR",#with BCX measures.
                    "BOPE","C_PNEU","M_PNEU",
                    "PCP",# up to now, all are bacteria or fungus.
                    "ECOL","PAER", # silver only pathogen  ########################## edited by Detian
                    "ADENO","CMV","COR_229","COR_43","COR_63","COR_HKU",
                    "FLU_C","HBOV","HMPV_A_B",
                    "FLU_A","FLU_B",
                    "PARA_1","PARA_2","PARA_3","PARA_4",
                    "PV_EV","RHINO","RSV")

Pathogen_anyorder   <- pathogen10new

## output cleaned PERCH data:
perch_data_clean_options <- list (case_def           =  "CXRFINCAT_5",
                                  case_def_val       =  1,
                                  X_strat            = c("newSITE"),
                                  X_strat_val        =  list(sitename),
                                  Pathogen           =  Pathogen_anyorder,
                                  extra_X            = c("ENRLDATE","patid","AGECAT","HIV"),#covariates besides case/control, and site
                                  RawMeasDir         = RawMeasDir,
                                  newSite_write_Dir  = "D:/Research_PERCH/working_version/package_test/PERCH_data_with_newSITE.csv",
                                  MeasDir            = "D:/Research_PERCH/working_version/package_test/PERCH_data_with_newSITE.csv",
                                  PathCatDir         = "D:/Research_PERCH/working_version/package_test/pathogen_category.csv")
cleaned_data  <- perch_data_clean(perch_data_clean_options)

## output specific data sets to be used in Bayesian fitting:
        Mobs  <- cleaned_data$Mobs
           Y  <- cleaned_data$Y
           X  <- cleaned_data$X
 Pathogen_cat <- cleaned_data$Pathogen_cat
          JSS <- cleaned_data$JSS
 Pathogen     <- cleaned_data$Pathogen_MSS_ordered

## BEGIN exploratory data analysis----------------------------------------------
eda_options <- list(             X_names     = list("site"),
                                 X_values    = list(sitename),# for plotting title on top of the figure.
                             total_positives = TRUE,
                                bubble_plot  = TRUE)
## END exploratory data analysis----------------------------------------------



## BEGIN preparation of data for nplcm_fit function-------------------
# to add: 1.stratifying functionality, not regression; 2. borrow TPR information
#         across the same category of pathogens
#       model_options1 <- list(M_use = c("BrS"), # has to be a subset of non-NA entries in Mobs.
#                             k_subclass = 1,
#                             TPR_prior  = c("noninformative"),#same length as M_use.
#                             Eti_prior  = "overall_uniform",
#                             allowed_list  = Pathogen,
#                             pathogen_list = Pathogen,
#                             FPR_regression = NULL,
#                             Eti_regression = NULL,
#                             pathogen_cat = Pathogen_cat)
#
#       model_options2 <- list(M_use = c("BrS","SS"), # has to be a subset of non-NA entries in Mobs.
#                              k_subclass = 1,
#                              TPR_prior  = c("noninformative","noninformative"),#same length as M_use.
#                              Eti_prior  = "overall_uniform",
#                              allowed_list  = Pathogen,
#                              pathogen_list = Pathogen,
#                              FPR_regression = NULL,
#                              Eti_regression = NULL,
#                              pathogen_cat = Pathogen_cat)
#
#       model_options3 <- list(M_use = c("BrS"), # has to be a subset of non-NA entries in Mobs.
#                              k_subclass = 2,
#                              TPR_prior  = c("noninformative"),#same length as M_use.
#                              Eti_prior  = "overall_uniform",
#                              allowed_list  = Pathogen,
#                              pathogen_list = Pathogen,
#                              FPR_regression = NULL,
#                              Eti_regression = NULL,
#                              pathogen_cat = Pathogen_cat)
#
#       model_options4 <- list(M_use = c("BrS","SS"), # has to be a subset of non-NA entries in Mobs.
#                              k_subclass = 2,
#                              TPR_prior  = c("noninformative","noninformative"),#same length as M_use.
#                              Eti_prior  = "overall_uniform",
#                              allowed_list  = Pathogen,
#                              pathogen_list = Pathogen,
#                              FPR_regression = NULL,
#                              Eti_regression = NULL,
#                              pathogen_cat = Pathogen_cat)

      model_options5 <- list(M_use = c("BrS","SS"), # has to be a subset of non-NA entries in Mobs.
                             k_subclass = 5,
                             TPR_prior  = c("informative","informative"),#same length as M_use.
                             Eti_prior  = "overall_uniform",
                             allowed_list  = Pathogen,
                             pathogen_list = Pathogen,
                             FPR_regression = NULL,
                             Eti_regression = NULL,
                             pathogen_cat = Pathogen_cat)

      model_options6 <- list(M_use = c("BrS","SS"), # has to be a subset of non-NA entries in Mobs.
                              k_subclass = 1,
                              TPR_prior  = c("informative","informative"),#same length as M_use.
                              Eti_prior  = "overall_uniform",
                              allowed_list  = Pathogen,
                              pathogen_list = Pathogen,
                              FPR_regression = NULL,
                              Eti_regression = NULL,
                              pathogen_cat = Pathogen_cat)

      Date = gsub("-", "_", Sys.Date())

      fname    <- paste0("D:\\Research_PERCH\\",Date,"_",sitename)
      dir.create(fname)
      fullname <- fname

      ## for finer scenarios
      result.folder <- fullname
      dir.create(result.folder)

      bugs.model.dir <- "D:\\Research_PERCH\\working_version\\winbugs_model_package\\"
      winbugs.dir    <- "D:\\WinBUGS14\\"

      mcmc_options <- list(debugstatus = !TRUE,
                                n.chains   = 1,
                                n.itermcmc = 20,#30000,
                                n.burnin   = 0,#10000,
                                n.thin     = 1,#20,
                           individual.pred = TRUE,
                                    ppd    = !TRUE,
                             result.folder = result.folder,
                             bugsmodel.dir = bugs.model.dir,
                             winbugs.dir   = winbugs.dir)
## END preparation of data for nplcm_fit function---------------------

model_options <- model_options5

## BEGIN recording settings of current analysis-----------------------------------
cat("Results stored in: ", result.folder,"\n")
#data clean options:
dput(perch_data_clean_options,paste0(mcmc_options$result.folder,
                                     "\\data_clean_options.txt"))
#model_options:
dput(model_options,paste0(mcmc_options$result.folder,"\\model_options.txt"))
#mcmc_options:
dput(mcmc_options,paste0(mcmc_options$result.folder,"\\mcmc_options.txt"))

#correlation structure:
pdf(paste0(mcmc_options$result.folder,"/logOR_bubble_plot.pdf"),height=20, width=20)
eda(Mobs, Y, X, eda_options,Pathogen)
dev.off()
## END recording settings of current analysis--------------------------------------

## WinBUGS fitting
gs <- nplcm_fit(Mobs,Y,X,model_options,mcmc_options)

# DN: 1.check simulation study compatibility
#     2.posterior-predictive checking procedure
#      separate: marginal, pairwise associations
#      combined:

## BEGIN results visualization ------------------------------------------------
## plot 1: three-panel plot:
pdf(paste0(result.folder,"\\",sitename,"_three_panel_plot.pdf"),width=11,height=10)
nplcm_plot_three_panel(DIR_NPLCM = result.folder,ss_upperlimit = .3,eti_upperlimit = .5)
dev.off()

## plot 2: individual diagnosis plot:
pdf(paste0(result.folder,"\\",sitename,"_individual_diagnosis.pdf"),width=16,height=16)
par(mfrow=c(4,4))
nplcm_plot_individual_diagnosis(DIR_NPLCM=result.folder,npat=16)
dev.off()
## END of results visualization -----------------------------------------------

## save workspace for future replot:
save.image(paste0(result.folder,"\\for_replot.RDATA"))




