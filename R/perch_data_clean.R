#' fucntion for cleaning PERCH data
#'
#' @param clean_options The list of options for cleaning PERCH data. The specific
#' elements are as follows.
#' \code{case_def}: the variable name that is used for case definition.
#' \code{case_def_val}: The value of the case-definition variable.
#' \code{X_strat}: A vector of strings, each defining the variables used to
#' stratify the data.
#' \code{X_strat_val}: A list of actual values that X_strat should take to
#' get stratified data sets.
#' \code{Pathogen}: The vector of pathogen names. It has to be a subset of
#' pathogen category information in the \code{PathCatDir}.
#' \code{extra_X}: A vector of strings, each being the covariate name that one
#' wants to have in the data set to be analyzed.
#' \code{RawMeasDir}: The file path to the raw data set
#' \code{newSite_write_Dir}: The file path where a new/cleaned data set for actual
#'  will be written
#' \code{MeasDir}: The file path to the cleaned data set
#' \code{PathCatDir}: The file path to the pathogen category list (.csv)
#'
#' @return         list(Mobs,Y,X,JSS,Pathogen_MSS_ordered,Pathogen_cat)
#'
#' @export
perch_data_clean <- function(clean_options){#         = "C:/package_test/pathogen_category.csv"){


        case_def     <- clean_options$case_def
        case_def_val <- clean_options$case_def_val
        X_strat      <- clean_options$X_strat
        X_strat_val  <- clean_options$X_strat_val
        Pathogen     <- clean_options$Pathogen
        extra_X      <- clean_options$extra_X
        RawMeasDir   <- clean_options$RawMeasDir
        newSite_write_Dir <- clean_options$newSite_write_Dir
        MeasDir           <- clean_options$MeasDir
        PathCatDir        <- clean_options$PathCatDir

        # combine two sub-sites:
        # 06NTH and 07STH  --> THA,
        # 08MBA and 09DBA  --> BAN:

        PERCH_data_with_newSITE <- combine_subsites(RawMeasDir,
                                    subsites_list = list(c("06NTH","07STH"),
                                                         c("08MBA","09DBA")),
                                    newsites_vec  = c("THA","BAN"))

        write.csv(PERCH_data_with_newSITE,newSite_write_Dir)

        # list the pathogen categories:
        Pathogen_cat <- read.csv(PathCatDir)

        Specimen  <- c("NP","B")
        Test      <- c("PCR","CX")

        #X_strat        <- c("newSITE","CASECONT")

        #newALLSITES <- c("01KEN","02GAM","03MAL","04ZAM","05SAF","THA","BAN")

        #sitename    <- newALLSITES[curr.site]
        #Xval.not.disease.status   <- sitename #stratifying variable

        datacase <- extract_data_raw(Pathogen,Specimen,Test,
                                     c(X_strat,case_def),append(X_strat_val,case_def_val),
                                     extra_covariates = extra_X,
                                     MeasDir,PathCatDir,silent=TRUE)
        #write.csv(datacase,"C:/package_test/datacase.csv")

        datactrl <- extract_data_raw(Pathogen,Specimen,Test,
                                     c(X_strat,"CASECONT"),append(X_strat_val,2),
                                     extra_covariates = extra_X,
                                     MeasDir,PathCatDir,silent=TRUE)
        #write.csv(datactrl,"C:/package_test/datactrl.csv")

        # deal with incompleteness of measurement data:

        # get pathogens that have many BcX measurements:
        SS_index  <- which(colMeans(is.na(datacase$BCX))<.9)
        cat("Pathogens that have blood culture measurements:","\n",
            Pathogen[SS_index],"\n")
        JSS       <- length(SS_index)
        Jfull     <- length(Pathogen)
        MSS_avail_index <- c(SS_index,(1:Jfull)[-SS_index])

        # get case/control indices who have complete observations on NPPCR&BCX(cases),
        # or NPPCR(ctrls):
        complete_case_index <- which(rowMeans(is.na(datacase$NPPCR))==0 &
                                       rowMeans(is.na(datacase$BCX[,SS_index]))==0)
        complete_ctrl_index <- which(rowMeans(is.na(datactrl$NPPCR))==0)

        # actual complete numbers of cases/controls:
        Nd  <- length(complete_case_index)
        Nu  <- length(complete_ctrl_index)

        M_NPPCR    <- rbind(datacase$NPPCR[complete_case_index,],
                            datactrl$NPPCR[complete_ctrl_index,])
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

        # create the list of measurements
        # look for pathogens that have more than one positives in BCX:
        cat("Total positives in MSS:","\n")
        print(table(rowSums(M_BCX[1:Nd,SS_index])))
        BCX_more_than_one_index <- which(rowSums(M_BCX[1:Nd,SS_index])>1)
        if (length(BCX_more_than_one_index)>0){
          cat("Removed case(s) who have more than one positive in BCX:","\n")
          print(datobs$patid[BCX_more_than_one_index])
          print(M_BCX[BCX_more_than_one_index,SS_index])

          Mobs   <- list(MBS = M_NPPCR[-BCX_more_than_one_index,MSS_avail_index],
                         MSS = M_BCX[-BCX_more_than_one_index,MSS_avail_index],
                         MGS = NA)
          Y      <- datobs$Y[-BCX_more_than_one_index]
          X      <- datobs[-BCX_more_than_one_index,c("raw.date","AGECAT","HIV")]
          Nd     <- Nd - length(BCX_more_than_one_index)
        } else {
          Mobs   <- list(MBS = M_NPPCR[,MSS_avail_index],
                         MSS = M_BCX[,MSS_avail_index],
                         MGS = NA)
          Y      <- datobs$Y
          X      <- datobs[,c("raw.date","AGECAT","HIV")]
        }

        #order the pathogens so that those have SS data comes first:
        Pathogen_MSS_ordered = Pathogen[MSS_avail_index]
        #get the list of pathogen categories for each of the above ordered pathogens:
        path_cat_ind <- sapply(1:Jfull,function(i)
                                   which(Pathogen_cat$X==Pathogen_MSS_ordered[i]))
        list(Mobs=Mobs,Y=Y,X=X,
             JSS=JSS,
             Pathogen_MSS_ordered = Pathogen_MSS_ordered,
             Pathogen_cat = Pathogen_cat[path_cat_ind,])
}
