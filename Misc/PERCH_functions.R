###############################################################################
## Functions for conducting both PERCH data analysis
## and simulation studies
##
## Zhenke Wu Nov 26th, 2013
## zhenkewu@gmail.com
###############################################################################
##------------------------------------------------------------------------##
## function to call bugs()
##------------------------------------------------------------------------##

call.bugs <- function(data, inits, parameters, m.file,
                      nitermcmc=n.itermcmc, nburnin=n.burnin, 
                      nthin=n.thin, nchains=n.chains,
                      dic=FALSE, is.debug=debugstatus, workd=result.folder,
                      bugsmodel.dir= "C:\\Users\\Detian Deng\\Desktop\\PERCH\\winbugs_model\\",
                      winbugs.dir  = "D:\\WinBUGS14\\" ,...) {
  require(R2WinBUGS)
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

bugs_CI_nppcr_noGS <- function(...) {
  inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                psiBS = rbeta(Jfull,1,1))};
  data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template","MBS");
  parameters <- c("thetaBS","psiBS","pEti","Icat");
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_CI_nppcr_noGS_from_R <- function(...) {
  inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                psiBS = rbeta(Jfull,1,1))};
  data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template","MBS",
                  "alphaB","betaB");
  parameters <- c("thetaBS","psiBS","pEti","Icat");
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_CI_nppcr_realGS <- function(...) {
  inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                psiBS = rbeta(Jfull,1,1))};
  data       <- c("Nd","Nu","NdGS","Jfull","Jallowed","alpha","template",
                  "MBS","known_Icat");
  parameters <- c("thetaBS","psiBS","pEti","Icat");
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_CI_realGS_no_nppcr <- function(...) {
  inits      <- function(){list(pEti = rep(1/Jfull,Jfull))};
  data       <- c("Nd","NdGS","Jallowed","alpha","known_Icat");
  parameters <- c("pEti","Icat");
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}


bugs_CI_norealGS_no_nppcr <- function(...) {
  inits      <- function(){list(pEti = rep(1/Jfull,Jfull))};
  data       <- c("Nd","Jallowed","alpha");
  parameters <- c("pEti","Icat");
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_CI_nppcr_bcx <- function(...) {
  inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                psiBS = rbeta(Jfull,1,1))};
  data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template",
                  "MBS","JGS","MGS");
  parameters <- c("thetaBS","psiBS","pEti","thetaGS","Icat");
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

## specify sensitivity/specificities from R directly
# bugs_CI_nppcr_bcx_from_R <- function(...) {
#   inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
#                                 psiBS = rbeta(Jfull,1,1))};
#   data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template",
#                   "MBS","JGS","MGS",
#                   "alphaB","betaB","alphaG","betaG");
#   #parameters <- c("thetaBS","psiBS","pEti","thetaGS","Icat");
#   parameters <- c("thetaBS","psiBS","pEti","thetaGS","prod")
#   rst.bugs   <- call.bugs(data, inits, parameters,...);
#   rst.bugs;
# }

bugs_CI_nppcr_bcx_from_R <- function(...) {
  inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                psiBS = rbeta(Jfull,1,1))};
  data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template",
                  "MBS","JGS","MGS",
                  "alphaB","betaB","alphaG","betaG");
  #parameters <- c("thetaBS","psiBS","pEti","thetaGS","Icat");
  parameters <- c("thetaBS","psiBS","pEti","thetaGS","Icat","MBS.new");
  #parameters <- c("thetaBS","psiBS","pEti","thetaGS","prod")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}




bugs_CI_nppcr_bcx_from_R_ABX <- function(...) {
  inits      <- function(){list(thetaBS = matrix(rbeta(Jfull*2,1,1),ncol=2),
                                psiBS = rbeta(Jfull,1,1))};
  data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template",
                  "MBS","JGS","MGS","PAB",
                  "alphaB","betaB","alphaG","betaG");
  #parameters <- c("thetaBS","psiBS","pEti","thetaGS","Icat");
  parameters <- c("thetaBS","psiBS","pEti","thetaGS","Icat","MBS.new");
  #parameters <- c("thetaBS","psiBS","pEti","thetaGS","prod")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

## specify sensitivity/specificities from R directly
bugs_CI_nppcr_bcx_from_R_TPRsame <- function(...) {
  inits      <- function(){list(thetaBS = rbeta(1,1,1),
                                psiBS = rbeta(Jfull,1,1))};
  data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template",
                  "MBS","JGS","MGS",
                  "alphaB","betaB","alphaG","betaG");
  #parameters <- c("thetaBS","psiBS","pEti","thetaGS","Icat");
  #parameters <- c("thetaBS","psiBS","pEti","thetaGS","Icat","MBS.new");
  #parameters <- c("thetaBS","psiBS","pEti","thetaGS","prod")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}


## specify sensitivity/specificities from R directly
bugs_CI_nppcr_bcx_from_R_TPRfix <- function(...) {
  inits      <- function(){list(psiBS = rbeta(Jfull,1,1))};
  data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template",
                  "MBS","JGS","MGS",
                  "alphaG","betaG","thetaBS");
  parameters <- c("psiBS","pEti","thetaGS","Icat");
  #parameters <- c("thetaBS","psiBS","pEti","thetaGS","prod")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}



## specify sensitivity/specificities from R directly
bugs_CI_nppcr_IS_from_R <- function(...) {
  inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                psiBS = rbeta(Jfull,1,1))};
  data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template",
                  "MBS","MGS",
                  "alphaB","betaB",
                  "alphaG.theta","betaG.theta",
                  "alphaG.psi","betaG.psi");
  parameters <- c("thetaBS","psiBS","pEti","thetaGS","psiGS","Icat");
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_pspline_fit <- function(...) {
  inits<-function(){list(beta=c(0,0),b=inits.b,taub=0.01)}
  parameters<-list("taub","sigmab","beta","b")
  data<-list("response","X","Z","n","num.knots")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_seasonality_tprs <- function(...) {
  inits<-function(){list(beta=matrix(0,nrow=Jfull,ncol=2),
                         b=matrix(inits.b,nrow=Jfull,ncol=num.knots,byrow=TRUE),
                         taub=rep(0.01,Jfull),
                         betaEti=rbind(matrix(0,nrow=Jfull-1,ncol=2),c(NA,NA)),
                         bEti=rbind(matrix(inits.bEti,
                                           nrow=Jfull-1,ncol=num.knots.case,
                                           byrow=TRUE),
                                    rep(NA,num.knots.case)),
                         taubEti=c(rep(0.01,Jfull-1)))}
  parameters<-list("taub","sigmab","beta","b",
                   "taubEti","sigmabEti","betaEti","bEti",
                   #"Delta","sigma_Delta",
                   "taubcase","sigmabcase","betacase","bcase",
                   "pEti","Icat",
                   "tpr","mu","fpr.case","fpr.ctrl")
  data<-list("MBS","Nd","Nu","Jfull","template",
              "X","Z","num.knots",#"invsqrt",
              "num.knots.case","num.knots.fprcase",
              #"Xfp","Zfp",
             "X.fprcase","Z.fprcase",
             "X.case","Z.case")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_seasonality_ncs <- function(...) {
  inits<-function(){list(beta=matrix(0,nrow=Jfull,ncol=2),
                         b=matrix(inits.b,nrow=Jfull,ncol=num.knots,byrow=TRUE),
                         taub=rep(0.01,Jfull),
                         betaEti=rbind(matrix(0,nrow=Jfull-1,ncol=num.knots.case+1),
                                       rep(NA,num.knots.case+1)))}
  #tauEti = rep(0.01,num.knots.case+1)
  parameters<-list("taub","sigmab","beta","b",
                   #"tauEti","sigmaEti",
                   "betaEti",
                   "pEti","Icat",
                   "tpr","mu","fpr.case","fpr.ctrl")
  data<-list("MBS","Nd","Nu","Jfull","template",
             "X.ctrl","Z.ctrl","num.knots","invsqrt",
             "Z.case","num.knots.case",
             "Xfp","Zfp")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_seasonality_ncs_comb <- function(...) {
  inits<-function(){list(beta=matrix(0,nrow=Jfull,ncol=2),
                         b=matrix(inits.b,nrow=Jfull,ncol=num.knots,byrow=TRUE),
                         taub=rep(0.01,Jfull),
                         betaEti=rbind(matrix(0,nrow=Jallowed-1,ncol=num.knots.case+1),
                                       rep(NA,num.knots.case+1)))}
  #tauEti = rep(0.01,num.knots.case+1)
  parameters<-list("taub","sigmab","beta","b",
                   #"tauEti","sigmaEti",
                   "betaEti",
                   "pEti","Icat",
                   "tpr","mu","fpr.case","fpr.ctrl")
  data<-list("MBS","Nd","Nu","Jfull","Jallowed","template",
             "X.ctrl","Z.ctrl","num.knots","invsqrt",
             "Z.case","num.knots.case",
             "Xfp","Zfp")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_seasonality_ncs_Deltafpr <- function(...) {
  inits<-function(){list(beta=matrix(0,nrow=Jfull,ncol=2),
                         b=matrix(inits.b,nrow=Jfull,ncol=num.knots,byrow=TRUE),
                         taub=rep(0.01,Jfull),
                         betaEti=rbind(matrix(0,nrow=Jfull-1,ncol=num.knots.case+1),
                                       rep(NA,num.knots.case+1)),
                         Delta = rep(0,Jfull),
                         tau_Delta = rep(0.01,Jfull))}
  #tauEti = rep(0.01,num.knots.case+1)
  parameters<-list("taub","sigmab","beta","b",
                   "Delta","tau_Delta","sigma_Delta",
                   "betaEti",
                   "pEti","Icat",
                   "tpr","mu","fpr.case","fpr.ctrl")
  data<-list("MBS","Nd","Nu","Jfull","template",
             "X.ctrl","Z.ctrl","num.knots","invsqrt",
             "Z.case","num.knots.case",
             "Xfp","Zfp","g")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_seasonality_ncs_tpr <- function(...) {
  inits<-function(){list(beta=matrix(0,nrow=Jfull,ncol=2),
                         b=matrix(inits.b,nrow=Jfull,ncol=num.knots,byrow=TRUE),
                         taub=rep(0.01,Jfull),
                         beta.tpr=matrix(0,nrow=Jfull,ncol=2),
                         b.tpr=matrix(inits.b.tpr,nrow=Jfull,ncol=num.knots.tpr,byrow=TRUE),
                         taub.tpr=rep(0.01,Jfull),
                         betaEti=rbind(matrix(0,nrow=Jfull-1,ncol=num.knots.case+1),
                                       rep(NA,num.knots.case+1)))}
  parameters<-list("taub","sigmab","beta","b",
                   #"tauEti","sigmaEti",
                   "betaEti",
                   "pEti","Icat",
                   "beta.tpr","b.tpr","sigmab.tpr","taub.tpr",
                   "tpr","mu","fpr.case","fpr.ctrl")
  data<-list("MBS","Nd","Nu","Jfull","template",
             "X.ctrl","Z.ctrl","num.knots","invsqrt",
             "Z.case","num.knots.case",
             "Z.tpr","X.tpr","num.knots.tpr",
             "Xfp","Zfp")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}



bugs_seasonality_ncs_re <- function(...) {
  inits<-function(){list(beta=matrix(0,nrow=Jfull,ncol=2),
                         b=matrix(inits.b,nrow=Jfull,ncol=num.knots,byrow=TRUE),
                         taub=rep(0.01,Jfull),
                         betaEti=rbind(matrix(0,nrow=Jfull-1,ncol=num.knots.case+1),
                                       rep(NA,num.knots.case+1)),
                         tau0=c(0.1,0.1))}
  parameters<-list("taub","sigmab","beta","b",
                   "betaEti",
                   "tau0","sigma0","b0",
                   "pEti","Icat",
                   "mu","logit.fpr.case","fpr.ctrl","logit.tpr",
                   "mustar")
  data<-list("MBS","Nd","Nu","Jfull","template",
             "X.ctrl","Z.ctrl","num.knots","invsqrt",
             "Z.case","num.knots.case",
             "Xfp","Zfp")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

bugs_seasonality_ncs_priorshape <- function(...) {
  inits<-function(){list(betaEti=rbind(matrix(0,nrow=Jfull-1,ncol=num.knots.case+1),
                                       rep(NA,num.knots.case+1)))}
  parameters<-list("pEti")
  data<-list("Nd","Jfull",
             "Z.case","num.knots.case")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}


bugs_seasonality_ncs_flex <- function(...) {
  inits<-function(){list(beta=matrix(0,nrow=Jfull,ncol=2),
                         b=matrix(inits.b,nrow=Jfull,ncol=num.knots,byrow=TRUE),
                         taub=rep(0.01,Jfull),
                         betaEti=rbind(matrix(0,nrow=Jfull,ncol=num.knots.case+1)))}
  #tauEti = rep(0.01,num.knots.case+1)
  parameters<-list("taub","sigmab","beta","b",
                   #"tauEti","sigmaEti",
                   "ind","ZEti",
                   "betaEti",
                   "tpr","mu","fpr.case","fpr.ctrl")
  data<-list("MBS","Nd","Nu","Jfull",
             "X.ctrl","Z.ctrl","num.knots","invsqrt",
             "Z.case","num.knots.case",
             "Xfp","Zfp")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}


bugs_seasonality_ncs_flex2 <- function(...) {
  inits<-function(){list(beta=matrix(0,nrow=Jfull,ncol=2),
                         b=matrix(inits.b,nrow=Jfull,ncol=num.knots,byrow=TRUE),
                         taub=rep(0.01,Jfull),
                         #tauV  = 0.1,
                         betaEti=array(0,c(Jfull,num.knots.case+1)))}
  #tauEti = rep(0.01,num.knots.case+1)
  parameters<-list(#"taub","sigmab","beta","b",
                   #"ind",
                   "Z",
                   "V",
                   "sigmaV",
                   #"betaEti",
                   #"tpr","mu",
                   "fpr.case","fpr.ctrl"
                   #"mutemp"
                   )
  data<-list("MBS","Nd","Nu","Jfull",
             "X.ctrl","Z.ctrl","num.knots","invsqrt",
             "Z.case","num.knots.case",
             "Xfp","Zfp","g")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}

#"alphaZ","betaZ",
bugs_seasonality_ncs_flex2_prior <- function(...) {
  inits<-function(){list(betaEti=array(0,c(Jfull,num.knots.case+1)))}
  #tauEti = rep(0.01,num.knots.case+1)
  parameters<-list("ind",
                   "Z",
                   "V")
  data<-list("Nd","Jfull",
             "Z.case","num.knots.case","alphaZ","betaZ","g")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}




bugs_seasonality_ncs_flex_feedback <- function(...) {
  inits<-function(){list(beta=matrix(0,nrow=Jfull,ncol=2),
                         b=matrix(inits.b,nrow=Jfull,ncol=num.knots,byrow=TRUE),
                         taub=rep(0.01,Jfull),
                         betaEti=array(0,c(Jfull,num.knots.case+1)),
                         theta_feed=0)}
  parameters<-list("taub","sigmab","beta","b",
    "ind","totalpathogen",
    "betaEti",
    "fpr.case","fpr.ctrl","theta_feed","tpr","mu","pEti_lim","pcov",
    "mu.rep","pEti_lim.rep","totalpathogen.rep","ind.rep")
  data<-list("MBS","Nd","Nu","Jfull",
             "X.ctrl","Z.ctrl","num.knots","invsqrt",
             "Z.case","num.knots.case",
             "Xfp","Zfp")
  rst.bugs   <- call.bugs(data, inits, parameters,...);
  rst.bugs;
}




###############################################################################
###############################################################################
###############################################################################
############   END OF BUGS CALL FUNCTIONS                ######################
###############################################################################
###############################################################################
###############################################################################

##------------------------------------------------------------------------##
## Extract data on pathogens from specimen and test. 
## It handles data before combination.
## MeasDir tells the function where the data sets (.csv) are.
## PathCatDir tells the function where the pathogen category .csv is.
##------------------------------------------------------------------------##
extract_data_raw = function(Pathogen,Specimen,Test,X,Xval,MeasDir,PathCatDir){
  pathogen_type = read.csv(PathCatDir)
  rownames(pathogen_type) = pathogen_type[,1]
  typeOrder = order(pathogen_type[Pathogen,2])
  Pathogen = Pathogen[typeOrder]
  #pathogen_type[Pathogen,]

  datraw = read.csv(MeasDir)
  #clean column names if the column names start with "_":
  delete_start_with = function(s,vec){
    ind = grep(s,substring(vec,1,nchar(s)))
    old = vec[ind]
    vec[ind] = substring(old,nchar(s)+1)
    return(vec)
  }
  cleanName=delete_start_with("X_",names(datraw))
  dat0 = datraw
  colnames(dat0) = cleanName
  
  indX = 1:nrow(dat0)
  for (j in 1:length(X)){
      indX = indX[which(dat0[indX,X[j]]==Xval[j] & !is.na(dat0[indX,X[j]]))]
  }
  dat = dat0[indX,]
  
  pstGrid = apply(expand.grid(paste(Pathogen,"_",sep=""),Specimen,Test),
                           1,paste,collapse="")
  
  # the output will be for each (specimen, test) pair, with pathogens aligned
  stGrid = apply(expand.grid(Specimen,Test),1,paste,collapse="")
  
  pstTable = matrix(NA,length(Pathogen),length(stGrid))
  colnames(pstTable) = stGrid
  rownames(pstTable) = Pathogen
   
  for (i in 1:length(stGrid)){
    for (j in 1:length(Pathogen)){
      tempName = paste(Pathogen[j],stGrid[i],sep="_")
      if (tempName%in%cleanName){
        if (sum(is.na(dat[,tempName]))!=nrow(dat)){
          pstTable[j,i] = TRUE
        }
      }
    }
  }

  
  notindata = which(rowSums(!is.na(pstTable))==0)
  if (length(notindata)>0) {cat(Pathogen[notindata],"can't be found in the dataset!")}
  
  naColumns = which(colSums(is.na(pstTable))==length(Pathogen))
  if (length(naColumns)==length(stGrid)) {
    cat("No test has available results on selected pathogens! Try other pathogens.")
    return(-1)
  } else{
    if (length(naColumns)==0){
      actualpstTable = pstTable
    } else{
      actualpstTable = pstTable[,-naColumns] 
    }
    
  resdat = list()
  
  for (j in 1:ncol(actualpstTable)){
      if (sum(is.na(actualpstTable[,j]))==0){
        tempnm = paste(Pathogen,colnames(actualpstTable)[j],sep="_")
        resdat[[j]]=dat[,tempnm] 
      } else{
        dftemp = as.data.frame(matrix(NA,nrow=nrow(dat),ncol=nrow(actualpstTable)))
        colnames(dftemp) = paste(Pathogen,colnames(actualpstTable)[j],sep="_")
        dftemp[,!is.na(actualpstTable[,j])] = 
          dat[,paste(Pathogen[!is.na(actualpstTable[,j])],colnames(actualpstTable)[j],sep="_")]
        resdat[[j]]=dftemp
      }
  }
    resdat[[ncol(actualpstTable)+1]]=dat[,"patid"]
    resdat[[ncol(actualpstTable)+2]]=dat[,"ENRLDATE"]
    resdat[[ncol(actualpstTable)+3]]=dat[,"AGECAT"]
    resdat[[ncol(actualpstTable)+4]]=dat[,"HIV"]
    #resdat[[ncol(actualpstTable)+5]]=dat[,"PRABXBC"]
    #names(resdat) = c(colnames(actualpstTable),"patid","ENRLDATE","AGECAT","HIV","PRABXBC")
     names(resdat) = c(colnames(actualpstTable),"patid","ENRLDATE","AGECAT","HIV")
 
  return(resdat)
  }
}

#####################################
##  example of usage
####################################
# Pathogen = c("PNEU","PV_EV","SASP","PARA1","HINF","INFLUENZA_A","LEGIO")
# Specimen = c("NP","B")
# Test = c("PCR","CX")
# 
# X = c("SITE","CASECONT","AGECAT")
# Xval = c("01KEN",1,1)
# 
# # need to build in pathogen category to sort measurement
# # need to build in for DataType = "raw" or "com"
# MeasDir = ifelse(as.numeric(Sys.info()["sysname"]=="Windows"),
#                  paste0("C:/Users/Administrator/",
#                         "Dropbox/ZW/working_projects/PERCH/data/EC2013LondonDataNov/PQ.csv"))
# PathCatDir = ifelse(as.numeric(Sys.info()["sysname"]=="Windows"),
#                     paste0("C:/Users/Administrator/",
#                            "Dropbox/ZW/working_projects/PERCH/data/EC2013LondonDataNov/pathogen_category_LH.csv"))
# 
# mydat = extract_data_raw(Pathogen,Specimen,Test,X,Xval,MeasDir,PathCatDir)


##------------------------------------------------------------------------##
## This function extracts GS/BS measure for each specified pathogen.
## Each element of the list will from BS or GS, with columns being pathogens.
## Also, another part of the list will be sens/spec for BS or GS
##------------------------------------------------------------------------##

extract_data_comb = function(Pathogen,Test,X,Xval,MeasDir,PathCatDir){
  pathogen_type = read.csv(PathCatDir)
  rownames(pathogen_type) = pathogen_type[,1]
  typeOrder = order(pathogen_type[Pathogen,2])
  Pathogen = Pathogen[typeOrder]
  pathogen_type[Pathogen,]
  
  datraw = read.csv(MeasDir)
  
  #clean column names if the column names start with "_":
  delete_start_with = function(s,vec){
    ind = grep(s,substring(vec,1,nchar(s)))
    old = vec[ind]
    vec[ind] = substring(old,nchar(s)+1)
    return(vec)
  }
  cleanName=delete_start_with("X_",names(datraw))
  dat0 = datraw
  colnames(dat0) = cleanName
  
  indX = 1:nrow(dat0)
  for (j in 1:length(X)){
    indX = indX[which(dat0[indX,X[j]]==Xval[j])]
  }
  dat = dat0[indX,]
  
  ptGrid = apply(expand.grid(paste(Pathogen,"_",sep=""),Test),
                 1,paste,collapse="")
  
  ptTable = matrix(NA,length(Pathogen),length(Test))
  colnames(ptTable) = Test
  rownames(ptTable) = Pathogen
  
  for (i in 1:length(Test)){
    for (j in 1:length(Pathogen)){
      tempName = paste(Pathogen[j],Test[i],sep="_")
      if (tempName%in%cleanName){
        if (sum(is.na(dat[,tempName]))!=nrow(dat)){
          ptTable[j,i] = TRUE
        }
      }
    }
  }
  
  notindata = which(rowSums(!is.na(ptTable))==0)
  if (length(notindata)>0) {cat(Pathogen[notindata],"can't be found in the dataset!")}
  
  naColumns = which(colSums(is.na(ptTable))==length(Pathogen))
  if (length(naColumns)==length(Test)) {
    cat("No test has available results on selected pathogens! Try other pathogens.")
    return(-1)
  } else{
    if (length(naColumns)==0){
      actualptTable = ptTable
    } else{
      actualptTable = ptTable[,-naColumns] 
    }
  
    resdat = list()
    
    for (j in 1:ncol(actualptTable)){
      if (sum(is.na(actualptTable[,j]))==0){
        tempnm = paste(Pathogen,colnames(actualptTable)[j],sep="_")
        resdat[[j]]=dat[,tempnm] 
      } else{
        dftemp = as.data.frame(matrix(NA,nrow=nrow(dat),ncol=nrow(actualptTable)))
        colnames(dftemp) = paste(Pathogen,colnames(actualptTable)[j],sep="_")
        dftemp[,!is.na(actualptTable[,j])] = 
          dat[,paste(Pathogen[!is.na(actualptTable[,j])],colnames(actualptTable)[j],sep="_")]
        resdat[[j]]=dftemp
      }
    }
    names(resdat) = colnames(actualptTable)
    return(resdat)
  }
}

#####################################
##  example of usage
####################################
# Pathogen = c("PNEU","PV_EV","SASP","PARA1","HINF","INFLUENZA_A","LEGIO")
# Test = c("BS","GS","BS_n","GS_n","BS_SEN","BS_SEN_MIN","BS_SEN_MAX","BS_SPE","BS_SPE_MIN","BS_SPE_MAX",
#          "GS_SEN","GS_SEN_MIN","GS_SEN_MAX","GS_SPE","GS_SPE_MIN","GS_SPE_MAX")
# 
# X = c("SITE","CASECONT","AGECAT")
# Xval = c("01KEN",1,1)
# 
# MeasDir= ifelse(as.numeric(Sys.info()["sysname"]=="Windows"),
#                 paste0("C:/Users/Administrator/",
#                        "Dropbox/ZW/working_projects/PERCH/data/EC2013LondonDataNov/PQ\ datasets/Senspec_before/PQ_EC_sero_24NOV13.csv"))
# PathCatDir = ifelse(as.numeric(Sys.info()["sysname"]=="Windows"),
#                     paste0("C:/Users/Administrator/",
#                            "Dropbox/ZW/working_projects/PERCH/data/EC2013LondonDataNov/pathogen_category_LH.csv"))
# 
# extract_data_comb(Pathogen,Test,X,Xval,MeasDir,PathCatDir)


##------------------------------------------------------------------------##
# Display pathogen and test pairs, no. of observations stratified by covariates,
# and cleaned dataset
##------------------------------------------------------------------------##

display_pathogen_test = function(ecdataraw,specimenName = c("B","PF","LA","NP","IS"),
                                 testName = c("CX","CX2","PCR")){
  #clean column names
  delete_start_with = function(s,vec){
    ind = grep(s,substring(vec,1,nchar(s)))
    old = vec[ind]
    vec[ind] = substring(old,nchar(s)+1)
    return(vec)
  }
  cleanName=delete_start_with("X_",names(ecdataraw))
  ecdata = ecdataraw
  colnames(ecdata) = cleanName
  
  #display all pathogens measure from a specific specimen-test pair
  display_pathogen = function(SpecTest,vec){
    PathTest = vec[grep(SpecTest,
                        substring(vec,nchar(vec)-nchar(SpecTest)+1,nchar(vec)))]
    res=substring(PathTest,1,nchar(PathTest)-nchar(SpecTest)-1)
    return(res)
  }
  
  #obtain pathogen-specimen-test table
  
  SpecimenTestGrid = apply(expand.grid(specimenName,testName),
                           1,paste,collapse="")
  
  pathDisplayList=sapply(SpecimenTestGrid,display_pathogen,cleanName)
  
  pathogens = sort(unique(unlist(pathDisplayList)))
  
  pathSpecTestTable = matrix(NA,length(pathogens),length(SpecimenTestGrid))
  colnames(pathSpecTestTable) = SpecimenTestGrid
  rownames(pathSpecTestTable) = pathogens
  for (i in 1:length(SpecimenTestGrid)){
    st = SpecimenTestGrid[i]
    temp = pathDisplayList[[st]]
    matchind = sapply(temp,function(v) which(pathogens==v))
    if (length(matchind)!=0){
      pathSpecTestTable[matchind,i] = "Yes"
    }
  }
  
  #reorder data by pathogen types, has to be specified for program input
  pathogen_type = read.csv(PathCatDir)
  rownames(pathogen_type) = pathogen_type[,1]
  typeOrder = order(pathogen_type[pathogens,2])
  pathSpecTestTable = pathSpecTestTable[typeOrder,]
  pathogens = pathogens[typeOrder]
  

  #up to now, all that is done is about finding path-spec-test triples that are in the dataset
  #allCaseType=sort(unique(ecdata$CASECONT))
  #caseInd = sapply(allCaseType,function(v) which(ecdata$CASECONT==v))

  # obsSpecTestName contains those (specimen,test) pairs that BOTH 1) have the name
  # appearing in the columns of raw data AND 2)the data below these names in the columns
  # are not all NA; obsSpecTestName1 contains those (specimen,test) pairs that 
  # just satisfy criteria 1)
  # obsSpecimenTestGrid: (pathogen,specimen_test) triples that satisfy criteria 1 and 2
  
    ecdatatemp = ecdata
    
    #get rid of SpecTest that either 1. not in the dataset or 2.in the dataset
    # but all are NA valued:
    notindata = sapply(SpecimenTestGrid,
                       function(v) sum(!is.na(ecdatatemp[,grep(v,cleanName)])))
    #can be many columns in the raw PQ.csv data. Pathogen_SpecTest.
    
    naColumns = which(colSums(is.na(pathSpecTestTable))==length(pathogens) | notindata==0 )
    
    obsSpecTestName = SpecimenTestGrid[-naColumns]
    
    pstTableNA = cbind(as.character(pathogen_type[pathogens,2]),pathSpecTestTable[,-naColumns])
    colnames(pstTableNA)= c("type",obsSpecTestName)
    pstTable = pstTableNA
    
    #write.csv(pstTable,"pstTable.csv")
    pstTable[is.na(pstTable)]="."
    pstTableNoType = pstTable[,-1]
    # write the table that shows (specimen,test) availability for each pathogen
    
    ## remove those columns that never appears in the columns of raw data:
    naColumns1 = which(colSums(is.na(pathSpecTestTable))==length(pathogens))
    
    obsSpecTestName1 = SpecimenTestGrid[-naColumns1]  
    pstTableNA1 = cbind(as.character(pathogen_type[typeOrder,2]),pathSpecTestTable[,-naColumns1])
    colnames(pstTableNA1)= c("type",obsSpecTestName1)
    pstTable1 = pstTableNA1                   
    pstTable1[is.na(pstTable1)]="."
    pstTableNoType1 = pstTable1[,-1] 
    obsSpecimenTestGrid1 = outer(pathogens,obsSpecTestName1,paste,sep="_")[(pstTableNoType1!=".")]
    
    #generate individual measurement profiles:
    obsSpecimenTestGrid = outer(pathogens,obsSpecTestName,paste,sep="_")[(pstTableNoType!=".")]
    M = array(NA,c(length(pathogens),length(obsSpecTestName),nrow(ecdatatemp)),
              dimnames=list(pathogens,obsSpecTestName,1:nrow(ecdatatemp)))
    
    progbar = txtProgressBar(min = 0, max = nrow(ecdatatemp), style = 3)
    for (i in 1:nrow(ecdatatemp)){
      tempmat = matrix(NA,nrow=length(pathogens),ncol=length(obsSpecTestName))
      tempmat[pstTableNoType=="Yes"] = as.numeric(ecdatatemp[i,obsSpecimenTestGrid])
      M[,,i] = tempmat 
      if(i==1){cat("-------Conversion to Individual Measurement Profiles--------","\n")
               cat("Percentage of Conversion:","\n")}
      setTxtProgressBar(progbar, i)  
    }
    close(progbar)
    
    MView = M
    MView[is.na(M)]="."
    
    
    #get proportion of data available for each path-spec-test triple
    propObs = round(100*apply(M,c(1,2),function(v) 1-sum(is.na(v))/nrow(ecdatatemp)),2)
    pstTableNoType[pstTableNoType=="Yes"]=paste(pstTableNoType[pstTableNoType=="Yes"], propObs[pstTableNoType=="Yes"],sep="|")
    pstTable[,-1] = pstTableNoType
    
    xNames = setdiff(cleanName,c(obsSpecimenTestGrid1))
    ecdatatempClean = ecdatatemp[,c(xNames,obsSpecimenTestGrid)]
 
    res = list(nrow(ecdatatemp),pstTable,ecdatatempClean)
    names(res) = c("nobs","pstTable","cleanedData")
    return(res)
}

# # ## example of usage
#  MeasDir = ifelse(as.numeric(Sys.info()["sysname"]=="Windows"),
#                   paste0("C:/Users/Administrator/",
#                          "Dropbox/ZW/working_projects/PERCH/data/EC2013LondonDataNov/PQ.csv"))
#  PathCatDir = ifelse(as.numeric(Sys.info()["sysname"]=="Windows"),
#                      paste0("C:/Users/Administrator/",
#                             "Dropbox/ZW/working_projects/PERCH/data/EC2013LondonDataNov/pathogen_category_LH.csv"))  
# # 
# ecdataraw      = read.csv(MeasDir)
# # display_pathogen_test(ecdataraw)
# disp = display_pathogen_test(subset(ecdataraw,SITE=="01KEN" & X_AGECAT==1 & X_CASECONT==1))
# disp[1:2]



##---------------------------------------------------------------------------##
## simple functions
##---------------------------------------------------------------------------##
make.foldername=function(parent.path,parameter.names,parameter.vals){
  subfolder = paste(parameter.names,parameter.vals,collapse="_",sep="=")
  res = paste(parent.path,subfolder,sep="\\")
  return(res)
}

make.filename = function(parameter.names,parameter.vals,format){
  res1=paste(parameter.names,parameter.vals,collapse="_",sep="=")
  res = paste(res1,format,sep=".")
  return(res)
}

# write ParameterList into .txt file
fnlist <- function(x, fil){ z <- deparse(substitute(x))
                            cat(z, "\n", file=fil)
                            nams=names(x) 
                            for (i in seq_along(x) ){ cat(nams[i], "\t",  x[[i]], "\n", 
                                                          file=fil, append=TRUE) }
}

logit = function(x) log(x)-log(1-x)
expit = function(x) 1/(1+exp(-x))


LOR<- function(x,y){
  require(epitools)
  res = oddsratio(x,y)
  res.or = res$measure[2,]
  res.p = res$p.value[2,"fisher.exact"]
  resl = list(res.or,res.p)
  names(resl) = c("OR","pval")
  resl
}


rvbern<-function(p){ 
  # fast sample multiple bernoulli variables with different means
  U<-runif(length(p),0,1)
  res<-(U<p)+0
  return(res)
}


bin2dec <- function(binaryvector) { 
  ##convert 0/1 binary coded sequence into decimal digits
  sum(2^(which(rev(binaryvector)==TRUE)-1))
}

symb2I = function(instance,pathlist){
  #this function convert names of pathogen/combinations into 0/1 coding
  #instance = i_L
  J = length(pathlist)
  splited = strsplit(instance,split="+",fixed=TRUE)
  deploy = function(inst,J) {
    if (inst[1] == "NoA"){
      rep(0,J)
    }else{
      sapply(inst,grep,x=pathlist)
    }
  }
  res = lapply(splited,deploy,J=J)
  nc = length(res)
  matres = t(sapply(1:nc,function(i) {
    tempres = rep(0,J)
    tempres[res[[i]]]=1
    tempres}))
  return(matres)
}
#symb2I("A+B",c("A","B","C"))
#symb2I("NoA",c("A","B","C"))

I2symb = function(code,pathlist){
  #this functio convert 0/1 coding to pathogen/combinations
  ind = grep("1",strsplit(code,split="")[[1]])
  res = ifelse(length(ind)==0,"NoA",paste(pathlist[ind],collapse="+"))
  return(res)
}

Imat2cat= function(mat,allowedlist,pathlist){# for case only
  known_code = apply(mat,1,function(v) paste(v,collapse=""))
  known_symb = sapply(known_code,I2symb,pathlist)
  known_Icat = sapply(known_symb,function(s) which(allowedlist==s))
  return(known_Icat)
}



##---------------------------------------------------------------------------##
## simulate data from PQ model specification, without covariates
## 
##---------------------------------------------------------------------------##
data_simu_nppcr = function(ParameterList){
  require(mvtnorm)
  require(e1071)
  FullPathogenName = ParameterList$FullPathogenName
  AllowedList      = ParameterList$AllowedList # including combinations
  Jfull = length(FullPathogenName)
  Jallowed = length(AllowedList)
  
  BSmodel = ParameterList$BSmodel
  
  psiBS = ParameterList$psiBS  
  thetaBS = ParameterList$thetaBS
  psiGS = ParameterList$psiGS  
  thetaGS = ParameterList$thetaGS
  
  Sigma0d = ParameterList$Sigma0d
  Sigma0u = ParameterList$Sigma0u
  Sigma1d = ParameterList$Sigma1d
  Sigma1u = ParameterList$Sigma1u
  
  Nd      =     ParameterList$Nd
  Nu      =     ParameterList$Nu
  
  #plist0  = ParameterList$plist0
  plist1  = ParameterList$plist1
  GSfraction = ParameterList$GSfraction
  etiology = ParameterList$etiology
  
  iLcat = rep(NA,Nd)
  iLall = matrix(NA,nrow=Nd+Nu,ncol=Jfull)
  etiologyMat = matrix(NA,nrow=Nd,ncol=Jallowed)
  for (k in 1:Nd){
    etiologyMat[k,] = etiology
    iLcat[k]   = sample(AllowedList,1,prob = etiologyMat[k,])
  }
  iL = symb2I(iLcat,FullPathogenName)
  iLall = rbind(iL,matrix(0,nrow=Nu,ncol=Jfull))
  iLcatAllnumeric=c(Imat2cat(iL,AllowedList,FullPathogenName),
                    rep(Jallowed+1,Nu))
  
  
  if (BSmodel=="FC1"){
    ## bronze standard (single, NP PCR)
    # random effects by partition
    b0d = rnorm(Nd,0,Sigma0d)
    b0u = rnorm(Nu,0,Sigma0u)
    b1d = rmvnorm(Nd,rep(0,length(plist1)),Sigma1d)
    b1u = rmvnorm(Nu,rep(0,length(plist1)),Sigma1u)
    
    MpdBS = matrix(NA,nrow=Nd,ncol=Jfull)# probabilities
    MpuBS = matrix(NA,nrow=Nu,ncol=Jfull)# probabilities
    MdBS = matrix(NA,nrow=Nd,ncol=Jfull) #binaries
    MuBS = matrix(NA,nrow=Nu,ncol=Jfull) #binaries
    
    
    for (k in 1:Nd){
      for (j in 1:Jfull){
        MpdBS[k,j] = expit(b0d[k]+b1d[k,vnu(j,plist1,FullPathogenName)]+
                             psiBS[j]*(1-iL[k,j])+thetaBS[j]*iL[k,j])
      }
      MdBS[k,]  = rvbern(MpdBS[k,])
    }
    
    for (k in 1:Nu){
      for (j in 1:Jfull){
        MpuBS[k,j] = expit(b0u[k]+b1u[k,vnu(j,plist1,FullPathogenName)]+
                             psiBS[j])
      }
      MuBS[k,] = rvbern(MpuBS[k,])
    }
  } else if (BSmodel=="CI"){
    ## assume conditional independence
    assign_NP = function(I,theta,psi){
      J =length(I)
      res = rep(NA,J)
      for (j in 1:J){
        res[j] = rbinom(1,size=1,ifelse(I[j]==1,theta[j],psi[j]))
      }
      return(res)
    }
    MdBS = t(apply(iL,1,assign_NP,theta=thetaBS,psi=psiBS))
    MuBS = t(apply(iLall[-(1:Nd),],1,assign_NP,theta=thetaBS,psi=psiBS))
    
  } else if (BSmodel=="FC0"){
    ## only frailty intercept
    b0d = rnorm(Nd,0,Sigma0d)
    b0u = rnorm(Nu,0,Sigma0u)
    
    
    MpdBS = matrix(NA,nrow=Nd,ncol=Jfull)# probabilities
    MpuBS = matrix(NA,nrow=Nu,ncol=Jfull)# probabilities
    MdBS = matrix(NA,nrow=Nd,ncol=Jfull) #binaries
    MuBS = matrix(NA,nrow=Nu,ncol=Jfull) #binaries
    for (k in 1:Nd){
      for (j in 1:Jfull){
        MpdBS[k,j] = expit(b0d[k]+psiBS[j]*(1-iL[k,j])+thetaBS[j]*iL[k,j])
      }
      MdBS[k,]  = rvbern(MpdBS[k,])
    }
    
    for (k in 1:Nu){
      for (j in 1:Jfull){
        MpuBS[k,j] = expit(b0u[k]+psiBS[j])
      }
      MuBS[k,] = rvbern(MpuBS[k,])
    }
  }
  
  
  
  ## gold standard data
  NdGS = floor(GSfraction*Nd)
  MpdGScomplete = matrix(NA,nrow=Nd,ncol=Jfull)
  MdGScomplete = matrix(NA,nrow=Nd,ncol=Jfull)
  
  if ( NdGS != 0){
    for (k in 1:Nd){
      MpdGScomplete[k,] = psiGS*(1-iL[k,])+thetaGS*iL[k,]
      MdGScomplete[k,]  = rvbern(MpdGScomplete[k,])
    }
    MdGSobs=MdGScomplete[1:NdGS,]
    MGSMiss = rbind(MdGSobs,matrix(NA,nrow=Nd+Nu-NdGS,ncol=Jfull))
  } else{
    MGSMiss = matrix(NA,nrow=Nu+Nd,ncol=Jfull)
  }
  
  
  ## organize case/control status, iL, BS, GS data into dataframes
  resCompleteGS = data.frame(Y = c(rep(1,Nd),rep(0,Nu)),
                             iLcat = iLcatAllnumeric,
                             iL = iLall,
                             MBS = rbind(MdBS,MuBS),
                             MGS = rbind(MdGScomplete,matrix(NA,nrow=Nu,ncol=Jfull)))
  resMissGS      = data.frame(Y = c(rep(1,Nd),rep(0,Nu)),
                              iLcat = iLcatAllnumeric,
                              iL = iLall,
                              MBS = rbind(MdBS,MuBS),
                              MGS = MGSMiss)
  datacolnames    = c("Y","iLcat",
                      paste("iL",FullPathogenName,sep="_"),
                      paste("MBS",FullPathogenName,sep="_"),
                      paste("MGS",FullPathogenName,sep="_"))
  colnames(resCompleteGS) = datacolnames
  colnames(resMissGS) = datacolnames
  template = as.matrix(rbind(symb2I(AllowedList,FullPathogenName),rep(0,Jfull)))
  colnames(template) = FullPathogenName
  rownames(template) = c(AllowedList,"control")
  return(list(BSmodel=BSmodel,template = template,NdGS=NdGS,
              dataComplete=resCompleteGS,dataMiss = resMissGS))
}


##---------------------------------------------------------------------------##
##simulate data with both blood culture and NP PCR
##---------------------------------------------------------------------------##
data_simu_bcx = function(ParameterList){
  require(mvtnorm)
  require(e1071)
  FullPathogenName = ParameterList$FullPathogenName
  AllowedList      = ParameterList$AllowedList # including combinations
  Jfull = length(FullPathogenName)
  Jallowed = length(AllowedList)
  
  BSmodel = ParameterList$BSmodel
  
  psiBS = ParameterList$psiBS  
  thetaBS = ParameterList$thetaBS
  psiGS = ParameterList$psiGS  
  thetaGS = ParameterList$thetaGS
  
  Nd      =     ParameterList$Nd
  Nu      =     ParameterList$Nu
  
  bcxList = ParameterList$bcxList
  JGS = length(bcxList)
  
  GSfraction = ParameterList$GSfraction
  
  etiology = ParameterList$etiology
  
  iLcat = rep(NA,Nd)
  iLall = matrix(NA,nrow=Nd+Nu,ncol=Jfull)
  etiologyMat = matrix(NA,nrow=Nd,ncol=Jallowed)
  for (k in 1:Nd){
    etiologyMat[k,] = etiology
    iLcat[k]   = sample(AllowedList,1,prob = etiologyMat[k,])
  }
  iL = symb2I(iLcat,FullPathogenName)
  iLall = rbind(iL,matrix(0,nrow=Nu,ncol=Jfull))
  iLcatAllnumeric=c(Imat2cat(iL,AllowedList,FullPathogenName),
                    rep(Jallowed+1,Nu))
  
  ## assume conditional independence
  assign_NP = function(I,theta,psi){
    J =length(I)
    res = rep(NA,J)
    for (j in 1:J){
      res[j] = rbinom(1,size=1,ifelse(I[j]==1,theta[j],psi[j]))
    }
    return(res)
  }
  MdBS = t(apply(iL,1,assign_NP,theta=thetaBS,psi=psiBS))
  MuBS = t(apply(iLall[-(1:Nd),],1,assign_NP,theta=thetaBS,psi=psiBS))
  
  ## gold standard data
  NdGS = floor(GSfraction*Nd)
  MpdGScomplete = matrix(NA,nrow=Nd,ncol=Jfull)
  MdGScomplete = matrix(NA,nrow=Nd,ncol=Jfull)
  
  if ( NdGS != 0){
    for (k in 1:Nd){
      MpdGScomplete[k,1:JGS] = psiGS*(1-iL[k,1:JGS])+thetaGS*iL[k,1:JGS]
      MdGScomplete[k,1:JGS]  = rvbern(MpdGScomplete[k,1:JGS])
    }
    MdGSobs=MdGScomplete[1:NdGS,]
    MGSMiss = rbind(MdGSobs,matrix(NA,nrow=Nd+Nu-NdGS,ncol=Jfull))
  } else{
    MGSMiss = matrix(NA,nrow=Nu+Nd,ncol=Jfull)
  }
  
  MGSMiss[,-match(bcxList,FullPathogenName)] = rep(NA,Nd+Nu)
  
  
  ## organize case/control status, iL, BS, GS data into dataframes
  resCompleteGS = data.frame(Y = c(rep(1,Nd),rep(0,Nu)),
                             iLcat = iLcatAllnumeric,
                             iL = iLall,
                             MBS = rbind(MdBS,MuBS),
                             MGS = rbind(MdGScomplete,matrix(NA,nrow=Nu,ncol=Jfull)))
  resMissGS      = data.frame(Y = c(rep(1,Nd),rep(0,Nu)),
                              iLcat = iLcatAllnumeric,
                              iL = iLall,
                              MBS = rbind(MdBS,MuBS),
                              MGS = MGSMiss)
  datacolnames    = c("Y","iLcat",
                      paste("iL",FullPathogenName,sep="_"),
                      paste("MBS",FullPathogenName,sep="_"),
                      paste("MGS",FullPathogenName,sep="_"))
  colnames(resCompleteGS) = datacolnames
  colnames(resMissGS) = datacolnames
  template = as.matrix(rbind(symb2I(AllowedList,FullPathogenName),rep(0,Jfull)))
  colnames(template) = FullPathogenName
  rownames(template) = c(AllowedList,"control")
  return(list(BSmodel=BSmodel,template = template,NdGS=NdGS,
              dataComplete=resCompleteGS,dataMiss = resMissGS))
}
##----------------------------------------------------------------------------##
## obtain beta parameters
##----------------------------------------------------------------------------##
beta.parms.from.quantiles <- function(q, p=c(0.025,0.975),
                                      precision=0.001, derivative.epsilon=1e-3,
                                      start.with.normal.approx=T, start=c(1, 1), plot=F)
{
  # Version 1.2.2 (December 2012)
  #
  # Function developed by 
  # Lawrence Joseph and Patrick Belisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@clinepi.mcgill.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html
  #
  # Please refer to our webpage for details on each argument.
  
  f <- function(x, theta){dbeta(x, shape1=theta[1], shape2=theta[2])}
  F.inv <- function(x, theta){qbeta(x, shape1=theta[1], shape2=theta[2])}
  f.cum <- function(x, theta){pbeta(x, shape1=theta[1], shape2=theta[2])}
  f.mode <- function(theta){a <- theta[1]; b <- theta[2]; mode <- ifelse(a>1, (a-1)/(a+b-2), NA); mode}
  theta.from.moments <- function(m, v){a <- m*m*(1-m)/v-m; b <- a*(1/m-1); c(a, b)}
  plot.xlim <- c(0, 1)
  
  dens.label <- 'dbeta'
  parms.names <- c('a', 'b')
  
  if (length(p) != 2) stop("Vector of probabilities p must be of length 2.")
  if (length(q) != 2) stop("Vector of quantiles q must be of length 2.")
  p <- sort(p); q <- sort(q)
  
  #_____________________________________________________________________________________________________
  
  print.area.text <- function(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim, M=30, M0=50)
  {
    par.usr <- par('usr')
    par.din <- par('din')
    
    p.string <- as.character(round(c(0,1) + c(1,-1)*p.check, digits=4))
    str.width <- strwidth(p.string, cex=cex)
    str.height <- strheight("0", cex=cex)
    
    J <- matrix(1, nrow=M0, ncol=1)
    
    x.units.1in <- diff(par.usr[c(1,2)])/par.din[1]
    y.units.1in <- diff(par.usr[c(3,4)])/par.din[2]
    aspect.ratio <- y.units.1in/x.units.1in
    
    # --- left area  -----------------------------------------------------------
    
    scatter.xlim <- c(max(plot.xlim[1], par.usr[1]), q[1])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=0, to=p[1], length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-1]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the right of the mode, if any
    w <- which(x>mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[1]+str.width[1]) <= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y, fromLast=T))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[1], adj=c(1,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[1], mean(par.usr[c(3,4)]), labels=p.string[1], col='gray', cex=cex, srt=90, adj=c(1,0))
    }
    
    # --- right area  ----------------------------------------------------------
    
    scatter.xlim <- c(q[2], plot.xlim[2])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=p[2], to=f.cum(plot.xlim[2], theta), length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-length(h)]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the left of the mode, if any
    w <- which(x<mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[2]-str.width[2]) >= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[2], adj=c(0,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[2], mean(par.usr[c(3,4)]), labels=p.string[2], col='gray', cex=cex, srt=-90, adj=c(1,0))
    }
  }
  
  # ......................................................................................................................................
  
  Newton.Raphson <- function(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start)
  {
    Hessian <- matrix(NA, 2, 2)
    
    if (start.with.normal.approx)
    {
      # Probably not a very good universal choice, but proved good in most cases in practice
      m <-  diff(q)/diff(p)*(0.5-p[1]) + q[1]
      v <- (diff(q)/diff(qnorm(p)))^2
      theta <- theta.from.moments(m, v)
    }
    else theta <- start
    
    
    change <- precision + 1
    niter <- 0
    # Newton-Raphson multivariate algorithm
    while (max(abs(change)) > precision)
    {
      Hessian[,1] <- (f.cum(q, theta) - f.cum(q, theta - c(derivative.epsilon, 0))) / derivative.epsilon
      Hessian[,2] <- (f.cum(q, theta) - f.cum(q, theta - c(0, derivative.epsilon))) / derivative.epsilon
      
      f <- f.cum(q, theta) - p
      change <- solve(Hessian) %*% f
      last.theta <- theta
      theta <- last.theta - change
      
      # If we step out of limits, reduce change
      
      if (any(theta<0))
      {
        k <- min(last.theta/change)
        theta <- last.theta - k/2*change
      }
      
      niter <- niter + 1
    }
    
    list(theta=as.vector(theta), niter=niter, last.change=as.vector(change))
  }
  
  # ...............................................................................................................
  
  plot.density <- function(p, q, f, f.cum, F.inv, mode, theta, plot.xlim, dens.label, parms.names, cex)
  {
    if (length(plot.xlim) == 0)
    {
      plot.xlim <- F.inv(c(0, 1), theta)
      
      if (is.infinite(plot.xlim[1]))
      {
        tmp <- min(c(0.001, p[1]/10))
        plot.xlim[1] <- F.inv(tmp, theta)
      }  
      
      if (is.infinite(plot.xlim[2]))
      {
        tmp <- max(c(0.999, 1 - (1-p[2])/10))
        plot.xlim[2] <- F.inv(tmp, theta)
      }
    }
    plot.xlim <- sort(plot.xlim)
    
    
    x <- seq(from=min(plot.xlim), to=max(plot.xlim), length=1000)
    h <- f(x, theta)
    x0 <- x; f0 <- h
    ylab <- paste(c(dens.label, '(x, ', parms.names[1], ' = ', round(theta[1], digits=5), ', ', parms.names[2], ' = ', round(theta[2], digits=5), ')'), collapse='')
    plot(x, h, type='l', ylab=ylab)
    
    # fill in area on the left side of the distribution
    x <- seq(from=plot.xlim[1], to=q[1], length=1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from=max(plot.xlim), to=q[2], length=1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # draw distrn again
    points(x0, f0, type='l')
    h <- f(q, theta)
    points(rep(q[1], 2), c(0, h[1]), type='l', col='orange')
    points(rep(q[2], 2), c(0, h[2]), type='l', col='orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)  
    
    xaxp <- par("xaxp")
    x.ticks <- seq(from=xaxp[1], to=xaxp[2], length=xaxp[3]+1)
    q2print <- as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    mtext(q2print, side=1, col='orange', at=q2print, cex=0.6, line=2.1)
    points(q, rep(par('usr')[3]+0.15*par('cxy')[2], 2), pch=17, col='orange')
  }
  
  #________________________________________________________________________________________________________________
  
  
  parms <- Newton.Raphson(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start=start)
  p.check <- f.cum(q, parms$theta)
  
  if (plot) plot.density(p, q, f, f.cum, F.inv, f.mode(parms$theta), parms$theta, plot.xlim, dens.label, parms.names, 0.8)
  
  list(a=parms$theta[1], b=parms$theta[2], last.change=parms$last.change, niter=parms$niter, q=q, p=p, p.check=p.check)
}



betaplot=function(a,b){
  x= seq(0,1,by=0.001)
  y = dbeta(x,a,b)
  plot(x,y,type="l",main=paste0("a=",a,",b=",b))
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#####  VISUALIZATION FUNCTIONS :
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


##----------------------------------------------------------------------------##
## Function to extract fitted etiology, i.e. \hat{\pi}_k(x).
## Note that, when Icat is extracted, we determine whether to include all people
## or only a subset without realGS (known infection status) according to 
## whether ParameterList has bcxList or not. When bcxList is absent, we always
## extract those people without realGS, even though model could be fitted 
## without using realGS data and hence predicting those people's infection status.
##----------------------------------------------------------------------------##

bugs_extraction = function(fullname,VarName,ParameterList){
  simsMat = read.coda(paste(fullname,"coda1.txt",sep="\\"),
                      paste(fullname,"codaIndex.txt",sep="\\"),quiet=TRUE)
  Jallowed = length(ParameterList$AllowedList)
  
  if ( VarName =="pgridPred"){
    grep(VarName,colnames(simsMat))
    Ngrid = length(ParameterList$ERcovPred)
    res = array(NA,c(nrow(simsMat),Jallowed,Ngrid))
    for (k in 1:Ngrid){
      SubVarName = rep(NA,Jallowed)
      for (j in 1:Jallowed){
        SubVarName[j] = paste(VarName,"[",k,",",j,"]",sep="")
      }
      res[,,k] = simsMat[,SubVarName]
    } 
    return(res)
  } else if(VarName=="pEti" ){
    grep(VarName,colnames(simsMat))
    res = matrix(NA,nrow=nrow(simsMat),ncol=Jallowed)
    SubVarName = rep(NA,Jallowed)
    for (j in 1:Jallowed){
      SubVarName[j] = paste(VarName,"[",j,"]",sep="")
    }
    res = simsMat[,SubVarName]
    return(res)
  } else if(VarName=="Icat"){
    grep(VarName,colnames(simsMat))
    if (!is.null(ParameterList$bcxList)){
      Npred = ParameterList$Nd
      resnm = 1:Nd
    } else {
      Npred = ParameterList$Nd-floor(ParameterList$GSfraction*ParameterList$Nd)
      resnm = tail(1:Nd,n=Npred)
    }
    res = array(NA,c(nrow(simsMat),Jallowed,Npred))
    
    multinom.template = rbind(diag(Jallowed),rep(0,Jallowed))
    for (k in 1:Npred){
      if (!is.null(ParameterList$bcxList)){
        kstar = k
      } else{
        kstar = k+floor(ParameterList$GSfraction*ParameterList$Nd)
      }
      SubVarName = paste(VarName,"[",kstar,"]",sep="")
      SubCat = simsMat[,SubVarName]
      res[,,k] = multinom.template[SubCat,]
    } 
    return(list(res=res,resnm=resnm))  
  } 
}

##----------------------------------------------------------------------------##
## visualize population pie at 2d, so there is no covariate involved
## (currently BS+GS) with three optional fit. Default only use nppcr data.
##----------------------------------------------------------------------------##

visualize_popPie_2d = function(ParameterList,directoryName,
                               nppcr.no.realGS.fit = TRUE,
                               nppcr.realGS.fit    = FALSE,
                               realGS.no.nppcr.fit = FALSE,
                               no.realGS.no.nppcr.fit =FALSE,
                               nppcr.bcx.fit       =  FALSE,
                               nppcr.IS.fit        =  FALSE,
                               pEtiName="pEti",ksFrac=1,
                               mcex=2,pcex = 2,
                               pcol = c("dodgerblue2","lightseagreen",
                                        "gold","gray","orange","brown"),
                               levellabel=c(5)){
  ##Ngrid: is the length of grid over which fitted etiology \hat{\pi}_k(x)
  ##will be plotted
  ##directoryName: the folder name that stores coda.txt file
  
  require(MASS)
  require(robCompositions)
  require(ks)
  require(compositions)
  
  Jallowed = length(ParameterList$AllowedList)
  truth = ParameterList$etiology
  AllowedList = ParameterList$AllowedList
  
  aug.eti.samp = rbind(truth)
  col.grp = "red"#c("tan4")
  pch.grp = c(10)
  cex.grp = c(5)
  
  #ternaryplot(aug.eti.samp,col=col.grp,pch=pch.grp)
  ternaryDiag(
    aug.eti.samp,
    pch = pch.grp,
    col = col.grp,
    cex=cex.grp,name=rep(NA,3),
    lwd=c(5),
    main = "Posterior Distribtion of Population Etiology"
  )
  bool =  which(c(nppcr.no.realGS.fit,nppcr.realGS.fit,realGS.no.nppcr.fit,no.realGS.no.nppcr.fit,
                  nppcr.bcx.fit,nppcr.IS.fit)==TRUE)
  legend("topright",
         c("nppcr.only.fit","nppcr.realGS.fit","realGS.no.nppcr.fit","no.realGS.no.nppcr.fit",
           "nppcr.bcx.fit","nppcr.IS.fit")[bool],
         col=pcol[bool],
         lty=c(2,1,2,2,1,1)[bool],
         lwd=3)
  
  name=AllowedList
  mtext(name[1], side = 1, line = -1, at = -0.1, cex = mcex)
  mtext(name[2], side = 1, line = -1, at = 1.1, cex = mcex)
  mtext(name[3], side = 3,  line=-1, at=0.5 ,cex = mcex)
  
  ##########################################
  ## add plots for different sources of data
  ##########################################
  if (nppcr.no.realGS.fit){
    pEti = bugs_extraction(paste0(directoryName,"\\nppcr.no.realGS.fit"),
                           pEtiName,ParameterList)
    postMean = colMeans(pEti)
    
    temp = compositions::ilr(pEti)[1:floor(ksFrac*nrow(pEti)),]
    
    coord1 = seq(0.001,0.999,by=0.01)
    coord2 = sqrt(0.75)*coord1
    
    coord.grid = expand.grid(coord1=coord1,coord2=coord2)
    prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                                 y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                                 x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                                   1/sqrt(3)*coord.grid$coord2)
    prob.grid = rev(prob.grid.temp)
    
    ilr.grid    = compositions::ilr(prob.grid)
    
    fhat1       = kde(x=temp,H=Hpi(temp),compute.cont=T)
    fhat2       = kde(x=temp,H=Hpi(temp),eval.points=ilr.grid)
    
    tri.z       = matrix(fhat2$estimate,nr=length(coord1),nc=length(coord2))
    
    f = function(x,y){
      (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
    }
    
    for (i in 1:length(coord1)){
      for (j in 1:length(coord2)){
        tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA) 
      }  
    }
    
    
    num.levels      = levellabel#c(5,10,20,50)
    cont.levels     = paste(num.levels,"%",sep="")
    label.levels     = paste(100-num.levels,"%",sep="")
    
    points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
           cex=pcex,col=pcol[1],pch=20)
    contour(coord1,coord2,tri.z,
            levels =(fhat1$cont)[cont.levels],
            labels = label.levels,add=TRUE,lwd=3,col=c(pcol[1]),lty=2)
  }
  
  
  if (nppcr.realGS.fit){
    pEti = bugs_extraction(paste0(directoryName,"\\nppcr.realGS.fit"),
                           pEtiName,ParameterList)
    postMean = colMeans(pEti)
    
    temp = compositions::ilr(pEti)[1:floor(ksFrac*nrow(pEti)),]
    
    coord1 = seq(0.001,0.999,by=0.01)
    coord2 = sqrt(0.75)*coord1
    
    coord.grid = expand.grid(coord1=coord1,coord2=coord2)
    prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                                 y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                                 x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                                   1/sqrt(3)*coord.grid$coord2)
    prob.grid = rev(prob.grid.temp)
    
    ilr.grid    = compositions::ilr(prob.grid)
    
    fhat1       = kde(x=temp,H=Hpi(temp),compute.cont=T)
    fhat2       = kde(x=temp,H=Hpi(temp),eval.points=ilr.grid)
    
    tri.z       = matrix(fhat2$estimate,nr=length(coord1),nc=length(coord2))
    
    f = function(x,y){
      (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
    }
    
    for (i in 1:length(coord1)){
      for (j in 1:length(coord2)){
        tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA) 
      }  
    }
    
    
    num.levels      = levellabel#c(5,10,20,50)
    cont.levels     = paste(num.levels,"%",sep="")
    label.levels     = paste(100-num.levels,"%",sep="")
    
    points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
           cex=pcex,col=pcol[2],pch=20)
    contour(coord1,coord2,tri.z,
            levels =(fhat1$cont)[cont.levels],
            labels = label.levels,add=TRUE,lwd=3,col=c(pcol[2]),lty=1)
  }
  
  
  if (realGS.no.nppcr.fit){
    pEti = bugs_extraction(paste0(directoryName,"\\realGS.no.nppcr.fit"),
                           pEtiName,ParameterList)
    postMean = colMeans(pEti)
    
    temp = compositions::ilr(pEti)[1:floor(ksFrac*nrow(pEti)),]
    
    coord1 = seq(0.001,0.999,by=0.01)
    coord2 = sqrt(0.75)*coord1
    
    coord.grid = expand.grid(coord1=coord1,coord2=coord2)
    prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                                 y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                                 x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                                   1/sqrt(3)*coord.grid$coord2)
    prob.grid = rev(prob.grid.temp)
    
    ilr.grid    = compositions::ilr(prob.grid)
    
    fhat1       = kde(x=temp,H=Hpi(temp),compute.cont=T)
    fhat2       = kde(x=temp,H=Hpi(temp),eval.points=ilr.grid)
    
    tri.z       = matrix(fhat2$estimate,nr=length(coord1),nc=length(coord2))
    
    f = function(x,y){
      (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
    }
    
    for (i in 1:length(coord1)){
      for (j in 1:length(coord2)){
        tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA) 
      }  
    }
    
    
    num.levels      = levellabel#c(5,10,20,50)
    cont.levels     = paste(num.levels,"%",sep="")
    label.levels     = paste(100-num.levels,"%",sep="")
    
    points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
           cex=pcex,col=pcol[3],pch=20)
    contour(coord1,coord2,tri.z,
            levels =(fhat1$cont)[cont.levels],
            labels = label.levels,add=TRUE,lwd=3,col=c(pcol[3]),lty=2)
  }
  
  if (no.realGS.no.nppcr.fit){
    pEti = bugs_extraction(paste0(directoryName,"\\no.realGS.no.nppcr.fit"),
                           pEtiName,ParameterList)
    postMean = colMeans(pEti)
    
    temp = compositions::ilr(pEti)[1:floor(ksFrac*nrow(pEti)),]
    
    coord1 = seq(0.001,0.999,by=0.01)
    coord2 = sqrt(0.75)*coord1
    
    coord.grid = expand.grid(coord1=coord1,coord2=coord2)
    prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                                 y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                                 x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                                   1/sqrt(3)*coord.grid$coord2)
    prob.grid = rev(prob.grid.temp)
    
    ilr.grid    = compositions::ilr(prob.grid)
    
    fhat1       = kde(x=temp,H=Hpi(temp),compute.cont=T)
    fhat2       = kde(x=temp,H=Hpi(temp),eval.points=ilr.grid)
    
    tri.z       = matrix(fhat2$estimate,nr=length(coord1),nc=length(coord2))
    
    f = function(x,y){
      (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
    }
    
    for (i in 1:length(coord1)){
      for (j in 1:length(coord2)){
        tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA) 
      }  
    }
    
    
    num.levels      = levellabel#c(5,10,20,50)
    cont.levels     = paste(num.levels,"%",sep="")
    label.levels     = paste(100-num.levels,"%",sep="")
    
    points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
           cex=pcex,col=pcol[4],pch=20)
    contour(coord1,coord2,tri.z,
            levels =(fhat1$cont)[cont.levels],
            labels = label.levels,add=TRUE,lwd=3,col=c(pcol[4]),lty=2)
  }
  
  if (nppcr.bcx.fit){
    pEti = bugs_extraction(paste0(directoryName,"\\nppcr.bcx.fit"),
                           pEtiName,ParameterList)
    postMean = colMeans(pEti)
    
    temp = compositions::ilr(pEti)[1:floor(ksFrac*nrow(pEti)),]
    
    coord1 = seq(0.001,0.999,by=0.01)
    coord2 = sqrt(0.75)*coord1
    
    coord.grid = expand.grid(coord1=coord1,coord2=coord2)
    prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                                 y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                                 x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                                   1/sqrt(3)*coord.grid$coord2)
    prob.grid = rev(prob.grid.temp)
    
    ilr.grid    = compositions::ilr(prob.grid)
    
    fhat1       = kde(x=temp,H=Hpi(temp),compute.cont=T)
    fhat2       = kde(x=temp,H=Hpi(temp),eval.points=ilr.grid)
    
    tri.z       = matrix(fhat2$estimate,nr=length(coord1),nc=length(coord2))
    
    f = function(x,y){
      (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
    }
    
    for (i in 1:length(coord1)){
      for (j in 1:length(coord2)){
        tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA) 
      }  
    }
    
    
    num.levels      = levellabel#c(5,10,20,50)
    cont.levels     = paste(num.levels,"%",sep="")
    label.levels     = paste(100-num.levels,"%",sep="")
    
    points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
           cex=pcex,col=pcol[5],pch=20)
    contour(coord1,coord2,tri.z,
            levels =(fhat1$cont)[cont.levels],
            labels = label.levels,add=TRUE,lwd=3,col=c(pcol[5]),lty=1)
  }
  
  
  if (nppcr.IS.fit){
    pEti = bugs_extraction(paste0(directoryName,"\\nppcr.IS.fit"),
                           pEtiName,ParameterList)
    postMean = colMeans(pEti)
    
    temp = compositions::ilr(pEti)[1:floor(ksFrac*nrow(pEti)),]
    
    coord1 = seq(0.001,0.999,by=0.01)
    coord2 = sqrt(0.75)*coord1
    
    coord.grid = expand.grid(coord1=coord1,coord2=coord2)
    prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                                 y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                                 x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                                   1/sqrt(3)*coord.grid$coord2)
    prob.grid = rev(prob.grid.temp)
    
    ilr.grid    = compositions::ilr(prob.grid)
    
    fhat1       = kde(x=temp,H=Hpi(temp),compute.cont=T)
    fhat2       = kde(x=temp,H=Hpi(temp),eval.points=ilr.grid)
    
    tri.z       = matrix(fhat2$estimate,nr=length(coord1),nc=length(coord2))
    
    f = function(x,y){
      (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
    }
    
    for (i in 1:length(coord1)){
      for (j in 1:length(coord2)){
        tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA) 
      }  
    }
    
    
    num.levels      = levellabel#c(5,10,20,50)
    cont.levels     = paste(num.levels,"%",sep="")
    label.levels     = paste(100-num.levels,"%",sep="")
    
    points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
           cex=pcex,col=pcol[6],pch=20)
    contour(coord1,coord2,tri.z,
            levels =(fhat1$cont)[cont.levels],
            labels = label.levels,add=TRUE,lwd=3,col=c(pcol[6]),lty=1)
  }
  
  
  #res = list(colMeans(pEti),matrix(truth,nrow=1))
  #names(res) = c("fitted mean", "true etiology")
  #res
}

#visualize_popPie_2d(ParameterList,fullname,T,T,T,ksFrac=1)

##----------------------------------------------------------------------------##
## visualize individual pie at 2d, so there is no covariate involved
## (currently BS+GS). Current implementation have nppcr +bcx prediction (black solid line)
## andnppcr only prediction (blue dotted line). The actual individual data
## are plotted at the top of the figure. 
##----------------------------------------------------------------------------##


visualize_indPie_2d = function(ParameterList, individual.index,
                               superdat = data,
                               directoryName=fullname,
                               nppcr.bcx.fit=TRUE,
                               nppcr.no.realGS.fit = TRUE,
                               nppcr.bcx.fit.dir = "\\nppcr.bcx.fit",
                               nppcr.no.realGS.fit.dir = "\\nppcr.no.realGS.fit",
                               IcatName="Icat"){
  #make sure that bacteria subset are in the first several slots  
  if (nppcr.bcx.fit){
    #read MCMC outputs from result.folder
    tmp = bugs_extraction(paste0(directoryName,nppcr.bcx.fit.dir),
                                 IcatName,ParameterList)
    iL_allowed = tmp$res
    allowed.individual.index = tmp$resnm
    if (!all(individual.index %in% allowed.individual.index)){
      cat("Your input individual index does not need his/her etiology prediction!
          Please change an individual index.")
      return(NULL)
    }
    Jallowed = length(ParameterList$AllowedList)
    Nd       = ParameterList$Nd
    #NdGS     = floor(ParameterList$GSfraction*ParameterList$Nd)
    indEtiMCMC = t(apply(iL_allowed,c(2,3),mean))
    colnames(indEtiMCMC)=paste("indpie_",ParameterList$AllowedList)
    indpie = cbind(superdat[1:Nd,],indEtiMCMC)
    
    #set up grid of all possible measurement patterns
    tmpBS = list()
    for (l in 1:length(ParameterList$FullPathogenName)){tmpBS[[l]] = c(1,0)}
    BSpatternGrid = as.matrix(expand.grid(tmpBS))
    colnames(BSpatternGrid) = paste("BS",ParameterList$FullPathogenName,sep="_")
    NBGrid = nrow(BSpatternGrid)
    
    tmpGS = list()
    for (l in 1:length(ParameterList$bcxList)){tmpGS[[l]] = c(1,0)}
    GSpatternGrid = as.matrix(expand.grid(tmpGS))
    colnames(GSpatternGrid) = paste("GS",ParameterList$bcxList,sep="_")
    NGGrid = nrow(GSpatternGrid)
    
    MeasGrid = list()
    count = 0
    for (i in 1:NBGrid){
      for (j in 1:NGGrid){
        count = count+1
        MeasGrid[[count]] = list(BS=BSpatternGrid[i,],GS=GSpatternGrid[j,]) 
      }
    }
    
    topprob=1
    
    for (i in individual.index){
      #o<-par(mar=c(7,4,2,2),mfrow=c(1,2))
      #B+G
      ##dev.off()
      op <- par(mar = c(0,4.1,3,2.1))
      layout(matrix(c(1,2),2,1),widths = c(8),heights = c(2,6))
      plot(1:Jallowed,rep(0,Jallowed),ylim=c(0,topprob),pch="",xaxt="n",yaxt="n",
           ylab="",main=paste0(i,"th case individual pie"))
      axis(2,at = topprob*c(0.3,0.8),c("GS","BS"),las=2)
      
      for (s in 1:Jallowed){
        colindex = grep(paste0("MBS_",ParameterList$AllowedList[s]),colnames(indpie))
        text(s,topprob*0.8,indpie[i,colindex],
             col="dodgerblue2")
      }
      for (s in 1:length(ParameterList$bcxList)){
        colindex = grep(paste0("MGS_",ParameterList$AllowedList[s]),colnames(indpie))
        text(s,topprob*0.3,indpie[i,colindex],
             col="gold")
      }
      par(op)
      par(mar = c(5.1,4.1,0.5,2.1))
      plot(1:Jallowed,indEtiMCMC[i,],
           xlab="pathogens",ylab="probability",xaxt="n",
           type="h",ylim=c(0,topprob),lwd=5)
      axis(1,at = 1:Jallowed, labels=ParameterList$AllowedList,las=2)
      
      if (nppcr.no.realGS.fit){
        iL_allowed.nppcr.only = bugs_extraction(paste0(directoryName,nppcr.no.realGS.fit.dir),
                                                IcatName,ParameterList)$res
        indEtiMCMC.nppcr.only = t(apply(iL_allowed.nppcr.only,c(2,3),mean))
        points((1:Jallowed)+0.1,indEtiMCMC.nppcr.only[i,],type="h",col="dodgerblue2",
               lwd=5,lty=1)
      } 
    }
  } else {
    cat("nppcr.bcx.fit has been set to FALSE. Please check whether the nppcr+bcx model is
        fitted.")
    return(NULL)
  } 
}


#visualize_indPie_2d(ParameterList)

##---------------------------------------------------------------------------##
## Infograph visualization for population pie
##---------------------------------------------------------------------------##
infograph_popPie = function(directoryName,
                            alphaE=alpha,
                            alphaB = rep(1,Jfull),
                            betaB = rep(1,Jfull),
                            alphaG = rep(1,JGS),
                            betaG  = rep(1,JGS),
                            Bdat=MBS,Gdat=MGS,
                            ncase = Nd,J.GS = JGS,
                            nppcr.no.realGS.fit=FALSE,
                            nppcr.no.realGS.fit.dir=prior.folder,#"\\nppcr.no.realGS.fit",
                            nppcr.bcx.fit=TRUE,
                            nppcr.bcx.fit.dir=prior.folder,#"\\nppcr.bcx.fit",
                            pEtiName="pEti",pathogens = Pathogen){
  
#   directoryName= fullname
#   alphaE = alpha
#   alphaB = alphaB
#   betaB = betaB
#   alphaG = alphaG
#   betaG  = betaG
#   nppcr.bcx.fit=TRUE
#   nppcr.no.realGS.fit = FALSE
#   pathogens=Pathogen
#   ncase = Nd
#   J.GS=JGS
#   pEtiName="pEti"
#   nppcr.bcx.fit.dir=prior.folder
#   Bdat=MBS;Gdat=MGS
  
  
  require(R2WinBUGS)
  #directoryName = fullname
  cexval=1
  srtval=0
  J.full = length(pathogens)
  if (nppcr.bcx.fit){
    layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), 
           widths=c(2,2,4), heights=c(8))
    #layout.show(n=3)
    curr.folder = paste0(directoryName,nppcr.bcx.fit.dir)
    simsMat = read.coda(paste(curr.folder,"coda1.txt",sep="\\"),
                        paste(curr.folder,"codaIndex.txt",sep="\\"),quiet=TRUE)  
    #grep(pEtiName,colnames(simsMat))
    res = matrix(NA,nrow=nrow(simsMat),ncol=J.full)
    SubVarName = rep(NA,J.full)
    for (j in 1:J.full){
      SubVarName[j] = paste(pEtiName,"[",j,"]",sep="")
    }
    pmatBG = simsMat[,SubVarName]
    pmeanBG = colMeans(pmatBG)
    
    pq1BG   = apply(pmatBG,2,quantile,probs=0.025)
    pq2BG   = apply(pmatBG,2,quantile,probs=0.975)
    p.innerq1BG   = apply(pmatBG,2,quantile,probs=0.25)
    p.innerq2BG   = apply(pmatBG,2,quantile,probs=0.75)

    ######################################################
    ##subplot 1: case-control rate comparisons
    ######################################################
    require(binom)
    nctrl = nrow(Bdat)-ncase
    tmp.case = binom.confint(colSums(Bdat[1:ncase,]), ncase, conf.level = 0.95, methods = "ac")
    tmp.ctrl = binom.confint(colSums(Bdat[-(1:ncase),]), nctrl, conf.level = 0.95, methods = "ac")
    Bcomp = rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5))
    Bcompq1 = rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")])
    Bcompq2 = rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")])
    
    thetamatBG = simsMat[,grep("thetaBS",colnames(simsMat))]
    thetameanBG = colMeans(thetamatBG)
    
    psimatBG = simsMat[,grep("psiBS",colnames(simsMat))]
    psimeanBG = colMeans(psimatBG)  
    
    fittedmeanBGcase = sapply(1:J.full, function(s) mean(pmatBG[,s]*thetamatBG[,s]+(1-pmatBG[,s])*psimatBG[,s]))
    fittedmeanBGcontrol = psimeanBG
    
    top2 = 1.3
    op<-par(mar=c(5.1,4.1,4.1,0))
    plotat = seq(0.5,J.full+0.5,by=1/4)[-(c(1,(1:J.full)*4+1))]
    plot(c(rbind(fittedmeanBGcase,Bcomp)),plotat,yaxt="n",xlim=c(0,top2),
         ylim=c(0.5,J.full+.5),xaxt="n",
         ylab="",xlab="probability",
         pch = c(rbind(rep(2,J.full),rep(20,J.full),rep(20,J.full))),
         col=c(rbind(rbind(rep(1,J.full),rep("dodgerblue2",J.full),rep("dodgerblue2",J.full)))),
         cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
    axis(2,at = plotat,labels=rep(c("","case","ctrl"),J.full),las=2)
    axis(1,at = c(0,0.2,0.4,0.6,0.8,1),labels= c(0,0.2,0.4,0.6,0.8,1),las=1)
    points(c(rbind(thetameanBG,Bcompq2)),plotat,
           pch=c(rbind(rep("+",J.full),rep("|",J.full),rep("|",J.full))),
           cex=c(rbind(rep(2,J.full),rep(1,J.full),rep(1,J.full))),
           col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
    points(c(rbind(psimeanBG,Bcompq1)),plotat,
           pch=c(rbind(rep("*",J.full),rep("|",J.full),rep("|",J.full))),
           cex=c(rbind(rep(2,J.full),rep(1,J.full),rep(1,J.full))),
           col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
    abline(h=seq(1.5,J.full-.5,by=1),lty=2,lwd=0.5,col="blue")
    abline(v=1,lty=2,lwd=.5)
    for (s in 1:J.full){
      if (any(Bcomp[,s]==0)){
        cat("Pathogen ",pathogens[s]," have zero cells cross-classifying Y and MBS.")
        text(top2*.9,s,NA)
      } else{
        
        tmp = LOR(MBS[,s],rep(c(1,0),times=c(ncase,nrow(MBS)-ncase)))
        L=round(tmp$OR[2],1)
        C=round(tmp$OR[1],1)
        R=round(tmp$OR[3],1)
        text(top2-0.1,s+1/(2*Jfull),C)
        text(top2-0.1,s-.2,paste(c(L,"   ",R),collapse=" "),cex=0.8)
      }
    }
    
    counter = 0
    for (s in 1:(3*J.full)){
      segments(y0=plotat[s],x0=c(rbind(psimeanBG,Bcompq1))[s],
               y1=plotat[s],x1=c(rbind(thetameanBG,Bcompq2))[s],col="black",
               lty=ifelse((s-1)%%3<1,4,1))
      if ((s-1)%%3>=1){
        counter=counter+1
        tmp.hpos = ifelse(c(Bcompq2)[counter]+0.15>0.95,c(Bcompq1)[counter]-0.2,c(Bcompq2)[counter]+0.15 )
        text(tmp.hpos,plotat[s],paste0(round(100*c(Bcomp),1)[counter],"%"),srt=srtval,cex=cexval)
      }
    }
    
    for (s in 1:(J.full)){
      segments(y0=plotat[3*s-1],x0=c(rbind(fittedmeanBGcase,Bcomp))[3*s-1],
               y1=plotat[3*s],x1=c(rbind(fittedmeanBGcase,Bcomp))[3*s],col="black")
    }
    # put prior shapes on bronze sensitivity
    for (s in 1:J.full){
      tmp = rbeta(10000,alphaB[s],betaB[s])
      boxplot(tmp,at = s-0.45, boxwex=1/10 , col="gray",add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
      tmp.post = as.matrix(thetamatBG)[,s]
      boxplot(tmp.post,at = s-0.35,boxwex=1/10,add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
    }
    axis(2,at=(1:J.full)-.45,label=rep("",J.full),las=2,cex.axis=.5)
    axis(2,at=(1:J.full)-.35,label=rep("",J.full),las=2,cex.axis=.5)
    
    #mtext("o:""observed rate in case/contrl",line=3,col="dodgerblue2",cex=.8)     
    #legend("top",c("obs","fitted","TPR","FPR"),pch=c(20,2,3,8),horiz=TRUE,bty="n",cex=1.1,
    #       col=c("blue","black","black","black"))
    
#     mtext(text=c(expression(dot("")),":observed rate","+",":TPR","*",":FPR"),
#           adj=c(0,.1,0.6,0.7,0.85,1), 
#           col=c("dodgerblue2","black","black","black","black","black"), 
#           side=3, line=.5,cex=c(3,1,1,1,1,1))
    mtext(expression(underline("BrS")),line=1,cex=1.8)
    par(op)
    
    ######################################################
    ##subplot 2: gold-standard data comparisons
    ######################################################
    
    Gcomplete.index = which(rowSums(is.na(Gdat[,1:J.GS]))==0)
    cat("No. of Complete Gold-standard data on ",pathogens[1:J.GS],":",
        length(Gcomplete.index),"\n")
    cat(ncase-length(Gcomplete.index)," cases deleted for plotting third subgraph.","\n")
    Gdat = Gdat[Gcomplete.index,]
    tmpG.case = binom.confint(colSums(Gdat[,1:J.GS]), nrow(Gdat), conf.level = 0.95, methods = "ac")
    Gcomp = rbind(round(tmpG.case$mean,5),rep(NA,J.GS))
    Gcompq1 = rbind(tmpG.case[,c("lower")],rep(NA,J.GS))
    Gcompq2 = rbind(tmpG.case[,c("upper")],rep(NA,J.GS))
    
    thetamatG = simsMat[,grep("thetaGS",colnames(simsMat))]
    thetameanG = colMeans(thetamatG)
    
    thetamatGq1=apply(thetamatG,2,quantile,0.025)
    thetamatGq2=apply(thetamatG,2,quantile,0.975)
    
    fittedmeanGpos = sapply(1:J.GS, function(s) mean(pmatBG[,s]*thetamatG[,s]))
    
    
    top3 = .25
    par(mar=c(5.1,0,4.1,3))
    #op<-par(mar=c(5.1,0,4.1,2.1),xpd=TRUE)
    plotat = seq(0.5,J.full+0.5,by=1/4)[-(c(1,(1:J.full)*4+1))]
    plotat.short = plotat[1:length(c(rbind(thetameanG,Gcomp)))]
    plot(c(rbind(thetameanG,Gcomp)),plotat.short,yaxt="n",xlim=c(0,top3),
         ylim=c(0.5,J.full+.5),xaxt="n",
         ylab="",xlab="probability",
         pch = c(rbind(rep(20,J.full),rep(20,J.full),rep(20,J.full))),
         col=c(rbind(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full)))),
         cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
    axis(1,at = seq(0,0.25,len=6),labels= seq(0,0.25,len=6),las=1)
    points(c(rbind(fittedmeanGpos,matrix("",nrow=2,ncol=J.GS))),plotat.short,yaxt="n",xlim=c(0,top3),
         ylim=c(0.5,J.full+.5),xaxt="n",
         ylab="",#xlab="Gold Positive Rate",
         pch = c(rbind(rep(2,J.full),rep(NA,J.full),rep(NA,J.full))),
         col=c(rbind(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full)))),
         cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
    
    #axis(2,at = plotat.short,labels=rep(c("TPR-G","case","ctrl"),J.GS),las=2)
    points(c(rbind(thetamatGq2,Gcompq2)),plotat.short,
           pch=c(rbind(rep("|",J.full),rep("|",J.full),rep("|",J.full))),
           cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
           col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
    points(c(rbind(thetamatGq1,Gcompq1)),plotat.short,
           pch=c(rbind(rep("|",J.full),rep("|",J.full),rep("|",J.full))),
           cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
           col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
    
    #inner 25%-75%
    thetamatGq.inner1=apply(thetamatG,2,quantile,0.25)
    thetamatGq.inner2=apply(thetamatG,2,quantile,0.75)
    points(c(rbind(thetamatGq.inner1,Gcompq1)),plotat.short,
           pch=c(rbind(rep("[",J.full),rep("|",J.full),rep("|",J.full))),
           cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
           col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
    points(c(rbind(thetamatGq.inner2,Gcompq1)),plotat.short,
           pch=c(rbind(rep("]",J.full),rep("|",J.full),rep("|",J.full))),
           cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
           col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
    counter = 0
    for (s in 1:length(plotat.short)){
     segments(y0=plotat.short[s],x0=c(rbind(thetamatGq.inner1,Gcompq1))[s],
              y1=plotat.short[s],x1=c(rbind(thetamatGq.inner2,Gcompq2))[s],
              col="black",
              lwd=2)
    }
    
    
    abline(h=seq(1.5,J.full-.5,by=1)[1:J.GS],lty=2,lwd=0.5,col="blue")
    
    
    counter = 0
    for (s in 1:length(plotat.short)){
      segments(y0=plotat.short[s],x0=c(rbind(thetamatGq1,Gcompq1))[s],
               y1=plotat.short[s],x1=c(rbind(thetamatGq2,Gcompq2))[s],col="black",
               lty=ifelse((s-1)%%3<2,1,1))
      if ((s-1)%%3>=1){
        counter=counter+1
        text(c(Gcomp)[counter]+0.03,plotat[s],paste0(round(100*c(Gcomp),1)[counter],"%"),srt=srtval,cex=cexval)
      }
    }
    
    
     for (s in 1:J.GS){
#       segments(y0=plotat.short[3*s-1],x0=c(rbind(thetameanG,Gcomp))[3*s-1],
#                y1=plotat.short[3*s],x1=c(rbind(thetameanG,Gcomp))[3*s],col="black")
       #text on estimated gold sensitivity
       text(thetameanG[s],plotat.short[3*s-2]+.125,paste(round(100*thetameanG[s],2),"%"))
     }
    
    # put prior shapes on gold sensitivity
    for (s in 1:J.GS){
      tmp = rbeta(10000,alphaG[s],betaG[s])
      boxplot(tmp,at = s-0.45, boxwex=1/8 ,col="gray", add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
      
    }
    #axis(2,at=(1:J.GS)-.45,label=rep("prior",J.GS),las=2)
    
    #mtext(text=c(expression(dot("")),":observed rate",
    #             expression(symbol(Delta)),":fitted pos(%)"),
    #      adj=c(0,.1,0.6,1), 
    #      col=c("gold","black","green","black"), 
    #      side=3, line=.5,cex=c(3,1,1,1))
    #legend(.1,5.5,c("observed rate","fitted rate"),pch=c(16,2),cex=1.2,bty="n")
    mtext(expression(underline("GnS")),line=1,cex=1.8)
    
    par(op)
    
    
    ######################################################
    ##subplot 3: etiology bars with prior shapes
    ######################################################
    top=0.4#max(pq2BG)
    dotcolor = "black"
    op <- par(mar=c(5.1,6,4.1,1.1))
    plot(c(pmeanBG),1:(J.full),
         yaxt="n",#xaxt="n",
         xlim=c(0,top),ylim=c(0.5,J.full+0.5),col=c("black"),
         ylab="",xlab="probability",
         pch=c(20),cex=2)
    axis(2,at=1:J.full,labels=pathogens,las=2,cex.axis=1.5)
    abline(h=seq(1.5,J.full-.5,by=1),lty=2,lwd=0.5,col="blue")
    #draw axis within plot:
    for (s in 1:(J.full-1)){
      axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=1,#labels=rep("",length(seq(0,1,by=.2))),
           pos = seq(1.5,J.full-.5,by=1)[s], cex.axis = 0.8,lty=2,col="blue")
      # axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=0,#labels=rep("",length(seq(0,1,by=.2))),
      #      pos = seq(1.5,J.full-.5,by=1)[s]+0.3, cex.axis = 0.8,lty=2,col="blue")
    }
    points(c(pq1BG),1:(J.full),pch="|",cex=1)
    points(c(pq2BG),1:(J.full),pch="|",cex=1)
    points(c(p.innerq1BG),1:(J.full),pch="[",cex=1)
    points(c(p.innerq2BG),1:(J.full),pch="]",cex=1)
    
    mtext(expression(underline(hat(pi))),line=1,cex=1.8)
    #mtext(c(expression(bold("--")),":prior","-",":posterior"),col=c("gray","black","black","black"),
    #      adj=c(0,0.1,0.3,0.4),line=.8,cex=.8,lwd=2)
    legend("topright",c("prior","posterior"),lty=c(2,1),col=c("gray","black"),
           lwd = 4,horiz=TRUE,cex=1.5,bty="n")
    pgrid = seq(0,1,by=0.01)
    
    for (s in 1:(J.full)){
      segments(y0=s,x0=c(pq1BG)[s],y1=s,x1=c(pq2BG)[s],col=dotcolor)
      segments(y0=s,x0=c(p.innerq1BG)[s],y1=s,x1=c(p.innerq2BG)[s],col=dotcolor,lwd=2)
      #text(pmeanBG[s],s+0.25,paste0(round(100*c(pmeanBG),2)[s],"%"),srt=srtval,cex=cexval)
      text(.35,s+0.25,paste0("=",paste0(round(100*c(pmeanBG),1)[s],"%")),srt=srtval,cex=2)
      text(.3,s+0.25,bquote(hat(pi)[.(s)]),srt=srtval,cex=2)
      #text(top-0.05,s,pathogens[s],cex=.8,srt=srtval)
      tmp.density = dbeta(pgrid,alphaE[s],sum(alphaE[-s]))
      points(pgrid,tmp.density/(3*max(tmp.density))+s-0.45,type="l",col="gray",lwd=4,lty=2)
      ##posterior density
      tmp.post.density = density(pmatBG[,s],from=0,to=1)
      tmp.x = tmp.post.density$x
      tmp.y = tmp.post.density$y
      points(tmp.x,tmp.y/(3*max(tmp.y))+s-0.45,col="black",lwd=4,type="l")
      
    }
    
    par(op)
    
    
    #######################################################
    #######################################################
    #######################################################
    #######################################################
    #######################################################
  } 
  if (nppcr.no.realGS.fit){  
    layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), 
           widths=c(3,2,2), heights=c(8))
    #layout.show(n=3)
    curr.folder = paste0(directoryName,nppcr.no.realGS.fit.dir)
    simsMat = read.coda(paste(curr.folder,"coda1.txt",sep="\\"),
                        paste(curr.folder,"codaIndex.txt",sep="\\"),quiet=TRUE)  
    #grep(pEtiName,colnames(simsMat))
    res = matrix(NA,nrow=nrow(simsMat),ncol=J.full)
    SubVarName = rep(NA,J.full)
    for (j in 1:J.full){
      SubVarName[j] = paste(pEtiName,"[",j,"]",sep="")
    }
    pmatBG = simsMat[,SubVarName]
    pmeanBG = colMeans(pmatBG)
    pq1BG   = apply(pmatBG,2,quantile,probs=0.025)
    pq2BG   = apply(pmatBG,2,quantile,probs=0.975)
    p.innerq1BG   = apply(pmatBG,2,quantile,probs=0.25)
    p.innerq2BG   = apply(pmatBG,2,quantile,probs=0.75)
    ######################################################
    ##subplot 1: etiology bars with prior shapes
    ######################################################
    top=1#max(pq2BG)
    dotcolor = "black"
    op <- par(mar=c(5.1,6,4.1,0))
    plot(c(pmeanBG),1:(J.full),
         yaxt="n",xaxt="n",
         xlim=c(0,top),ylim=c(0.5,J.full+0.5),col=c("black"),
         ylab="",xlab="probability",
         pch=c(20),cex=2)
    axis(2,at=1:J.full,labels=pathogens,las=2)
    abline(h=seq(1.5,J.full-.5,by=1),lty=2,lwd=0.5,col="blue")
    #draw axis within plot:
    for (s in 1:(J.full-1)){
      axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=1,#labels=rep("",length(seq(0,1,by=.2))),
           pos = seq(1.5,J.full-.5,by=1)[s], cex.axis = 0.8,lty=2,col="blue")
      #axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=0,#labels=rep("",length(seq(0,1,by=.2))),
      #     pos = seq(1.5,J.full-.5,by=1)[s]+0.3, cex.axis = 0.8,lty=2,col="blue")
    }
    points(c(pq1BG),1:(J.full),pch="|",cex=1)
    points(c(pq2BG),1:(J.full),pch="|",cex=1)
    points(c(p.innerq1BG),1:(J.full),pch="[",cex=1)
    points(c(p.innerq2BG),1:(J.full),pch="]",cex=1)
    
    mtext(expression(underline("Population Etiology")),line=2)
    mtext(c("--",":prior","-",":posterior"),col=c("gray","black","black","black"),
          adj=c(0,0.1,0.3,0.4),line=.8,cex=.8)
    pgrid = seq(0,1,by=0.01)
    
    for (s in 1:(J.full)){
      segments(y0=s,x0=c(pq1BG)[s],y1=s,x1=c(pq2BG)[s],col=dotcolor)
      segments(y0=s,x0=c(p.innerq1BG)[s],y1=s,x1=c(p.innerq2BG)[s],col=dotcolor,lwd=2)
      text(pmeanBG[s],s+0.25,paste0(round(100*c(pmeanBG),1)[s],"%"),srt=srtval,cex=cexval)
      #text(top-0.05,s,pathogens[s],cex=.8,srt=srtval)
      tmp.density = dbeta(pgrid,alphaE[s],sum(alphaE[-s]))
      points(pgrid,tmp.density/(3*max(tmp.density))+s-0.45,type="l",col="gray",lwd=2,lty=2)
      ##posterior density
      tmp.post.density = density(pmatBG[,s],from=0,to=1)
      tmp.x = tmp.post.density$x
      tmp.y = tmp.post.density$y
      points(tmp.x,tmp.y/(3*max(tmp.y))+s-0.45,col="black",lwd=2,type="l")
      
    }
    
    par(op)
    
    ######################################################
    ##subplot 2: case-control rate comparisons
    ######################################################
    require(binom)
    nctrl = nrow(Bdat)-ncase
    tmp.case = binom.confint(colSums(Bdat[1:ncase,]), ncase, conf.level = 0.95, methods = "ac",)
    tmp.ctrl = binom.confint(colSums(Bdat[-(1:ncase),]), nctrl, conf.level = 0.95, methods = "ac",)
    Bcomp = rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5))
    Bcompq1 = rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")])
    Bcompq2 = rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")])
    
    thetamatBG = simsMat[,grep("thetaBS",colnames(simsMat))]
    thetameanBG = colMeans(thetamatBG)
    
    psimatBG = simsMat[,grep("psiBS",colnames(simsMat))]
    psimeanBG = colMeans(psimatBG)  
    
    fittedmeanBGcase = sapply(1:J.full, function(s) mean(pmatBG[,s]*thetamatBG[,s]+(1-pmatBG[,s])*psimatBG[,s]))
    fittedmeanBGcontrol = psimeanBG
    
    top2 = 1.3
    op<-par(mar=c(5.1,4.1,4.1,0))
    plotat = seq(0.5,J.full+0.5,by=1/4)[-(c(1,(1:J.full)*4+1))]
    plot(c(rbind(fittedmeanBGcase,Bcomp)),plotat,yaxt="n",xlim=c(0,top2),
         ylim=c(0.5,J.full+.5),
         ylab="",xlab="Bronze Positive Rate",
         pch = c(rbind(rep(2,J.full),rep(20,J.full),rep(20,J.full))),
         col=c(rbind(rbind(rep(3,J.full),rep("dodgerblue2",J.full),rep("dodgerblue2",J.full)))),
         cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
    axis(2,at = plotat,labels=rep(c("B","case","ctrl"),J.full),las=2)
    points(c(rbind(thetameanBG,Bcompq2)),plotat,
           pch=c(rbind(rep("+",J.full),rep("|",J.full),rep("|",J.full))),
           cex=c(rbind(rep(2,J.full),rep(1,J.full),rep(1,J.full))),
           col=c(rbind(rep(3,J.full),rep(1,J.full),rep(1,J.full))))
    points(c(rbind(psimeanBG,Bcompq1)),plotat,
           pch=c(rbind(rep("*",J.full),rep("|",J.full),rep("|",J.full))),
           cex=c(rbind(rep(2,J.full),rep(1,J.full),rep(1,J.full))),
           col=c(rbind(rep(3,J.full),rep(1,J.full),rep(1,J.full))))
    abline(h=seq(1.5,J.full-.5,by=1),lty=2,lwd=0.5,col="blue")
    abline(v=1,lty=2,lwd=.5)
    for (s in 1:J.full){
      if (any(Bcomp[,s]==0)){
        cat("Pathogen ",pathogens[s]," have zero cells cross-classifying Y and MBS.")
        text(top2*.9,s,NA)
      } else{
        tmp = LOR(MBS[,s],rep(c(1,0),times=c(ncase,nrow(MBS)-ncase)))
        text(top2-0.15,s+1/(2*Jfull),paste0("OR:",round(tmp$OR[1],1)))
        text(top2-0.15,s-.35,paste0("pval:",round(tmp$pval,3)))
      }
    }
    
    counter = 0
    for (s in 1:(3*J.full)){
      segments(y0=plotat[s],x0=c(rbind(psimeanBG,Bcompq1))[s],
               y1=plotat[s],x1=c(rbind(thetameanBG,Bcompq2))[s],col="black",
               lty=ifelse((s-1)%%3<2,4,1))
      if ((s-1)%%3>=1){
        counter=counter+1
        tmp.hpos = ifelse(c(Bcompq2)[counter]+0.15>0.95,c(Bcompq1)[counter]-0.2,c(Bcompq2)[counter]+0.15 )
        text(tmp.hpos,plotat[s],paste0(round(100*c(Bcomp),2)[counter],"%"),srt=srtval,cex=cexval)
        
      }
    }
    
    for (s in 1:(J.full)){
      segments(y0=plotat[3*s-1],x0=c(rbind(fittedmeanBGcase,Bcomp))[3*s-1],
               y1=plotat[3*s],x1=c(rbind(fittedmeanBGcase,Bcomp))[3*s],col="black")
    }
    # put prior shapes on bronze sensitivity
    for (s in 1:J.full){
      tmp = rbeta(10000,alphaB[s],betaB[s])
      boxplot(tmp,at = s-0.45, boxwex=1/10 , col="gray",add=TRUE,horizontal=TRUE,outline=FALSE)
      tmp.post = as.matrix(thetamatBG)[,s]
      boxplot(tmp.post,at = s-0.35,boxwex=1/10,add=TRUE,horizontal=TRUE,outline=FALSE)
    }
    axis(2,at=(1:J.full)-.45,label=rep("TPR-B-prior",J.full),las=2,cex.axis=.5)
    axis(2,at=(1:J.full)-.35,label=rep("TPR-B-post",J.full),las=2,cex.axis=.5)
    
    #mtext("o:""observed rate in case/contrl",line=3,col="dodgerblue2",cex=.8)      
    mtext(text=c(expression(dot("")),":observed rate","+",":TPR","*",":FPR"),
          adj=c(0,.1,0.6,0.7,0.85,1), 
          col=c("dodgerblue2","black","green","black","green","black"), 
          side=3, line=.5,cex=c(3,1,1,1,1,1))
    mtext(expression(underline("Bronze-standard data")),line=2)
    par(op)
    
    
    ########################################################################
    ########################################################################
    ########################################################################
    ########################################################################
    ########################################################################
    ########################################################################
    
  }
  #else{
  #  cat("NPPCR+BCX data were not fitted together! Please make sure nppcr.bcx.fit
  #      is TRUE!")
  #  return(NULL)
  #}
}
#infograph_popPie(fullname)


##----------------------------------------------------------------------------##
## seasonal prevalence 
##----------------------------------------------------------------------------##
##smoothed:
seasonal_smoothed_prev = function(path,DATACASE = datcase.season,
                                  DATACTRL = datctrl.season){
  datcase.path = data.frame(DATACASE[,c(colnames(DATACASE)[pathind[path]],"NumDateCentered")])
  datctrl.path = data.frame(DATACTRL[,c(colnames(DATACTRL)[pathind[path]],"NumDateCentered")])
  colnames(datcase.path)[1] = "M"
  colnames(datctrl.path)[1] = "M"
  
  fitcase=gam(M~s(NumDateCentered,k=8),
              data=datcase.path,
              family=binomial(logit))
  
  EDAcase = predict(fitcase,type="response")[(1:nrow(datcase.path))]
  
  fitctrl=gam(M~s(NumDateCentered,k=8),
              data=datctrl.path,
              family=binomial(logit))
  
  EDActrl = predict(fitctrl,type="response")[(1:nrow(datctrl.path))]
  res = list(EDAcase,EDActrl)
  names(res) = c("EDAcase","EDActrl")
  return(res)
}

## discrete:
seasonal_prev = function(pathind,interval = 50,fig = FALSE,DATACASE = datacase,DATACTRL=datactrl){
  require(lubridate)
  require(binom)
  datcase.season = cbind(DATACASE$NPPCR,date= datacase$ENRLDATE)
  datctrl.season = cbind(DATACTRL$NPPCR,date= datactrl$ENRLDATE)
  
  date.case = as.Date(datcase.season$date, "%d-%B-%y")
  date.ctrl = as.Date(datctrl.season$date, "%d-%B-%y")
  
  
#   if (interval =="monthly"){
#     uniq.month.case=unique(paste(month(date.case),year(date.case),sep="-"))
#     uniq.month.ctrl=unique(paste(month(date.ctrl),year(date.ctrl),sep="-"))
#     if (!all(uniq.month.case==uniq.month.ctrl)){
#       cat("Cases and controls have different unique enrollment months!","\n")
#     } else{
#       uniq.month = uniq.month.case
#       m.y.grid = strsplit(uniq.month, split="-")
#       p.case = matrix(NA,nrow=length(Pathogen),ncol=length(uniq.month))
#       rownames(p.case) = Pathogen
#       colnames(p.case) = uniq.month
#       p.ctrl = p.case
#       
#       n.case = rep(NA,length=length(uniq.month))
#       names(n.case) = uniq.month
#       n.ctrl = n.case
#       
#       np.case = p.case
#       np.ctrl = p.case
#       up.case = p.case
#       lo.case = p.case
#       up.ctrl = p.case
#       lo.ctrl = p.case
#       
#       for (i in 1:length(m.y.grid)){
#         curr.m.y = m.y.grid[[i]]
#         m.tmp = curr.m.y[1]
#         y.tmp = curr.m.y[2]
#         ind.case = which(month(date.case)==m.tmp & year(date.case)==y.tmp)
#         ind.ctrl = which(month(date.ctrl)==m.tmp & year(date.ctrl)==y.tmp)
#         n.case[i] = length(ind.case)
#         n.ctrl[i] = length(ind.ctrl)
#         for (j in 1:length(Pathogen)){
#           p.case[j,i] = mean(datcase.season[ind.case,j])
#           p.ctrl[j,i] = mean(datctrl.season[ind.ctrl,j])
#           np.case[j,i] = sum(datcase.season[ind.case,j])
#           np.ctrl[j,i] = sum(datctrl.season[ind.ctrl,j])
#           
#           tmp.case = binom.confint(np.case[j,i], n.case[i], 
#                                    conf.level = 0.95, methods = "bayes")
#           tmp.ctrl = binom.confint(np.ctrl[j,i], n.ctrl[i], 
#                                    conf.level = 0.95, methods = "bayes")
#           
#           lo.case[j,i] = tmp.case$lower
#           up.case[j,i] = tmp.case$upper                  
#           lo.ctrl[j,i] = tmp.ctrl$lower
#           up.ctrl[j,i] = tmp.ctrl$upper
#         }
#       }
#       
#       op0<-par(mfrow=c(2,1),oma=par("oma")+c(3,1,3,5),xpd=TRUE)
#       op <-par(mar=c(0,5,3,10))
#       
#       matplot(1:length(uniq.month),t(p.case[c(pathind),]),type="l",ylab="prevalence",
#               xaxt = "n",xlab="",lty=pathind,col=pathind,pch=pathind,las=2,
#               ylim=c(0,1),lwd=3)
#       
#       for (s in pathind){
#         polygon(c(1:length(uniq.month),rev(1:length(uniq.month))),
#                 c(lo.case[s,],rev(up.case[s,])),lty=2,border=s)
#       }
#       
#       
#       text(i+4,.9,"cases",cex=3)
#       for (i in 1:length(uniq.month)){
#         text(i,0,n.case[i])
#       }
#       text(i+2,0,"#cases")
#       par(op)
#       #axis(1,at = 1:length(uniq.month),labels=uniq.month,las=2)
#       op<-par(mar=c(6,5,0,10),xpd=T)
#       matplot(t(p.ctrl[c(pathind),]),type="l",ylab="prevalence",
#               xlab="",xaxt = "n",lty=pathind,col=pathind,
#               pch=pathind,las=2,ylim=c(0,1),lwd=3)
#       
#       for (s in pathind){
#         polygon(c(1:length(uniq.month),rev(1:length(uniq.month))),
#                 c(lo.ctrl[s,],rev(up.ctrl[s,])),lty=2,border=s)
#       }
#       
#       
#       text(i+4,0,"ctrls",cex=3,las=2)
#       mtext("month-year",1,line=5,cex=1,las=1)
#       axis(1,at = 1:length(uniq.month),labels=uniq.month,las=2)
#       mtext(paste0("Site: ", sitename),outer=TRUE,cex=2)
#       legend(length(uniq.month)+1,0.9,
#              c(Pathogen[pathind]),lty=pathind,col=pathind,
#              bty="n",lwd=rep(3,length(pathind)))
#       for (i in 1:length(uniq.month)){
#         text(i,1,n.ctrl[i])
#       }
#       text(i+2,1,"#control")
#     }
#     par(op0)
#   } else {
    
    intervals = seq(min(c(date.case,date.ctrl)),
                    max(c(date.case,date.ctrl)),by=interval)
    
    if (max(intervals)!=max(c(date.case,date.ctrl))){
      intervals = c(intervals,max(c(date.case,date.ctrl))+1)
    }
    
    mid.dates = as.Date(rep(NA,length(intervals)-1))
    for (s in 1:(length(intervals)-1)){
      mid.dates[s] = intervals[s]+interval/2
      if (s == (length(intervals)-1) ){
        mid.dates[s] = intervals[s]+(intervals[s+1]-intervals[s])/2
      }
    }
    
    p.case = matrix(NA,nrow=length(Pathogen),ncol=length(intervals)-1)
    rownames(p.case) = Pathogen
    colnames(p.case) = as.character(mid.dates)
    p.ctrl = p.case
    
    n.case = rep(NA,length=length(intervals)-1)
    names(n.case) = as.character(mid.dates)
    n.ctrl = n.case
    
    ORseason = p.case
    
    np.case = p.case
    np.ctrl = p.case
    up.case = p.case
    lo.case = p.case
    up.ctrl = p.case
    lo.ctrl = p.case
    
    for (i in 1:(length(intervals)-1)){
      ind.case = which(date.case>=intervals[i] & date.case<intervals[i+1])
      ind.ctrl = which(date.ctrl>=intervals[i] & date.ctrl<intervals[i+1])
      n.case[i] = length(ind.case)
      n.ctrl[i] = length(ind.ctrl)
      for (j in 1:length(Pathogen)){
        p.case[j,i] = mean(datcase.season[ind.case,j])
        p.ctrl[j,i] = mean(datctrl.season[ind.ctrl,j])
        np.case[j,i] = sum(datcase.season[ind.case,j])
        np.ctrl[j,i] = sum(datctrl.season[ind.ctrl,j])
        
        tmp.case = binom.confint(np.case[j,i], n.case[i], conf.level = 0.95,
                                 methods = "bayes",tol=1e-11)
        tmp.ctrl = binom.confint(np.ctrl[j,i], n.ctrl[i], conf.level = 0.95, 
                                 methods = "bayes",tol=1e-11)
        
        lo.case[j,i] = tmp.case$lower
        up.case[j,i] = tmp.case$upper                  
        lo.ctrl[j,i] = tmp.ctrl$lower
        up.ctrl[j,i] = tmp.ctrl$upper
        
        ORseason[j,i] = ((np.case[j,i]+0.5)/(n.case[i]-np.case[j,i]+0.5))/((np.ctrl[j,i]+0.5)/(n.ctrl[i]-np.ctrl[j,i]+0.5))
      }
    }
    if (fig==TRUE){
          op0<-par(mfrow=c(2,1),oma=par("oma")+c(3,1,3,5),xpd=TRUE)
          op <-par(mar=c(0,5,3,10))
          
          plot(mid.dates,p.case[c(pathind),],type="l",ylab="prevalence",xlab="",
               lty=pathind,col=pathind,pch=pathind,las=2,xaxt="n",
                  ylim=c(0,1),lwd=3)
          
          for (s in pathind){
            polygon(c(mid.dates,rev(mid.dates)),
                    c(lo.case[s,],rev(up.case[s,])),lty=2,border=s)
          }
          
          
          text(tail(intervals,1)+25,.9,"cases",cex=3)
          for (i in 1:length(intervals[-1])){
            text(mid.dates[i],0,n.case[i])
          }
          text(tail(intervals,1)+25+2,0,"#cases")
          par(op)
          #axis(1,at = 1:length(uniq.month),labels=uniq.month,las=2)
          op<-par(mar=c(6,5,0,10),xpd=T)
          plot(mid.dates,p.ctrl[c(pathind),],type="l",ylab="prevalence",
                  xlab="",xaxt = "n",lty=pathind,col=pathind,
                  pch=pathind,las=2,ylim=c(0,1),lwd=3)
          
          for (s in pathind){
            polygon(c(mid.dates,rev(mid.dates)),
                    c(lo.ctrl[s,],rev(up.ctrl[s,])),lty=2,border=s)
          }
          
          
          text(tail(intervals,1)+25+2,0,"ctrls",cex=3,las=2)
          mtext("month-year",1,line=5,cex=1,las=1)
          axis(1,at = mid.dates,labels=mid.dates,las=2)
          mtext(paste0("Site: ", sitename),outer=TRUE,cex=2)
          legend(length(intervals[-1])+1,0.9,
                 c(Pathogen[pathind]),lty=pathind,col=pathind,
                 bty="n",lwd=rep(3,length(pathind)))
          for (i in 1:length(intervals[-1])){
            text(mid.dates[i],1,n.ctrl[i])
            for (j in pathind){
              text(mid.dates[i],0.8*(j/length(Pathogen)),round(ORseason[j,i],2),
                   col=j)
            }
          }
          text(tail(intervals,1)+25+2,1,"#control")
          par(op0)
          par(op)
    }
  #}
  res = list(p.case,p.ctrl,n.case,n.ctrl,np.case,np.ctrl,mid.dates,intervals)
  names(res) = c("p.case","p.ctrl","n.case","n.ctrl","np.case","np.ctrl","mid.dates","intervals")
  return(res)
}

infograph_popPie_ord = function(directoryName,
                                alphaE=alpha,
                                alphaB = rep(1,Jfull),
                                betaB = rep(1,Jfull),
                                alphaG = rep(1,JGS),
                                betaG  = rep(1,JGS),
                                Bdat=MBS,Gdat=MGS,
                                ncase = Nd,J.GS = JGS,
                                nppcr.no.realGS.fit=FALSE,
                                nppcr.no.realGS.fit.dir=prior.folder,#"\\nppcr.no.realGS.fit",
                                nppcr.bcx.fit=TRUE,
                                nppcr.bcx.fit.dir=prior.folder,#"\\nppcr.bcx.fit",
                                pEtiName="pEti",pathogens = Pathogen,
                                top3=0.25){
  
#   
#     directoryName= fullname
#     alphaE = alpha
#     alphaB = alphaB
#     betaB = betaB
#     alphaG = alphaG
#     betaG  = betaG
#     nppcr.bcx.fit=TRUE
#     nppcr.no.realGS.fit = FALSE
#     pathogens=Pathogen
#     ncase = Nd
#     J.GS=JGS
#     pEtiName="pEti"
#     nppcr.bcx.fit.dir=prior.folder
#     Bdat=MBS;Gdat=MGS
#   
  require(R2WinBUGS)
  cexval=1
  srtval=0
  J.full = length(pathogens)
  
  layout(matrix(c(1,2,3),1,3,byrow = TRUE), 
         widths=c(3,2,3),heights=c(8))
  #layout.show(n=3)
  
  curr.folder = paste0(directoryName,nppcr.bcx.fit.dir)
  simsMat = read.coda(paste(curr.folder,"coda1.txt",sep="\\"),
                      paste(curr.folder,"codaIndex.txt",sep="\\"),quiet=TRUE)  
  
  res = matrix(NA,nrow=nrow(simsMat),ncol=J.full)
  SubVarName = rep(NA,J.full)
  for (j in 1:J.full){
    SubVarName[j] = paste(pEtiName,"[",j,"]",sep="")
  }
  pmatBG = simsMat[,SubVarName]
  pmeanBG = colMeans(pmatBG)
  pmeanBG0 = colMeans(pmatBG)
  
  ord = order(pmeanBG)
  #### from now on permute!
  pathogens.ord = pathogens[ord]
  pmeanBG.ord = pmeanBG[ord]
  pmatBG.ord = pmatBG[,ord]
  
  ## quantiles for etiology: outer is 97.5% CI, inner is 50% CI
  pq1BG   = apply(pmatBG,2,quantile,probs=0.025)[ord]
  pq2BG   = apply(pmatBG,2,quantile,probs=0.975)[ord]
  p.innerq1BG   = apply(pmatBG,2,quantile,probs=0.25)[ord]
  p.innerq2BG   = apply(pmatBG,2,quantile,probs=0.75)[ord]
  
  ######################################################
  ##subplot 1: case-control rate comparisons
  ######################################################
  require(binom)
  nctrl = nrow(Bdat)-ncase
  tmp.case = binom.confint(colSums(Bdat[1:ncase,]), ncase, conf.level = 0.95, methods = "ac")[ord,]
  tmp.ctrl = binom.confint(colSums(Bdat[-(1:ncase),]), nctrl, conf.level = 0.95, methods = "ac")[ord,]
  
  ## case and control positive rate, lower and upper limit
  Bcomp = rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5))
  Bcompq1 = rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")])
  Bcompq2 = rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")])
  
  ## posterior distribution of TPR
  thetamatBG = (simsMat[,grep("thetaBS",colnames(simsMat))])[,ord]
  thetameanBG = colMeans(thetamatBG)
  
  ## posterior distribution of FPR
  psimatBG = (simsMat[,grep("psiBS",colnames(simsMat))])[,ord]
  psimeanBG = colMeans(psimatBG)  
  
  ## model fitted postive rate for each pathogen
  fittedmeanBGcase = sapply(1:J.full, 
                            function(s) mean(pmatBG.ord[,s]*thetamatBG[,s]+(1-pmatBG.ord[,s])*psimatBG[,s]))
  fittedmeanBGcontrol = psimeanBG
  
  top2 = 1.3
  op<-par(mar=c(5.1,4.1,4.1,0))
  plotat = seq(0.5,J.full+0.5,by=1/4)[-(c(1,(1:J.full)*4+1))]
  plot(c(rbind(fittedmeanBGcase,Bcomp)),plotat,yaxt="n",xlim=c(0,top2),
       ylim=c(0.5,J.full+.5),xaxt="n",
       ylab="",xlab="probability",
       pch = c(rbind(rep(2,J.full),rep(20,J.full),rep(20,J.full))),
       col=c(rbind(rbind(rep(1,J.full),rep("dodgerblue2",J.full),rep("dodgerblue2",J.full)))),
       cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
  axis(2,at = plotat,labels=rep(c("","case","ctrl"),J.full),las=2)
  axis(1,at = c(0,0.2,0.4,0.6,0.8,1),labels= c(0,0.2,0.4,0.6,0.8,1),las=1)
  points(c(rbind(thetameanBG,Bcompq2)),plotat,
         pch=c(rbind(rep("+",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(2,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  points(c(rbind(psimeanBG,Bcompq1)),plotat,
         pch=c(rbind(rep("*",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(2,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  abline(h=seq(1.5,J.full-.5,by=1),lty=2,lwd=0.5,col="blue")
  abline(v=1,lty=2,lwd=.5)
  
  ## conditional odds ratios
  COR = function(brs.data,nd,pathogens){
    y = rep(c(1,0),times=c(nd,nrow(brs.data)-nd))
    #X = matrix(NA,nrow=nrow(brs.data),ncol=length(pathogens))
    x.nm   = paste(pathogens,"NPPCR",sep="_")
    #colnames(X) = x.nm
    #for (j in 1:length(pathogens)){
    #  X[,j] = brs.data[,x.nm[j]]
    #}
    dat.reg = as.data.frame(cbind(y,brs.data))
    form = as.formula(paste0("y~",paste(x.nm,collapse="+")))
    fit = glm(form,data=dat.reg,family=binomial)
    res = cbind(exp(confint(fit)),exp(fit$coef))[-1,]
  }
  
  tmp = COR(MBS,ncase,pathogens.ord)
    for (s in (1:J.full)){
        L=round(tmp[s,1],1)
        C=round(tmp[s,3],1)
        R=round(tmp[s,2],1)
        text(top2-0.12,s+1/(2*Jfull),C,cex=1.5)
        text(top2-0.12,s-.2,paste(c(L,"   ",R),collapse=" "),cex=1.2)
    }
  legend("topright","conditional OR",bty="n")

#   ## marginal odds ratios
#   for (s in (1:J.full)){
#     if (any(Bcomp[,s]==0)){
#       cat("Pathogen ",pathogens[ord[s]]," have zero cells cross-classifying Y and MBS.")
#       text(top2*.9,s,NA)
#     } else{
#       
#       tmp = LOR(Bdat[,ord[s]],rep(c(1,0),times=c(ncase,nrow(Bdat)-ncase)))
#       L=round(tmp$OR[2],1)
#       C=round(tmp$OR[1],1)
#       R=round(tmp$OR[3],1)
#       text(top2-0.12-0.5,s+1/(2*Jfull),C,cex=1.5)
#       text(top2-0.12-0.5,s-.2,paste(c(L,"   ",R),collapse=" "),cex=1.2)
#     }
#   }
  
  counter = 0
  for (s in 1:(3*J.full)){
    segments(y0=plotat[s],x0=c(rbind(psimeanBG,Bcompq1))[s],
             y1=plotat[s],x1=c(rbind(thetameanBG,Bcompq2))[s],col="black",
             lty=ifelse((s-1)%%3<1,4,1))
    if ((s-1)%%3>=1){
      counter=counter+1
      tmp.hpos = ifelse(c(Bcompq2)[counter]+0.15>0.95,c(Bcompq1)[counter]-0.2,c(Bcompq2)[counter]+0.15 )
      text(tmp.hpos,plotat[s],paste0(round(100*c(Bcomp),1)[counter],"%"),srt=srtval,cex=cexval)
    }
  }
  
  for (s in 1:(J.full)){
    segments(y0=plotat[3*s-1],x0=c(rbind(fittedmeanBGcase,Bcomp))[3*s-1],
             y1=plotat[3*s],x1=c(rbind(fittedmeanBGcase,Bcomp))[3*s],col="black")
  }
  # put prior shapes on bronze sensitivity
  for (s in 1:J.full){
    tmp = rbeta(10000,alphaB[ord[s]],betaB[ord[s]])
    boxplot(tmp,at = s-0.45, boxwex=1/10 , col="gray",add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
    tmp.post = as.matrix(thetamatBG)[,s]
    boxplot(tmp.post,at = s-0.35,boxwex=1/10,add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
  }
  axis(2,at=(1:J.full)-.45,label=rep("",J.full),las=2,cex.axis=.5)
  axis(2,at=(1:J.full)-.35,label=rep("",J.full),las=2,cex.axis=.5)
  
  mtext(expression(underline("BrS")),line=1,cex=1.8)
  par(op)
  
  ######################################################
  ##subplot 2: silver standard data comparisons
  ######################################################
  
  Gcomplete.index = which(rowSums(is.na(Gdat[,1:J.GS]))==0)
  cat("No. of Complete Silver-standard data on ",pathogens[1:J.GS],":",
      length(Gcomplete.index),"\n")
  cat(ncase-length(Gcomplete.index)," cases deleted for plotting third subgraph.","\n")
  Gdat = Gdat[Gcomplete.index,]
  ord.gs = order(pmeanBG0[1:J.GS])
  ind.gs = rep(NA,J.GS) # tells where the the gs row should go
  for (j in 1:J.GS){
    ind.gs[j] = which(ord==j)
  }
  
  tmpG.case = binom.confint(colSums(Gdat[,1:J.GS]), nrow(Gdat), conf.level = 0.95, methods = "ac")
  Gcomp = rbind(round(tmpG.case$mean,5),rep(NA,J.GS))
  Gcompq1 = rbind(tmpG.case[,c("lower")],rep(NA,J.GS))
  Gcompq2 = rbind(tmpG.case[,c("upper")],rep(NA,J.GS))
  
  thetamatG = (simsMat[,grep("thetaGS",colnames(simsMat))])########### stopped here
  thetameanG = colMeans(thetamatG)
  
  thetamatGq1=apply(thetamatG,2,quantile,0.025)
  thetamatGq2=apply(thetamatG,2,quantile,0.975)
  
  fittedmeanGpos = sapply(1:J.GS, function(s) mean(pmatBG.ord[,ind.gs[s]]*thetamatG[,s]))
  
  
  #top3 = .25
  par(mar=c(5.1,0,4.1,0))
  #op<-par(mar=c(5.1,0,4.1,2.1),xpd=TRUE)
  
  plotat = seq(0.5,J.full+0.5,by=1/4)[-(c(1,(1:J.full)*4+1))]
  #plotat.short = plotat[1:length(c(rbind(thetameanG,Gcomp)))]
  plotat.calc = function(j) {c(3*j-2,3*j-1,3*j)}
  plotat.short = rep(NA,J.GS*3)
  for (j in 1:J.GS){
    plotat.short[c(3*j-2,3*j-1,3*j)] = plotat[plotat.calc(ind.gs[j])]
  }
  
  plot(c(rbind(thetameanG,Gcomp)),plotat.short,yaxt="n",xlim=c(0,top3),
       ylim=c(0.5,J.full+.5),#xaxt="n",
       ylab="",xlab="probability",
       pch = c(rbind(rep(20,J.full),rep(20,J.full),rep(20,J.full))),
       col=c(rbind(rbind(rep(1,J.full),rep("blue",J.full),rep(1,J.full)))),
       cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
  #axis(1,at = seq(0,top3,len=10),labels= seq(0,top3,len=10),las=1)
  points(c(rbind(fittedmeanGpos,matrix("",nrow=2,ncol=J.GS))),plotat.short,yaxt="n",xlim=c(0,top3),
         ylim=c(0.5,J.full+.5),xaxt="n",
         ylab="",#xlab="Gold Positive Rate",
         pch = c(rbind(rep(2,J.full),rep(NA,J.full),rep(NA,J.full))),
         col=c(rbind(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full)))),
         cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
  
  #axis(2,at = plotat.short,labels=rep(c("TPR-G","case","ctrl"),J.GS),las=2)
  points(c(rbind(thetamatGq2,Gcompq2)),plotat.short,
         pch=c(rbind(rep("|",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  points(c(rbind(thetamatGq1,Gcompq1)),plotat.short,
         pch=c(rbind(rep("|",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  
  #inner 25%-75%
  thetamatGq.inner1=apply(thetamatG,2,quantile,0.25)
  thetamatGq.inner2=apply(thetamatG,2,quantile,0.75)
  points(c(rbind(thetamatGq.inner1,Gcompq1)),plotat.short,
         pch=c(rbind(rep("[",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  points(c(rbind(thetamatGq.inner2,Gcompq1)),plotat.short,
         pch=c(rbind(rep("]",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  counter = 0
  for (s in 1:length(plotat.short)){
    segments(y0=plotat.short[s],x0=c(rbind(thetamatGq.inner1,Gcompq1))[s],
             y1=plotat.short[s],x1=c(rbind(thetamatGq.inner2,Gcompq2))[s],
             col="black",
             lwd=1)
  }
  
  # row separation lines
  abline(h=seq(1.5,J.full-.5,by=1)[ind.gs],lty=2,lwd=0.5,col="blue")
  abline(h=seq(1.5,J.full-.5,by=1)[ind.gs]-1,lty=2,lwd=0.5,col="blue")
  
  
  counter = 0
  for (s in 1:length(plotat.short)){
    segments(y0=plotat.short[s],x0=c(rbind(thetamatGq1,Gcompq1))[s],
             y1=plotat.short[s],x1=c(rbind(thetamatGq2,Gcompq2))[s],col="black",
             lty=ifelse((s-1)%%3<2,1,1))
    if ((s-1)%%3>=1){
      counter=counter+1
      text(c(Gcomp)[counter],plotat.short[s]+0.125,paste0(round(100*c(Gcomp),1)[counter],"%"),srt=srtval,cex=cexval)
    }
  }
  
  
  for (s in 1:J.GS){
    
    text(thetameanG[s],plotat.short[3*s-2]+.125,paste(round(100*thetameanG[s],2),"%"))
    
    # put prior shapes on gold sensitivity
    tmp = rbeta(10000,alphaG[s],betaG[s])
    boxplot(tmp,at = ind.gs[s]-0.45, boxwex=1/8 ,col="gray", 
            add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
    
  }
  
  mtext(expression(underline("SS")),line=1,cex=1.8)
  
  par(op)
  
  
  ######################################################
  ##subplot 3: etiology bars with prior shapes
  ######################################################
  top=0.4#max(pq2BG)
  dotcolor = "black"
  #op <- par(mar=c(5.1,6,4.1,1.1))
  op <- par(mar=c(5.1,0,4.1,9))
  plot(c(pmeanBG.ord),1:(J.full),
       yaxt="n",#xaxt="n",
       xlim=c(0,top),ylim=c(0.5,J.full+0.5),col=c("black"),
       ylab="",xlab="probability",
       pch=c(20),cex=2)
  axis(4,at=1:J.full,labels=pathogens.ord,las=2,cex.axis=1.5)
  abline(h=seq(1.5,J.full-.5,by=1),lty=2,lwd=0.5,col="blue")
  #draw axis within plot:
  for (s in 1:(J.full-1)){
    axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=1,#labels=rep("",length(seq(0,1,by=.2))),
         pos = seq(1.5,J.full-.5,by=1)[s], cex.axis = 0.8,lty=2,col="blue")
    # axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=0,#labels=rep("",length(seq(0,1,by=.2))),
    #      pos = seq(1.5,J.full-.5,by=1)[s]+0.3, cex.axis = 0.8,lty=2,col="blue")
  }
  points(c(pq1BG),1:(J.full),pch="|",cex=1)
  points(c(pq2BG),1:(J.full),pch="|",cex=1)
  points(c(p.innerq1BG),1:(J.full),pch="[",cex=1)
  points(c(p.innerq2BG),1:(J.full),pch="]",cex=1)
  
  mtext(expression(underline(hat(pi))),line=1,cex=1.8)
  #mtext(c(expression(bold("--")),":prior","-",":posterior"),col=c("gray","black","black","black"),
  #      adj=c(0,0.1,0.3,0.4),line=.8,cex=.8,lwd=2)
  legend("topright",c("prior","posterior"),lty=c(2,1),col=c("gray","black"),
         lwd = 4,horiz=TRUE,cex=1.5,bty="n")
  pgrid = seq(0,1,by=0.01)
  
  for (s in 1:(J.full)){
    segments(y0=s,x0=c(pq1BG)[s],y1=s,x1=c(pq2BG)[s],col=dotcolor)
    segments(y0=s,x0=c(p.innerq1BG)[s],y1=s,x1=c(p.innerq2BG)[s],col=dotcolor,lwd=2)
    #text(pmeanBG[s],s+0.25,paste0(round(100*c(pmeanBG),2)[s],"%"),srt=srtval,cex=cexval)
    text(.35,s+0.25,paste0("=",paste0(round(100*c(pmeanBG.ord),1)[s],"%")),srt=srtval,cex=2)
    text(.27,s+0.25,bquote(hat(pi)[.(ord[s])]),srt=srtval,cex=2)
    #text(top-0.05,s,pathogens[s],cex=.8,srt=srtval)
    tmp.density = dbeta(pgrid,alphaE[ord[s]],sum(alphaE[-ord[s]]))
    points(pgrid,tmp.density/(3*max(tmp.density))+s-0.45,type="l",col="gray",lwd=4,lty=2)
    ##posterior density
    tmp.post.density = density(pmatBG.ord[,s],from=0,to=1)
    tmp.x = tmp.post.density$x
    tmp.y = tmp.post.density$y
    points(tmp.x,tmp.y/(3*max(tmp.y))+s-0.45,col="black",lwd=4,type="l")
    
  }
  
  par(op)
  
  
  #######################################################
  #######################################################
  #######################################################
  #######################################################
  #######################################################
  
}

visualize_poppie_pathogen_category = function(directoryName,
                                              alphaE=alpha,
                                              alphaB = rep(1,Jfull),
                                              betaB = rep(1,Jfull),
                                              alphaG = rep(1,JGS),
                                              betaG  = rep(1,JGS),
                                              Bdat=MBS,Gdat=MGS,
                                              ncase = Nd,J.GS = JGS,
                                              nppcr.no.realGS.fit=FALSE,
                                              nppcr.no.realGS.fit.dir=prior.folder,#"\\nppcr.no.realGS.fit",
                                              nppcr.bcx.fit=TRUE,
                                              nppcr.bcx.fit.dir=prior.folder,#"\\nppcr.bcx.fit",
                                              pEtiName="pEti",pathogens = Pathogen,
                                              ksFrac=0.1){
  
  ## load packages
  pckg = try(require(MASS))
  if(!pckg) {
    cat("Installing 'MASS' from CRAN\n")
    getPckg("MASS")
    require("MASS")
  }
  
  pckg = try(require(robCompositions))
  if(!pckg) {
    cat("Installing 'robCompositions' from CRAN\n")
    getPckg("robCompositions")
    require("robCompositions")
  }
  
  pckg = try(require(ks))
  if(!pckg) {
    cat("Installing 'ks' from CRAN\n")
    getPckg("ks")
    require("ks")
  }
  
  pckg = try(require(compositions))
  if(!pckg) {
    cat("Installing 'compositions' from CRAN\n")
    getPckg("compositions")
    require("compositions")
  }
  
  pckg = try(require(R2WinBUGS))
  if(!pckg) {
    cat("Installing 'compositions' from CRAN\n")
    getPckg("R2WinBUGS")
    require("R2WinBUGS")
  }
    
#       directoryName= fullname
#       alphaE = alpha
#       alphaB = alphaB
#       betaB = betaB
#       alphaG = alphaG
#       betaG  = betaG
#       nppcr.bcx.fit=TRUE
#       nppcr.no.realGS.fit = FALSE
#       pathogens=Pathogen
#       ncase = Nd
#       J.GS=JGS
#       pEtiName="pEti"
#       nppcr.bcx.fit.dir=prior.folder
#       Bdat=MBS;Gdat=MGS
  
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
         widths=c(4,4), heights=c(5,3))
  #layout.show(n=3)
  cexval=1
  srtval=0
  J.full = length(pathogens)
  
  #layout.show(n=3)
  
  curr.folder = paste0(directoryName,nppcr.bcx.fit.dir)
  simsMat = read.coda(paste(curr.folder,"coda1.txt",sep="\\"),
                      paste(curr.folder,"codaIndex.txt",sep="\\"),quiet=TRUE)  
  
  res = matrix(NA,nrow=nrow(simsMat),ncol=J.full)
  SubVarName = rep(NA,J.full)
  for (j in 1:J.full){
    SubVarName[j] = paste(pEtiName,"[",j,"]",sep="")
  }
  pmatBG = simsMat[,SubVarName]
  pmeanBG = colMeans(pmatBG)
  pmeanBG0 = colMeans(pmatBG)
  
  ord = order(pmeanBG)
  #### from now on permute!
  pathogens.ord = pathogens[ord]
  pmeanBG.ord = pmeanBG[ord]
  pmatBG.ord = pmatBG[,ord]
  ##############################################
  ## plot 1: pr (virus is a cause)
  ##############################################
  op<-par()
  par(mai=c(0.5,3,1.09333,3))
  pr.v = rowSums(pmatBG[,ord[which(ord>J.GS)]])
  plot(density(pr.v),type="l",xlim=c(0,1),xlab="",bty="n",
       ylab="",yaxt="n",lwd=5,main="",
       cex.lab=2,cex.axis=2)
  tmp.rnd = rbeta(100000,sum(alphaE[-(1:J.GS)]),sum(alphaE[(1:J.GS)]))
  points(density(tmp.rnd),type="l",lty=2,col="grey",lwd=5)
  mtext("Probability of Viral Cause",side=3,line=0,cex=2)
  
  op


#################################################
## plot 2: pr(B1, B2, others | cause is a Bacteria)
################################################
par(mar = c(6.1,0,1,4.1),mai=c(.5,1,.5,.5),xpd=TRUE)
pEti0 = pmatBG[,ord[which(ord<J.GS+1)]]
pEti1  = data.frame(others=rowSums(pEti0[,-c(ncol(pEti0)-1,ncol(pEti0))]),
                    B2 = pEti0[,ncol(pEti0)-1],
                    B1 = pEti0[,ncol(pEti0)] )
pEti = pEti1/rowSums(pEti1)
names(pEti) = c("B-rest",
                rev(pathogens.ord[which(ord<J.GS+1)])[2],
                rev(pathogens.ord[which(ord<J.GS+1)])[1])
#ksFrac = 1
levellabel = c(5,20,50)
postMean = colMeans(pEti)

ternaryDiag(
  as.data.frame(matrix(postMean,nrow=1,ncol=3)),
  pch = "",
  col = "black",
  cex=3,name=rep(NA,3),
  lwd=3)

ind.temp = floor(seq(1,nrow(pEti),len=floor(ksFrac*nrow(pEti))))
temp = compositions::ilr(pEti)[ind.temp,]

coord1 = seq(0.001,0.999,by=0.01)
coord2 = sqrt(0.75)*coord1

coord.grid = expand.grid(coord1=coord1,coord2=coord2)
prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                             y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                             x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                               1/sqrt(3)*coord.grid$coord2)
prob.grid = rev(prob.grid.temp)

ilr.grid    = compositions::ilr(prob.grid)

fhat1       = kde(x=temp,H=Hpi(temp),compute.cont=T)
fhat2       = kde(x=temp,H=Hpi(temp),eval.points=ilr.grid)

tri.z       = matrix(fhat2$estimate,nr=length(coord1),nc=length(coord2))

f = function(x,y){
  (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
}

for (i in 1:length(coord1)){
  for (j in 1:length(coord2)){
    tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA) 
  }  
}


num.levels      = levellabel#c(5,10,20,50)
cont.levels     = paste(num.levels,"%",sep="")
label.levels     = paste(100-num.levels,"%",sep="")

points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
       cex=2,col="red",pch=15)
contour(coord1,coord2,tri.z,
        levels =(fhat1$cont)[cont.levels],
        labels = label.levels,add=TRUE,lwd=4,col=c("blue"),lty=1)
mcex=2
name = names(pEti)
mtext(name[1], side = 1, line = 1, at = -0.1, cex = mcex)
mtext(name[2], side = 1, line = 1, at = 1.1, cex = mcex)
mtext(name[3], side = 3,  line= -1, at=0.5 ,cex = mcex)

  #################################################
  ## plot 3: pr(V1, V2, others | cause is a virus)
  ################################################
  op<-par()
  par(mar = c(6.1,4.1,1,2.1),mai=c(.5,.5,.5,1))
  pEti0 = pmatBG[,ord[which(ord>J.GS)]]
  pEti1  = data.frame(V2 = pEti0[,ncol(pEti0)-1],
                      others=rowSums(pEti0[,-c(ncol(pEti0)-1,ncol(pEti0))]),
                      V1 = pEti0[,ncol(pEti0)])
  pEti = pEti1/rowSums(pEti1)
  names(pEti) = c(rev(pathogens.ord[which(ord>J.GS)])[2],"V-rest",
                  rev(pathogens.ord[which(ord>J.GS)])[1])
  #ksFrac = 1
  levellabel = c(5,20,50)
  postMean = colMeans(pEti)
  
  ternaryDiag(
    as.data.frame(matrix(postMean,nrow=1,ncol=3)),
    pch = "",
    col = "black",
    cex=3,name=rep(NA,3),
    lwd=3)
  
  ind.temp = floor(seq(1,nrow(pEti),len=floor(ksFrac*nrow(pEti))))
  temp = compositions::ilr(pEti)[ind.temp,]
  
  coord1 = seq(0.001,0.999,by=0.01)
  coord2 = sqrt(0.75)*coord1
  
  coord.grid = expand.grid(coord1=coord1,coord2=coord2)
  prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                               y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                               x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                                 1/sqrt(3)*coord.grid$coord2)
  prob.grid = rev(prob.grid.temp)
  
  ilr.grid    = compositions::ilr(prob.grid)
  
  fhat1       = kde(x=temp,H=Hpi(temp),compute.cont=T)
  fhat2       = kde(x=temp,H=Hpi(temp),eval.points=ilr.grid)
  
  tri.z       = matrix(fhat2$estimate,nr=length(coord1),nc=length(coord2))
  
  f = function(x,y){
    (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
  }
  
  for (i in 1:length(coord1)){
    for (j in 1:length(coord2)){
      tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA) 
    }  
  }
  
  
  num.levels      = levellabel#c(5,10,20,50)
  cont.levels     = paste(num.levels,"%",sep="")
  label.levels     = paste(100-num.levels,"%",sep="")
  
  points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
         cex=2,col="red",pch=15)
  contour(coord1,coord2,tri.z,
          levels =(fhat1$cont)[cont.levels],
          labels = label.levels,add=TRUE,lwd=4,col=c("blue"),lty=1)
  mcex=2
  name = names(pEti)
  mtext(name[1], side = 1, line = 1, at = -0.1, cex = mcex)
  mtext(name[2], side = 1, line = 1, at = 1.1, cex = mcex)
  mtext(name[3], side = 3,  line= -1, at=0.5 ,cex = mcex)
  op
  
  
}





infograph_popPie_ord_fixTPR = function(directoryName,
                                alphaE=alpha,
                                thetaBS = thetaBS,
                                alphaG = rep(1,JGS),
                                betaG  = rep(1,JGS),
                                Bdat=MBS,Gdat=MGS,
                                ncase = Nd,J.GS = JGS,
                                nppcr.no.realGS.fit=FALSE,
                                nppcr.no.realGS.fit.dir=prior.folder,#"\\nppcr.no.realGS.fit",
                                nppcr.bcx.fit=TRUE,
                                nppcr.bcx.fit.dir=prior.folder,#"\\nppcr.bcx.fit",
                                pEtiName="pEti",pathogens = Pathogen,
                                top3=0.25){
  
#     
#       directoryName= fullname
#       alphaE = alpha
#       alphaB = alphaB
#       betaB = betaB
#       alphaG = alphaG
#       betaG  = betaG
#       nppcr.bcx.fit=TRUE
#       nppcr.no.realGS.fit = FALSE
#       pathogens=Pathogen
#       ncase = Nd
#       J.GS=JGS
#       pEtiName="pEti"
#       nppcr.bcx.fit.dir=prior.folder
#       Bdat=MBS;Gdat=MGS
    
  require(R2WinBUGS)
  cexval=1
  srtval=0
  J.full = length(pathogens)
  
  layout(matrix(c(1,2,3),1,3,byrow = TRUE), 
         widths=c(3,2,3),heights=c(8))
  #layout.show(n=3)
  
  curr.folder = paste0(directoryName,nppcr.bcx.fit.dir)
  simsMat = read.coda(paste(curr.folder,"coda1.txt",sep="\\"),
                      paste(curr.folder,"codaIndex.txt",sep="\\"),quiet=TRUE)  
  
  res = matrix(NA,nrow=nrow(simsMat),ncol=J.full)
  SubVarName = rep(NA,J.full)
  for (j in 1:J.full){
    SubVarName[j] = paste(pEtiName,"[",j,"]",sep="")
  }
  pmatBG = simsMat[,SubVarName]
  pmeanBG = colMeans(pmatBG)
  pmeanBG0 = colMeans(pmatBG)
  
  ord = order(pmeanBG)
  #### from now on permute!
  pathogens.ord = pathogens[ord]
  pmeanBG.ord = pmeanBG[ord]
  pmatBG.ord = pmatBG[,ord]
  
  ## quantiles for etiology: outer is 97.5% CI, inner is 50% CI
  pq1BG   = apply(pmatBG,2,quantile,probs=0.025)[ord]
  pq2BG   = apply(pmatBG,2,quantile,probs=0.975)[ord]
  p.innerq1BG   = apply(pmatBG,2,quantile,probs=0.25)[ord]
  p.innerq2BG   = apply(pmatBG,2,quantile,probs=0.75)[ord]
  
  ######################################################
  ##subplot 1: case-control rate comparisons
  ######################################################
  require(binom)
  nctrl = nrow(Bdat)-ncase
  tmp.case = binom.confint(colSums(Bdat[1:ncase,]), ncase, conf.level = 0.95, methods = "ac")[ord,]
  tmp.ctrl = binom.confint(colSums(Bdat[-(1:ncase),]), nctrl, conf.level = 0.95, methods = "ac")[ord,]
  
  ## case and control positive rate, lower and upper limit
  Bcomp = rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5))
  Bcompq1 = rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")])
  Bcompq2 = rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")])
  
  ## posterior distribution of TPR
  thetamatBG = matrix(thetaBS,nrow=nrow(pmatBG),ncol=J.full)#(simsMat[,grep("thetaBS",colnames(simsMat))])[,ord]
  thetameanBG = rep(thetaBS,J.full)
  
  ## posterior distribution of FPR
  psimatBG = (simsMat[,grep("psiBS",colnames(simsMat))])[,ord]
  psimeanBG = colMeans(psimatBG)  
  
  ## model fitted postive rate for each pathogen
  fittedmeanBGcase = sapply(1:J.full, 
                            function(s) mean(pmatBG.ord[,s]*thetamatBG[,s]+(1-pmatBG.ord[,s])*psimatBG[,s]))
  fittedmeanBGcontrol = psimeanBG
  
  top2 = 1.3
  op<-par(mar=c(5.1,4.1,4.1,0))
  plotat = seq(0.5,J.full+0.5,by=1/4)[-(c(1,(1:J.full)*4+1))]
  plot(c(rbind(fittedmeanBGcase,Bcomp)),plotat,yaxt="n",xlim=c(0,top2),
       ylim=c(0.5,J.full+.5),xaxt="n",
       ylab="",xlab="probability",
       pch = c(rbind(rep(2,J.full),rep(20,J.full),rep(20,J.full))),
       col=c(rbind(rbind(rep(1,J.full),rep("dodgerblue2",J.full),rep("dodgerblue2",J.full)))),
       cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
  axis(2,at = plotat,labels=rep(c("","case","ctrl"),J.full),las=2)
  axis(1,at = c(0,0.2,0.4,0.6,0.8,1),labels= c(0,0.2,0.4,0.6,0.8,1),las=1)
  points(c(rbind(thetameanBG,Bcompq2)),plotat,
         pch=c(rbind(rep("+",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(2,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  points(c(rbind(psimeanBG,Bcompq1)),plotat,
         pch=c(rbind(rep("*",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(2,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  abline(h=seq(1.5,J.full-.5,by=1),lty=2,lwd=0.5,col="blue")
  abline(v=1,lty=2,lwd=.5)
  
  ## conditional odds ratios
  COR = function(brs.data,nd,pathogens){
    y = rep(c(1,0),times=c(nd,nrow(brs.data)-nd))
    #X = matrix(NA,nrow=nrow(brs.data),ncol=length(pathogens))
    x.nm   = paste(pathogens,"NPPCR",sep="_")
    #colnames(X) = x.nm
    #for (j in 1:length(pathogens)){
    #  X[,j] = brs.data[,x.nm[j]]
    #}
    dat.reg = as.data.frame(cbind(y,brs.data))
    form = as.formula(paste0("y~",paste(x.nm,collapse="+")))
    fit = glm(form,data=dat.reg,family=binomial)
    res = cbind(exp(confint(fit)),exp(fit$coef))[-1,]
  }
  
  tmp = COR(MBS,ncase,pathogens.ord)
  for (s in (1:J.full)){
    L=round(tmp[s,1],1)
    C=round(tmp[s,3],1)
    R=round(tmp[s,2],1)
    text(top2-0.12,s+1/(2*Jfull),C,cex=1.5)
    text(top2-0.12,s-.2,paste(c(L,"   ",R),collapse=" "),cex=1.2)
  }
  legend("topright","conditional OR",bty="n")
  
  #   ## marginal odds ratios
  #   for (s in (1:J.full)){
  #     if (any(Bcomp[,s]==0)){
  #       cat("Pathogen ",pathogens[ord[s]]," have zero cells cross-classifying Y and MBS.")
  #       text(top2*.9,s,NA)
  #     } else{
  #       
  #       tmp = LOR(Bdat[,ord[s]],rep(c(1,0),times=c(ncase,nrow(Bdat)-ncase)))
  #       L=round(tmp$OR[2],1)
  #       C=round(tmp$OR[1],1)
  #       R=round(tmp$OR[3],1)
  #       text(top2-0.12-0.5,s+1/(2*Jfull),C,cex=1.5)
  #       text(top2-0.12-0.5,s-.2,paste(c(L,"   ",R),collapse=" "),cex=1.2)
  #     }
  #   }
  
  counter = 0
  for (s in 1:(3*J.full)){
    segments(y0=plotat[s],x0=c(rbind(psimeanBG,Bcompq1))[s],
             y1=plotat[s],x1=c(rbind(thetameanBG,Bcompq2))[s],col="black",
             lty=ifelse((s-1)%%3<1,4,1))
    if ((s-1)%%3>=1){
      counter=counter+1
      tmp.hpos = ifelse(c(Bcompq2)[counter]+0.15>0.95,c(Bcompq1)[counter]-0.2,c(Bcompq2)[counter]+0.15 )
      text(tmp.hpos,plotat[s],paste0(round(100*c(Bcomp),1)[counter],"%"),srt=srtval,cex=cexval)
    }
  }
  
  for (s in 1:(J.full)){
    segments(y0=plotat[3*s-1],x0=c(rbind(fittedmeanBGcase,Bcomp))[3*s-1],
             y1=plotat[3*s],x1=c(rbind(fittedmeanBGcase,Bcomp))[3*s],col="black")
  }
  # put prior shapes on bronze sensitivity
  for (s in 1:J.full){
    tmp = rep(thetaBS,100)
    boxplot(tmp,at = s-0.45, boxwex=1/10 , col="gray",add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
    tmp.post = as.matrix(thetamatBG)[,s]
    boxplot(tmp.post,at = s-0.35,boxwex=1/10,add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
  }
  axis(2,at=(1:J.full)-.45,label=rep("",J.full),las=2,cex.axis=.5)
  axis(2,at=(1:J.full)-.35,label=rep("",J.full),las=2,cex.axis=.5)
  
  mtext(expression(underline("BrS")),line=1,cex=1.8)
  par(op)
  
  ######################################################
  ##subplot 2: silver standard data comparisons
  ######################################################
  
  Gcomplete.index = which(rowSums(is.na(Gdat[,1:J.GS]))==0)
  cat("No. of Complete Silver-standard data on ",pathogens[1:J.GS],":",
      length(Gcomplete.index),"\n")
  cat(ncase-length(Gcomplete.index)," cases deleted for plotting third subgraph.","\n")
  Gdat = Gdat[Gcomplete.index,]
  ord.gs = order(pmeanBG0[1:J.GS])
  ind.gs = rep(NA,J.GS) # tells where the the gs row should go
  for (j in 1:J.GS){
    ind.gs[j] = which(ord==j)
  }
  
  tmpG.case = binom.confint(colSums(Gdat[,1:J.GS]), nrow(Gdat), conf.level = 0.95, methods = "ac")
  Gcomp = rbind(round(tmpG.case$mean,5),rep(NA,J.GS))
  Gcompq1 = rbind(tmpG.case[,c("lower")],rep(NA,J.GS))
  Gcompq2 = rbind(tmpG.case[,c("upper")],rep(NA,J.GS))
  
  thetamatG = (simsMat[,grep("thetaGS",colnames(simsMat))])########### stopped here
  thetameanG = colMeans(thetamatG)
  
  thetamatGq1=apply(thetamatG,2,quantile,0.025)
  thetamatGq2=apply(thetamatG,2,quantile,0.975)
  
  fittedmeanGpos = sapply(1:J.GS, function(s) mean(pmatBG.ord[,ind.gs[s]]*thetamatG[,s]))
  
  
  #top3 = .25
  par(mar=c(5.1,0,4.1,0))
  #op<-par(mar=c(5.1,0,4.1,2.1),xpd=TRUE)
  
  plotat = seq(0.5,J.full+0.5,by=1/4)[-(c(1,(1:J.full)*4+1))]
  #plotat.short = plotat[1:length(c(rbind(thetameanG,Gcomp)))]
  plotat.calc = function(j) {c(3*j-2,3*j-1,3*j)}
  plotat.short = rep(NA,J.GS*3)
  for (j in 1:J.GS){
    plotat.short[c(3*j-2,3*j-1,3*j)] = plotat[plotat.calc(ind.gs[j])]
  }
  
  plot(c(rbind(thetameanG,Gcomp)),plotat.short,yaxt="n",xlim=c(0,top3),
       ylim=c(0.5,J.full+.5),#xaxt="n",
       ylab="",xlab="probability",
       pch = c(rbind(rep(20,J.full),rep(20,J.full),rep(20,J.full))),
       col=c(rbind(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full)))),
       cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
  #axis(1,at = seq(0,top3,len=10),labels= seq(0,top3,len=10),las=1)
  points(c(rbind(fittedmeanGpos,matrix("",nrow=2,ncol=J.GS))),plotat.short,yaxt="n",xlim=c(0,top3),
         ylim=c(0.5,J.full+.5),xaxt="n",
         ylab="",#xlab="Gold Positive Rate",
         pch = c(rbind(rep(2,J.full),rep(NA,J.full),rep(NA,J.full))),
         col=c(rbind(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full)))),
         cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
  
  #axis(2,at = plotat.short,labels=rep(c("TPR-G","case","ctrl"),J.GS),las=2)
  points(c(rbind(thetamatGq2,Gcompq2)),plotat.short,
         pch=c(rbind(rep("|",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  points(c(rbind(thetamatGq1,Gcompq1)),plotat.short,
         pch=c(rbind(rep("|",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  
  #inner 25%-75%
  thetamatGq.inner1=apply(thetamatG,2,quantile,0.25)
  thetamatGq.inner2=apply(thetamatG,2,quantile,0.75)
  points(c(rbind(thetamatGq.inner1,Gcompq1)),plotat.short,
         pch=c(rbind(rep("[",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  points(c(rbind(thetamatGq.inner2,Gcompq1)),plotat.short,
         pch=c(rbind(rep("]",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  counter = 0
  for (s in 1:length(plotat.short)){
    segments(y0=plotat.short[s],x0=c(rbind(thetamatGq.inner1,Gcompq1))[s],
             y1=plotat.short[s],x1=c(rbind(thetamatGq.inner2,Gcompq2))[s],
             col="black",
             lwd=2)
  }
  
  # row separation lines
  abline(h=seq(1.5,J.full-.5,by=1)[ind.gs],lty=2,lwd=0.5,col="blue")
  abline(h=seq(1.5,J.full-.5,by=1)[ind.gs]-1,lty=2,lwd=0.5,col="blue")
  
  
  counter = 0
  for (s in 1:length(plotat.short)){
    segments(y0=plotat.short[s],x0=c(rbind(thetamatGq1,Gcompq1))[s],
             y1=plotat.short[s],x1=c(rbind(thetamatGq2,Gcompq2))[s],col="black",
             lty=ifelse((s-1)%%3<2,1,1))
    if ((s-1)%%3>=1){
      counter=counter+1
      text(c(Gcomp)[counter],plotat.short[s]+0.125,paste0(round(100*c(Gcomp),1)[counter],"%"),srt=srtval,cex=cexval)
    }
  }
  
  
  for (s in 1:J.GS){
    
    text(thetameanG[s],plotat.short[3*s-2]+.125,paste(round(100*thetameanG[s],2),"%"))
    
    # put prior shapes on gold sensitivity
    tmp = rbeta(10000,alphaG[s],betaG[s])
    boxplot(tmp,at = ind.gs[s]-0.45, boxwex=1/8 ,col="gray", 
            add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
    
  }
  
  mtext(expression(underline("SS")),line=1,cex=1.8)
  
  par(op)
  
  
  ######################################################
  ##subplot 3: etiology bars with prior shapes
  ######################################################
  top=0.4#max(pq2BG)
  dotcolor = "black"
  #op <- par(mar=c(5.1,6,4.1,1.1))
  op <- par(mar=c(5.1,0,4.1,9))
  plot(c(pmeanBG.ord),1:(J.full),
       yaxt="n",#xaxt="n",
       xlim=c(0,top),ylim=c(0.5,J.full+0.5),col=c("black"),
       ylab="",xlab="probability",
       pch=c(20),cex=2)
  axis(4,at=1:J.full,labels=pathogens.ord,las=2,cex.axis=1.5)
  abline(h=seq(1.5,J.full-.5,by=1),lty=2,lwd=0.5,col="blue")
  #draw axis within plot:
  for (s in 1:(J.full-1)){
    axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=1,#labels=rep("",length(seq(0,1,by=.2))),
         pos = seq(1.5,J.full-.5,by=1)[s], cex.axis = 0.8,lty=2,col="blue")
    # axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=0,#labels=rep("",length(seq(0,1,by=.2))),
    #      pos = seq(1.5,J.full-.5,by=1)[s]+0.3, cex.axis = 0.8,lty=2,col="blue")
  }
  points(c(pq1BG),1:(J.full),pch="|",cex=1)
  points(c(pq2BG),1:(J.full),pch="|",cex=1)
  points(c(p.innerq1BG),1:(J.full),pch="[",cex=1)
  points(c(p.innerq2BG),1:(J.full),pch="]",cex=1)
  
  mtext(expression(underline(hat(pi))),line=1,cex=1.8)
  #mtext(c(expression(bold("--")),":prior","-",":posterior"),col=c("gray","black","black","black"),
  #      adj=c(0,0.1,0.3,0.4),line=.8,cex=.8,lwd=2)
  legend("topright",c("prior","posterior"),lty=c(2,1),col=c("gray","black"),
         lwd = 4,horiz=TRUE,cex=1.5,bty="n")
  pgrid = seq(0,1,by=0.01)
  
  for (s in 1:(J.full)){
    segments(y0=s,x0=c(pq1BG)[s],y1=s,x1=c(pq2BG)[s],col=dotcolor)
    segments(y0=s,x0=c(p.innerq1BG)[s],y1=s,x1=c(p.innerq2BG)[s],col=dotcolor,lwd=2)
    #text(pmeanBG[s],s+0.25,paste0(round(100*c(pmeanBG),2)[s],"%"),srt=srtval,cex=cexval)
    text(.35,s+0.25,paste0("=",paste0(round(100*c(pmeanBG.ord),1)[s],"%")),srt=srtval,cex=2)
    text(.3,s+0.25,bquote(hat(pi)[.(ord[s])]),srt=srtval,cex=2)
    #text(top-0.05,s,pathogens[s],cex=.8,srt=srtval)
    tmp.density = dbeta(pgrid,alphaE[ord[s]],sum(alphaE[-ord[s]]))
    points(pgrid,tmp.density/(3*max(tmp.density))+s-0.45,type="l",col="gray",lwd=4,lty=2)
    ##posterior density
    tmp.post.density = density(pmatBG.ord[,s],from=0,to=1)
    tmp.x = tmp.post.density$x
    tmp.y = tmp.post.density$y
    points(tmp.x,tmp.y/(3*max(tmp.y))+s-0.45,col="black",lwd=4,type="l")
    
  }
  
  par(op)
  
  
  #######################################################
  #######################################################
  #######################################################
  #######################################################
  #######################################################
  
}



infograph_popPie_ord_sameTPR = function(directoryName,
                                alphaE=alpha,
                                alphaB = 1,
                                betaB = 1,
                                alphaG = rep(1,JGS),
                                betaG  = rep(1,JGS),
                                Bdat=MBS,Gdat=MGS,
                                ncase = Nd,J.GS = JGS,
                                nppcr.no.realGS.fit=FALSE,
                                nppcr.no.realGS.fit.dir=prior.folder,#"\\nppcr.no.realGS.fit",
                                nppcr.bcx.fit=TRUE,
                                nppcr.bcx.fit.dir=prior.folder,#"\\nppcr.bcx.fit",
                                pEtiName="pEti",pathogens = Pathogen,
                                top3=0.25){
  
    
#       directoryName= fullname
#       alphaE = alpha
#       alphaB = alphaB
#       betaB = betaB
#       alphaG = alphaG
#       betaG  = betaG
#       nppcr.bcx.fit=TRUE
#       nppcr.no.realGS.fit = FALSE
#       pathogens=Pathogen
#       ncase = Nd
#       J.GS=JGS
#       pEtiName="pEti"
#       nppcr.bcx.fit.dir=prior.folder
#       Bdat=MBS;Gdat=MGS
    
  require(R2WinBUGS)
  cexval=1
  srtval=0
  J.full = length(pathogens)
  
  layout(matrix(c(1,2,3),1,3,byrow = TRUE), 
         widths=c(3,2,3),heights=c(8))
  #layout.show(n=3)
  
  curr.folder = paste0(directoryName,nppcr.bcx.fit.dir)
  simsMat = read.coda(paste(curr.folder,"coda1.txt",sep="\\"),
                      paste(curr.folder,"codaIndex.txt",sep="\\"),quiet=TRUE)  
  
  res = matrix(NA,nrow=nrow(simsMat),ncol=J.full)
  SubVarName = rep(NA,J.full)
  for (j in 1:J.full){
    SubVarName[j] = paste(pEtiName,"[",j,"]",sep="")
  }
  pmatBG = simsMat[,SubVarName]
  pmeanBG = colMeans(pmatBG)
  pmeanBG0 = colMeans(pmatBG)
  
  ord = order(pmeanBG)
  #### from now on permute!
  pathogens.ord = pathogens[ord]
  pmeanBG.ord = pmeanBG[ord]
  pmatBG.ord = pmatBG[,ord]
  
  ## quantiles for etiology: outer is 97.5% CI, inner is 50% CI
  pq1BG   = apply(pmatBG,2,quantile,probs=0.025)[ord]
  pq2BG   = apply(pmatBG,2,quantile,probs=0.975)[ord]
  p.innerq1BG   = apply(pmatBG,2,quantile,probs=0.25)[ord]
  p.innerq2BG   = apply(pmatBG,2,quantile,probs=0.75)[ord]
  
  ######################################################
  ##subplot 1: case-control rate comparisons
  ######################################################
  require(binom)
  nctrl = nrow(Bdat)-ncase
  tmp.case = binom.confint(colSums(Bdat[1:ncase,]), ncase, conf.level = 0.95, methods = "ac")[ord,]
  tmp.ctrl = binom.confint(colSums(Bdat[-(1:ncase),]), nctrl, conf.level = 0.95, methods = "ac")[ord,]
  
  ## case and control positive rate, lower and upper limit
  Bcomp = rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5))
  Bcompq1 = rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")])
  Bcompq2 = rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")])
  
  ## posterior distribution of TPR
  thetamatBG = matrix(simsMat[,grep("thetaBS",colnames(simsMat))],nrow=nrow(pmatBG),
                      ncol=J.full)
  
  thetameanBG = colMeans(thetamatBG)
  
  ## posterior distribution of FPR
  psimatBG = (simsMat[,grep("psiBS",colnames(simsMat))])[,ord]
  psimeanBG = colMeans(psimatBG)  
  
  ## model fitted postive rate for each pathogen
  fittedmeanBGcase = sapply(1:J.full, 
                            function(s) mean(pmatBG.ord[,s]*thetamatBG[,s]+(1-pmatBG.ord[,s])*psimatBG[,s]))
  fittedmeanBGcontrol = psimeanBG
  
  top2 = 1.3
  op<-par(mar=c(5.1,4.1,4.1,0))
  plotat = seq(0.5,J.full+0.5,by=1/4)[-(c(1,(1:J.full)*4+1))]
  plot(c(rbind(fittedmeanBGcase,Bcomp)),plotat,yaxt="n",xlim=c(0,top2),
       ylim=c(0.5,J.full+.5),xaxt="n",
       ylab="",xlab="probability",
       pch = c(rbind(rep(2,J.full),rep(20,J.full),rep(20,J.full))),
       col=c(rbind(rbind(rep(1,J.full),rep("dodgerblue2",J.full),rep("dodgerblue2",J.full)))),
       cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
  axis(2,at = plotat,labels=rep(c("","case","ctrl"),J.full),las=2)
  axis(1,at = c(0,0.2,0.4,0.6,0.8,1),labels= c(0,0.2,0.4,0.6,0.8,1),las=1)
  points(c(rbind(thetameanBG,Bcompq2)),plotat,
         pch=c(rbind(rep("+",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(2,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  points(c(rbind(psimeanBG,Bcompq1)),plotat,
         pch=c(rbind(rep("*",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(2,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  abline(h=seq(1.5,J.full-.5,by=1),lty=2,lwd=0.5,col="blue")
  abline(v=1,lty=2,lwd=.5)
  
  ## conditional odds ratios
  COR = function(brs.data,nd,pathogens){
    y = rep(c(1,0),times=c(nd,nrow(brs.data)-nd))
    #X = matrix(NA,nrow=nrow(brs.data),ncol=length(pathogens))
    x.nm   = paste(pathogens,"NPPCR",sep="_")
    #colnames(X) = x.nm
    #for (j in 1:length(pathogens)){
    #  X[,j] = brs.data[,x.nm[j]]
    #}
    dat.reg = as.data.frame(cbind(y,brs.data))
    form = as.formula(paste0("y~",paste(x.nm,collapse="+")))
    fit = glm(form,data=dat.reg,family=binomial)
    res = cbind(exp(confint(fit)),exp(fit$coef))[-1,]
  }
  
  tmp = COR(MBS,ncase,pathogens.ord)
  for (s in (1:J.full)){
    L=round(tmp[s,1],1)
    C=round(tmp[s,3],1)
    R=round(tmp[s,2],1)
    text(top2-0.12,s+1/(2*Jfull),C,cex=1.5)
    text(top2-0.12,s-.2,paste(c(L,"   ",R),collapse=" "),cex=1.2)
  }
  legend("topright","conditional OR",bty="n")
  
  #   ## marginal odds ratios
  #   for (s in (1:J.full)){
  #     if (any(Bcomp[,s]==0)){
  #       cat("Pathogen ",pathogens[ord[s]]," have zero cells cross-classifying Y and MBS.")
  #       text(top2*.9,s,NA)
  #     } else{
  #       
  #       tmp = LOR(Bdat[,ord[s]],rep(c(1,0),times=c(ncase,nrow(Bdat)-ncase)))
  #       L=round(tmp$OR[2],1)
  #       C=round(tmp$OR[1],1)
  #       R=round(tmp$OR[3],1)
  #       text(top2-0.12-0.5,s+1/(2*Jfull),C,cex=1.5)
  #       text(top2-0.12-0.5,s-.2,paste(c(L,"   ",R),collapse=" "),cex=1.2)
  #     }
  #   }
  
  counter = 0
  for (s in 1:(3*J.full)){
    segments(y0=plotat[s],x0=c(rbind(psimeanBG,Bcompq1))[s],
             y1=plotat[s],x1=c(rbind(thetameanBG,Bcompq2))[s],col="black",
             lty=ifelse((s-1)%%3<1,4,1))
    if ((s-1)%%3>=1){
      counter=counter+1
      tmp.hpos = ifelse(c(Bcompq2)[counter]+0.15>0.95,c(Bcompq1)[counter]-0.2,c(Bcompq2)[counter]+0.15 )
      text(tmp.hpos,plotat[s],paste0(round(100*c(Bcomp),1)[counter],"%"),srt=srtval,cex=cexval)
    }
  }
  
  for (s in 1:(J.full)){
    segments(y0=plotat[3*s-1],x0=c(rbind(fittedmeanBGcase,Bcomp))[3*s-1],
             y1=plotat[3*s],x1=c(rbind(fittedmeanBGcase,Bcomp))[3*s],col="black")
  }
  # put prior shapes on bronze sensitivity
  for (s in 1:J.full){
    tmp = rbeta(10000,alphaB,betaB)
    boxplot(tmp,at = s-0.45, boxwex=1/10 , col="gray",add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
    tmp.post = as.matrix(thetamatBG)[,s]
    boxplot(tmp.post,at = s-0.35,boxwex=1/10,add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
  }
  axis(2,at=(1:J.full)-.45,label=rep("",J.full),las=2,cex.axis=.5)
  axis(2,at=(1:J.full)-.35,label=rep("",J.full),las=2,cex.axis=.5)
  
  mtext(expression(underline("BrS")),line=1,cex=1.8)
  par(op)
  
  ######################################################
  ##subplot 2: silver standard data comparisons
  ######################################################
  
  Gcomplete.index = which(rowSums(is.na(Gdat[,1:J.GS]))==0)
  cat("No. of Complete Silver-standard data on ",pathogens[1:J.GS],":",
      length(Gcomplete.index),"\n")
  cat(ncase-length(Gcomplete.index)," cases deleted for plotting third subgraph.","\n")
  Gdat = Gdat[Gcomplete.index,]
  ord.gs = order(pmeanBG0[1:J.GS])
  ind.gs = rep(NA,J.GS) # tells where the the gs row should go
  for (j in 1:J.GS){
    ind.gs[j] = which(ord==j)
  }
  
  tmpG.case = binom.confint(colSums(Gdat[,1:J.GS]), nrow(Gdat), conf.level = 0.95, methods = "ac")
  Gcomp = rbind(round(tmpG.case$mean,5),rep(NA,J.GS))
  Gcompq1 = rbind(tmpG.case[,c("lower")],rep(NA,J.GS))
  Gcompq2 = rbind(tmpG.case[,c("upper")],rep(NA,J.GS))
  
  thetamatG = (simsMat[,grep("thetaGS",colnames(simsMat))])########### stopped here
  thetameanG = colMeans(thetamatG)
  
  thetamatGq1=apply(thetamatG,2,quantile,0.025)
  thetamatGq2=apply(thetamatG,2,quantile,0.975)
  
  fittedmeanGpos = sapply(1:J.GS, function(s) mean(pmatBG.ord[,ind.gs[s]]*thetamatG[,s]))
  
  
  #top3 = .25
  par(mar=c(5.1,0,4.1,0))
  #op<-par(mar=c(5.1,0,4.1,2.1),xpd=TRUE)
  
  plotat = seq(0.5,J.full+0.5,by=1/4)[-(c(1,(1:J.full)*4+1))]
  #plotat.short = plotat[1:length(c(rbind(thetameanG,Gcomp)))]
  plotat.calc = function(j) {c(3*j-2,3*j-1,3*j)}
  plotat.short = rep(NA,J.GS*3)
  for (j in 1:J.GS){
    plotat.short[c(3*j-2,3*j-1,3*j)] = plotat[plotat.calc(ind.gs[j])]
  }
  
  plot(c(rbind(thetameanG,Gcomp)),plotat.short,yaxt="n",xlim=c(0,top3),
       ylim=c(0.5,J.full+.5),#xaxt="n",
       ylab="",xlab="probability",
       pch = c(rbind(rep(20,J.full),rep(20,J.full),rep(20,J.full))),
       col=c(rbind(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full)))),
       cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
  #axis(1,at = seq(0,top3,len=10),labels= seq(0,top3,len=10),las=1)
  points(c(rbind(fittedmeanGpos,matrix("",nrow=2,ncol=J.GS))),plotat.short,yaxt="n",xlim=c(0,top3),
         ylim=c(0.5,J.full+.5),xaxt="n",
         ylab="",#xlab="Gold Positive Rate",
         pch = c(rbind(rep(2,J.full),rep(NA,J.full),rep(NA,J.full))),
         col=c(rbind(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full)))),
         cex = c(rbind(rep(1,J.full),rep(2,J.full),rep(2,J.full))))
  
  #axis(2,at = plotat.short,labels=rep(c("TPR-G","case","ctrl"),J.GS),las=2)
  points(c(rbind(thetamatGq2,Gcompq2)),plotat.short,
         pch=c(rbind(rep("|",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  points(c(rbind(thetamatGq1,Gcompq1)),plotat.short,
         pch=c(rbind(rep("|",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  
  #inner 25%-75%
  thetamatGq.inner1=apply(thetamatG,2,quantile,0.25)
  thetamatGq.inner2=apply(thetamatG,2,quantile,0.75)
  points(c(rbind(thetamatGq.inner1,Gcompq1)),plotat.short,
         pch=c(rbind(rep("[",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  points(c(rbind(thetamatGq.inner2,Gcompq1)),plotat.short,
         pch=c(rbind(rep("]",J.full),rep("|",J.full),rep("|",J.full))),
         cex=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))),
         col=c(rbind(rep(1,J.full),rep(1,J.full),rep(1,J.full))))
  counter = 0
  for (s in 1:length(plotat.short)){
    segments(y0=plotat.short[s],x0=c(rbind(thetamatGq.inner1,Gcompq1))[s],
             y1=plotat.short[s],x1=c(rbind(thetamatGq.inner2,Gcompq2))[s],
             col="black",
             lwd=2)
  }
  
  # row separation lines
  abline(h=seq(1.5,J.full-.5,by=1)[ind.gs],lty=2,lwd=0.5,col="blue")
  abline(h=seq(1.5,J.full-.5,by=1)[ind.gs]-1,lty=2,lwd=0.5,col="blue")
  
  
  counter = 0
  for (s in 1:length(plotat.short)){
    segments(y0=plotat.short[s],x0=c(rbind(thetamatGq1,Gcompq1))[s],
             y1=plotat.short[s],x1=c(rbind(thetamatGq2,Gcompq2))[s],col="black",
             lty=ifelse((s-1)%%3<2,1,1))
    if ((s-1)%%3>=1){
      counter=counter+1
      text(c(Gcomp)[counter],plotat.short[s]+0.125,paste0(round(100*c(Gcomp),1)[counter],"%"),srt=srtval,cex=cexval)
    }
  }
  
  
  for (s in 1:J.GS){
    
    text(thetameanG[s],plotat.short[3*s-2]+.125,paste(round(100*thetameanG[s],2),"%"))
    
    # put prior shapes on gold sensitivity
    tmp = rbeta(10000,alphaG[s],betaG[s])
    boxplot(tmp,at = ind.gs[s]-0.45, boxwex=1/8 ,col="gray", 
            add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
    
  }
  
  mtext(expression(underline("SS")),line=1,cex=1.8)
  
  par(op)
  
  
  ######################################################
  ##subplot 3: etiology bars with prior shapes
  ######################################################
  top=0.4#max(pq2BG)
  dotcolor = "black"
  #op <- par(mar=c(5.1,6,4.1,1.1))
  op <- par(mar=c(5.1,0,4.1,9))
  plot(c(pmeanBG.ord),1:(J.full),
       yaxt="n",#xaxt="n",
       xlim=c(0,top),ylim=c(0.5,J.full+0.5),col=c("black"),
       ylab="",xlab="probability",
       pch=c(20),cex=2)
  axis(4,at=1:J.full,labels=pathogens.ord,las=2,cex.axis=1.5)
  abline(h=seq(1.5,J.full-.5,by=1),lty=2,lwd=0.5,col="blue")
  #draw axis within plot:
  for (s in 1:(J.full-1)){
    axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=1,#labels=rep("",length(seq(0,1,by=.2))),
         pos = seq(1.5,J.full-.5,by=1)[s], cex.axis = 0.8,lty=2,col="blue")
    # axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=0,#labels=rep("",length(seq(0,1,by=.2))),
    #      pos = seq(1.5,J.full-.5,by=1)[s]+0.3, cex.axis = 0.8,lty=2,col="blue")
  }
  points(c(pq1BG),1:(J.full),pch="|",cex=1)
  points(c(pq2BG),1:(J.full),pch="|",cex=1)
  points(c(p.innerq1BG),1:(J.full),pch="[",cex=1)
  points(c(p.innerq2BG),1:(J.full),pch="]",cex=1)
  
  mtext(expression(underline(hat(pi))),line=1,cex=1.8)
  #mtext(c(expression(bold("--")),":prior","-",":posterior"),col=c("gray","black","black","black"),
  #      adj=c(0,0.1,0.3,0.4),line=.8,cex=.8,lwd=2)
  legend("topright",c("prior","posterior"),lty=c(2,1),col=c("gray","black"),
         lwd = 4,horiz=TRUE,cex=1.5,bty="n")
  pgrid = seq(0,1,by=0.01)
  
  for (s in 1:(J.full)){
    segments(y0=s,x0=c(pq1BG)[s],y1=s,x1=c(pq2BG)[s],col=dotcolor)
    segments(y0=s,x0=c(p.innerq1BG)[s],y1=s,x1=c(p.innerq2BG)[s],col=dotcolor,lwd=2)
    #text(pmeanBG[s],s+0.25,paste0(round(100*c(pmeanBG),2)[s],"%"),srt=srtval,cex=cexval)
    text(.35,s+0.25,paste0("=",paste0(round(100*c(pmeanBG.ord),1)[s],"%")),srt=srtval,cex=2)
    text(.3,s+0.25,bquote(hat(pi)[.(ord[s])]),srt=srtval,cex=2)
    #text(top-0.05,s,pathogens[s],cex=.8,srt=srtval)
    tmp.density = dbeta(pgrid,alphaE[ord[s]],sum(alphaE[-ord[s]]))
    points(pgrid,tmp.density/(3*max(tmp.density))+s-0.45,type="l",col="gray",lwd=4,lty=2)
    ##posterior density
    tmp.post.density = density(pmatBG.ord[,s],from=0,to=1)
    tmp.x = tmp.post.density$x
    tmp.y = tmp.post.density$y
    points(tmp.x,tmp.y/(3*max(tmp.y))+s-0.45,col="black",lwd=4,type="l")
    
  }
  
  par(op)
  
  
  #######################################################
  #######################################################
  #######################################################
  #######################################################
  #######################################################
  
}

