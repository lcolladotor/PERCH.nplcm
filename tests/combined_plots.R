library(R2WinBUGS)
library(coda)
library(nplcm)

##Kenya:
rm(list=ls())
parent_folder = "C:/Users/Administrator/Dropbox/PERCH Core/Data Analysis/PQ Analysis/Meeting Materials/2014 10 03/20141002_01KEN"
DIR_NPLCM = paste0(parent_folder, "/20141002_01KEN_CD")
DIR_PLCM  = paste0(parent_folder, "/20141002_01KEN_CI")

pdf(paste0(parent_folder,"/combined_individual_diagnosis.pdf"),
    width=16,height=16)
par(mfrow=c(4,4))
combined_visualization(DIR_NPLCM,DIR_PLCM,npat = 16)
dev.off()


# rm(list=ls())
# parent_folder = "C:/Users/Administrator/Dropbox/PERCH Core/Data Analysis/PQ Analysis/Meeting Materials/2014 10 03/20141002_01KEN"
# DIR_NPLCM = paste0(parent_folder, "/20141002_01KEN_CD")

# pdf(paste0(DIR_NPLCM,"/individual_diagnosis.pdf"),
    # width=16,height=16)
# par(mfrow=c(4,4))
# nplcm_plot_individual_diagnosis(DIR_NPLCM,npat = 16)
# dev.off()



## Gambia:
rm(list=ls())
parent_folder = "C:/Users/Administrator/Dropbox/PERCH Core/Data Analysis/PQ Analysis/Meeting Materials/2014 10 03/20141002_02GAM"
DIR_NPLCM = paste0(parent_folder, "/20141002_02GAM_CD")
DIR_PLCM  = paste0(parent_folder, "/20141002_02GAM_CI")

pdf(paste0(parent_folder,"/combined_individual_diagnosis.pdf"),
    width=16,height=16)
par(mfrow=c(4,4))
combined_visualization(DIR_NPLCM,DIR_PLCM,npat = 16)
dev.off()


library(R2WinBUGS)
library(coda)
library(nplcm)

rm(list=ls())
parent_folder = "C:/Users/Administrator/Dropbox/PERCH Core/Data Analysis/PQ Analysis/Meeting Materials/2014 10 03/20141002_02GAM"
DIR_NPLCM = paste0(parent_folder, "/20141002_02GAM_CD")

nplcm_plot_three_panel(DIR_NPLCM,Mobs,Y,X,model_options)

pdf(paste0(DIR_NPLCM,"/individual_diagnosis.pdf"),
    width=16,height=16)
par(mfrow=c(4,4))
nplcm_plot_individual_diagnosis(DIR_NPLCM,npat = 16)
dev.off()
