
#panels 2 and 3 are in this script
source('SST_MPA_Fishing.R')

#panel 1 in this script
source('SST_MPA_Fishing_June2024.R')


ggarrange(p3, p1, p2, nrow=1, ncol=3, common.legend = TRUE) 
