
#panels 2 and 3 are in this script
source('SST_MPA_Fishing.R')

#panel 1 in this script
source('SST_MPA_Fishing_Relative_to_Baseline.R')


figure <-ggarrange(p3 +rremove("xlab"), p1 +rremove("xlab"), p2+rremove("xlab"), nrow=1, ncol=3, common.legend = TRUE) 

annotate_figure(
  figure,
  bottom = text_grob("Fishing effort", size = 20, vjust = 0)
)
