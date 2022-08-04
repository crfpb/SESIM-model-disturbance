#############################################
#             output data frames            #
#  SESIM model - Journal of Animal Ecology  #
#############################################

# experimental design
experiment <- data.frame(treatment=rep(1:(length(rate)*length(eP)*length(eH)*length(scenario)*length(disturb))),
                         dispersal=NA, environment=NA, heterog=NA, interact=NA, disturbance=NA)

# sampled abundances
X_save <- array(data=NA, dim=c(numCom, nSp, length(sampleV)))

# community level results
results_comm <- data.frame(treatment=as.numeric(), time=as.numeric(),
                           replicate=as.numeric(), patch=as.numeric(),
                           richness=as.numeric())

# metacommunity level results
results_meta <- data.frame(treatment=as.numeric(), time=as.numeric(), replicate=as.numeric(), 
                           richness=as.numeric())

