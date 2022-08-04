#############################################
#             input parameters              #
#  SESIM model - Journal of Animal Ecology  #
#############################################

# spatial network parameters #
# number of patches
numCom <- 5
# connectivity matrix
d <- matrix(c(1.0), nrow=numCom, ncol=numCom)
diag(d) <- 0
# dispersal rates
#rate <- a <- c(0, 0.0001, 0.001, 0.01, 0.1, 1)
rate <- a <- 0.01
# dispersal kernel
k <- 1
# dispersal ability
d_exp <- exp(-k*d) - diag(nrow(d))
dispersal_m <- apply(d_exp, 1, function(x) x/sum(x))
dispersal_m[is.na(dispersal_m)] <- 0

# species attributes parameters #
# (presented in the following order: Pa Pb S C T V E B)
# species-specific dispersal ability
disp <- c(0.64, 0.62, 0.53, 0.90, 0.62, 0.44, 0.71, 1.00)
# species-specific growth rate
r.high <- c(1.32, 0.92, 1.64,  0.11, 2.63,  0.15, 1.78, 1.31)/10
r.base <- c(1.47, 0.38, 1.64,  0.70, 1.61,  0.42, 1.49, 1.32)/10
# species-specific cell mass
mass <- c(4.3e-8, 1.96e-7, 3.76e-6, 1.52e-8, 4.77e-9, 6.9e-8, 8.05e-8, 2.27e-7)

# food web parameters #
# number of species
nSp <- 8
# species interactions
#scenario <- c("species","competitive")
interact <- scenario <- "species"
# constant
o = 0.0007
# positive interactions strenghts
source("./pos_BB.R")
# negative interactions strenghts
source("./neg_BB.R")
# for species interactions scenario:
spp_matrices <- function(interact){
  if(interact=="species"){
    # positive interactions matrix
    BB_pos <- pos <- fullFW_pos
    # competitive interactions matrix
    comp_BB <- comp <- matrix(c(-0.1, rep(-0.05, nSp)), nSp, nSp)/10
    # positive interactions matrix
    BB_neg <- neg <- fullFW_neg + comp_BB
  } 
  # for purely competitive scenario:
  if(interact=="competitive"){
    # positive interactions matrix
    BB_pos <- pos <- fullFW_pos
    BB_pos[BB_pos > 0] <- 0
    # competitive interactions matrix
    comp_BB <- neg <- comp <- matrix(c(-0.1, rep(-0.05, nSp)), nSp, nSp)/10
    # positive interactions matrix
    BB_neg <- comp_BB 
  }
  return(array(c(BB_pos, BB_neg), dim=c(8,8,2))) 
}

# environment parameters #
# number of environmental variables
nEnv <- 1
# amplitude of environmental fluctuations
eA <- 1
# strength of environmental effect
enveff <- 1
# period of environmental fluctuation
eP <- 50
# environmental heterogeneity
eH <- c(0, 1)

# disturbance parameters #
# disturbance gradient
disturb <- pdistr <- c(0, 20, 40, 60, 80, 100)
# disturbance time
disturb_time <- 500


