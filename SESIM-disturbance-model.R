#############################################
#               model dynamics              #
#  SESIM model - Journal of Animal Ecology  #
#############################################

# input parameters
# note: ignore warning
source("./input.R")

# simulation parameters
source("./simulation.R")

# output parameters
source("./output.R")

# model #
# experimental treatment count
treatment_id = 1
# loop count
loop_id = 1
# for loop of experimental design
for(interact in scenario){
  for(a in rate){
    for(eHeterog in eH){
      for(pdistr in disturb){
        for(ePeriod in eP){
          for(r in 1:reps){
            
            # sampling count          
            sampling = 1
            # set seed          
            set.seed(xseed[r])
            
            # initial abundances
            high <- base <- rep(100, nSp)
            X <- cbind(high, base, high, high, base)
            Xd <- t(X)
            
            # growth rate matrix
            C <- cbind(r.high, r.base, r.high, r.high, r.base)
            
            # species environmental optima
            Env_Opt <- matrix(runif(nSp, 0, 1), nSp, nEnv)
            
            # environmental heterogeneity states
            if(eHeterog==0){
              envt.v <- c(1, 1, 1, 1, 1)
              env1 <- c(1, 1, 1, 1, 1)
              env2 <- c(0, 0, 0, 0, 0)
            }
            if(eHeterog==1){
              envt.v <- c(1, 0, 1, 0, 1)
              env1 <- c(1, 0, 1, 0, 1)
              env2 <- c(0, 1, 0, 1, 0)
            }
            
            # environmental fluctuation
            if(ePeriod==0) fluc.time <- 0
            if(ePeriod!=0) fluc.time <- seq(0, Tmax, by=ePeriod)
            
            # interaction strengths
            # positive interactions
            pos <- BB_pos <- spp_matrices(interact)[,,1]
            # negative interactions
            neg <- BB_neg <- spp_matrices(interact)[,,2]
            
            # interaction strengths matrices
            # positive interactions
            for(j in 1:ncol(BB_pos)){
              for(i in 1:nrow(BB_pos)){
                set.seed(xseed[r])
                pos[i,j] <- runif(1, min=min(BB_pos[i,j]/2, BB_pos[i,j]), max=max(BB_pos[i,j]/2, BB_pos[i,j]))
              }
            }
            pos[pos > 0] <- 0
            # negative interactions
            for(j in 1:ncol(BB_neg)){
              for(i in 1:nrow(BB_neg)){
                set.seed(xseed[r])
                neg[i,j] <- runif(1, min=min(BB_neg[i,j]/2, BB_neg[i,j]), max=max(BB_neg[i,j]/2, BB_neg[i,j]))
              }
            }
            
            # for loop of species dynamics
            for(l in 1:Tmax){
              
              # environment
              Env <- matrix(rep(envt.v, each=nSp), nSp, numCom)
              
              # environmental match
              enviro <- t(apply(1-abs(Env_Opt[,1] - Env), 1, function(x) x/enveff))
              
              # negative interactions
              neg_interactions <- neg %*% (enviro*X)    
              # positive interactions
              pos_interactions <- (pos %*% X) * enviro  
              # total interactions
              interactions <- pos_interactions + neg_interactions
              
              # growth
              growth <- C * enviro
              
              # migration
              migrants <- apply(X, 2, function(x) x*(a*disp))
              
              # immigration
              immigrants <- matrix(NA, nSp, numCom)
              for(i in 1:nSp){
                immigrants[i,] <- dispersal_m%*%X[i,]*(a*disp[i])
              }
              
              # generalized lotka-volterra equation
              Xt <- X * exp(growth + interactions) + immigrants - migrants
              
              # extinctions
              Xt[Xt < a] <- 0
              
              # disturbance
              if(any(l==disturb_time) & pdistr!=0){
                if(pdistr==20) Xt[,c(1)] <- 0
                if(pdistr==40) Xt[,c(1,4)] <- 0
                if(pdistr==60) Xt[,c(1,3,4)] <- 0
                if(pdistr==80) Xt[,c(1,3,4,5)] <- 0
                if(pdistr==100) Xt[,c(1,2,3,4,5)] <- 0
              }
              
              # abundance matrix
              X <- Xt
              Xd <- t(X)
              
              # sample
              if(any(sampleV==l) && l<=Tmax){
                X_save[,,sampling] <- Xd
                sampling <- sampling+1
              }
              
              # environmental fluctuation
              if(ePeriod!=0 && any(l==fluc.time)){
                if(envt.v[1]==1){
                  envt.v <- env2
                }else{
                  envt.v <- env1
                }
              }
              
            } # end of for loop for species dynamics
            
            # output data frames
            
            # experimental design
            experiment$treatment[treatment_id] <- treatment_id
            experiment$environment[treatment_id] <- ePeriod
            experiment$dispersal[treatment_id] <- a
            experiment$heterog[treatment_id] <- eHeterog
            experiment$sppinteract[treatment_id] <- interact
            experiment$disturbance[treatment_id] <- pdistr
            
            # calculate community results
            comm_res <- data.frame(treatment=rep(treatment_id, 5*Tmax), time=rep(1:Tmax, each=5),
                                   replicate=rep(r, 5*Tmax), patch=c(1:5),
                                   richness=as.vector(apply(X_save, MARGIN=3, FUN=specnumber)))
            # combine results
            results_comm <- rbind(results_comm, comm_res)
            
            # metacommunity level results
            meta_res <- data.frame(treatment=rep(treatment_id, Tmax), time=rep(1:Tmax), replicate=rep(r, Tmax),
                                   richness=apply(apply(X_save, c(2,3), sum), 2, specnumber))
            # combine results
            results_meta <- rbind(results_meta, meta_res)
            
            # loop count
            loop_id <- loop_id+1
            
          } # end of for loop for replicates
          
          # count treatment          
          treatment_id  <- treatment_id+1
          
        }
      }
    } 
  }
}

# end

# combine experiment and results using "treatment"
community_results <- merge(experiment, results_comm, by="treatment")
metacommunity_results <- merge(experiment, results_meta, by="treatment")

# save results
#write.csv(community_results, file=paste(a, scenario, "comm", sep="-"))
#write.csv(metacommunity_results, file=paste(a, scenario, "meta", sep="-"))
