###############
## LOADING FUNCTIONS AND PACKAGES ##
###############
source("Farming_plant_cooperation_R_converted_analytical_model.R") # file path needs to be adjusted
source("Farming_plant_cooperation_fecundity.R") # file path needs to be adjusted
library(sf)
library(Deriv)

###############
## DEFINITION OF THE SIMULATION FUNCTION ##
###############


##'@param xopt = optimal plant height for individual seed production
##'@param a = light interception ability = angle at the summit of the triangles
##'@param nG = number of generations to be simulated
##'@param nF = number of fields to be simulated
##'@param nP = number of individuals simulated in each field
##'@param mut_rate = mutation rate
##'@param mut_effect = mutation effect (absolute value)
##'@param nS = number of seeds sampled from generation "g" to repopulate a field in generation "g+1"
##'@param m = seed migration rate
##'@param alpha = weight given to the among-field competition in the selective process. For positive values, founders coming from the pooled harvest of all fields have a much higher probability to come from the most productive field than it would if field productivity alone was condidered.For grp_selec_intentisty = 0, only field productivity is considered when chosing founders in the pooled harvest of all fields
##'@param init = initial phenotypic value for the monomorphic meta-population
##'@param d = physical distance between plants
##'@param seed = random number generator parameter, to ensure reproductibility


simul <- function(xopt, a, nG, nF, nP, mut_rate, mut_effect, nS, m, alpha, init, d, seed) {
  
  #########################################################################################
  ##### COMPUTING RELATEDNESS AND EXPECTED ESS #######
  
  # relatedness #
  gamma <- (1-mut_rate)^2
  
  if(alpha==0) {
    phi <- 0
  } else {
    phi <- 1
  }
  
  # ESS #
  dx <- Deriv(Ws,"x")
  dy <- Deriv(Ws,"y")
  
  equ <- function (z) {
    dx(z, z, z) + relatedness(nS, m, gamma, phi)*dy(z, z, z)
  }
  
  ess <- c()
  for (j in c(10E-3:200)) {
    tryCatch({root <- uniroot(equ, c(j,j+1))$root
    ess <- c(ess,root)}, error=function(e) paste("Oops"))
  }
  
  
  #########################################################################################
  ##### SIMULATION OF THE FIRST GENERATION #######
  farm <- data.frame(IDinit=NULL,IDfield=NULL, IDgenotype=NULL, height=NULL) # empty object to store all fields of the first generation
  
  for (s in 1:nF) { # each field is populated consecutively
    
    set.seed(s) # setting the seed, a parameter for random number sampling. The seed enables to reproduce a given random sampling by setting the same seed parameter before using "sample()" function
    
    field <- data.frame(IDinit=init, IDfield=rep(s,nP),IDgenotype=c(1:nP), height= phenotypes[[init]][sample(length(phenotypes[[init]]), nP, replace=T)]) # simulation of the field : phenotypes are sampled following a uniform probability density within the range of phenotypic values "i".
    
    farm <- rbind(farm, field) # combining the different fields with the different initial phenotypes into a single object "farm"
    remove(field) # removing temporary objects
    
  }
  
  
  rec_init <- data.frame(IDinit=init, mean_height=c(by(farm$height,farm$IDinit, FUN=mean)), sd_height= c(by(farm$height, farm$IDinit, FUN=sd))) # creation of a dataframe to store initial mean height and initial height standard deviation for the different initial states (across the different fields)
  ######################################################################################################
  
  #########################################################################################
  #### SIMULATION OF THE SUBSEQUENT GENERATIONS ########
  set.seed(seed) ## setting the seed parameter to ensure reproductibility
  phenotype_record <- data.frame(seed=NULL, init=NULL,distance=NULL,generation=NULL, mean_height=NULL, sd_height=NULL, relatedness=NULL, ess1=NULL, ess2=NULL, ess3=NULL)## creating an empty object to store phenotypic evolution accross generations (every ten generation)
  
  
  farm_i <- farm[which(farm$IDinit==init),] 
  farm_i_fec <- list()
  
  for (f in 1:nF) {
    farm_i_f <- farm_i[which(farm_i$IDfield==f),]
    farm_i_f$fecundity <- rpois(nP,fec_extended(farm_i_f$IDgenotype, farm_i_f$height, xopt, d, a))
    farm_i_fec[[f]]  <- farm_i_f
    print(f)
  }
  
  farm_i_fec <- do.call(rbind, farm_i_fec)
  
  for (g in 1:nG) { ## simulation of the the consecutive generations

    generation <- data.frame(IDgenotype=NULL, height=NULL, IDfield=NULL) # creation of an empty data-frame to store each generation
    
    fld_fec <- aggregate(fecundity~IDfield, data=farm_i_fec, FUN=sum) # computing field fecundity as the sum of individual fecundity
    fld_fec$rel_fecundity <- fld_fec$fecundity - max(fld_fec$fecundity)
    
    fld_fec$weights <- exp(alpha*fld_fec$rel_fecundity)/mean(exp(alpha*fld_fec$rel_fecundity)) # computing the weight given to among-field competition using "grp_selec_intensity" parameter
    
    farm_i_fec <- merge(farm_i_fec, fld_fec[,c("IDfield","weights")], by="IDfield", all.x=T) # reporting the weights in the global data.frame
    
    if (sum(fld_fec$fecundity)==0 | (any(fld_fec$fecundity==0) & m==0)) {
      break
    } else {
    
    for (f in 1:nF) { # within one generation, simulations are performed consecutively
      
      farm_i_f <- farm_i_fec[which(farm_i_fec$IDfield==f),] # subseting individuals that are in the good field
      
      #### PICKING FOUNDERS ####        
      
      founders <- vector() ## empty vector
      nf <- 0 ## iteration variable
      
      while (nf < nS) {
        
        pool <- sample(c("mig","res"), 1, prob=c(m, (1-m))) ## deciding wether the founder is picked in the pooled harvest of all fields (mig) with probability "mig_rate" or in the harvest of the current field (res) with probability "1-mig_rate"
        if (pool == "res") {
          
          founder <- sample(farm_i_f$height, 1, prob=farm_i_f$fecundity/sum(farm_i_f$fecundity), replace=T) ## in the current field, the probablity for an individual to be picked up as a founder equals its relative fitness in the field
          
        } else {
          
          founder <- sample(farm_i_fec$height, 1, prob=(farm_i_fec$fecundity/sum(farm_i_fec$fecundity))*farm_i_fec$weights, replace=T) ## in the pooled harvest, the probability for an individual to be picked up as a founder equals its relative fitness in the pooled harvest times the weight given to among-field competition
          
        }
        
        nf <- nf+1
        founders <- c(founders, founder)
        
      }
      
      #### MULTIPLYING FOUNDERS ###
      one_field_generation <- data.frame(IDfield=rep(f, nP), height=rep(founders, each=nP/nS))
      
      ### APLLYING MUTATION ###          
      one_field_generation$height <- apply(one_field_generation,1, function(x) {sample(c(x["height"], x["height"]+mut_effect, x["height"]-mut_effect),1,prob=c(1-mut_rate, mut_rate/2, mut_rate/2))})
      one_field_generation$height[one_field_generation$height<0] <- 0 # replacing negative phenotypic values with 0  
      
      ### RANDOMIZING INDIVIDUALS ###
      one_field_generation <- one_field_generation[sample(1:nrow(one_field_generation)),]
      one_field_generation$IDgenotype <- c(1:nP)
      
      one_field_generation$fecundity<- rpois(nP,fec_extended(one_field_generation$IDgenotype, one_field_generation$height, xopt, d, a))
      generation <- rbind(generation, one_field_generation) ## combining the different fields of the same generation in a single data-frame
      
    }
    
    
    if(g %in% seq(0,nG,by=10)) { ## Every ten generation, phenotypic informations are recorded
      
      pheno <- data.frame(seed=seed,init=init,distance=d,generation=g, mean_height=mean(generation$heigh), sd_height=sd(generation$height), relatedness=relatedness(nS, m, gamma, phi), ess1=ess[1], ess2=ess[2], ess3=ess[3])
      phenotype_record <- rbind(phenotype_record, pheno)
      
    }
    
    farm_i_fec <- generation ## progenies are then now considered as parents for the next generation
    print(g)
    
  }
  }
  
  return(list(rec_init, phenotype_record))
  
}


############################
# CHOSING  PARAMETERS #
############################
nF <- 5
nP <- 200 
mut_rate <- 0.1
mut_effect <- 5
xopt <- 1/(2*0.006)
a <- pi/50
d <- 50
phenotypes <- list(25,50,100,125) 
init <- 2
nG <- 3
nS <- 2
m <- 1
alpha <- 0
seed <- 2019

# name of the output file
output_name <- paste(nS,"f",m,"m",alpha,"comp", nF,"fld",nP, "ind", nG, "g",mut_rate,"mut_r",mut_effect,"mut_e",d,"dist",a,"angle",seed, "seed", sep="_")



################
## RUNNING SIMULATION ##
###############

result <- simul(xopt, a, nG, nF , nP, mut_rate, mut_effect, nS, m, alpha, init, d, seed)
rec_init <- result[[1]]
phenotype_record <- result[[2]]

################
## OUTPUTS ##
###############

write.csv(phenotype_record, file=paste(output_name, "pheno_rec.csv", sep=""), row.names=F) ## .csv output with phenotypic records across generations

