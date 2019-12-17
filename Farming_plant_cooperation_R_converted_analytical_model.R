#### This code is based on the Mathematica notebook "Farming_plant_cooperation_analytical_model.nb". It was converted into R code to be able to source some key functions such as relatedness or selection gradients in the "Farming_plant_cooperation_ind_based_simulation.R" code. #####

### PARAMETERS ###

##'@param x = focal plant height
##'@param xopt = optimal plant height for individual seed production
##'@param yl = left neighbor plant height
##'@param yr = right neighbor plant height
##'@param d = distance between focal plant and neighbor plants
##'@param a = light interception ability, = angle at the summit of the triangles
##'@param z = mean plant height in the meta-population
##'@param nP = number of plants per field
##'@param nF = number of fields in the meta-population
##'@param m = seed migration rate
##'@param alpha = group selection intensity
##'@param phi = parameter that controls the probability that two seeds in the migrant pool come from the same population. When phi=0, all fields contribute to the migrant pool equally, so that two seeds sampled in the migrant pool come from the same population with probability 1/nF. When phi=1, all seeds from the migrant pool come from the same field, so that two seeds sampled in the migrant pool come from the same population with probability 1
##'@param mu = mutation rate




#############################
#### I. FOCAL PLANT FECUNDITY ####
#############################


## "f" function = focal plant fecundity.We assumed a single focal plant is surrounded by two neighbors, on the left and the right sides. The two neighbors are positionned at the same distance from the focal plant,d


f <-
  function (x, xopt, yl, yr, d, a) {
    # arguments : (y: focal plant height), (v : "height cost" parameter), (z_l : left neighbor plant height), (z_r : right neighbor plant height), (d: physical distance between focal and neighbor plants), (a: shading angle (radians))
    
    if (x > 2 * xopt) {
      # focal plant fitness when focal plant height is > 1/v
      
      0
      
    } else {
      if (yl <= (d / tan(a) - x) &
          yr <= (d / tan(a) - x)) {
        # focal plant fitness when both neighbors are too far away from the focal plant
        
        x * (1 - x / (2 * xopt))
        
      } else {
        if (yl > (d / tan(a) + x) |
            yr > (d / tan(a) + x)) {
          # focal plant fitness when one of the neighbors completly shades focal plant
          
          0
          
        } else {
          if (yl <= (d / tan(a) - x) &
              yr > (d / tan(a) - x)) {
            # focal plant fitness when left neighbor is too far away from the focal plant but right neighbor is shading the focal plant
            
            x * (1 - x / (2 * xopt)) * (tan(a) * x ^ 2 - (((
              tan(a) * (x + yr) - d
            ) ^ 2) / (4 * tan(a)))) / (tan(a) * x ^ 2)
            
          } else {
            if (yl > (d / tan(a) - x) &
                yr <= (d / tan(a) - x)) {
              # focal plant fitness when right neighbor is too far away from the focal plant but left neighbor is shading the focal plant
              
              x * (1 - x / (2 * xopt)) * (tan(a) * x ^ 2 - (((
                tan(a) * (x + yl) - d
              ) ^ 2) / (4 * tan(a)))) / (tan(a) * x ^ 2)
              
            } else {
              if (yl > (d / tan(a) - x) &
                  yr > (d / tan(a) - x) &
                  tan(a) * (yl + yr) <= (2 * d)) {
                # focal plant fitness when both neighbors are shading the focal plant with no overlap in their shades
                
                x * (1 - x / (2 * xopt)) * (tan(a) * x ^ 2 - (((
                  tan(a) * (x + yl) - d
                ) ^ 2) / (4 * tan(a))) - (((
                  tan(a) * (x + yr) - d
                ) ^ 2) / (4 * tan(a)))) / (tan(a) * x ^ 2)
                
              } else {
                if (yl > (d / tan(a) - x) &
                    yr > (d / tan(a) - x) &
                    tan(a) * (yl + yr) > (2 * d)) {
                  # focal plant fitness when both neighbors are shading the focal plant with an overlap in their shades
                  
                  x * (1 - x / (2 * xopt)) * (tan(a) * x ^ 2 - (((
                    tan(a) * (x + yl) - d
                  ) ^ 2) / (4 * tan(a))) - (((
                    tan(a) * (x + yr) - d
                  ) ^ 2) / (4 * tan(a))) + (((
                    tan(a) * (yl + yr) - 2 * d
                  ) ^ 2) / (4 * tan(a)))) / (tan(a) * x ^ 2)
                  
                }
              }
            }
          }
        }
      }
    }
  }










#############################
#### II. FOCAL PLANT FITNESS ####
#############################

## We assumed that plant evolved in a field-structured meta-population and that harvested seeds are used for sowing in the next season. We account for the fact that farmers can have different ways to chose seeds at sowing. Seeds could either be picked in the pooled harvest of all fields, in the harvest of the best field alone or each field can be sown with seeds picked in their own harvest.

###
fec <- function(x, y) {
  f(x, xopt, y, y, d, a)
}


### "W" = focal plant fitness ###
W <- function (x, y, z) {
  (1 - m) * W_philo(x, y, z) + m * W_allo(x, y, z)
}

### "W_philo" = Philopatric fitness component (seeds of the focal plant that stay in the parental field) ###
W_philo <- function (x, y, z) {
  fec(x, y) / (comp_philo_loc(x, y) + comp_immig(x, y, z))
}

### "comp_philo_loc" = competition between the seeds of the focal plant that stay in their parental field with the seeds of other individuals in the focal field (resident seeds) ###
comp_philo_loc <- function (x, y) {
  (1 - m) * fec(yR(x, y), yR(x, y))
}

### "yR" = mean phenotype in the patch of the focal, including the focal ###
yR <- function (x, y) {
  (1 / nP) * x + ((nP - 1) / nP) * y
}

### "comp_immig " = competition between the seeds of the focal plant that stay in their parental field with the seeds coming from the migrant pool. There might be seeds from the focal plant in the migrant pool that finally return to their parental field with probability 1/nF ###
comp_immig <- function (x, y, z) {
  m * (1 / nF * (gp_sl_int(yR(x, y), z) * fec(yR(x, y), yR(x, y))) + (nF - 1) / nF * (gp_sl_int(z, z) *
                                                                                        fec(z, z)))
}

### "grp_sel_int" = group selection intensity. Seeds from the most productive fields have more chance to be picked in the migrant pool. We wanted to modulate this phenomenon by introducing paramter alpha. When alpha >>0, only seeds from the most productive field can be picked in the migrant pool ###

gp_sl_int <- function(y, z) {
  exp(alpha * (fec(y, y) - fec(z, z)))
}

### "W_allo" = allopatric fitness component (seeds of the focal plant that go to the migrant pool) ###

W_allo <- function (x, y, z) {
  W_allo_return(x, y, z) + W_allo_true(x, y, z)
}

### "W_allo_return" = when the seeds of the focal plant go to the migrant pool, they can ultimately return to their parental field (with probability 1/nF) ###

W_allo_return <- function (x, y, z) {
  1 / nF * (fec(x, y) * gp_sl_int(yR(x, y), z)) / (comp_philo_loc(x, y) + comp_immig(x, y, z))
}

### "W_allo_true" = when the seeds of the focal plant go to the migrant pool and then migrate to a different population from their parental population, they compete with the local non-migrated seeds and the migrants from other populations ###

W_allo_true <- function (x, y, z) {
  (nF - 1) / nF * (fec(x, y) * gp_sl_int(yR(x, y), z)) / (comp_allo_loc(z) + comp_immig(x, y, z))
}

### "comp_allo_loc" = competition between the seeds of the focal plant that migrate in a different population and the local non-migrated seeds ###
comp_allo_loc <- function (z) {
  (1 - m) * fec(z, z)
}



#############################
#### III. FOCAL PLANT FITNESS APPROXIMATION WHEN NF >>> 0 ###
#############################

### When NF >> 0, terms with 1/nF become negligible compared to termes with (nF-1)/nF. As a results, both comp_immig function and W_allo function are simplified.

Ws <- function (x, y, z) {
  (1 - m) * fec(x, y) / ((1 - m) * fec(y, y) + m * fec(z, z)) + m * fec(x, y) * gp_sl_int(y, z) / fec(z, z)
}




#############################
#### IV. RELATEDNESS ####
#############################
# We define relatedness as R = (Qw-Qb)/(1-Qb) where Qw and Qb are the probabilities of identity of two distinct plants sampled from the same or different fields, repesctively.

R <- expression((Qw - Qb) / (1 - Qb))

# We follow Cockerham and Wei (1987) to compute Qw and Qb.

## A is the probability that two founding seeds from one population at time t were produced in the same population at t-1 ##
Pw <- function (m, phi, nF) {
  (1 - m) ^ 2 + m ^ 2 * (phi + (1 - phi) / nF)
}

## B is the probability that two founding seeds from distinct population at time t were produced in the same population at t-1 ##
Pa <- function(m, phi, nF) {
  (2 * m * (1 - m)) / (nF + m ^ 2 * (phi + (1 - phi) / nF))
}

## From Pw and Pa expressions, we can write the following equations for Qw and Qb

Qw <- function(m, phi, nF, gamma, nS, nP) {
  (gamma * (
    -nF * (-1 + nP + nS) * (-1 + gamma) + 2 * m * (nF * (-1 + nS) * (-1 + gamma) +
                                                     (-1 + nP + nS) * gamma) + m ^ 2 * (
                                                       -1 - nF + nS + nF * nS + 2 * gamma + nF * gamma - nP * gamma - 2 * nS *
                                                         gamma - nF * nS * gamma + (-1 + nF) * (-1 + nS + nP * gamma) * phi
                                                     )
  )) / (nF * ((-1 + gamma) * (-nP * nS + (-1 + m) ^ 2 * (-1 + nP) * (-1 +
                                                                       nS) * gamma) + m ^ 2 * (-1 + nP + nS) * gamma * phi
  ) + m * gamma * (
    2 * nP * nS - 2 * (-1 + nP) * (-1 + nS) * gamma + m * nP * (1 + 2 * nS *
                                                                  (-1 + gamma) - 2 * gamma - phi) - m * (-1 + nS) * (-1 + 2 * gamma + phi)
  ))
}

Qb <- function(m, phi, nF, gamma, nS, nP) {
  (m * gamma * (nS + (-1 + nP) * gamma) * (2 + m * (-1 + (-1 + nF) * phi))) /
    (nF * ((-1 + gamma) * (-nP * nS + (-1 + m) ^ 2 * (-1 + nP) * (-1 + nS) *
                             gamma) + m ^ 2 * (-1 + nP + nS) * gamma * phi
    ) + m * gamma * (
      2 * nP * nS - 2 * (-1 + nP) * (-1 + nS) * gamma + m * nP * (1 + 2 * nS *
                                                                    (-1 + gamma) - 2 * gamma - phi) - m * (-1 + nS) * (-1 + 2 * gamma + phi)
    ))
}

##### ***********************************************
### System resolution
##### ************************************************

relatedness_complete <- function (nF, nP, nS, m, gamma, phi) {
  (Qw(m, phi, nF, gamma, nS, nP) - Qb(m, phi, nF, gamma, nS, nP)) / (1 - Qb(m, phi, nF, gamma, nS, nP))
}

##### ***********************************************
### approximation when nF and nP are large
##### ************************************************

relatedness <- function (nS, m, gamma, phi) {
  -(gamma / (nS * (-1 + ((
    -1 + m
  ) ^ 2) * gamma) - gamma * (1 + m * (-2 + m + m * phi))))
}


#############################
#### V. ESS EVALUATION ####
#############################
#library(Deriv)
############################
# CHOSING  PARAMETERS #
############################
# nF <- 20
# nP <- 200
# mut_rate <- 0.1
# mut_effect <- 5
# xopt <- 1 / (2 * 0.006)
# a <- pi / 50
# d <- 10
# phenotypes <- list(25, 50, 100, 140)
# init <- 2
# nG <- 5000
# nS <- 10
# m <- 1
# alpha <- 100
# gamma <- (1 - mut_rate) ^ 2
# if (alpha == 0) {
#   phi <- 0
# } else {
#   phi <- 1
# }
# 
# 
# dx <- Deriv(Ws, "x")
# dy <- Deriv(Ws, "y")
# 
# 
# equ <- function (z) {
#   dx(z, z, z) + relatedness(nS, m, gamma, phi) * dy(z, z, z)
# }
# 
# 
# ess <- c()
# for (j in c(10E-3:165)) {
#   tryCatch({
#     root <- uniroot(equ, c(j, j + 1))$root
#     ess <- c(ess, root)
#   }, error = function(e)
#     paste("Oops"))
# }
# 
# ess
