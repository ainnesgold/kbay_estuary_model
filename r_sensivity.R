#Sensitivity test on r_H, r_P

library(tidyverse)

H_mort <- 0
P_mort <- 0

Vm = 2    # max phtyoplankton growth rate
Rm = 1.1  # max fish ingestion rate, 1.1
ks = 1 # half-saturation constant for phytoplankton N uptake 
m = 0.1      # phytoplankton mortality
kg = 12      # grazing half saturation or (1/Ivlev constant)
y = 0.1  # unassimilated phytoplankton fraction

e <- 0.2
r_H <- seq(0, 2, by = 0.1)
r_P <- seq(0, 2, by = 0.1)
b <- 0.3
c <- 0.5
d <- 100
K_P <- 50
K_H <- 100

parameter_grid <- expand.grid(r_H = r_H,
                              r_P = r_P)
number_patches = 1
#outputs
outcome_N <- matrix(0, ncol = number_patches, nrow = nrow(parameter_grid))
outcome_PH <- matrix(0, ncol = number_patches, nrow = nrow(parameter_grid))
outcome_H <- matrix(0, ncol = number_patches, nrow = nrow(parameter_grid))
outcome_P <- matrix(0, ncol = number_patches, nrow = nrow(parameter_grid))

ingestion <- function(predation_type, kg, P){
  if (predation_type == 'quadratic'){
    I = P^2/(kg^2+P^2)
  }   
  else if (predation_type == 'Ivlev'){
    I = 1-exp(-kg^P) 
  }
  else if (predation_type == 'M-P'){
    I = kg^-1*P*(1-exp(-kg^-1*P))
  }
  else {
    I = P/(kg+P) 
  }
  return(I)
}
predation_type = 'M-P'


dt = 0.01
timesteps = 365*100

prod <- array(NA, dim = c(timesteps))
phmort <- array(NA, dim = c(timesteps))
grazing <- array(NA, dim = c(timesteps))
PH <- array(NA, dim = c(timesteps))
N <- array(NA, dim = c(timesteps))
H <- array(NA, dim = c(timesteps))
P <- array(NA, dim = c(timesteps))


#LOOP
#Rm (max fish ingestion rate) and Vm (max phyto growth rate) are daily rates i think
#growth rate of of Predator is annual
#conversion of phytoplankton to herbivore biomass
n = nrow(parameter_grid)
outcome_list <- vector("list", length = n)


for (iter in 1:n) {
  
  PH[1] <- 1
  N[1] <- 1
  H[1] <- 20
  P[1] <- 5
  
  for(t in 2:timesteps){
    prod[t-1] = PH[t-1] * ((Vm * N[t-1])/ ks + N[t-1])
    phmort[t-1] = PH[t-1] * m
    grazing[t-1] = Rm * ingestion(predation_type, kg, PH[t-1]) * H[t-1]
    
    PH[t] = PH[t-1] + (prod[t-1] - phmort[t-1] - grazing[t-1]) *dt
    N[t] = N[t-1] + ((-prod[t-1] + phmort[t-1] + y*grazing[t-1] + (H_mort * H[t-1]) + 
                        (P_mort * P[t-1])) *dt)
    
    H[t] = H[t-1] + ((parameter_grid[['r_H']][[iter]] * H[t-1]) * (1 - H[t-1] / K_H) + ((1-y) * grazing[t-1]) - (H_mort*H[t-1]) -
                       (c * H[t-1] * P[t-1]) / (d + H[t-1])
    ) * dt
    
    P[t] = P[t-1] + ((parameter_grid[['r_P']][[iter]] * P[t-1]) * (1 - P[t-1] / K_P) - (P_mort*P[t-1]) +
                       (b * H[t-1] * P[t-1]) / (d + H[t-1])
    ) * dt
  }
  
  outcome_list[[iter]] <- list(PH = PH, N = N, H = H, P = P)
  
  outcome_N[iter] <- N[timesteps]
  outcome_PH[iter] <- PH[timesteps]
  outcome_H[iter] <- H[timesteps]
  outcome_P[iter] <- P[timesteps]
}


#Individual trajectories
plot(1:timesteps, outcome_list[[1]]$N)
plot(1:timesteps, outcome_list[[1]]$PH)
plot(1:timesteps, outcome_list[[1]]$H)
plot(1:timesteps, outcome_list[[1]]$P)

#End value plots
outcome <- cbind(parameter_grid, outcome_N, outcome_PH, outcome_H, outcome_P)

ggplot(outcome, aes(x=r_H, y=outcome_N, col = as.factor(r_P))) +
  geom_point() +
  geom_line()

ggplot(outcome, aes(x=r_H, y=outcome_PH, col = as.factor(r_P))) +
  geom_point() +
  geom_line()

ggplot(outcome, aes(x=r_H, y=outcome_H, col = as.factor(r_P))) +
  geom_point() +
  geom_line()

ggplot(outcome, aes(x=r_H, y=outcome_P, col = as.factor(r_P))) +
  geom_point() +
  geom_line()


