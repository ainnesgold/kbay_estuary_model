library(tidyverse)


H_mort <- 0.1
P_mort <- 0.1
Vm = 2    # max phtyoplankton growth rate
Rm = 1.1  # max fish ingestion rate, 1.1
ks = 1 # half-saturation constant for phytoplankton N uptake 
m = 0.1      # phytoplankton mortality
kg = 12      # grazing half saturation or (1/Ivlev constant)
y = 0.1  # unassimilated phytoplankton fraction
#g = 0.001      # fish mortality
e <- 0.2
r_H <- 0.3
r_P <- 0.3
b <- 0.3
c <- 0.5
d <- 100
K_P <- 50
K_H <- 100


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
pmort <- array(NA, dim = c(timesteps))
grazing <- array(NA, dim = c(timesteps))
PH <- array(NA, dim = c(timesteps))
N <- array(NA, dim = c(timesteps))
H <- array(NA, dim = c(timesteps))
P <- array(NA, dim = c(timesteps))

PH[1] <- 1
N[1] <- 1
H[1] <- 20
P[1] <- 5

#LOOP
#Rm (max fish ingestion rate) and Vm (max phyto growth rate) are daily rates i think
#growth rate of of Predator is annual
#conversion of phytoplankton to herbivore biomass

for(t in 2:timesteps){
prod[t-1] = PH[t-1] * ((Vm * N[t-1])/ ks + N[t-1])
pmort[t-1] = PH[t-1] * m
grazing[t-1] = Rm * ingestion(predation_type, kg, PH[t-1]) * H[t-1]

PH[t] = PH[t-1] + (prod[t-1] - pmort[t-1] - grazing[t-1]) *dt
N[t] = N[t-1] + ((-prod[t-1] + pmort[t-1] + y*grazing[t-1] + (H_mort * H[t-1]) + 
                    (P_mort * P[t-1])) *dt)

H[t] = H[t-1] + ((r_H * H[t-1]) * (1 - H[t-1] / K_H) + ((1-y) * grazing[t-1]) - (H_mort*H[t-1]) -
                   (c * H[t-1] * P[t-1]) / (d + H[t-1])
) * dt

P[t] = P[t-1] + (r_P * P[t-1] * (1 - P[t-1] / K_P) - (P_mort*P[t-1]) +
                   (b * H[t-1] * P[t-1]) / (d + H[t-1])
                ) * dt
}


df <- data.frame(N, PH, H, P)
time <- c(1:timesteps)
df <- cbind(df, time)
df_long <- df %>% pivot_longer(N:P, names_to = "variable", values_to = "value")

ggplot(data=df_long, aes(x=time, y = value, col=variable)) + geom_line()














##old
#H[t] = H[t-1] + (((1-y) * grazing[t-1]) * (1 - H[t-1] / K_H) - (H_mort*H[t-1]) - 
#                  (c * H[t-1] * P[t-1] / (d + H[t-1]))) * dt

#H[t] = H[t-1] + (((1-y) * grazing[t-1]) * (1 - H[t-1] / K_H) - (H_mort*H[t-1]) -
#                  (c * H[t-1] * P[t-1]) / (d + H[t-1])
#               ) * dt

#P[t] = P[t-1] + (e * c * H[t-1] * P[t-1] / (d + H[t-1]) - (P_mort*P[t-1])) * dt