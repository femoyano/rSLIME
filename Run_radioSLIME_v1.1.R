
library(deSolve) #needed for solver
library(lhs) #needed for latin hypercube sampler
library(gdata)#for reading excel files
library(plotrix)#for rescaling the lhs
library(RColorBrewer) # for color palettes
library(hydroGOF) # for error metrics

#╔════════════════════════════════════════════════════════════════════════════╗
#║              load conversion (and/or additional) functions                         
#╚════════════════════════════════════════════════════════════════════════════╝
source("Conversion_functions.R")

#╔════════════════════════════════════════════════════════════════════════════╗
#║                                 load radioSLIME
#╚════════════════════════════════════════════════════════════════════════════╝
source("radioSLIME_v1.1.R")




#╔════════════════════════════════════════════════════════════════════════════╗
#║                    Running the model for testing purposes
#╚════════════════════════════════════════════════════════════════════════════╝

Background14C_data = (read.xls ("Dataset_corrected.xlsx", sheet = 6, header = TRUE))
Background_14C<-Background14C_data$X14C.atmosphere

# define initial conditions for F and H pools
inits <-  c(init_C=30,  # initial C content
            init_C14=aF_to_m14C(30, Background_14C[1]/100)) # initial 14C content, can be a fraction of the above
# define decay and initialization parameters
pars <- c(init_ratio=0.7,  #ratio between the two litter layers
          kF = 0.02, #decay of OF
          kH = 0.15, #decay of OH
          kFH = 0.15, #rate of transformation of OF in OH
          If = 0.55) #input rate for tot C in units of C per unit of time
#run the model
out <- rSLIME(pars = pars, inits = inits, pmc14C_atm=Background14C_data$X14C.atmosphere) #"length" is length of the simulation



##plotting the results
par(mfrow=c(2,1))
plot(out$time, out$OF, main = "SLIME simulation, tot C", ylab = "C units",  xlab = "time", type = "l", col="darkgreen", ylim=c(0,50))
lines(out$time, out$OH, type = "l", col="red")
lines(out$time, out$OH+out$OF, type = "l", col="black")
plot(out$time, out$pmc_F_14C, main = "SLIME simulation, 14C", ylab = "pmC",  xlab = "time", type = "l", col="darkgreen", ylim=c(80,150))
lines(out$time, out$pmc_H_14C, type = "l", col="red")

#╔════════════════════════════════════════════════════════════════════════════╗
#║                                end of test run
#╚════════════════════════════════════════════════════════════════════════════╝




#╔════════════════════════════════════════════════════════════════════════════╗
#║              Preparing the data and simulating missing ones                
#║
#║               STANDARDIZED BOOTSTRAPPING NUMBER = 1000
#╚════════════════════════════════════════════════════════════════════════════╝
# since the dataset has irregular shape, the samples will be standardized to 100
# by bootstrapping. Each case will be reshaped by choosing an appropriate 
# distribution. After that data will be downsampled to 50 to keep it manageable

N.BOOT = 1000

Litter14C_data = na.omit(read.xls ("Dataset_corrected.xlsx", sheet = 2, header = TRUE))
LitterC_data = na.omit(read.xls ("Dataset_corrected.xlsx", sheet = 3, header = TRUE))
Background14C_data = (read.xls ("Dataset_corrected.xlsx", sheet = 6, header = TRUE))

Background_14C<-Background14C_data$X14C.atmosphere


LitterC_data_mod<-cbind(LitterC_data, LitterC_data$mean_C-LitterC_data$variance_C, LitterC_data$mean_C+LitterC_data$variance_C)
colnames(LitterC_data_mod)[12:13]<-c("minC", "maxC")

Litter14C_data_mod<-cbind(Litter14C_data, Litter14C_data$F14C.-Litter14C_data$F14C_err_procent, Litter14C_data$F14C.+Litter14C_data$F14C_err_procent)
colnames(Litter14C_data_mod)[12:13]<-c("min14C", "max14C")

Litter14C_data[Litter14C_data$Site== "K1",]
Litter14C_data[Litter14C_data$Site== "K1" & Litter14C_data$Type== "AoF" ,]


# needed data
# F_floor_depth
# F_floor_BD
# F_floor_C

#K1
K1_F_floor_depth = 2 # cm
K1_F_floor_BD = 0.95 # g/cm3
K1_H_floor_depth = 1.5 # cm
K1_H_floor_BD = 0.95 # g/cm3

K1_F_floor_C = sample(runif(min=min(LitterC_data_mod[LitterC_data_mod$Site== "K1" & LitterC_data_mod$Type== "AoF" ,]$minC),
                      max=max(LitterC_data_mod[LitterC_data_mod$Site== "K1" & LitterC_data_mod$Type== "AoF" ,]$maxC),
                      n=N.BOOT), 50)
K1_F_C <- K1_F_floor_C*K1_F_floor_depth*K1_F_floor_BD

K1_F_floor_14C = sample(runif(min=min(Litter14C_data_mod[Litter14C_data_mod$Site== "K1" & Litter14C_data_mod$Type== "AoF" ,]$min14C),
                              max=max(Litter14C_data_mod[Litter14C_data_mod$Site== "K1" & Litter14C_data_mod$Type== "AoF" ,]$max14C),
                              n=N.BOOT), 50)
K1_F_14C <- K1_F_floor_14C*K1_F_floor_depth*K1_F_floor_BD


K1_H_floor_C = sample(runif(min=min(LitterC_data_mod[LitterC_data_mod$Site== "K1" & LitterC_data_mod$Type== "AoH" ,]$minC),
                            max=max(LitterC_data_mod[LitterC_data_mod$Site== "K1" & LitterC_data_mod$Type== "AoH" ,]$maxC),
                            n=N.BOOT), 50)
K1_H_C <- K1_F_floor_C*K1_F_floor_depth*K1_F_floor_BD

K1_H_floor_14C = sample(runif(min=min(Litter14C_data_mod[Litter14C_data_mod$Site== "K1" & Litter14C_data_mod$Type== "AoH" ,]$min14C),
                              max=max(Litter14C_data_mod[Litter14C_data_mod$Site== "K1" & Litter14C_data_mod$Type== "AoH" ,]$max14C),
                              n=N.BOOT), 50)
K1_H_14C <- K1_H_floor_14C*K1_H_floor_depth*K1_H_floor_BD


#SP3
SP3_F_floor_depth = 2 # cm
SP3_F_floor_BD = 0.95 # g/cm3
SP3_H_floor_depth = 1.5 # cm
SP3_H_floor_BD = 0.95 # g/cm3

SP3_F_floor_C = sample(runif(min=min(LitterC_data_mod[LitterC_data_mod$Site== "SP3" & LitterC_data_mod$Type== "AoF" ,]$minC),
                            max=max(LitterC_data_mod[LitterC_data_mod$Site== "SP3" & LitterC_data_mod$Type== "AoF" ,]$maxC),
                            n=N.BOOT), 50)
SP3_F_C <- SP3_F_floor_C*SP3_F_floor_depth*SP3_F_floor_BD

SP3_F_floor_14C = sample(runif(min=min(Litter14C_data_mod[Litter14C_data_mod$Site== "SP3" & Litter14C_data_mod$Type== "AoF" ,]$min14C),
                              max=max(Litter14C_data_mod[Litter14C_data_mod$Site== "SP3" & Litter14C_data_mod$Type== "AoF" ,]$max14C),
                              n=N.BOOT), 50)
SP3_F_14C <- SP3_F_floor_14C*SP3_F_floor_depth*SP3_F_floor_BD


SP3_H_floor_C = sample(runif(min=min(LitterC_data_mod[LitterC_data_mod$Site== "SP3" & LitterC_data_mod$Type== "AoH" ,]$minC),
                            max=max(LitterC_data_mod[LitterC_data_mod$Site== "SP3" & LitterC_data_mod$Type== "AoH" ,]$maxC),
                            n=N.BOOT), 50)
SP3_H_C <- SP3_F_floor_C*SP3_F_floor_depth*SP3_F_floor_BD

SP3_H_floor_14C = sample(runif(min=min(Litter14C_data_mod[Litter14C_data_mod$Site== "SP3" & Litter14C_data_mod$Type== "AoH" ,]$min14C),
                              max=max(Litter14C_data_mod[Litter14C_data_mod$Site== "SP3" & Litter14C_data_mod$Type== "AoH" ,]$max14C),
                              n=N.BOOT), 50)
SP3_H_14C <- SP3_H_floor_14C*SP3_H_floor_depth*SP3_H_floor_BD


#SP4
SP4_F_floor_depth = 2 # cm
SP4_F_floor_BD = 0.95 # g/cm3
SP4_H_floor_depth = 1.5 # cm
SP4_H_floor_BD = 0.95 # g/cm3

SP4_F_floor_C = sample(runif(min=min(LitterC_data_mod[LitterC_data_mod$Site== "SP4" & LitterC_data_mod$Type== "AoF" ,]$minC),
                            max=max(LitterC_data_mod[LitterC_data_mod$Site== "SP4" & LitterC_data_mod$Type== "AoF" ,]$maxC),
                            n=N.BOOT), 50)
SP4_F_C <- SP4_F_floor_C*SP4_F_floor_depth*SP4_F_floor_BD

SP4_F_floor_14C = sample(runif(min=min(Litter14C_data_mod[Litter14C_data_mod$Site== "SP4" & Litter14C_data_mod$Type== "AoF" ,]$min14C),
                              max=max(Litter14C_data_mod[Litter14C_data_mod$Site== "SP4" & Litter14C_data_mod$Type== "AoF" ,]$max14C),
                              n=N.BOOT), 50)
SP4_F_14C <- SP4_F_floor_14C*SP4_F_floor_depth*SP4_F_floor_BD


SP4_H_floor_C = sample(runif(min=min(LitterC_data_mod[LitterC_data_mod$Site== "SP4" & LitterC_data_mod$Type== "AoH" ,]$minC),
                            max=max(LitterC_data_mod[LitterC_data_mod$Site== "SP4" & LitterC_data_mod$Type== "AoH" ,]$maxC),
                            n=N.BOOT), 50)
SP4_H_C <- SP4_F_floor_C*SP4_F_floor_depth*SP4_F_floor_BD

SP4_H_floor_14C = sample(runif(min=min(Litter14C_data_mod[Litter14C_data_mod$Site== "SP4" & Litter14C_data_mod$Type== "AoH" ,]$min14C),
                              max=max(Litter14C_data_mod[Litter14C_data_mod$Site== "SP4" & Litter14C_data_mod$Type== "AoH" ,]$max14C),
                              n=N.BOOT), 50)
SP4_H_14C <- SP4_H_floor_14C*SP4_H_floor_depth*SP4_H_floor_BD









#╔════════════════════════════════════════════════════════════════════════════╗
#║              Preparing the parameter table, Latin Hypercube sampling               
#║
#║                              NUMBER OF RUNS = N.RUNS
#╚════════════════════════════════════════════════════════════════════════════╝

N.RUNS = 5000
params_names_list<-  c("kF",
                       "kH",
                       "kFH",
                       "If",
                       "init_ratio",
                       "initial_C",
                       "initial_pMC")
n.params=length(params_names_list)

params_table <- as.data.frame(mat.or.vec(N.RUNS, n.params+12))
names(params_table)<-c(params_names_list, 
                       "K1_fit_C_F","K1_fit_C_H", "K1_fit_14C_F", "K1_fit_14C_H",
                       "SP3_fit_C_F","SP3_fit_C_H", "SP3_fit_14C_F", "SP3_fit_14C_H",
                       "SP4_fit_C_F","SP4_fit_C_H", "SP4_fit_14C_F", "SP4_fit_14C_H")
LHS_table<-randomLHS(N.RUNS,n.params)

#declare min and max values for the parameters
params_priors_table<-as.data.frame(mat.or.vec(2,n.params))
names(params_priors_table)<-params_names_list
params_priors_table[,1]<-c(0, 100)
params_priors_table[,2]<-c(0, 50)
params_priors_table[,3]<-c(0, 0.8)
params_priors_table[,4]<-c(0,140) # tons C per hectare
params_priors_table[,5]<-c(0,1)
params_priors_table[,6]<-c(0, 100) # tons C per hectare
params_priors_table[,7]<-c(min(Background_14C[0:5]),max(Background_14C[0:5])) 

params_table[,1]<- rescale(LHS_table[,1], c(params_priors_table[1,1], params_priors_table[2,1]))
params_table[,2]<- rescale(LHS_table[,2], c(params_priors_table[1,2], params_priors_table[2,2]))
params_table[,3]<- rescale(LHS_table[,3], c(params_priors_table[1,3], params_priors_table[2,3]))
params_table[,4]<- rescale(LHS_table[,4], c(params_priors_table[1,4], params_priors_table[2,4]))
params_table[,5]<- rescale(LHS_table[,5], c(params_priors_table[1,5], params_priors_table[2,5]))
params_table[,6]<- rescale(LHS_table[,6], c(params_priors_table[1,6], params_priors_table[2,6]))
params_table[,7]<- rescale(LHS_table[,7], c(params_priors_table[1,7], params_priors_table[2,7]))





#╔════════════════════════════════════════════════════════════════════════════╗
#║ 
#║              Running the model and calculating the fitness!!!!               
#║ 
#╚════════════════════════════════════════════════════════════════════════════╝


#                   ╔═════════════════════════════════════════════╗
#                   ║               Model run loop...               
start.time <- Sys.time()

for (i in 1: N.RUNS)
{

  # define decay and initialization parameters
  pars <- c(kF = params_table[i,1], #decay of OF
            kH = params_table[i,2], #decay of OH
            kFH = params_table[i,3], #rate of transformation of OF in OH
            If = params_table[i,4],
            init_ratio=params_table[i,5]) 
  
  # define initial conditions for F and H pools
  inits <-  c(init_C= params_table[i,6],  # initial C content
              init_C14=aF_to_m14C(30, params_table[i,7]/100)) # initial 14C content, can be a fraction of the above
  
  #run the model
  out <- rSLIME(pars = pars, inits = inits, pmc14C_atm=Background14C_data$X14C.atmosphere) #"length" is length of the simulation

  
#                         ╔═════════════════════════════╗
#                         ║    Fitness calculations...               

  
  params_table[i,]$K1_fit_C_F    <- pbias(tail(out, n=1)$OF,mean(K1_F_floor_C))
  params_table[i,]$K1_fit_C_H    <- pbias(tail(out, n=1)$OH,mean(K1_H_floor_C))
  params_table[i,]$K1_fit_14C_F  <- pbias(tail(out, n=1)$pmc_F_14C,mean(K1_F_floor_14C))
  params_table[i,]$K1_fit_14C_H  <- pbias(tail(out, n=1)$pmc_H_14C,mean(K1_H_floor_14C))
  
  params_table[i,]$SP3_fit_C_F    <- pbias(tail(out, n=1)$OF,mean(SP3_F_floor_C))
  params_table[i,]$SP3_fit_C_H    <- pbias(tail(out, n=1)$OH,mean(SP3_H_floor_C))
  params_table[i,]$SP3_fit_14C_F  <- pbias(tail(out, n=1)$pmc_F_14C,mean(SP3_F_floor_14C))
  params_table[i,]$SP3_fit_14C_H  <- pbias(tail(out, n=1)$pmc_H_14C,mean(SP3_H_floor_14C))
  
  params_table[i,]$SP4_fit_C_F    <- pbias(tail(out, n=1)$OF,mean(SP4_F_floor_C))
  params_table[i,]$SP4_fit_C_H    <- pbias(tail(out, n=1)$OH,mean(SP4_H_floor_C))
  params_table[i,]$SP4_fit_14C_F  <- pbias(tail(out, n=1)$pmc_F_14C,mean(SP4_F_floor_14C))
  params_table[i,]$SP4_fit_14C_H  <- pbias(tail(out, n=1)$pmc_H_14C,mean(SP4_H_floor_14C))
  
  
#                         ╚═════════════════════════════╝
    
  
}
end.time <- Sys.time()
time.taken <- end.time - start.time
#                   ║               End of run loop...               
#                   ╚═════════════════════════════════════════════╝

time.taken




#╔════════════════════════════════════════════════════════════════════════════╗
#║ 
#║                      Calibration metrics and plotting               
#║ 
#╚════════════════════════════════════════════════════════════════════════════╝
parameter_palette<-brewer.pal(10,"Set3")
parameter_palette_alpha<-add.alpha(parameter_palette, 0.6)
site_palette<-brewer.pal(4,"Dark2")[2:4]
site_palette_alpha<-add.alpha(site_palette, 0.6)

#   ╔═══════════════════════════════════════════════════════════════════╗
#   ║                  parameter density plots, priors            
#   ╚═══════════════════════════════════════════════════════════════════╝

par_densities_priors <- list()

for(i in 1:n.params){
  par_densities_priors[[i]]<-density(params_table[,i])
}

svg("Parameter_density_plot_priors.svg", width = 12, height = 12, pointsize = 12)
par(mfrow=c(3,3))
for(i in 1:n.params){
  plot(par_densities_priors[[i]], main=names(params_table)[i])
  polygon(par_densities_priors[[i]], col=parameter_palette_alpha[i])
}
dev.off()


#   ╔═══════════════════════════════════════════════════════════════════╗
#   ║        parameter density plots, posteriors, best criterium           
#   ╚═══════════════════════════════════════════════════════════════════╝


params_table_posteriors_best_K1   <-params_table[order(rowMeans(params_table[,n.params:(n.params+4)])),][1:100,]
params_table_posteriors_best_SP3  <-params_table[order(rowMeans(params_table[,(n.params+5):(n.params+8)])),][1:100,]
params_table_posteriors_best_SP4  <-params_table[order(rowMeans(params_table[,(n.params+9):(n.params+12)])),][1:100,]

par_densities_posteriors_best_K1  <- list()
par_densities_posteriors_best_SP3 <- list()
par_densities_posteriors_best_SP4 <- list()


for(i in 1:n.params){
  par_densities_posteriors_best_K1[[i]]<-density(params_table_posteriors_best_K1[,i])
}
for(i in 1:n.params){
  par_densities_posteriors_best_SP3[[i]]<-density(params_table_posteriors_best_SP3[,i])
}
for(i in 1:n.params){
  par_densities_posteriors_best_SP4[[i]]<-density(params_table_posteriors_best_SP4[,i])
}

svg("Parameter_density_plot_posteriors_best.svg", width = 12, height = 12, pointsize = 12)
par(mfrow=c(3,3))
for(i in 1:n.params){
  plot(par_densities_posteriors_best_K1[[i]], main=names(params_table)[i], 
       xlim=c(min(c(par_densities_posteriors_best_K1[[i]]$x,par_densities_posteriors_best_SP3[[i]]$x,par_densities_posteriors_best_SP4[[i]]$x)),
              max(c(par_densities_posteriors_best_K1[[i]]$x,par_densities_posteriors_best_SP3[[i]]$x,par_densities_posteriors_best_SP4[[i]]$x))),
       ylim=c(0,
              1.2*max(c(par_densities_posteriors_best_K1[[i]]$y,par_densities_posteriors_best_SP3[[i]]$y,par_densities_posteriors_best_SP4[[i]]$y))
                     ))
  polygon(par_densities_posteriors_best_K1[[i]], col=site_palette_alpha[1])
  polygon(par_densities_posteriors_best_SP3[[i]], col=site_palette_alpha[2])
  polygon(par_densities_posteriors_best_SP4[[i]], col=site_palette_alpha[3])
  legend("topleft", c("K1", "SP3", "SP4"), pch=16, col=site_palette_alpha, bty="n")
  
  }
dev.off()



#   ╔═══════════════════════════════════════════════════════════════════╗
#   ║        parameter density plots, posteriors, limits criterium           
#   ╚═══════════════════════════════════════════════════════════════════╝


rowMeans(params_table[,n.params:(n.params+4)])
