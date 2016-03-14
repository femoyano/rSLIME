
#▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄
#
#
# ║║╔║║╔╗ ║
# ╠╣╠║║║║ ║
# ║║╚╚╚╚╝ O
#
# Model name(short):  SLIME
# Model name:         Successional LItter ModEl 
# Author:             Lorenzo Menichetti
# Author email:       ilmenichetti@gmail.com
# ★★★please contact the author before using it publicly (only to let me know...)★★★
#
# require(deSolve)
# the package is not loaded in order to save on execution time, needs to be loaded outside
#
# input parameters (list "pars to pass to the command)
# init_ratio :ratio between the two litter layers
# kF         :decay of OF
# kH         :decay of OH
# kFH        :rate of transformation of OF in OH
# If         :input rate for tot C in units of C per unit of time
#
# input inits (list "inits" to pass to the command)
# init_C initial C content
# init_C14= aF_to_m14C(init_C, initial pMC/100)) initial 14C content. It's a mass to save on execution time by skipping the conversion
#
#▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄▀▄



##############################################################################
############ defining the SLIME  system of differential functions ############
##############################################################################

rSLIME <- function (pars, inits, pmc14C_atm)  {
  with (as.list(c(pars, inits)),{ #opens a "with" environment in order to use vectors for inputs
        
  ########################### radiocarbon conversions ########################### 
  aF_D14C  <- pmc14C_atm/100 # from pmc to Delta units
  If14_ts<- aF_to_m14C(aF=aF_D14C, C_sample= If )
  If14_fun <- splinefun(If14_ts, method="monoH.FC")
  ########################### radiocarbon conversions ########################### 
  
    
  ##finding the derivatives. Refer to the model schema
    derivs <- function(time, y, pars) {#y is the vector of the state variables
      with (as.list(c(pars, y)), {
          
        #tot C module
          dOF = If-kF*OF-kFH*OF
          dIH = kFH*OF
          dOH = dIH-kH*OH
          
        #14C module
        #CONV# radiocarbon signature of inputs in C units
          If14  <- If14_fun(time)
          dOF14 = If14-kF*OF14-kFH*OF14
          dIH14 = kFH*OF14
          dOH14 = dIH14-kH*OH14
          
          return(list(c(dOF, dOH, dOF14, dOH14), If14 = If14))
        } )
      }
  ##end of derivatives part  
  
  #initial conditions
  OF_0 = init_C*init_ratio
  OH_0 = init_C*(1-init_ratio)
  OF14_0 = init_C14*init_ratio
  OH14_0 = init_C14*(1-init_ratio)
  
  y <- c(OF=OF_0, OH=OH_0, OF14=OF14_0, OH14=OH14_0)

  times <- seq(0,length(pmc14C_atm)+3,1)

  out <- as.data.frame(ode(y=y, parms = pars, times = times, func= derivs))
  
  F_14C<-out[,4] #14C of F litter
  H_14C<-out[,5] #14C of H litter
  aF_F_14C<- m14C_to_aF(m14C= F_14C, C_sample = out[,2])
  aF_H_14C<- m14C_to_aF(m14C= H_14C, C_sample = out[,3])
  pmc_F_14C<-(aF_F_14C*100)
  pmc_H_14C<-(aF_H_14C*100)
  
  return(cbind(out, pmc_F_14C, pmc_H_14C))
  })  # !! this closes the "with" environment
  
}

##############################################################################
############            end of differential functions             ############
##############################################################################



