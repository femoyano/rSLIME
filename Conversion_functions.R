
#===============================================================================
### Conversion of aF to mass of 14C, m14C
aF_to_m14C <- function
#===============================================================================
(aF             ##<< absolute Fraction modern = A_SN/A_abs
 ,C_sample            ##<< mass of the total sample
 ,A_abs = 1.176e-12*exp((1950-2015)/8267)   ##<< 0.95 times the specific activity of NBS Oxalic Acid I (SRM 4990B)
 ## normalized to delta13CVPDB=-19 per mil absolute. Is equal to an
 ## activity (AD 1950) of 1.176 ? 0.010e-12
 
){
  m14C <- aF*A_abs*C_sample
  m14C
  ### mass of 14C in [M]
}

# Author: B. Ahrens, modified by me for decay
#===============================================================================


#===============================================================================
### Conversion of mass of 14C, m14C, to Af
m14C_to_aF <- function
#===============================================================================
(m14C                 ##<< mass of 14C in the sample
 ,C_sample            ##<< mass of the C total in the  sample
 ,A_abs = 1.176e-12*exp((1950-2015)/8267)   ##<< 0.95 times the specific activity of NBS Oxalic Acid I (SRM 4990B)
 ## normalized to delta13CVPDB=-19 per mil absolute. Is equal to an
 ## activity (AD 1950) of 1.176 ? 0.010e-12
 
){
  aF <- m14C/(A_abs*C_sample)
  aF
  ### mass of 14C in [M]
}
#===============================================================================



#===============================================================================
### Add an alpha value to a color
add.alpha <- function(col, alpha=1){
#===============================================================================
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

