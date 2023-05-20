## These simulations were performed on a Server Dell PowerEdge R650xs equipped 
## with 2x Intel Xeon Silver 4309Y processor at 2.8 GHz and 256 GB of RAM DDR4.
## The total simulation time was 33 minutes and 50 seconds.

## reading R sources
source("invar_LW.R")

## loading 'alabama' (required to constrained optimization)
library(alabama)

## creating Table 1... (it will take long time)
z05 <- summary.LW(Nsize = 10000, coef = c(1,1), sigma = 1, k = 5)
z02 <- summary.LW(Nsize = 10000, coef = c(1,1), sigma = 1, k = 2)
o02 <- summary.LW(Nsize = 10000, coef = c(1,1), sigma = 1, k = -2)
o05 <- summary.LW(Nsize = 10000, coef = c(1,1), sigma = 1, k = -5)

out <- rbind(z05$output, z02$output, o02$output, o05$output) / 100

print(out, digits = 2)
# Output:
#    Wald.A Wald.B  BF.A  BF.B    BF    LM     D
#25   0.058  0.246 0.058 0.184 0.058 0.058 0.058
#50   0.049  0.216 0.049 0.154 0.049 0.049 0.049
#100  0.049  0.155 0.049 0.108 0.049 0.049 0.049
#500  0.051  0.084 0.051 0.065 0.051 0.051 0.051
#25   0.058  0.127 0.058 0.203 0.216 0.058 0.058
#50   0.049  0.104 0.049 0.073 0.049 0.049 0.049
#100  0.049  0.074 0.049 0.065 0.049 0.049 0.049
#500  0.051  0.058 0.051 0.055 0.051 0.051 0.051
#25   0.058  0.155 0.058 0.137 0.080 0.058 0.058
#50   0.049  0.135 0.049 0.085 0.036 0.049 0.049
#100  0.049  0.103 0.049 0.066 0.042 0.049 0.049
#500  0.051  0.058 0.051 0.050 0.051 0.051 0.051
#25   0.058  0.234 0.058 0.163 0.058 0.058 0.058
#50   0.049  0.216 0.049 0.143 0.049 0.049 0.049
#100  0.049  0.171 0.049 0.108 0.049 0.049 0.049
#500  0.051  0.096 0.051 0.060 0.051 0.051 0.051
