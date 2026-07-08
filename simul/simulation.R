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

# Output:
out
#     Wald.A Wald.B   BF.A   BF.B     BF     LM      D
# 25  0.0580 0.2461 0.0580 0.1837 0.0580 0.0580 0.0580
# 50  0.0492 0.2163 0.0492 0.1538 0.0492 0.0492 0.0492
# 100 0.0489 0.1553 0.0489 0.1080 0.0489 0.0489 0.0489
# 500 0.0509 0.0838 0.0509 0.0653 0.0509 0.0509 0.0509
# 25  0.0580 0.1274 0.0580 0.2026 0.2161 0.0580 0.0580
# 50  0.0492 0.1042 0.0492 0.0729 0.0491 0.0492 0.0492
# 100 0.0489 0.0744 0.0489 0.0646 0.0489 0.0489 0.0489
# 500 0.0509 0.0580 0.0509 0.0553 0.0509 0.0509 0.0509
# 25  0.0580 0.1548 0.0580 0.1373 0.0800 0.0580 0.0580
# 50  0.0492 0.1345 0.0492 0.0850 0.0364 0.0492 0.0492
# 100 0.0489 0.1034 0.0489 0.0657 0.0421 0.0489 0.0489
# 500 0.0509 0.0581 0.0509 0.0496 0.0509 0.0509 0.0509
# 25  0.0580 0.2342 0.0580 0.1626 0.0580 0.0580 0.0580
# 50  0.0492 0.2165 0.0492 0.1430 0.0492 0.0492 0.0492
# 100 0.0489 0.1710 0.0489 0.1083 0.0489 0.0489 0.0489
# 500 0.0509 0.0961 0.0509 0.0597 0.0509 0.0509 0.0509
