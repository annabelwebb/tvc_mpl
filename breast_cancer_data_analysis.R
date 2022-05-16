
bc_longformat = read.table("bc_longformat.txt", header = TRUE, row.names = NULL)
bc_baseline = read.table("bc_baseline.txt", header = TRUE, row.names = NULL)
bc_left = read.table("bc_left.txt", header = TRUE, row.names = NULL)

source("functions_diffN.R")
source("fitting_diffN.R")
library(splines2)
library(dplyr)

ctrl = tvc_mpl_control(lambda = 0, iter = c(10,3000), n_knots = 6, par_initial = c(0,0,1), 
                       range = c(0.1,0.9), line_search = c(1, 1, 1), reg_conv = 1e-5)
test = tvc_fit((bc_longformat), (bc_baseline), (bc_left), ctrl)

test$parameters
test$se_Q

2*(1 - pnorm( abs(c(test$parameters$beta, test$parameters$gamma)/test$se_Q[1:6])))







