# The simulation code of PEER

## Requirement

* Pre-installed R packages : rrpack, glmnet, BeSS, MASS, secure, doParallel, foreach

* Installation of peer package: R command "install_github("DongRuipeng/peer-code/peer")" in "devtools" R packages 

## Simulation script

* All simulations results are save as RData files firstly then formatted as corresponding tables and figures. 

* Tables 1 and 2 are generated by "sim-table-default.R" with n.fix = 100, 200.

* Table 3 is generated by "sim-table-secure.R" with default settings. 

* Table 4 is generated by "sim-real.R" with default settings. 

* Figure 2 is generated by "sim-miss.R" with default settings. 

* After getting the corresponding RData file, scripts in folder format-script will generate the final tables and figures. 

## Supplement

* We employ clustering computing to speed up the computation. The parameter "cl.num" is the thread number, which should be suitable to the current computing device. 
