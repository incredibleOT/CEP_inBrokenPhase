L_s                       64
L_t                       128
antiperiodic_L_t          0
	
# m0Squared                   0.35

## uses kappa instead of m0Squared. If set, setting for m0Squared are ignored
use_kappa                  1
kappa                   0.1173062

lambda                    -0.390

lambda_6                    0.10

## 175/246=0.711382114    500/246=2.032520325
yukawa_t                  0.711382114

##y_t/y_b
yukawa_ratio               1


field_min                  0.04
field_max                  0.2
field_step                 0.001

use_listOfFermContr        0
##a file withe v   U_f(v) 
# listOfFermContr            filename.txt


## default: N_f=1, rho=1, r=0.5
# N_f                       1
# rho                       1.0
# r                         0.5

## default: absolut_tolerance_for_minimization=1.0e-1, relative_tolerance_for_minimization=1.0e-7, tolerance_for_HiggsMassSquared=1.0e-5
absolut_tolerance_for_minimization    1.0e-7
relative_tolerance_for_minimization   1.0e-7
tolerance_for_HiggsMassSquared        1.0e-5

##
# max_numer_of_iterations_minimizer            100
# max_numer_of_iterations_HiggsMassSquared     100
	
	
#1 gsl_min_fminimizer_goldensection, 2 gsl_min_fminimizer_brent, 3 gsl_min_fminimizer_quad_golden
minimization_algorithm         2

# before the minimization starts, a scan is performed to find the initial guess
testvalue_min                    1.0e-20 
testvalue_max                    5.05 
testvalue_step                    0.05


#may use: [Ls]->Lxx [Lt]->Txx [m0Sq_k]->m0Sq_xxxx or k_xxxx [l]->l_xxxx [l6]->l6_xxxx [yt]->yt_xxxx [yr]->yr_xxxxfor the scanables: e.g. k_xx_yy with max and min)
outputfile                  output/potential_[Ls][Lt]_[yt]_[yr]_[l]_[l6]_[m0Sq_k]_[field].txt

