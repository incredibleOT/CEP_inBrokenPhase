L_s                       16
L_t                       32
antiperiodic_L_t          0
	
# scan_m0Squared            1
m0Squared                   0.15
# m0Squared_min               -0.5
# m0Squared_max               +0.5
# m0Squared_step              +0.01

## uses kappa instead of m0Squared. If set, setting for m0Squared are ignored
use_kappa                  0
# scan_kappa            1
# kappa                   0.12
# kappa_min               -0.1220
# kappa_max               +0.1235
# kappa_step              +0.00001

# scan_lambda             1
lambda                    0.03
# lambda_min                0.0
# lambda_max                0.5
# lambda_step               0.01

# scan_lambda_6             1
lambda_6                    0.03
# lambda_6_min                0.0
# lambda_6_max                0.5
# lambda_6_step               0.01

## 175/246=0.711382114    500/246=2.032520325
# scan_yukawa_t           1
yukawa_t                  0.711382114
# yukawa_t_min            0.5
# yukawa_t_max            3.0
# yukawa_t_step           0.1

##y_t/y_b
# scan_yukawa_ratio         0
yukawa_ratio               1
# yukawa_ratio_min         
# yukawa_ratio_max
# yukawa_ratio_step

#if set, the contributions from the goldstones will be ignored (basically setting simlpy the propagator sums to zero))
exclude_goldstones         0

use_listOfFermContr        0
##a file withe v   U_f(v) 
# listOfFermContr            filename.txt

## default: N_f=1, rho=1, r=0.5
# N_f                       1
# rho                       1.0
# r                         0.5

## default: absolut_tolerance_for_minimization=1.0e-7, relative_tolerance_for_minimization=1.0e-7, tolerance_for_HiggsMassSquared=1.0e-5
absolut_tolerance_for_minimization    1.0e-7
relative_tolerance_for_minimization   1.0e-7
tolerance_for_HiggsMassSquared        1.0e-5

##
# max_numer_of_iterations_minimizer            100
# max_numer_of_iterations_HiggsMassSquared     100
	
	
#1 gsl_min_fminimizer_goldensection, 2 gsl_min_fminimizer_brent, 3 gsl_min_fminimizer_quad_golden
minimization_algorithm         2

# before the minimization starts, a scan is performed to find the initial guess
testvalue_min                    1.0e-7 
testvalue_max                    3.05 
testvalue_step                    0.05


#may use: [Ls]->Lxx [Lt]->Txx [m0Sq_k]->m0Sq_xxxx or k_xxxx [l]->l_xxxx [l6]->l6_xxxx [yt]->yt_xxxx [yr]->yr_xxxxfor the scanables: e.g. k_xx_yy with max and min)
outputfile                  output/output_[Ls][Lt]_[yt]_[yr]_[l]_[l6]_[m0Sq_k].txt

