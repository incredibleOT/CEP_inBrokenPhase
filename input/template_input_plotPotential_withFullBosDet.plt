L_s                       4
L_t                       4
antiperiodic_L_t          0
	
m0Squared                   0.15

## uses kappa instead of m0Squared. If set, setting for m0Squared are ignored
use_kappa                  0
# kappa                   0.15

lambda                    -0.008

lambda_6                    0.001

## 175/246=0.711382114    500/246=2.032520325
yukawa_t                  0.711382114

##y_t/y_b
yukawa_ratio               1


field_min                  0.01
field_max                  1.00
field_step                 0.01


## default: N_f=1, rho=1, r=0.5
# N_f                       1
# rho                       1.0
# r                         0.5

#may use: [Ls]->Lxx [Lt]->Txx [m0Sq_k]->m0Sq_xxxx or k_xxxx [l]->l_xxxx [l6]->l6_xxxx [yt]->yt_xxxx [yr]->yr_xxxxfor the scanables: e.g. k_xx_yy with max and min)
outputfile                  output/potential_prodprod_[Ls][Lt]_[yt]_[yr]_[l]_[l6]_[m0Sq_k]_[field].txt

