*! mex_tour.do -- easi on a real household survey with a sampling design
*!
*! mex_bench.dta is a reduced extract of Mexico's ENIGH 2014 (cereals):
*! 2,477 households, 703 primary sampling units, 28 strata, already -svyset-
*! with the expansion factor sweight.  Three goods: corn (w1), wheat (w2) and
*! the composite of everything else (w3); log prices lp1 lp2 lp3 and log
*! per-capita expenditure lx, centred on their means (the uncentred logs are
*! lx_raw lp1_raw lp2_raw); demographics z1 (household size), z2 (urban),
*! age, isMale, and an exhaustive set of composition dummies nocup0-nocup4.
*!
*! Run it line by line.  Every step is a few seconds.

clear all
set more off
use mex_bench.dta, clear
describe, short
svyset

* ---- 1. look before you estimate ------------------------------------------
* The diagnostic reads the design at the Stone-index starting point and does
* not estimate anything.  On the raw logs it flags the scaling; on the
* centred ones it predicts the number of iterations.
easidiag w1 w2 w3, lnprices(lp1_raw lp2_raw lp3) lnexpenditure(lx_raw)	///
	demographics(z1 z2) power(3)
easidiag w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx)		///
	demographics(z1 z2) power(3)

* An exhaustive set of dummies: the estimator would run, pinning one
* coefficient to zero without a word.  The diagnostic names the group.
easidiag w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx)		///
	demographics(nocup0 nocup1 nocup2 nocup3 nocup4) power(3)

* ---- 2. the same model under four variance estimators ---------------------
* Unweighted, robust: the households are treated as independent draws.
easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx)			///
	demographics(z1 z2) power(3) py nolog
estimates store srs

* Sampling weights: the elasticities become population quantities.
easi w1 w2 w3 [pw=sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx)	///
	demographics(z1 z2) power(3) py nolog
estimates store pw

* Clusters: scores summed within PSU.
easi w1 w2 w3 [pw=sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx)	///
	demographics(z1 z2) power(3) py vce(cluster pu) nolog
estimates store cl

* The design as declared: PSUs within strata, weight taken from svyset.
* This is the estimator the note validates against a design bootstrap
* (Table 3): ignoring the clusters understates the standard errors by 7.5%
* on average and 36% at worst on these data.
easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx)			///
	demographics(z1 z2) power(3) py vce(svy) nolog
estimates store svy
di as txt "PSUs " e(N_psu) ", strata " e(N_strata) ", weight " e(svywvar)

* Standard errors of the expenditure elasticities, side by side
foreach m in srs pw cl svy {
	estimates restore `m'
	matrix S = e(elast_exp_se)
	di as txt "  `m'" _col(8) as res %8.4f S[1,1] "  " %8.4f S[1,2] "  " %8.4f S[1,3]
}

* ---- 3. what the note calls the ratio of two estimators -------------------
* The expenditure elasticity of a good with a 2% share divides by an
* estimated mean share.  The standard error above already carries it; the
* full influence function is what makes vce(svy) agree with a design
* bootstrap on goods 1 and 2 (Table 4 of the note).
estimates restore svy
matrix list e(elast_exp)
matrix list e(elast_exp_se)

* ---- 4. everything else works as on any e-class command -------------------
predict double shat*, shares
su shat*
predict double yhat, y
estat engel, n(60)
easi, compensated checks			// replay with more tables
matrix list e(elast_price_c)

* ---- 5. the design bootstrap, if you want to see for yourself -------------
* 50 replications take under a minute; the note uses 400.  Resampling PSUs
* within strata is what makes this an oracle for vce(svy): a bootstrap over
* households would share its blind spot.
capture program drop _eb
program _eb, eclass
	version 14.2
	qui easi w1 w2 w3 [pw=sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) ///
		demographics(z1 z2) power(3) py nolog noelastse
end
set seed 2026
bootstrap _b, reps(50) strata(st) cluster(pu) idcluster(_pu) nowarn nodots: _eb
matrix Vb = e(V)
estimates restore svy
matrix Va = e(V)
mata: r = sort(sqrt(diagonal(st_matrix("Vb"))) :/ sqrt(diagonal(st_matrix("Va"))), 1); st_numscalar("ratio", r[ceil(rows(r)/2)])
di as txt "median(bootstrap SE / vce(svy) SE) over the coefficients = " as res %6.3f ratio
