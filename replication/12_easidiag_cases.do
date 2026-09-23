*! 12_easidiag_cases.do -- section 6.1 of the note: the four easidiag cases
*!
*!   A. the Mexican bench with prices and expenditure in raw levels
*!   A2. the same, centred
*!   B. an exhaustive set of household-composition dummies
*!   C. a rare binary demographic under pz
*!   D. hixdata, the reference specification
*!
*! Produces out/easidiag_cases.log, the reports quoted in the note.
*! Run from the replication/ directory:  do 12_easidiag_cases.do

clear all
set more off
set linesize 100
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"
local OUT "`ROOT'/replication/out"

log using "`OUT'/easidiag_cases.log", replace text name(d)

use "`ROOT'/examples/mex_bench.dta", clear
di as txt _n "==================== A. raw levels"
easidiag w1 w2 w3, lnprices(lp1_raw lp2_raw lp3) lnexpenditure(lx_raw)	///
	demographics(hhsize isMale) power(3)

di as txt _n "==================== A2. centred"
easidiag w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx)		///
	demographics(hhsize isMale) power(3)

di as txt _n "==================== B. exhaustive dummies"
easidiag w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) power(3)	///
	demographics(nocup0 nocup1 nocup2 nocup3 nocup4)
* the estimator itself converges on this specification without a word
qui easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) power(3)	///
	demographics(nocup0 nocup1 nocup2 nocup3 nocup4) nolog noelastse
di as txt "   easi on the same specification: converged = " e(converged)

* Table 5 of the note: the dummy set perturbed away from exact collinearity.
* Exact and 1e-12 are benign (the generalised inverse pins one column at
* zero); 1e-9 leaves the matrix full rank and the estimate explodes.
set seed 11
qui gen double u = runiform() - 0.5
di as txt _n "   perturbation      max |coef|     max std. err."
foreach eps in 0 1e-12 1e-9 {
	foreach v in nocup0 nocup1 nocup2 nocup3 nocup4 {
		capture drop p_`v'
		qui gen double p_`v' = `v' + `eps' * u
	}
	qui easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) power(3)	///
		demographics(p_nocup0 p_nocup1 p_nocup2 p_nocup3 p_nocup4)	///
		nolog noelastse
	mata: st_numscalar("bmax", max(abs(st_matrix("e(b)"))))
	mata: st_numscalar("smax", max(sqrt(diagonal(st_matrix("e(V)")))))
	di as txt "   " %-14s "`eps'" as res %12.4f bmax %16.4f smax
}
drop p_nocup* u

di as txt _n "==================== C. rare category under pz"
set seed 7
qui gen byte rare = runiform() < 0.012
easidiag w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) power(3)	///
	demographics(hhsize isMale rare) pz

di as txt _n "==================== D. hixdata, reference specification"
use "`ROOT'/examples/hixdata.dta", clear
easidiag sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,	///
	lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers)	///
	lnexpenditure(log_y) demographics(age hsex carown time tran)		///
	power(5) py pz zy
local pred = r(iter)
qui easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,	///
	lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers)	///
	lnexpenditure(log_y) demographics(age hsex carown time tran)		///
	power(5) py pz zy nolog noelastse
di as txt "   iterations predicted " `pred' ", observed " e(iter)

log close d
di as res _n "12_easidiag_cases: done"
