*! Table10.do -- Table 10 of the technical note and the numbers of Section
*! 6.9: the households that do not buy, on the reduced Mexican survey
*! (2,477 households; 9% buy no corn, 15% no wheat), vce(svy).
*!   1. easidiag before estimating: the probits, VIF(delta) at the start,
*!      with no variable of the probit only and with the age and the sex of
*!      the head in the probits only (selvars(age isMale));
*!   2. easi without and with the correction: delta, VIF(delta) at the
*!      converged y, the completion of y, the expenditure and own-price
*!      elasticities of the households and of the market;
*!   3. the linearized survey standard errors against vce(bootstrap, svy),
*!      500 replications of the whole procedure (PSUs within strata).
*! About fifteen minutes.

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table10.do
* Nothing needs to be edited.  The check below stops with a clear message
* when the working directory is not replication/.
capture confirm file "master.do"
if _rc {
	di as error "run this script from the replication/ directory:"
	di as error "    cd <path-to-repository>/replication"
	exit 601
}
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"
local OUT "`ROOT'/replication/out"
capture log close t10
log using "`OUT'/selection_mex.log", replace text name(t10)

use "`ROOT'/examples/mex_bench.dta", clear
local M "w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) nolog"
local SV "selvars(age isMale)"

*------------------------------------------------ 1. before estimating
di as txt _n "1. easidiag: the probits and VIF(delta) at the start"
quietly easidiag w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) ///
	demographics(z1 z2) power(3) selection
tempname D0 D1
matrix `D0' = r(sel_diag)
quietly easidiag w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) ///
	demographics(z1 z2) power(3) `SV'
matrix `D1' = r(sel_diag)
di as txt "good" _col(10) "buyers %" _col(21) "pseudo-R2" _col(32) "VIF, no exclusion" _col(52) "VIF, selvars(age isMale)"
forvalues j = 1/2 {
	di as txt "w`j'" _col(10) as res %7.1f `D1'[`j', 1] _col(21) %8.3f `D1'[`j', 2] ///
		_col(34) %10.1f `D0'[`j', 4] _col(58) %10.1f `D1'[`j', 4]
}

*------------------------------------------------ 2. without and with the correction
di as txt _n "2. easi, vce(svy): without and with the correction"
quietly easi `M' vce(svy) notable compensated
tempname E0 P0 E0m P0m
matrix `E0'  = e(elast_exp)
matrix `P0'  = vecdiag(e(elast_price_nc))
matrix `E0m' = e(elast_exp_mkt)
matrix `P0m' = vecdiag(e(elast_price_nc_mkt))
local t0 = e(time)
easi `M' vce(svy) `SV' compensated
tempname E1 E1s P1 P1s E1m P1m G B1 V1 VS1 ES1 ES1m PS1 PS1m
matrix `E1'  = e(elast_exp)
matrix `E1s' = e(elast_exp_se)
matrix `P1'  = vecdiag(e(elast_price_nc))
matrix `P1s' = vecdiag(e(elast_price_nc_se))
matrix `E1m' = e(elast_exp_mkt)
matrix `P1m' = vecdiag(e(elast_price_nc_mkt))
matrix `G'   = e(sel_diag)
matrix `B1'  = e(b)
matrix `VS1' = e(V_sel)
matrix `ES1' = e(elast_exp_se)
matrix `ES1m' = e(elast_exp_mkt_se)
matrix `PS1'  = e(elast_price_nc_se)
matrix `PS1m' = e(elast_price_nc_mkt_se)
local t1 = e(time)
di as txt _n "Table 10."
di as txt "good" _col(9) "delta (z)" _col(28) "VIF" _col(36) "expenditure: none   corrected (se)" ///
	_col(72) "own price: none   corrected (se)"
forvalues j = 1/3 {
	local q : word `j' of w1 w2 w3
	local dl ""
	if `j' < 3 {
		local d = _b[`q':delta]
		local z = `d' / _se[`q':delta]
		local dl : di %7.4f `d' " (" %4.2f `z' ")"
		local vf : di %6.1f `G'[`j', 4]
	}
	else local vf ""
	di as txt "`q'" _col(8) as res "`dl'" _col(26) "`vf'" _col(38) %7.4f `E0'[1, `j'] ///
		_col(49) %7.4f `E1'[1, `j'] " (" %5.3f `E1s'[1, `j'] ")" ///
		_col(76) %7.4f `P0'[1, `j'] _col(87) %7.4f `P1'[1, `j'] " (" %5.3f `P1s'[1, `j'] ")"
}
di as txt "market:"
forvalues j = 1/3 {
	local q : word `j' of w1 w2 w3
	di as txt "`q'" _col(38) as res %7.4f `E0m'[1, `j'] _col(49) %7.4f `E1m'[1, `j'] ///
		_col(76) %7.4f `P0m'[1, `j'] _col(87) %7.4f `P1m'[1, `j']
}
di as txt "completion of y: " as res %4.1f e(sel_ycomp_pct) as txt "% of the households, mean |change| " ///
	as res %6.4f e(sel_ycomp_mean) as txt " = " as res %5.3f e(sel_ycomp_mean) / e(sel_ycomp_sd) as txt " sd of y"
di as txt "execution time: without " as res %5.2f `t0' as txt " s, with " as res %5.2f `t1' as txt " s"

*------------------------------------------------ 3. the survey variance against the design bootstrap
di as txt _n "3. vce(svy) against vce(bootstrap, svy), 500 replications"
timer clear 1
timer on 1
quietly easi `M' `SV' vce(bootstrap, reps(500) seed(20260927) svy) compensated notable
timer off 1
quietly timer list 1
di as txt "bootstrap: " e(N_reps_ok) " of " e(N_reps) " replications, " %6.0f r(t1) " s, resampling `e(boot_design)'"
mata:
real scalar _med(real colvector x)
{
	real colvector y
	real scalar m
	y = sort(x, 1); m = rows(y)
	return(mod(m, 2) ? y[(m + 1) / 2] : (y[m / 2] + y[m / 2 + 1]) / 2)
}
void _cmp(string scalar nm, real rowvector a, real rowvector b)
{
	real rowvector r
	r = a :/ b
	printf("{txt}  %-34s {res}%8.3f %8.3f %8.3f %5.0f\n", nm, _med(r'), min(r), max(r), cols(r))
}
cn  = st_matrixcolstripe(st_local("B1"))[, 1] :+ ":" :+ st_matrixcolstripe(st_local("B1"))[, 2]
K   = rows(cn)
dl  = selectindex(strpos(cn', ":delta") :> 0)
sb  = J(1, K, 1); sb[dl] = J(1, cols(dl), 0)
sa  = sqrt(diagonal(st_matrix(st_local("VS1"))))'
sbo = sqrt(diagonal(st_matrix("e(V_sel)")))'
printf("\n{txt}SE vce(svy) / SE bootstrap (svy), by block\n")
printf("{txt}  %-34s %8s %8s %8s %5s\n", "", "median", "min", "max", "cells")
_cmp("latent coefficients", sa[selectindex(sb)], sbo[selectindex(sb)])
_cmp("delta", sa[dl], sbo[dl])
_cmp("probit coefficients", sa[| K + 1 \ cols(sa) |], sbo[| K + 1 \ cols(sbo) |])
_cmp("expenditure, households", st_matrix(st_local("ES1")), st_matrix("e(elast_exp_se)"))
_cmp("expenditure, market", st_matrix(st_local("ES1m")), st_matrix("e(elast_exp_mkt_se)"))
_cmp("uncompensated, households", vec(st_matrix(st_local("PS1")))', vec(st_matrix("e(elast_price_nc_se)"))')
_cmp("uncompensated, market", vec(st_matrix(st_local("PS1m")))', vec(st_matrix("e(elast_price_nc_mkt_se)"))')
end
log close t10
di as res _n "Table10: done"
