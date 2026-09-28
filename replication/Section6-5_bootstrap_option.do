*! test_step25.do -- easi 2.0.0, phase C: vce(bootstrap) of the whole
*! procedure (the design of equaids).
*!   1. the point estimates are those of the analytic estimation; e(V), the
*!      standard errors of the elasticities (households and market) and, with
*!      selection, e(V_sel) are replaced by those of the replications;
*!   2. the bootstrap standard errors are of the order of the analytic ones
*!      (a loose band: B = 20 only);
*!   3. svy: the PSUs within strata; a seed reproduces the draws.
clear all
set more off
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Section6-5_bootstrap_option.do
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

global NFAIL 0
capture program drop _chk
program define _chk
	args label value tol
	if (`value' <= `tol') di as txt "  ok    " as res %10.3g `value' as txt "  `label'"
	else {
		di as err "  FAIL  " %10.3g `value' "  `label' (tolerance `tol')"
		global NFAIL = ${NFAIL} + 1
	}
end

use "`ROOT'/examples/mex_bench.dta", clear
local M "w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) nolog notable"

foreach sel in "" selection {
	di as txt _n "1. `=cond("`sel'" == "", "without", "with")' selection"
	quietly easi `M' vce(cluster psu) `sel'
	tempname B V E Es
	matrix `B' = e(b)
	matrix `V' = e(V)
	matrix `E' = e(elast_exp)
	matrix `Es' = e(elast_exp_se)
	quietly easi `M' `sel' vce(bootstrap, reps(20) seed(11) psu(psu))
	_chk "the same point estimates" `=mreldif(e(b), `B')' 1e-12
	_chk "the same elasticities" `=mreldif(e(elast_exp), `E')' 1e-12
	_chk "e(vce) = bootstrap, 20 replications, all converged, the seed" ///
		`=("`e(vce)'" != "bootstrap") + (e(N_reps) != 20) + (e(N_reps_ok) != 20) + ("`e(boot_seed)'" != "11")' 0
	_chk "e(V) replaced" `=(mreldif(e(V), `V') < 1e-6)' 0
	mata: st_numscalar("__r", max(abs(ln(st_matrix("e(elast_exp_se)") :/ st_matrix("`Es'")))))
	_chk "bootstrap / analytic SE of the expenditure elasticities within a factor 2" `=(__r > ln(2))' 0
	capture confirm matrix e(elast_exp_mkt_se)
	_chk "market SEs present" `=_rc != 0' 0
	if "`sel'" != "" {
		capture confirm matrix e(V_sel)
		_chk "e(V_sel) present" `=_rc != 0' 0
	}
}

di as txt _n "3. svy and the seed"
quietly easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) nolog notable ///
	vce(bootstrap, reps(10) seed(3) svy)
tempname V1
matrix `V1' = e(V)
_chk "svy: resampling PSUs within strata" `=("`e(boot_design)'" != "PSUs psu within strata strata")' 0
quietly easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) nolog notable ///
	vce(bootstrap, reps(10) seed(3) svy)
_chk "the same seed, the same variance" `=mreldif(e(V), `V1')' 0
capture easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) nolog vce(bootstrap, reps(1))
_chk "reps(1) refused (rc `=_rc')" `=abs(_rc - 198)' 0

di as txt _n "{hline 60}"
if ${NFAIL} == 0 di as res "test_step25: all checks passed"
else di as err "test_step25: ${NFAIL} check(s) failed"
