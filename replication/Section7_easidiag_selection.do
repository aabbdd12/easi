*! test_step26.do -- easi 2.0.0, phase C: the selection diagnostics of
*! easidiag, before any estimation.
*!   1. the probits are those of the estimation: buyers %, pseudo-R2 and
*!      separation equal e(sel_diag) of easi on the same specification;
*!   2. oracle: VIF(delta) at y = Stone index is 1/(1 - R2), R2 the
*!      uncentred R2 of the regression of phi on the columns Phi x(y) (regress,
*!      noconstant), Phi and phi from Stata's probit;
*!   3. selvars(): the probit-only variables are counted and change the
*!      probits as in the estimation; the refusals of the estimation hold.
clear all
set more off
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Section7_easidiag_selection.do
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
local M "w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3)"

di as txt _n "1. The probits of the estimation"
easidiag `M' selection
tempname D E
matrix `D' = r(sel_diag)
matrix list `D'
_chk "r(selgoods) = w1 w2" `=("`r(selgoods)'" != "w1 w2")' 0
quietly easi `M' selection vce(robust) nolog notable
matrix `E' = e(sel_diag)
matrix __d = `D'[1..2, 1..3]
matrix __e = `E'[1..2, 1..3]
_chk "buyers %, pseudo-R2, separation = e(sel_diag)" `=mreldif(__d, __e)' 1e-10
* at the start and at the converged y: the same order of magnitude
mata: st_numscalar("__r", max(abs(ln(st_matrix("`D'")[., 4] :/ st_matrix("`E'")[., 4]))))
_chk "VIF(delta) at the start within a factor 3 of the converged one" `=(__r > ln(3))' 0

di as txt _n "2. Oracle: VIF(delta) = 1/(1 - R2)"
quietly gen double ystone = lx - (w1 * lp1 + w2 * lp2 + w3 * lp3)
quietly gen double y2 = ystone^2
quietly gen double y3 = ystone^3
quietly gen double np1 = lp1 - lp3
quietly gen double np2 = lp2 - lp3
forvalues j = 1/3 {
	quietly gen double s`j' = lp`j' - lx
}
forvalues i = 1/2 {
	quietly gen byte d`i' = w`i' > 0
	quietly probit d`i' s1 s2 s3 z1 z2 [pw = sweight]
	quietly predict double xb`i', xb
	quietly gen double Ph`i' = normal(xb`i')
	quietly gen double ph`i' = normalden(xb`i')
	local cols ""
	foreach v in one ystone y2 y3 z1 z2 np1 np2 {
		if "`v'" == "one" quietly gen double P`i'_one = Ph`i'
		else quietly gen double P`i'_`v' = Ph`i' * `v'
		local cols `cols' P`i'_`v'
	}
	quietly regress ph`i' `cols' [aw = sweight], noconstant
	local vif = 1 / (1 - e(r2))
	_chk "w`i': VIF(delta) = 1/(1 - R2) (`: di %8.2f `vif'')" `=reldif(`D'[`i', 4], `vif')' 1e-6
}

di as txt _n "3. selvars() and the refusals"
set seed 26
quietly gen double q = runiform()
easidiag `M' selvars(w1: q)
matrix `D' = r(sel_diag)
* selvars() adds variables to the probits; the goods corrected are still
* those with zero shares
_chk "selvars(w1: q): w1 and w2 still corrected" `=("`r(selgoods)'" != "w1 w2")' 0
quietly easi `M' selvars(w1: q) vce(robust) nolog notable
matrix __d = `D'[1, 1..3]
matrix __e = e(sel_diag)
matrix __e = __e[1, 1..3]
_chk "the probit with q = the estimation's" `=mreldif(__d, __e)' 1e-10
capture easidiag `M' selgoods(w3)
_chk "selgoods(w3), the last good, refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture easidiag `M' selvars(z1)
_chk "selvars(z1), a demographic, refused (rc `=_rc')" `=abs(_rc - 198)' 0
quietly easidiag `M'
_chk "no selection: no r(sel_diag)" `=("`r(selgoods)'" != "")' 0

di as txt _n "{hline 60}"
if ${NFAIL} == 0 di as res "test_step26: all checks passed"
else di as err "test_step26: ${NFAIL} check(s) failed"
