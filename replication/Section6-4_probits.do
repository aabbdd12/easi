*! test_step22.do -- easi 2.0.0, phase C step 2: the probits of buying.
*!   1. oracle: e(sel_alpha) = Stata's probit of buying on (ln p - ln x, z)
*!      [pw], good by good, with and without a variable of the probit only;
*!   2. the diagnostics (percent of buyers) and the default goods corrected;
*!   3. the system is corrected: one delta per corrected good;
*!   4. refusals (the rules of duvm and equaids).
clear all
set more off
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Section6-4_probits.do
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
local M "w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) nolog notable vce(cluster psu)"
quietly gen byte d1 = w1 > 0
quietly gen byte d2 = w2 > 0
forvalues j = 1/3 {
	quietly gen double s`j' = lp`j' - lx
}

di as txt _n "1. The probits against Stata's probit"
quietly easi `M'
tempname B0
matrix `B0' = e(b)
quietly easi `M' selection
tempname A D
matrix `A' = e(sel_alpha)
matrix `D' = e(sel_diag)
_chk "default goods corrected: w1 w2 (e(selgoods) = `e(selgoods)')" `=("`e(selgoods)'" != "w1 w2")' 0
_chk "e(selection) = shonkwiler-yen" `=("`e(selection)'" != "shonkwiler-yen")' 0
forvalues g = 1/2 {
	quietly probit d`g' s1 s2 s3 z1 z2 [pw = sweight]
	tempname P
	matrix `P' = e(b)
	* easi's order: _cons, the prices, the demographics
	* Stata's order: s1 s2 s3 z1 z2 _cons
	matrix `P' = (`P'[1, 6], `P'[1, 1..5])'
	matrix __a = `A'[1..6, `g']
	_chk "good `g': probit coefficients = Stata's probit" `=mreldif(__a, `P')' 1e-7
	quietly summarize d`g' [aw = sweight]
	_chk "good `g': percent of buyers" `=reldif(`D'[`g', 1], 100 * r(mean))' 1e-12
}
* s3 = -lx: the probit is homogeneous of degree zero in prices and expenditure

di as txt _n "2. A variable of the probit only (rururb, good 1)"
quietly easi `M' selection selvars(w1: rururb)
matrix `A' = e(sel_alpha)
quietly probit d1 s1 s2 s3 z1 z2 rururb [pw = sweight]
tempname P
matrix `P' = e(b)
matrix `P' = (`P'[1, 7], `P'[1, 1..6])'
matrix __a = `A'[1..7, 1]
_chk "good 1 with rururb: = Stata's probit" `=mreldif(__a, `P')' 1e-7
_chk "good 2 has no rururb (missing coefficient)" `=!missing(`A'[7, 2])' 0
quietly easi `M' selgoods(w2)
_chk "selgoods(w2): only good 2 corrected" `=("`e(selgoods)'" != "w2")' 0
_chk "good 1 not corrected (missing column)" `=!missing(e(sel_alpha)[1, 1])' 0

di as txt _n "3. The system is corrected (step 3): one delta per corrected good"
quietly easi `M' selection
_chk "e(b) has the two delta beyond the latent coefficients" `=colsof(e(b)) - colsof(`B0') - 2' 0

di as txt _n "4. Refusals"
capture easi `M' selgoods(w3)
_chk "selgoods(): the last good refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture easi `M' selgoods(w9)
_chk "selgoods(): an unknown good refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture easi `M' selvars(w1: z1)
_chk "selvars(): a model variable refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture easi `M' selvars(w1 w2: rururb)
_chk "selvars(): two goods in a segment refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture easi `M' selgoods(w2) selvars(w1: rururb)
_chk "selvars(): a good not corrected refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture easi `M' selvars(w1: nosuchvar)
_chk "selvars(): an unknown variable refused (rc `=_rc')" `=abs(_rc - 111)' 0
capture easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) nolog vce(conventional) selection
_chk "selection with vce(conventional) refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) nolog compat selection
_chk "selection with compat refused (rc `=_rc')" `=abs(_rc - 198)' 0
* the last good with zeros
quietly gen double w3z = w3
quietly replace w3z = 0 in 1
capture easi w1 w2 w3z [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) nolog selection
_chk "the last good with zeros refused (rc `=_rc')" `=abs(_rc - 198)' 0

di as txt _n "{hline 60}"
if ${NFAIL} == 0 di as res "test_step22: all checks passed"
else di as err "test_step22: ${NFAIL} check(s) failed"
