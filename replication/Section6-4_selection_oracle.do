*! test_step24.do -- easi 2.0.0, phase C step 3b: the selection of the
*! buyers (Shonkwiler and Yen), estimation.
*!   1. oracle: two goods, power 1, no interaction, legacy(instr) (the
*!      instrument is y_tilda): one equation, exactly identified, so the GMM
*!      is Stata's ivregress 2sls of w1 on (Phi, Phi y, Phi z, Phi np, phi)
*!      with instruments (Phi, Phi z, Phi np, Phi y_tilda, phi), Phi and phi
*!      from Stata's probit, y the completed implicit utility (closed form
*!      with power 1: y = [y_stone - (1-d) np (b0 + b_z z + a np - delta
*!      phi/(1-Phi)) + a np^2/2] / [1 + (1-d) np b_y]);
*!   2. three goods: symmetry holds, delta estimated, selgoods() changes only
*!      the corrected equation's form.
clear all
set more off
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Section6-4_selection_oracle.do
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
quietly gen double w23 = w2 + w3

di as txt _n "1. Oracle: two goods, power 1, against ivregress 2sls"
quietly easi w1 w23 [pw = sweight], lnprices(lp1 lp3) lnexpenditure(lx) demographics(z1 z2) ///
	power(1) nolog notable legacy(instr) selection tolerance(1e-13) iterate(1000)
tempname B
matrix `B' = e(b)
local b0 = `B'[1, 1]
local by = `B'[1, 2]
local bz1 = `B'[1, 3]
local bz2 = `B'[1, 4]
local a  = `B'[1, 5]
local dl = `B'[1, 6]
_chk "converged" `=1 - e(converged)' 0
* the probit, Phi and phi
quietly gen byte d = w1 > 0
quietly gen double s1 = lp1 - lx
quietly gen double s2 = lp3 - lx
quietly probit d s1 s2 z1 z2 [pw = sweight]
quietly predict double xb, xb
quietly gen double Ph = normal(xb)
quietly gen double ph = normalden(xb)
* y: the completed implicit utility, closed form at the estimate
quietly gen double np = lp1 - lp3
quietly gen double ystone = lx - (w1 * lp1 + w23 * lp3)
quietly gen double mills = ph / normal(-xb)
quietly gen double y = (ystone - (1 - d) * np * (`b0' + `bz1' * z1 + `bz2' * z2 + `a' * np - `dl' * mills) ///
	+ 0.5 * `a' * np^2) / (1 + (1 - d) * np * `by')
* y_tilda: ln x deflated by the Stone index at the weighted mean shares
quietly summarize w1 [aw = sweight]
local m1 = r(mean)
quietly gen double ytil = lx - (`m1' * lp1 + (1 - `m1') * lp3)
foreach v in y z1 z2 np ytil {
	quietly gen double P_`v' = Ph * `v'
}
quietly ivregress 2sls w1 Ph P_z1 P_z2 P_np ph (P_y = P_ytil) [pw = sweight], noconstant
tempname I
matrix `I' = e(b)
* ivregress: P_y Ph P_z1 P_z2 P_np ph ; easi: _cons y z1 z2 np delta
matrix `I' = (`I'[1, 2], `I'[1, 1], `I'[1, 3], `I'[1, 4], `I'[1, 5], `I'[1, 6])
_chk "coefficients and delta = ivregress 2sls" `=mreldif(`I', `B')' 1e-7

di as txt _n "1b. predict after selection, against the same closed forms"
* ivregress has replaced e(): the same easi again
quietly easi w1 w23 [pw = sweight], lnprices(lp1 lp3) lnexpenditure(lx) demographics(z1 z2) ///
	power(1) nolog notable legacy(instr) selection tolerance(1e-13) iterate(1000)
quietly predict double __y, y
quietly gen double __d = reldif(__y, y)
quietly summarize __d
_chk "predict, y: the completed implicit utility" r(max) 1e-9
quietly gen double f1 = `b0' + `by' * y + `bz1' * z1 + `bz2' * z2 + `a' * np
quietly gen double E1 = Ph * f1 + `dl' * ph
quietly predict double __e*, shares
quietly replace __d = reldif(__e1, E1)
quietly summarize __d
_chk "predict, shares: E[w_1] = Phi f + delta phi" r(max) 1e-9
quietly replace __d = abs(__e1 + __e2 - 1)
quietly summarize __d
_chk "the expected shares sum to one" r(max) 1e-12
quietly predict double __f*, shares latent
quietly replace __d = reldif(__f1, f1)
quietly summarize __d
_chk "predict, latent: f" r(max) 1e-9
quietly predict double __r*, residuals
quietly replace __d = reldif(__r1, w1 - E1)
quietly summarize __d
_chk "predict, residuals: w - E[w]" r(max) 1e-9

di as txt _n "1c. The completion of y for the non-buyers (decision DC1)"
* without the completion (the observed zeros), y = y_stone + a np^2/2 with
* power 1 and no interaction; the buyers are not completed
quietly gen double ync = ystone + 0.5 * `a' * np^2
quietly gen double ady = abs(y - ync) if !d
quietly summarize ady [aw = sweight]
local mdy = r(mean)
quietly summarize d [aw = sweight]
local pct = 100 * (1 - r(mean))
_chk "e(sel_ycomp_pct): weighted percent of the non-buyers (`: di %5.2f `pct'')" `=reldif(e(sel_ycomp_pct), `pct')' 1e-10
_chk "e(sel_ycomp_mean): their mean |y - y of the zeros| (`: di %8.5f `mdy'')" `=reldif(e(sel_ycomp_mean), `mdy')' 1e-8
quietly summarize y [aw = sweight]
_chk "e(sel_ycomp_sd): the sd of y (divisor N - 1, as summarize)" `=reldif(e(sel_ycomp_sd), r(sd))' 1e-10

di as txt _n "2. Three goods"
quietly easi w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) ///
	power(3) nolog notable vce(cluster psu) selection
matrix `B' = e(b)
_chk "symmetry: w1:lp2 = w2:lp1" `=abs(_b[w1:lp2] - _b[w2:lp1])' 1e-12
_chk "two delta" `=missing(_b[w1:delta]) + missing(_b[w2:delta])' 0
_chk "the SE of delta positive" `=(_se[w1:delta] <= 0) + (_se[w2:delta] <= 0)' 0
quietly easi w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) ///
	power(3) nolog notable vce(cluster psu) selgoods(w2)
capture local x = _b[w1:delta]
_chk "selgoods(w2): no delta for w1" `=(_rc == 0)' 0
_chk "selgoods(w2): delta for w2" `=missing(_b[w2:delta])' 0
_chk "symmetry still" `=abs(_b[w1:lp2] - _b[w2:lp1])' 1e-12
di as txt _n "3. estat engel after selection"
quietly easi w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) ///
	power(3) nolog notable vce(cluster psu) selection
tempname VS
matrix `VS' = e(V_sel)
_chk "e(V_sel): the coefficients and the probits (18 + 2 x 6)" `=abs(rowsof(`VS') - 30)' 0
matrix __v = `VS'[1..18, 1..18]
_chk "e(V_sel) contains e(V)" `=mreldif(__v, e(V))' 1e-12
capture noisily estat engel, n(12) nodraw data("`c(tmpdir)'/__eng_sel")
_chk "estat engel at the means (rc `=_rc')" `=_rc != 0' 0
preserve
quietly use "`c(tmpdir)'/__eng_sel", clear
quietly ds
di as txt "  curve variables: `r(varlist)'"
restore
quietly predict double __s*, shares
capture noisily estat engel, asobserved n(12) nodraw data("`c(tmpdir)'/__eng_sel2")
_chk "estat engel as observed (rc `=_rc')" `=_rc != 0' 0

quietly easi w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) ///
	power(3) nolog notable vce(cluster psu)
capture predict double __l*, shares latent
_chk "predict, latent refused without selection (rc `=_rc')" `=abs(_rc - 198)' 0

di as txt _n "{hline 60}"
if ${NFAIL} == 0 di as res "test_step24: all checks passed"
else di as err "test_step24: ${NFAIL} check(s) failed"
