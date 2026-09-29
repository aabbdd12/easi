*! test_step21.do -- easi 2.0.0, phase C step 1: missing prices, pimpute().
*!   1. oracle: the prices filled by hand (weighted mean of the log prices of
*!      the same PSU, then of the same stratum) give the same point estimates
*!      as pimpute(psu strata): e(b), elasticities, N (not the same variance:
*!      the analytic standard errors of pimpute() include the imputation);
*!   1b. the influence of the imputation sums to zero within each group: with
*!      pimpute(psu) under vce(cluster psu) it cancels, and e(V) and every
*!      standard error equal those of the prices filled by hand at the PSU
*!      level (iterated, the households of PSUs without a donor leaving);
*!   1c. under vce(robust) it does not cancel: e(V) differs;
*!   2. predict and estat engel after pimpute() see the same prices as the
*!      estimation (against the same commands on the data filled by hand);
*!   3. easidiag with pimpute() = easidiag on the data filled by hand;
*!   4. with no missing price, pimpute() changes nothing; the weight of the
*!      imputation is that of the estimation (hhsize() included).
clear all
set more off
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Section6-1_pimpute.do
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
capture program drop _same
program define _same
	args lab items
	local d 0
	foreach m of local items {
		local d = max(`d', mreldif(e(`m'), __R_`m'))
	}
	_chk "`lab'" `d' 1e-10
end
local IT b V elast_exp elast_exp_se elast_price_nc elast_price_nc_se elast_demo elast_exp_mkt elast_price_nc_mkt elast_price_nc_mkt_se
local ITp b elast_exp elast_price_nc elast_demo elast_exp_mkt elast_price_nc_mkt

* the Mexican cereals; the prices of good 1 missing for 30% of the households
* (all of some PSUs, so that the stratum level is used too), of good 2 for 10%
use "`ROOT'/examples/mex_bench.dta", clear
set seed 2109
quietly bysort psu: gen double __u = runiform() if _n == 1
quietly bysort psu: replace __u = __u[1]
quietly replace lp1 = . if __u < 0.08 | runiform() < 0.25
quietly replace lp2 = . if runiform() < 0.10
drop __u
quietly count if missing(lp1)
local m1 = r(N)
quietly count if missing(lp2)
di as txt "missing prices: good 1 " `m1' ", good 2 " r(N)
tempfile raw
quietly save "`raw'"

* the prices filled by hand: PSU first, then stratum, weighted by sweight;
* the donors are the households whose price is OBSERVED, at every level
foreach v in lp1 lp2 {
	tempvar o
	quietly gen double `o' = `v'
	foreach g in psu strata {
		tempvar num den
		quietly egen double `num' = total(cond(!missing(`o'), sweight * `o', .)), by(`g')
		quietly egen double `den' = total(cond(!missing(`o'), sweight, .)), by(`g')
		quietly replace `v' = `num' / `den' if missing(`v') & `den' > 0 & !missing(`den')
		drop `num' `den'
	}
	drop `o'
}
quietly count if missing(lp1) | missing(lp2)
_chk "every price filled by hand (no household left)" r(N) 0
tempfile filled
quietly save "`filled'"

local M "lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) py nolog notable"

di as txt _n "1. Oracle: pimpute(psu strata) = the prices filled by hand"
use "`filled'", clear
quietly easi w1 w2 w3 [pw = sweight], `M' vce(cluster psu)
foreach m of local IT {
	matrix __R_`m' = e(`m')
}
local N0 = e(N)
quietly predict double __f*, shares
quietly estat engel, n(15) nodraw data("`c(tmpdir)'/__eng0")
keep __f1 __f2 __f3
quietly gen long __id = _n
tempfile f0
quietly save "`f0'"

use "`raw'", clear
quietly easi w1 w2 w3 [pw = sweight], `M' vce(cluster psu) pimpute(psu strata)
_same "e(b) and the elasticities (point estimates)" "`ITp'"
_chk "the same N (`=e(N)' vs `N0')" `=abs(e(N) - `N0')' 0
_chk "e(pimpute) stored" `=("`e(pimpute)'" != "psu strata")' 0

di as txt _n "2. predict and estat engel after pimpute()"
quietly predict double __g*, shares
quietly gen long __id = _n
quietly merge 1:1 __id using "`f0'", nogenerate
quietly gen double __d = max(reldif(__f1, __g1), reldif(__f2, __g2), reldif(__f3, __g3))
quietly summarize __d if e(sample)
_chk "predict, shares: the same fitted shares on e(sample)" r(max) 1e-10
quietly count if e(sample) & missing(__g1)
_chk "no missing prediction on e(sample)" r(N) 0
quietly estat engel, n(15) nodraw data("`c(tmpdir)'/__eng1")
preserve
quietly use "`c(tmpdir)'/__eng1", clear
* the curves; their bands (_se, _lo, _hi) carry the analytic standard errors,
* which include the imputation under pimpute() and not with the prices filled
* by hand
quietly ds _se* _lo* _hi*, not
local cv `r(varlist)'
mkmat `cv', matrix(__E1)
quietly use "`c(tmpdir)'/__eng0", clear
mkmat `cv', matrix(__E0)
restore
_chk "estat engel: the same curves" `=mreldif(__E1, __E0)' 1e-10

di as txt _n "3. easidiag with pimpute()"
use "`filled'", clear
quietly easidiag w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3)
local SC : r(scalars)
foreach x of local SC {
	local R_`x' = r(`x')
}
use "`raw'", clear
quietly easidiag w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) pimpute(psu strata)
local d 0
foreach x of local SC {
	local d = max(`d', reldif(r(`x'), `R_`x''))
}
_chk "easidiag pimpute() = easidiag on the prices filled by hand" `d' 1e-10

di as txt _n "1b. The influence of the imputation cancels within the clusters of its first level"
* the prices filled by hand at the PSU level only, iterated as pimpute(psu):
* the households still without a price leave, and the means are computed
* again on the households that stay
use "`raw'", clear
tempvar o1 o2
quietly gen double `o1' = lp1
quietly gen double `o2' = lp2
local more 1
while `more' {
	quietly replace lp1 = `o1'
	quietly replace lp2 = `o2'
	foreach k in 1 2 {
		tempvar num den
		quietly egen double `num' = total(cond(!missing(`o`k''), sweight * `o`k'', .)), by(psu)
		quietly egen double `den' = total(cond(!missing(`o`k''), sweight, .)), by(psu)
		quietly replace lp`k' = `num' / `den' if missing(lp`k') & `den' > 0 & !missing(`den')
		drop `num' `den'
	}
	quietly count if missing(lp1) | missing(lp2)
	local more = (r(N) > 0)
	quietly drop if missing(lp1) | missing(lp2)
}
quietly easi w1 w2 w3 [pw = sweight], `M' vce(cluster psu)
foreach m of local IT {
	matrix __R_`m' = e(`m')
}
local N0 = e(N)
use "`raw'", clear
quietly easi w1 w2 w3 [pw = sweight], `M' vce(cluster psu) pimpute(psu)
_same "pimpute(psu), vce(cluster psu): e(b), e(V), elasticities and SEs = filled by hand" "`IT'"
_chk "the same N (`=e(N)' vs `N0')" `=abs(e(N) - `N0')' 0

di as txt _n "1c. Under vce(robust) the influence of the imputation is there"
quietly easi w1 w2 w3 [pw = sweight], `M' vce(robust) pimpute(psu strata)
matrix __V1 = e(V)
use "`filled'", clear
quietly easi w1 w2 w3 [pw = sweight], `M' vce(robust)
* scale-free: max |V1 - V0| / sqrt(V0_ii V0_jj)
mata: V0 = st_matrix("e(V)"); V1 = st_matrix("__V1"); d = sqrt(diagonal(V0))
mata: st_numscalar("__d21", max(abs(V1 - V0) :/ (d * d')))
di as txt "  vce(robust): the imputation moves e(V) by " as res %6.4f __d21 as txt " (scale-free)"
_chk "vce(robust): e(V) differs from that of the prices filled by hand" `=(__d21 < 1e-4)' 0

di as txt _n "4. No missing price; the weight of the imputation"
use "`ROOT'/examples/mex_bench.dta", clear
quietly easi w1 w2 w3, `M'
foreach m of local IT {
	matrix __R_`m' = e(`m')
}
quietly easi w1 w2 w3, `M' pimpute(psu)
_same "no missing price: pimpute() changes nothing" "`IT'"
* individuals: the imputation is weighted by the weight x hhsize()
use "`raw'", clear
quietly easi w1 w2 w3 [pw = sweight], `M' vce(cluster psu) pimpute(psu strata) hhsize(hhsize)
foreach m in b elast_exp {
	matrix __R_`m' = e(`m')
}
foreach v in lp1 lp2 {
	tempvar o
	quietly gen double `o' = `v'
	foreach g in psu strata {
		tempvar num den
		quietly egen double `num' = total(cond(!missing(`o'), sweight * hhsize * `o', .)), by(`g')
		quietly egen double `den' = total(cond(!missing(`o'), sweight * hhsize, .)), by(`g')
		quietly replace `v' = `num' / `den' if missing(`v') & `den' > 0 & !missing(`den')
		drop `num' `den'
	}
	drop `o'
}
quietly easi w1 w2 w3 [pw = sweight], `M' vce(cluster psu) hhsize(hhsize)
_same "individuals: filled with the weight x hhsize()" "b elast_exp"

di as txt _n "{hline 60}"
if ${NFAIL} == 0 di as res "test_step21: all checks passed"
else di as err "test_step21: ${NFAIL} check(s) failed"
