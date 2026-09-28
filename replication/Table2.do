*! Table2.do -- Table 2 of the technical note: the three types of
*! elasticities -- the mean over the households of the elasticities of each
*! household, weighted by the weight (households), the weight times the size
*! of the household (individuals) or the weight times its expenditure
*! (market) -- on the reduced Mexican survey, vce(svy).
*! Also the numbers of Section 4.3: how far the reference household moves
*! with the mean taken (ln of the mean of x against the mean of ln x), and
*! what the compensated elasticities change when the mean of the products of
*! the shares replaces the product of their means (cov(w_j, w_k)/mean w_j).

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table2.do
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
capture log close tt
log using "`OUT'/types.log", replace text name(tt)

*------------------------------------------------ the three types, mex_bench
use "`ROOT'/examples/mex_bench.dta", clear
local M "w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) py vce(svy) nolog notable"
quietly easi `M'
tempname EH EHs PH PHs EM EMs PM PMs EI EIs PI PIs
matrix `EH'  = e(elast_exp)
matrix `EHs' = e(elast_exp_se)
matrix `PH'  = vecdiag(e(elast_price_nc))
matrix `PHs' = vecdiag(e(elast_price_nc_se))
matrix `EM'  = e(elast_exp_mkt)
matrix `EMs' = e(elast_exp_mkt_se)
matrix `PM'  = vecdiag(e(elast_price_nc_mkt))
matrix `PMs' = vecdiag(e(elast_price_nc_mkt_se))
local chk1 = e(chk_engel) + e(chk_cournot) + e(chk_engel_mkt) + e(chk_cournot_mkt)
quietly easi `M' hhsize(hhsize)
matrix `EI'  = e(elast_exp)
matrix `EIs' = e(elast_exp_se)
matrix `PI'  = vecdiag(e(elast_price_nc))
matrix `PIs' = vecdiag(e(elast_price_nc_se))
local chk2 = e(chk_engel) + e(chk_cournot)

di as txt _n "Table 2. Expenditure and own-price elasticities of the three types,"
di as txt "reduced Mexican survey (2,477 households), power(3) py, vce(svy);"
di as txt "standard errors in parentheses."
di as txt "{hline 78}"
di as txt "good" _col(12) "type" _col(26) "expenditure" _col(46) "own price"
di as txt "{hline 78}"
local gn "corn wheat other"
forvalues j = 1/3 {
	local g : word `j' of `gn'
	foreach t in H I M {
		local lab = cond("`t'" == "H", "households", cond("`t'" == "I", "individuals", "market"))
		di as txt "`g'" _col(12) "`lab'" _col(24) as res %8.4f `E`t''[1, `j'] ///
			as txt " (" as res %6.4f `E`t's'[1, `j'] as txt ")" _col(44) ///
			as res %8.4f `P`t''[1, `j'] as txt " (" as res %6.4f `P`t's'[1, `j'] as txt ")"
		local g ""
	}
}
di as txt "{hline 78}"
di as txt "Engel and Cournot residuals, all types: " as res %9.2e `chk1' + `chk2'

*------------------------------------------------ the reference household
di as txt _n "ln(mean x) - mean(ln x): how far the reference household moves"
quietly gen double __x = exp(lx)
quietly summarize lx [aw = sweight]
local ml = r(mean)
quietly summarize __x [aw = sweight]
local d = ln(r(mean)) - `ml'
di as txt "  mex_bench (weighted):  " as res %6.3f `d' as txt "  (the mean of x is " ///
	as res %4.1f 100 * (exp(`d') - 1) as txt "% above exp(mean ln x))"
use "`ROOT'/examples/hixdata.dta", clear
quietly gen double __x = exp(log_y)
quietly summarize log_y
local ml = r(mean)
quietly summarize __x
local d = ln(r(mean)) - `ml'
di as txt "  hixdata:               " as res %6.3f `d' as txt "  (the mean of x is " ///
	as res %4.1f 100 * (exp(`d') - 1) as txt "% above exp(mean ln x))"

*------------------------------------------------ mean of products, product of means
* the compensated elasticity of the households is (mean Psi + mean(w_j w_k))
* / mean w_j - 1{j = k}; the product of the means changes it by
* cov(w_j, w_k) / mean w_j
local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
mata: W = st_data(., tokens("`SH'")); m = mean(W); C = quadvariance(W) :* ((rows(W) - 1) / rows(W))
mata: D = C :/ m'
mata: st_numscalar("__omin", min(diagonal(D))); st_numscalar("__omax", max(diagonal(D)))
mata: st_numscalar("__cmax", max(abs(D - diag(diagonal(D)))))
di as txt _n "cov(w_j, w_k) / mean w_j, hixdata: own-price " as res %6.3f __omin ///
	as txt " to " as res %6.3f __omax as txt ", cross-price up to " as res %6.3f __cmax
use "`ROOT'/examples/mex_bench.dta", clear
mata: W = st_data(., ("w1", "w2", "w3")); ww = st_data(., "sweight"); ww = ww :/ mean(ww)
mata: m = mean(W, ww); C = quadcross(W :- m, ww, W :- m) :/ rows(W); D = C :/ m'
mata: st_numscalar("__omin", min(diagonal(D))); st_numscalar("__omax", max(diagonal(D)))
di as txt "cov(w_j, w_k) / mean w_j, mex_bench: own-price " as res %6.3f __omin ///
	as txt " to " as res %6.3f __omax
log close tt
di as res _n "Table2: done"
