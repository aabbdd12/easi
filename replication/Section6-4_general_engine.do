*! test_step23.do -- easi 2.0.0, phase C step 3a: the general engine.
*! Without selection, the general engine (every equation with its own
*! regressors and instruments, GMM form of the restricted 3SLS, Jacobians by
*! central differences) must give the usual engine's results: the same
*! coefficients and Sigma (the same algebra), the same variance and
*! elasticity standard errors to the precision of the differences.  The
*! developer scalar _easi_geng switches it on.
clear all
set more off
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Section6-4_general_engine.do
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

local HSH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local HPR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local HZ  age hsex carown time tran
local S1 "`HSH', lnprices(`HPR') lnexpenditure(log_y) demographics(`HZ') power(5) py zy pz interpz(`HZ') nolog notable compensated"
local S2 "`HSH', lnprices(`HPR') lnexpenditure(log_y) demographics(age hsex carown) power(3) nolog notable compensated"
local S4 "w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 isMale) power(3) nolog notable vce(svy) compensated"
local S5 "w1 w2 w3 [pw=sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 isMale) power(3) nolog notable vce(cluster psu) compensated"
local S6 "w1 w2 w3 [pw=sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 isMale) power(3) py zy pz nolog notable vce(cluster psu) compensated"
local EX b Sigma
local EV V elast_exp_se elast_price_nc_se elast_price_c_se elast_demo_se elast_exp_mkt_se elast_price_nc_mkt_se
local EP elast_exp elast_price_nc elast_price_c elast_demo elast_exp_mkt

foreach k in 2 4 5 6 1 {
	if inlist(`k', 1, 2) use "`ROOT'/examples/hixdata.dta", clear
	else use "`ROOT'/examples/mex_bench.dta", clear
	capture scalar drop _easi_geng
	quietly easi `S`k''
	local it0 = e(iter)
	foreach m in `EX' `EV' `EP' {
		matrix __R_`m' = e(`m')
	}
	scalar _easi_geng = 1
	quietly easi `S`k''
	scalar drop _easi_geng
	di as txt _n "specification `k' (`e(vce)'), `=e(iter)' iterations (usual engine: `it0')"
	_chk "same iterations" `=abs(e(iter) - `it0')' 0
	local d 0
	foreach m of local EX {
		local d = max(`d', mreldif(e(`m'), __R_`m'))
	}
	_chk "coefficients and Sigma: the same algebra" `d' 1e-9
	local d 0
	foreach m of local EP {
		local d = max(`d', mreldif(e(`m'), __R_`m'))
	}
	_chk "elasticities" `d' 1e-9
	local d 0
	foreach m of local EV {
		local d = max(`d', mreldif(e(`m'), __R_`m'))
	}
	_chk "variance and standard errors (central differences)" `d' 1e-6
}

di as txt _n "{hline 60}"
if ${NFAIL} == 0 di as res "test_step23: all checks passed"
else di as err "test_step23: ${NFAIL} check(s) failed"
