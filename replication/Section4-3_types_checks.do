*! test_step20.do -- easi 2.0.0, phase D: the types of elasticities.
*!   1. market: the elasticities of the market are the household formulas with
*!      every mean weighted by the weight times the total expenditure. Oracle:
*!      the point formulas (_easi_epoint) re-run from outside on e(b), the
*!      implicit utility of predict and the weights w x exp(ln x), against
*!      e(elast_*_mkt); the same route with the weights w gives e(elast_*),
*!      which checks the route itself. Specifications: hixdata with all the
*!      interactions (power 5, py zy pz), mex_bench [pw] vce(svy).
*!   2. the self-test of the influence function (numerators rebuilt = point
*!      estimates) and the aggregation identities, for both types;
*!   3. replay: elasticities(market) after a households estimation, and back;
*!      refusals (market under compat, an unknown type, individuals on replay).
clear all
set more off
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Section4-3_types_checks.do
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
run "`ROOT'/replication/dgp_easi.do"

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

* The Mata functions of an ado-file are private to it: compile a copy of the
* Mata block of easi.ado (from the line "mata:" to the end) for the oracle.
tempfile mcopy
mata: _L = cat("`ROOT'/src/easi.ado"); _i = selectindex(strtrim(_L) :== "mata:")[1]
mata: _fh = fopen(st_local("mcopy"), "w"); for (_j = _i; _j <= rows(_L); _j++) fput(_fh, _L[_j]); fclose(_fh)
quietly run "`mcopy'"

* the point elasticities of e() re-computed from outside, with the weights wv
capture program drop _oracle
program define _oracle
	args shares lnp demo wv sfx
	tempvar yh es
	quietly predict double `yh' if e(sample), y
	quietly gen byte `es' = e(sample)
	local ipz ""
	if e(pz) {
		foreach v in `e(interpz)' {
			local ipz `ipz' `: list posof "`v'" in demo'
		}
	}
	mata: _orc("`shares'", "`lnp'", "`demo'", "`yh'", "`wv'", "`ipz'", "`sfx'", "`es'")
end
mata:
void _orc(string scalar sv, string scalar pv, string scalar zv, string scalar yv,
	string scalar wv, string scalar sipz, string scalar sfx, string scalar es)
{
	real matrix s, P, z, EPRICE, EP, EPS, EPQ, EZ
	real colvector y, wt, b, Dn2
	real rowvector EI, ER, ws, ipz
	real scalar J
	s = st_data(., sv, es)
	P = st_data(., pv, es)
	z = st_data(., zv, es)
	y = st_data(., yv, es)
	wt = st_data(., wv, es)
	b = st_matrix("e(b)")'
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))
	J = cols(s)
	_easi_epoint(b, s, P, z, y, wt, st_numscalar("e(power)"),
		st_numscalar("e(py)"), st_numscalar("e(zy)"), st_numscalar("e(pz)"),
		ipz, 0, 0, EI, ER, EPRICE, EP, EPS, EPQ, EZ, ws, Dn2)
	st_numscalar("__d1", mreldif(EI, st_matrix("e(elast_exp" + sfx + ")")))
	st_numscalar("__d2", mreldif(EPRICE', st_matrix("e(elast_price_nc" + sfx + ")")))
	st_numscalar("__d3", mreldif(EPQ - I(J), st_matrix("e(elast_price_c" + sfx + ")")))
	st_numscalar("__d4", mreldif(EZ, st_matrix("e(elast_demo" + sfx + ")")))
}
end

local HSH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local HPR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local HZ  age hsex carown time tran

di as txt _n "1. market = the household formulas weighted by w x exp(ln x)"
forvalues k = 1/2 {
	if `k' == 1 {
		use "`ROOT'/examples/hixdata.dta", clear
		local SH `HSH'
		local PR `HPR'
		local Z  `HZ'
		local lx log_y
		quietly easi `SH', lnprices(`PR') lnexpenditure(`lx') demographics(`Z') ///
			power(5) py zy pz nolog tolerance(1e-13) iterate(500) notable
		quietly gen double __w = 1
	}
	else {
		use "`ROOT'/examples/mex_bench.dta", clear
		local SH w1 w2 w3
		local PR lp1 lp2 lp3
		local Z  z1 isMale
		local lx lx
		quietly easi `SH', lnprices(`PR') lnexpenditure(`lx') demographics(`Z') ///
			power(3) py nolog vce(svy) tolerance(1e-13) iterate(500) notable
		quietly summarize sweight if e(sample), meanonly
		quietly gen double __w = sweight * r(N) / r(sum)
	}
	local lab = cond(`k' == 1, "hixdata, all interactions", "mex_bench, vce(svy)")
	_chk "`lab': converged" `=1 - e(converged)' 0
	_oracle "`SH'" "`PR'" "`Z'" __w ""
	_chk "`lab': route check, households = e(elast_exp)" __d1 1e-8
	_chk "`lab': route check, households = e(elast_price_nc)" __d2 1e-8
	_chk "`lab': route check, households = e(elast_price_c)" __d3 1e-8
	_chk "`lab': route check, households = e(elast_demo)" __d4 1e-8
	quietly gen double __wx = __w * exp(`lx')
	_oracle "`SH'" "`PR'" "`Z'" __wx "_mkt"
	_chk "`lab': market expenditure elasticities" __d1 1e-8
	_chk "`lab': market uncompensated price elasticities" __d2 1e-8
	_chk "`lab': market compensated price elasticities" __d3 1e-8
	_chk "`lab': market demographic semi-elasticities" __d4 1e-8
	* the weights matter: market differs from households
	_chk "`lab': market differs from households" `=(mreldif(e(elast_exp), e(elast_exp_mkt)) < 1e-4)' 0
}

di as txt _n "2. Self-tests and identities, both types"
use "`ROOT'/examples/hixdata.dta", clear
quietly easi `HSH', lnprices(`HPR') lnexpenditure(log_y) demographics(age hsex carown) power(3) nolog compensated notable
_chk "numerators of the influence function = estimates, households" scalar(_easi_esampchk) 1e-12
_chk "numerators of the influence function = estimates, market" scalar(_easi_esampchk_mkt) 1e-12
_chk "Engel, households" e(chk_engel) 1e-8
_chk "Cournot, households" e(chk_cournot) 1e-8
_chk "Engel, market" e(chk_engel_mkt) 1e-8
_chk "Cournot, market" e(chk_cournot_mkt) 1e-8
_chk "e(elasticities) = households by default" `=("`e(elasticities)'" != "households")' 0
foreach m in elast_exp_mkt elast_exp_mkt_se elast_price_nc_mkt elast_price_nc_mkt_se elast_price_c_mkt elast_price_c_mkt_se elast_demo_mkt elast_demo_mkt_se {
	capture confirm matrix e(`m')
	_chk "e(`m') stored" `=_rc != 0' 0
}

di as txt _n "3. Replay and refusals"
tempname M H
matrix `M' = e(elast_exp_mkt)
matrix `H' = e(elast_exp)
capture noisily easi, elasticities(market) notable
_chk "replay elasticities(market) (rc `=_rc')" `=_rc != 0' 0
_chk "replay leaves e() unchanged" `=mreldif(`M', e(elast_exp_mkt)) + mreldif(`H', e(elast_exp))' 0
quietly easi `HSH', lnprices(`HPR') lnexpenditure(log_y) demographics(age hsex carown) power(3) nolog elasticities(mkt) notable
_chk "estimation with elasticities(mkt): e(elasticities) = market" `=("`e(elasticities)'" != "market")' 0
_chk "the same numbers as the households estimation" `=mreldif(`M', e(elast_exp_mkt)) + mreldif(`H', e(elast_exp))' 1e-12
capture easi, elasticities(households) notable
_chk "replay back to households (rc `=_rc')" `=_rc != 0' 0
capture easi `HSH', lnprices(`HPR') lnexpenditure(log_y) demographics(age hsex carown) power(3) nolog elasticities(market) compat
_chk "market under compat refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture easi `HSH', lnprices(`HPR') lnexpenditure(log_y) demographics(age hsex carown) power(3) nolog legacy(elast) elasticities(market)
_chk "market under legacy(elast) refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture easi `HSH', lnprices(`HPR') lnexpenditure(log_y) demographics(age hsex carown) power(3) nolog elasticities(average)
_chk "unknown type refused (rc `=_rc')" `=abs(_rc - 198)' 0
quietly easi `HSH', lnprices(`HPR') lnexpenditure(log_y) demographics(age hsex carown) power(3) nolog compat notable
capture confirm matrix e(elast_exp_mkt)
_chk "compat stores no market elasticities" `=_rc == 0' 0
capture easi, elasticities(market)
_chk "replay market after compat refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture easi, elasticities(individuals)
_chk "replay individuals after a households estimation refused (rc `=_rc')" `=abs(_rc - 198)' 0

di as txt _n "4. individuals = households with the weight x hhsize() in the whole estimation"
local IT b V Sigma elast_exp elast_exp_se elast_price_nc elast_price_nc_se elast_price_c elast_price_c_se elast_demo elast_demo_se slutsky
capture program drop _same
program define _same
	args lab items
	local d 0
	foreach m of local items {
		local d = max(`d', mreldif(e(`m'), __R_`m'))
	}
	_chk "`lab'" `d' 1e-10
end
* hixdata, no weight: hhsize(n) = [aw = n]
use "`ROOT'/examples/hixdata.dta", clear
quietly gen byte n = 1 + mod(_n, 5)
local H3 "lnprices(`HPR') lnexpenditure(log_y) demographics(age hsex carown) power(3) nolog compensated notable"
quietly easi `HSH' [aw = n], `H3'
foreach m of local IT {
	matrix __R_`m' = e(`m')
}
quietly easi `HSH', `H3' hhsize(n)
_chk "hhsize() alone gives individuals" `=("`e(elasticities)'" != "individuals")' 0
_same "hixdata: hhsize(n) = [aw = n], households (all of e())" "`IT'"
_chk "e(wtype) e(wexp) = aweight = n" `=("`e(wtype)'`e(wexp)'" != "aweight= n")' 0
_chk "e(hhsize) = n" `=("`e(hhsize)'" != "n")' 0
capture confirm matrix e(elast_exp_mkt)
_chk "no market elasticities after individuals" `=_rc == 0' 0
quietly easi `HSH', `H3' elasticities(individuals) hhsize(n)
_same "elasticities(individuals) hhsize(n): the same" "`IT'"
capture noisily easi, elasticities(market)
_chk "replay market after individuals refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture noisily easi, elasticities(households)
_chk "replay households after individuals refused (rc `=_rc')" `=abs(_rc - 198)' 0
capture noisily easi, elasticities(individuals) notable
_chk "replay individuals after individuals (rc `=_rc')" `=_rc != 0' 0
capture noisily estat engel, n(20)
_chk "estat engel after individuals (rc `=_rc')" `=_rc != 0' 0
* an explicit type wins: hhsize() not used
quietly easi `HSH', `H3'
foreach m of local IT {
	matrix __R_`m' = e(`m')
}
quietly easi `HSH', `H3' elasticities(households) hhsize(n)
_same "elasticities(households) hhsize(n) = no hhsize()" "`IT'"
* refusals
capture easi `HSH', `H3' elasticities(individuals)
_chk "individuals without hhsize() refused (rc `=_rc')" `=abs(_rc - 198)' 0
quietly replace n = 0 in 3
capture easi `HSH', `H3' hhsize(n)
_chk "hhsize() with a zero refused (rc `=_rc')" `=abs(_rc - 411)' 0
quietly replace n = 2 in 3
capture easi `HSH', `H3' hhsize(n) compat
_chk "individuals under compat refused (rc `=_rc')" `=abs(_rc - 198)' 0
* mex_bench, vce(svy): the svyset weight x hhsize() = [pw = sweight x n]
use "`ROOT'/examples/mex_bench.dta", clear
quietly gen byte n = 1 + mod(_n, 4)
quietly gen double swn = sweight * n
local M3 "lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 isMale) power(3) py nolog vce(svy) compensated notable"
quietly easi w1 w2 w3 [pw = swn], `M3'
foreach m of local IT {
	matrix __R_`m' = e(`m')
}
quietly easi w1 w2 w3, `M3' hhsize(n)
_same "mex_bench vce(svy): hhsize(n) = [pw = sweight x n] (all of e())" "`IT'"
_chk "e(wexp) = (sweight) * (n)" `=("`e(wexp)'" != "= (sweight) * (n)")' 0

di as txt _n "5. easidiag: weights, and the weights of individuals"
* 1.0.0 stopped on every weighted diagnostic (abs(= w)): fixed in 2.0.0
use "`ROOT'/examples/mex_bench.dta", clear
quietly gen byte n = 1 + mod(_n, 4)
local D "lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 isMale) power(3)"
capture noisily easidiag w1 w2 w3 [pw = sweight], `D'
_chk "easidiag [pw] runs (rc `=_rc')" `=_rc != 0' 0
quietly easidiag w1 w2 w3 [aw = n], `D'
local SC : r(scalars)
foreach x of local SC {
	local R_`x' = r(`x')
}
quietly easidiag w1 w2 w3, `D' hhsize(n)
local d 0
foreach x of local SC {
	local d = max(`d', reldif(r(`x'), `R_`x''))
}
_chk "easidiag hhsize(n) = easidiag [aw = n] (all r() scalars)" `d' 1e-10
quietly easidiag w1 w2 w3, `D' hhsize(n) elasticities(households)
local d 0
quietly easidiag w1 w2 w3, `D'
foreach x of local SC {
	local R_`x' = r(`x')
}
quietly easidiag w1 w2 w3, `D' hhsize(n) elasticities(households)
foreach x of local SC {
	local d = max(`d', reldif(r(`x'), `R_`x''))
}
_chk "easidiag elasticities(households) hhsize(n) = no weight" `d' 1e-10
capture easidiag w1 w2 w3, `D' elasticities(individuals)
_chk "easidiag individuals without hhsize() refused (rc `=_rc')" `=abs(_rc - 198)' 0

di as txt _n "6. The y term of the elasticity Jacobian: analytic total = central differences with y(beta)"
* y = (y_stone + p'A(z)p/2)/(1 - p'Bp/2) depends on the coefficients; the
* reported Jacobian is _easi_Gan (y fixed) + _easi_Gy (the y term), left in
* memory by the developer scalar _easi_jaccheck
mata:
real colvector _t6_e(real colvector b, real matrix s, real matrix P, real matrix z,
	real colvector lnx, real colvector wt, real scalar R, real scalar py,
	real scalar zy, real scalar pz, real rowvector ipz)
{
	real matrix np, EPRICE, EP, EPS, EPQ, EZ
	real colvector pAp, pBp, y, Dn2
	real rowvector EI, ER, ws
	np = P[, 1::(cols(s) - 1)] :- P[, cols(s)]
	_easi_quadf(b, np, z, R, py, pz, zy, ipz, 0, pAp, pBp)
	y = (lnx - rowsum(s :* P) :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)
	_easi_epoint(b, s, P, z, y, wt, R, py, zy, pz, ipz, 0, 0,
		EI, ER, EPRICE, EP, EPS, EPQ, EZ, ws, Dn2)
	return(_easi_estack(EI, EPRICE, EZ, EPQ, 1))
}
void _t6(string scalar es, string scalar sipz)
{
	real matrix s, P, z, Gt, Ga
	real colvector b, lnx, wt, bp, bm
	real rowvector ipz
	real scalar R, py, pz, zy, j, h
	s = st_data(., "w1 w2 w3", es); P = st_data(., "lp1 lp2 lp3", es)
	z = st_data(., "z1 z2", es); lnx = st_data(., "lx", es)
	wt = J(rows(s), 1, 1); b = st_matrix("e(b)")'
	R = st_numscalar("e(power)"); py = st_numscalar("e(py)")
	pz = st_numscalar("e(pz)"); zy = st_numscalar("e(zy)")
	ipz = strtoreal(tokens(sipz))
	Gt = J(3 + 2 * 9 + 6, rows(b), .)
	for (j = 1; j <= rows(b); j++) {
		h = 1e-5 * max((1, abs(b[j])))
		bp = b; bp[j] = bp[j] + h
		bm = b; bm[j] = bm[j] - h
		Gt[, j] = (_t6_e(bp, s, P, z, lnx, wt, R, py, zy, pz, ipz)
			- _t6_e(bm, s, P, z, lnx, wt, R, py, zy, pz, ipz)) :/ (2 * h)
	}
	Ga = st_matrix("_easi_Gan") + st_matrix("_easi_Gy")
	st_numscalar("__d6", max(abs(Ga - Gt)) / max(abs(Gt)))
	st_numscalar("__d6y", max(abs(st_matrix("_easi_Gy"))) / max(abs(Gt)))
}
end
easi_dgp, n(500) seed(5) power(3)
scalar _easi_jaccheck = 1
quietly easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) ///
	power(3) py zy pz nolog compensated notable tolerance(1e-13) iterate(500)
scalar drop _easi_jaccheck
gen byte es = e(sample)
mata: _t6("es", "1 2")
_chk "power 3, py zy pz: analytic total Jacobian = central differences" __d6 1e-6
_chk "the y term is not zero here" `=(__d6y < 1e-6)' 0

di as txt _n "{hline 60}"
if ${NFAIL} == 0 di as res "test_step20: all checks passed"
else di as err "test_step20: ${NFAIL} check(s) failed"
