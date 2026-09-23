*! test_step9.do -- do the reported standard errors match the bootstrap?
*!                  Simple SRS case, no weights, no survey design.
*!
*! test_step8 answered this for the coefficients: the homoskedastic 3SLS
*! standard error is 23% too small on the terms in y, and -vce(robust)- -- now
*! the default -- closes the gap.  The generated-regressor term that Pendakur's
*! comment points at is real but worth under half a percent (test_step11).
*!
*! Here we ask it of what users actually report: the ELASTICITIES.  Their
*! reported standard error is the R package's heuristic -- the median across
*! households of the pointwise delta-method standard error of the corresponding
*! SEMI-elasticity, divided by the mean budget share.  For the price elasticity
*! that heuristic ignores the B(C+D+G) term entirely, so it is not obvious that
*! it should agree with anything.
*!
*! The bootstrap re-runs the whole procedure, elasticity computation included,
*! so it is the reference.

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table5.do
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

log using "`ROOT'/replication/out/step9.log", replace text name(t9)

use "`ROOT'/examples/hixdata.dta", clear
local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local SPEC lnprices(`PR') lnexpenditure(log_y)				///
	   demographics(age hsex carown) power(3) nolog
local J 9

di ""
di as txt "{hline 78}"
di as txt "  Elasticites : ecart-type rapporte contre bootstrap  (SRS, sans poids)"
di as txt "{hline 78}"

*------------------------------------------------ reported
qui easi `SH', `SPEC'
matrix EI  = e(elast_income)
matrix EIS = e(elast_income_se)
matrix EP  = e(elast_price)
matrix EPS = e(elast_price_se)

*------------------------------------------------ bootstrap
* bootstrap cannot subscript an e() matrix inside its expression list, so the
* quantities of interest are returned by a small rclass wrapper -- the standard
* route.  easi exits r(430) if it fails to converge, so the wrapper propagates
* the failure and bootstrap counts it.
global SH_  "`SH'"
global SPEC_ "`SPEC'"

capture program drop _elboot
program define _elboot, rclass
	version 14.2
	tempname EI EP
	qui easi $SH_, $SPEC_ noelastse
	matrix `EI' = e(elast_income)
	matrix `EP' = e(elast_price)
	forvalues j = 1/9 {
		return scalar ei`j' = `EI'[1,`j']
		return scalar ep`j' = `EP'[`j',`j']
	}
end

local el
forvalues j = 1/`J' {
	local el `el' (ei`j': r(ei`j'))
}
forvalues j = 1/`J' {
	local el `el' (ep`j': r(ep`j'))
}

local REPS 400
set seed 20260921
qui bootstrap `el', reps(`REPS') nodots nowarn: _elboot
matrix BB = e(b)
matrix BV = e(V)
local nfail = e(N_misreps)
if "`nfail'" == "" local nfail 0
di as txt "  replications echouees : " as res `nfail' as txt " sur `REPS'"

*------------------------------------------------ comparison
preserve
	clear
	qui set obs `=2*`J''
	gen str12 nm   = ""
	gen str10 kind = ""
	gen double est = .
	gen double ser = .
	gen double seb = .
	forvalues j = 1/`J' {
		local k = `J' + `j'
		qui replace nm   = "depense `j'"     in `j'
		qui replace kind = "depense"         in `j'
		qui replace est  = EI[1,`j']         in `j'
		qui replace ser  = abs(EIS[1,`j'])   in `j'
		qui replace seb  = sqrt(BV[`j',`j']) in `j'

		qui replace nm   = "prix `j'"        in `k'
		qui replace kind = "prix"            in `k'
		qui replace est  = EP[`j',`j']       in `k'
		qui replace ser  = abs(EPS[`j',`j']) in `k'
		qui replace seb  = sqrt(BV[`k',`k']) in `k'
	}
	gen double ratio = ser / seb

	di ""
	di as txt "  Rapport  se(delta, methode corrigee) / se(bootstrap)"
	di as txt "    (1 = concordance ; < 1 = l'ecart-type rapporte SOUS-estime)"
	local worst = 0
	foreach t in depense prix {
		qui su ratio if kind == "`t'", detail
		di as txt "    elasticites-`t'" _col(26) "mediane " as res %6.3f r(p50) ///
		   as txt "   min " as res %6.3f r(min)				///
		   as txt "   max " as res %6.3f r(max)
		local worst = max(`worst', abs(r(p50) - 1))
	}
	di ""
	di as txt "    ecart median maximal au bootstrap " as res %6.3f `worst'	///
	   _col(52) as res cond(`worst' < 0.10, "ok", "ECHEC")

	di ""
	di as txt "  Detail :"
	format est ser seb ratio %9.4f
	list nm est ser seb ratio, noobs sepby(kind) abbreviate(12)

	qui save "`ROOT'/replication/out/elast_se_compare.dta", replace
restore

di as txt "{hline 78}"

log close t9
