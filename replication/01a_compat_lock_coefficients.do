*! test_step1.do -- does the Mata estimator reproduce the R package?
*! Reference: R_reference/out/R_coef.csv (hixdata, vignette specification)

clear all
set more off

local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"

log using "`ROOT'/replication/out/step1.log", replace text name(t1)

use "`ROOT'/examples/hixdata.dta", clear

timer clear 1
timer on 1
easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,	///
	lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers) ///
	lnexpenditure(log_y)						///
	demographics(age hsex carown time tran)				///
	power(5) py pz zy compat
timer off 1
timer list 1

matrix b = e(b)
local K = colsof(b)
di as txt "K = `K'   (R reference: 320)"

*---------------------------------------------------------------- compare
matrix V = e(V)

preserve
	import delimited using "`ROOT'/replication/R_reference/out/R_coef.csv", clear varnames(1) ///
		stringcols(_all)
	rename v1 rname
	gen double rb  = real(estimate)
	gen double rse = real(stderror)
	keep rname rb rse
	tempfile rcoef
	qui save `rcoef'
	qui count
	di as txt "R coefficients read: " r(N)
restore

preserve
	clear
	qui set obs `K'
	gen double sb  = .
	gen double sse = .
	forvalues i = 1/`K' {
		qui replace sb  = b[1,`i']          in `i'
		qui replace sse = sqrt(V[`i',`i'])  in `i'
	}
	qui merge 1:1 _n using `rcoef', nogen

	foreach q in b se {
		gen double ad_`q' = abs(s`q' - r`q')
		gen double rd_`q' = ad_`q' / max(abs(r`q'), 1e-10)
		qui su ad_`q'
		local amax_`q' = r(max)
		qui su rd_`q'
		local rmax_`q' = r(max)
	}
	qui count
	local nR = r(N)

	di ""
	di as txt "{hline 66}"
	di as txt "  Stata vs R -- `nR' coefficients (hixdata, vignette spec)"
	di as txt "{hline 66}"
	di as txt "               " _col(26) "max |abs|" _col(44) "max |rel|"
	di as txt "  coefficients " _col(26) as res %12.4e `amax_b'  _col(44) %12.4e `rmax_b'
	di as txt "  std. errors  " _col(26) as res %12.4e `amax_se' _col(44) %12.4e `rmax_se'
	di as txt "{hline 66}"

	local ok = (`amax_b' < 1e-8) & (`amax_se' < 1e-8)
	if `ok' di as result "  PASS -- verrou 1e-8 atteint sur b et se"
	if !`ok' {
		di as error "  FAIL"
		gen double ad = max(ad_b, ad_se)
		gsort -ad
		list rname sb rb ad_b sse rse ad_se in 1/10, noobs
	}
restore

log close t1
