*! test_step2.do -- elasticities and their standard errors against the R package
*! Reference: R_reference/out/R_elast_*.csv, R_semi_*.csv, R_slutsky.csv, R_compensated_q.csv

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do 01b_compat_lock_elasticities.do
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

log using "`ROOT'/replication/out/step2.log", replace text name(t2)

use "`ROOT'/examples/hixdata.dta", clear

qui easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,	///
	lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers) ///
	lnexpenditure(log_y)						///
	demographics(age hsex carown time tran)				///
	power(5) py pz zy compat nolog

local mats  elast_income elast_income_se elast_price elast_price_se	///
            elast_demo   elast_demo_se   semi_income semi_price		///
            slutsky      compensated_q
local files R_elast_income R_elast_income_se R_elast_price R_elast_price_se ///
            R_elast_demo   R_elast_demo_se   R_semi_income R_semi_price	    ///
            R_slutsky      R_compensated_q

tempname M
local nfail = 0
local nm : word count `mats'

di ""
di as txt "{hline 74}"
di as txt "  Elasticites : Stata vs R   (hixdata, spec vignette, mode compat)"
di as txt "{hline 74}"
di as txt "  matrice" _col(24) "dim" _col(34) "max |abs|" _col(50) "max |rel|" _col(66) "verdict"

forvalues m = 1/`nm' {
	local sm : word `m' of `mats'
	local rf : word `m' of `files'

	matrix `M' = e(`sm')
	local nr = rowsof(`M')
	local nc = colsof(`M')

	preserve
		import delimited using "`ROOT'/replication/R_reference/out/`rf'.csv", clear	///
			varnames(nonames) rowrange(2) stringcols(_all)
		qui count
		local nrcsv = r(N)
		local nccsv = c(k) - 1

		if (`nrcsv' != `nr') | (`nccsv' != `nc') {
			di as error "  `sm'" _col(24) "Stata `nr'x`nc' vs R `nrcsv'x`nccsv'" ///
			   _col(66) "DIM"
			local ++nfail
		}
		else {
			forvalues j = 1/`nc' {
				local jj = `j' + 1
				qui gen double x`j' = real(v`jj')
			}
			local amax = 0
			local rmax = 0
			forvalues i = 1/`nr' {
				forvalues j = 1/`nc' {
					local rv = x`j'[`i']
					local sv = `M'[`i',`j']
					local ad = abs(`sv' - `rv')
					local rd = `ad' / max(abs(`rv'), 1e-10)
					if `ad' > `amax' local amax = `ad'
					if `rd' > `rmax' local rmax = `rd'
				}
			}
			local tag = cond(`amax' < 1e-7, "ok", "FAIL")
			if "`tag'" == "FAIL" local ++nfail
			di as txt "  `sm'" _col(24) as txt "`nr'x`nc'"		///
			   _col(34) as res %12.4e `amax'			///
			   _col(50) as res %12.4e `rmax'			///
			   _col(66) as txt "`tag'"
		}
	restore
}

di as txt "{hline 74}"
if `nfail' == 0 di as result "  PASS -- les `nm' matrices concordent"
if `nfail' >  0 di as error  "  FAIL -- `nfail' matrice(s) en desaccord"

log close t2
