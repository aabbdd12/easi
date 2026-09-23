*! test_step3.do -- how much does each correction move the results?
*! Starts from -compat- (= the R package, bit for bit) and switches the fixes
*! on one at a time, on the reference example.

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do 02_effect_of_each_correction.do
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

log using "`ROOT'/replication/out/step3.log", replace text name(t3)

use "`ROOT'/examples/hixdata.dta", clear

local SPEC sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,	///
	lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers) ///
	lnexpenditure(log_y)						///
	demographics(age hsex carown time tran)				///
	power(5) py pz zy nolog

local OUT elast_income elast_price elast_demo elast_income_se elast_price_se elast_demo_se

*---------------------------------------------------------------- reference
qui easi `SPEC' legacy(all)
local Kref = colsof(e(b))
foreach m of local OUT {
	tempname R_`m'
	matrix `R_`m'' = e(`m')
}
matrix REF_ei = e(elast_income)

di ""
di as txt "{hline 92}"
di as txt "  Effet de chaque correction -- ecart max par rapport a -compat- (= paquet R)"
di as txt "  hixdata, 9 biens, power 5, py+pz+zy, 5 demographiques"
di as txt "{hline 92}"
di as txt "  correction activee" _col(26) "K"				///
   _col(32) "el.depense" _col(46) "el.prix" _col(60) "el.demo"		///
   _col(72) "se.dep" _col(82) "se.prix"
di as txt "{hline 92}"

local ALL interpz instr quant ez lastse

*------------------------------------------------------- one fix at a time
foreach f of local ALL {
	local keep : list ALL - f
	if "`keep'" == "" local keep none
	qui easi `SPEC' legacy(`keep')
	local K = colsof(e(b))

	local line
	foreach m of local OUT {
		tempname C
		matrix `C' = e(`m')
		mata: st_local("d", strofreal(max(abs(st_matrix("`C'") -		///
			st_matrix("`R_`m''")))))
		local d_`m' = `d'
	}
	di as txt "  +`f'" _col(26) as res `K'					///
	   _col(32) %11.3e `d_elast_income'					///
	   _col(46) %11.3e `d_elast_price'					///
	   _col(60) %11.3e `d_elast_demo'					///
	   _col(72) %9.2e  `d_elast_income_se'					///
	   _col(82) %9.2e  `d_elast_price_se'
}

*----------------------------------------------------------- all the fixes
qui easi `SPEC' legacy(none)
local K = colsof(e(b))
foreach m of local OUT {
	tempname C
	matrix `C' = e(`m')
	mata: st_local("d", strofreal(max(abs(st_matrix("`C'") -		///
		st_matrix("`R_`m''")))))
	local d_`m' = `d'
}
di as txt "{hline 92}"
di as txt "  TOUTES (mode corrige)" _col(26) as res `K'			///
   _col(32) %11.3e `d_elast_income'					///
   _col(46) %11.3e `d_elast_price'					///
   _col(60) %11.3e `d_elast_demo'					///
   _col(72) %9.2e  `d_elast_income_se'					///
   _col(82) %9.2e  `d_elast_price_se'
di as txt "{hline 92}"
di as txt "  (K = nombre de coefficients ; reference K = `Kref')"

*------------------------------------------- elasticites-depense cote a cote
matrix COR_ei = e(elast_income)
di ""
di as txt "  Elasticites-depense : compat (= R) vs mode corrige"
matrix CMP = REF_ei \ COR_ei \ (COR_ei - REF_ei)
matrix rownames CMP = compat corrige ecart
matlist CMP, format(%9.4f) twidth(12)

*--------------------------------------------- ecart-type du dernier article
qui easi `SPEC' legacy(all)
scalar se_last_compat = e(elast_income_se)[1, 9]
qui easi `SPEC' legacy(none)
scalar se_last_corr = e(elast_income_se)[1, 9]
di ""
di as txt "  Ecart-type de l'elasticite-depense du dernier bien (spers) :"
di as txt "      compat  = " as res %9.4f se_last_compat as txt "   <- signe negatif (bug 5)"
di as txt "      corrige = " as res %9.4f se_last_corr

log close t3
