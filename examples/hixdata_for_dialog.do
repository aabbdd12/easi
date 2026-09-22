*! hixdata_for_dialog.do -- prepare hixdata for the dialog box
*!
*! The dialog offers Prices and Expenditure in LEVELS: it emits prices() and
*! expenditure(), and easi takes logarithms internally.  hixdata is already in
*! logarithms -- pfoodh..ppers are log prices and log_y is log expenditure,
*! all of them centred and therefore often NEGATIVE -- so feeding them to the
*! dialog as they stand asks easi for the log of a negative number.
*!
*! This exponentiates them once so the dialog has something to work with, and
*! checks that the result is identical to the log-based command.

clear all
set more off

local ROOT "C:/Users/aabd/OneDrive/Desktop/EASI_project"
adopath ++ "`ROOT'/src"

use "`ROOT'/examples/hixdata.dta", clear

local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local LP pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers

di ""
di as txt "{hline 70}"
di as txt "  hixdata is in logarithms: how many negative values?"
di as txt "{hline 70}"
foreach v of local LP {
	qui count if `v' < 0
	di as txt "   `v'" _col(14) as res %6.0f r(N) as txt " negative of " ///
	   as res _N
}
qui count if log_y < 0
di as txt "   log_y" _col(14) as res %6.0f r(N) as txt " negative of " as res _N
di as txt "   -> prices() and expenditure() would give ln() of a negative number"

*------------------------------------------------- level variables
local i 0
local PLEV
foreach v of local LP {
	local ++i
	gen double p`i' = exp(`v')
	label var p`i' "price of good `i', level"
	local PLEV `PLEV' p`i'
}
gen double totexp = exp(log_y)
label var totexp "total expenditure, level"

di ""
di as txt "  created: " as res "`PLEV' totexp"

*------------------------------------------------- the two paths agree?
qui easi `SH', lnprices(`LP') lnexpenditure(log_y)			///
	demographics(age hsex carown) power(3) nolog noelastse
matrix A = e(elast_exp)
qui easi `SH', prices(`PLEV') expenditure(totexp)			///
	demographics(age hsex carown) power(3) nolog noelastse
matrix B = e(elast_exp)
mata: st_numscalar("d", max(abs(st_matrix("A") - st_matrix("B"))))
di ""
di as txt "  lnprices() path vs prices() path, max difference" _col(56)	///
   as res %11.3e d _col(70) as res cond(d < 1e-10, "ok", "ECHEC")

qui save "`ROOT'/examples/hixdata_levels.dta", replace
di ""
di as txt "  saved: examples/hixdata_levels.dta"
di as txt "{hline 70}"
