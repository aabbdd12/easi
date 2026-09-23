*! audit_homogeneity.do -- degree-zero homogeneity
*!
*! Multiplying every price AND total expenditure by the same factor cannot
*! change any elasticity: it is the same economy in different units.  The EASI
*! model satisfies this exactly, because adding up forces sum_j a_jk = 0 and
*! sum_j b_jk = 0, so both y and the budget shares are invariant.
*!
*! The elasticities reported by the R package are NOT invariant.  That is a
*! test no correct formula can fail, so it settles the question.

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do SectionA5_homogeneity.do
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

log using "`ROOT'/replication/out/homogeneity.log", replace text name(ho)

use "`ROOT'/examples/hixdata.dta", clear
local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local DE age hsex carown time tran

* shifted units: every log price and log expenditure moved by the same constant
local c = 1
local PR2
foreach v of local PR {
	qui gen double s_`v' = `v' + `c'
	local PR2 `PR2' s_`v'
}
qui gen double s_log_y = log_y + `c'

di ""
di as txt "{hline 76}"
di as txt "  HOMOGENEITE DE DEGRE ZERO   (prix et depense multiplies par e^`c')"
di as txt "{hline 76}"

*----------------------------------------------------------- shares invariant?
qui easi `SH', lnprices(`PR') lnexpenditure(log_y) demographics(`DE')	///
	power(5) py pz zy compat nolog
matrix EI0 = e(elast_income)
matrix EP0 = e(elast_price)
predict double y_a, y
predict double w_a*, shares

qui easi `SH', lnprices(`PR2') lnexpenditure(s_log_y) demographics(`DE')	///
	power(5) py pz zy compat nolog
matrix EI1 = e(elast_income)
matrix EP1 = e(elast_price)
predict double y_b, y
predict double w_b*, shares

qui gen double d_y = abs(y_a - y_b)
qui su d_y
* tolerance 1e-6: the shift multiplies the data's own |sum(shares)-1| slack
* (about 5e-8 in hixdata, single precision source) by the shift constant
di as txt "  1. l'utilite implicite y est invariante" _col(52) as res %11.3e r(max) ///
   _col(66) as txt cond(r(max) < 1e-6, "ok", "FAIL")

qui gen double d_w = 0
forvalues j = 1/9 {
	qui replace d_w = max(d_w, abs(w_a`j' - w_b`j'))
}
qui su d_w
di as txt "  2. les parts ajustees sont invariantes" _col(52) as res %11.3e r(max) ///
   _col(66) as txt cond(r(max) < 1e-6, "ok", "FAIL")

di as txt "     -> le MODELE est bien homogene de degre zero"
di ""

*--------------------------- same thing again, in corrected mode
qui easi `SH', lnprices(`PR') lnexpenditure(log_y) demographics(`DE')	///
	power(5) py pz zy nolog
matrix KI0 = e(elast_income)
matrix KP0 = e(elast_price)
qui easi `SH', lnprices(`PR2') lnexpenditure(s_log_y) demographics(`DE')	///
	power(5) py pz zy nolog
matrix KI1 = e(elast_income)
matrix KP1 = e(elast_price)

*------------------------------------------------- elasticities invariant?
mata: st_numscalar("dI", max(abs(st_matrix("EI0") - st_matrix("EI1"))))
mata: st_numscalar("dP", max(abs(st_matrix("EP0") - st_matrix("EP1"))))
di as txt "  3. elasticites-depense rapportees, ecart max" _col(52) as res %11.3e dI ///
   _col(66) as txt cond(dI < 1e-6, "ok", "FAIL")
di as txt "  4. elasticites-prix rapportees,    ecart max" _col(52) as res %11.3e dP ///
   _col(66) as txt cond(dP < 1e-6, "ok", "FAIL")

mata: st_numscalar("kI", max(abs(st_matrix("KI0") - st_matrix("KI1"))))
mata: st_numscalar("kP", max(abs(st_matrix("KP0") - st_matrix("KP1"))))
di ""
di as txt "  MODE CORRIGE :"
di as txt "  5. elasticites-depense, ecart max" _col(52) as res %11.3e kI	///
   _col(66) as txt cond(kI < 1e-6, "ok", "FAIL")
di as txt "  6. elasticites-prix,    ecart max" _col(52) as res %11.3e kP	///
   _col(66) as txt cond(kP < 1e-6, "ok", "FAIL")

di ""
di as txt "  Elasticites-depense dans les deux systemes d'unites :"
matrix CMP = EI0 \ EI1 \ (EI1 - EI0)
matrix rownames CMP = unites_1 unites_2 ecart
matlist CMP, format(%9.4f) twidth(11)

di as txt "{hline 76}"
di as txt "  Le modele est invariant (1 et 2), les elasticites rapportees ne le"
di as txt "  sont pas (3 et 4) : la formule est donc fausse, independamment de"
di as txt "  toute reference externe."
di as txt "{hline 76}"

log close ho
