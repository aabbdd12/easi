*! test_step10.do -- alternate which good is dropped
*!
*! The system drops one good and recovers it by adding up.  Estimating a second
*! time with a different good dropped gives two things:
*!
*!  (a) a REFERENCE for the dropped good.  In run A good 9 is recovered by
*!      adding up; in run B it is an equation of its own, with its own
*!      coefficients and its own standard error.  Comparing the two checks the
*!      adding-up recovery -- the quantity that was 14 times too large before
*!      the delta-method fix -- without going through the bootstrap.
*!
*!  (b) a test of INVARIANCE to the deletion choice.  For a singular demand
*!      system the estimator is invariant to the equation deleted only if it is
*!      the ML/iterated one (Barten 1969).  systemfit runs with maxiter = 1, so
*!      Sigma is formed once from the 2SLS residuals and never re-iterated; the
*!      loop we do iterate is the one on y.  Invariance is therefore NOT
*!      guaranteed, and how far it fails is worth knowing.

clear all
set more off

local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"

log using "`ROOT'/replication/out/step10.log", replace text name(ta)

use "`ROOT'/examples/hixdata.dta", clear

* order A: spers is last, so it is the good dropped
local SHA sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PRA pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
* order B: sfoodh and spers swapped, so sfoodh is dropped instead
local SHB spers  sfoodr srent soper sfurn scloth stranop srecr sfoodh
local PRB ppers  pfoodr prent poper pfurn pcloth ptranop precr pfoodh

local OPT lnexpenditure(log_y) demographics(age hsex carown) power(3) nolog

di ""
di as txt "{hline 76}"
di as txt "  Bien omis alterne : reference externe et test d'invariance"
di as txt "{hline 76}"

qui easi `SHA', lnprices(`PRA') `OPT'
matrix IA  = e(elast_income)
matrix ISA = e(elast_income_se)
matrix PA  = e(elast_price)
matrix PSA = e(elast_price_se)
local iterA = e(iter)

qui easi `SHB', lnprices(`PRB') `OPT'
matrix IB  = e(elast_income)
matrix ISB = e(elast_income_se)
matrix PB  = e(elast_price)
matrix PSB = e(elast_price_se)
local iterB = e(iter)

di as txt "  run A : spers  omis (recupere par adding-up), " as res `iterA' as txt " iterations"
di as txt "  run B : sfoodh omis, " as res `iterB' as txt " iterations"

* B is A with positions 1 and 9 exchanged
mata:
    pm  = (9, 2, 3, 4, 5, 6, 7, 8, 1)
    IBp = st_matrix("IB")[., pm]
    ISBp= st_matrix("ISB")[., pm]
    PBp = st_matrix("PB")[pm, pm]
    PSBp= st_matrix("PSB")[pm, pm]
    st_matrix("IBp", IBp);  st_matrix("ISBp", ISBp)
    st_matrix("PBp", PBp);  st_matrix("PSBp", PSBp)
    st_numscalar("dI", max(abs(st_matrix("IA") - IBp)))
    st_numscalar("dP", max(abs(st_matrix("PA") - PBp)))
end

di ""
di as txt "  (b) INVARIANCE au bien omis (apres permutation)"
di as txt "      elasticites-depense, ecart max" _col(48) as res %11.3e dI
di as txt "      elasticites-prix,    ecart max" _col(48) as res %11.3e dP
* 1e-6: what is left is the tolerance of the Sigma loop (1e-11 on b)
*       propagating; not a failure of invariance
di as txt "      -> l'estimateur " as res cond(dI < 1e-6 & dP < 1e-6,	///
   "EST invariant", "n'est PAS invariant") as txt " au choix du bien omis"

di ""
di as txt "  (a) Le bien omis, vu des deux cotes"
di as txt "      spers  : recupere par adding-up (A) contre estime (B)"
di as txt "      sfoodh : estime (A) contre recupere par adding-up (B)"
di ""
di as txt "      bien          elast.A    elast.B      se A       se B    se B/A"
foreach pair in "spers 9 1" "sfoodh 1 9" {
	gettoken g rest : pair
	gettoken ja jb  : rest
	local eA = IA[1,`ja']
	local eB = IB[1,`jb']
	local sA = ISA[1,`ja']
	local sB = ISB[1,`jb']
	di as txt "      `g'" _col(20) as res %9.5f `eA' _col(31) %9.5f `eB'	///
	   _col(42) %9.5f `sA' _col(53) %9.5f `sB' _col(64) %7.3f `sB'/`sA'
}

di ""
di as txt "  Elasticites-depense, les deux ordres :"
matrix CMP = IA \ IBp \ (IBp - IA)
matrix rownames CMP = ordre_A ordre_B ecart
matlist CMP, format(%9.5f) twidth(10)

di ""
di as txt "  Ecarts-types, les deux ordres :"
matrix CMS = ISA \ ISBp \ (ISBp - ISA)
matrix rownames CMS = ordre_A ordre_B ecart
matlist CMS, format(%9.5f) twidth(10)

di as txt "{hline 76}"

log close ta
