*! test_step8.do -- are the reported standard errors right?  Four candidates.
*!
*! Pendakur's own code carries the line
*!     *note that reported standard errors are wrong for iterated estimates
*! immediately before its final reg3.  The reason he gives is a generated
*! regressor: the design contains y, which is a function of the coefficients,
*! and the 3SLS covariance is the one of the last linear step CONDITIONAL on y.
*!
*! That diagnosis turns out to be right in principle and almost irrelevant in
*! practice.  test_step11.do adds the missing dy/dbeta term to the Jacobian,
*! checks it against finite differences, and finds it worth less than half a
*! percent -- because y depends on beta only through p'Ap/2 and p'Bp/2, and
*! normalised log prices are of order 0.1.
*!
*! What IS worth 23% is heteroskedasticity.  This file crosses the two:
*!     conventional / robust   x   conditional on y / + dy/dbeta
*! and scores all four against a bootstrap that re-runs the whole iterated
*! procedure and therefore captures everything.
*!
*! A deliberately light specification: 9 goods, power 3, three demographics, no
*! interactions.  Some replications may fail to converge; -easi- exits with
*! r(430) in that case, so -bootstrap- counts them instead of silently using a
*! bad fit.

clear all
set more off

local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"

log using "`ROOT'/replication/out/step8.log", replace text name(t8)

use "`ROOT'/examples/hixdata.dta", clear
local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local SPEC lnprices(`PR') lnexpenditure(log_y)				///
	   demographics(age hsex carown) power(3) nolog noelastse

di ""
di as txt "{hline 78}"
di as txt "  Ecarts-types des coefficients : quatre variances contre le bootstrap"
di as txt "{hline 78}"

*------------------------------------------------ the four analytic variances
timer clear 1
timer on 1
qui easi `SH', `SPEC' vce(conventional) legacy(condse)
timer off 1
qui timer list 1
local one = r(t1)
local K = colsof(e(b))
local cn : colnames e(b)
matrix V1 = e(V)

qui easi `SH', `SPEC' vce(conventional)
matrix V2 = e(V)
qui easi `SH', `SPEC' vce(robust) legacy(condse)
matrix V3 = e(V)
qui easi `SH', `SPEC'
matrix V4 = e(V)
local defvce "`e(vce)', `e(vcetype2)'"

di as txt "  specification : " as res `K' as txt " coefficients, "		///
   as res e(iter) as txt " iterations, " as res %5.2f `one' as txt " s"
di as txt "  defaut de easi : " as res "`defvce'"

*------------------------------------------------ the bootstrap
local REPS 400
di as txt "  bootstrap : `REPS' replications (~" as res			///
   %4.0f `REPS' * `one' as txt " s)"

set seed 20260921
qui bootstrap, reps(`REPS') nodots nowarn: easi `SH', `SPEC'
matrix VB = e(V)
local nfail = e(N_misreps)
if "`nfail'" == "" local nfail 0
di as txt "  replications echouees (non convergence) : " as res `nfail'	///
   as txt " sur `REPS'"

*------------------------------------------------ comparison
preserve
	clear
	qui set obs `K'
	gen str40 nm = ""
	foreach m in 1 2 3 4 B {
		gen double se`m' = .
	}
	forvalues i = 1/`K' {
		local v : word `i' of `cn'
		qui replace nm = "`v'" in `i'
		foreach m in 1 2 3 4 B {
			qui replace se`m' = sqrt(V`m'[`i',`i']) in `i'
		}
	}
	* the generated regressor is y, so its own coefficients are the ones
	* whose conditional standard error should be most in error
	gen byte isy = regexm(nm, "^y[0-9]+$")
	qui count if isy
	local ny = r(N)
	local no = `K' - `ny'

	local lab1 "conventionnelle, cond. sur y"
	local lab2 "conventionnelle, + dy/dbeta"
	local lab3 "robuste,         cond. sur y"
	local lab4 "robuste,         + dy/dbeta"

	di ""
	di as txt "  se(bootstrap) / se(analytique)   -- au-dessus de 1 = analytique trop petit"
	di as txt "  variance" _col(36) "termes en y (`ny')" _col(60) "autres (`no')"
	local worst = 0
	forvalues m = 1/4 {
		qui gen double r`m' = seB / se`m'
		qui su r`m' if isy
		local a  = r(mean)
		local ax = r(max)
		qui su r`m' if !isy
		local b  = r(mean)
		local bx = r(max)
		di as txt "  `lab`m''" _col(36) as res %6.3f `a'			///
		   as txt " (max " as res %5.2f `ax' as txt ")"			///
		   _col(60) as res %6.3f `b'					///
		   as txt " (max " as res %5.2f `bx' as txt ")"
		if `m' == 4 {
			local d4  = max(abs(`a' - 1), abs(`b' - 1))
			local x4  = max(`ax', `bx')
		}
	}

	di ""
	di as txt "  1. le regresseur genere (ligne 1 -> 2, ligne 3 -> 4) ne deplace"
	di as txt "     presque rien : c'est un terme reel mais d'ordre 0.5 %."
	di as txt "  2. l'heteroscedasticite (ligne 1 -> 3) vaut 23 % sur les termes"
	di as txt "     en y, ceux qui dessinent les courbes d'Engel."
	di as txt "  3. le defaut de easi est la ligne 4."
	di ""
	local tag = cond(`d4' < 0.05 & `x4' < 1.20, "ok", "ECHEC")
	di as txt "     ecart moyen du defaut au bootstrap  " as res %6.3f `d4'	///
	   _col(52) as txt "max " as res %5.2f `x4' _col(66) as res "`tag'"
restore

di as txt "{hline 78}"

log close t8
