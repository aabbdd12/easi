*! test_step14.do -- les ecarts-types sous pweight, contre le bootstrap
*!
*! Regle du projet : valider les ecarts-types sur SRS, PUIS sous pweight, PUIS
*! sous plan de sondage -- jamais dans l'autre sens.  Le SRS est verrouille par
*! test_step8 et test_step9 sur hixdata.  Ici on monte d'un cran.
*!
*! Le bootstrap doit rebattre l'UNITE DE TIRAGE du plan.  Sous pweight, c'est
*! l'observation AVEC son poids : le poids est une colonne, il voyage avec la
*! ligne, donc -bootstrap: easi ... [pw=w]- fait exactement cela.
*!
*! Donnees : le banc mexicain, qui porte un vrai poids d'expansion (et, pour
*! l'etape suivante, des strates et des PSU).  hixdata reste reserve aux
*! estimateurs sans poids.
*!
*! Ce que ce test NE valide PAS : la stratification et les grappes.  Sous
*! -vce(robust)-, chaque menage est sa propre unite ; le plan reel a 703 PSU
*! dans 28 strates.  C'est l'objet de vce(svy), l'etape suivante.

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table4.do
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

log using "`ROOT'/replication/out/step14.log", replace text name(ta)

use "`ROOT'/examples/mex_bench.dta", clear

local SH w1 w2 w3
local PR lp1 lp2 lp3
local SPEC lnprices(`PR') lnexpenditure(lx) demographics(z1 isMale)	///
	   power(3) nolog noelastse
local REPS 400

di ""
di as txt "{hline 78}"
di as txt "  Ecarts-types sous pweight : analytique contre bootstrap"
di as txt "{hline 78}"

* -bootstrap- refuses a weight on the prefixed command, so the weight goes
* INSIDE a wrapper.  That is not a trick: it is the only way to make the
* bootstrap resample the row together with its weight, which is the sampling
* unit under pweight.
global SH_   "`SH'"
global SPEC_ "`SPEC'"
global WT_   ""
capture program drop _eb
program define _eb, eclass
	version 14.2
	if "$WT_" == "" qui easi $SH_, $SPEC_
	else            qui easi $SH_ $WT_, $SPEC_
end

qui easi `SH', `SPEC'
local K = colsof(e(b))
local cn : colnames e(b)
di as txt "  banc mexicain : " as res _N as txt " menages, " as res `K'	///
   as txt " coefficients, " as res e(iter) as txt " iterations"
di as txt "  " as res `REPS' as txt " replications par configuration"

*======================================================= coefficients
tempname R
matrix `R' = J(2, 4, .)
local rn
local i 0
foreach cfg in "SRS" "pweight" {
	local ++i
	local rn `rn' `cfg'
	if "`cfg'" == "SRS" {
		local W
		global WT_ ""
	}
	else {
		local W [pw=sweight]
		global WT_ "[pw=sweight]"
	}

	qui easi `SH' `W', `SPEC'
	matrix VA = e(V)

	set seed 20260922
	qui bootstrap, reps(`REPS') nodots nowarn: _eb
	matrix VB = e(V)
	local nfail = e(N_misreps)
	if "`nfail'" == "" local nfail 0

	preserve
		clear
		qui set obs `K'
		gen str32 nm = ""
		gen double sa = .
		gen double sb = .
		forvalues j = 1/`K' {
			local v : word `j' of `cn'
			qui replace nm = "`v'"              in `j'
			qui replace sa = sqrt(VA[`j',`j'])  in `j'
			qui replace sb = sqrt(VB[`j',`j'])  in `j'
		}
		gen double r = sb / sa
		gen byte isy = regexm(nm, "^y[0-9]+$")
		qui su r
		matrix `R'[`i',1] = r(mean)
		matrix `R'[`i',2] = r(max)
		qui su r if isy
		matrix `R'[`i',3] = r(mean)
		qui su r if !isy
		matrix `R'[`i',4] = r(mean)
	restore
	di as txt "  `cfg' : " as res `nfail' as txt " replication(s) echouee(s)"
}
matrix rownames `R' = `rn'
matrix colnames `R' = moyenne max termes_y autres

di ""
di as txt "  se(bootstrap) / se(analytique), coefficients"
di as txt "  au-dessus de 1 = l'analytique sous-estime"
matlist `R', format(%9.3f) twidth(12) border(all) rowtitle("")

*======================================================= elasticites
* bootstrap ne sait pas indexer une matrice e() dans sa liste d'expressions,
* d'ou le relais rclass -- meme procede que test_step9.
capture program drop _elb
program define _elb, rclass
	version 14.2
	tempname EI EP
	if "$WT_" == "" qui easi $SH_, $SPEC_
	else            qui easi $SH_ $WT_, $SPEC_
	matrix `EI' = e(elast_exp)
	matrix `EP' = e(elast_price_nc)
	forvalues j = 1/3 {
		return scalar ei`j' = `EI'[1,`j']
		return scalar ep`j' = `EP'[`j',`j']
	}
end

local el
forvalues j = 1/3 {
	local el `el' (ei`j': r(ei`j'))
}
forvalues j = 1/3 {
	local el `el' (ep`j': r(ep`j'))
}

global WT_ "[pw=sweight]"
qui easi `SH' [pw=sweight], lnprices(`PR') lnexpenditure(lx)		///
	demographics(z1 isMale) power(3) nolog
matrix EI = e(elast_exp)
matrix ES = e(elast_exp_se)
matrix EP = e(elast_price_nc)
matrix PS = e(elast_price_nc_se)

set seed 20260922
qui bootstrap `el', reps(`REPS') nodots nowarn: _elb
matrix BV = e(V)

di ""
di as txt "  se(delta) / se(bootstrap), elasticites sous pweight"
di as txt "  bien" _col(16) "estime" _col(30) "se delta" _col(44)	///
   "se boot" _col(58) "rapport"
local worst 0
forvalues j = 1/3 {
	local e1 = EI[1,`j']
	local s1 = abs(ES[1,`j'])
	local b1 = sqrt(BV[`j',`j'])
	local r1 = `s1' / `b1'
	local worst = max(`worst', abs(`r1' - 1))
	di as txt "  depense `j'" _col(14) as res %10.4f `e1' _col(28)	///
	   %10.4f `s1' _col(42) %10.4f `b1' _col(56) %10.4f `r1'
}
forvalues j = 1/3 {
	local k = 3 + `j'
	local e1 = EP[`j',`j']
	local s1 = abs(PS[`j',`j'])
	local b1 = sqrt(BV[`k',`k'])
	local r1 = `s1' / `b1'
	local worst = max(`worst', abs(`r1' - 1))
	di as txt "  prix `j'" _col(14) as res %10.4f `e1' _col(28)	///
	   %10.4f `s1' _col(42) %10.4f `b1' _col(56) %10.4f `r1'
}

di ""
local cmax = max(abs(`R'[2,3] - 1), abs(`R'[2,4] - 1))
di as txt "  ecart max au bootstrap, coefficients pweight" _col(56)	///
   as res %9.3f `cmax' _col(68)						///
   as res cond(`cmax' < 0.12, "ok", "ECHEC")
di as txt "  ecart max au bootstrap, elasticites pweight" _col(56)	///
   as res %9.3f `worst' _col(68)					///
   as res cond(`worst' < 0.15, "ok", "ECHEC")

di ""
di as txt "  Rappel : ce test ne couvre ni la stratification ni les grappes."
di as txt "  Le banc a 703 PSU dans 28 strates ; vce(robust) traite chaque"
di as txt "  menage comme sa propre unite.  C'est l'objet de vce(svy)."
di as txt "{hline 78}"

log close ta
