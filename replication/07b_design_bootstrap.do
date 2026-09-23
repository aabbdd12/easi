*! test_step15.do -- ce que le plan de sondage coute, mesure contre le bon oracle
*!
*! test_step14 comparait vce(robust) a un bootstrap qui rebat des OBSERVATIONS.
*! Les deux ignorent la meme chose -- la stratification et les grappes -- donc
*! leur accord ne prouvait rien sur le plan : c'etait un miroir, pas un oracle.
*!
*! Le bon oracle rebat l'unite de tirage du plan : les PSU A L'INTERIEUR de
*! chaque strate, avec remise, m_h = n_h, en emportant le poids de chaque ligne.
*! -idcluster()- est obligatoire, faute de quoi une PSU tiree deux fois compte
*! pour une seule grappe et la variance sort trop petite.
*!
*! Le banc mexicain : 2 477 menages, 703 PSU, 28 strates, poids d'expansion.
*!
*! Ce que ce test mesure : l'ecart entre ce que easi sait faire aujourd'hui --
*! vce(robust), vce(cluster) -- et la verite du plan.  Cet ecart est le cahier
*! des charges de vce(svy).
*!
*! Reserve : le bootstrap a grappes de Stata ne re-echelonne pas les poids a
*! chaque replication (ce que ferait un bootstrap de Rao-Wu).  Pour un oracle
*! de validation c'est l'usage courant et suffisant.

clear all
set more off

local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"

log using "`ROOT'/replication/out/step15.log", replace text name(ta)

use "`ROOT'/examples/mex_bench.dta", clear

local SH w1 w2 w3
local PR lp1 lp2 lp3
local SPEC lnprices(`PR') lnexpenditure(lx) demographics(z1 isMale)	///
	   power(3) nolog noelastse
local REPS 400

di ""
di as txt "{hline 78}"
di as txt "  Le plan de sondage : ce que vce(robust) et vce(cluster) manquent"
di as txt "{hline 78}"

qui su st
local nst = r(max)
qui su pu
local npu = r(max)
di as txt "  banc : " as res _N as txt " menages, " as res `npu'		///
   as txt " PSU, " as res `nst' as txt " strates, poids d'expansion"

* Le poids ne peut pas etre passe a la commande prefixee par -bootstrap-,
* il doit vivre dans l'enveloppe.  C'est aussi ce qui garantit qu'il est
* rebattu avec sa ligne.
global SH_   "`SH'"
global SPEC_ "`SPEC'"
capture program drop _eb
program define _eb, eclass
	version 14.2
	qui easi $SH_ [pw=sweight], $SPEC_
end

*--------------------------------------------------- les variances analytiques
qui easi `SH' [pw=sweight], `SPEC' vce(robust)
matrix VR = e(V)
local K = colsof(e(b))
local cn : colnames e(b)
matrix BR = e(b)

qui easi `SH' [pw=sweight], `SPEC' vce(cluster pu)
matrix VC = e(V)
matrix BC = e(b)

qui easi `SH', `SPEC' vce(svy)
matrix VS = e(V)
local n_str = e(N_strata)
local n_psu2 = e(N_psu)

mata: st_numscalar("dpt", max(abs(st_matrix("BR") - st_matrix("BC"))))
di as txt "  les estimes ponctuels ne dependent pas du vce" _col(60)	///
   as res %11.3e dpt _col(72) as res cond(dpt < 1e-10, "ok", "ECHEC")

*--------------------------------------------------- l'oracle : bootstrap de plan
di as txt "  bootstrap de plan : " as res `REPS' as txt			///
   " replications, PSU tirees dans les strates"
set seed 20260922
qui bootstrap, reps(`REPS') strata(st) cluster(pu) idcluster(_newpu)	///
	nodots nowarn: _eb
matrix VB = e(V)
local nfail = e(N_misreps)
if "`nfail'" == "" local nfail 0
di as txt "  replications echouees" _col(60) as res %11.0f `nfail'

*--------------------------------------------------- comparaison
preserve
	clear
	qui set obs `K'
	gen str32 nm = ""
	gen double sr = .
	gen double sc = .
	gen double ss = .
	gen double sb = .
	forvalues j = 1/`K' {
		local v : word `j' of `cn'
		qui replace nm = "`v'"             in `j'
		qui replace sr = sqrt(VR[`j',`j']) in `j'
		qui replace sc = sqrt(VC[`j',`j']) in `j'
		qui replace ss = sqrt(VS[`j',`j']) in `j'
		qui replace sb = sqrt(VB[`j',`j']) in `j'
	}
	gen double r_rob = sb / sr
	gen double r_clu = sb / sc
	gen double r_svy = sb / ss

	di ""
	di as txt "  se(bootstrap de plan) / se(analytique)"
	di as txt "  au-dessus de 1 = l'analytique sous-estime"
	di as txt "  variante" _col(30) "moyenne" _col(44) "min"	_col(58) "max"
	foreach v in rob clu svy {
		local lbl "vce(robust)"
		if "`v'" == "clu" local lbl "vce(cluster pu)"
		if "`v'" == "svy" local lbl "vce(svy)"
		qui su r_`v'
		di as txt "  `lbl'" _col(28) as res %9.3f r(mean) _col(42)	///
		   %9.3f r(min) _col(56) %9.3f r(max)
		local m_`v' = r(mean)
		local x_`v' = r(max)
	}

	di ""
	di as txt "  Les coefficients les plus touches :"
	gsort -r_rob
	di as txt "    coefficient" _col(22) "robust" _col(34)		///
	   "cluster" _col(46) "svy" _col(58) "boot" _col(70) "boot/svy"
	forvalues j = 1/5 {
		di as txt "    " as res %-14s nm[`j'] _col(20) %9.4f sr[`j']	///
		   _col(32) %9.4f sc[`j'] _col(44) %9.4f ss[`j']		///
		   _col(56) %9.4f sb[`j'] _col(68) %9.3f r_svy[`j']
	}
restore

*=========================================== les elasticites sous plan
* Les coefficients sont valides ; les elasticites ne le sont pas encore.  Elles
* heritent du plan PAR les coefficients, puisque la methode delta s'applique a
* V.  Mais une elasticite s'ecrit 1 + aleph/wbar, et wbar est une moyenne de
* population elle aussi estimee, elle aussi soumise au plan : notre jacobien la
* traite comme fixe.  Ce bloc mesure ce que cette omission coute.
capture program drop _elb
program define _elb, rclass
	version 14.2
	tempname EI EP
	qui easi $SH_ [pw=sweight], $SPEC_
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

qui easi `SH', lnprices(`PR') lnexpenditure(lx) demographics(z1 isMale)	///
	power(3) nolog vce(svy)
matrix EI = e(elast_exp)
matrix ES = e(elast_exp_se)
matrix EP = e(elast_price_nc)
matrix PS = e(elast_price_nc_se)

set seed 20260922
qui bootstrap `el', reps(`REPS') strata(st) cluster(pu) idcluster(_newpu2)	///
	nodots nowarn: _elb
matrix BE = e(V)

di ""
di as txt "  Elasticites : se(delta, vce(svy)) / se(bootstrap de plan)"
di as txt "  bien" _col(14) "part" _col(26) "estime" _col(38) "se delta"	///
   _col(50) "se boot" _col(62) "rapport"
local wmax 0
local pmax 0
forvalues j = 1/3 {
	qui su w`j' [aw=sweight]
	local wj = r(mean)
	local r1 = abs(ES[1,`j']) / sqrt(BE[`j',`j'])
	local wmax = max(`wmax', abs(`r1' - 1))
	di as txt "  depense `j'" _col(12) as res %9.4f `wj' _col(24)	///
	   %9.4f EI[1,`j'] _col(36) %9.4f abs(ES[1,`j'])		///
	   _col(48) %9.4f sqrt(BE[`j',`j']) _col(60) %9.3f `r1'
}
forvalues j = 1/3 {
	local k = 3 + `j'
	qui su w`j' [aw=sweight]
	local wj = r(mean)
	local r1 = abs(PS[`j',`j']) / sqrt(BE[`k',`k'])
	local pmax = max(`pmax', abs(`r1' - 1))
	di as txt "  prix `j'" _col(12) as res %9.4f `wj' _col(24)	///
	   %9.4f EP[`j',`j'] _col(36) %9.4f abs(PS[`j',`j'])		///
	   _col(48) %9.4f sqrt(BE[`k',`k']) _col(60) %9.3f `r1'
}

*--------------------------------------------------- verdict
di ""
local gap_rob = abs(`m_rob' - 1)
local gap_clu = abs(`m_clu' - 1)
di as txt "  ecart moyen de vce(robust) au plan" _col(60)		///
   as res %9.3f `gap_rob'
di as txt "  ecart moyen de vce(cluster) au plan" _col(60)		///
   as res %9.3f `gap_clu'
local gap_svy = abs(`m_svy' - 1)
di as txt "  ecart moyen de vce(svy) au plan" _col(60)			///
   as res %9.3f `gap_svy' _col(72)					///
   as res cond(`gap_svy' < 0.05, "ok", "ECHEC")
di ""
if `gap_clu' < `gap_rob' {
	di as txt "  -> vce(cluster pu) recupere l'essentiel de l'effet de grappe."
	di as txt "     Ce que vce(svy) ajouterait est la stratification, qui joue"
	di as txt "     en sens inverse et reduit la variance."
}
else {
	di as txt "  -> le clustering ne suffit pas a lui seul."
}
di ""
di as txt "  ecart max, elasticites-depense" _col(60) as res %9.3f `wmax'
di as txt "  ecart max, elasticites-prix propres" _col(60)		///
   as res %9.3f `pmax' _col(72)						///
   as res cond(`pmax' < 0.15, "ok", "ECHEC")
di as txt "{hline 78}"

log close ta
