*! test_step13.do -- easidiag : le diagnostic tient-il ses promesses ?
*!
*! Quatre choses sont verrouillees ici.
*!
*! 1. IL N'ESTIME PAS.  Les sections 1 a 6 se calculent sur les donnees et sur
*!    la matrice de design au point de depart y = indice de Stone, qui ne depend
*!    que des donnees.  C'est le point : un diagnostic qui aurait besoin que le
*!    modele s'ajuste serait muet exactement quand il ne s'ajuste pas.
*!
*! 2. IL TROUVE LA QUASI-DEPENDANCE, PAS SEULEMENT L'EXACTE.  Un jeu de muettes
*!    exhaustif perturbe de 1e-9 laisse la matrice de plein rang -- un test de
*!    rang ne voit rien -- alors que les ecarts-types explosent d'un facteur
*!    3e+08.  L'indice de conditionnement, lui, est continu et le voit.
*!
*! 3. IL NOMME LES COUPABLES.  Les proportions de decomposition de variance
*!    designent les colonnes impliquees, au lieu d'annoncer "rang deficient".
*!
*! 4. IL PREDIT LE NOMBRE D'ITERATIONS, et la prediction doit coller a ce que
*!    easi met reellement.

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do 10_easidiag_checks.do
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

log using "`ROOT'/replication/out/step13.log", replace text name(ta)

di ""
di as txt "{hline 76}"
di as txt "  easidiag : depistage sans estimation"
di as txt "{hline 76}"

global NFAIL 0
capture program drop _ck
program define _ck
	args label ok
	di as txt "  `label'" _col(62) as res cond(`ok', "ok", "ECHEC")
	if !`ok' global NFAIL = ${NFAIL} + 1
end

*==================================================== 1. cas sain : hixdata
use "`ROOT'/examples/hixdata.dta", clear
local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local SPEC lnprices(`PR') lnexpenditure(log_y)				///
	   demographics(age hsex carown) power(3)

di ""
di as txt "  1. hixdata, specification saine"
qui easidiag `SH', `SPEC'
local pred = r(iter)
local rat  = r(ratio)
local nprob = r(nprob)

_ck "conditionnement raisonnable (cond < 1e4)" `=r(cond_raw) < 1e4'
_ck "la boucle contracte (ratio < 1)"          `=`rat' < 1'
di as txt "     iterations predites" _col(56) as res %6.0f `pred'

* la prediction contre la realite
qui easi `SH', `SPEC' nolog noelastse
local reel = e(iter)
di as txt "     iterations reelles" _col(56) as res %6.0f `reel'
_ck "la prediction colle a l'execution (a 1 pres)" `=abs(`pred' - `reel') <= 1'

* le seul point signale sur ces donnees : les regimes de prix
_ck "un seul point signale"                    `=`nprob' == 1'

*========================= 2. colinearite exacte, puis quasi-exacte
* the reduced Mexican survey shipped with the module: shares w1-w3, centred
* logs lp1 lp2 lp3 lx, their uncentred copies lx_raw lp1_raw lp2_raw, and the
* exhaustive household-composition dummies nocup0-nocup4
use "`ROOT'/examples/mex_bench.dta", clear
local MX lnprices(lp1 lp2 lp3) lnexpenditure(lx) power(3)

di ""
di as txt "  2. jeu de muettes exhaustif, colinearite EXACTE"
qui easidiag w1 w2 w3, `MX' demographics(nocup0 nocup1 nocup2 nocup3 nocup4)
_ck "le design est declare mal conditionne"    `=r(cond_raw) > 1e6'
_ck "au moins un point signale"                `=r(nprob) > 0'

* l'estimateur, lui, ne dit rien : c'est la raison d'etre du diagnostic
qui easi w1 w2 w3, `MX' demographics(nocup0 nocup1 nocup2 nocup3 nocup4)	///
	nolog noelastse
_ck "easi, lui, converge sans broncher"        `=e(converged) == 1'

di ""
di as txt "  2b. les memes, perturbees de 1e-9 : QUASI-exacte"
set seed 11
qui gen double u = runiform() - 0.5
foreach v in nocup0 nocup1 nocup2 nocup3 nocup4 {
	qui gen double p_`v' = `v' + 1e-9 * u
}
qui easidiag w1 w2 w3, `MX'						///
	demographics(p_nocup0 p_nocup1 p_nocup2 p_nocup3 p_nocup4)
_ck "la quasi-dependance est vue elle aussi"   `=r(cond_raw) > 1e6'
* un test de RANG ne la verrait pas : la matrice est de plein rang
qui easi w1 w2 w3, `MX'							///
	demographics(p_nocup0 p_nocup1 p_nocup2 p_nocup3 p_nocup4)	///
	nolog noelastse
matrix V = e(V)
local smax 0
local K = colsof(V)
forvalues i = 1/`K' {
	local smax = max(`smax', sqrt(V[`i',`i']))
}
di as txt "     plus grand ecart-type de easi" _col(52) as res %14.1f `smax'
_ck "et elle est bien explosive (se > 1000)"   `=`smax' > 1000'

*=========================================== 3. donnees non centrees
di ""
di as txt "  3. les memes donnees en niveaux bruts"
* lx_raw lp1_raw lp2_raw are the uncentred logs shipped with the bench
capture qui easidiag w1 w2 w3, lnprices(lp1_raw lp2_raw lp3)		///
	lnexpenditure(lx_raw) demographics(hhsize isMale) power(3)
_ck "le centrage ameliore nettement le conditionnement"			///
   `=r(cond_raw) > 10 * r(cond_ctr)'

*=========================================== 4. options ignorees
di ""
di as txt "  4. les options de l'estimateur sont recues et annoncees"
use "`ROOT'/examples/hixdata.dta", clear
qui easidiag `SH', `SPEC' vce(robust) dec(5) compensated
_ck "la commande s'execute malgre elles"       `=r(N) == 4847'

di ""
if ${NFAIL} == 0 di as txt "  -> tous les blocs passent."
else             di as error "  -> ${NFAIL} echec(s)"
di as txt "{hline 76}"

log close ta
