*! test_step16.do -- jacobien analytique des elasticites
*!
*! Trois choses sont verrouillees ici.
*!
*! 1. ANALYTIQUE == DIFFERENCES FINIES.  Les ecarts-types des elasticites sont
*!    un delta method, se = sqrt(g'Vg).  Le gradient g etait obtenu par
*!    differences finies sur _easi_epoint (K evaluations, 86 % du temps
*!    d'execution) ; il est maintenant analytique (_easi_ejac).  L'oracle est
*!    l'ancien jacobien numerique, que le scalaire _easi_jaccheck fait calculer
*!    a cote du nouveau.  Les deux doivent coincider au plancher des differences
*!    finies, famille par famille -- EI, EPRICE, EZ, EPQ -- et sur toutes les
*!    interactions (py, pz, zy), parce que c'est la ou les derivees a la main se
*!    trompent : la completion par adding-up de bjk et des coefficients de prix
*!    n'est pas symetrique, et un gradient ecrit en coordonnees libres qui
*!    l'ignore est faux de 1.6 % avec py (mesure, hixdata).  Le code n'ecrit
*!    aucune derivee en coordonnees libres : il derive en coordonnees completees
*!    et applique la completion comme une matrice, G = Gc * C.
*!
*! 2. CELLULES DES SE DEMOGRAPHIQUES.  Sans pz ni zy, EZ[t,i] est exactement le
*!    coefficient de z_t dans l'equation i, donc son ecart-type est exactement
*!    sqrt(V) de ce coefficient, cellule par cellule, et pour le bien omis
*!    sqrt(1'V1).  Cet oracle ne passe par aucun empilement.  Il attrape le bug
*!    de _easi_estack qui empilait vec(EZ') la ou le depilage attend vec(EZ) :
*!    avec T = 3 et 9 biens, la cellule (2,1) de e(elast_demo_se) affichait
*!    l'ecart-type de la cellule (1,2).
*!
*! 3. LES SE N'ONT PAS BOUGE.  Les tables d'ecarts-types de la specification de
*!    reference sont identiques a celles du jacobien numerique, a la precision
*!    des differences finies pres : le remplacement est transparent.
*!
*! Donnees : hixdata et le banc mexicain (celui-ci pour pz avec un plan
*! stratifie, et parce que ses prix sont centres : P_J n'y est pas nul par
*! observation, ce qui est exactement le terme qui manquait).

clear all
set more off

* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Section4-4_analytic_jacobian.do
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

log using "`ROOT'/replication/out/step16.log", replace text name(ta)

mata:
// ecart relatif max entre le jacobien analytique et le numerique, par famille
real rowvector _jac_gap(real scalar Jg, real scalar T, real scalar doCQ)
{
	real matrix a, f, af, ff
	real rowvector lo, hi, out
	real scalar i, nf
	a = st_matrix("_easi_Gan"); f = st_matrix("_easi_Gfd")
	lo = (1, Jg + 1, Jg + Jg^2 + 1, Jg + Jg^2 + T*Jg + 1)
	hi = (Jg, Jg + Jg^2, Jg + Jg^2 + T*Jg, Jg + 2*Jg^2 + T*Jg)
	nf = doCQ ? 4 : 3
	out = J(1, 4, .)
	for (i = 1; i <= nf; i++) {
		af = a[| lo[i], 1 \ hi[i], . |]
		ff = f[| lo[i], 1 \ hi[i], . |]
		out[i] = max(abs(af - ff)) / max((max(abs(ff)), 1e-12))
	}
	return(out)
}
end

* ======================================================================
* 1. analytique == differences finies, famille par famille
* ======================================================================
di ""
di as txt "{hline 70}"
di as txt "1. Analytic Jacobian against the finite-difference one"
di as txt "{hline 70}"
scalar _easi_jaccheck = 1

* les familles non lineaires (EI, EPRICE) sont bornees par le bruit des
* differences finies ; les lineaires (EZ, EPQ) par la precision machine
local TOL_NL 1e-5
local TOL_L  1e-9

use "`ROOT'/examples/mex_bench.dta", clear
local MSPEC lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 isMale) power(3) nolog
qui easi w1 w2 w3, `MSPEC'
mata: st_matrix("gap", _jac_gap(3, 2, 0))
di as txt "   mex_bench, no interactions      EI " %9.2e gap[1,1] "  EPRICE " %9.2e gap[1,2] "  EZ " %9.2e gap[1,3]
assert gap[1,1] < `TOL_NL' & gap[1,2] < `TOL_NL' & gap[1,3] < `TOL_L'

qui easi w1 w2 w3, `MSPEC' py pz zy compensated
mata: st_matrix("gap", _jac_gap(3, 2, 1))
di as txt "   mex_bench, py pz zy, compensated EI " %9.2e gap[1,1] "  EPRICE " %9.2e gap[1,2] "  EZ " %9.2e gap[1,3] "  EPQ " %9.2e gap[1,4]
assert gap[1,1] < `TOL_NL' & gap[1,2] < `TOL_NL' & gap[1,3] < `TOL_L' & gap[1,4] < `TOL_L'

qui easi w1 w2 w3 [pw=sweight], `MSPEC' py vce(svy)
mata: st_matrix("gap", _jac_gap(3, 2, 0))
di as txt "   mex_bench, pweight + vce(svy)   EI " %9.2e gap[1,1] "  EPRICE " %9.2e gap[1,2] "  EZ " %9.2e gap[1,3]
assert gap[1,1] < `TOL_NL' & gap[1,2] < `TOL_NL' & gap[1,3] < `TOL_L'

use "`ROOT'/examples/hixdata.dta", clear
local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local HSPEC lnprices(`PR') lnexpenditure(log_y)

qui easi `SH', `HSPEC' demographics(age hsex carown) power(3) py nolog
mata: st_matrix("gap", _jac_gap(9, 3, 0))
di as txt "   hixdata, py                     EI " %9.2e gap[1,1] "  EPRICE " %9.2e gap[1,2] "  EZ " %9.2e gap[1,3]
assert gap[1,1] < `TOL_NL' & gap[1,2] < `TOL_NL' & gap[1,3] < `TOL_L'

qui easi `SH', `HSPEC' demographics(age hsex carown time tran) power(5) py pz zy compensated nolog
mata: st_matrix("gap", _jac_gap(9, 5, 1))
di as txt "   hixdata, reference spec         EI " %9.2e gap[1,1] "  EPRICE " %9.2e gap[1,2] "  EZ " %9.2e gap[1,3] "  EPQ " %9.2e gap[1,4]
assert gap[1,1] < `TOL_NL' & gap[1,2] < `TOL_NL' & gap[1,3] < `TOL_L' & gap[1,4] < `TOL_L'

* les tables d'ecarts-types de la specification de reference, telles que
* rapportees, contre les memes tables recalculees avec le jacobien numerique
tempname Van Vfd Sfd
matrix `Van' = e(V)
mata: st_matrix("`Sfd'", sqrt(diagonal(st_matrix("_easi_Gfd") * st_matrix("`Van'") * st_matrix("_easi_Gfd")'))')
matrix PSE = e(elast_price_se)
matrix ZSE = e(elast_demo_se)
mata: pse = st_matrix("PSE"); zse = st_matrix("ZSE"); sfd = st_matrix("`Sfd'")
mata: st_numscalar("dP", max(abs(vec(pse) - sfd[| 10 \ 90 |]') :/ sfd[| 10 \ 90 |]'))
mata: st_numscalar("dZ", max(abs(vec(zse) - sfd[| 91 \ 135 |]') :/ sfd[| 91 \ 135 |]'))
di as txt "   reported SE vs numerical-Jacobian SE:  price " %9.2e dP "   demographic " %9.2e dZ
assert dP < 1e-5 & dZ < 1e-5

scalar drop _easi_jaccheck
di as res "   OK"

* ======================================================================
* 2. cellules des SE demographiques : oracle sans empilement
* ======================================================================
di ""
di as txt "{hline 70}"
di as txt "2. Demographic-elasticity standard errors, cell by cell"
di as txt "{hline 70}"
qui easi `SH', `HSPEC' demographics(age hsex carown) power(3) nolog
matrix ZSE = e(elast_demo_se)
matrix V   = e(V)
local T 3
local J 9
local k = colsof(V) / (`J' - 1)
tempname one
local worst 0
forvalues t = 1/`T' {
	forvalues i = 1/8 {
		local p = (`i' - 1) * `k' + 1 + 3 + `t'		// power(3): 1 + R + t
		local se = sqrt(V[`p', `p'])
		local d = abs(ZSE[`t', `i'] - `se') / `se'
		if `d' > `worst' local worst `d'
	}
	* le bien omis : -sum_i gjt[t,i], donc sqrt(1'V1) sur les 8 positions
	matrix `one' = J(1, colsof(V), 0)
	forvalues i = 1/8 {
		matrix `one'[1, (`i' - 1) * `k' + 1 + 3 + `t'] = 1
	}
	matrix VJ = `one' * V * `one''
	local se = sqrt(VJ[1, 1])
	local d = abs(ZSE[`t', `J'] - `se') / `se'
	if `d' > `worst' local worst `d'
}
di as txt "   worst relative deviation over " `T' " x " `J' " cells: " %9.2e `worst'
assert `worst' < 1e-8
di as res "   OK -- each cell of e(elast_demo_se) is the standard error of its own cell"

log close ta
di ""
di as res "test_step16: all assertions passed"
