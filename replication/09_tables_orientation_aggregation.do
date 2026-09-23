*! test_step12.do -- tables d'elasticites : orientation, compensees, agregation
*!
*! Trois choses sont verrouillees ici.
*!
*! 1. ORIENTATION.  e(elast_price) est indexee [prix, bien] et e(compensated_q)
*!    [bien, prix] : les deux matrices heritees du paquet R se lisent en sens
*!    OPPOSES, et comme leurs deux axes portent les memes noms de biens, rien
*!    ne signale au lecteur qu'il doit transposer.  Les noms courants
*!    e(elast_price_nc) et e(elast_price_c) sont tous deux [bien, prix] : la
*!    ligne j, colonne k est l'elasticite du bien j au prix du bien k.
*!
*! 2. COMPENSEES.  e(compensated_q) vaut eta^H + I, pas eta^H : la matrice de
*!    Slutsky en parts est EPS = EP + w'w - diag(w), donc diag(w)^-1 EPS =
*!    compensated_q - I.  Retrancher l'identite n'est pas cosmetique, c'est
*!    toute la diagonale, c'est-a-dire chaque elasticite-prix propre.
*!
*! 3. AGREGATION.  Engel sum_j w_j eta^x_j = 1 et Cournot sum_j w_j eta^k_j =
*!    -w_k sont exactes pour les formules en moyennes, parce que sum_j aleph_j
*!    = 0 et sum_j a^jk = 0.  Elles tournent sur les donnees de l'utilisateur :
*!    c'est un autotest du code d'elasticites, pas une statistique descriptive.
*!
*! Donnees : hixdata (pas de poids, pas de plan de sondage -- le banc mexicain
*! est reserve a la validation ponderee).

clear all
set more off

local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"

log using "`ROOT'/replication/out/step12.log", replace text name(ta)

use "`ROOT'/examples/hixdata.dta", clear
local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local SPEC lnprices(`PR') lnexpenditure(log_y)				///
	   demographics(age hsex carown) power(3) nolog
local J 9

di ""
di as txt "{hline 76}"
di as txt "  Tables d'elasticites : orientation, compensees, agregation"
di as txt "{hline 76}"

local fail 0
capture program drop _chk
program define _chk
	args label value tol
	local tag = cond(`value' <= `tol', "ok", "ECHEC")
	di as txt "  `label'" _col(52) as res %11.3e `value' _col(66) as res "`tag'"
	if `value' > `tol' {
		global NFAIL = ${NFAIL} + 1
	}
end
global NFAIL 0

*------------------------------------------------ 1. orientation et compensees
qui easi `SH', `SPEC' compensated
matrix NC  = e(elast_price_nc)
matrix CC  = e(elast_price_c)
matrix EX  = e(elast_exp)
matrix OLD = e(elast_price)
matrix CQ  = e(compensated_q)
matrix SL  = e(slutsky)

di ""
di as txt "  1. Orientation et definition"
mata:
    nc = st_matrix("NC"); cc = st_matrix("CC"); ex = st_matrix("EX")
    old = st_matrix("OLD"); cq = st_matrix("CQ"); sl = st_matrix("SL")
    Jg  = cols(ex)
    ws  = mean(st_data(., st_local("SH")))
    // elast_price_nc est bien la transposee de l'ancienne
    st_numscalar("d_t", max(abs(nc - old')))
    // elast_price_c est bien compensated_q - I
    st_numscalar("d_c", max(abs(cc - (cq - I(Jg)))))
    // ... et diag(w)^-1 * Slutsky
    st_numscalar("d_s", max(abs(cc - diag(1:/ws) * sl)))
    // homogeneite de la matrice compensee : chaque ligne somme a zero
    st_numscalar("d_h", max(abs(rowsum(cc))))
    // Slutsky : eta^M[j,k] = eta^H[j,k] - w_k eta^x_j, aux covariances pres
    D = J(Jg, Jg, 0)
    for (j = 1; j <= Jg; j++) {
        for (k = 1; k <= Jg; k++) D[j,k] = nc[j,k] - (cc[j,k] - ws[k]*ex[j])
    }
    st_numscalar("d_sl", max(abs(D)))
    // les compensees propres sont moins negatives que les non compensees
    st_numscalar("d_own", max(diagonal(nc) - diagonal(cc)))
end
_chk "elast_price_nc == transposee de elast_price" d_t 1e-14
_chk "elast_price_c  == compensated_q - I"         d_c 1e-14
_chk "elast_price_c  == diag(w)^-1 * slutsky"      d_s 1e-12
_chk "homogeneite : chaque ligne compensee somme a 0" d_h 1e-8
di as txt "  Slutsky eta^M = eta^H - w_k eta^x (residu = covariance)" ///
   _col(52) as res %11.3e d_sl
di as txt "  max(eta^M propre - eta^H propre) < 0 attendu" _col(52)	///
   as res %11.3e d_own _col(66) as res cond(d_own < 0, "ok", "ECHEC")
if d_own >= 0 global NFAIL = ${NFAIL} + 1

*------------------------------------------------ 2. identites d'agregation
di ""
di as txt "  2. Identites d'agregation, postees dans e()"
_chk "Engel   |sum_j w_j eta^x_j - 1|"       "e(chk_engel)"   1e-8
_chk "Cournot max_k |sum_j w_j eta^k_j + w_k|" "e(chk_cournot)" 1e-8

*------------------------------------------------ 3. options de reporting
di ""
di as txt "  3. Les tables additionnelles sont bien sous option"
qui easi `SH', `SPEC' noelastse
local r0 "`e(report)'"
qui easi `SH', `SPEC' noelastse compensated
local r1 "`e(report)'"
qui easi `SH', `SPEC' noelastse detail
local r2 "`e(report)'"
di as txt "     defaut      : report = " as res `""`r0'""'		///
   _col(52) as res cond("`r0'" == "", "ok", "ECHEC")
di as txt "     compensated : report = " as res `""`r1'""'		///
   _col(52) as res cond("`r1'" == "compensated", "ok", "ECHEC")
di as txt "     detail      : report = " as res `""`r2'""'		///
   _col(52) as res cond(strpos("`r2'","compensated") & strpos("`r2'","demoelast") ///
   & strpos("`r2'","checks"), "ok", "ECHEC")
if "`r0'" != ""            global NFAIL = ${NFAIL} + 1
if "`r1'" != "compensated" global NFAIL = ${NFAIL} + 1

* sans -compensated-, les SE des compensees ne sont pas calculees
matrix CS0 = e(elast_price_c_se)
qui easi `SH', `SPEC' compensated
matrix CS1 = e(elast_price_c_se)
local m0 = missing(CS0[1,1])
local m1 = missing(CS1[1,1])
di as txt "     SE compensees absentes par defaut" _col(52)		///
   as res cond(`m0', "ok", "ECHEC")
di as txt "     SE compensees presentes sous option" _col(52)		///
   as res cond(!`m1', "ok", "ECHEC")
if !`m0' global NFAIL = ${NFAIL} + 1
if `m1'  global NFAIL = ${NFAIL} + 1

*------------------------------------------------ 4. noms hérités intacts
di ""
di as txt "  4. Les noms herites gardent leur contenu"
qui easi `SH', `SPEC' noelastse
foreach m in elast_income elast_price compensated_q slutsky semi_income {
	capture confirm matrix e(`m')
	local ok = (_rc == 0)
	di as txt "     e(`m')" _col(52) as res cond(`ok', "present", "ABSENT")
	if !`ok' global NFAIL = ${NFAIL} + 1
}
matrix A1 = e(elast_income)
matrix A2 = e(elast_exp)
mata: st_numscalar("d_a", max(abs(st_matrix("A1") - st_matrix("A2"))))
_chk "elast_exp == elast_income (meme contenu)" d_a 0

di ""
if ${NFAIL} == 0 di as txt "  -> les 4 blocs passent."
else             di as error "  -> ${NFAIL} echec(s)"
di as txt "{hline 76}"

log close ta
