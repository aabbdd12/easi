*! easi_tour.do -- tour d'horizon du module easi sur hixdata
*!
*! Trois temps :
*!   1. un modele simple, sortie par defaut
*!   2. le meme, avec les tables additionnelles (reporting)
*!   3. post-estimation : predict, test/lincom, courbes d'Engel
*!
*! hixdata : 4 847 menages canadiens, 9 biens, prix et depense DEJA NORMALISES
*! autour de la periode de base (log_y de moyenne -0,11).  Cette normalisation
*! n'est pas cosmetique : en niveaux bruts, les colonnes 1, y, y^2, y^3 du
*! design sont correlees a plus de 0,999 et la matrice normale devient
*! singuliere.  Si vous arrivez avec vos propres donnees, centrez.
*!
*! Rien ici n'a besoin de poids ni de plan de sondage : c'est le jeu de donnees
*! des estimateurs.  Pour la validation ponderee, voir examples/mex_bench.dta.

clear all
set more off

local ROOT "C:/Users/aabd/OneDrive/Desktop/EASI_project"
adopath ++ "`ROOT'/src"
cd "`ROOT'/examples"

use "`ROOT'/examples/hixdata.dta", clear

local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local NM food_home food_rest rent operation furniture clothing transport ///
	 recreation personal

di ""
di as txt "{hline 78}"
di as txt "  hixdata : " _N " menages, 9 biens"
di as txt "{hline 78}"
tabstat `SH', stat(mean sd min max) columns(statistics) format(%9.4f)


*==========================================================================
di ""
di as txt "{hline 78}"
di as txt "  1.  MODELE SIMPLE"
di as txt "{hline 78}"
*
* Tout par defaut : ecarts-types robustes, Sigma iteree, methode delta sur les
* elasticites, terme du regresseur genere inclus.  -snames()- ne sert qu'a
* etiqueter les tables.
*==========================================================================

easi `SH', lnprices(`PR') lnexpenditure(log_y)				///
	demographics(age hsex carown) power(3)				///
	snames(`NM') dec(3)

di ""
di as txt "  Lecture : les elasticites-depense sous 1 sont des necessites"
di as txt "  (alimentation a domicile 0,52 ; logement 0,70), celles au-dessus"
di as txt "  des biens superieurs (restaurants 1,54 ; ameublement 1,98)."
di as txt ""
di as txt "  La Table 03 se lit LIGNE = BIEN, COLONNE = PRIX : la diagonale"
di as txt "  porte les elasticites-prix propres, toutes negatives ici."


*==========================================================================
di ""
di as txt "{hline 78}"
di as txt "  2.  LE MEME, AVEC LES TABLES ADDITIONNELLES"
di as txt "{hline 78}"
*
* -compensated- ajoute les elasticites de Hicks ET calcule leurs ecarts-types
* (elles sont hors defaut parce que la sortie est deja longue, et les laisser
* hors defaut evite aussi leur cout).  -demoelast- ajoute les demographiques,
* -checks- les identites d'agregation.  -detail- fait les trois.
*==========================================================================

easi `SH', lnprices(`PR') lnexpenditure(log_y)				///
	demographics(age hsex carown) power(3)				///
	snames(`NM') dec(3) detail

di ""
di as txt "  Verification de Slutsky, a l'oeil : les elasticites-prix propres"
di as txt "  COMPENSEES (Table 04) sont systematiquement MOINS negatives que"
di as txt "  les non compensees (Table 03), puisque eta^H = eta^M + w * eta^x"
di as txt "  et que w * eta^x > 0 pour un bien normal."
di ""
di as txt "  Les identites d'agregation sont un autotest du code d'elasticites"
di as txt "  execute sur VOS donnees.  Elles sont toujours calculees ; un echec"
di as txt "  serait signale meme sans l'option -checks-."
di as txt "     e(chk_engel)   = " as res %11.3e e(chk_engel)
di as txt "     e(chk_cournot) = " as res %11.3e e(chk_cournot)

di ""
di as txt "  Les memes options se donnent au rejeu, sans reestimer :"
di as txt "     . easi, compensated"

di ""
di as txt "  Matrices stockees (toutes en [bien, prix]) :"
di as txt "     e(elast_exp)       e(elast_exp_se)"
di as txt "     e(elast_price_nc)  e(elast_price_nc_se)   non compensees"
di as txt "     e(elast_price_c)   e(elast_price_c_se)    compensees"
di as txt "     e(elast_demo)      e(elast_demo_se)"
di as txt "     e(slutsky)  e(semi_exp)  e(semi_price)  e(Sigma)"
matrix list e(elast_price_c), format(%9.4f) title("  e(elast_price_c)")


*==========================================================================
di ""
di as txt "{hline 78}"
di as txt "  3.  POST-ESTIMATION"
di as txt "{hline 78}"
*==========================================================================

*---------------------------------------------------------------- predict
di ""
di as txt "  3a. predict : parts ajustees, indice d'utilite implicite, residus"

predict double wh*, shares
predict double yhat, y
predict double rr*, residuals

di ""
di as txt "     les parts ajustees somment-elles a 1 ?"
qui gen double wsum = wh1+wh2+wh3+wh4+wh5+wh6+wh7+wh8+wh9
qui su wsum
di as txt "       min " as res %12.10f r(min) as txt "   max " as res %12.10f r(max)

di ""
di as txt "     l'indice d'utilite implicite y :"
qui su yhat
di as txt "       moyenne " as res %8.4f r(mean) as txt "   ecart-type " ///
   as res %8.4f r(sd)
di as txt "     (y n'est PAS log_y : il corrige la depense par le terme de prix"
di as txt "      p'A(z)p/2, ce qui en fait une mesure de depense reelle)"
qui corr yhat log_y
di as txt "       correlation avec log_y = " as res %6.4f r(rho)

di ""
di as txt "     residus : observe moins ajuste, un par bien"
qui su rr1
di as txt "       rr1 : moyenne " as res %11.3e r(mean) as txt "   ecart-type " ///
   as res %7.4f r(sd)

*------------------------------------------------------ test et lincom
di ""
di as txt "  3b. easi est e-class : test, lincom, nlcom fonctionnent"
di ""
di as txt "     Les termes en y de l'equation du logement sont-ils tous nuls ?"
di as txt "     (c'est-a-dire : la courbe d'Engel du logement est-elle plate ?)"
test [srent]y1 [srent]y2 [srent]y3

*------------------------------------------------------ courbes d'Engel
di ""
di as txt "  3c. estat engel : les courbes d'Engel"
di ""
di as txt "     Par defaut (-atmeans-) : la part ajustee en fonction de la"
di as txt "     depense totale, demographiques et prix tenus a leur moyenne"
di as txt "     ponderee, le point fixe (w, y) resolu a chaque noeud.  C'est"
di as txt "     donc la fonction exacte : n() est une RESOLUTION, pas un"
di as txt "     parametre de lissage, et bwidth() n'a pas de sens ici."

estat engel, n(60) saving("engel_atmeans.gph", replace)
return list

qui graph use "engel_atmeans.gph"
qui graph export "engel_atmeans.png", replace width(1100)
di as txt "     -> engel_atmeans.png"

di ""
di as txt "     Avec -asobserved- : covariables laissees telles qu'observees."
di as txt "     C'est un nuage de points, donc il EST lisse, par une regression"
di as txt "     locale lineaire ; bwidth() s'applique alors, avec la regle du"
di as txt "     pouce des polynomes locaux par defaut (et non celle de"
di as txt "     Silverman, qui est une regle de DENSITE)."

estat engel, asobserved n(60) saving("engel_asobserved.gph", replace)
qui graph use "engel_asobserved.gph"
qui graph export "engel_asobserved.png", replace width(1100)
di as txt "     -> engel_asobserved.png"


*==========================================================================
di ""
di as txt "{hline 78}"
di as txt "  ANNEXE : le mode -compat- reproduit le paquet R, bugs compris"
di as txt "{hline 78}"
*==========================================================================

qui easi `SH', lnprices(`PR') lnexpenditure(log_y)			///
	demographics(age hsex carown) power(3) snames(`NM') nolog
matrix COR = e(elast_exp)
qui easi `SH', lnprices(`PR') lnexpenditure(log_y)			///
	demographics(age hsex carown) power(3) snames(`NM') nolog compat
matrix CMP = e(elast_income)

matrix CC = CMP \ COR \ (COR - CMP)
matrix rownames CC = paquet_R corrige ecart
di ""
di as txt "  Elasticites-depense : paquet R (compat) contre mode corrige"
matlist CC, format(%9.4f) twidth(12)
di ""
di as txt "  Le controle d'agregation detecte le defaut TOUT SEUL, et isole"
di as txt "  sa cause : c'est -legacy(elast)-, les formules d'elasticite, et"
di as txt "  rien d'autre.  Aucune difference finie, aucun test d'homogeneite,"
di as txt "  et cela marche sur n'importe quel jeu de donnees."
di ""
di as txt "     mode" _col(28) "Engel" _col(45) "Cournot"
foreach m in "" "compat" "legacy(elast)" "legacy(ez)" {
	qui easi `SH', lnprices(`PR') lnexpenditure(log_y)		///
		demographics(age hsex carown) power(3) nolog noelastse `m'
	local lbl = cond("`m'" == "", "corrige (defaut)", "`m'")
	di as txt "     `lbl'" _col(26) as res %11.3e e(chk_engel)	///
	   _col(43) %11.3e e(chk_cournot)
}
di ""
di as txt "  L'ecart vient des formules d'elasticite : la derivation du paquet"
di as txt "  R differencie les parts a parts CONSTANTES, donc elle manque le"
di as txt "  fait que y repond aux prix et a la depense.  Les formules"
di as txt "  corrigees collent aux differences finies a 7e-09, celles du"
di as txt "  paquet R a 1,9e-02."

di ""
di as txt "{hline 78}"
