*! easi_tour.do -- a tour of easi on hixdata
*!
*! Three steps:
*!   1. a simple model, the default output
*!   2. the same, with the additional tables
*!   3. after estimation: predict, test, Engel curves
*! and an annex: compat reproduces the R package, defects included.
*!
*! hixdata: 4,847 Canadian households, 9 goods, prices and expenditure already
*! normalized around the base period (mean log_y -0.11).  The normalization is
*! not cosmetic: in raw levels, the columns 1, y, y^2, y^3 of the design are
*! correlated above 0.999 and the normal matrix becomes singular.  With your
*! own data, centre them.
*!
*! Nothing here needs weights or a survey design: these are the data of the
*! estimators.  For a survey with its design, see mex_tour.do on mex_bench.dta.
*!
*! Run it from the folder where "ssc install easi, all" (or "net get easi")
*! copied it with hixdata.dta.  The graphs are written to that folder.

clear all
set more off

capture confirm file "hixdata.dta"
if _rc {
	di as err "hixdata.dta is not in the current folder: copy it with"
	di as err "  ssc install easi, all replace   (or net get easi)"
	exit 601
}
use hixdata, clear

local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local NM food_home food_rest rent operation furniture clothing transport ///
	 recreation personal

di ""
di as txt "{hline 78}"
di as txt "  hixdata: " _N " households, 9 goods"
di as txt "{hline 78}"
tabstat `SH', stat(mean sd min max) columns(statistics) format(%9.4f)


*==========================================================================
di ""
di as txt "{hline 78}"
di as txt "  1.  A SIMPLE MODEL"
di as txt "{hline 78}"
*
* Everything by default: robust standard errors, iterated Sigma, the
* elasticities with their influence-function standard errors, the term of the
* generated regressor included.  -snames()- only labels the tables.
*==========================================================================

easi `SH', lnprices(`PR') lnexpenditure(log_y)				///
	demographics(age hsex carown) power(3)				///
	snames(`NM') dec(3)

tempname EX
matrix `EX' = e(elast_exp)
di ""
di as txt "  Reading: expenditure elasticities below 1 are necessities (food at"
di as txt "  home " as res %4.2f el(`EX', 1, 1) as txt ", rent " as res %4.2f el(`EX', 1, 3) ///
   as txt "), those above 1 luxuries (restaurants " as res %4.2f el(`EX', 1, 2) ///
   as txt ", furniture " as res %4.2f el(`EX', 1, 5) as txt ")."
di as txt ""
di as txt "  The price-elasticity table reads ROW = GOOD, COLUMN = PRICE: the"
di as txt "  diagonal holds the own-price elasticities."


*==========================================================================
di ""
di as txt "{hline 78}"
di as txt "  2.  THE SAME, WITH THE ADDITIONAL TABLES"
di as txt "{hline 78}"
*
* -compensated- adds the Hicksian elasticities and their standard errors
* (not by default: the output is already long).  -demoelast- adds the
* demographic elasticities, -checks- the aggregation identities.  -detail-
* does all three.
*==========================================================================

easi `SH', lnprices(`PR') lnexpenditure(log_y)				///
	demographics(age hsex carown) power(3)				///
	snames(`NM') dec(3) detail

di ""
di as txt "  Slutsky at a glance: the COMPENSATED own-price elasticities are"
di as txt "  less negative than the uncompensated ones, since"
di as txt "  eta^H = eta^M + w * eta^x and w * eta^x > 0 for a normal good."
di ""
di as txt "  The aggregation identities are a self-test of the elasticity code, run"
di as txt "  on YOUR data.  They are always computed; a failure would be reported"
di as txt "  even without the option -checks-."
di as txt "     e(chk_engel)   = " as res %11.3e e(chk_engel)
di as txt "     e(chk_cournot) = " as res %11.3e e(chk_cournot)

di ""
di as txt "  The same options work on replay, without estimating again:"
di as txt "     . easi, compensated"

di ""
di as txt "  Stored matrices (all [good, price]):"
di as txt "     e(elast_exp)       e(elast_exp_se)"
di as txt "     e(elast_price_nc)  e(elast_price_nc_se)   uncompensated"
di as txt "     e(elast_price_c)   e(elast_price_c_se)    compensated"
di as txt "     e(elast_demo)      e(elast_demo_se)"
di as txt "     e(slutsky)  e(semi_exp)  e(semi_price)  e(Sigma)"
matrix list e(elast_price_c), format(%9.4f) title("  e(elast_price_c)")


*==========================================================================
di ""
di as txt "{hline 78}"
di as txt "  3.  AFTER ESTIMATION"
di as txt "{hline 78}"
*==========================================================================

*---------------------------------------------------------------- predict
di ""
di as txt "  3a. predict: fitted shares, the implicit utility index, residuals"

predict double wh*, shares
predict double yhat, y
predict double rr*, residuals

di ""
di as txt "     do the fitted shares add up to 1?"
qui gen double wsum = wh1+wh2+wh3+wh4+wh5+wh6+wh7+wh8+wh9
qui su wsum
di as txt "       min " as res %12.10f r(min) as txt "   max " as res %12.10f r(max)

di ""
di as txt "     the implicit utility index y:"
qui su yhat
di as txt "       mean " as res %8.4f r(mean) as txt "   standard deviation " ///
   as res %8.4f r(sd)
di as txt "     (y is NOT log_y: it corrects expenditure by the price term"
di as txt "      p'A(z)p/2, which makes it a measure of real expenditure)"
qui corr yhat log_y
di as txt "       correlation with log_y = " as res %6.4f r(rho)

di ""
di as txt "     residuals: observed minus fitted, one per good"
qui su rr1
di as txt "       rr1: mean " as res %11.3e r(mean) as txt "   standard deviation " ///
   as res %7.4f r(sd)

*------------------------------------------------------ test
di ""
di as txt "  3b. easi is e-class: test, lincom, nlcom work"
di ""
di as txt "     Are all the y terms of the rent equation zero?"
di as txt "     (that is: is the Engel curve of rent flat?)"
test [srent]y1 [srent]y2 [srent]y3

*------------------------------------------------------ Engel curves
di ""
di as txt "  3c. estat engel: the Engel curves"
di ""
di as txt "     By default (-atmeans-): the fitted share as a function of total"
di as txt "     expenditure, demographics and prices held at their weighted"
di as txt "     means, the fixed point (w, y) solved at each node.  It is the"
di as txt "     exact function: n() is a RESOLUTION, not a smoothing parameter,"
di as txt "     and bwidth() has no meaning here."

estat engel, n(60) saving("engel_atmeans.gph", replace)
return list

qui graph use "engel_atmeans.gph"
qui graph export "engel_atmeans.png", replace width(1100)
di as txt "     -> engel_atmeans.png"

di ""
di as txt "     With -asobserved-: the covariates as observed.  It is a scatter,"
di as txt "     so it IS smoothed, by a local linear regression; bwidth() then"
di as txt "     applies, with the rule of thumb of local polynomials by default"
di as txt "     (not Silverman's, which is a rule for DENSITIES)."

estat engel, asobserved n(60) saving("engel_asobserved.gph", replace)
qui graph use "engel_asobserved.gph"
qui graph export "engel_asobserved.png", replace width(1100)
di as txt "     -> engel_asobserved.png"


*==========================================================================
di ""
di as txt "{hline 78}"
di as txt "  ANNEX: -compat- reproduces the R package, defects included"
di as txt "{hline 78}"
*==========================================================================

qui easi `SH', lnprices(`PR') lnexpenditure(log_y)			///
	demographics(age hsex carown) power(3) snames(`NM') nolog
matrix COR = e(elast_exp)
local e0 = e(chk_engel)
local c0 = e(chk_cournot)
qui easi `SH', lnprices(`PR') lnexpenditure(log_y)			///
	demographics(age hsex carown) power(3) snames(`NM') nolog compat
matrix CMP = e(elast_income)
local e1 = e(chk_engel)
local c1 = e(chk_cournot)

matrix CC = CMP \ COR \ (COR - CMP)
matrix rownames CC = R_package corrected difference
di ""
di as txt "  Expenditure elasticities: the R package (compat) against the default"
matlist CC, format(%9.4f) twidth(12)
di ""
di as txt "  The aggregation checks detect the defect by themselves, on any data,"
di as txt "  without finite differences:"
di as txt "     mode" _col(28) "Engel" _col(45) "Cournot"
di as txt "     default" _col(26) as res %11.3e `e0' _col(43) %11.3e `c0'
di as txt "     compat"  _col(26) as res %11.3e `e1' _col(43) %11.3e `c1'
di ""
di as txt "  The gap comes from the elasticity formulas: the derivation of the R"
di as txt "  package differentiates the shares at CONSTANT y, so it misses that y"
di as txt "  responds to prices and to expenditure; the formulas of easi include it."

di ""
di as txt "{hline 78}"
