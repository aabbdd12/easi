*! Table13.do -- Table 13 (Section 8.7) of the note: the running times of the
*! note.
*!
*! The reference specification on hixdata (4,847 households, 9 goods, 576
*! parameters): compat mode (the like-for-like comparison with the R
*! package), the default mode without and with the elasticity standard
*! errors, with every reporting table, with the elasticities of the
*! individuals (the whole estimation weighted by a household size), and the
*! conventional variance.  The reduced Mexican survey (2,477 households,
*! 3 goods), vce(svy): without and with the correction for the non-buyers,
*! and the cost of one bootstrap replication of the whole procedure.  Times
*! depend on the machine; the ratios are what the note relies on.
*!
*! Run from the replication/ directory:  do Table13.do

clear all
set more off
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Table13.do
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

use "`ROOT'/examples/hixdata.dta", clear
local SH sfoodh sfoodr srent soper sfurn scloth stranop srecr spers
local PR pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers
local SPEC lnprices(`PR') lnexpenditure(log_y)				///
	   demographics(age hsex carown time tran) power(5) py pz zy nolog notable

qui easi `SH', `SPEC' compat noelastse
di as txt "hixdata, compat, no elasticity SE        " as res %7.2f e(time) " s"
qui easi `SH', `SPEC' noelastse
di as txt "hixdata, default, no elasticity SE       " as res %7.2f e(time) " s"
qui easi `SH', `SPEC'
di as txt "hixdata, default, with elasticity SE     " as res %7.2f e(time) " s"
qui easi `SH', `SPEC' detail
di as txt "hixdata, default, every reporting table  " as res %7.2f e(time) " s"
qui gen byte __n = 1 + mod(_n, 5)
qui easi `SH', `SPEC' hhsize(__n)
di as txt "hixdata, individuals (a household size)  " as res %7.2f e(time) " s"
qui easi `SH', `SPEC' noelastse vce(conventional)
di as txt "hixdata, conventional variance, no SE    " as res %7.2f e(time) " s"

use "`ROOT'/examples/mex_bench.dta", clear
local M "w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) vce(svy) nolog notable"
qui easi `M'
di as txt "mex_bench, vce(svy)                      " as res %7.2f e(time) " s"
qui easi `M' selvars(age isMale)
di as txt "mex_bench, vce(svy), selection           " as res %7.2f e(time) " s"
qui easi `M' selvars(age isMale) noelastse
di as txt "mex_bench, selection, no elasticity SE   " as res %7.2f e(time) " s"
timer clear 1
timer on 1
qui easi w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) ///
	selvars(age isMale) vce(bootstrap, reps(20) seed(1) svy) nolog notable
timer off 1
qui timer list 1
di as txt "mex_bench, selection, bootstrap: " as res %5.2f r(t1) / 21 as txt " s per estimation (20 replications and the full sample)"
di as res _n "Table13: done"
