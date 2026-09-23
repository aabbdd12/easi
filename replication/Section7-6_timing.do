*! 13_timing.do -- the running times quoted in section 7.6 of the note
*!
*! Reference specification on hixdata (4,847 households, 9 goods, 576
*! parameters): compat mode (the like-for-like comparison with the R package),
*! the default mode without and with the elasticity standard errors, and the
*! default mode with every reporting table.  Times depend on the machine;
*! the ratios are what the note relies on.
*!
*! Run from the replication/ directory:  do Section7-6_timing.do

clear all
set more off
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do Section7-6_timing.do
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
	   demographics(age hsex carown time tran) power(5) py pz zy nolog

qui easi `SH', `SPEC' compat noelastse
di as txt "compat, no elasticity SE            " as res %7.2f e(time) " s"
qui easi `SH', `SPEC' noelastse
di as txt "default, no elasticity SE           " as res %7.2f e(time) " s"
qui easi `SH', `SPEC'
di as txt "default, with elasticity SE         " as res %7.2f e(time) " s"
qui easi `SH', `SPEC' detail
di as txt "default, every reporting table      " as res %7.2f e(time) " s"
qui easi `SH', `SPEC' noelastse vce(conventional)
di as txt "conventional variance, no elast. SE " as res %7.2f e(time) " s"
di as res _n "13_timing: done"
