*! 11_hixdata_results.do -- section 7.7 of the note: the reference
*! specification on hixdata, its elasticity tables, and the Engel curves
*!
*! Produces
*!   out/hixdata_output.log        the Stata output quoted in the note
*!   out/elast_*.csv               e() matrices, read by make_tables.py
*!   out/engel_atmeans.pdf/.gph    Figure 1, top panel
*!   out/engel_asobserved.pdf/.gph Figure 1, bottom panel
*!
*! Run from the replication/ directory:  do 11_hixdata_results.do

clear all
set more off
set linesize 100
* ---- run this from the replication/ directory ----------------------------
* Every script locates the module (../src), the data (../examples) and the
* frozen R reference (R_reference/out) relative to the current directory:
*     cd <path-to-repository>/replication
*     do 11_hixdata_results.do
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
local OUT "`ROOT'/replication/out"

use "`ROOT'/examples/hixdata.dta", clear

log using "`OUT'/hixdata_output.log", replace text name(o)
easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,		///
     lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers)	///
     lnexpenditure(log_y) demographics(age hsex carown time tran)		///
     power(5) py pz zy snames(foodh foodr rent oper furn cloth tranop recr pers) ///
     compensated demoelast checks nolog
log close o

* ---- e() matrices to CSV, full double precision ------------------------
mata:
function wcsv(string scalar f, string scalar m) {
	real matrix M; string vector rn, cn; real scalar fh, i, j; string scalar s
	M = st_matrix(m)
	rn = st_matrixrowstripe(m)[,2]; cn = st_matrixcolstripe(m)[,2]
	fh = fopen(f, "w")
	fput(fh, "," + invtokens(cn', ","))
	for (i = 1; i <= rows(M); i++) {
		s = rn[i]
		for (j = 1; j <= cols(M); j++) s = s + "," + strofreal(M[i,j], "%18.0g")
		fput(fh, s)
	}
	fclose(fh)
}
d = "`OUT'/"
wcsv(d + "elast_exp.csv",         "e(elast_exp)")
wcsv(d + "elast_exp_se.csv",      "e(elast_exp_se)")
wcsv(d + "elast_price_nc.csv",    "e(elast_price_nc)")
wcsv(d + "elast_price_nc_se.csv", "e(elast_price_nc_se)")
wcsv(d + "elast_price_c.csv",     "e(elast_price_c)")
wcsv(d + "elast_price_c_se.csv",  "e(elast_price_c_se)")
wcsv(d + "elast_demo.csv",        "e(elast_demo)")
wcsv(d + "elast_demo_se.csv",     "e(elast_demo_se)")
end
di as txt "N = " e(N) ", iterations = " e(iter) ", time = " %6.2f e(time) " s"
di as txt "Engel residual " %9.3e e(chk_engel) ", Cournot residual " %9.3e e(chk_cournot)

* ---- Figure 1 ------------------------------------------------------------
estat engel, saving("`OUT'/engel_atmeans.gph", replace)
graph export "`OUT'/engel_atmeans.pdf", replace
estat engel, asobserved saving("`OUT'/engel_asobserved.gph", replace)
di as txt "asobserved bandwidth = " %7.4f r(bwidth)
graph export "`OUT'/engel_asobserved.pdf", replace

di as res _n "11_hixdata_results: done.  Then:  python make_tables.py out"
