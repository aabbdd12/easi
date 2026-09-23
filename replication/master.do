*! master.do -- reproduce every number, table and figure of the technical note
*!
*! Run from this directory:
*!     cd replication
*!     do master.do
*!
*! Data come from ../examples (hixdata.dta, mex_bench.dta), the module from
*! ../src, the frozen R reference from R_reference/out.  Everything written
*! goes to out/.  README.md maps each script to the tables and figures.
*!
*! The bootstrap scripts (Table2, Table3, Table4, Table5) run 400 replications each and
*! take the longest -- about an hour in all on a laptop.  Set BOOT to 0 to
*! skip them on a first pass.

clear all
set more off
if "$BOOT" == "" global BOOT 1		// global BOOT 0 before -do master- skips the bootstraps

local scripts							///
	SectionA1_compat_coefficients				///
	SectionA1_compat_elasticities				///
	SectionA_corrections_one_at_a_time			///
	Table1							///
	Section4-3_orientation_aggregation			///
	Section4-4_analytic_jacobian				///
	Section5-3_generated_regressor				///
	Table6							///
	Section6_easidiag_checks				///
	Section7-6_timing					///
	Table8_9_10_Figure1					///
	Table12							///
	SectionA5_homogeneity
if $BOOT {
	local scripts `scripts'					///
	Table2							///
	Table5							///
	Table4							///
	Table3
}

capture mkdir out
timer clear 1
timer on 1
foreach s of local scripts {
	di as txt _n "{hline 78}" _n "==> `s'" _n "{hline 78}"
	do `s'.do
}
timer off 1
timer list 1
di as res _n "master.do: all scripts ran.  Tables: python make_tables.py out"
