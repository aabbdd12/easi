*! master.do -- reproduce every number, table and figure of the technical note
*!
*! Run from this directory:
*!     cd replication
*!     do master.do
*!
*! Data come from ../examples (hixdata.dta, mex_bench.dta), the module from
*! ../src, the frozen R reference from R_reference/out, the data-generating
*! process of the Monte Carlo experiments from dgp_easi.do.  Everything
*! written goes to out/.  README.md maps each script to the tables and figures.
*!
*! Three groups, by cost:
*!   - the checks and the tables without resampling: a few minutes;
*!   - the bootstraps (Tables 3 to 6, 9 and 10), 500 replications each: about
*!     forty minutes; global BOOT 0 before -do master- skips them;
*!   - the brute force and the Monte Carlo experiments (Tables 7, 8 and 9):
*!     about an hour and a half, most of it the Monte Carlo under selection;
*!     global MC 0 before -do master- skips them.

clear all
set more off
if "$BOOT" == "" global BOOT 1
if "$MC"   == "" global MC 1

local scripts							///
	SectionA1_compat_coefficients				///
	SectionA1_compat_elasticities				///
	SectionA_corrections_one_at_a_time			///
	Table1							///
	Table2							///
	Section4-3_types_checks					///
	Section4-4_orientation_aggregation			///
	Section4-5_analytic_jacobian				///
	Section5-3_generated_regressor				///
	Section6-1_pimpute					///
	Section6-4_probits					///
	Section6-4_general_engine				///
	Section6-4_selection_oracle				///
	Section6-5_bootstrap_option				///
	Section6-6_selection_elasticities			///
	Table11							///
	Section7_easidiag_checks				///
	Section7_easidiag_selection				///
	Table13							///
	Table14_15_16_Figure1					///
	Table18							///
	SectionA5_homogeneity
if $BOOT {
	local scripts `scripts'					///
	Table3							///
	Table4							///
	Table5							///
	Table6							///
	Table9_bootstrap					///
	Table10
}
if $MC {
	local scripts `scripts'					///
	Table7							///
	Table8							///
	Table9_bruteforce					///
	Table9_montecarlo
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
di as res _n "master.do: all scripts ran.  Tables 14 to 16: python make_tables.py out"
