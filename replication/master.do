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
*! The bootstrap scripts (05a, 05b, 07a, 07b) run 400 replications each and
*! take the longest -- about an hour in all on a laptop.  Set BOOT to 0 to
*! skip them on a first pass.

clear all
set more off
if "$BOOT" == "" global BOOT 1		// global BOOT 0 before -do master- skips the bootstraps

local scripts							///
	01a_compat_lock_coefficients				///
	01b_compat_lock_elasticities				///
	02_effect_of_each_correction				///
	03_invariance_to_dropped_good				///
	04a_formulas_finite_differences				///
	04b_formulas_homogeneity				///
	06_generated_regressor					///
	08_analytic_jacobian					///
	09_tables_orientation_aggregation			///
	10_easidiag_checks					///
	11_hixdata_results					///
	12_easidiag_cases					///
	13_timing
if $BOOT {
	local scripts `scripts'					///
	05a_coefficient_se_bootstrap				///
	05b_elasticity_se_bootstrap				///
	07a_pweight_bootstrap					///
	07b_design_bootstrap
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
