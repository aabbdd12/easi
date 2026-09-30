*! version 2.0.1  30sep2026
*! License: GPL-3.0-or-later (see LICENSE); https://github.com/aabbdd12/easi
*! EASI demand system -- pure Stata/Mata, no R: iterated restricted 3SLS,
*! elasticities of three types, influence-function variance, the households
*! that do not buy (Shonkwiler-Yen), easidiag, predict, estat engel
*! Araar Abdelkrim

program easi, eclass
	version 14.2

	* hidden subcommand used by easi_p.ado; see -program Predict- below
	gettoken sub rest : 0
	if `"`sub'"' == "_predict" {
		Predict `rest'
		exit
	}
	if `"`sub'"' == "_engel" {
		Engel `rest'
		exit
	}
	* hidden subcommand used by easidiag.ado
	if `"`sub'"' == "_diag" {
		Diagnose `rest'
		exit
	}

	if replay() {
		if "`e(cmd)'" != "easi" error 301
		Display `0'
		exit
	}
	* vce(bootstrap): the whole procedure on resamples of the households
	local 0save `"`0'"'
	capture syntax anything [if] [in] [aweight fweight pweight iweight], VCE(string asis) [*]
	if !_rc & regexm(`"`vce'"', "^ *boot") {
		_easi_boot `0save'
		exit
	}
	local 0 `"`0save'"'
	Estimate `0'
end

*--------------------------------------------------------------------------
* vce(bootstrap [, reps(#) seed(#) strata(varname) psu(varname) svy]): the
* point estimates of the full sample, and variances from the replications of
* the whole procedure (missing prices filled, probits and system estimated
* again on each resample of households, or of PSUs within strata); the
* households keep their weights.  Replaces the variance of the coefficients,
* the standard errors of the elasticities (households and market) and, with
* selection, e(V_sel).  The design of equaids.
program _easi_boot, eclass
	version 14.2
	syntax varlist(min=2 numeric) [if] [in] [aweight fweight pweight iweight],	///
		VCE(string asis) [ STars SAVEres(string) noTABle noLOg noELASTse * ]
	local vl `varlist'
	local ifin `if' `in'
	local wt0 `weight'
	local ex0 `"`exp'"'
	local opts `"`options'"'
	local rstars `stars'
	local rsave `"`saveres'"'
	local rtable `table'
	gettoken vt vo : vce, parse(",")
	gettoken comma vo : vo, parse(",")
	local 0 `", `vo'"'
	syntax [, Reps(integer 200) SEED(string) STRata(varname) PSU(varname) SVY]
	if `reps' < 2 {
		di as error "{bf:vce(bootstrap)}: reps() must be at least 2"
		exit 198
	}
	local weight `wt0'
	local exp `"`ex0'"'
	local wgt = cond("`weight'" != "", "[`weight'`exp']", "")
	if "`svy'" != "" {
		quietly svyset
		if "`r(settings)'" == ", clear" | "`r(settings)'" == "" {
			di as error "{bf:vce(bootstrap, svy)}: the data are not svyset"
			exit 119
		}
		if "`psu'" == "" & !inlist("`r(su1)'", "", "_n") local psu `r(su1)'
		if "`strata'" == "" local strata `r(strata1)'
		if "`r(wvar)'" != "" {
			if "`weight'" != "" {
				di as error "{bf:vce(bootstrap, svy)}: the weight comes from svyset; do not specify weights"
				exit 198
			}
			local wgt "[pweight = `r(wvar)']"
		}
	}
	local bso ""
	if "`strata'" != "" local bso `bso' strata(`strata')
	if "`psu'" != "" local bso `bso' cluster(`psu')
	local design = trim(cond("`psu'" != "", "PSUs `psu'", "households") +	///
		cond("`strata'" != "", " within strata `strata'", ""))
	if "`seed'" != "" set seed `seed'

	* the full sample: the point estimates
	di as txt "(estimating the full sample)"
	quietly Estimate `vl' `ifin' `wgt', `opts' vce(robust) notable nolog
	tempvar es
	quietly gen byte `es' = e(sample)
	local sel = ("`e(selection)'" != "")
	local mkt = 0
	capture confirm matrix e(elast_exp_mkt)
	if !_rc local mkt 1
	tempname full
	_estimates hold `full', copy restore

	* the replications
	di as txt "Bootstrap replications (" as res `reps' as txt ", resampling `design')"
	mata: __ea_bs = J(0, 0, .)
	local ok 0
	forvalues r = 1/`reps' {
		preserve
		quietly keep if `es'
		quietly bsample, `bso'
		capture quietly Estimate `vl' `wgt', `opts' vce(robust) noelastse notable nolog
		* a replication that does not converge is dropped and counted
		local rc = _rc
		if !`rc' & e(converged) != 1 local rc 430
		if !`rc' {
			mata: __ea_bs = (rows(__ea_bs) ? __ea_bs \ _easi_bvec(`mkt', `sel') : _easi_bvec(`mkt', `sel'))
			local ++ok
			di as txt "." _c
		}
		else di as error "x" _c
		restore
		if mod(`r', 50) == 0 di as txt " " `r'
	}
	if mod(`reps', 50) di ""
	_estimates unhold `full'
	if `ok' < 2 {
		di as error "{bf:vce(bootstrap)}: fewer than 2 successful replications"
		exit 430
	}

	* the variances from the replications, under the names of the analytic ones
	tempname T
	mata: _easi_bpost(__ea_bs, `mkt', `sel')
	mata: mata drop __ea_bs
	matrix `T' = e(V)
	mata: st_replacematrix("`T'", st_matrix("_easi_bV"))
	ereturn repost V = `T'
	local sfxl ""
	if `mkt' local sfxl "_mkt"
	foreach sf in "" `sfxl' {
		foreach m in exp price_nc price_c demo {
			capture confirm matrix e(elast_`m'`sf'_se)
			if !_rc {
				matrix `T' = e(elast_`m'`sf'_se)
				mata: st_replacematrix("`T'", st_matrix("_easi_bse_`m'`sf'"))
				ereturn matrix elast_`m'`sf'_se = `T'
			}
		}
	}
	* the legacy names, kept consistent
	matrix `T' = e(elast_exp_se)
	ereturn matrix elast_income_se = `T'
	matrix `T' = e(elast_price_nc_se)'
	ereturn matrix elast_price_se = `T'
	if `sel' {
		matrix `T' = e(V_sel)
		mata: st_replacematrix("`T'", st_matrix("_easi_bVs"))
		ereturn matrix V_sel = `T'
	}
	foreach m in bV bVs bse_exp bse_price_nc bse_price_c bse_demo bse_exp_mkt	///
		bse_price_nc_mkt bse_price_c_mkt bse_demo_mkt {
		capture matrix drop _easi_`m'
	}
	ereturn local  vce       "bootstrap"
	ereturn local  vcetype   "Bootstrap"
	ereturn scalar N_reps    = `reps'
	ereturn scalar N_reps_ok = `ok'
	ereturn local  boot_design "`design'"
	if "`seed'" != "" ereturn local boot_seed "`seed'"
	local sv ""
	if `"`rsave'"' != "" local sv `"saveres(`rsave')"'
	Display, `rstars' `sv' `rtable'
end

*--------------------------------------------------------------------------
program Estimate, eclass
	version 14.2

	* Stata's timers are one global pool of 100 and a running timer cannot be
	* read back, so none of them can be saved and restored.  Number 100 is the
	* least likely to be in use; the help file says so.
	timer clear 100
	timer on 100

	syntax varlist(min=2 numeric) [if] [in]				///
		[aweight fweight pweight iweight] ,			///
		[ PRices(varlist numeric)				///
		  LNPRices(varlist numeric)				///
		  EXPenditure(varname numeric)				///
		  LNEXPenditure(varname numeric)			///
		  DEMOgraphics(varlist numeric)				///
		  POWer(integer 5)					///
		  PY PZ ZY						///
		  INPY(integer 0) INPZ(integer 0) INZY(integer 0)	///
		  INTERPZ(varlist numeric)				///
		  COMPat						///
		  LEGacy(string)					///
		  VCE(string)						///
		  SNames(string)					///
		  DEC(integer 4)					///
		  DISLAS(integer 1)					///
		  DREGRES(integer 0)					///
		  TOLerance(real 1e-6)					///
		  ITERate(integer 100)					///
		  noELASTse						///
		  COMPENSated						///
		  DEMOELast						///
		  CHECKS						///
		  DETail						///
		  STars							///
		  SAVEres(string)					///
		  noTABle						///
		  ELASticities(string)					///
		  HHsize(varname numeric)				///
		  PIMpute(varlist)					///
		  SELection SELGoods(namelist) SELVars(string)		///
		  noLOg ]

	* fweights are refused (since 2.0.0): they are almost always sampling
	* weights in disguise, such as fw = int(pw*10000), and a frequency weight
	* would count each household that many times
	_easi_nofw `weight'

	* -detail- is the three reporting tables at once
	if "`detail'" != "" {
		local compensated compensated
		local demoelast   demoelast
		local checks      checks
	}
	local docq = ("`compensated'" != "")

	local shares `varlist'
	local J : word count `shares'
	local neq = `J' - 1

	if ("`prices'" != "") & ("`lnprices'" != "") {
		di as error "specify either {bf:prices()} or {bf:lnprices()}, not both"
		exit 198
	}
	if ("`prices'`lnprices'" == "") {
		di as error "must specify {bf:prices()} or {bf:lnprices()}"
		exit 198
	}
	local inprice `prices'`lnprices'
	local npr : word count `inprice'
	if `npr' != `J' {
		di as error "`npr' price variable(s) for `J' budget share(s)"
		exit 198
	}

	if ("`expenditure'" != "") & ("`lnexpenditure'" != "") {
		di as error "specify either {bf:expenditure()} or {bf:lnexpenditure()}, not both"
		exit 198
	}
	if ("`expenditure'`lnexpenditure'" == "") {
		di as error "must specify {bf:expenditure()} or {bf:lnexpenditure()}"
		exit 198
	}

	local T : word count `demographics'
	if `T' == 0 {
		di as error "{bf:demographics()} is required (at least one variable)"
		exit 198
	}
	if `power' < 1 {
		di as error "{bf:power()} must be a positive integer"
		exit 198
	}

	* py/pz/zy are the documented flags; inpy()/inpz()/inzy() are the legacy
	* sr_easi spellings, accepted so that old do-files and the dialog keep working
	local pyi = ("`py'" != "") | (`inpy' != 0)
	local pzi = ("`pz'" != "") | (`inpz' != 0)
	local zyi = ("`zy'" != "") | (`inzy' != 0)
	local cmp = ("`compat'" != "")
	* the lists that name goods and variables, checked together
	_easi_chkvars, shares(`shares') prices(`prices' `lnprices')		///
		expend(`expenditure' `lnexpenditure') demo(`demographics')	///
		interpz(`interpz') snames(`snames')
	* interpz() names the demographics crossed with prices: it implies pz
	if "`interpz'" != "" & !`pzi' {
		di as txt "(interpz() implies pz: the prices are crossed with `interpz')"
		local pzi 1
	}

	* ---- legacy behaviours ------------------------------------------------
	* -compat- switches all of them on at once, which is what reproduces the R
	* package bit for bit.  -legacy()- is an undocumented developer switch that
	* turns them on one at a time, so that the effect of each fix can be
	* measured separately (see tests/test_step3.do).  The names are documented
	* in tests/README.md.
	local LGALL interpz instr quant ez lastse elast elastse sigma1 condse
	foreach f of local LGALL {
		local lg_`f' = `cmp'
	}
	if "`legacy'" != "" {
		foreach f of local LGALL {
			local lg_`f' = 0
		}
		foreach f of local legacy {
			if "`f'" == "all" {
				foreach g of local LGALL {
					local lg_`g' = 1
				}
			}
			else if "`f'" != "none" {
				local ok : list posof "`f'" in LGALL
				if `ok' == 0 {
					di as error "{bf:legacy()}: unknown item `f'"
					exit 198
				}
				local lg_`f' = 1
			}
		}
	}
	local lgon
	foreach f of local LGALL {
		if `lg_`f'' local lgon `lgon' `f'
	}
	local allleg = ("`lgon'" == "`LGALL'")

	* ---- which elasticities ------------------------------------------------
	* The mean over the households of the household elasticities, weighted
	* by the weight (households) or by the weight times the total expenditure
	* (market, the elasticities of the total demand); both come from the same
	* estimation and both are stored.  The legacy formulas reproduce the R
	* package, which has the first only.
	* individuals: each household counts for its weight times its size,
	* hhsize(), in the whole estimation; the default when hhsize() is given.
	if `"`elasticities'"' == "" & "`hhsize'" != "" local etype individuals
	else {
		_easi_etype `elasticities'
		local etype `r(etype)'
	}
	local indiv = ("`etype'" == "individuals")
	if `indiv' & "`hhsize'" == "" {
		di as error "{bf:elasticities(individuals)} requires {bf:hhsize(}{it:varname}{bf:)}, the size of the household"
		exit 198
	}
	local lgelas = `lg_elast' | `lg_ez' | `lg_elastse'
	if "`etype'" != "households" & `lgelas' {
		di as error "{bf:elasticities(`etype')} needs the corrected elasticity formulas;"
		di as error "{bf:compat} and {bf:legacy()} reproduce the R package, which has the elasticities of the households only"
		exit 198
	}
	if `indiv' & `"`elasticities'"' == "" {
		di as txt "(hhsize() given: the elasticities of the individuals, each household counting"
		di as txt " for its weight x `hhsize' in the whole estimation)"
	}
	if !`indiv' & "`hhsize'" != "" {
		di as txt "(hhsize() is used by elasticities(individuals) only: `hhsize' not used)"
	}
	local domkt = !`lgelas' & !`indiv'

	* ---- selection of the buyers (Shonkwiler and Yen 1999) -----------------
	* The goods corrected (selgoods(); default: those with zero shares, the
	* last good excepted, which closes the system) and the variables of the
	* probits only (selvars(): segments separated by ";", "vars" for every
	* corrected good or "good: vars"); the names and the rules of duvm and
	* equaids.
	if "`selgoods'" != "" | `"`selvars'"' != "" local selection selection
	local sel = ("`selection'" != "")
	local zall ""
	local lastgood : word `J' of `shares'
	if `sel' {
		if `lgelas' {
			di as error "{bf:selection} needs the corrected elasticity formulas (not {bf:compat} or {bf:legacy()})"
			exit 198
		}
		gettoken v1x : vce
		if inlist("`v1x'", "conventional", "oim") {
			di as error "{bf:selection}: use vce(robust), vce(cluster) or vce(svy); the conventional variance ignores the probits"
			exit 198
		}
		_easi_selparse, shares(`shares') selgoods(`selgoods')			///
			selvars(`"`selvars'"') xin(`shares' `prices' `lnprices'	///
			`expenditure' `lnexpenditure' `demographics')
		local selgoods `r(selgoods)'
		local zall `r(zall)'
		local zin `r(zin)'
		forvalues k = 1/`J' {
			local g : word `k' of `shares'
			local z_`g' `r(z_`k')'
		}
	}

	* ---- which demographics interact with prices ------------------------
	* corrected: those listed in interpz(), all of them by default -- the
	* documented intent.  compat: reproduce the R bug, where
	*   interpz <- ifelse(length(interpz) > 1, interpz, 1:nsoc)
	* collapses the vector to a single element (see tests/README.md).
	local ipz
	if `pzi' {
		if "`interpz'" == "" local interpz `demographics'
		foreach v of local interpz {
			local k : list posof "`v'" in demographics
			if `k' == 0 {
				di as error "{bf:interpz()}: `v' is not in {bf:demographics()}"
				exit 198
			}
			local ipz `ipz' `k'
		}
		if `lg_interpz' {
			local nn : word count `ipz'
			if `nn' > 1 local ipz : word 1 of `ipz'
			else        local ipz 1
		}
	}
	local nipz : word count `ipz'

	* ---- estimation sample ----------------------------------------------
	marksample touse
	* pimpute(): a missing price is filled further down, so it does not take
	* the household out of the sample here (a nonpositive price level does)
	local mkp = cond("`pimpute'" == "", "`inprice'", "")
	markout `touse' `shares' `mkp' `expenditure' `lnexpenditure' `demographics'
	if `sel' & "`zall'" != "" {
		qui count if `touse'
		local n0 = r(N)
		markout `touse' `zall'
		qui count if `touse'
		if r(N) < `n0' di as txt "(selvars(): `=string(`n0' - r(N), "%12.0fc")' households dropped for missing values of `zall')"
	}

	tempvar lnx
	if ("`expenditure'" != "") qui gen double `lnx' = ln(`expenditure') if `touse'
	else                       qui gen double `lnx' = `lnexpenditure'   if `touse'

	local lnplist
	foreach v of local inprice {
		tempvar lp
		if ("`prices'" != "") {
			if "`pimpute'" != "" qui replace `touse' = 0 if `v' <= 0
			qui gen double `lp' = ln(`v') if `touse'
		}
		else qui gen double `lp' = `v' if `touse'
		local lnplist `lnplist' `lp'
	}
	markout `touse' `lnx' `=cond("`pimpute'" == "", "`lnplist'", "")'
	if `indiv' {
		markout `touse' `hhsize'
		qui count if `touse' & `hhsize' <= 0
		if r(N) {
			di as error "{bf:hhsize()}: `hhsize' has nonpositive values"
			exit 411
		}
	}

	* ---- survey design -----------------------------------------------------
	* Read before the weight block, so that the svyset weight can be adopted
	* when the user did not give one.  -svyset- stores its design in the data
	* characteristics; the strata variable is often a string and the PSU
	* identifier is usually unique only within its stratum, so both are turned
	* into dense integer codes further down, once the sample is final.
	qui count if `touse'

	local svy_on 0
	if "`vce'" != "" {
		gettoken v1x vrestx : vce
		if "`v1x'" == "svy" {
			local svy_on 1
			local svyw : char _dta[_svy_wvar]
			local svyt : char _dta[_svy_wtype]
			local svyp : char _dta[_svy_su1]
			local svys : char _dta[_svy_strata1]
			if "`svyp'`svys'" == "" {
				di as error "{bf:vce(svy)}: the data are not {bf:svyset}"
				exit 119
			}
			local svyf : char _dta[_svy_fpc1]
			if "`svyp'" != "" markout `touse' `svyp', strok
			if "`svys'" != "" markout `touse' `svys', strok
			if "`svyf'" != "" markout `touse' `svyf'
			if "`weight'" == "" & "`svyw'" != "" {
				local weight pweight
				local exp "= `svyw'"
			}
			qui count if `touse'
			local N = r(N)
			if `N' == 0 error 2000
		}
	}
	local N = r(N)
	if `N' == 0 error 2000

	local noisy = ("`log'" == "")
	* the delta-method elasticity standard errors need a numerical Jacobian,
	* costing K re-evaluations.  Inside -bootstrap- they are pure waste.
	local doese = ("`elastse'" != "noelastse")

	* ---- weights ----------------------------------------------------------
	* individuals: the weight times hhsize() (an aweight of hhsize() alone
	* without a weight), for the whole estimation; e(wexp) keeps it, so that
	* estat rebuilds the same.
	if `indiv' {
		if "`weight'" != "" {
			local exp "= (`=trim(subinstr("`exp'", "=", "", 1))') * (`hhsize')"
		}
		else {
			local weight aweight
			local exp "= `hhsize'"
		}
	}
	* aw/pw are normalised to sum to N (Stata's convention); fw/iw are used as
	* they stand.  The R package has no weights at all, so nothing here can be
	* validated against it -- see tests/test_step4.do for the internal checks.
	tempvar wvar
	if "`weight'" != "" {
		qui gen double `wvar' `exp' if `touse'
		qui replace `touse' = 0 if missing(`wvar') | `wvar' <= 0
		qui count if `touse'
		local N = r(N)
		if `N' == 0 error 2000
		qui su `wvar' if `touse', meanonly
		if inlist("`weight'", "aweight", "pweight") {
			qui replace `wvar' = `wvar' * r(N) / r(sum) if `touse'
			local neff = r(N)
		}
		else local neff = r(sum)
	}
	else {
		qui gen double `wvar' = 1 if `touse'
		local neff = `N'
	}

	* ---- variance estimator ------------------------------------------------
	* Default: ROBUST.  Budget-share equations are heteroskedastic, and a 400-
	* replication bootstrap puts the conventional standard error 23% below the
	* truth on the coefficients in y -- the ones that shape the Engel curves --
	* while the robust one lands within 1% everywhere (tests/test_step8.do).
	* That gap is heteroskedasticity, not the generated regressor, which
	* test_step11.do measures at under half a percent.  -compat- (equivalently
	* -legacy(all)-) and -vce(conventional)- give the homoskedastic 3SLS
	* covariance that systemfit, reg3 and the R package report.
	local clvar
	local vtype = cond(`allleg', 0, 1)
	if "`vce'" == "" & "`weight'" == "pweight" local vce robust
	if "`vce'" != "" {
		gettoken v1 vrest : vce
		if inlist("`v1'", "conventional", "oim") local vtype 0
		else if "`v1'" == "robust" local vtype 1
		else if "`v1'" == "cluster" {
			local vtype 2
			local clvar = trim("`vrest'")
			cap confirm variable `clvar'
			if _rc {
				di as error "{bf:vce(cluster ...)}: `clvar' is not a variable"
				exit 198
			}
			markout `touse' `clvar', strok
			qui count if `touse'
			local N = r(N)
		}
		else if "`v1'" == "svy" {
			local vtype 3
		}
		else {
			di as error "{bf:vce()}: must be conventional, robust, "	///
			   "cluster {it:varname} or svy"
			exit 198
		}
	}

	* ---- missing prices: pimpute() -----------------------------------------
	* The weighted mean of the log prices of the households of the same group
	* that have one, the groups of pimpute() in order (_easi_pimpute, the rule
	* of equaids); households still without a price leave the sample.  Done
	* once the sample is otherwise final, with the weight of the estimation,
	* so that the donors are the households estimated on (estat and predict
	* redo it on e(sample)).  The analytic standard errors include the
	* imputation (each donor moves the mean of its group, hence the prices it
	* fills): the log prices before filling and the groups, as integers, go
	* to Mata (_easi_impmap).
	local lp0list ""
	local impglist ""
	if "`pimpute'" != "" {
		local k0 0
		foreach v of local lnplist {
			local ++k0
			tempvar lp0`k0'
			quietly gen double `lp0`k0'' = `v' if `touse'
			local lp0list `lp0list' `lp0`k0''
		}
	}
	if "`pimpute'" != "" {
		_easi_pimpute `lnplist', touse(`touse') wt(`wvar') groups(`pimpute') ///
			names(`shares')
		local k0 0
		foreach g of local pimpute {
			local ++k0
			tempvar ig`k0'
			quietly egen long `ig`k0'' = group(`g') if `touse'
			local impglist `impglist' `ig`k0''
		}
		local pinotes `"`r(notes)'"'
		foreach l of local pinotes {
			di as txt `"`l'"'
		}
		qui count if `touse'
		local N = r(N)
		if `N' == 0 error 2000
		* the weights of the final sample
		if inlist("`weight'", "aweight", "pweight") {
			qui su `wvar' if `touse', meanonly
			qui replace `wvar' = `wvar' * r(N) / r(sum) if `touse'
			local neff = r(N)
		}
		else if "`weight'" != "" {
			qui su `wvar' if `touse', meanonly
			local neff = r(sum)
		}
		else local neff = `N'
	}

	* ---- selection: the goods corrected, the probits ------------------------
	* One probit of buying per corrected good, on the final sample, with the
	* weight of the estimation: its Phi, phi and the influence functions of
	* its coefficients are left in temporary variables for the system.
	local selgstr ""
	local zmstr ""
	local phv ""
	local phiv ""
	local ifv ""
	if `sel' {
		_easi_selgoods, shares(`shares') touse(`touse')			///
			selgoods(`selgoods') zall(`zall') zin(`zin')
		local selgoods `r(selgoods)'
		local selgstr `r(selgstr)'
		local zmstr `r(zmstr)'
		foreach g of local selgoods {
			tempvar ph phi
			qui gen double `ph' = .
			qui gen double `phi' = .
			local phv `phv' `ph'
			local phiv `phiv' `phi'
			local nzg : word count `z_`g''
			local np = 1 + `J' + `T' + `nzg'
			forvalues c = 1/`np' {
				tempvar ifc
				qui gen double `ifc' = .
				local ifv `ifv' `ifc'
			}
		}
		tempname selalpha seldiag selas
		mata: _easi_selprobits("`shares'", "`lnplist'", "`lnx'",		///
			"`demographics'", "`zall'", "`touse'", "`wvar'", "`selgstr'",	///
			"`zmstr'", "`phv'", "`phiv'", "`ifv'", "`selalpha'", "`seldiag'", "`selas'")
		local rn _cons `inprice' `demographics' `zall'
		matrix rownames `selalpha' = `rn'
		matrix colnames `selalpha' = `: list shares - lastgood'
		matrix rownames `seldiag' = `: list shares - lastgood'
	}

	* ---- survey identifiers ------------------------------------------------
	tempvar svystr svypsu
	if `svy_on' {
		if "`svys'" != "" qui egen long `svystr' = group(`svys') if `touse'
		else              qui gen  long `svystr' = 1             if `touse'
		* a PSU label is usually recycled across strata, so the cluster
		* identifier has to be the pair
		if "`svyp'" != "" {
			qui egen long `svypsu' = group(`svys' `svyp') if `touse'
		}
		else qui gen long `svypsu' = _n if `touse'

		tempvar tagpsu tagstr
		qui bysort `touse' `svypsu': gen byte `tagpsu' = (_n == 1) & `touse'
		qui bysort `touse' `svystr': gen byte `tagstr' = (_n == 1) & `touse'
		qui count if `tagpsu'
		local n_psu = r(N)
		qui count if `tagstr'
		local n_str = r(N)
		* a stratum with one PSU contributes nothing and, left unsaid, would
		* quietly shrink the variance
		tempvar npsu
		qui bysort `touse' `svystr': egen long `npsu' = total(`tagpsu')
		qui count if `tagstr' & `npsu' < 2
		local n_single = r(N)
		if `n_single' > 0 {
			di as txt "note: `n_single' stratum(s) with a single PSU "	///
			   "contribute nothing to the variance"
		}
	}
	else {
		qui gen long `svystr' = .
		qui gen long `svypsu' = .
	}

	* ---- Mata -------------------------------------------------------------
	tempname b V Sig it crit conv
	* the general engine (every equation with its own regressors and
	* instruments) with selection; the developer scalar _easi_geng runs it
	* without selection too, where it must give the usual engine's results
	capture confirm scalar _easi_geng
	local geng = `sel' | (_rc == 0)
	if `geng' {
		mata: _easi_g_run("`shares'", "`lnplist'", "`lnx'", "`demographics'",	///
			"`touse'", "`wvar'", `neff', `vtype', "`clvar'",		///
			"`=cond(`svy_on', "`svystr'", "")'",				///
			"`=cond(`svy_on', "`svypsu'", "")'",				///
			"`=cond(`svy_on', "`svyf'", "")'",				///
			`power', `pyi', `pzi', `zyi', "`ipz'",			///
			`lg_instr', `lg_quant', `lg_ez', `lg_lastse', `lg_elast',	///
			`lg_elastse', `doese', `docq', `lg_condse',			///
			`=cond(`lg_sigma1',1,100)',					///
			`tolerance', `iterate', `noisy',				///
			"`b'", "`V'", "`Sig'", "`it'", "`crit'", "`conv'", `domkt',	///
			"`selgstr'", "`zall'", "`zmstr'", "`selas'", "`ifv'",	///
			"`lp0list'", "`impglist'", `=("`pimpute'" != "")')
	}
	else mata: _easi_run("`shares'", "`lnplist'", "`lnx'", "`demographics'",	///
		"`touse'", "`wvar'", `neff', `vtype', "`clvar'",		///
		"`=cond(`svy_on', "`svystr'", "")'",				///
		"`=cond(`svy_on', "`svypsu'", "")'",				///
		"`=cond(`svy_on', "`svyf'", "")'",				///
		`power', `pyi', `pzi', `zyi', "`ipz'",			///
		`lg_instr', `lg_quant', `lg_ez', `lg_lastse', `lg_elast',	///
		`lg_elastse', `doese', `docq', `lg_condse',			///
		`=cond(`lg_sigma1',1,100)',					///
		`tolerance', `iterate', `noisy',				///
		"`b'", "`V'", "`Sig'", "`it'", "`crit'", "`conv'", `domkt',	///
		"`lp0list'", "`impglist'", `=("`pimpute'" != "")')

	* ---- coefficient names -------------------------------------------------
	local cn
	forvalues i = 1/`neq' {
		local eqn : word `i' of `shares'
		local cn `cn' `eqn':_cons
		forvalues r = 1/`power' {
			local cn `cn' `eqn':y`r'
		}
		foreach z of local demographics {
			local cn `cn' `eqn':`z'
		}
		if `zyi' {
			foreach z of local demographics {
				local cn `cn' `eqn':y_`z'
			}
		}
		forvalues j = 1/`neq' {
			local pj : word `j' of `inprice'
			local cn `cn' `eqn':`pj'
		}
		if `pyi' {
			forvalues j = 1/`neq' {
				local pj : word `j' of `inprice'
				local cn `cn' `eqn':y_`pj'
			}
		}
		if `pzi' {
			foreach t of local ipz {
				local zt : word `t' of `demographics'
				forvalues j = 1/`neq' {
					local pj : word `j' of `inprice'
					local cn `cn' `eqn':`zt'_`pj'
				}
			}
		}
		if `sel' {
			if `: word `i' of `selgstr'' local cn `cn' `eqn':delta
		}
	}
	matrix colnames `b' = `cn'
	matrix colnames `V' = `cn'
	matrix rownames `V' = `cn'

	local K = colsof(`b')
	timer off 100
	qui timer list 100
	local etime = r(t100)
	timer clear 100

	ereturn post `b' `V', esample(`touse') obs(`N')

	ereturn local  cmd          "easi"
	ereturn local  predict      "easi_p"
	ereturn local  estat_cmd    "easi_estat"
	ereturn local  shares       "`shares'"
	ereturn local  prices       "`inprice'"
	ereturn local  expvar       "`expenditure'`lnexpenditure'"
	ereturn scalar prlog      = ("`prices'" != "")
	ereturn scalar explog     = ("`expenditure'" != "")
	ereturn local  demographics "`demographics'"
	if `pzi' {
		local ipznames
		foreach t of local ipz {
			local zt : word `t' of `demographics'
			local ipznames `ipznames' `zt'
		}
		ereturn local interpz "`ipznames'"
	}
	ereturn local  vce       = cond(`vtype'==0, "conventional",		///
		cond(`vtype'==1, "robust", cond(`vtype'==2, "cluster", "svy")))
	if `vtype' == 2 ereturn local clustvar "`clvar'"
	if `vtype' == 3 {
		ereturn local  svyfpc    "`svyf'"
		ereturn local  svyunit   "`svyp'"
		ereturn local  svystrata "`svys'"
		ereturn local  svywvar   "`svyw'"
		ereturn scalar N_psu     = `n_psu'
		ereturn scalar N_strata  = `n_str'
	}
	ereturn local  vcetype2  = cond(`lg_condse', "conditional on y", "generated regressor")
	if "`weight'" != "" {
		ereturn local wtype "`weight'"
		ereturn local wexp  "`exp'"
	}
	ereturn scalar time      = `etime'
	ereturn scalar N_eff     = `neff'
	ereturn scalar dec       = `dec'
	ereturn scalar dislas    = `dislas'
	ereturn scalar dregres   = `dregres'
	ereturn local  legacy    "`lgon'"
	ereturn local  mode      = cond(`allleg', "compat", cond("`lgon'" == "", "corrected", "mixed"))

	* ---- elasticities ------------------------------------------------------
	* short names for the elasticity tables; default = the share variables
	local gn
	forvalues i = 1/`J' {
		local gn `gn' `: word `i' of `shares''
	}
	if "`snames'" != "" {
		local nsn : word count `snames'
		if `nsn' != `J' {
			di as error "{bf:snames()}: `nsn' name(s) for `J' good(s)"
			exit 198
		}
		local gn `snames'
		ereturn local snames "`snames'"
	}
	foreach m in ei eise er {
		matrix colnames _easi_`m' = `gn'
	}
	matrix rownames _easi_ei   = elasticity
	matrix rownames _easi_eise = std_err
	matrix rownames _easi_er   = semielast
	foreach m in ep epse semip slut cq epnc epncse epc epcse {
		matrix colnames _easi_`m' = `gn'
		matrix rownames _easi_`m' = `gn'
	}
	foreach m in ez ezse {
		matrix colnames _easi_`m' = `gn'
		matrix rownames _easi_`m' = `demographics'
	}
	* Legacy names, kept bit for bit so that existing do-files keep working.
	* elast_price is [price, good]; compensated_q is eta^H + I, [good, price].
	* Neither is what the command displays -- see the elast_*_nc / _c pair.
	ereturn matrix compensated_q  = _easi_cq
	ereturn matrix slutsky        = _easi_slut
	ereturn matrix semi_price     = _easi_semip
	ereturn matrix semi_income    = _easi_er, copy
	ereturn matrix elast_demo_se  = _easi_ezse
	ereturn matrix elast_demo     = _easi_ez
	ereturn matrix elast_price_se = _easi_epse
	ereturn matrix elast_price    = _easi_ep
	ereturn matrix elast_income_se = _easi_eise, copy
	ereturn matrix elast_income   = _easi_ei, copy

	* Current names.  Every price matrix here is [good, price]: row j, column k
	* is the elasticity of good j with respect to the price of good k.
	ereturn matrix elast_price_c_se  = _easi_epcse
	ereturn matrix elast_price_c     = _easi_epc
	ereturn matrix elast_price_nc_se = _easi_epncse
	ereturn matrix elast_price_nc    = _easi_epnc
	ereturn matrix elast_exp_se      = _easi_eise
	ereturn matrix elast_exp         = _easi_ei
	ereturn matrix semi_exp          = _easi_er

	* Aggregation identities, computed on the estimation sample.  Engel:
	* sum_j w_j eta^x_j = 1.  Cournot: sum_j w_j eta^k_j = -w_k.  Both are
	* exact for the mean-based formulas, so a departure means a defect in the
	* elasticity code, not in the data.
	ereturn scalar chk_engel   = _easi_chkE
	ereturn scalar chk_cournot = _easi_chkC
	ereturn local  report      = trim("`compensated' `demoelast' `checks'")

	* The elasticities of the market, from the same estimation: the same
	* names with _mkt (e(elast_exp_mkt), e(elast_exp_mkt_se)...).
	if `domkt' {
		foreach m in ei eise {
			matrix colnames _easi_`m'_mkt = `gn'
		}
		matrix rownames _easi_ei_mkt   = elasticity
		matrix rownames _easi_eise_mkt = std_err
		foreach m in epnc epncse epc epcse {
			matrix colnames _easi_`m'_mkt = `gn'
			matrix rownames _easi_`m'_mkt = `gn'
		}
		foreach m in ez ezse {
			matrix colnames _easi_`m'_mkt = `gn'
			matrix rownames _easi_`m'_mkt = `demographics'
		}
		ereturn matrix elast_price_c_mkt_se  = _easi_epcse_mkt
		ereturn matrix elast_price_c_mkt     = _easi_epc_mkt
		ereturn matrix elast_price_nc_mkt_se = _easi_epncse_mkt
		ereturn matrix elast_price_nc_mkt    = _easi_epnc_mkt
		ereturn matrix elast_demo_mkt_se     = _easi_ezse_mkt
		ereturn matrix elast_demo_mkt        = _easi_ez_mkt
		ereturn matrix elast_exp_mkt_se      = _easi_eise_mkt
		ereturn matrix elast_exp_mkt         = _easi_ei_mkt
		ereturn scalar chk_engel_mkt   = _easi_chkE_mkt
		ereturn scalar chk_cournot_mkt = _easi_chkC_mkt
		foreach m in ep epse er semip slut cq {
			capture matrix drop _easi_`m'_mkt
		}
	}
	ereturn local elasticities "`etype'"
	if "`hhsize'" != "" ereturn local hhsize "`hhsize'"
	if "`pimpute'" != "" ereturn local pimpute "`pimpute'"
	if `sel' {
		matrix `seldiag' = `seldiag', _easi_selvif'
		matrix drop _easi_selvif
		matrix colnames `seldiag' = pct_buyers pseudo_r2 n_separated vif_delta
		ereturn scalar sel_ycomp_pct  = _easi_ycomp[1, 1]
		ereturn scalar sel_ycomp_mean = _easi_ycomp[1, 2]
		ereturn scalar sel_ycomp_sd   = _easi_ycomp[1, 3]
		matrix drop _easi_ycomp
		ereturn local  selection "shonkwiler-yen"
		ereturn local  selgoods  "`selgoods'"
		ereturn local  selvars   `"`selvars'"'
		ereturn matrix sel_alpha = `selalpha'
		ereturn matrix sel_diag  = `seldiag'
		foreach g of local selgoods {
			if "`z_`g''" != "" ereturn local sel_z_`g' "`z_`g''"
		}
		* for predict and estat: the goods corrected (one flag per equation),
		* the probit-only variables and which good uses which
		ereturn local sel_flags "`selgstr'"
		ereturn local sel_zall  "`zall'"
		ereturn local sel_zmask "`zmstr'"
		* the joint variance of the coefficients and of the probits
		local an ""
		foreach g of local selgoods {
			local an `an' sel_`g':_cons
			foreach v of local inprice {
				local an `an' sel_`g':`v'
			}
			foreach v of local demographics {
				local an `an' sel_`g':`v'
			}
			foreach v of local z_`g' {
				local an `an' sel_`g':`v'
			}
		}
		matrix colnames _easi_Vsel = `cn' `an'
		matrix rownames _easi_Vsel = `cn' `an'
		ereturn matrix V_sel = _easi_Vsel
	}

	ereturn matrix Sigma     = `Sig'
	ereturn scalar N         = `N'
	ereturn scalar neq       = `neq'
	ereturn scalar ngoods    = `J'
	ereturn scalar nsoc      = `T'
	ereturn scalar power     = `power'
	ereturn scalar py        = `pyi'
	ereturn scalar pz        = `pzi'
	ereturn scalar zy        = `zyi'
	ereturn scalar k_eq      = `K' / `neq'
	ereturn scalar iter      = `it'
	ereturn scalar crit      = `crit'
	ereturn scalar converged = `conv'

	* A non-converged fit must FAIL, not be returned with a warning: under
	* -bootstrap- a silently bad replication is worse than a lost one.  Stata
	* counts the failures and reports them.
	if !`conv' {
		di as error "convergence not achieved in `iterate' iterations"
		di as error "criterion = " %10.3e `crit' ", tolerance = " %10.3e `tolerance'
		exit 430
	}

	if `"`saveres'"' != "" local saveres `"saveres(`saveres')"'
	Display, `stars' `saveres' `table'
end

*--------------------------------------------------------------------------
* Pre-estimation diagnostic.  Reached through easidiag.ado.
program Diagnose, rclass
	version 14.2

	syntax varlist(min=2 numeric) [if] [in]				///
		[aweight fweight pweight iweight] ,			///
		[ PRices(varlist numeric)				///
		  LNPRices(varlist numeric)				///
		  EXPenditure(varname numeric)				///
		  LNEXPenditure(varname numeric)			///
		  DEMOgraphics(varlist numeric)				///
		  POWer(integer 5)					///
		  PY PZ ZY						///
		  INPY(integer 0) INPZ(integer 0) INZY(integer 0)	///
		  INTERPZ(varlist numeric)				///
		  SNames(string)					///
		  TOLerance(real 1e-6)					///
		  ELASticities(string)					///
		  HHsize(varname numeric)				///
		  PIMpute(varlist)					///
		  SELection SELGoods(namelist) SELVars(string)		///
		  * ]

	* fweights are refused (since 2.0.0): they are almost always sampling
	* weights in disguise, such as fw = int(pw*10000), and a frequency weight
	* would count each household that many times
	_easi_nofw `weight'

	local shares `varlist'
	local J : word count `shares'
	local neq = `J' - 1
	local lastgood : word `J' of `shares'

	if ("`prices'" != "") & ("`lnprices'" != "") {
		di as error "specify either {bf:prices()} or {bf:lnprices()}, not both"
		exit 198
	}
	if ("`prices'" == "") & ("`lnprices'" == "") {
		di as error "{bf:prices()} or {bf:lnprices()} required"
		exit 198
	}
	if ("`expenditure'" == "") & ("`lnexpenditure'" == "") {
		di as error "{bf:expenditure()} or {bf:lnexpenditure()} required"
		exit 198
	}
	if "`demographics'" == "" {
		di as error "{bf:demographics()} required"
		exit 198
	}

	local pyi = ("`py'" != "") | (`inpy' != 0)
	local pzi = ("`pz'" != "") | (`inpz' != 0)
	local zyi = ("`zy'" != "") | (`inzy' != 0)
	_easi_chkvars, shares(`shares') prices(`prices' `lnprices')		///
		expend(`expenditure' `lnexpenditure') demo(`demographics')	///
		interpz(`interpz') snames(`snames')
	if "`interpz'" != "" & !`pzi' {
		di as txt "(interpz() implies pz: the prices are crossed with `interpz')"
		local pzi 1
	}
	foreach v of local interpz {
		if !`: list v in demographics' {
			di as error "{bf:interpz()}: `v' is not in {bf:demographics()}"
			exit 198
		}
	}

	local ipz
	if `pzi' {
		if "`interpz'" == "" local interpz `demographics'
		foreach v of local interpz {
			local kk : list posof "`v'" in demographics
			local ipz `ipz' `kk'
		}
	}

	* the selection of the buyers: the probits the estimation would run
	if "`selgoods'" != "" | `"`selvars'"' != "" local selection selection
	local sel = ("`selection'" != "")
	local zall ""
	if `sel' {
		_easi_selparse, shares(`shares') selgoods(`selgoods')			///
			selvars(`"`selvars'"') xin(`shares' `prices' `lnprices'	///
			`expenditure' `lnexpenditure' `demographics')
		local selgoods `r(selgoods)'
		local zall `r(zall)'
		local zin `r(zin)'
		local zcnt ""
		forvalues k = 1/`neq' {
			local zcnt `zcnt' `: word count `r(z_`k')''
		}
	}

	* the weights of the estimation it prepares: elasticities(individuals),
	* the default with hhsize(), weights each household by its weight times
	* hhsize() in the whole estimation, so the diagnostic does the same
	if `"`elasticities'"' == "" & "`hhsize'" != "" local etype individuals
	else {
		_easi_etype `elasticities'
		local etype `r(etype)'
	}
	local indiv = ("`etype'" == "individuals")
	if `indiv' & "`hhsize'" == "" {
		di as error "{bf:elasticities(individuals)} requires {bf:hhsize(}{it:varname}{bf:)}, the size of the household"
		exit 198
	}
	if `indiv' {
		if "`weight'" != "" {
			local exp "= (`=trim(subinstr("`exp'", "=", "", 1))') * (`hhsize')"
		}
		else {
			local weight aweight
			local exp "= `hhsize'"
		}
		di as txt "(the weights of elasticities(individuals): `weight' `exp')"
	}
	local nipz : word count `ipz'
	local T : word count `demographics'

	local inprice `prices' `lnprices'
	marksample touse
	local mkp = cond("`pimpute'" == "", "`inprice'", "")
	markout `touse' `shares' `mkp' `expenditure' `lnexpenditure' `demographics'
	if `sel' & "`zall'" != "" {
		qui count if `touse'
		local n0 = r(N)
		markout `touse' `zall'
		qui count if `touse'
		if r(N) < `n0' di as txt "(selvars(): `=string(`n0' - r(N), "%12.0fc")' households dropped for missing values of `zall')"
	}
	if `indiv' {
		markout `touse' `hhsize'
		qui count if `touse' & `hhsize' <= 0
		if r(N) {
			di as error "{bf:hhsize()}: `hhsize' has nonpositive values"
			exit 411
		}
	}
	qui count if `touse'
	local N = r(N)
	qui count
	local Nall = r(N)
	if `N' == 0 {
		di as error "no observations"
		exit 2000
	}

	tempvar lnx wvar
	if ("`expenditure'" != "") qui gen double `lnx' = ln(`expenditure') if `touse'
	else                       qui gen double `lnx' = `lnexpenditure'   if `touse'

	local lnplist
	if ("`prices'" != "") {
		foreach v of local prices {
			tempvar lp
			if "`pimpute'" != "" qui replace `touse' = 0 if `v' <= 0
			qui gen double `lp' = ln(`v') if `touse'
			local lnplist `lnplist' `lp'
		}
	}
	else if "`pimpute'" != "" {
		* copies: the imputation writes into them
		foreach v of local lnprices {
			tempvar lp
			qui gen double `lp' = `v' if `touse'
			local lnplist `lnplist' `lp'
		}
	}
	else local lnplist `lnprices'

	if "`weight'" != "" {
		* as the estimation: the households with a missing or nonpositive
		* weight are left out (1.0.0 wrote abs(`exp'), which is abs(= w):
		* every weighted diagnostic stopped on "invalid name")
		qui gen double `wvar' `exp' if `touse'
		qui replace `touse' = 0 if missing(`wvar') | `wvar' <= 0
		qui count if `touse'
		local N = r(N)
		if `N' == 0 {
			di as error "no observations"
			exit 2000
		}
		if inlist("`weight'", "aweight", "pweight") {
			qui su `wvar' if `touse', meanonly
			qui replace `wvar' = `wvar' * `N' / r(sum) if `touse'
			local neff = `N'
		}
		else {
			qui su `wvar' if `touse', meanonly
			local neff = r(sum)
		}
	}
	else {
		qui gen double `wvar' = 1 if `touse'
		local neff = `N'
	}

	* missing prices, as the estimation would fill them (pimpute())
	if "`pimpute'" != "" {
		_easi_pimpute `lnplist', touse(`touse') wt(`wvar') groups(`pimpute') ///
			names(`shares')
		local pinotes `"`r(notes)'"'
		foreach l of local pinotes {
			di as txt `"`l'"'
		}
		qui count if `touse'
		local N = r(N)
		if `N' == 0 {
			di as error "no observations"
			exit 2000
		}
		if inlist("`weight'", "aweight", "pweight") {
			qui su `wvar' if `touse', meanonly
			qui replace `wvar' = `wvar' * `N' / r(sum) if `touse'
			local neff = `N'
		}
		else if "`weight'" != "" {
			qui su `wvar' if `touse', meanonly
			local neff = r(sum)
		}
		else local neff = `N'
	}

	* short names for the goods
	local gn
	forvalues i = 1/`J' {
		local v : word `i' of `shares'
		local gn `gn' `v'
	}
	if "`snames'" != "" {
		local nsn : word count `snames'
		if `nsn' == `J' local gn `snames'
	}

	* names of the design columns, in the order _easi_design() builds them
	local dn _cons
	forvalues r = 1/`power' {
		local dn `dn' y`r'
	}
	local dn `dn' `demographics'
	if `zyi' {
		foreach v of local demographics {
			local dn `dn' `v'#y
		}
	}
	forvalues i = 1/`neq' {
		local g : word `i' of `gn'
		local dn `dn' np_`g'
	}
	if `pyi' {
		forvalues i = 1/`neq' {
			local g : word `i' of `gn'
			local dn `dn' np_`g'#y
		}
	}
	if `pzi' {
		foreach t of local ipz {
			local zt : word `t' of `demographics'
			forvalues i = 1/`neq' {
				local g : word `i' of `gn'
				local dn `dn' np_`g'#`zt'
			}
		}
	}

	* The dialog emits one superset of options for both commands, and a future
	* version may use more of them, so the catch-all stays.  But an option that
	* is received and dropped without a word is the failure mode this module
	* exists to avoid, so say what was ignored.
	if `"`options'"' != "" {
		di as txt "note: not used by {bf:easidiag}, ignored: "		///
		   as res `"`options'"'
	}

	mata: _easi_diagm("`shares'", "`lnplist'", "`lnx'", "`demographics'",	///
		"`touse'", "`wvar'", `neff', `power', `pyi', `pzi', `zyi',	///
		"`ipz'", `tolerance')

	* the probits of the estimation (same sample, same weights), and the
	* identification of each delta at the starting point
	if `sel' {
		_easi_selgoods, shares(`shares') touse(`touse')			///
			selgoods(`selgoods') zall(`zall') zin(`zin')
		local selgoods `r(selgoods)'
		local selgstr `r(selgstr)'
		local zmstr `r(zmstr)'
		tempname selalpha seldiag selas
		mata: _easi_selprobits("`shares'", "`lnplist'", "`lnx'",		///
			"`demographics'", "`zall'", "`touse'", "`wvar'", "`selgstr'",	///
			"`zmstr'", "", "", "", "`selalpha'", "`seldiag'", "`selas'")
		mata: _easi_d_selvif("`shares'", "`lnplist'", "`lnx'",		///
			"`demographics'", "`zall'", "`touse'", "`wvar'", "`selgstr'",	///
			"`zmstr'", "`selas'", `power', `pyi', `pzi', `zyi', "`ipz'")
		matrix `seldiag' = `seldiag', _d_selvif'
		matrix drop _d_selvif
		matrix colnames `seldiag' = pct_buyers pseudo_r2 n_separated vif_delta
		matrix rownames `seldiag' = `: list shares - lastgood'
	}

	*---------------------------------------------------------------- display
	* nblock counts what makes the design unfittable -- zero variance, a
	* near-dependency, ill-conditioning, rank deficiency.  Those, and only
	* those, make the fixed-point section meaningless, because it needs a
	* design that can be solved.  Few price regimes or a weak first stage are
	* worth saying but do not stop the estimator.
	local nprob 0
	local nblock 0
	local advice

	di ""
	di as txt "{hline 78}"
	di as txt "EASI model diagnostic" _col(50) as txt			///
	   "Observations   = " as res %10.0fc `N'
	di as txt "Goods          = " as res %-10.0g `J' _col(50) as txt	///
	   "Coef. per eq.  = " as res %10.0g _d_k
	di as txt "y power        = " as res %-10.0g `power' _col(50) as txt	///
	   "Coefficients   = " as res %10.0g _d_K
	local ia = cond(`pyi',"py ","") + cond(`pzi',"pz ","") + cond(`zyi',"zy ","")
	if "`ia'" == "" local ia "none"
	di as txt "Interactions   = " as res %-10s "`ia'" _col(50) as txt	///
	   "Instruments/eq = " as res %10.0g _d_nz
	di as txt "{hline 78}"

	*---- 1
	di ""
	di as txt "1. Data admissibility"
	local t1 = cond(_d_sum1 < 1e-7, "ok", "FAIL")
	di as txt "   budget shares sum to 1, max |sum-1|" _col(56)		///
	   as res %11.3e _d_sum1 _col(70) as res "`t1'"
	if _d_sum1 >= 1e-7 {
		local ++nprob
		local advice `"`advice' "the budget shares do not sum to 1" "'
		local ++nblock
	}
	di as txt "   share values below 0" _col(56) as res %11.0f _d_neg	///
	   _col(70) as res cond(_d_neg == 0, "ok", "note")
	di as txt "   share values above 1" _col(56) as res %11.0f _d_gt1	///
	   _col(70) as res cond(_d_gt1 == 0, "ok", "note")
	di as txt "   observations dropped as missing" _col(56)			///
	   as res %11.0f `=`Nall' - `N'' _col(70)				///
	   as res cond(`Nall' == `N', "ok", "note")

	*---- 2  depistage des colonnes du design, SANS estimation
	di ""
	di as txt "2. Screening of the design columns"
	di as txt "   (evaluated at the starting point y = Stone index; no estimation)"
	local kk = colsof(_d_scr)
	local kk = rowsof(_d_scr)
	local nflat 0
	local nvif  0
	di as txt "   column" _col(32) "sd" _col(45) "VIF"
	* the constant is not a regressor to screen, and its VIF is undefined.
	* Guard every comparison with missing(): in Stata, . > 30 is TRUE.
	forvalues i = 2/`kk' {
		local nm : word `i' of `dn'
		local tag
		if _d_scr[`i',2] < 1e-10 {
			local tag "  no variation"
			local ++nflat
		}
		else if !missing(_d_scr[`i',3]) & _d_scr[`i',3] > 30 {
			local tag "  collinear"
			local ++nvif
		}
		else if missing(_d_scr[`i',3]) {
			local tag "  collinear (VIF undefined)"
			local ++nvif
		}
		local show = ("`tag'" != "")
		if !missing(_d_scr[`i',3]) & _d_scr[`i',3] > 10 local show 1
		if `show' {
			local vshow : di %11.1f _d_scr[`i',3]
			if missing(_d_scr[`i',3]) local vshow "          ."
			di as txt "   `nm'" _col(26) as res %11.4f _d_scr[`i',2]	///
			   _col(40) "`vshow'" as txt "`tag'"
		}
	}
	if `nflat' + `nvif' == 0 {
		di as txt "   (no column with zero variance or VIF above 30)"
	}
	if `nflat' > 0 {
		local ++nprob
		local advice `"`advice' "`nflat' design column(s) have no variation" "'
		local ++nblock
	}
	* a high VIF is a warning: the design still solves.  What stops the
	* fixed point being measured is an exact dependency, a flat column, a
	* rank deficiency or cond(X'X) beyond 1e6, each counted where found.
	if `nvif' > 0 {
		local ++nprob
		local vifadv 1
	}

	* near-dependencies, named.  A rank test would miss the dangerous case:
	* perturbing an exact collinearity by 1e-9 leaves the matrix full rank
	* while multiplying the coefficients by 2800 and the standard errors by
	* 3e+08.  The condition index measures the continuum, so it fires there.
	local nci 0
	forvalues c = 1/`kk' {
		if _d_ci[`c',1] > 30 | missing(_d_ci[`c',1]) {
			local grp
			local ng 0
			forvalues i = 1/`kk' {
				if _d_vp[`c',`i'] > 0.5 {
					local nm : word `i' of `dn'
					local grp `grp' `nm'
					local ++ng
				}
			}
			if `ng' >= 2 {
				local ++nci
				local ci : di %9.0f _d_ci[`c',1]
				local ci = trim("`ci'")
				* 1e8 is the floor imposed on a numerically zero
				* eigenvalue, so reaching it means exact, not 1e8
				if _d_ci[`c',1] >= 9.9e7 local ci "effectively infinite"
				di ""
				di as txt "   near-dependency, condition index "	///
				   as res "`ci'" as txt " :"
				di as txt "      " as res "`grp'"
				* an exact dependency cannot be solved through;
				* a finite condition index is a warning, and the
				* estimator will still fit (hixdata's reference
				* specification carries an index of 52 and
				* converges in four iterations)
				if _d_ci[`c',1] >= 9.9e7 local nblockci 1
				* the constant plus a set of BINARY demographics
				* is the classic exhaustive set of category
				* dummies; the constant plus y, y^2, y^3 is not,
				* it is uncentred data (section 4)
				if strpos(" `grp' ", " _cons ") {
					local alldum 1
					foreach nm of local grp {
						if "`nm'" == "_cons" continue
						local isz : list nm in demographics
						if !`isz' {
							local alldum 0
							continue
						}
						qui count if `touse' & `nm' != 0 & `nm' != 1
						if r(N) > 0 local alldum 0
					}
					if `alldum' {
						di as txt "      (this group sums to a "	///
						   "constant: an exhaustive set of dummies)"
						local catadv 1
					}
				}
			}
		}
	}
	if `nci' > 0 {
		local ++nprob
		if "`nblockci'" != "" local ++nblock
		local advice `"`advice' "`nci' near-dependency(ies) among the design columns, named above" "'
	}
	* a high VIF and a rank deficiency are two more views of the same
	* dependency; listing them again would turn one problem into four
	if "`vifadv'" != "" & `nci' == 0 {
		local advice `"`advice' "`nvif' design column(s) have VIF above 30" "'
	}
	if "`catadv'" != "" {
		local advice `"`advice' "drop one category from the exhaustive dummy set" "'
	}

	*---- 3  les variables categorielles, cause classique et silencieuse
	*
	* Trois pieges, tous invisibles dans une sortie ordinaire.
	*
	*  - un jeu de muettes EXHAUSTIF somme a la constante.  L'estimateur ne
	*    plante pas: invsym renvoie un inverse generalise qui epingle une
	*    colonne a zero, les elasticites restent justes, mais un coefficient
	*    sort avec un ecart-type comme s'il avait ete estime.  Et si la somme
	*    ne vaut 1 qu'a 1e-9 pres -- parts derivees stockees en float -- rien
	*    n'est epingle et les ecarts-types explosent d'un facteur 3e+08.
	*
	*  - une modalite rare: ses interactions prix x z sous -pz- sont non nulles
	*    pour une poignee de menages et n'apportent que du bruit.
	*
	*  - une ordinale entree telle quelle impose un effet cardinal: passer de
	*    la modalite 1 a la 2 aurait le meme effet que de la 10 a la 11.
	di ""
	di as txt "3. Categorical and discrete variables"
	di as txt "   variable" _col(26) "levels" _col(38) "storage"		///
	   _col(50) "smallest cell"
	local nrare 0
	local ncard 0
	local nfloat 0
	foreach v of local demographics {
		* -levelsof- packs every distinct value into a macro, which blows
		* past the macro length limit on a continuous variable: 15,079
		* distinct values were enough to stop the command.  Counting in
		* Mata has no limit and leaves the sort order alone.
		mata: st_numscalar("_d_nl",					///
			rows(uniqrows(st_data(., "`v'", "`touse'"))))
		local nl = _d_nl
		local ty : type `v'
		local note
		local cell "---"

		capture assert `v' == int(`v') if `touse'
		local isint = (_rc == 0)

		if `nl' == 2 {
			qui su `v' if `touse', meanonly
			local mn = r(mean)
			qui su `v' if `touse'
			local lo = r(min)
			local sh = (`mn' - `lo') / (r(max) - `lo')
			local sh = min(`sh', 1 - `sh')
			local cell : di %9.4f `sh'
			local cell = trim("`cell'")
			if `sh' < 0.02 {
				local note "rare, merge it"
				local ++nrare
			}
			}
		else if `nl' <= 12 & `isint' {
			local note "entered linearly (cardinal)"
			local ++ncard
		}
		if "`ty'" == "float" & !`isint' {
			local note "`note'  float, use double"
			local ++nfloat
		}
		di as txt "   `v'" _col(26) as res %5.0f `nl' _col(37)		///
		   as res %-10s "`ty'" _col(50) as res %13s "`cell'"		///
		   as txt "   `note'"
	}

	* An exhaustive dummy set is named by the condition index in section 2,
	* which handles any subset and not just the case where every binary
	* demographic happens to belong to it, so there is nothing to add here.
	if `nrare' > 0 {
		local ++nprob
		local advice `"`advice' "`nrare' binary demographic(s) have a cell below 2%; merge categories" "'
	}
	if `ncard' > 0 {
		local advice `"`advice' "`ncard' variable(s) look categorical but enter linearly; use dummies if the effect is not cardinal" "'
		local ++nprob
	}
	if `nfloat' > 0 {
		local ++nprob
		local advice `"`advice' "`nfloat' non-integer regressor(s) stored as float; float carries 6e-08 relative error, use double" "'
	}
	if `nrare' + `ncard' + `nfloat' == 0 {
		di as txt "   (nothing to flag)"
	}

	*---- 4
	di ""
	di as txt "4. Normalisation and conditioning"
	di as txt "   variable" _col(30) "mean" _col(43) "sd" _col(56) "|mean|/sd"
	local nunc 0
	local vn "lnexpenditure"
	forvalues i = 1/`=`J'+1' {
		if `i' > 1 {
			local g : word `=`i'-1' of `gn'
			local vn "price of `g'"
		}
		local flat = cond(_d_var[`i',2] < 1e-12, "  flat", "")
		if _d_var[`i',3] > 2 & _d_var[`i',2] > 1e-12 {
			local flat "  not centred"
			local ++nunc
		}
		di as txt "   `vn'" _col(26) as res %10.4f _d_var[`i',1]		///
		   _col(38) %10.4f _d_var[`i',2] _col(53) %10.2f _d_var[`i',3]	///
		   as txt "`flat'"
	}
	* A variable whose mean sits more than two standard deviations from zero
	* is not centred.  This is the check that matters, and it is far more
	* legible than a condition number: EASI is written for log prices and log
	* expenditure measured around a base period, and fed raw levels the columns
	* 1, y, y^2, y^3 of the design become collinear to four decimal places.
	di as txt "   variables not centred (|mean|/sd > 2)" _col(56)		///
	   as res %11.0f `nunc' _col(70)					///
	   as res cond(`nunc' == 0, "ok", "CENTRE")
	if `nunc' > 0 {
		local ++nprob
		local advice `"`advice' "`nunc' variable(s) are not centred; subtract their means" "'
	}
	di as txt "   cond(X'X) as given" _col(56) as res %11.3e _d_condR	///
	   _col(70) as res cond(_d_condR < 1e6, "ok", "ILL-COND")
	di as txt "   cond(X'X) if centred" _col(56) as res %11.3e _d_condC
	if _d_condR >= 1e6 {
		local ++nprob
		local gain = round(_d_condR / _d_condC)
		local ++nblock
		if `nci' > 0 {
			local advice `"`advice' "ill-conditioned design; the cause is the near-dependency named above, not the scaling" "'
		}
		else {
			local advice `"`advice' "ill-conditioned design; centring divides cond(X'X) by `gain'" "'
		}
	}

	*---- 3b  causes structurelles propres au modele
	di ""
	di as txt "5. Structural sources of a flat objective"
	di as txt "   polynomial basis 1,y,...,y^R  cond(Gram) by power:"
	local pline "     "
	forvalues r = 1/`power' {
		local v : di %9.2e _d_poly[`r',1]
		local pline "`pline'  R=`r' `=trim("`v'")'"
	}
	di as txt "`pline'"
	local pbad = (_d_poly[`power',1] > 1e6)
	di as txt "   cond at the chosen power(`power')" _col(56)		///
	   as res %11.3e _d_poly[`power',1] _col(70)				///
	   as res cond(`pbad', "TOO HIGH", "ok")
	if `pbad' {
		local ++nprob
		local advice `"`advice' "the y polynomial is ill-conditioned at power(`power'); centre y or lower the power" "'
	}
	di as txt "   distinct price vectors" _col(50) as res %6.0f _d_nprice	///
	   as txt " / " as res %-6.0f `N' _col(70)				///
	   as res cond(_d_nprice > 2 * _d_npar, "ok", "FEW")
	di as txt "   price parameters in the symmetric A block" _col(56)	///
	   as res %11.0f _d_npar
	if _d_nprice <= 2 * _d_npar {
		local ++nprob
		local advice `"`advice' "only `=_d_nprice' distinct price vectors for `=_d_npar' price parameters" "'
	}
	* Rank below half the goods, or a single component carrying more than 90%,
	* means prices essentially move together and A is identified in only a few
	* directions.  Rank 6 of 8 with a 66% first component, as on the Canadian
	* data, is concentrated but not degenerate.
	local rkbad = (_d_rkp < ceil(`neq'/2)) | (_d_pc1 > 0.9)
	di as txt "   effective rank of var(np)" _col(50) as res %6.0f _d_rkp	///
	   as txt " / " as res %-6.0f `neq' _col(70)				///
	   as res cond(`rkbad', "LOW", "ok")
	if `rkbad' {
		local ++nprob
		local advice `"`advice' "prices move almost proportionally; the price block is weakly identified" "'
	}
	di as txt "   share of the first price component" _col(56)		///
	   as res %11.4f _d_pc1
	di as txt "   first-stage R2 of y on the instruments" _col(56)		///
	   as res %11.4f _d_fs _col(70)						///
	   as res cond(_d_fs > 0.5, "ok", "WEAK")
	if _d_fs <= 0.5 {
		local ++nprob
		local advice `"`advice' "the first stage for y is weak (R2 = `=string(_d_fs,"%5.3f")')" "'
	}
	di as txt "   cond(var(estimated shares))" _col(56)			///
	   as res %11.3e _d_conds _col(70)					///
	   as res cond(_d_conds < 1e4, "ok", "HIGH")

	*---- 5
	di ""
	di as txt "6. Scale and variability"
	di as txt "   good" _col(26) "mean share" _col(40) "sd share"		///
	   _col(54) "sd of norm. price"
	forvalues i = 1/`J' {
		local g : word `i' of `gn'
		local sm = cond(_d_good[`i',1] < 0.01, "  small", "")
		if `i' == `J' {
			* the last good is the numeraire: np_J = p_J - p_J = 0
			di as txt "   `g'" _col(26) as res %10.4f _d_good[`i',1]	///
			   _col(40) %10.4f _d_good[`i',2] _col(58)		///
			   as txt "(base)" as txt "`sm'"
		}
		else {
			di as txt "   `g'" _col(26) as res %10.4f _d_good[`i',1]	///
			   _col(40) %10.4f _d_good[`i',2] _col(56)		///
			   as res %10.4f _d_good[`i',3] as txt "`sm'"
		}
	}
	di as txt "   smallest mean budget share" _col(56) as res %11.4f _d_wmin	///
	   _col(70) as res cond(_d_wmin >= 0.01, "ok", "SMALL")
	if _d_wmin < 0.01 {
		local ++nprob
		local advice `"`advice' "a mean budget share below 0.01 scales its semi-elasticities by 1/w" "'
	}
	di as txt "   smallest sd of a normalised log price" _col(56)		///
	   as res %11.4f _d_npmin _col(70)					///
	   as res cond(_d_npmin > 1e-6, "ok", "FLAT")
	if _d_npmin <= 1e-6 {
		local ++nprob
		local advice `"`advice' "a normalised log price does not vary; its price effects are unidentified" "'
	}

	*---- 6  le point fixe, en dernier et court-circuite
	di ""
	di as txt "7. Convergence of the fixed point"
	if `nblock' > 0 {
		di as txt "   skipped: the screening above found `nblock' problem(s) "	///
		   "that stop the design being fitted."
		di as txt "   Fix those first; this step needs a design that "	///
		   "can be solved."
	}
	else if _d_fail {
		di as error "   the first restricted GLS step returned missing:"
		di as error "   the normal matrix is singular.  See section 2."
		local ++nprob
	}
	else {
		di as txt "   Phi = 1 - p'Bp/2" _col(30) "min " as res %8.4f _d_phimin ///
		   as txt "   mean " as res %8.4f _d_phimean				///
		   as txt "   max " as res %8.4f _d_phimax
		di as txt "   households with Phi <= 0" _col(56)			///
		   as res %11.0f _d_phineg _col(70)					///
		   as res cond(_d_phineg == 0, "ok", "OUT OF RANGE")
		if _d_phineg > 0 {
			local ++nprob
			local advice `"`advice' "Phi = 1 - p'Bp/2 is non-positive somewhere; the cost function is irregular there" "'
		}
		di as txt "   first step |y1 - y0|" _col(56) as res %11.3e _d_step1
		di as txt "   contraction factor |y2-y1|/|y1-y0|" _col(56)		///
		   as res %11.4f _d_ratio _col(70)					///
		   as res cond(_d_ratio < 1, "ok", "NO CONTRACTION")
		if _d_ratio >= 1 {
			local ++nprob
			local advice `"`advice' "the iteration on y is not contracting; it may fail to converge" "'
		}
		else {
			local tolfmt : di %7.1e `tolerance'
			di as txt "   iterations predicted to reach `=trim("`tolfmt'")'" ///
			   _col(56) as res %11.0f _d_iter
		}
	}

	*---- 7
	di ""
	di as txt "8. Identification"
	di as txt "   rank of Z'WZ" _col(50) as res %6.0f _d_rz as txt " / "	///
	   as res %-6.0f _d_nz _col(70)						///
	   as res cond(_d_rz == _d_nz, "ok", "DEFICIENT")
	di as txt "   rank of Xhat'W Xhat" _col(50) as res %6.0f _d_rx		///
	   as txt " / " as res %-6.0f _d_k _col(70)				///
	   as res cond(_d_rx == _d_k, "ok", "DEFICIENT")
	if (_d_rz != _d_nz) | (_d_rx != _d_k) {
		local ++nprob
		local ++nblock
		if `nci' == 0 {
			local advice `"`advice' "rank deficient design or instruments; a regressor is collinear" "'
		}
	}
	di as txt "   instruments vs coefficients per equation" _col(50)	///
	   as res %6.0f _d_nz as txt " / " as res %-6.0f _d_k _col(70)		///
	   as res cond(_d_nz >= _d_k, "ok", "UNDER")

	*---- 9  the selection of the buyers
	* delta, the coefficient of phi_i in E[w_i] = Phi_i f_i + delta_i phi_i,
	* is identified through the variables of the probit only (selvars()) or,
	* failing them, through the nonlinearity of the normal distribution
	* alone: phi_i is then close to a combination of the Phi_i x columns.
	* VIF(delta) = [(X'WX)^-1]_dd (phi'W phi), X = (Phi_i x(y), phi_i),
	* measures it here at the starting point y = Stone index, before any
	* estimation (the estimation reports it again at the converged y).
	if `sel' {
		di ""
		di as txt "9. Selection of the buyers (Shonkwiler-Yen probits)"
		di as txt "   (the probits the estimation runs; VIF(delta) at y = Stone index)"
		di as txt "   good" _col(16) "buyers %" _col(26) "pseudo-R2"		///
		   _col(37) "separated" _col(48) "probit-only" _col(60) "VIF(delta)"
		local wweak ""
		local wsep ""
		forvalues k = 1/`neq' {
			if !`: word `k' of `selgstr'' continue
			local g : word `k' of `gn'
			local v : word `k' of `shares'
			local nzk : word `k' of `zcnt'
			local vif = `seldiag'[`k', 4]
			local tag ""
			if missing(`vif') | `vif' > 10 {
				local tag "  WEAK"
				local wweak `wweak' `v'
			}
			if `seldiag'[`k', 3] > 0 {
				local tag "`tag' SEP"
				local wsep `wsep' `v'
			}
			local vshow : di %10.1f `vif'
			if missing(`vif') local vshow "         ."
			di as txt "   `g'" _col(16) as res %8.1f `seldiag'[`k', 1]	///
			   _col(26) %9.4f `seldiag'[`k', 2] _col(37) %9.0f		///
			   `seldiag'[`k', 3] _col(48) %11.0f `nzk' _col(60)		///
			   "`vshow'" as txt "`tag'"
		}
		if "`wweak'" != "" {
			local ++nprob
			local advice `"`advice' "delta weakly identified for `wweak' (VIF(delta) above 10); add a variable of the probit only in selvars(), an exclusion restriction, or leave the good uncorrected with selgoods()" "'
		}
		if "`wsep'" != "" {
			local ++nprob
			local advice `"`advice' "the probit of `wsep' predicts households with probability 0 or 1 (separation)" "'
		}
		if "`wweak'`wsep'" == "" di as txt "   (nothing to flag)"
	}

	*---- verdict
	di ""
	di as txt "{hline 78}"
	local nadv : list sizeof advice
	if `nprob' == 0 & `nadv' == 0 {
		di as txt "Verdict: nothing to flag.  The specification should estimate cleanly"
		di as txt "         in about " as res _d_iter as txt " iterations."
	}
	else {
		di as txt "Verdict: `nadv' point(s) to look at."
		local i 0
		foreach a of local advice {
			local ++i
			di as txt "  `i'. `a'"
		}
	}
	di as txt "{hline 78}"

	return scalar N        = `N'
	return scalar nprob    = `nprob'
	return scalar cond_raw = _d_condR
	return scalar cond_ctr = _d_condC
	return scalar ratio    = _d_ratio
	return scalar iter     = _d_iter
	return scalar phimin   = _d_phimin
	return scalar wmin     = _d_wmin
	return matrix var      = _d_var
	return matrix good     = _d_good
	if `sel' {
		return local  selgoods "`selgoods'"
		return matrix sel_diag = `seldiag'
	}
end

*--------------------------------------------------------------------------
* The selection options of easi and easidiag, parsed once for both: the goods
* to correct (selgoods()) and the variables of the probits only (selvars():
* segments separated by ";", "vars" for every corrected good or "good:
* vars"); the names and the rules of duvm and equaids.  Returns the goods
* asked for, the probit-only variables r(zall), those of the k-th good in
* r(z_k) and the same as flags, r(zin): J-1 rows of 0/1 over r(zall).
program _easi_selparse, rclass
	syntax , SHares(varlist) [ SELGoods(string) SELVars(string) XIN(string) ]
	local J : word count `shares'
	local neq = `J' - 1
	local lastgood : word `J' of `shares'
	local zall ""
	foreach g of local selgoods {
		if !`: list g in shares' {
			di as error "{bf:selgoods()}: `g' is not one of the shares (`shares')"
			exit 198
		}
		if "`g'" == "`lastgood'" {
			di as error "{bf:selgoods()}: `g' is the last good, which closes the system and is not corrected"
			exit 198
		}
	}
	local selgoods : list uniq selgoods
	local zcom ""
	foreach g of local shares {
		local z_`g' ""
	}
	local rest `"`selvars'"'
	while `"`rest'"' != "" {
		gettoken seg rest : rest, parse(";")
		local seg = trim(`"`seg'"')
		if `"`seg'"' == ";" | `"`seg'"' == "" continue
		if strpos(`"`seg'"', ":") {
			local gg = trim(substr(`"`seg'"', 1, strpos(`"`seg'"', ":") - 1))
			local vv = trim(substr(`"`seg'"', strpos(`"`seg'"', ":") + 1, .))
			if `"`gg'"' == "" {
				di as error `"{bf:selvars()}: a good name is missing before ":" in "`seg'""'
				exit 198
			}
			if `: word count `gg'' > 1 {
				di as error `"{bf:selvars()}: one good per segment, "`gg'" is not one good; write "good1: vars ; good2: vars""'
				exit 198
			}
			if strpos(`"`vv'"', ":") {
				di as error `"{bf:selvars()}: one colon per segment in "`seg'"; separate the goods by ";""'
				exit 198
			}
			if `"`vv'"' == "" {
				di as error "{bf:selvars()}: no variable after `gg':"
				exit 198
			}
			if !`: list gg in shares' {
				di as error "{bf:selvars()}: `gg' is not one of the shares (`shares')"
				exit 198
			}
			if "`selgoods'" != "" & !`: list gg in selgoods' {
				di as error "{bf:selvars()}: `gg' is not corrected; add it to selgoods() (`selgoods')"
				exit 198
			}
			capture unab vv : `vv'
			if _rc {
				di as error "{bf:selvars()}, `gg': variable not found in `vv'"
				exit 111
			}
			confirm numeric variable `vv'
			local z_`gg' `z_`gg'' `vv'
		}
		else {
			capture unab vv : `seg'
			if _rc {
				di as error `"{bf:selvars()}: variable not found in "`seg'""'
				exit 111
			}
			confirm numeric variable `vv'
			local zcom `zcom' `vv'
		}
	}
	* a probit-only variable is an exclusion restriction: not in the model
	foreach g of local shares {
		local z_`g' : list zcom | z_`g'
		local zall : list zall | z_`g'
	}
	local bad : list zall & xin
	if "`bad'" != "" {
		di as error "{bf:selvars()}: `bad' already in the model (share, price, expenditure or demographic);"
		di as error "a variable of the probit only must be excluded from the share equations"
		exit 198
	}
	local zin ""
	forvalues k = 1/`neq' {
		local g : word `k' of `shares'
		foreach v of local zall {
			local inz : list v in z_`g'
			local zin `zin' `inz'
		}
	}
	return local selgoods "`selgoods'"
	return local zall "`zall'"
	return local zin "`zin'"
	forvalues k = 1/`J' {
		local g : word `k' of `shares'
		return local z_`k' "`z_`g''"
	}
end

*--------------------------------------------------------------------------
* The goods corrected, on the final sample: those asked for, or by default
* those with zero shares; the last good, which closes the system, must be
* bought by everybody.  Returns them in the order of the shares, a flag per
* equation r(selgstr) and the mask r(zmstr) of the probit-only variables of
* each equation.
program _easi_selgoods, rclass
	syntax , SHares(varlist) TOUSE(varname) [ SELGoods(string) ZALL(string) ZIN(string) ]
	local J : word count `shares'
	local neq = `J' - 1
	local lastgood : word `J' of `shares'
	local nz : word count `zall'
	local selgstr ""
	local zmstr ""
	qui count if `touse' & `lastgood' <= 0
	if r(N) {
		di as error "{bf:selection}: the last good (`lastgood') has zero shares; place last a good that every household buys"
		exit 198
	}
	if "`selgoods'" == "" {
		forvalues k = 1/`neq' {
			local g : word `k' of `shares'
			qui count if `touse' & `g' <= 0
			if r(N) local selgoods `selgoods' `g'
		}
		if "`selgoods'" == "" {
			di as error "{bf:selection}: no good has zero shares; nothing to correct"
			exit 459
		}
	}
	else {
		local sg ""
		foreach g of local selgoods {
			qui count if `touse' & `g' <= 0
			if r(N) == 0 di as txt "(selection: every household buys `g'; it is left uncorrected)"
			else local sg `sg' `g'
		}
		if "`sg'" == "" {
			di as error "{bf:selection}: every household buys each of the goods to correct; nothing to correct"
			exit 459
		}
		local selgoods `sg'
	}
	* in the order of the shares, as the parameters delta
	local sgo ""
	forvalues k = 1/`neq' {
		local g : word `k' of `shares'
		local on : list g in selgoods
		local selgstr `selgstr' `on'
		if `on' local sgo `sgo' `g'
		forvalues j = 1/`nz' {
			local inz : word `=(`k' - 1) * `nz' + `j'' of `zin'
			local zmstr `zmstr' `=`on' & `inz''
		}
	}
	local selgoods `sgo'
	return local selgoods "`selgoods'"
	return local selgstr "`selgstr'"
	return local zmstr "`zmstr'"
end

*--------------------------------------------------------------------------
* The selection of the buyers: one probit per corrected good, and the
* coefficient delta of its density term in the expected share
program _easi_seltable
	tempname G
	matrix `G' = e(sel_diag)
	local nm "`e(snames)'"
	if "`nm'" == "" local nm "`e(shares)'"
	local sh "`e(shares)'"
	di ""
	di as txt "Selection of the buyers (Shonkwiler and Yen 1999)"
	di as txt "  E[w_i] = Phi_i f_i + delta_i phi_i, one probit per corrected good"
	di as txt "{hline 13}{c TT}{hline 64}"
	di as txt _col(14) "{c |}" _col(17) "Buyers %" _col(28) "Pseudo R2" _col(39) "Perfect"	///
		_col(48) "VIF(delta)" _col(61) "delta" _col(70) "z"
	di as txt "{hline 13}{c +}{hline 64}"
	local wvif ""
	local wper 0
	forvalues i = 1/`=rowsof(`G')' {
		if missing(`G'[`i', 4]) continue
		local g : word `i' of `nm'
		local q : word `i' of `sh'
		local d = _b[`q':delta]
		local z = `d' / _se[`q':delta]
		di as txt %12s abbrev("`g'", 12) " {c |}" as res _col(16) %9.1f `G'[`i', 1]	///
			_col(28) %9.3f `G'[`i', 2] _col(39) %7.0f `G'[`i', 3] _col(48) %10.1f `G'[`i', 4] ///
			_col(59) %8.4f `d' _col(68) %6.2f `z'
		if `G'[`i', 4] > 10 local wvif "`wvif' `g'"
		if `G'[`i', 3] > 0 local wper 1
	}
	di as txt "{hline 13}{c BT}{hline 64}"
	local zz ""
	foreach q of local sh {
		if "`e(sel_z_`q')'" != "" local zz "`zz' `q': `e(sel_z_`q')';"
	}
	if "`zz'" != "" di as txt "Probit-only variables:" as res "`zz'"
	else {
		di as txt "Probits on ln p - ln x and the demographics, with no variable of the"
		di as txt "probit only (no exclusion restriction)"
	}
	di as txt "VIF(delta): variance inflation of delta due to the collinearity of phi with"
	di as txt "the other regressors of its equation; above 10, delta rests on the curvature"
	di as txt "of the probit"
	if !missing(e(sel_ycomp_pct)) {
		di as txt "y completed for the non-buyers (the households with a zero share of a"
		di as txt "corrected good): " as res %4.1f e(sel_ycomp_pct) as txt "% of the households,"	///
		   as txt " mean |change| " as res %6.4f e(sel_ycomp_mean) as txt ", "	///
		   as res %5.3f e(sel_ycomp_mean) / e(sel_ycomp_sd) as txt " sd of y"
	}
	if "`wvif'" != "" {
		* under vce(bootstrap) the remedy is already applied: a note, not a
		* warning
		if "`e(vce)'" == "bootstrap" {
			di as txt "note:`wvif': delta is weakly identified (VIF(delta) above 10); the"
			di as txt "  bootstrap standard errors take it into account"
		}
		else {
			di as err "warning:`wvif': delta is weakly identified; leave the good uncorrected"
			di as err "  with selgoods(), add a variable of the probit only in selvars(), or use"
			di as err "  vce(bootstrap)"
		}
	}
	if `wper' di as err "warning: some probits predict purchase perfectly for some households"
end

*--------------------------------------------------------------------------
* elasticities(): households (the default), individuals or market; any
* abbreviation, and mkt
program _easi_etype, rclass
	local t = lower(trim(`"`0'"'))
	if "`t'" == "" local t households
	else if strpos("households", "`t'") == 1 local t households
	else if strpos("individuals", "`t'") == 1 local t individuals
	else if strpos("market", "`t'") == 1 | "`t'" == "mkt" local t market
	else {
		di as error `"{bf:elasticities(`0')}: households, individuals or market"'
		exit 198
	}
	return local etype "`t'"
end

*--------------------------------------------------------------------------
* the lists that name goods and variables: each variable once, no variable in
* two roles, one short name per good and no name twice
program _easi_chkvars
	syntax , shares(string) prices(string) expend(string) demo(string)	///
		[ interpz(string) snames(string) ]
	foreach L in shares prices demo interpz snames {
		local d : list dups `L'
		if "`d'" != "" {
			local o = cond("`L'" == "shares", "the budget shares",		///
				cond("`L'" == "demo", "{bf:demographics()}", "{bf:`L'()}"))
			di as error "`o': `d' given more than once"
			exit 198
		}
	}
	foreach a in shares prices expend demo {
		foreach b in shares prices expend demo {
			if "`a'" >= "`b'" continue
			local both : list `a' & `b'
			if "`both'" != "" {
				di as error "`both': in two roles (`a' and `b'); a variable is either a share, a price, the expenditure or a demographic"
				exit 198
			}
		}
	}
end

*--------------------------------------------------------------------------
* fweights: refused with the way out
program _easi_nofw
	if "`0'" == "fweight" {
		di as error "fweights are not allowed: sampling weights go in {bf:[pweight=]} (or {bf:[aweight=]});"
		di as error "a frequency weight built from a sampling weight, such as int(pw*10000), is a pweight"
		exit 101
	}
end

*--------------------------------------------------------------------------
program Display
	syntax [, DEC(integer -1) DISLAS(integer -1) DREGRES(integer -1)		///
		  COMPENSated DEMOELast CHECKS DETail STars SAVEres(string)	///
		  noTABle ELASticities(string) ]

	* which elasticities: asked for now, else those of the estimation (an
	* estimation of easi 1.0.0 has the households only)
	local etype "`e(elasticities)'"
	if "`etype'" == "" local etype households
	local est "`etype'"
	if `"`elasticities'"' != "" {
		_easi_etype `elasticities'
		local etype `r(etype)'
	}
	if ("`etype'" == "individuals") != ("`est'" == "individuals") {
		if "`est'" == "individuals" {
			di as error "the estimation is for individuals (weight x `e(hhsize)'): only {bf:elasticities(individuals)} can be shown"
		}
		else di as error "{bf:elasticities(individuals)} is an estimation, with {bf:hhsize()}: run {bf:easi} again"
		exit 198
	}
	local sfx = cond("`etype'" == "market", "_mkt", "")
	if "`sfx'" != "" {
		capture confirm matrix e(elast_exp_mkt)
		if _rc {
			di as error "the market elasticities are not in e(): they need the corrected formulas (not {bf:compat} or {bf:legacy()})"
			exit 198
		}
	}

	if `dec'     == -1 local dec     = e(dec)
	if `dislas'  == -1 local dislas  = e(dislas)
	if `dregres' == -1 local dregres = e(dregres)

	* what to show: what was asked for at estimation, plus anything asked for
	* now.  -detail- is the three at once.
	local rep "`e(report)'"
	if "`detail'" != "" local rep "compensated demoelast checks"
	foreach o in compensated demoelast checks {
		if "``o''" != "" local rep "`rep' `o'"
	}
	local shcq : list posof "compensated" in rep
	local shdz : list posof "demoelast"   in rep
	local shck : list posof "checks"      in rep

	local ia = cond(e(py),"py ","") + cond(e(pz),"pz ","") + cond(e(zy),"zy ","")
	if "`ia'" == "" local ia "none"
	local J   = e(ngoods)
	local shw = cond(`dislas', `J', `J' - 1)
	local f   = "%12.`dec'f"

	di ""
	di as txt "{hline 78}"
	di as txt "EASI demand system"  _col(50) as txt "Number of obs  = " as res %10.0fc e(N)
	di as txt "Goods          = " as res %-10.0g e(ngoods)			///
	   _col(50) as txt "Coef. per eq.  = " as res %10.0g e(k_eq)
	di as txt "y power        = " as res %-10.0g e(power)			///
	   _col(50) as txt "Coefficients   = " as res %10.0g colsof(e(b))
	di as txt "Interactions   = " as res %-10s "`ia'"			///
	   _col(50) as txt "Mode           = " as res %10s "`=e(mode)'"
	if e(pz) di as txt "  prices x     : " as res "`=e(interpz)'"
	if "`e(wtype)'" != "" {
		di as txt "Weights        = " as res "`e(wtype)' `e(wexp)'"
	}
	if "`etype'" == "households"  local etxt "each household counts for its weight"
	if "`etype'" == "market"      local etxt "the total demand: weight x total expenditure"
	if "`etype'" == "individuals" local etxt "each household counts for its weight x its size"
	di as txt "Elasticities   = " as res "`etype'" as txt " (`etxt')"
	if "`e(pimpute)'" != "" {
		di as txt "Missing prices = " as res "filled, pimpute(`e(pimpute)')"
	}
	local cvtxt = cond("`e(clustvar)'" != "", " (`e(clustvar)')", "")
	local gvtxt = cond("`e(vcetype2)'" == "conditional on y", ", conditional on y", "")
	di as txt "Variance       = " as res "`=e(vce)'`cvtxt'`gvtxt'"
	if e(converged) {
		di as txt "Converged in " as res e(iter)				///
		   as txt " iteration(s), criterion = " as res %10.3e e(crit)	///
		   _col(50) as txt "  Execution time = " as res			///
		   %8.2f e(time) as txt " s"
	}
	else di as error "warning: did not converge (criterion = " %10.3e e(crit) ")"
	di as txt "{hline 78}"

	if "`e(selection)'" != "" & "`table'" == "" _easi_seltable

	*---------------------------------------------------- Table 01: coefficients
	if `dregres' & "`table'" == "" {
		di ""
		di as txt "Table 01: Estimated coefficients"
		ereturn display
	}

	*------------------------------------------------- Tables 02 to 05
	* through the private copy of tabstars: on the screen, in a file
	* (saveres(): .docx .tex .xlsx .csv .md), or both
	if "`table'" == "" | `"`saveres'"' != "" {
		local fq ""
		if `"`saveres'"' != "" local fq `""`saveres'""'
		_easi_tables `fq', dec(`dec') shw(`shw') shcq(`shcq') shdz(`shdz')	///
			`stars' screen(`="`table'" == ""') sfx(`sfx') etype(`etype')
	}

	*------------------------------------------------ aggregation identities
	* Always computed; shown on request, but a failure speaks up regardless --
	* it would mean the elasticity code is wrong on these data, which is the
	* user's problem to know about whether or not they asked.
	local tolagg 1e-6
	local bad = (e(chk_engel`sfx') > `tolagg') | (e(chk_cournot`sfx') > `tolagg')
	if (`shck' & "`table'" == "") | `bad' {
		di ""
		di as txt "Aggregation identities (exact by construction)"
		di as txt "  Engel   |sum_j w_j eta^x_j - 1|" _col(46)		///
		   as res %11.3e e(chk_engel`sfx')
		di as txt "  Cournot max_k |sum_j w_j eta^k_j + w_k|" _col(46)	///
		   as res %11.3e e(chk_cournot`sfx')
	}
	if `bad' {
		* Under -compat- / -legacy(elast)- the failure is the point: the R
		* package's elasticity formulas differentiate the shares holding the
		* shares fixed, which breaks Cournot aggregation by about 3e-04 on
		* the reference data.  The check finds it on its own, without finite
		* differences and on any dataset -- so say what it means rather than
		* asking for a bug report about a defect we reproduce deliberately.
		local leg "`e(legacy)'"
		local lgel : list posof "elast" in leg
		if `lgel' {
			di as txt "  The legacy elasticity formulas do not satisfy " ///
			   "Cournot aggregation."
			di as txt "  This is expected under {bf:compat} and is one of " ///
			   "the defects it reproduces."
		}
		else {
			di as error "  warning: an aggregation identity fails by " ///
			   "more than `tolagg'."
			di as error "  The reported elasticities should not be " ///
			   "trusted; please report this."
		}
	}

	if `"`saveres'"' != "" di as txt "Tables written to " as res `"`saveres'"'
	if "`table'" != "" exit
	if `dislas' == 0 {
		di as txt "  (last good omitted: {bf:dislas(1)} to show it)"
	}
	capture assert missing(e(elast_price_nc_se)[1,1])
	if _rc {
		di ""
		if "`e(legacy)'" == "" {
			if "`e(vce)'" == "bootstrap" {
				di as txt "Standard errors: bootstrap of the whole procedure (" as res e(N_reps_ok)	///
				   as txt " replications, resampling `e(boot_design)')."
				if e(N_reps_ok) < e(N_reps) {
					di as txt "(" as res e(N_reps) - e(N_reps_ok) as txt	///
					   " replication(s) not converged or failed, dropped)"
				}
			}
			else {
			di as txt "Standard errors: the influence function of each reported mean, the households"
			di as txt "sampled and the coefficients estimated (" as res "`e(vce)'"	///
			   as txt " covariance" cond("`e(vcetype2)'" == "conditional on y", "", ", generated regressor included") ")."
			}
		}
		else {
			di as txt "Standard errors are the delta method applied to the "	///
			   "quantity reported, from"
			di as txt "the " as res "`e(vce)'" as txt " coefficient covariance"	///
			   cond("`e(vcetype2)'" == "conditional on y", "", " with the generated-regressor term")
			di as txt "included."
		}
		di as txt "{bf:noelastse} skips them."
	}
	if !`shcq' {
		di as txt "  {bf:compensated} adds the Hicksian table, " ///
		   "{bf:demoelast} the demographic one, {bf:detail} both."
	}
end

*--------------------------------------------------------------------------
* Tables 02 to 05 through _easi_tabstars, the private copy of tabstars shipped
* with easi (dev/make_easi_tabstars.py): the same numbers as e(), the standard
* errors in a companion table (02-b, 03-b...), stars on request (with t on
* the design degrees of freedom under vce(svy)); the price tables in one block
* however many goods (nowrap), as -table- did on Stata 17 and later.
program _easi_tables
	syntax [anything(name=fname)] , DEC(integer) SHW(integer) SHCQ(integer)	///
		SHDZ(integer) [STars SCREEN(integer 1) SFX(string) ETYPE(string)]
	local file `fname'
	local show = cond(`screen', "", "nodisplay")
	local st   = cond("`stars'" != "", "", "nostars")
	local df   = cond("`e(vce)'" == "svy", e(N_psu) - e(N_strata), 0)
	tempname T S
	if "`file'" != "" _easi_tabstars export begin using "`file'", replace
	capture noisily {
		* Table 02: expenditure elasticities, one row per good
		matrix `T' = e(elast_exp`sfx')
		matrix `T' = (`T'[1, 1..`shw'])'
		matrix colnames `T' = Elasticity
		local tt "Expenditure elasticities (`etype'), d ln q / d ln x"
		matrix `S' = e(elast_exp`sfx'_se)
		if !missing(`S'[1,1]) {
			matrix `S' = (`S'[1, 1..`shw'])'
			matrix colnames `S' = Elasticity
			_easi_tab, num(02) title("`tt'") est(`T') se(`S') dec(`dec')	///
				rowtitle("Good") df(`df') `st' `show'
		}
		else _easi_tab, num(02) title("`tt'") est(`T') dec(`dec')		///
			rowtitle("Good") `show'
		* Tables 03 (uncompensated) and 04 (compensated): rows goods,
		* columns prices
		foreach w in nc c {
			if "`w'" == "c" & !`shcq' continue
			if "`w'" == "nc" {
				local num 03
				local tt "Uncompensated (Marshallian) price elasticities (`etype'), d ln q / d ln p (rows: goods, columns: prices)"
			}
			else {
				local num 04
				local tt "Compensated (Hicksian) price elasticities (`etype'), d ln q / d ln p at constant utility (rows: goods, columns: prices)"
			}
			matrix `T' = e(elast_price_`w'`sfx')
			matrix `T' = `T'[1..`shw', 1..`shw']
			matrix `S' = e(elast_price_`w'`sfx'_se)
			if !missing(`S'[1,1]) {
				matrix `S' = `S'[1..`shw', 1..`shw']
				_easi_tab, num(`num') title("`tt'") est(`T') se(`S')	///
					dec(`dec') rowtitle("Good") df(`df') `st' `show' nowrap
			}
			else _easi_tab, num(`num') title("`tt'") est(`T') dec(`dec')	///
				rowtitle("Good") `show' nowrap
		}
		* Table 05: demographic elasticities, rows demographics, columns goods
		if `shdz' {
			local tt "Demographic elasticities (`etype'), d w / d z, share semi-elasticities (rows: demographics, columns: goods)"
			matrix `T' = e(elast_demo`sfx')
			matrix `T' = `T'[1..., 1..`shw']
			matrix `S' = e(elast_demo`sfx'_se)
			if !missing(`S'[1,1]) {
				matrix `S' = `S'[1..., 1..`shw']
				_easi_tab, num(05) title("`tt'") est(`T') se(`S')		///
					dec(`dec') rowtitle("Demographic") df(`df') `st'	///
					`show' nowrap
			}
			else _easi_tab, num(05) title("`tt'") est(`T') dec(`dec')	///
				rowtitle("Demographic") `show' nowrap
		}
	}
	if _rc {
		local rc = _rc
		capture _easi_tabstars export clear
		if "`file'" != "" di as err "saveres(): the tables could not be written to `file'"
		exit `rc'
	}
	if "`file'" != "" _easi_tabstars export end
end

* one table (and its standard errors) through the private tabstars
program _easi_tab
	syntax, num(string) title(string) est(name) [se(name) dec(integer 4)	///
		noSTars rowtitle(string) NODISplay NOWrap DF(integer 0)]
	if "`nodisplay'" == "" di ""
	if "`se'" != "" _easi_tabstars `est', se(`se') dec(`dec') title(`"`title'"')	///
		tabnumber(`num') rowtitle(`"`rowtitle'"') df(`df') `nodisplay'	///
		`stars' `nowrap'
	else _easi_tabstars `est', dec(`dec') title(`"`title'"') tabnumber(`num')	///
		rowtitle(`"`rowtitle'"') `nodisplay' `nowrap'
end

*--------------------------------------------------------------------------
* Engel curves.  Lives here, beside the Mata routines, for the same reason as
* -program Predict-.  Reached through -estat engel-.
program Engel, rclass
	version 14.2

	syntax [if] [in] , [ ATMeans ASOBserved N(integer 100) BWidth(real 0)	///
			     TRIM(real 1) Level(cilevel) noCI			///
			     DATA(string asis) SAVing(string asis) noDRAW * ]

	* Two different objects, and only the second one involves smoothing.
	*
	*  atmeans (default) -- the Engel curve proper: the model's fitted share as
	*    a function of total expenditure, with the demographics and the prices
	*    held at their weighted means.  It is the exact function evaluated on a
	*    grid, so n() is a RESOLUTION, not a smoothing parameter, and a
	*    bandwidth would be meaningless here.
	*
	*  asobserved -- the sample profile: the fitted shares as they actually vary
	*    with expenditure, demographics and prices left as observed.  That is a
	*    scatter and does need smoothing, so bwidth() applies here.  We use a
	*    local linear smoother rather than bins, because bins are a rectangular
	*    kernel with arbitrary breakpoints.  The default is the rule of thumb
	*    for local polynomial REGRESSION, which -lpoly- computes; Silverman's
	*    rule is a density rule and is not the right one for a smoother.
	if "`atmeans'" != "" & "`asobserved'" != "" {
		di as error "specify either {bf:atmeans} or {bf:asobserved}, not both"
		exit 198
	}
	local asobs = ("`asobserved'" != "")
	if `n' < 5 {
		di as error "{bf:n()} must be at least 5"
		exit 198
	}
	if !`asobs' & `n' > c(matsize) {
		di as error "{bf:n()} may not exceed matsize (" c(matsize) ")"
		exit 198
	}
	if `bwidth' < 0 {
		di as error "{bf:bwidth()} must be positive"
		exit 198
	}
	if `bwidth' > 0 & !`asobs' {
		di as error "{bf:bwidth()} applies only with {bf:asobserved}: the curve at the means is exact and is not smoothed"
		exit 198
	}
	if `trim' < 0 | `trim' >= 50 {
		di as error "{bf:trim()} must be in [0, 50)"
		exit 198
	}

	tempvar touse
	mark `touse' `if' `in'
	local shv `e(shares)'
	local prv `e(prices)'
	local dev `e(demographics)'
	local exv `e(expvar)'
	local J = e(ngoods)
	local pim "`e(pimpute)'"
	markout `touse' `shv' `=cond("`pim'" == "", "`prv'", "")' `exv' `dev'

	tempvar lnx es
	qui gen byte `es' = e(sample)
	if e(explog) qui gen double `lnx' = ln(`exv') if `touse'
	else         qui gen double `lnx' = `exv'     if `touse'
	local lnplist
	foreach v of local prv {
		tempvar lp
		if e(prlog) qui gen double `lp' = ln(`v') if `touse' | (`es' & "`pim'" != "")
		else        qui gen double `lp' = `v'      if `touse' | (`es' & "`pim'" != "")
		local lnplist `lnplist' `lp'
	}
	_easi_repimpute `lnplist', es(`es') touse(`touse')
	markout `touse' `lnx' `lnplist'

	tempvar wt
	if "`e(wtype)'" != "" qui gen double `wt' `e(wexp)' if `touse'
	else                  qui gen double `wt' = 1       if `touse'
	qui replace `touse' = 0 if missing(`wt') | `wt' <= 0

	local ipzn `e(interpz)'
	local ipz
	foreach v of local ipzn {
		local k : list posof "`v'" in dev
		local ipz `ipz' `k'
	}
	local lgq = (strpos("`e(legacy)'", "quant") > 0)

	local nm `e(snames)'
	if "`nm'" == "" local nm `shv'
	local zc = invnormal(1 - (100 - `level') / 200)

	*------------------------------------------------------------------
	* with selection: the expected shares E[w], their standard errors with
	* the variance of the probits (e(V_sel))
	local sel = ("`e(selection)'" != "")
	if `sel' markout `touse' `e(sel_zall)'
	if !`asobs' & `sel' {
		tempname EG
		mata: _easi_g_engelgrid("`shv'", "`lnplist'", "`lnx'", "`dev'",	///
			"`e(sel_zall)'", "`touse'", "`wt'", `n', `trim', `e(power)',	///
			`e(py)', `e(pz)', `e(zy)', "`ipz'", `lgq', "`EG'",		///
			"`e(sel_flags)'", "`e(sel_zmask)'")
		local what "at the mean of the demographics and of prices"
	}
	else if !`asobs' {
		tempname EG
		mata: _easi_engelgrid("`shv'", "`lnplist'", "`lnx'", "`dev'",	///
			"`touse'", "`wt'", `n', `trim', `e(power)', `e(py)',	///
			`e(pz)', `e(zy)', "`ipz'", `lgq', "`EG'")
		local what "at the mean of the demographics and of prices"
	}
	else {
		local wv
		local sv
		forvalues j = 1/`J' {
			tempvar w`j' s`j'
			qui gen double `w`j'' = . if `touse'
			qui gen double `s`j'' = . if `touse'
			local wv `wv' `w`j''
			local sv `sv' `s`j''
		}
		if `sel' {
			mata: _easi_g_engel("`shv'", "`lnplist'", "`lnx'", "`dev'",	///
				"`e(sel_zall)'", "`touse'", "`wv'", "`sv'", `e(power)',	///
				`e(py)', `e(pz)', `e(zy)', "`ipz'", `lgq',		///
				"`e(sel_flags)'", "`e(sel_zmask)'")
		}
		else mata: _easi_engel("`shv'", "`lnplist'", "`lnx'", "`dev'",	///
			"`touse'", "`wv'", "`sv'", `e(power)', `e(py)',		///
			`e(pz)', `e(zy)', "`ipz'", `lgq')
	}

	*------------------------------------------------------------------
	preserve
	if !`asobs' {
		qui drop _all
		qui svmat double `EG', names(_eg)
		qui rename _eg1 pctile
		qui rename _eg2 lnexp
		qui rename _eg3 y
		forvalues j = 1/`J' {
			local c1 = 3 + `j'
			local c2 = 3 + `J' + `j'
			qui rename _eg`c1' _w`j'
			qui rename _eg`c2' _se`j'
		}
	}
	else {
		qui keep if `touse'
		forvalues j = 1/`J' {
			qui rename `w`j'' _raw`j'
			qui rename `s`j'' _rawse`j'
		}

		* evaluation grid: trimmed percentiles of total expenditure
		tempvar gx gp
		qui gen double `gx' = .
		qui gen double `gp' = .
		forvalues i = 1/`n' {
			local q = `trim' + (100 - 2 * `trim') * (`i' - 0.5) / `n'
			qui _pctile `lnx' [aw=`wt'], percentiles(`q')
			qui replace `gx' = r(r1) in `i'
			qui replace `gp' = `q'   in `i'
		}

		* one bandwidth for every panel, so that the panels are comparable
		if `bwidth' == 0 {
			qui lpoly _raw1 `lnx' [aw=`wt'], degree(1) nograph
			local bwidth = r(bwidth)
		}
		forvalues j = 1/`J' {
			qui lpoly _raw`j' `lnx' [aw=`wt'], degree(1)		///
				bwidth(`bwidth') at(`gx') nograph generate(_w`j')
			qui lpoly _rawse`j' `lnx' [aw=`wt'], degree(1)		///
				bwidth(`bwidth') at(`gx') nograph generate(_se`j')
		}
		qui keep in 1/`n'
		qui gen double pctile = `gp'
		qui rename `gx' lnexp
		local bws = string(`bwidth', "%6.4f")
		local what "as observed, local linear, bandwidth `bws'"
	}

	forvalues j = 1/`J' {
		qui gen double _lo`j' = _w`j' - `zc' * _se`j'
		qui gen double _hi`j' = _w`j' + `zc' * _se`j'
	}
	label variable pctile "Percentiles of total expenditure"
	keep pctile lnexp _w* _se* _lo* _hi*
	order pctile lnexp

	if `"`data'"' != "" {
		qui save `data', replace
		di as txt `"curve data saved to {bf:`data'}"'
	}

	if "`draw'" != "nodraw" {
		local plots
		forvalues j = 1/`J' {
			local ttl : word `j' of `nm'
			local band
			if "`ci'" != "noci" {
				local band (rarea _lo`j' _hi`j' pctile,		///
					color(navy%25) lwidth(none))
			}
			tempname g`j'
			twoway `band'						///
			       (line _w`j' pctile, lcolor(navy) lwidth(medthick)), ///
				title("`ttl'", size(medsmall))			///
				ytitle("Budget share", size(vsmall))		///
				xtitle("") ylabel(, labsize(vsmall) angle(0))	///
				xlabel(0(20)100, labsize(vsmall))		///
				legend(off) graphregion(color(white))		///
				name(`g`j'', replace) nodraw
			local plots `plots' `g`j''
		}
		* two lines, so that the note fits a 2x2 canvas as well as a 4x4
		local nt `""Fitted budget shares, `what'""'
		local nt2
		if `trim' > 0       local nt2 "tails trimmed at `trim'%"
		if "`ci'" != "noci" {
			if "`nt2'" != "" local nt2 "`nt2'; "
			local nt2 "`nt2'`level'% confidence band"
		}
		if "`nt2'" != "" local nt `"`nt' "`nt2'""'
		* a near-square grid whatever the number of goods: 3 -> 2x2,
		* 4 -> 2x2, 9 -> 3x3, 12 -> 4x3, 16 -> 4x4; the canvas grows
		* with the grid so that the panels keep their size.  A cols()
		* or rows() given by the user wins.
		local ncol = ceil(sqrt(`J'))
		local nrow = ceil(`J' / `ncol')
		local grid
		if !strpos(`"`options'"', "cols(") & !strpos(`"`options'"', "rows(") {
			local grid cols(`ncol')
		}
		local gsize
		if !strpos(`"`options'"', "xsize(") & !strpos(`"`options'"', "ysize(") {
			local gsize xsize(`=min(2.3 * `ncol', 12)') ysize(`=min(1.9 * `nrow' + 1, 12)')
		}
		graph combine `plots', `grid' `gsize'				///
			title("Engel curves, EASI estimation")			///
			b1title("Percentiles of total expenditure", size(small)) ///
			note(`nt', size(vsmall))				///
			graphregion(color(white)) `options'
		if `"`saving'"' != "" graph save `saving'
		graph drop `plots'
	}
	restore

	return scalar n = `n'
	if `asobs' return scalar bwidth = `bwidth'
end

*--------------------------------------------------------------------------
* -predict- lives here, not in easi_p.ado, because the Mata routines below
* are private to this ado-file and cannot be reached from another one.
* easi_p.ado is a two-line relay that calls back into -easi _predict-.
program Predict
	version 14.2

	syntax [anything(name=vlist)] [if] [in] , [ SHares Y RESiduals LATent ]

	local nstat = ("`shares'" != "") + ("`y'" != "") + ("`residuals'" != "")
	if `nstat' > 1 {
		di as error "specify only one of {bf:shares}, {bf:y} or {bf:residuals}"
		exit 198
	}
	if `nstat' == 0 local shares shares
	local sel = ("`e(selection)'" != "")
	if "`latent'" != "" {
		if !`sel' {
			di as error "{bf:latent} applies after an estimation with {bf:selection}"
			exit 198
		}
		if "`y'`residuals'" != "" {
			di as error "{bf:latent} goes with {bf:shares} only"
			exit 198
		}
	}

	local what = cond("`y'" != "", 2, cond("`residuals'" != "", 3, 1))
	local nv   = cond(`what' == 2, 1, e(ngoods))

	if `"`vlist'"' == "" {
		di as error "new variable name required"
		exit 100
	}
	if index(`"`vlist'"', "*") {
		_stubstar2names `vlist', nvars(`nv')
		local varlist `s(varlist)'
		local typlist `s(typlist)'
		confirm new variable `varlist'
	}
	else {
		local 0 `"`vlist'"'
		syntax newvarlist(min=`nv' max=`nv')
	}

	tempvar touse
	mark `touse' `if' `in'
	local shv  `e(shares)'
	local prv  `e(prices)'
	local dev  `e(demographics)'
	local exv  `e(expvar)'
	local pim "`e(pimpute)'"
	markout `touse' `=cond("`pim'" == "", "`prv'", "")' `exv' `dev' `e(sel_zall)'
	if `what' != 2 markout `touse' `shv'
	* y itself needs the observed shares (it is a Stone-type index)
	markout `touse' `shv'

	tempvar lnx
	if e(explog) qui gen double `lnx' = ln(`exv') if `touse'
	else         qui gen double `lnx' = `exv'     if `touse'

	local lnplist
	tempvar es
	qui gen byte `es' = e(sample)
	foreach v of local prv {
		tempvar lp
		if e(prlog) qui gen double `lp' = ln(`v') if `touse' | (`es' & "`pim'" != "")
		else        qui gen double `lp' = `v'      if `touse' | (`es' & "`pim'" != "")
		local lnplist `lnplist' `lp'
	}
	* the missing prices of the estimation sample, filled as at estimation;
	* outside it a missing price gives a missing prediction
	_easi_repimpute `lnplist', es(`es') touse(`touse') keepout
	markout `touse' `lnx' `lnplist'

	* positions of the demographics that interact with prices
	local ipzn `e(interpz)'
	local ipz
	foreach v of local ipzn {
		local k : list posof "`v'" in dev
		local ipz `ipz' `k'
	}
	local lgq = (strpos("`e(legacy)'", "quant") > 0)

	local i = 0
	foreach v of local varlist {
		local ++i
		local ty : word `i' of `typlist'
		if "`ty'" == "" local ty double
		qui gen `ty' `v' = . if `touse'
	}

	if `sel' {
		mata: _easi_g_pred("`shv'", "`lnplist'", "`lnx'", "`dev'",	///
			"`e(sel_zall)'", "`touse'", "`varlist'", `e(power)', `e(py)',	///
			`e(pz)', `e(zy)', "`ipz'", `lgq', `what', "`latent'" != "",	///
			"`e(sel_flags)'", "`e(sel_zmask)'")
	}
	else mata: _easi_pred("`shv'", "`lnplist'", "`lnx'", "`dev'", "`touse'",	///
		"`varlist'", `e(power)', `e(py)', `e(pz)', `e(zy)', "`ipz'",	///
		`lgq', `what')

	* labels
	local i = 0
	foreach v of local varlist {
		local ++i
		if `what' == 2 label variable `v' "Implicit utility y"
		else {
			local g : word `i' of `shv'
			local lab = cond(`what' == 1, "Fitted share, `g'",	///
					 "Residual, `g'")
			if `sel' & `what' == 1 {
				local lab = cond("`latent'" != "", "Latent share, `g'",	///
					"Expected share E[w], `g'")
			}
			if `sel' & `what' == 3 local lab "Residual w - E[w], `g'"
			label variable `v' "`lab'"
		}
	}
end

*--------------------------------------------------------------------------
* After an estimation with pimpute(): fill the missing log prices (copies,
* modified in place) as the estimation did -- the donors are the estimation
* sample, weighted by the weight of the estimation -- so that estat engel and
* predict see the same prices.  Without keepout, the working sample is
* restricted to the estimation sample (estat engel); with it, households
* outside keep their own prices (predict).  Nothing to do without pimpute().
program _easi_repimpute
	syntax varlist, ES(varname) TOUSE(varname) [KEEPOUT]
	if "`e(pimpute)'" == "" exit
	tempvar wpi es2
	if "`e(wtype)'" != "" qui gen double `wpi' `e(wexp)' if `es'
	else                  qui gen double `wpi' = 1       if `es'
	qui gen byte `es2' = `es'
	quietly _easi_pimpute `varlist', touse(`es2') wt(`wpi') groups(`e(pimpute)')
	if "`keepout'" == "" qui replace `touse' = 0 if !`es'
end

*==========================================================================
version 14.2
mata:
mata set matastrict on


// R's round() breaks ties to even; Stata's round() breaks them away from zero,
// which would make round(0.5) = 1 and shift p'Bp by 1e-6 when it should be 0.
real colvector _easi_rhalfeven(real colvector x)
{
	real colvector f, r
	f = floor(x)
	r = x :- f
	return(f :+ (r :> 0.5) :+ ((r :== 0.5) :* (mod(f, 2) :!= 0)))
}

// the quantisation the R package applies to p'Ap and p'Bp
real colvector _easi_quant(real colvector x)
{
	return(_easi_rhalfeven(1e6 :* x :+ 0.5) :/ 1e6)
}

// design matrix, in the exact column order the R package produces
real matrix _easi_design(real colvector y, real matrix z, real matrix np,
	real scalar R, real scalar py, real scalar zy, real scalar pz,
	real rowvector ipz)
{
	real matrix X
	real scalar i, t, neq

	neq = cols(np)
	X   = J(rows(y), 1, 1)
	for (i = 1; i <= R; i++) X = X, y:^i
	X = X, z
	if (zy) X = X, (y :* z)
	X = X, np
	if (py) X = X, (y :* np)
	if (pz) {
		for (t = 1; t <= cols(ipz); t++) {
			for (i = 1; i <= neq; i++) X = X, (np[, i] :* z[, ipz[t]])
		}
	}
	return(X)
}

// instruments.
//
// In the R package the estimation loop (its second `while`) restarts from the
// pristine data frame and -- unlike its first loop -- never refreshes the
// instrument columns.  So the instruments that produce the reported estimates
// are frozen at y_tilda: the powers y_inst^r AND the interactions y_inst*z,
// y_inst*np.  compat mode reproduces that (Z is then constant across
// iterations); corrected mode refreshes them with the current y_inst.
real matrix _easi_instr(real colvector yp, real matrix z,
	real matrix np, real scalar R, real scalar py, real scalar zy,
	real scalar pz, real rowvector ipz)
{
	real matrix Z
	real scalar i, t, neq

	neq = cols(np)
	Z   = J(rows(yp), 1, 1), z, np
	if (pz) {
		for (t = 1; t <= cols(ipz); t++) {
			for (i = 1; i <= neq; i++) Z = Z, (np[, i] :* z[, ipz[t]])
		}
	}
	for (i = 1; i <= R; i++) Z = Z, yp:^i
	if (zy) Z = Z, (yp :* z)
	if (py) Z = Z, (yp :* np)
	return(Z)
}

// p'A(z)p and p'Bp, straight from the coefficients.  Used by the estimation
// loop and by the diagnostic; keeping one copy is what stops them drifting.
void _easi_quadf(real colvector b, real matrix np, real matrix z,
	real scalar R, real scalar py, real scalar pz, real scalar zy,
	real rowvector ipz, real scalar lgQuant,
	real colvector pAp, real colvector pBp)
{
	real matrix B, Aco
	real scalar n, neq, T, nipz, onp, oynp, onpz, i, j, t

	n = rows(np); neq = cols(np); T = cols(z); nipz = cols(ipz)
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	B    = rowshape(b', neq)'

	pAp = J(n, 1, 0)
	pBp = J(n, 1, 0)

	Aco = B[| onp + 1, 1 \ onp + neq, neq |]'
	for (i = 1; i <= neq; i++) {
		for (j = 1; j <= neq; j++) {
			pAp = pAp + Aco[i, j] :* np[, i] :* np[, j]
		}
	}
	if (pz) {
		for (t = 1; t <= nipz; t++) {
			Aco = B[| onpz + (t-1)*neq + 1, 1 \ onpz + t*neq, neq |]'
			for (i = 1; i <= neq; i++) {
				for (j = 1; j <= neq; j++) {
					pAp = pAp + Aco[i, j] :*
						np[, i] :* np[, j] :* z[, ipz[t]]
				}
			}
		}
	}
	if (py) {
		Aco = B[| oynp + 1, 1 \ oynp + neq, neq |]'
		for (i = 1; i <= neq; i++) {
			for (j = 1; j <= neq; j++) {
				pBp = pBp + Aco[i, j] :* np[, i] :* np[, j]
			}
		}
	}
	if (lgQuant) {
		pAp = _easi_quant(pAp)
		pBp = _easi_quant(pBp)
	}
}

// cross-equation Slutsky symmetry:  eq_i_<block_j> - eq_j_<block_i> = 0
real matrix _easi_restrict(real scalar neq, real scalar k, real scalar onp,
	real scalar py, real scalar oynp, real scalar pz, real scalar onpz,
	real scalar nipz)
{
	real matrix Rr
	real rowvector row
	real scalar i, j, t, K

	K  = neq * k
	Rr = J(0, K, 0)
	for (i = 1; i <= neq - 1; i++) {
		for (j = i + 1; j <= neq; j++) {
			row = J(1, K, 0)
			row[(i - 1) * k + onp + j] =  1
			row[(j - 1) * k + onp + i] = -1
			Rr = Rr \ row
		}
	}
	if (py) {
		for (i = 1; i <= neq - 1; i++) {
			for (j = i + 1; j <= neq; j++) {
				row = J(1, K, 0)
				row[(i - 1) * k + oynp + j] =  1
				row[(j - 1) * k + oynp + i] = -1
				Rr = Rr \ row
			}
		}
	}
	if (pz) {
		for (t = 1; t <= nipz; t++) {
			for (i = 1; i <= neq - 1; i++) {
				for (j = i + 1; j <= neq; j++) {
					row = J(1, K, 0)
					row[(i-1)*k + onpz + (t-1)*neq + j] =  1
					row[(j-1)*k + onpz + (t-1)*neq + i] = -1
					Rr = Rr \ row
				}
			}
		}
	}
	return(Rr)
}

// Restricted 3SLS, reproducing systemfit's defaults exactly:
//   maxiter = 1            -> restricted 2SLS, then ONE GLS step
//   methodResidCov =        "geomean", computed from the 2SLS residuals
//   residCovRestricted =    TRUE, method3sls = "GLS"
// Every equation shares X and Z, so the stacked normal matrix is I # A and we
// never form the (n*neq) stack.
// The KKT system of the restricted GLS step borders the normal block by the
// symmetry restrictions.  The block has entries of order Sigma^-1 * A, which
// with small budget shares runs to 1e7 and beyond, while the restriction rows
// are +-1 -- fourteen orders of magnitude inside one matrix, and LU loses all
// of them.  Measured on a 3-good cereals system with mean shares 0.02:
// cond(Sigma) = 1.5, cond(Sigma^-1 # A) = 343, cond(bordered) = 2.7e14, and
// lusolve() returns missing.  It does not bite on data with shares near 0.1,
// which is why it went unnoticed.
//
// Scaling the restriction rows fixes it and changes nothing: the constraint
// set is the same, and the block of the inverse we use as the covariance,
// P = A^-1 - A^-1 R'(R A^-1 R')^-1 R A^-1, is invariant to R -> cR.  The
// right-hand side of those rows is zero, so it needs no rescaling either.
real matrix _easi_border(real matrix M, real matrix Rr)
{
	real scalar sc
	sc = trace(M) / rows(M)
	if (sc <= 0 | sc >= .) sc = 1
	return((M, sc :* Rr' \ sc :* Rr, J(rows(Rr), rows(Rr), 0)))
}
real colvector _easi_ksolve(real matrix W, real colvector rhs)
{
	real colvector sol
	sol = lusolve(W, rhs)
	if (hasmissing(sol)) sol = qrsolve(W, rhs)
	return(sol)
}
real matrix _easi_kinv(real matrix W, real scalar K)
{
	real matrix Wi
	Wi = luinv(W)
	if (hasmissing(Wi)) Wi = pinv(W)
	return(Wi[| 1, 1 \ K, K |])
}

void _easi_3sls(real matrix X, real matrix Z, real matrix Y, real matrix Rr,
	real colvector wt, real scalar neff, real scalar vtype,
	real colvector clid, real colvector strid, real colvector fpc,
	real scalar maxsig, real colvector b, real matrix V, real matrix Sigma,
	real matrix Pout)
{
	real matrix Xhat, A, Ry, W, Si, B, E, ZZiZX, P, Bmat
	real colvector rhs, sol, bold
	real scalar k, neq, K, r, it

	k = cols(X); neq = cols(Y); K = neq * k; r = rows(Rr)

	// weighted first stage: Xhat = Z (Z'WZ)^-1 Z'WX
	ZZiZX = cholsolve(quadcross(Z, wt, Z), quadcross(Z, wt, X))
	if (hasmissing(ZZiZX)) ZZiZX = qrsolve(quadcross(Z, wt, Z), quadcross(Z, wt, X))
	Xhat = Z * ZZiZX
	A    = quadcross(Xhat, wt, Xhat)
	Ry   = quadcross(Xhat, wt, Y)

	// step 1 -- restricted 2SLS
	W   = _easi_border(I(neq) # A, Rr)
	rhs = vec(Ry) \ J(r, 1, 0)
	sol = _easi_ksolve(W, rhs)
	b   = sol[1::K]

	// steps 2-3 -- Sigma from the current residuals (geomean divisor: k_i = k
	// for all i, because the symmetry restrictions are cross-equation and so
	// do not reduce any single equation's coefficient count), then a GLS step.
	//
	// maxsig = 1 is what systemfit and reg3 do, and what -legacy(sigma1)-
	// reproduces.  Iterating Sigma to convergence gives the ML estimator of
	// the singular system which, unlike the one-step version, is invariant to
	// which good is dropped (Barten 1969); see tests/test_step10.do.
	for (it = 1; it <= maxsig; it++) {
		bold  = b
		B     = rowshape(b', neq)'
		E     = Y - X * B
		Sigma = quadcross(E, wt, E) :/ (neff - k)
		Si    = invsym(Sigma)
		W     = _easi_border(Si # A, Rr)
		rhs   = vec(Ry * Si) \ J(r, 1, 0)
		sol   = _easi_ksolve(W, rhs)
		b     = sol[1::K]
		if (it > 1 & max(abs(b - bold)) < 1e-11) break
	}

	// P plays the role of A^-1 under the restrictions (P A P = P)
	P = _easi_kinv(W, K)
	P = (P + P') / 2
	Pout = P

	if (vtype == 0) {
		V = P
		return
	}

	// sandwich: the estimator is linear in the moment vector Xhat'Omega^-1 y,
	// so V = P * Var(score) * P.
	E    = Y - X * rowshape(b', neq)'
	Bmat = _easi_bmat(Xhat, E, Si, wt, vtype, clid, strid, fpc, k, neq)
	V    = P * Bmat * P
	V    = (V + V') / 2
}

// Design-based variance of a sum of per-observation contributions.  The input
// rows must already be the contributions themselves -- scores for a
// coefficient variance, influence functions for anything else -- so the same
// three designs serve both.
real matrix _easi_agg(real matrix S, real scalar vtype, real colvector clid,
	real colvector strid, real colvector fpc)
{
	real matrix G, D, info, sinfo, Gh, B, Ss
	real colvector perm, sg, st, fg, fv
	real rowvector mh
	real scalar c, h, nc, ns, nh, K, fh

	K = cols(S)
	// robust: the small-sample factor n/(n-1), the one cluster (nc/(nc-1))
	// and svy (nh/(nh-1)) apply below, so that robust is the cluster
	// variance with one observation per cluster and the svy variance of a
	// simple random sample (easi 2.0.0; 1.0.0 had no factor).  The legacy
	// paths (vtype 0) keep the plain cross-product.
	if (vtype == 1) return((rows(S) / (rows(S) - 1)) * quadcross(S, S))
	if (vtype <= 0) return(quadcross(S, S))

	if (vtype == 2) {
		perm = order(clid, 1)
		Ss   = S[perm, .]
		info = panelsetup(clid[perm], 1)
		nc   = rows(info)
		G    = J(nc, K, 0)
		for (c = 1; c <= nc; c++) {
			G[c, .] = colsum(panelsubmatrix(Ss, c, info))
		}
		return((nc / (nc - 1)) * quadcross(G, G))
	}

	perm = order((strid, clid), (1, 2))
	Ss   = S[perm, .]
	st   = strid[perm]
	fv   = (rows(fpc) == rows(S) ? fpc[perm] : J(rows(S), 1, 0))
	info = panelsetup(clid[perm], 1)
	nc   = rows(info)
	G    = J(nc, K, 0)
	sg   = J(nc, 1, 0)
	fg   = J(nc, 1, 0)
	for (c = 1; c <= nc; c++) {
		G[c, .] = colsum(panelsubmatrix(Ss, c, info))
		sg[c]   = st[info[c, 1]]
		fg[c]   = fv[info[c, 1]]
	}
	sinfo = panelsetup(sg, 1)
	ns    = rows(sinfo)
	B     = J(K, K, 0)
	for (h = 1; h <= ns; h++) {
		Gh = panelsubmatrix(G, h, sinfo)
		nh = rows(Gh)
		if (nh > 1) {
			mh = colsum(Gh) / nh
			D  = Gh :- mh
			// Stata's convention: a value at or below one is a
			// sampling RATE, anything larger is a population SIZE.
			// Only the first stage matters -- that is exactly what
			// makes the ultimate-cluster approximation legitimate.
			fh = fg[sinfo[h, 1]]
			if (fh > 1) fh = nh / fh
			if (fh <= 0 | fh >= 1 | missing(fh)) fh = 0
			B = B + (1 - fh) * (nh / (nh - 1)) * quadcross(D, D)
		}
	}
	return(B)
}

// Cross-product of the scores.  The score of observation t is
// w_t * (u_t # xhat_t') with u_t = Sigma^-1 e_t'; -vce(cluster)- sums them
// within cluster first and applies the usual nc/(nc-1) correction.
real matrix _easi_bmat(real matrix Xhat, real matrix E, real matrix Si,
	real colvector wt, real scalar vtype, real colvector clid,
	real colvector strid, real colvector fpc, real scalar k, real scalar neq)
{
	return(_easi_agg(_easi_scores(Xhat, E, Si, wt, k, neq), vtype, clid,
		strid, fpc))
}

// Per-observation score contributions, wt_t * (u_t # xhat_t') with
// u_t = Sigma^-1 e_t'.  Shared by the coefficient variance and by the
// influence function of beta.
real matrix _easi_scores(real matrix Xhat, real matrix E, real matrix Si,
	real colvector wt, real scalar k, real scalar neq)
{
	real matrix S, U
	real scalar j, n, K

	n = rows(Xhat); K = neq * k
	U = E * Si
	S = J(n, K, 0)
	for (j = 1; j <= neq; j++) {
		S[| 1, (j - 1) * k + 1 \ n, j * k |] = U[, j] :* Xhat
	}
	return(wt :* S)
}

// Influence function of beta-hat, one row per observation: IF_t = Pu * s_t,
// where s_t is the score contribution and Pu the projection that produced the
// reported covariance.  Summing the rows and aggregating them by the design
// reproduces V exactly, which is what makes this the right object to combine
// with the influence of a sample mean.
void _easi_ifb(real matrix X, real matrix Z, real matrix Y, real colvector wt,
	real matrix Sigma, real colvector b, real matrix Pu, real matrix IFb)
{
	real matrix ZZiZX, Xhat, E, Si
	real scalar k, neq

	k = cols(X); neq = cols(Y)
	ZZiZX = cholsolve(quadcross(Z, wt, Z), quadcross(Z, wt, X))
	if (hasmissing(ZZiZX)) ZZiZX = qrsolve(quadcross(Z, wt, Z), quadcross(Z, wt, X))
	Xhat = Z * ZZiZX
	E    = Y - X * rowshape(b', neq)'
	Si   = invsym(Sigma)
	IFb  = _easi_scores(Xhat, E, Si, wt, k, neq) * Pu'
}

// ------------------------------------------- generated-regressor variance
//
// The design contains y, which is itself a function of beta, so the 3SLS
// covariance -- the one systemfit, reg3 and the R package report -- is that of
// the last linear step CONDITIONAL on y.  Pendakur's code says as much, twice,
// immediately before its final reg3:
//     *note that reported standard errors are wrong for iterated estimates
// A bootstrap audit puts the shortfall at 23% on average for the coefficients
// on y, and at nothing at all for the other 96 (tests/AUDIT.md, section I).
//
// Write the estimator as GMM in beta with fixed instruments and residual
// e_im(beta) = w_im - x_i(y(beta))'b_m.  The moment vector has blocks
//     m_j(beta) = sum_i wt_i xhat_i sum_m sigma^{jm} e_im(beta)
// and, since dx_i/dy = d_i and b_m'd_i = aleph_im (the semi-elasticity of good
// m, already computed for the elasticities),
//     de_im/db_l' = -[ x_i' 1{m=l} + aleph_im * (dy_i/db_l') ] .
// Hence  -dm/dbeta' = G = (Sigma^-1 # A) + C  with
//     C row-block j = sum_i wt_i v_ij xhat_i q_i',
//     v_ij = sum_m sigma^{jm} aleph_im,        q_i = dy_i/dbeta .
// The conditional covariance is the special case C = 0.
//
// dy/dbeta is closed form.  y = (y_stone + p'Ap/2)/(1 - p'Bp/2) with p'Ap and
// p'Bp both LINEAR in beta, so with D_i = 1 - p'Bp_i/2,
//     q_i = (1/(2 D_i)) * [ d(p'Ap)/dbeta + y_i d(p'Bp)/dbeta ]
//         = (1/(2 D_i)) * (np_i # g_i),
// g_i being nonzero only on the price positions of the equation's block:
// np_ip at onp+p, y_i*np_ip at oynp+p, np_ip*z_is at onpz+(s-1)neq+p.
//
// -legacy(condse)- keeps the conditional covariance; -compat- turns it on.
void _easi_gencorr(real matrix X, real matrix Z, real matrix Y, real matrix Rr,
	real colvector wt, real matrix np, real matrix z, real colvector y,
	real colvector pBp, real matrix Sigma, real colvector b, real scalar R,
	real scalar py, real scalar zy, real scalar pz, real rowvector ipz,
	real scalar vtype, real colvector clid, real colvector strid,
	real colvector fpc, real matrix V, real matrix Pout)
{
	real matrix Xhat, ZZiZX, A, Bm, E, Dx, Al, Vv, Q, C, Gj, W, PG, Bv, Si
	real scalar n, k, T, neq, K, nipz, r, i, onp, oynp, onpz

	n = rows(X); k = cols(X); neq = cols(Y); K = neq * k
	T = cols(z); nipz = cols(ipz); r = rows(Rr)
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)

	ZZiZX = cholsolve(quadcross(Z, wt, Z), quadcross(Z, wt, X))
	if (hasmissing(ZZiZX)) ZZiZX = qrsolve(quadcross(Z, wt, Z), quadcross(Z, wt, X))
	Xhat = Z * ZZiZX
	A    = quadcross(Xhat, wt, Xhat)
	Si   = invsym(Sigma)

	Bm = rowshape(b', neq)'				// k x neq
	E  = Y - X * Bm

	// d_i = dx_i/dy, then the semi-elasticities aleph and v_ij
	Dx = J(n, k, 0)
	for (i = 1; i <= R; i++) Dx[, 1 + i] = i :* y:^(i - 1)
	if (zy) Dx[| 1, 1 + R + T + 1 \ n, 1 + R + 2 * T |] = z
	if (py) Dx[| 1, oynp + 1 \ n, oynp + neq |] = np
	Al = Dx * Bm					// n x neq
	Vv = Al * Si					// n x neq

	// q_i = dy_i/dbeta
	Q = _easi_dydb(np, z, y, pBp, R, py, pz, zy, ipz)

	C = J(K, K, 0)
	for (i = 1; i <= neq; i++) {
		C[| (i - 1) * k + 1, 1 \ i * k, K |] = quadcross(Xhat, wt :* Vv[, i], Q)
	}

	Gj = (Si # A) + C
	W  = _easi_border(Gj, Rr)
	PG = _easi_kinv(W, K)
	Pout = PG

	Bv = (vtype == 0 ? (Si # A) :
		_easi_bmat(Xhat, E, Si, wt, vtype, clid, strid, fpc, k, neq))
	V  = PG * Bv * PG'
	V  = (V + V') / 2
}

// dy_i/dbeta, one row per observation, in free coordinates: with D_i =
// 1 - p'Bp_i/2, q_i = (1/(2 D_i)) * (np_i # g_i), g_i nonzero only on the
// price positions of the equation's block (np_ip at onp+p, y_i*np_ip at
// oynp+p, np_ip*z_is at onpz+(s-1)neq+p).  Shared by the generated-regressor
// term of the coefficient variance and by the y term of the elasticity
// Jacobian (_easi_ejacy).
real matrix _easi_dydb(real matrix np, real matrix z, real colvector y,
	real colvector pBp, real scalar R, real scalar py, real scalar pz,
	real scalar zy, real rowvector ipz)
{
	real matrix g, Q
	real colvector fac
	real scalar n, neq, T, nipz, k, onp, oynp, onpz, e, t

	n = rows(np); neq = cols(np); T = cols(z); nipz = cols(ipz)
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)

	g = J(n, k, 0)
	g[| 1, onp + 1 \ n, onp + neq |] = np
	if (py) g[| 1, oynp + 1 \ n, oynp + neq |] = y :* np
	if (pz) {
		for (t = 1; t <= nipz; t++) {
			g[| 1, onpz + (t - 1) * neq + 1 \ n, onpz + t * neq |] =
				np :* z[, ipz[t]]
		}
	}
	fac = 0.5 :/ (1 :- 0.5 :* pBp)

	Q = J(n, neq * k, 0)
	for (e = 1; e <= neq; e++) {
		Q[| 1, (e - 1) * k + 1 \ n, e * k |] = (fac :* np[, e]) :* g
	}
	return(Q)
}

// Condition number of a cross-product matrix, scaled by its own diagonal.
// The raw condition number of X'X mixes the scale of the columns with their
// collinearity; scaling by sqrt(diag) leaves the collinearity, which is what
// we want to report.
real scalar _easi_scond(real matrix G)
{
	real colvector d
	d = sqrt(diagonal(G))
	d = d :+ (d :== 0)
	return(cond(G :/ (d * d')))
}

// weighted standard deviation of a column
real scalar _easi_wsd(real colvector x, real colvector wt)
{
	real scalar m, v
	m = _easi_wm(x, wt)
	v = _easi_wm((x :- m):^2, wt)
	return(v > 0 ? sqrt(v * rows(x) / max((1, rows(x) - 1))) : 0)
}

// Screening of the design columns, with NO estimation anywhere.
//
// The two things that stop the EASI fixed point converging are a regressor
// with almost no variance and a near-linear dependency among regressors, and
// both live in the design matrix alone.  The design is evaluated at the
// STARTING point y = y_stone, which is a function of the data only, so this
// runs even when the model cannot be fitted at all -- which is exactly when it
// is needed.
//
// Two instruments, because they answer different questions.  The VIF of column
// j, the j-th diagonal element of the inverse correlation matrix, says how much
// of that column is explained by the others.  It flags a variable but not the
// company it keeps.  The condition indexes of Belsley, Kuh and Welsch do name
// the company: scale every column to unit length, take the eigenvalues of the
// cross-product, and a small eigenvalue is a near-dependency whose condition
// index sqrt(lmax/l) measures its severity.  The variance-decomposition
// proportions then say which columns load on it.  A component with an index
// above 30 on which two or more columns carry more than half their variance is
// the textbook signature of harmful collinearity.
void _easi_screen(real matrix X, real colvector wt)
{
	real matrix G, Gs, Rc, Ri, V, VP, CI, OUT
	real colvector d, m, sd, lam, vif, tot
	real scalar k, j, i, lmax, flo

	k = cols(X)

	m  = J(k, 1, 0)
	sd = J(k, 1, 0)
	for (j = 1; j <= k; j++) {
		m[j]  = _easi_wm(X[, j], wt)
		sd[j] = _easi_wsd(X[, j], wt)
	}

	// ---- VIF, from the correlation matrix of the non-constant columns
	vif = J(k, 1, .)
	Rc  = J(k, k, 0)
	for (i = 1; i <= k; i++) {
		for (j = 1; j <= k; j++) {
			if (sd[i] > 0 & sd[j] > 0) {
				Rc[i, j] = (_easi_wm(X[, i] :* X[, j], wt)
					- m[i] * m[j]) / (sd[i] * sd[j])
			}
		}
	}
	for (j = 1; j <= k; j++) if (sd[j] == 0) Rc[j, j] = 1
	Ri = invsym(Rc)
	for (j = 1; j <= k; j++) {
		if (sd[j] > 0) vif[j] = (Ri[j, j] != 0 ? Ri[j, j] : .)
	}
	st_matrix("_d_scr", (m, sd, vif))

	// ---- condition indexes and variance proportions
	G  = quadcross(X, wt, X)
	d  = sqrt(diagonal(G))
	d  = d :+ (d :== 0)
	Gs = G :/ (d * d')
	symeigensystem(Gs, V, lam)
	lam  = lam'
	lmax = max(lam)

	// An exactly collinear group gives an eigenvalue that comes out as zero or
	// slightly NEGATIVE in floating point (-5.5e-17 was observed).  Guarding
	// with lam > 0 would zero out precisely the component that names the
	// culprits, so the eigenvalue is floored at a tiny multiple of the largest
	// instead: the proportions then load on the offending group as they should,
	// and the condition index comes out large and finite rather than missing.
	flo = 1e-16 * lmax
	lam = (lam :< flo) :* flo + (lam :>= flo) :* lam

	VP = J(k, k, 0)
	for (j = 1; j <= k; j++) {
		for (i = 1; i <= k; i++) {
			VP[i, j] = V[j, i]^2 / lam[i]
		}
	}
	tot = colsum(VP)'
	for (j = 1; j <= k; j++) {
		if (tot[j] > 0) VP[, j] = VP[, j] :/ tot[j]
	}
	CI = sqrt(lmax :/ lam)
	st_matrix("_d_ci", CI)
	st_matrix("_d_vp", VP)
	st_numscalar("_d_cimax", max(CI))
}

// Structural sources of a flat objective, specific to EASI, all computable
// from the data alone.
//
//  a. THE POLYNOMIAL BASIS.  The design carries 1, y, y^2, ..., y^R, a
//     Vandermonde basis whose conditioning explodes when the dispersion of y is
//     small relative to its level.  Measured: centred, each extra power costs a
//     factor of about 4; uncentred, a factor of about 100.  With power(5), the
//     default, that is 2.1e+02 on centred data and 1.1e+11 on raw levels -- past
//     what double precision can carry.  Reporting the number per power tells the
//     user exactly how far they can push R on their own data.
//
//  b. DISTINCT PRICE REGIMES.  When prices are regional or temporal indexes,
//     every household in a cell shares one price vector, and the price
//     parameters are identified off the number of cells, not the sample size.
//     On the Canadian reference data that is 48 vectors for 4,847 households,
//     against 36 parameters in the symmetric A block alone.
//
//  c. RANK OF THE PRICE VARIATION.  If prices move proportionally, var(np) is
//     near rank one and A is identified in only a few directions.
//
//  d. FIRST-STAGE STRENGTH.  y is instrumented by the mean-share Stone index.
//     A weak first stage flattens the objective in the direction of the
//     coefficients on y, which are the ones that shape the Engel curves.
//
//  e. CONDITIONING OF THE RESIDUAL COVARIANCE, proxied before estimation by the
//     covariance of the estimated shares themselves.
void _easi_struct(real colvector ystone, real colvector ytil, real matrix np,
	real matrix s, real colvector wt, real scalar R)
{
	real matrix Yp, Gy, C, Vp, Z, S
	real colvector d, lam, res, PC
	real scalar n, neq, r

	n = rows(np); neq = cols(np)

	PC = J(R, 1, .)
	Yp = J(n, 1, 1)
	for (r = 1; r <= R; r++) {
		Yp = Yp, ystone:^r
		Gy = quadcross(Yp, wt, Yp)
		d  = sqrt(diagonal(Gy)); d = d :+ (d :== 0)
		PC[r] = cond(Gy :/ (d * d'))
	}
	st_matrix("_d_poly", PC)

	st_numscalar("_d_nprice", rows(uniqrows(np)))
	st_numscalar("_d_npar",   neq * (neq + 1) / 2)

	C = variance(np)
	symeigensystem(C, Vp, lam)
	lam = lam'
	st_numscalar("_d_rkp", sum((lam :/ max(lam)) :> 0.01))
	st_numscalar("_d_pc1", max(lam) / sum(lam))

	Z   = J(n, 1, 1), np, ytil
	res = ystone - Z * qrsolve(quadcross(Z, wt, Z), quadcross(Z, wt, ystone))
	st_numscalar("_d_fs", 1 - variance(res) / variance(ystone))

	S = variance(s[, 1::neq])
	d = sqrt(diagonal(S)); d = d :+ (d :== 0)
	st_numscalar("_d_conds", cond(S :/ (d * d')))
}

// Standard error of the expenditure elasticities, by the FULL delta method.
//
//     eta_j = 1 + alephbar_j / wbar_j
//
// alephbar_j is a weighted mean of a per-household quantity that depends on
// beta; wbar_j is a weighted mean of the observed shares.  Both are estimators,
// so both must enter -- gradient times variance, cross term included.  Taking
// the Jacobian over beta alone, with wbar held fixed, keeps one term of three.
//
// The omission is worth +18% and +15% on the two small shares of a real survey,
// and it is not a missing variance: corr(alephbar, wbar) = -0.6 while alephbar
// is negative, so the cross term is NEGATIVE.  Numerator and denominator move
// together, the ratio is steadier than either, and ignoring that OVERSTATES the
// standard error -- the opposite of the usual intuition.
//
// Closed form throughout.  With alpha_hj = N_hj / Phi_h,
//     N_hj  = sum_r r b_jr y^(r-1) + sum_t h_jt z_t + sum_k b_jk np_k,
//     Phi_h = 1 - p'Bp_h / 2,
// both linear in beta,
//     d alpha_hj/d beta = (1/Phi) dN/dbeta - (N/Phi^2) dPhi/dbeta,
//     dN_hj/dbeta  = e_j (x) Dx[h,.],
//     dPhi_h/dbeta = -(1/2) (np_h (x) gB_h),
// with gB_h zero outside the y*np block, where it holds np_h.  Without the
// price-expenditure interaction Phi is one and the second term disappears.
//
// The influence function per household is then
//     IF_h(eta_j) = wt_h/sum(wt) * [ (alpha_hj - alephbar_j)/wbar_j
//                                    - (alephbar_j/wbar_j^2)(s_hj - wbar_j) ]
//                   + (1/wbar_j) (d alephbar_j/d beta) IF_h(beta),
// and its variance goes through the same aggregator as everything else, so
// robust, cluster and svy need no special case.
void _easi_eise(real matrix s, real matrix np, real colvector pJ, real matrix z,
	real colvector y, real colvector wt, real matrix IFb, real colvector Dn2,
	real rowvector ws, real colvector b, real scalar R, real scalar py,
	real scalar zy, real scalar vtype, real colvector clid,
	real colvector strid, real colvector fpc, real scalar lgB,
	real rowvector EIse)
{
	real matrix Gb, Tc, M, Vh, Dx, Dh, Bm, Al, DxP
	real colvector alp, c, Phi, Phi2, sw, tmp, v
	real rowvector gmean, gmeanP
	real scalar n, Jg, neq, k, K, T_, oynp, j, e, i, sww, ah, gJ, lo, hi

	n = rows(s); Jg = cols(s); neq = Jg - 1
	T_ = cols(z); k = cols(IFb) / neq; K = neq * k
	oynp = 1 + R + T_ + (zy ? T_ : 0) + neq
	sww  = quadsum(wt)
	Phi  = (lgB ? J(n, 1, 1) : Dn2)
	Phi2 = Phi:^2

	// dx/dy, the same object the semi-elasticities use
	Dx = J(n, k, 0)
	for (i = 1; i <= R; i++) Dx[, 1 + i] = i :* y:^(i - 1)
	if (zy) Dx[| 1, 1 + R + T_ + 1 \ n, 1 + R + 2 * T_ |] = z
	if (py) Dx[| 1, oynp + 1 \ n, oynp + neq |] = np

	Bm = rowshape(b', neq)'
	Al = Dx * Bm
	Al = Al, (-rowsum(Al))		// the dropped good, by adding up

	// d Dn2 / d beta.  The y-power and z*y coefficients enter Dn2 through
	// (C_j + D_j) np_j, so their block e is np_e * Dx[h,.]: Dh is Dx with
	// the y*np block zeroed, and that block is handled separately below.
	// Treating Dn2 as 1 - p'Bp/2 -- which is the denominator of y, a
	// different object -- gets both the sign and the missing terms wrong.
	Dh = Dx
	if (py & !lgB) Dh[| 1, oynp + 1 \ n, oynp + neq |] = J(n, neq, 0)

	DxP   = Dx :/ Phi
	gmean = quadcross(wt, DxP) :/ sww

	// The point estimate takes the y*np term of good i as P * bjk[, i]
	// over ALL J goods with bjk completed by adding up, and that
	// completion is not symmetric: the free entry b[a,e] enters good e
	// with P_a and good a with -P_J.  Al = Dx * Bm above uses np_a = P_a -
	// P_J instead; same value once symmetry holds at the estimate, but a
	// different entry-wise derivative.  gmeanP carries P_a on the own
	// block, gJ the -P_J on entry (e, j) of every block e.  The two agree
	// on the diagonal, and their difference lies along the symmetry
	// constraints -- so G V G' is unaffected, but the Jacobian is not
	// equal to the finite-difference one unless it is written this way.
	gmeanP = gmean
	gJ = 0
	if (py) {
		gmeanP[| 1, oynp + 1 \ 1, oynp + neq |] =
			quadcross(wt, (np :+ pJ) :/ Phi) :/ sww
		gJ = quadcross(wt, pJ :/ Phi) / sww
	}

	Gb = J(Jg, K, 0)
	Tc = J(n, Jg, 0)
	sw = wt :/ sww

	for (j = 1; j <= Jg; j++) {
		alp = Al[, j] :/ Phi
		ah  = quadcross(wt, alp) / sww

		if (j <= neq) {
			Gb[| j, (j-1)*k + 1 \ j, j*k |] = gmeanP
			if (py) {
				for (e = 1; e <= neq; e++) {
					Gb[j, (e-1)*k + oynp + j] = Gb[j, (e-1)*k + oynp + j] - gJ
				}
			}
		}
		else {
			for (e = 1; e <= neq; e++) {
				Gb[| j, (e-1)*k + 1 \ j, e*k |] = -gmean
			}
		}
		if (!lgB) {
			c = wt :* (Al[, j] :/ Phi2)
			M = quadcross(np, c, Dh) :/ sww
			if (py) {
				// The y*np coefficients enter Dn2 through the quadratic
				// form (1/2) sum_ab bjk[a,b] P_a P_b over ALL J goods,
				// with bjk completed by adding up.  That completion is
				// not symmetric: the free entry b[j,e] also sits at
				// (j,J), (J,j) and (J,J), so
				//     d Dn2 / d b[j,e] = (P_j P_e - 2 P_j P_J + P_J^2) / 2
				//                      = np_j np_e / 2 + P_J (np_e - np_j) / 2.
				// The second piece is antisymmetric in (j,e) -- it lies
				// along the Slutsky-symmetry constraints, so it cannot
				// move G V G' -- but the entry-wise Jacobian must carry
				// it to equal the finite-difference one.  Verified to
				// 1e-10 per entry against the coded Dn2.
				v = 0.5 :* quadcross(np :* pJ, c)
				M[| 1, oynp + 1 \ neq, oynp + neq |] =
					(0.5 :* quadcross(np, c, np)
					 + v * J(1, neq, 1) - J(neq, 1, 1) * v') :/ sww
			}
			for (e = 1; e <= neq; e++) {
				lo = (e-1)*k + 1
				hi = e*k
				Gb[| j, lo \ j, hi |] = Gb[| j, lo \ j, hi |] - M[e, .]
			}
		}
		Gb[j, .] = Gb[j, .] :/ ws[j]

		tmp = (alp :- ah) :/ ws[j] :- (ah / (ws[j]^2)) :* (s[, j] :- ws[j])
		Tc[, j] = sw :* tmp
	}

	Tc = Tc + IFb * Gb'
	Vh = _easi_agg(Tc, vtype, clid, strid, fpc)
	EIse = sqrt(diagonal(Vh))'
}

// ----------------------------------------------------------- diagnostic
//
// Everything here is cheap: two restricted GLS steps and a handful of
// cross-products.  It answers the four questions that decide whether an EASI
// specification will behave, in the order in which they bite.
//
//  1. Are the data admissible at all -- shares summing to one, prices usable.
//
//  2. Is the model NORMALISED.  EASI is written for log prices and log
//     expenditure centred around a base period.  Fed raw levels, the columns
//     1, y, y^2, y^3 of the design correlate above 0.999 and the normal matrix
//     goes singular; on the Mexican survey cond(X'X) was 4.5e+10 in levels and
//     4.5e+01 centred.  Reporting both numbers turns a baffling r(504) into
//     one line of advice.
//
//  3. Is anything too small or too flat.  A price that does not vary carries
//     no information about A(z); a mean budget share below about 0.01
//     multiplies every semi-elasticity by 1/w and degrades the conditioning of
//     the restricted system.
//
//  4. Will the fixed point CONVERGE, and how fast.  The iteration is
//     y <- f(b(y)).  Running it twice from the 2SLS start gives an empirical
//     contraction factor r = |y2-y1| / |y1-y0|, and if r < 1 the number of
//     iterations to reach tol follows from log(tol/|y1-y0|)/log(r).  The
//     companion quantity is Phi = 1 - p'Bp/2, the denominator of y: as it
//     approaches zero the map stops contracting, and a negative Phi means the
//     cost function has left its regular region.
void _easi_diagm(string scalar svars, string scalar lpvars, string scalar lxvar,
	string scalar zvars, string scalar touse, string scalar wvar,
	real scalar neff, real scalar R, real scalar py, real scalar pz,
	real scalar zy, string scalar sipz, real scalar tol)
{
	real matrix s, p, z, np, X, Xc, Xh, Z, Y, Rr, V, Sigma, G, A, VT, GT, Pu
	real colvector lnx, ystone, ytil, wt, b, pAp, pBp, y1, y2, Phi
	real rowvector ms, ipz
	real scalar n, J, neq, T, k, nipz, onp, oynp, onpz, i, j
	real scalar d0, d1, ratio, nit

	s   = st_data(., svars,  touse)
	p   = st_data(., lpvars, touse)
	lnx = st_data(., lxvar,  touse)
	z   = st_data(., zvars,  touse)
	wt  = st_data(., wvar,   touse)
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))

	n = rows(s); J = cols(s); neq = J - 1; T = cols(z); nipz = cols(ipz)

	np     = p[, 1::neq] :- p[, J]
	ystone = lnx - rowsum(s :* p)
	ms     = quadcross(wt, s) :/ quadsum(wt)
	ytil   = lnx - p * ms'
	Y      = s[, 1::neq]

	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	Rr   = _easi_restrict(neq, k, onp, py, oynp, pz, onpz, nipz)

	// ---- 1. admissibility
	st_numscalar("_d_sum1", max(abs(rowsum(s) :- 1)))
	st_numscalar("_d_neg",  sum(s :< 0))
	st_numscalar("_d_gt1",  sum(s :> 1))
	st_numscalar("_d_n",    n)
	st_numscalar("_d_k",    k)
	st_numscalar("_d_K",    neq * k)

	// ---- 2. screening of the design, no estimation involved
	X = _easi_design(ystone, z, np, R, py, zy, pz, ipz)
	_easi_screen(X, wt)
	_easi_struct(ystone, ytil, np, s, wt, R)

	// ---- 3. normalisation and conditioning
	G = quadcross(X, wt, X)
	st_numscalar("_d_condR", _easi_scond(G))

	Xc = X
	for (i = 2; i <= cols(Xc); i++) {
		Xc[, i] = Xc[, i] :- _easi_wm(Xc[, i], wt)
	}
	st_numscalar("_d_condC", _easi_scond(quadcross(Xc, wt, Xc)))

	// mean, sd, |mean|/sd for log expenditure and every log price
	VT = J(1 + J, 3, 0)
	VT[1, 1] = _easi_wm(lnx, wt); VT[1, 2] = _easi_wsd(lnx, wt)
	for (j = 1; j <= J; j++) {
		VT[1 + j, 1] = _easi_wm(p[, j], wt)
		VT[1 + j, 2] = _easi_wsd(p[, j], wt)
	}
	VT[, 3] = abs(VT[, 1]) :/ (VT[, 2] :+ (VT[, 2] :== 0))
	st_matrix("_d_var", VT)

	// ---- 3. scale and variability, good by good
	GT = J(J, 3, 0)
	for (j = 1; j <= J; j++) {
		GT[j, 1] = _easi_wm(s[, j], wt)
		GT[j, 2] = _easi_wsd(s[, j], wt)
		GT[j, 3] = (j <= neq ? _easi_wsd(np[, j], wt) : 0)
	}
	st_matrix("_d_good", GT)
	st_numscalar("_d_wmin", min(GT[, 1]))
	st_numscalar("_d_npmin", min(GT[| 1, 3 \ neq, 3 |]))

	// ---- 5. identification (computed here, reported last)
	Z  = _easi_instr(ytil, z, np, R, py, zy, pz, ipz)
	Xh = Z * qrsolve(quadcross(Z, wt, Z), quadcross(Z, wt, X))
	A  = quadcross(Xh, wt, Xh)
	st_numscalar("_d_nz",   cols(Z))
	st_numscalar("_d_rz",   rank(quadcross(Z, wt, Z)))
	st_numscalar("_d_rx",   rank(A))
	st_numscalar("_d_conda", _easi_scond(A))

	// ---- 4. convergence of the fixed point, two steps
	_easi_3sls(X, Z, Y, Rr, wt, neff, 0, J(0,1,0), J(0,1,0), J(0,1,0), 1, b, V, Sigma, Pu)
	if (hasmissing(b)) {
		st_numscalar("_d_fail", 1)
		st_numscalar("_d_ratio", .); st_numscalar("_d_iter", .)
		st_numscalar("_d_phimin", .); st_numscalar("_d_phimax", .)
		st_numscalar("_d_phimean", .); st_numscalar("_d_phineg", .)
		return
	}
	st_numscalar("_d_fail", 0)
	_easi_quadf(b, np, z, R, py, pz, zy, ipz, 0, pAp, pBp)
	y1 = (ystone :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)

	X = _easi_design(y1, z, np, R, py, zy, pz, ipz)
	_easi_3sls(X, Z, Y, Rr, wt, neff, 0, J(0,1,0), J(0,1,0), J(0,1,0), 1, b, V, Sigma, Pu)
	_easi_quadf(b, np, z, R, py, pz, zy, ipz, 0, pAp, pBp)
	y2 = (ystone :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)

	Phi = 1 :- 0.5 :* pBp
	st_numscalar("_d_phimin",  min(Phi))
	st_numscalar("_d_phimax",  max(Phi))
	st_numscalar("_d_phimean", _easi_wm(Phi, wt))
	st_numscalar("_d_phineg",  sum(Phi :<= 0))

	d0 = max(abs(y1 - ystone))
	d1 = max(abs(y2 - y1))
	ratio = (d0 > 0 ? d1 / d0 : 0)
	st_numscalar("_d_ratio", ratio)
	st_numscalar("_d_step1", d0)
	nit = 2
	if (ratio >= 1) nit = .
	else if (ratio > 0) nit = ceil(log(tol / d1) / log(ratio)) + 2
	st_numscalar("_d_iter", nit)
}

// Weighted median.  With unit weights this is exactly R's median() -- sort,
// middle element, or the mean of the two middle ones -- and with integer
// frequency weights it is exactly the median of the physically expanded
// sample, which is what makes -fweight- equivalent to duplicating rows.
real scalar _easi_wmedian(real colvector x, real colvector wt)
{
	real matrix XW
	real colvector cw
	real scalar n, i, half

	n = rows(x)
	if (n == 0) return(.)
	XW   = sort((x, wt), 1)
	cw   = runningsum(XW[, 2])
	half = cw[n] / 2
	for (i = 1; i <= n; i++) {
		if (cw[i] == half & i < n) return((XW[i, 1] + XW[i + 1, 1]) / 2)
		if (cw[i] >  half)        return(XW[i, 1])
	}
	return(XW[n, 1])
}

// median over households of the pointwise delta-method standard error,
// median( sqrt( diag( MAT * DD * MAT' ) ) ), without forming the n x n matrix
real scalar _easi_dse(real matrix MAT, real matrix DD, real colvector wt)
{
	return(_easi_wmedian(sqrt(rowsum((MAT * DD) :* MAT)), wt))
}

// weighted mean of a column, and of each column of a matrix
real scalar _easi_wm(real colvector x, real colvector wt)
{
	return(quadcross(wt, x) / quadsum(wt))
}
real rowvector _easi_wms(real matrix X, real colvector wt)
{
	return(quadcross(wt, X) :/ quadsum(wt))
}

// ------------------------------------------------------------ elasticities
//
// Transcription of intermediate.blocs() + elastic() from the R package.
// Lines flagged "R bug" are reproduced verbatim under -compat- and fixed
// otherwise; each one is documented in tests/README.md.
// stack the reported elasticities into one vector, for the Jacobian
// Compensated elasticities ride along only when asked: the numerical Jacobian
// loops over the K coefficients, so extra OUTPUTS cost almost nothing, but the
// table is off by default and there is no point computing what is not shown.
real colvector _easi_estack(real rowvector EI, real matrix EPRICE, real matrix EZ,
	real matrix EPQ, real scalar doCQ)
{
	// vec(EZ), not vec(EZ'): the unstacking in _easi_elast is
	// rowshape(v', Jg)', which is the inverse of vec() for every block.
	// Stacking the transpose put the demographic standard errors in the
	// wrong cells whenever T > 1 -- with T = 3 and 9 goods, cell (2,1)
	// showed the standard error of cell (1,2).
	if (doCQ) return(EI' \ vec(EPRICE) \ vec(EZ) \ vec(EPQ))
	return(EI' \ vec(EPRICE) \ vec(EZ))
}

// The reported elasticities as a FUNCTION of the coefficient vector.  Split out
// of _easi_elast so that the same code can be re-evaluated at perturbed b, which
// is what the delta-method standard errors need.
void _easi_epoint(real colvector b, real matrix s, real matrix P, real matrix z,
	real colvector y, real colvector wt, real scalar R, real scalar py,
	real scalar zy, real scalar pz, real rowvector ipz,
	real scalar lgEZ, real scalar lgB,
	real rowvector EI, real rowvector ER, real matrix EPRICE, real matrix EP,
	real matrix EPS, real matrix EPQ, real matrix EZ,
	real rowvector ws, real colvector Dn2)
{
	pointer(real matrix) rowvector A
	real matrix B, Zc, bjk, bjr, gjt, hjt, AZ, Ypow, M2
	real rowvector mAZ, mHS
	real colvector tot2, Ci, Di, Gi, Ui, Ei, Fi, Hv, cum, den1
	real scalar n, Jg, neq, T, nipz, k, onp, oynp, onpz
	real scalar i, j, q, t, r, d, den, mG, mCDG, Bsc, Hm, acc, my, mU, sIdx

	n = rows(s); Jg = cols(s); neq = Jg - 1; T = cols(z); nipz = cols(ipz)
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)

	B  = rowshape(b', neq)'			// k x neq
	Zc = J(n, 1, 1), z			// n x (T+1)
	Ypow = J(n, R, 1)
	for (r = 2; r <= R; r++) Ypow[, r] = y:^(r - 1)
	ws = _easi_wms(s, wt)			// 1 x Jg (weighted)

	// ---- coefficient blocks, with the adding-up completion ------------
	bjr = J(R, neq + 1, 0)
	for (i = 1; i <= neq; i++) {
		for (r = 1; r <= R; r++) bjr[r, i] = B[1 + r, i]
	}
	for (r = 1; r <= R; r++) bjr[r, neq + 1] = -sum(bjr[r, 1::neq])

	gjt = J(T, neq + 1, 0)
	for (i = 1; i <= neq; i++) {
		for (t = 1; t <= T; t++) gjt[t, i] = B[1 + R + t, i]
	}
	for (t = 1; t <= T; t++) gjt[t, neq + 1] = -sum(gjt[t, 1::neq])

	hjt = J(T, neq + 1, 0)
	if (zy) {
		for (i = 1; i <= neq; i++) {
			for (t = 1; t <= T; t++) hjt[t, i] = B[1 + R + T + t, i]
		}
		for (t = 1; t <= T; t++) hjt[t, neq + 1] = -sum(hjt[t, 1::neq])
	}

	bjk = J(neq + 1, neq + 1, 0)
	if (py) {
		for (i = 1; i <= neq; i++) {
			for (j = 1; j <= neq; j++) bjk[j, i] = B[oynp + j, i]
		}
		for (j = 1; j <= neq; j++) bjk[j, neq + 1] = -sum(bjk[j, 1::neq])
		for (j = 1; j <= neq; j++) bjk[neq + 1, j] = bjk[j, neq + 1]
		bjk[neq + 1, neq + 1] = -sum(bjk[neq + 1, 1::neq])
	}

	// a[t, k, i] -- one (T+1) x (neq+1) matrix per good
	A = J(1, neq + 1, NULL)
	for (i = 1; i <= neq + 1; i++) A[i] = &(J(T + 1, neq + 1, 0))
	for (i = 1; i <= neq; i++) {
		for (j = 1; j <= neq; j++) (*A[i])[1, j] = B[onp + j, i]
	}
	for (j = 1; j <= neq; j++) {
		acc = 0
		for (i = 1; i <= neq; i++) acc = acc + (*A[i])[1, j]
		(*A[neq + 1])[1, j] = -acc
	}
	for (i = 1; i <= neq; i++) (*A[i])[1, neq + 1] = (*A[neq + 1])[1, i]
	(*A[neq + 1])[1, neq + 1] = -sum((*A[neq + 1])[1, 1::neq])

	if (pz) {
		for (i = 1; i <= neq; i++) {
			for (sIdx = 1; sIdx <= nipz; sIdx++) {
				for (q = 1; q <= neq; q++) {
					(*A[i])[ipz[sIdx] + 1, q] =
						B[onpz + (sIdx - 1) * neq + q, i]
				}
			}
		}
		for (sIdx = 1; sIdx <= nipz; sIdx++) {
			t = ipz[sIdx] + 1
			for (i = 1; i <= neq; i++) {
				acc = 0
				for (j = 1; j <= neq; j++) acc = acc + (*A[j])[t, i]
				(*A[neq + 1])[t, i] = -acc
			}
			for (i = 1; i <= neq; i++) (*A[i])[t, neq + 1] = (*A[neq + 1])[t, i]
			acc = 0
			for (i = 1; i <= neq; i++) acc = acc + (*A[i])[t, neq + 1]
			(*A[neq + 1])[t, neq + 1] = -acc
		}
	}

	tot2 = J(n, 1, 0)
	if (py) {
		for (j = 1; j <= neq; j++) {
			for (q = 1; q <= neq; q++) {
				tot2 = tot2 + bjk[j, q] :* P[, j] :* P[, q]
			}
		}
	}
	den1 = 1 :- 0.5 :* tot2
	den  = _easi_wm(den1, wt)

	// Corrected derivation (see tests/AUDIT.md).  Differentiating the cost
	// identity ln x = ln C(p, y, z, eps) at fixed x, with d lnC/d l_i = w_i
	// (Shephard's lemma) and d lnC/dy = Dn2, gives
	//     dy/dl_i = -w_i / Dn2        dy/dlnx = 1 / Dn2
	// with Dn2 = 1 + sum_j (C_j + D_j) l_j + 1/2 sum_jk b_jk l_j l_k summed
	// over ALL J goods.  The R package instead differentiates the y formula
	// holding the budget shares fixed, which leaves in terms that ought to
	// cancel and gets the denominator wrong; -legacy(elast)- reproduces that.
	Dn2 = J(n, 1, 1)
	if (!lgB) {
		for (i = 1; i <= neq + 1; i++) {
			Ci = Ypow * ((1::R) :* bjr[, i])
			Di = zy ? z * hjt[, i] : J(n, 1, 0)
			Dn2 = Dn2 + (Ci + Di) :* P[, i]
		}
		Dn2 = Dn2 + 0.5 :* rowsum((P * bjk) :* P)
	}

	// ---- expenditure (income) elasticities -----------------------------
	ER = J(1, neq + 1, 0); EI = J(1, neq + 1, 0)
	EPRICE = J(neq + 1, neq + 1, 0)
	EP     = J(neq + 1, neq + 1, 0)

	for (i = 1; i <= neq + 1; i++) {
		Ci = Ypow * ((1::R) :* bjr[, i])
		Di = zy ? z * hjt[, i] : J(n, 1, 0)
		Gi = py ? P * bjk[, i] : J(n, 1, 0)
		// ER is the semi-elasticity dw/dy of the reference paper: unchanged
		ER[i] = _easi_wm(Ci + Di + Gi, wt)
		if (lgB) EI[i] = 1 + ER[i] / ws[i]
		else     EI[i] = 1 + _easi_wm((Ci + Di + Gi) :/ Dn2, wt) / ws[i]

		AZ = Zc * (*A[i])			// n x (neq+1), one matmul
		Ui = rowsum(AZ :* P)
		mG   = _easi_wm(Gi, wt)
		mCDG = _easi_wm(Ci + Di + Gi, wt)

		// all prices at once: the per-household vectors enter only through
		// weighted means, so each is one vector-matrix product
		my  = _easi_wm(y, wt)
		mAZ = quadcross(wt, AZ) :/ quadsum(wt)		// 1 x (neq+1)
		if (!lgB) {
			mHS = -(quadcross(wt, ((Ci + Di + Gi) :/ Dn2) :* s)) :/ quadsum(wt)
		}
		else mU = _easi_wm(Ui, wt)
		for (q = 1; q <= neq + 1; q++) {
			if (lgB) {
				Bsc = -(ws[q] + mU) / den - my / den * mG
				Hm  = Bsc * mCDG + mAZ[q] + (py ? bjk[q, i] * my : 0)
			}
			else Hm = mHS[q] + mAZ[q] + (py ? bjk[q, i] * my : 0)
			EPRICE[q, i] = Hm / ws[i] - (i == q)
			EP[q, i]     = mAZ[q] + bjk[q, i] * my
		}
	}

	// The Slutsky matrix and the compensated elasticities are means over the
	// households, like the others: the compensated elasticity of household h
	// is (EP_h[q,i] + w_hq w_hi) / w_hq - (q == i), and its mean weighted by
	// the share w_hq carries the MEAN of the share products, not the product
	// of the mean shares.  The two differ by Cov(w_q, w_i) / wbar_q (0.03 to
	// 0.07 on the own-price elasticities of hixdata).  The mean of the
	// household Slutsky matrices is also negative semidefinite whenever each
	// of them is.  -legacy(elast)- keeps the R package's product of means.
	M2  = (lgB ? ws' * ws : quadcross(s, wt, s) :/ quadsum(wt))
	EPS = EP + M2 - diag(ws)
	EPQ = invsym(diag(ws)) * (EP + M2)

	// ---- demographic elasticities ---------------------------------------
	// R bug 1: the running sum `tempo4` is initialised outside the t loop, so
	//          it accumulates across demographics.
	// R bug 2: it indexes a[t, k, i] where the first dimension is
	//          (constant, z_1, ..., z_T) -- so t = 1 picks the CONSTANT block
	//          instead of demographic 1.  Correct index is t + 1.
	EZ = J(T, neq + 1, 0)
	for (i = 1; i <= neq + 1; i++) {
		cum = J(n, 1, 0)
		for (t = 1; t <= T; t++) {
			if (!lgEZ) cum = J(n, 1, 0)
			if ((py | pz | zy) & anyof(ipz, t)) {
				for (q = 1; q <= neq + 1; q++) {
					cum = cum + (*A[i])[lgEZ ? t : t + 1, q] :* P[, q]
				}
			}
			EZ[t, i] = _easi_wm(gjt[t, i] :+ hjt[t, i] :* y + cum, wt)
		}
	}

}

// ------------------------------------------------- analytic elasticity Jacobian
//
// The reported elasticities are written in terms of the COMPLETED coefficient
// arrays -- bjr, gjt, hjt, bjk, A[i] -- where the dropped good's entries are
// filled in by adding up.  Two of those completions are not symmetric: the
// free y*np entry b[q,i] also sits at (q,J), (J,q) and (J,J) of bjk, and the
// free price entry A[i][t,q] also sits at A[J][t,q], A[q][t,J] and A[J][t,J].
// Differentiating by hand in free coordinates is where the sign errors live.
//
// So we do not.  The completion is a linear map, b -> C b, with C a sparse
// 0/+-1 matrix built from the SAME rules as _easi_epoint.  Every elasticity is
// then differentiated in completed coordinates, where each formula is a plain
// weighted mean and the derivatives are one quadcross each, and the chain
// rule closes it: G = Gc * C.  One place for the adding-up, no hand-derived
// adjoint, and C can be checked by reproducing the completed arrays.
//
// The implicit utility y is held fixed here, as in the finite-difference
// Jacobian this replaced.  That is not the whole derivative: the
// generated-regressor term carries d y / d beta into the variance of the
// COEFFICIENTS, not into the elasticities, which depend on y household by
// household.  _easi_ejacy adds that part (easi 2.0.0).

// place the four adding-up images of price entry (i, tt, q) into column col
void _easi_cA(real matrix C, real scalar oA, real scalar na, real scalar Jg,
	real scalar i, real scalar tt, real scalar q, real scalar col)
{
	C[oA + ((i  - 1) * na + tt - 1) * Jg + q,  col] =  1
	C[oA + ((Jg - 1) * na + tt - 1) * Jg + q,  col] = -1
	C[oA + ((q  - 1) * na + tt - 1) * Jg + Jg, col] = -1
	C[oA + ((Jg - 1) * na + tt - 1) * Jg + Jg, col] =  1
}

// The completion map C (Kc x K).  Completed layout, column-major throughout:
//   bjr (r,i) at (i-1)R + r         gjt, hjt (t,i) at (i-1)T + t
//   bjk (q,i) at (i-1)Jg + q        A (i,tt,q) at ((i-1)na + tt-1)Jg + q
// with tt = 1 the np row and tt = s+1 the interaction with z[ipz[s]].
real matrix _easi_ecomp(real scalar neq, real scalar R, real scalar T,
	real scalar zy, real scalar py, real scalar pz, real rowvector ipz,
	real scalar k)
{
	real matrix C
	real scalar Jg, nipz, na, Kc, i, q, r, t, s, p
	real scalar onp, oynp, onpz, obr, ogt, oht, obk, oA

	Jg = neq + 1; nipz = cols(ipz); na = nipz + 1
	onp = 1 + R + T + (zy ? T : 0); oynp = onp + neq; onpz = oynp + (py ? neq : 0)
	obr = 0
	ogt = obr + R * Jg
	oht = ogt + T * Jg
	obk = oht + T * Jg
	oA  = obk + Jg * Jg
	Kc  = oA + na * Jg * Jg
	C = J(Kc, neq * k, 0)
	for (i = 1; i <= neq; i++) {
		p = (i - 1) * k
		for (r = 1; r <= R; r++) {
			C[obr + (i  - 1) * R + r, p + 1 + r] =  1
			C[obr + (Jg - 1) * R + r, p + 1 + r] = -1
		}
		for (t = 1; t <= T; t++) {
			C[ogt + (i  - 1) * T + t, p + 1 + R + t] =  1
			C[ogt + (Jg - 1) * T + t, p + 1 + R + t] = -1
			if (zy) {
				C[oht + (i  - 1) * T + t, p + 1 + R + T + t] =  1
				C[oht + (Jg - 1) * T + t, p + 1 + R + T + t] = -1
			}
		}
		for (q = 1; q <= neq; q++) {
			_easi_cA(C, oA, na, Jg, i, 1, q, p + onp + q)
			if (py) {
				C[obk + (i  - 1) * Jg + q,  p + oynp + q] =  1	// (q, i)
				C[obk + (Jg - 1) * Jg + q,  p + oynp + q] = -1	// (q, J)
				C[obk + (q  - 1) * Jg + Jg, p + oynp + q] = -1	// (J, q)
				C[obk + (Jg - 1) * Jg + Jg, p + oynp + q] =  1	// (J, J)
			}
			if (pz) {
				for (s = 1; s <= nipz; s++) {
					_easi_cA(C, oA, na, Jg, i, s + 1, q,
						p + onpz + (s - 1) * neq + q)
				}
			}
		}
	}
	return(C)
}

// Jacobian of the stacked elasticities (EI \ vec(EPRICE) \ vec(EZ) [\ vec(EPQ)])
// with respect to b, analytic.  Non-legacy formulas only; the legacy modes
// keep the finite-difference path.
real matrix _easi_ejac(real colvector b, real matrix s, real matrix P,
	real matrix z, real colvector y, real colvector wt, real colvector Dn2,
	real rowvector ws, real scalar R, real scalar py, real scalar zy,
	real scalar pz, real rowvector ipz, real scalar doCQ)
{
	real matrix C, Gc, Yr, Zp, Al, bjr, hjt, bjk, blk
	real colvector cb, u, v, uq
	real rowvector mZp, mP, gY, gZ, gP
	real scalar n, Jg, neq, T, nipz, na, k, K, Kc, nE
	real scalar onp, oynp, onpz, obr, ogt, oht, obk, oA
	real scalar i, q, r, t, sx, tt, sww, my, row, rEP, rEZ, rEQ

	n = rows(s); Jg = cols(s); neq = Jg - 1; T = cols(z)
	nipz = cols(ipz); na = nipz + 1
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	K    = neq * k
	obr = 0; ogt = obr + R * Jg; oht = ogt + T * Jg; obk = oht + T * Jg
	oA  = obk + Jg * Jg
	Kc  = oA + na * Jg * Jg

	C  = _easi_ecomp(neq, R, T, zy, py, pz, ipz, k)
	cb = C * b
	bjr = colshape(cb[| obr + 1 \ obr + R * Jg |], R)'	// R x Jg
	hjt = (T ? colshape(cb[| oht + 1 \ oht + T * Jg |], T)' : J(0, Jg, 0))
	bjk = colshape(cb[| obk + 1 \ obk + Jg * Jg |], Jg)'	// Jg x Jg, (q, i)

	Yr = J(n, R, 1)
	for (r = 2; r <= R; r++) Yr[, r] = r :* y:^(r - 1)
	Zp = J(n, 1, 1)
	if (nipz) Zp = Zp, z[, ipz]
	Al = Yr * bjr
	if (zy) Al = Al + z * hjt
	if (py) Al = Al + P * bjk
	sww = quadsum(wt)
	my  = quadcross(wt, y) / sww
	mZp = quadcross(wt, Zp) :/ sww
	mP  = quadcross(wt, P)  :/ sww
	gY  = quadcross(wt, Yr :/ Dn2) :/ sww
	gZ  = quadcross(wt, z  :/ Dn2) :/ sww
	gP  = quadcross(wt, P  :/ Dn2) :/ sww

	rEP = Jg
	rEZ = rEP + Jg * Jg
	rEQ = rEZ + T * Jg
	nE  = rEQ + (doCQ ? Jg * Jg : 0)
	Gc  = J(nE, Kc, 0)

	for (i = 1; i <= Jg; i++) {
		// ---- expenditure elasticity: EI_i = 1 + wm(Al_i / Dn2) / ws_i
		u = wt :* Al[, i] :/ (Dn2:^2)
		blk = -quadcross(Yr, u, P) :/ sww
		blk[, i] = blk[, i] + gY'
		Gc[| i, obr + 1 \ i, obr + R * Jg |] = vec(blk)' :/ ws[i]
		if (zy) {
			blk = -quadcross(z, u, P) :/ sww
			blk[, i] = blk[, i] + gZ'
			Gc[| i, oht + 1 \ i, oht + T * Jg |] = vec(blk)' :/ ws[i]
		}
		if (py) {
			blk = -0.5 :* quadcross(P, u, P) :/ sww
			blk[, i] = blk[, i] + gP'
			Gc[| i, obk + 1 \ i, obk + Jg * Jg |] = vec(blk)' :/ ws[i]
		}

		// ---- price elasticities: EPRICE[q,i] = Hm / ws_i - (i == q),
		//      Hm = -wm(s_q Al_i / Dn2) + wm(Zp A_i[., q]) + bjk[q,i] my
		for (q = 1; q <= Jg; q++) {
			row = rEP + (i - 1) * Jg + q
			v   = wt :* s[, q] :/ Dn2
			uq  = u :* s[, q]
			blk = quadcross(Yr, uq, P) :/ sww
			blk[, i] = blk[, i] - (quadcross(v, Yr) :/ sww)'
			Gc[| row, obr + 1 \ row, obr + R * Jg |] = vec(blk)' :/ ws[i]
			if (zy) {
				blk = quadcross(z, uq, P) :/ sww
				blk[, i] = blk[, i] - (quadcross(v, z) :/ sww)'
				Gc[| row, oht + 1 \ row, oht + T * Jg |] = vec(blk)' :/ ws[i]
			}
			if (py) {
				blk = 0.5 :* quadcross(P, uq, P) :/ sww
				blk[, i] = blk[, i] - (quadcross(v, P) :/ sww)'
				blk[q, i] = blk[q, i] + my
				Gc[| row, obk + 1 \ row, obk + Jg * Jg |] = vec(blk)' :/ ws[i]
			}
			for (tt = 1; tt <= na; tt++) {
				Gc[row, oA + ((i - 1) * na + tt - 1) * Jg + q] = mZp[tt] / ws[i]
			}
		}

		// ---- demographic elasticities: EZ[t,i] = wm(gjt + hjt y + sum_q A_i[s+1,q] P_q)
		for (t = 1; t <= T; t++) {
			row = rEZ + (i - 1) * T + t
			Gc[row, ogt + (i - 1) * T + t] = 1
			if (zy) Gc[row, oht + (i - 1) * T + t] = my
			if (pz & anyof(ipz, t)) {
				sx = selectindex(ipz :== t)[1]
				for (q = 1; q <= Jg; q++) {
					Gc[row, oA + ((i - 1) * na + sx) * Jg + q] = mP[q]
				}
			}
		}

		// ---- compensated quantity: EPQ[q,i] = (EP[q,i] + ws_q ws_i) / ws_q
		if (doCQ) {
			for (q = 1; q <= Jg; q++) {
				row = rEQ + (i - 1) * Jg + q
				for (tt = 1; tt <= na; tt++) {
					Gc[row, oA + ((i - 1) * na + tt - 1) * Jg + q] = mZp[tt] / ws[q]
				}
				if (py) Gc[row, obk + (i - 1) * Jg + q] = my / ws[q]
			}
		}
	}
	return(Gc * C)
}

// Sampling term of the influence function of the reported elasticities.
//
// Every reported elasticity is a ratio of two weighted means over the
// households, Nbar/Dbar (minus the identity for the price elasticities), or a
// weighted mean (the demographic semi-elasticities).  With sw_h = wt_h/sum(wt),
//     IF_h(Nbar/Dbar) = sw_h [ (N_h - Nbar)/Dbar - (Nbar/Dbar^2)(D_h - Dbar) ],
//     IF_h(Nbar)      = sw_h (N_h - Nbar).
// This is the variability of the households the means are taken over; the
// estimation of the coefficients adds IF(beta) times the Jacobian.  For good i
// and price q (stacked as _easi_estack does):
//     expenditure        N = a_hi = dw_hi/dln x         D = w_hi
//     uncompensated      N = dw_hi/dln p_q (y implicit) D = w_hi
//     demographic        N = dw_hi/dz_t                 (a mean)
//     compensated        N = EP_h[q,i] + w_hq w_hi      D = w_hq
// The numerators are those of _easi_epoint, rebuilt from the completed
// coefficients; est returns the ratios, which must equal the point estimates.
real colvector _easi_rif(real colvector N, real colvector D,
	real colvector sw, real scalar est)
{
	real scalar Nb, Db
	Nb  = quadcross(sw, N)
	Db  = quadcross(sw, D)
	est = Nb / Db
	return(sw :* ((N :- Nb) :/ Db :- (Nb / Db^2) :* (D :- Db)))
}

// ---------------------------------------------------------------------------
// pimpute(): the influence of the imputation of the prices on the analytic
// standard errors.  A household that gives its price to the mean m_gj of its
// group moves it by w_d (lp_dj - m_gj) / W_gj when its weight moves, and with
// it the prices it fills; its influence function gains that change times the
// derivative, in the filled prices of the group, of the estimating equations
// (the moments of the system, the scores of the probits) and of the sample
// terms of the elasticities.  The term sums to zero within each group: it
// cancels under vce(cluster) or vce(svy) at the level of the first grouping
// variable.  Checked against a brute force (tests/test_step28.do): analytic /
// brute-force standard errors 0.998-1.001 (without selection) and 0.998-1.000
// (with selection, where omitting the term overstated standard errors of
// elasticities by up to 9%).  The map and its application are those of
// equaids 1.2.0 (_eq_impmap, _eq_impapply):
// GI (n x J, the group of each filled price, 0 otherwise) and DT (one row
// per donor and group: good, row, group, w_d (lp_dj - m_gj) / W_gj).
// ---------------------------------------------------------------------------
// the rows of X summed by key (1..n), n x cols(X)
real matrix _easi_impsumby(real colvector key, real matrix X, real scalar n)
{
	real colvector o, ko
	real matrix    info, R
	R = J(n, cols(X), 0)
	if (rows(key) == 0) return(R)
	o  = order(key, 1)
	ko = key[o]
	info = panelsetup(ko, 1)
	R[ko[info[., 1]], .] = panelsum(X[o, .], info)
	return(R)
}

// the map of the imputation (LP0: the log prices before filling; GR: the
// groups of pimpute(), as integers); the largest gap between a filled price
// and the mean that should have filled it goes to _easi_imp_mapchk
void _easi_impmap(real matrix LP, real matrix LP0, real matrix GR,
	real colvector om, real matrix GI, real matrix DT)
{
	real colvector isI, isD, key, done, ix, id
	real rowvector vals
	real scalar    N, M, j, l, k, v, W1, m1, ng, chk
	N = rows(LP) ; M = cols(LP)
	GI = J(N, M, 0) ; DT = J(0, 4, .)
	ng = 0 ; chk = 0
	for (j = 1; j <= M; j++) {
		isI = (LP0[., j] :>= .)
		if (!sum(isI)) continue
		done = J(N, 1, 0)
		for (l = 1; l <= cols(GR); l++) {
			key = GR[., l]
			isD = (LP0[., j] :< .) :& (key :< .)
			if (!sum(isD)) continue
			vals = uniqrows(select(key, isD))'
			for (k = 1; k <= cols(vals); k++) {
				v  = vals[k]
				ix = selectindex(isI :& !done :& (key :== v))
				if (rows(ix) == 0) continue
				id = selectindex(isD :& (key :== v))
				W1 = sum(om[id])
				if (W1 <= 0) continue
				m1 = sum(om[id] :* LP0[id, j]) / W1
				chk = max((chk, max(abs(LP[ix, j] :- m1))))
				ng++
				GI[ix, j] = J(rows(ix), 1, ng)
				DT = DT \ (J(rows(id), 1, j), id, J(rows(id), 1, ng),
					om[id] :* (LP0[id, j] :- m1) :/ W1)
				done[ix] = J(rows(ix), 1, 1)
			}
		}
	}
	st_numscalar("_easi_imp_mapchk", chk)
}

// Dp[j]: n x q derivatives of the contributions in lp_j (NULL: none);
// returns the n x q influence of the imputation
real matrix _easi_impapply(pointer(real matrix) rowvector Dp, real scalar N,
	real scalar q, real matrix GI, real matrix DT)
{
	real matrix    X, SG
	real colvector key, r
	real scalar    j
	if (rows(DT) == 0) return(J(N, q, 0))
	X = J(0, q, .) ; key = J(0, 1, .)
	for (j = 1; j <= cols(Dp); j++) {
		if (Dp[j] == NULL) continue
		r = selectindex(GI[., j] :> 0)
		if (rows(r) == 0) continue
		X   = X \ (*Dp[j])[r, .]
		key = key \ GI[r, j]
	}
	SG = _easi_impsumby(key, X, max(GI))
	return(_easi_impsumby(DT[., 2], DT[., 4] :* SG[DT[., 3], .], N))
}

// the scores of the system (those of _easi_ifb) at the log prices P, the
// coefficients b fixed: y from the cost identity, the instrument from the
// coefficients b1 of the first pass (Pendakur's scheme), the first stage
// ZZiZX and Sigma fixed
real matrix _easi_impS(real colvector b, real colvector b1, real scalar nphase,
	real matrix Si, real matrix ZZiZX, real matrix s, real matrix P,
	real colvector lnx, real matrix z, real colvector wt, real rowvector ms,
	real scalar R, real scalar py, real scalar zy, real scalar pz,
	real rowvector ipz, real scalar lgQuant, real scalar k)
{
	real matrix np, X, Zi, Xhat, E
	real colvector ystone, pAp, pBp, y, ytil, yinst, pAp1, pBp1
	real scalar Jg, neq
	Jg = cols(s); neq = Jg - 1
	np = P[, 1::neq] :- P[, Jg]
	ystone = lnx - rowsum(s :* P)
	pAp = .; pBp = .
	_easi_quadf(b, np, z, R, py, pz, zy, ipz, lgQuant, pAp, pBp)
	y = (ystone :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)
	X = _easi_design(y, z, np, R, py, zy, pz, ipz)
	ytil = lnx - P * ms'
	if (nphase > 1) {
		pAp1 = .; pBp1 = .
		_easi_quadf(b1, np, z, R, py, pz, zy, ipz, lgQuant, pAp1, pBp1)
		yinst = (ytil :+ 0.5 :* pAp1) :/ (1 :- 0.5 :* pBp1)
	}
	else yinst = ytil
	Zi   = _easi_instr(yinst, z, np, R, py, zy, pz, ipz)
	Xhat = Zi * ZZiZX
	E    = s[, 1::neq] - X * rowshape(b', neq)'
	return(_easi_scores(Xhat, E, Si, wt, k, neq))
}

// the numerators of the elasticities, household by household (those of
// _easi_esamp), at the log prices P, the coefficients fixed
real matrix _easi_impNN(real colvector b, real matrix s, real matrix P,
	real matrix z, real colvector lnx, real colvector wt, real scalar R,
	real scalar py, real scalar zy, real scalar pz, real rowvector ipz,
	real scalar lgQuant, real scalar lgEZ, real scalar lgB, real scalar doCQ)
{
	real matrix np, EPRICE, EP, EPS, EPQ, EZ, NN, DD, Ts
	real rowvector EI, ER, ws
	real colvector ystone, pAp, pBp, y, Dn2, est
	real scalar Jg, neq
	Jg = cols(s); neq = Jg - 1
	np = P[, 1::neq] :- P[, Jg]
	ystone = lnx - rowsum(s :* P)
	pAp = .; pBp = .
	_easi_quadf(b, np, z, R, py, pz, zy, ipz, lgQuant, pAp, pBp)
	y = (ystone :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)
	EI = .; ER = .; EPRICE = .; EP = .; EPS = .; EPQ = .; EZ = .; ws = .; Dn2 = .
	_easi_epoint(b, s, P, z, y, wt, R, py, zy, pz, ipz, lgEZ, lgB, EI, ER,
		EPRICE, EP, EPS, EPQ, EZ, ws, Dn2)
	est = .; NN = .; DD = .
	Ts = _easi_esamp(b, s, P, z, y, wt, Dn2, R, py, zy, pz, ipz, doCQ, est,
		NN, DD)
	return(NN)
}


// pimpute(), with selection: the moments of the engine, household by
// household (those of _easi_g_run), at the log prices P, at fixed theta and
// alpha: Phi, phi and the Mills ratios from the probit designs at P, the
// completed y re-solved from y0, the instrument from the first pass b1
real matrix _easi_impgGh(real colvector th, real rowvector al, real colvector b1,
	real scalar nphase, real matrix P, real colvector lnx, real matrix z,
	real matrix Q, real matrix zm, real rowvector cix, real rowvector ao,
	real matrix s, real matrix D0, real colvector y0, real colvector wt,
	real rowvector ms, real rowvector ki, real rowvector ko,
	real rowvector bidx, real rowvector didx, real scalar R, real scalar py,
	real scalar zy, real scalar pz, real rowvector ipz, real scalar lgQuant)
{
	real matrix np, PH, PHI, M0, X, Z0, Xs, Zs, Gh
	real colvector ystone, y, ytil, yinst, pAp1, pBp1, b
	real rowvector dl, rg
	real scalar Jg, neq, i
	Jg = cols(s); neq = Jg - 1
	np = P[, 1::neq] :- P[, Jg]
	ystone = lnx - rowsum(s :* P)
	PH = .; PHI = .; M0 = .
	_easi_g_probs(al, P, lnx, z, Q, zm, cix, ao, PH, PHI, M0)
	b  = th[bidx]
	dl = J(1, cols(ki), 0)
	if (cols(didx)) dl[selectindex(cix :> 0)] = th[didx]'
	y  = _easi_g_y(y0, ystone, np, z, D0, M0, cix, b, dl, R, py, zy, pz, ipz,
		lgQuant)
	X  = _easi_design(y, z, np, R, py, zy, pz, ipz)
	ytil = lnx - P * ms'
	if (nphase > 1) {
		pAp1 = .; pBp1 = .
		_easi_quadf(b1, np, z, R, py, pz, zy, ipz, lgQuant, pAp1, pBp1)
		yinst = (ytil :+ 0.5 :* pAp1) :/ (1 :- 0.5 :* pBp1)
	}
	else yinst = ytil
	Z0 = _easi_instr(yinst, z, np, R, py, zy, pz, ipz)
	Xs = .; Zs = .
	_easi_g_blocks(X, Z0, PH, PHI, cix, ki, ko, Xs, Zs)
	Gh = J(rows(s), sum(ki), .)
	for (i = 1; i <= neq; i++) {
		rg = (ko[i] + 1)..(ko[i] + ki[i])
		Gh[., rg] = Zs[., rg] :* (wt :* (s[., i] - Xs[., rg] * th[rg]))
	}
	return(Gh)
}

// the scores of the probits of the corrected goods (w lambda s, the order
// of alpha) at the log prices P, alpha fixed
real matrix _easi_imppscore(real rowvector al, real matrix P, real colvector lnx,
	real matrix z, real matrix Q, real matrix zm, real rowvector cix,
	real rowvector ao, real matrix s, real colvector wt)
{
	real matrix Sc, S
	real colvector xb, a, r
	real scalar i, c
	Sc = J(rows(P), 0, .)
	for (i = 1; i <= cols(cix); i++) {
		if (!cix[i]) continue
		c  = cix[i]
		S  = _easi_seldesign(P, lnx, z, Q, (cols(Q) ? zm[i, .] : J(1, 0, .)))
		xb = S * al[| ao[c] + 1 \ ao[c] + cols(S) |]'
		a  = 2 :* (s[., i] :> 0) :- 1
		r  = a :* normalden(a :* xb) :/ normal(a :* xb)
		Sc = Sc, (wt :* r) :* S
	}
	return(Sc)
}

// the inverses of the observed Hessians of the probits, block-diagonal in
// the order of alpha (those of _easi_probit)
real matrix _easi_imppIi(real rowvector al, real matrix P, real colvector lnx,
	real matrix z, real matrix Q, real matrix zm, real rowvector cix,
	real rowvector ao, real matrix s, real colvector wt)
{
	real matrix B, S
	real colvector xb, a, r
	real rowvector rg
	real scalar i, c
	B = J(cols(al), cols(al), 0)
	for (i = 1; i <= cols(cix); i++) {
		if (!cix[i]) continue
		c  = cix[i]
		S  = _easi_seldesign(P, lnx, z, Q, (cols(Q) ? zm[i, .] : J(1, 0, .)))
		rg = (ao[c] + 1)..(ao[c] + cols(S))
		xb = S * al[rg]'
		a  = 2 :* (s[., i] :> 0) :- 1
		r  = a :* normalden(a :* xb) :/ normal(a :* xb)
		B[rg, rg] = invsym(quadcross(S, wt :* (r :* (r :+ xb)), S))
	}
	return(B)
}


real matrix _easi_esamp(real colvector b, real matrix s, real matrix P,
	real matrix z, real colvector y, real colvector wt, real colvector Dn2,
	real scalar R, real scalar py, real scalar zy, real scalar pz,
	real rowvector ipz, real scalar doCQ, real colvector est,
	| real matrix NN, real matrix DD)
{
	real matrix C, Yr, Zp, Al, alp, bjr, gjt, hjt, bjk, Ai, ZA, Ts
	real colvector cb, sw, Nh, one
	real scalar n, Jg, neq, T, nipz, na, k, nE, i, q, t, r, sx, e
	real scalar obr, ogt, oht, obk, oA, onp, oynp, onpz, rEP, rEZ, rEQ

	n = rows(s); Jg = cols(s); neq = Jg - 1; T = cols(z)
	nipz = cols(ipz); na = nipz + 1
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	obr = 0; ogt = obr + R * Jg; oht = ogt + T * Jg; obk = oht + T * Jg
	oA  = obk + Jg * Jg

	C   = _easi_ecomp(neq, R, T, zy, py, pz, ipz, k)
	cb  = C * b
	bjr = colshape(cb[| obr + 1 \ obr + R * Jg |], R)'	// R x Jg
	gjt = colshape(cb[| ogt + 1 \ ogt + T * Jg |], T)'	// T x Jg
	hjt = colshape(cb[| oht + 1 \ oht + T * Jg |], T)'	// T x Jg, 0 without zy
	bjk = colshape(cb[| obk + 1 \ obk + Jg * Jg |], Jg)'	// (q, i), 0 without py

	Yr = J(n, R, 1)
	for (r = 2; r <= R; r++) Yr[, r] = r :* y:^(r - 1)
	Zp = J(n, 1, 1)
	if (nipz) Zp = Zp, z[, ipz]
	Al = Yr * bjr
	if (zy) Al = Al + z * hjt
	if (py) Al = Al + P * bjk
	alp = Al :/ Dn2				// dw/dln x, household by household
	sw  = wt :/ quadsum(wt)
	one = J(n, 1, 1)

	rEP = Jg; rEZ = rEP + Jg * Jg; rEQ = rEZ + T * Jg
	nE  = rEQ + (doCQ ? Jg * Jg : 0)
	Ts  = J(n, nE, 0)
	est = J(nE, 1, .)
	if (args() > 14) {
		NN = J(n, nE, 0)
		DD = J(n, nE, 0)
	}

	for (i = 1; i <= Jg; i++) {
		Ts[, i] = _easi_rif(alp[, i], s[, i], sw, e)
		est[i]  = 1 + e
		if (args() > 14) {
			NN[, i] = alp[, i]
			DD[, i] = s[, i]
		}
		// A_i[tt, q]: tt = 1 the price term, tt = s + 1 its interaction
		// with the demographic ipz[s]
		Ai = colshape(cb[| oA + (i - 1) * na * Jg + 1 \ oA + i * na * Jg |], Jg)
		ZA = Zp * Ai				// n x Jg
		for (q = 1; q <= Jg; q++) {
			Nh = ZA[, q] + bjk[q, i] :* y - s[, q] :* alp[, i]
			Ts[, rEP + (i - 1) * Jg + q] = _easi_rif(Nh, s[, i], sw, e)
			est[rEP + (i - 1) * Jg + q]  = e - (i == q)
			if (args() > 14) {
				NN[, rEP + (i - 1) * Jg + q] = Nh
				DD[, rEP + (i - 1) * Jg + q] = s[, i]
			}
			if (doCQ) {
				Nh = ZA[, q] + bjk[q, i] :* y + s[, q] :* s[, i]
				Ts[, rEQ + (i - 1) * Jg + q] = _easi_rif(Nh, s[, q], sw, e)
				est[rEQ + (i - 1) * Jg + q]  = e
				if (args() > 14) {
					NN[, rEQ + (i - 1) * Jg + q] = Nh
					DD[, rEQ + (i - 1) * Jg + q] = s[, q]
				}
			}
		}
		for (t = 1; t <= T; t++) {
			Nh = gjt[t, i] :+ hjt[t, i] :* y
			if (pz & anyof(ipz, t)) {
				sx = selectindex(ipz :== t)[1]
				Nh = Nh + P * Ai[sx + 1, .]'
			}
			Ts[, rEZ + (i - 1) * T + t] = _easi_rif(Nh, one, sw, e)
			est[rEZ + (i - 1) * T + t]  = e
			if (args() > 14) {
				NN[, rEZ + (i - 1) * T + t] = Nh
				DD[, rEZ + (i - 1) * T + t] = one
			}
		}
	}
	return(Ts)
}

// The y term of the Jacobian of the elasticities.
//
// The implicit utility is a function of the coefficients,
//     y = (y_stone + p'A(z)p/2) / (1 - p'Bp/2),
// and the elasticities depend on y household by household: through
// a_hi = dw_hi/dln x (the powers of y in its numerator, and Dn2) and through
// the y x price terms.  _easi_ejac holds y fixed; this is the rest of the
// total derivative, element by element,
//     (1/Dbar) sum_h sw_h (dN_h/dy) (dy_h/dbeta),
// with, for good i and price q,
//     expenditure     dN/dy = d a_hi/dy,              D = w_hi
//     uncompensated   dN/dy = b_qi - w_hq d a_hi/dy,  D = w_hi
//     demographic     dN/dy = h_ti                    (a mean)
//     compensated     dN/dy = b_qi,                   D = w_hq
// and d a_hi/dy = (dAl_hi Dn2_h - Al_hi dDn2_h) / Dn2_h^2, where Al_hi is
// the numerator of a_hi and dDn2_h = sum_j p_hj dAl_hj.  In free
// coordinates, as dy/dbeta (_easi_dydb).  It weighs at most 0.2% of the
// standard errors on the examples; it is part of the derivative all the
// same (tests/test_step20.do checks the total against central differences
// with y recomputed from the coefficients).
real matrix _easi_ejacy(real colvector b, real matrix s, real matrix P,
	real matrix z, real colvector y, real colvector wt, real colvector Dn2,
	real scalar R, real scalar py, real scalar zy, real scalar pz,
	real rowvector ipz, real scalar doCQ)
{
	real matrix C, np, Q, bjr, hjt, bjk, Yr, Yd, Al, dAl, dal, Gy
	real colvector cb, sw, pAp, pBp, dDn2, v
	real rowvector ws, mQ
	real scalar n, Jg, neq, T, nipz, k, K, nE, i, q, t, r
	real scalar obr, ogt, oht, obk, onp, oynp, onpz, rEP, rEZ, rEQ

	n = rows(s); Jg = cols(s); neq = Jg - 1; T = cols(z); nipz = cols(ipz)
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	K    = neq * k
	obr = 0; ogt = obr + R * Jg; oht = ogt + T * Jg; obk = oht + T * Jg

	C   = _easi_ecomp(neq, R, T, zy, py, pz, ipz, k)
	cb  = C * b
	bjr = colshape(cb[| obr + 1 \ obr + R * Jg |], R)'	// R x Jg
	hjt = colshape(cb[| oht + 1 \ oht + T * Jg |], T)'	// T x Jg
	bjk = colshape(cb[| obk + 1 \ obk + Jg * Jg |], Jg)'	// (q, i)

	np = P[, 1::neq] :- P[, Jg]
	_easi_quadf(b, np, z, R, py, pz, zy, ipz, 0, pAp, pBp)
	Q  = _easi_dydb(np, z, y, pBp, R, py, pz, zy, ipz)	// n x K
	sw = wt :/ quadsum(wt)
	ws = quadcross(sw, s)
	mQ = quadcross(sw, Q)

	Yr = J(n, R, 1)
	Yd = J(n, R, 0)
	for (r = 2; r <= R; r++) {
		Yr[, r] = r :* y:^(r - 1)
		Yd[, r] = r * (r - 1) :* y:^(r - 2)
	}
	Al = Yr * bjr
	if (zy) Al = Al + z * hjt
	if (py) Al = Al + P * bjk
	dAl  = Yd * bjr
	dDn2 = rowsum(dAl :* P)
	dal  = (dAl :* Dn2 :- Al :* dDn2) :/ (Dn2:^2)

	rEP = Jg; rEZ = rEP + Jg * Jg; rEQ = rEZ + T * Jg
	nE  = rEQ + (doCQ ? Jg * Jg : 0)
	Gy  = J(nE, K, 0)
	for (i = 1; i <= Jg; i++) {
		v = sw :* dal[, i]
		Gy[i, .] = quadcross(v, Q) :/ ws[i]
		for (q = 1; q <= Jg; q++) {
			Gy[rEP + (i - 1) * Jg + q, .] =
				(bjk[q, i] :* mQ - quadcross(v :* s[, q], Q)) :/ ws[i]
			if (doCQ) {
				Gy[rEQ + (i - 1) * Jg + q, .] = bjk[q, i] :* mQ :/ ws[q]
			}
		}
		for (t = 1; t <= T; t++) Gy[rEZ + (i - 1) * T + t, .] = hjt[t, i] :* mQ
	}
	return(Gy)
}

// ---------------------------------------------------------------- selection
//
// Selection of the buyers (option selection): Shonkwiler and Yen (1999).
// For each corrected good i, a probit of buying on s = (1, ln p_1 - ln x,
// ..., ln p_J - ln x, z, the probit-only variables of the good): homogeneous
// of degree zero in prices and expenditure.  The probit engine is that of
// duvm and equaids.

// weighted probit of dd (0/1) on S by Fisher scoring (concave log
// likelihood); Phi and phi at the estimate, IF the household influence
// functions of the coefficients (w r s') I^-1.  ok = 0 when everybody or
// nobody buys.
real colvector _easi_probit(real colvector dd, real matrix S, real colvector w,
	real colvector Ph, real colvector ph, real matrix IF, real scalar ok)
{
	real colvector g, xb, P, f, r, a, step
	real matrix Ii
	real scalar it, n, n1

	n  = rows(dd)
	g  = J(cols(S), 1, 0)
	n1 = sum(dd)
	ok = (n1 > 0 & n1 < n)
	if (!ok) {
		Ph = J(n, 1, 1); ph = J(n, 1, 0); IF = J(n, cols(S), 0)
		return(g)
	}
	for (it = 1; it <= 200; it++) {
		xb = S * g
		P  = rowmin((rowmax((normal(xb), J(n, 1, 1e-15))), J(n, 1, 1 - 1e-15)))
		f  = normalden(xb)
		r  = (dd - P) :* f :/ (P :* (1 :- P))
		a  = (f:^2) :/ (P :* (1 :- P))
		Ii = invsym(quadcross(S, w :* a, S))
		step = Ii * quadcross(S, w :* r)
		g  = g + step
		if (max(abs(step)) < 1e-11) break
	}
	xb = S * g
	P  = rowmin((rowmax((normal(xb), J(n, 1, 1e-15))), J(n, 1, 1 - 1e-15)))
	f  = normalden(xb)
	r  = (dd - P) :* f :/ (P :* (1 :- P))
	Ph = normal(xb)
	ph = f
	// The influence function of the estimator uses the OBSERVED Hessian,
	// sum w lambda (lambda + xb) s s' with lambda = q phi(q xb)/Phi(q xb),
	// q = 2d - 1 (Greene): the derivative of the score the estimator
	// actually solves.  The expected (Fisher) information that the scoring
	// steps use is equal to it only in expectation under a correct probit;
	// on survey data the two differ at first order (a brute-force influence
	// function measured the difference at up to 4.6% of the standard
	// errors of the system).
	a  = (2 :* dd :- 1)
	r  = a :* normalden(a :* xb) :/ normal(a :* xb)
	Ii = invsym(quadcross(S, w :* (r :* (r :+ xb)), S))
	IF = ((w :* r) :* S) * Ii
	return(g)
}

// the probit design of good i: (1, ln p - ln x, z, q[mask_i])
real matrix _easi_seldesign(real matrix LP, real colvector LX, real matrix Z,
	real matrix Q, real rowvector mask)
{
	real matrix S
	S = (J(rows(LP), 1, 1), LP :- LX, Z)
	if (cols(Q)) {
		if (any(mask)) S = S, Q[., selectindex(mask)]
	}
	return(S)
}

// the probits of the corrected goods, before the system: Phi, phi and the
// influence functions of the coefficients stored in the variables phv, phiv,
// ifv (one set per corrected good, in the order of the shares); the
// coefficients (rows: _cons, the prices, the demographics, the probit-only
// variables; columns: the J-1 goods, missing when not corrected) and the
// diagnostics (percent of buyers, McFadden pseudo-R2, households predicted
// with probability 0 or 1) in Stata matrices.
void _easi_selprobits(string scalar svars, string scalar lpvars,
	string scalar lxvar, string scalar zvars, string scalar qvars,
	string scalar touse, string scalar wvar, string scalar sgstr,
	string scalar zmstr, string scalar phv, string scalar phiv,
	string scalar ifv, string scalar alphan, string scalar diagn,
	string scalar alsn)
{
	real matrix W, LP, Z, Q, zm, S, IFg, AL, DG, PH, PHI, IFA
	real colvector LX, om, dd, Ph, ph, g, xb, pr
	real rowvector sg, als
	real scalar N, M, nq, pb, pz, i, okp, ll, ll0, pbar

	W  = st_data(., svars, touse)
	LP = st_data(., lpvars, touse)
	LX = st_data(., lxvar, touse)
	Z  = st_data(., zvars, touse)
	Q  = (qvars == "" ? J(rows(W), 0, .) : st_data(., qvars, touse))
	om = st_data(., wvar, touse)
	N = rows(W); M = cols(W); nq = cols(Q)
	sg = strtoreal(tokens(sgstr))
	zm = (nq ? rowshape(strtoreal(tokens(zmstr)), M - 1) : J(M - 1, 0, .))
	pb = 1 + M + cols(Z); pz = pb + nq
	AL = J(pz, M - 1, .); DG = J(M - 1, 3, .)
	PH = J(N, 0, .); PHI = J(N, 0, .); IFA = J(N, 0, .); als = J(1, 0, .)
	for (i = 1; i < M; i++) {
		if (!sg[i]) continue
		S  = _easi_seldesign(LP, LX, Z, Q, (nq ? zm[i, .] : J(1, 0, .)))
		dd = W[., i] :> 0
		g  = _easi_probit(dd, S, om, Ph, ph, IFg, okp)
		if (!okp) {
			errprintf("easi: selection: every household (or none) buys good %g\n", i)
			exit(459)
		}
		PH = PH, Ph; PHI = PHI, ph; IFA = IFA, IFg; als = als, g'
		if (nq) {
			if (any(zm[i, .])) AL[(1::pb) \ (pb :+ selectindex(zm[i, .])'), i] = g
			else AL[1::pb, i] = g
		}
		else AL[1::pb, i] = g
		pbar = quadsum(om :* dd) / quadsum(om)
		xb  = S * g
		ll  = quadsum(om :* (dd :* lnnormal(xb) + (1 :- dd) :* lnnormal(-xb)))
		ll0 = quadsum(om) * (pbar * ln(pbar) + (1 - pbar) * ln(1 - pbar))
		pr  = normal(xb)
		DG[i, 1] = 100 * pbar
		DG[i, 2] = 1 - ll / ll0
		DG[i, 3] = sum((pr :< 1e-10) :| (pr :> 1 - 1e-10))
	}
	if (phv != "") {
		st_store(., tokens(phv), touse, PH)
		st_store(., tokens(phiv), touse, PHI)
		st_store(., tokens(ifv), touse, IFA)
	}
	st_matrix(alphan, AL)
	st_matrix(diagn, DG)
	st_matrix(alsn, als)
}

// easidiag, selection: VIF(delta) of each corrected good at the starting
// point y = Stone index -- [(X'WX)^-1]_dd (phi'W phi) with X = (Phi_i x(y),
// phi_i), the regressors of its equation in the estimation (_easi_g_blocks),
// Phi and phi from the probits (stacked coefficients in alsn).  Missing for
// the goods not corrected.
void _easi_d_selvif(string scalar svars, string scalar lpvars,
	string scalar lxvar, string scalar zvars, string scalar qvars,
	string scalar touse, string scalar wvar, string scalar sgstr,
	string scalar zmstr, string scalar alsn, real scalar R, real scalar py,
	real scalar pz, real scalar zy, string scalar sipz)
{
	real matrix s, p, z, Q, zm, np, X, PH, PHI, M0, Xi
	real colvector lnx, wt, ystone
	real rowvector cix, ao, al, ipz, vif
	real scalar n, J, neq, C, i, c, nq, kk

	s   = st_data(., svars,  touse)
	p   = st_data(., lpvars, touse)
	lnx = st_data(., lxvar,  touse)
	z   = st_data(., zvars,  touse)
	wt  = st_data(., wvar,   touse)
	Q   = (qvars == "" ? J(rows(s), 0, .) : st_data(., qvars, touse))
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))
	n = rows(s); J = cols(s); neq = J - 1; nq = cols(Q)
	np     = p[, 1::neq] :- p[, J]
	ystone = lnx - rowsum(s :* p)
	zm = (nq ? rowshape(strtoreal(tokens(zmstr)), neq) : J(neq, 0, .))
	cix = strtoreal(tokens(sgstr))
	C = 0
	for (i = 1; i <= neq; i++) {
		if (cix[i]) {
			C++
			cix[i] = C
		}
	}
	ao = J(1, C, 0)
	for (i = 1; i <= neq; i++) {
		if (!cix[i]) continue
		c = cix[i]
		if (c < C) ao[c + 1] = ao[c] + 1 + J + cols(z) + (nq ? sum(zm[i, .]) : 0)
	}
	al = st_matrix(alsn)
	_easi_g_probs(al, p, lnx, z, Q, zm, cix, ao, PH, PHI, M0)
	X = _easi_design(ystone, z, np, R, py, zy, pz, ipz)
	vif = J(1, neq, .)
	for (i = 1; i <= neq; i++) {
		if (!cix[i]) continue
		c  = cix[i]
		Xi = (PH[., c] :* X, PHI[., c])
		kk = cols(Xi)
		vif[i] = invsym(quadcross(Xi, wt, Xi))[kk, kk] *
			quadcross(PHI[., c], wt, PHI[., c])
	}
	st_matrix("_d_selvif", vif)
}

// ===========================================================================
// The general engine: every equation with its own regressors and
// instruments.  Built for the selection of the buyers (Shonkwiler and Yen
// 1999): the regressors of a corrected good i are (Phi_i x(y), phi_i), its
// instruments (Phi_i z, phi_i); the other equations keep x(y) and z.
//
// It is the GMM form of the restricted 3SLS: the moments g_i = Z_i' W e_i,
// stacked, weighted by Omega^-1 with Omega_ij = sigma_ij Z_i' W Z_j, under
// the symmetry restrictions (bordered).  With H = M Omega^-1 M',
// M = diag(X_i' W Z_i), the step is H theta = M Omega^-1 (Z_i' W w_i).  With
// the same X and Z in every equation, H = Sigma^-1 # A and this IS
// _easi_3sls.  The "second stage on the fitted X-hat" form would not be
// convergent here: the residual of the projection on z is not orthogonal to
// Phi_i z.
//
// Variance: theta-hat - theta = P L (sum_h g_h + D_a (alpha-hat - alpha)),
// L = M Omega^-1, P the restricted inverse of L D_t, D_t = -d sum g/d theta
// with y recomputed from theta (the generated regressor), D_a = d sum g/d alpha
// (the probits: Phi, phi and y move with alpha).  The two Jacobians are taken
// by central differences.
// ===========================================================================

// the probit indices, Phi, phi and the inverse Mills ratio of the non-buyers
// phi/(1 - Phi) of the corrected goods, at alpha (stacked, good by good)
void _easi_g_probs(real rowvector al, real matrix P, real colvector lnx,
	real matrix z, real matrix Q, real matrix zm, real rowvector cix,
	real rowvector ao, real matrix PH, real matrix PHI, real matrix M0,
	| real matrix XB)
{
	real matrix S
	real colvector xb
	real scalar i, c, n

	n = rows(P)
	PH = J(n, 0, .); PHI = J(n, 0, .); M0 = J(n, 0, .); XB = J(n, 0, .)
	for (i = 1; i <= cols(cix); i++) {
		if (!cix[i]) continue
		c  = cix[i]
		S  = _easi_seldesign(P, lnx, z, Q, (cols(Q) ? zm[i, .] : J(1, 0, .)))
		xb = S * al[| ao[c] + 1 \ ao[c] + cols(S) |]'
		PH  = PH, normal(xb)
		PHI = PHI, normalden(xb)
		M0  = M0, normalden(xb) :/ normal(-xb)
		XB  = XB, xb
	}
}

// the implicit utility y (decision DC1): the Stone index of the LATENT
// shares -- those of the non-buyers of a corrected good replaced by
// E[w*_i | not bought] = f_i - delta_i phi_i / (1 - Phi_i), the last good
// absorbing so that they still sum to one, hence the relative price np_i --
// a fixed point in y.  Without a corrected good it is the usual y.
real colvector _easi_g_y(real colvector y0, real colvector ystone,
	real matrix np, real matrix z, real matrix D0, real matrix M0,
	real rowvector cix, real colvector b, real rowvector dl, real scalar R,
	real scalar py, real scalar zy, real scalar pz, real rowvector ipz,
	real scalar lgQuant)
{
	real matrix B, F
	real colvector pAp, pBp, y, yold, comp
	real scalar neq, i, it

	neq = cols(np)
	_easi_quadf(b, np, z, R, py, pz, zy, ipz, lgQuant, pAp, pBp)
	if (!any(cix)) return((ystone :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp))
	B = rowshape(b', neq)'
	y = y0
	for (it = 1; it <= 500; it++) {
		yold = y
		F    = _easi_design(y, z, np, R, py, zy, pz, ipz) * B
		comp = J(rows(y), 1, 0)
		for (i = 1; i <= neq; i++) {
			if (!cix[i]) continue
			comp = comp + D0[, i] :* np[, i] :* (F[, i] - dl[i] :* M0[, cix[i]])
		}
		y = (ystone - comp :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)
		if (max(abs(y - yold)) < 1e-14) break
	}
	return(y)
}

// the regressors and the instruments of every equation, side by side
void _easi_g_blocks(real matrix X, real matrix Z0, real matrix PH,
	real matrix PHI, real rowvector cix, real rowvector ki, real rowvector ko,
	real matrix Xs, real matrix Zs)
{
	real scalar i, k, c, n

	n = rows(X); k = cols(X)
	Xs = J(n, sum(ki), .); Zs = J(n, sum(ki), .)
	for (i = 1; i <= cols(ki); i++) {
		if (cix[i]) {
			c = cix[i]
			Xs[| 1, ko[i] + 1 \ n, ko[i] + k |] = PH[, c] :* X
			Zs[| 1, ko[i] + 1 \ n, ko[i] + k |] = PH[, c] :* Z0
			Xs[., ko[i] + k + 1] = PHI[, c]
			Zs[., ko[i] + k + 1] = PHI[, c]
		}
		else {
			Xs[| 1, ko[i] + 1 \ n, ko[i] + k |] = X
			Zs[| 1, ko[i] + 1 \ n, ko[i] + k |] = Z0
		}
	}
}

// the restricted GMM step (see above); Sigma iterated maxsig times
void _easi_g_solve(real matrix Xs, real matrix Zs, real matrix Y,
	real matrix Rt, real colvector wt, real scalar neff, real rowvector ki,
	real rowvector ko, real rowvector eqi, real scalar maxsig,
	real colvector th, real matrix Sigma, real matrix L, real matrix H)
{
	real matrix ZZ, M, Om, Oi, E
	real colvector zy, sol, thold
	real rowvector rg
	real scalar neq, Kt, r, i, j, it

	neq = cols(Y); Kt = sum(ki); r = rows(Rt)
	ZZ = quadcross(Zs, wt, Zs)
	M  = J(Kt, Kt, 0); zy = J(Kt, 1, 0); Om = J(Kt, Kt, 0)
	for (i = 1; i <= neq; i++) {
		rg = (ko[i] + 1)..(ko[i] + ki[i])
		M[rg, rg]  = quadcross(Xs[., rg], wt, Zs[., rg])
		zy[rg]     = quadcross(Zs[., rg], wt, Y[., i])
		Om[rg, rg] = ZZ[rg, rg]
	}
	// step 1: Sigma = I, the restricted 2SLS
	Oi  = invsym(Om)
	L   = M * Oi
	H   = L * M'
	sol = _easi_ksolve(_easi_border(H, Rt), (L * zy) \ J(r, 1, 0))
	th  = sol[1::Kt]
	E   = J(rows(Y), neq, .)
	for (it = 1; it <= maxsig; it++) {
		thold = th
		for (i = 1; i <= neq; i++) {
			rg = (ko[i] + 1)..(ko[i] + ki[i])
			E[., i] = Y[., i] - Xs[., rg] * th[rg]
		}
		Sigma = J(neq, neq, .)
		for (i = 1; i <= neq; i++) {
			for (j = 1; j <= neq; j++) {
				Sigma[i, j] = quadcross(E[., i], wt, E[., j]) /
					sqrt((neff - ki[i]) * (neff - ki[j]))
			}
		}
		Om  = Sigma[eqi, eqi] :* ZZ
		Oi  = invsym(Om)
		L   = M * Oi
		H   = L * M'
		sol = _easi_ksolve(_easi_border(H, Rt), (L * zy) \ J(r, 1, 0))
		th  = sol[1::Kt]
		if (it > 1 & max(abs(th - thold)) < 1e-11) break
	}
}

// the moments summed over the households, at theta (y recomputed from theta
// unless ycond) and at the probits given by PH, PHI, M0
real colvector _easi_g_gsum(real colvector th, real colvector yc,
	real colvector ystone, real matrix np, real matrix z, real matrix D0,
	real matrix PH, real matrix PHI, real matrix M0, real matrix Z0,
	real matrix Y, real colvector wt, real rowvector cix, real rowvector ki,
	real rowvector ko, real rowvector bidx, real rowvector didx, real scalar R,
	real scalar py, real scalar zy, real scalar pz, real rowvector ipz,
	real scalar lgQuant, real scalar ycond)
{
	real matrix X, Xs, Zs
	real colvector y, g
	real rowvector dl, rg
	real scalar i

	dl = J(1, cols(ki), 0)
	if (cols(didx)) dl[selectindex(cix :> 0)] = th[didx]'
	y  = (ycond ? yc : _easi_g_y(yc, ystone, np, z, D0, M0, cix, th[bidx],
		dl, R, py, zy, pz, ipz, lgQuant))
	X  = _easi_design(y, z, np, R, py, zy, pz, ipz)
	_easi_g_blocks(X, Z0, PH, PHI, cix, ki, ko, Xs, Zs)
	g  = J(sum(ki), 1, 0)
	for (i = 1; i <= cols(ki); i++) {
		rg = (ko[i] + 1)..(ko[i] + ki[i])
		g[rg] = quadcross(Zs[., rg], wt, Y[., i] - Xs[., rg] * th[rg])
	}
	return(g)
}

void _easi_g_run(string scalar svars, string scalar lpvars, string scalar lxvar,
	string scalar zvars, string scalar touse, string scalar wvar,
	real scalar neff, real scalar vtype, string scalar clvar,
	string scalar strvar, string scalar psuvar, string scalar fpcvar,
	real scalar R,
	real scalar py, real scalar pz, real scalar zy, string scalar sipz,
	real scalar lgInstr, real scalar lgQuant, real scalar lgEZ, real scalar lgSE,
	real scalar lgB, real scalar lgESE, real scalar doESE, real scalar doCQ,
	real scalar lgCondSE, real scalar maxsig,
	real scalar tol, real scalar maxit, real scalar noisy,
	string scalar bnm, string scalar Vnm, string scalar Snm,
	string scalar itnm, string scalar crnm, string scalar cvnm,
	real scalar doMkt, string scalar sgstr, string scalar qvars,
	string scalar zmstr, string scalar alnm, string scalar ifv,
	string scalar lp0vars, string scalar impgvars, real scalar doImp)
{
	real matrix s, p, z, np, X, Z0, Y, Rr, Rt, V, Sigma, L, H, Xs, Zs
	real matrix GI, DT, KDg, KDsc, BIi, Pup, Pdn, Gh0
	real matrix Nu, Du, Nd, Dd, Alu, Ald, HSu, HSd
	real rowvector offu, offd
	real colvector b1, rr, yu, yd
	pointer(real matrix) rowvector DGh, DSc, DNg, DDg
	real scalar jj, eps
	real matrix Q, zm, PH, PHI, M0, D0, IFA, Gh, Dt, Da, Gam, Pm, IFt, E
	real matrix PH2, PHI2, M02, Vb
	real colvector lnx, ystone, ytil, y, yold, yinst, th, thold, pAp, pBp
	real colvector wt, clid, strid, fpcv, b, gp, gm, thp, thm, yz, tch
	real rowvector ipz, ms, cix, ki, ko, eqi, bidx, didx, ao, al, alp, alm, dl, rg, vif
	real scalar n, J, neq, T, k, K, Kt, C, i, c, it, crit, conv, nq
	real scalar onp, oynp, onpz, nipz, interact, ph, nphase, j, h, na

	s   = st_data(., svars,  touse)
	p   = st_data(., lpvars, touse)
	lnx = st_data(., lxvar,  touse)
	z   = st_data(., zvars,  touse)
	wt  = st_data(., wvar,   touse)
	clid  = (clvar  == "" ? J(0, 1, 0) : st_data(., clvar,  touse))
	strid = (strvar == "" ? J(0, 1, 0) : st_data(., strvar, touse))
	fpcv  = (fpcvar == "" ? J(0, 1, 0) : st_data(., fpcvar, touse))
	if (psuvar != "") clid = st_data(., psuvar, touse)
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))

	n = rows(s); J = cols(s); neq = J - 1; T = cols(z); nipz = cols(ipz)
	np     = p[, 1::neq] :- p[, J]
	ystone = lnx - rowsum(s :* p)
	ms     = quadcross(wt, s) :/ quadsum(wt)
	ytil   = lnx - p * ms'
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	K    = neq * k
	Rr   = _easi_restrict(neq, k, onp, py, oynp, pz, onpz, nipz)
	Y    = s[, 1::neq]
	interact = (py | pz | zy)

	// the corrected goods, the layout of theta = (b_1, delta_1, b_2, ...)
	cix = J(1, neq, 0)
	if (sgstr != "") cix = strtoreal(tokens(sgstr))
	C = 0
	for (i = 1; i <= neq; i++) {
		if (cix[i]) {
			C++
			cix[i] = C
		}
	}
	ki = J(1, neq, k) + (cix :> 0)
	ko = J(1, neq, 0)
	for (i = 2; i <= neq; i++) ko[i] = ko[i - 1] + ki[i - 1]
	Kt = sum(ki)
	eqi = J(1, Kt, .); bidx = J(1, 0, .); didx = J(1, 0, .)
	for (i = 1; i <= neq; i++) {
		eqi[(ko[i] + 1)..(ko[i] + ki[i])] = J(1, ki[i], i)
		bidx = bidx, (ko[i] + 1)..(ko[i] + k)
		if (cix[i]) didx = didx, ko[i] + k + 1
	}
	Rt = J(rows(Rr), Kt, 0)
	if (rows(Rr)) Rt[., bidx] = Rr

	// the probits of the corrected goods
	Q  = (qvars == "" ? J(n, 0, .) : st_data(., qvars, touse))
	nq = cols(Q)
	zm = (nq ? rowshape(strtoreal(tokens(zmstr)), neq) : J(neq, 0, .))
	ao = J(1, C, 0); al = J(1, 0, .); IFA = J(n, 0, .)
	D0 = J(n, neq, 0)
	if (C) {
		al  = st_matrix(alnm)
		IFA = st_data(., ifv, touse)
		for (i = 1; i <= neq; i++) {
			if (!cix[i]) continue
			c = cix[i]
			if (c < C) ao[c + 1] = ao[c] + 1 + J + T + (nq ? sum(zm[i, .]) : 0)
			D0[., i] = (s[., i] :<= 0)
		}
	}
	_easi_g_probs(al, p, lnx, z, Q, zm, cix, ao, PH, PHI, M0)

	nphase = lgInstr ? 1 : 2
	yinst  = ytil
	pAp = J(n, 1, 0); pBp = J(n, 1, 0)
	if (noisy) printf("\n{txt}iterating on the implicit utility index y...\n")
	for (ph = 1; ph <= nphase; ph++) {
		y  = ystone
		th = J(Kt, 1, 0)
		crit = 1; conv = 0; it = 0
		Z0 = _easi_instr(yinst, z, np, R, py, zy, pz, ipz)
		while (it < maxit) {
			it++
			thold = th
			yold  = y
			X = _easi_design(y, z, np, R, py, zy, pz, ipz)
			_easi_g_blocks(X, Z0, PH, PHI, cix, ki, ko, Xs, Zs)
			_easi_g_solve(Xs, Zs, Y, Rt, wt, neff, ki, ko, eqi, maxsig, th, Sigma, L, H)
			b  = th[bidx]
			dl = J(1, neq, 0)
			if (C) dl[selectindex(cix :> 0)] = th[didx]'
			_easi_quadf(b, np, z, R, py, pz, zy, ipz, lgQuant, pAp, pBp)
			y = _easi_g_y(y, ystone, np, z, D0, M0, cix, b, dl, R, py, zy,
				pz, ipz, lgQuant)
			if (interact) {
				if (it > 1) crit = sum((th - thold):^2)
			}
			else crit = max(abs(y - yold))
			if (noisy) printf("{txt}  iteration %3.0f    criterion = %12.5e\n", it, crit)
			if (crit <= tol) {
				conv = 1
				break
			}
		}
		if (ph == 1) b1 = b
		if (ph < nphase) yinst = (ytil :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)
	}

	// final refit at the converged y
	X = _easi_design(y, z, np, R, py, zy, pz, ipz)
	_easi_g_blocks(X, Z0, PH, PHI, cix, ki, ko, Xs, Zs)
	_easi_g_solve(Xs, Zs, Y, Rt, wt, neff, ki, ko, eqi, maxsig, th, Sigma, L, H)

	// the identification of delta: its variance inflation due to the
	// collinearity of phi_i with the other regressors of its equation,
	// Phi_i x -- [(X'WX)^-1]_dd (phi'W phi); above 10, delta rests on
	// little more than the curvature of the probit
	if (C) {
		vif = J(1, neq, .)
		for (i = 1; i <= neq; i++) {
			if (!cix[i]) continue
			rg = (ko[i] + 1)..(ko[i] + ki[i])
			vif[i] = invsym(quadcross(Xs[., rg], wt, Xs[., rg]))[ki[i], ki[i]] *
				quadcross(Xs[., ko[i] + ki[i]], wt, Xs[., ko[i] + ki[i]])
		}
		st_matrix("_easi_selvif", vif)

		// the completion of y for the non-buyers (decision DC1): the
		// households with a zero share of a corrected good (weighted
		// percent), the mean absolute change of their y against the y of
		// the observed zeros, and the standard deviation of y
		yz  = _easi_g_y(y, ystone, np, z, J(n, neq, 0), M0, cix, b, dl, R,
			py, zy, pz, ipz, lgQuant)
		tch = rowsum(D0) :> 0
		st_matrix("_easi_ycomp", (100 * quadsum(wt :* tch) / quadsum(wt),
			(quadsum(wt :* tch) > 0 ?
			quadsum(wt :* tch :* abs(y - yz)) / quadsum(wt :* tch) : 0),
			sqrt(quadvariance(y, wt))))
	}

	// ---- variance: the influence functions of theta ----------------------
	E = J(n, neq, .)
	Gh = J(n, Kt, .)
	for (i = 1; i <= neq; i++) {
		rg = (ko[i] + 1)..(ko[i] + ki[i])
		E[., i] = Y[., i] - Xs[., rg] * th[rg]
		Gh[., rg] = Zs[., rg] :* (wt :* E[., i])
	}

	// pimpute(): the influence of the imputation of the prices -- the
	// moments of the system (Phi, phi, the completed y and the instrument
	// move with the filled prices, at fixed theta, alpha, Sigma and first
	// pass), the scores of the probits (through their observed Hessian),
	// and the numerators and denominators of the elasticities; central
	// differences, good by good, every filled price of the good at once
	DNg = J(1, J, NULL); DDg = J(1, J, NULL); GI = J(0, 0, .); DT = J(0, 4, .)
	KDg = J(n, Kt, 0)
	if (doImp & vtype) {
		_easi_impmap(p, st_data(., lp0vars, touse), st_data(., impgvars, touse),
			wt, GI, DT)
	}
	if (rows(DT)) {
		Gh0 = _easi_impgGh(th, al, b1, nphase, p, lnx, z, Q, zm, cix, ao, s, D0, y,
			wt, ms, ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant)
		st_numscalar("_easi_imp_gchk", mreldif(Gh0, Gh))
		DGh = J(1, J, NULL); DSc = J(1, J, NULL)
		eps = 1e-5
		for (jj = 1; jj <= J; jj++) {
			rr = selectindex(GI[., jj] :> 0)
			if (rows(rr) == 0) continue
			Pup = p; Pup[rr, jj] = Pup[rr, jj] :+ eps
			Pdn = p; Pdn[rr, jj] = Pdn[rr, jj] :- eps
			DGh[jj] = &((_easi_impgGh(th, al, b1, nphase, Pup, lnx, z, Q, zm, cix,
				ao, s, D0, y, wt, ms, ki, ko, bidx, didx, R, py, zy, pz, ipz,
				lgQuant) - _easi_impgGh(th, al, b1, nphase, Pdn, lnx, z, Q, zm,
				cix, ao, s, D0, y, wt, ms, ki, ko, bidx, didx, R, py, zy, pz,
				ipz, lgQuant)) :/ (2 * eps))
			if (C) DSc[jj] = &((_easi_imppscore(al, Pup, lnx, z, Q, zm, cix, ao, s,
				wt) - _easi_imppscore(al, Pdn, lnx, z, Q, zm, cix, ao, s, wt)) :/
				(2 * eps))
			if (C) {
				_easi_g_parts(th, al, y, s, Pup, z, lnx, Q, zm, cix, ao,
					lnx - rowsum(s :* Pup), Pup[, 1::neq] :- Pup[, J], D0, ki,
					ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, 1, Nu, Du,
					offu, Alu, HSu, yu)
				_easi_g_parts(th, al, y, s, Pdn, z, lnx, Q, zm, cix, ao,
					lnx - rowsum(s :* Pdn), Pdn[, 1::neq] :- Pdn[, J], D0, ki,
					ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, 1, Nd, Dd,
					offd, Ald, HSd, yd)
				DNg[jj] = &((Nu - Nd) :/ (2 * eps))
				DDg[jj] = &((Du - Dd) :/ (2 * eps))
			}
		}
		KDg = _easi_impapply(DGh, n, Kt, GI, DT)
		if (C) {
			BIi = _easi_imppIi(al, p, lnx, z, Q, zm, cix, ao, s, wt)
			st_numscalar("_easi_imp_pchk", mreldif(_easi_imppscore(al, p, lnx, z, Q,
				zm, cix, ao, s, wt) * BIi, IFA))
			KDsc = _easi_impapply(DSc, n, cols(al), GI, DT)
			IFA  = IFA + KDsc * BIi
		}
	}

	Dt = J(Kt, Kt, .)
	for (j = 1; j <= Kt; j++) {
		h = 1e-6 * max((1, abs(th[j])))
		thp = th; thp[j] = thp[j] + h
		thm = th; thm[j] = thm[j] - h
		gp = _easi_g_gsum(thp, y, ystone, np, z, D0, PH, PHI, M0, Z0, Y, wt,
			cix, ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, lgCondSE)
		gm = _easi_g_gsum(thm, y, ystone, np, z, D0, PH, PHI, M0, Z0, Y, wt,
			cix, ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, lgCondSE)
		Dt[., j] = -(gp - gm) :/ (2 * h)
	}
	Gam = L * Dt
	Pm  = _easi_kinv(_easi_border(Gam, Rt), Kt)
	IFt = Gh + KDg
	na  = cols(al)
	if (C) {
		// the probits: Phi, phi, the Mills ratios and y move with alpha
		Da = J(Kt, na, .)
		for (j = 1; j <= na; j++) {
			h = 1e-6 * max((1, abs(al[j])))
			alp = al; alp[j] = alp[j] + h
			alm = al; alm[j] = alm[j] - h
			_easi_g_probs(alp, p, lnx, z, Q, zm, cix, ao, PH2, PHI2, M02)
			gp = _easi_g_gsum(th, y, ystone, np, z, D0, PH2, PHI2, M02, Z0, Y,
				wt, cix, ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, 0)
			_easi_g_probs(alm, p, lnx, z, Q, zm, cix, ao, PH2, PHI2, M02)
			gm = _easi_g_gsum(th, y, ystone, np, z, D0, PH2, PHI2, M02, Z0, Y,
				wt, cix, ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, 0)
			Da[., j] = (gp - gm) :/ (2 * h)
		}
		IFt = IFt + IFA * Da'
	}
	IFt = IFt * L' * Pm'
	if (rows(st_numscalar("_easi_ifkeep"))) st_matrix("_easi_ifkeep_b", IFt)
	if (vtype == 0) V = Pm * H * Pm'
	else            V = _easi_agg(IFt, vtype, clid, strid, fpcv)
	V = (V + V') / 2
	// the joint variance of the coefficients and of the probits (estat engel)
	if (C & vtype) st_matrix("_easi_Vsel", _easi_agg((IFt, IFA), vtype, clid, strid, fpcv))
	if (rows(st_numscalar("_easi_geng"))) st_matrix("_easi_g_IFt", IFt)

	// ---- elasticities ------------------------------------------------------
	b  = th[bidx]
	Vb = V[bidx, bidx]
	if (!C) {
		_easi_elast(b, Vb, s, p, z, y, wt, R, py, zy, pz, ipz, lgEZ, lgSE, lgB,
			lgESE, doESE, doCQ, IFt[., bidx], vtype, clid, strid, fpcv, "")
		if (doMkt) {
			_easi_elast(b, Vb, s, p, z, y, wt :* exp(lnx :- max(lnx)), R, py,
				zy, pz, ipz, lgEZ, lgSE, lgB, lgESE, doESE, doCQ, IFt[., bidx],
				vtype, clid, strid, fpcv, "_mkt")
		}
	}
	else {
		// the elasticities on the expected shares
		_easi_g_elast(th, al, y, s, p, z, lnx, Q, zm, cix, ao, ystone, np, D0,
			ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, doESE, doCQ, wt,
			IFt, IFA, vtype, clid, strid, fpcv, "", DNg, DDg, GI, DT)
		if (doMkt) {
			_easi_g_elast(th, al, y, s, p, z, lnx, Q, zm, cix, ao, ystone, np,
				D0, ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, doESE, doCQ,
				wt :* exp(lnx :- max(lnx)), IFt, IFA, vtype, clid, strid, fpcv,
				"_mkt", DNg, DDg, GI, DT)
		}
	}

	st_matrix(bnm, th')
	st_matrix(Vnm, V)
	st_matrix(Snm, Sigma)
	st_numscalar(itnm, it)
	st_numscalar(crnm, crit)
	st_numscalar(cvnm, conv)
}

// ---------------------------------------------------------------------------
// The elasticities under selection (step 3c): on the EXPECTED shares
//     E[w_i] = Phi_i f_i + delta_i phi_i    (a corrected good),  f_i otherwise,
// the last good by adding-up.  Household by household (xb = s'alpha_i):
//     dE_i/dln x  = Phi_i df_i/dln x  + phi_i (f_i - delta_i xb)(- sum_j a_pj)
//     dE_i/dln p_q = Phi_i df_i/dln p_q + phi_i (f_i - delta_i xb) a_pq
//     dE_i/dz_t    = Phi_i df_i/dz_t    + phi_i (f_i - delta_i xb) a_zt
// with the latent derivatives of EASI, y responding through the Stone index
// of the latent shares completed for the non-buyers (dy/dln p_q = -w~_q/Dn2).
// The means over the households of the types (households, individuals,
// market) as in phase D, the denominators the EXPECTED shares (decision DC2:
// Engel and Cournot aggregation exact); the compensated elasticities by the
// Slutsky equation household by household (the expected shares have no
// symmetry to rely on).  Variance: the sampling term, and the Jacobians in
// theta and alpha (central differences, y recomputed) times their influence
// functions.
// ---------------------------------------------------------------------------

// per-household numerators N and denominators D of the stacked elasticities
// (the order of _easi_estack), and the offsets (1 for the expenditure
// elasticities, -1 on the diagonal of the uncompensated ones); ER and EP (the
// latent semi-elasticities, means of dw/dy and of dw/dln p at fixed y) and
// the Slutsky matrix of the expected shares ride along
void _easi_g_parts(real colvector th, real rowvector al, real colvector y0,
	real matrix s, real matrix P, real matrix z, real colvector lnx,
	real matrix Q, real matrix zm, real rowvector cix, real rowvector ao,
	real colvector ystone, real matrix np, real matrix D0, real rowvector ki,
	real rowvector ko, real rowvector bidx, real rowvector didx, real scalar R,
	real scalar py, real scalar zy, real scalar pz, real rowvector ipz,
	real scalar lgQuant, real scalar doCQ, real matrix N, real matrix D,
	real rowvector off, real matrix Al, real matrix HS, real colvector y)
{
	real matrix C, PH, PHI, M0, XB, B, F, Wt, Yr, CD, Zp, Ai, ZA, E, DX, DZ
	real matrix DPl, DP, S
	real colvector cb, b, Dn2, kk, sa
	real rowvector dl, al_c
	real scalar n, J, neq, T, nipz, na, k, i, q, t, c, sx, r, e, rEP, rEZ, rEQ, nE
	real scalar obr, ogt, oht, obk, oA, onp, oynp, onpz
	real matrix bjr, gjt, hjt, bjk

	n = rows(s); J = cols(s); neq = J - 1; T = cols(z); nipz = cols(ipz)
	na = nipz + 1
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	obr = 0; ogt = obr + R * J; oht = ogt + T * J; obk = oht + T * J
	oA  = obk + J * J

	b  = th[bidx]
	dl = J(1, neq, 0)
	if (cols(didx)) dl[selectindex(cix :> 0)] = th[didx]'
	_easi_g_probs(al, P, lnx, z, Q, zm, cix, ao, PH, PHI, M0, XB)
	y  = _easi_g_y(y0, ystone, np, z, D0, M0, cix, b, dl, R, py, zy, pz, ipz,
		lgQuant)

	C   = _easi_ecomp(neq, R, T, zy, py, pz, ipz, k)
	cb  = C * b
	bjr = colshape(cb[| obr + 1 \ obr + R * J |], R)'
	gjt = colshape(cb[| ogt + 1 \ ogt + T * J |], T)'
	hjt = colshape(cb[| oht + 1 \ oht + T * J |], T)'
	bjk = colshape(cb[| obk + 1 \ obk + J * J |], J)'

	// the latent shares and the completed shares
	B = rowshape(b', neq)'
	F = _easi_design(y, z, np, R, py, zy, pz, ipz) * B
	F = F, (1 :- rowsum(F))
	Wt = s
	for (i = 1; i <= neq; i++) {
		if (!cix[i]) continue
		Wt[., i] = s[., i] :* (1 :- D0[., i]) + D0[., i] :* (F[., i] - dl[i] :* M0[., cix[i]])
	}
	Wt[., J] = 1 :- rowsum(Wt[| 1, 1 \ n, neq |])

	// the latent derivatives
	Yr = J(n, R, 1)
	for (r = 2; r <= R; r++) Yr[, r] = r :* y:^(r - 1)
	CD = Yr * bjr
	if (zy) CD = CD + z * hjt
	Al = CD
	if (py) Al = Al + P * bjk
	Dn2 = 1 :+ rowsum(CD :* P) :+ 0.5 :* rowsum((P * bjk) :* P)
	Zp = J(n, 1, 1)
	if (nipz) Zp = Zp, z[, ipz]
	DX  = Al :/ Dn2
	DPl = J(n, J * J, .)
	HS  = J(n, J * J, .)
	DZ  = J(n, T * J, .)
	for (i = 1; i <= J; i++) {
		Ai = colshape(cb[| oA + (i - 1) * na * J + 1 \ oA + i * na * J |], J)
		ZA = Zp * Ai
		for (q = 1; q <= J; q++) {
			HS[., (i - 1) * J + q]  = ZA[., q] + bjk[q, i] :* y
			DPl[., (i - 1) * J + q] = HS[., (i - 1) * J + q] - Wt[., q] :* DX[., i]
		}
		for (t = 1; t <= T; t++) {
			DZ[., (i - 1) * T + t] = gjt[t, i] :+ hjt[t, i] :* y
			if (pz & anyof(ipz, t)) {
				sx = selectindex(ipz :== t)[1]
				DZ[., (i - 1) * T + t] = DZ[., (i - 1) * T + t] + P * Ai[sx + 1, .]'
			}
		}
	}

	// the expected shares and their derivatives
	E  = F
	DP = DPl
	for (i = 1; i <= neq; i++) {
		if (!cix[i]) continue
		c    = cix[i]
		al_c = al[| ao[c] + 1 \ ao[c] + 1 + J + T |]
		kk   = PHI[., c] :* (F[., i] - dl[i] :* XB[., c])
		E[., i]  = PH[., c] :* F[., i] + dl[i] :* PHI[., c]
		DX[., i] = PH[., c] :* DX[., i] - kk :* sum(al_c[2..(J + 1)])
		for (q = 1; q <= J; q++) {
			DP[., (i - 1) * J + q] = PH[., c] :* DPl[., (i - 1) * J + q] + kk :* al_c[1 + q]
		}
		for (t = 1; t <= T; t++) {
			DZ[., (i - 1) * T + t] = PH[., c] :* DZ[., (i - 1) * T + t] + kk :* al_c[1 + J + t]
		}
	}
	// the last good closes the system
	E[., J]  = 1 :- rowsum(E[| 1, 1 \ n, neq |])
	DX[., J] = -rowsum(DX[| 1, 1 \ n, neq |])
	for (q = 1; q <= J; q++) {
		sa = J(n, 1, 0)
		for (i = 1; i <= neq; i++) sa = sa + DP[., (i - 1) * J + q]
		DP[., (J - 1) * J + q] = -sa
	}
	for (t = 1; t <= T; t++) {
		sa = J(n, 1, 0)
		for (i = 1; i <= neq; i++) sa = sa + DZ[., (i - 1) * T + t]
		DZ[., (J - 1) * T + t] = -sa
	}

	// stacked as _easi_estack: EI, vec(EPRICE) ([q, i] at (i-1)J + q),
	// vec(EZ) ([t, i] at (i-1)T + t), vec(EPQ) ([q, i]: the compensated
	// elasticity of good q with respect to the price of good i, plus one on
	// the diagonal)
	rEP = J; rEZ = rEP + J * J; rEQ = rEZ + T * J
	nE  = rEQ + (doCQ ? J * J : 0)
	N = J(n, nE, .); D = J(n, nE, 1); off = J(1, nE, 0)
	for (i = 1; i <= J; i++) {
		N[., i] = DX[., i]; D[., i] = E[., i]; off[i] = 1
		for (q = 1; q <= J; q++) {
			e = rEP + (i - 1) * J + q
			N[., e] = DP[., (i - 1) * J + q]; D[., e] = E[., i]
			off[e] = -(i == q)
			if (doCQ) {
				e = rEQ + (i - 1) * J + q
				N[., e] = DP[., (q - 1) * J + i] + E[., q] :* E[., i] + DX[., q] :* E[., i]
				D[., e] = E[., q]
			}
		}
		for (t = 1; t <= T; t++) N[., rEZ + (i - 1) * T + t] = DZ[., (i - 1) * T + t]
	}
	// the Slutsky matrix of the expected shares, [q, i], rides in HS
	S = J(n, J * J, .)
	for (i = 1; i <= J; i++) {
		for (q = 1; q <= J; q++) {
			S[., (i - 1) * J + q] = DP[., (i - 1) * J + q] + DX[., i] :* E[., q] +
				E[., i] :* E[., q] - (i == q) :* E[., i]
		}
	}
	HS = HS, S
}

// the weighted means, and the estimates
real rowvector _easi_g_est(real matrix N, real matrix D, real rowvector off,
	real colvector sw)
{
	return(quadcross(sw, N) :/ quadcross(sw, D) :+ off)
}

void _easi_g_elast(real colvector th, real rowvector al, real colvector y0,
	real matrix s, real matrix P, real matrix z, real colvector lnx,
	real matrix Q, real matrix zm, real rowvector cix, real rowvector ao,
	real colvector ystone, real matrix np, real matrix D0, real rowvector ki,
	real rowvector ko, real rowvector bidx, real rowvector didx, real scalar R,
	real scalar py, real scalar zy, real scalar pz, real rowvector ipz,
	real scalar lgQuant, real scalar doESE, real scalar doCQ,
	real colvector wt, real matrix IFt, real matrix IFA, real scalar vtype,
	real colvector clid, real colvector strid, real colvector fpc,
	string scalar sfx, | pointer(real matrix) rowvector DNg,
	pointer(real matrix) rowvector DDg, real matrix GI, real matrix DT)
{
	real matrix KDe
	real rowvector Nb, Db
	pointer(real matrix) rowvector Pk
	real scalar jj
	real matrix N, D, Np, Dp, Nm, Dm, Al, HS, Gt, Ga, Ts, IFe, Ve, EPRICE, EZ
	real matrix EPQ, EP, EPS, EPRse, EZse, EPQse
	real colvector sw, yy, thp, thm, sev
	real rowvector off, est, estp, estm, alp, alm, ws, EI, EIse, ER
	real scalar J, T, nE, j, h, e, rEP, rEZ, rEQ, K, na

	J = cols(s); T = cols(z)
	sw = wt :/ quadsum(wt)
	// the compensated elasticities are always computed, as by the usual
	// engine; compensated decides only whether their standard errors are
	_easi_g_parts(th, al, y0, s, P, z, lnx, Q, zm, cix, ao, ystone, np, D0,
		ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, 1, N, D, off, Al,
		HS, yy)
	est = _easi_g_est(N, D, off, sw)
	nE  = cols(est)
	if (rows(st_numscalar("_easi_ifkeep"))) st_matrix("_easi_ifkeep_e0" + sfx, est)
	rEP = J; rEZ = rEP + J * J; rEQ = rEZ + T * J

	EI     = est[1..J]
	EPRICE = rowshape(est[(rEP + 1)..(rEP + J * J)], J)'
	EZ     = rowshape(est[(rEZ + 1)..(rEZ + T * J)], J)'
	ws     = quadcross(sw, D[., 1..J])
	ER     = quadcross(sw, Al)
	EP     = rowshape(quadcross(sw, HS[., 1..(J * J)]), J)'
	EPS    = rowshape(quadcross(sw, HS[., (J * J + 1)..(2 * J * J)]), J)'
	EPQ = rowshape(est[(rEQ + 1)..nE], J)'

	EIse = J(1, J, .); EPRse = J(J, J, .); EZse = J(T, J, .); EPQse = J(J, J, .)
	if (doESE) {
		// the sampling term
		Ts = J(rows(N), nE, .)
		for (e = 1; e <= nE; e++) {
			Ts[., e] = _easi_rif(N[., e], D[., e], sw, h)
		}
		// pimpute(): the filled prices move N and D (expected shares) of
		// their households: d(Nbar/Dbar) = sw (dN Dbar - Nbar dD) / Dbar^2
		if (args() > 34) {
			if (rows(DT)) {
				Nb = quadcross(sw, N); Db = quadcross(sw, D)
				Pk = J(1, cols(DNg), NULL)
				for (jj = 1; jj <= cols(DNg); jj++) {
					if (DNg[jj] != NULL) Pk[jj] = &(sw :* ((*DNg[jj]) :* Db -
						Nb :* (*DDg[jj])) :/ (Db:^2))
				}
				KDe = _easi_impapply(Pk, rows(N), nE, GI, DT)
				Ts  = Ts + KDe
			}
		}
		// the Jacobian in theta (y recomputed from theta)
		K  = rows(th)
		Gt = J(nE, K, .)
		for (j = 1; j <= K; j++) {
			h = 1e-6 * max((1, abs(th[j])))
			thp = th; thp[j] = thp[j] + h
			thm = th; thm[j] = thm[j] - h
			_easi_g_parts(thp, al, yy, s, P, z, lnx, Q, zm, cix, ao, ystone, np,
				D0, ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, 1, Np, Dp,
				off, Al, HS, yy)
			estp = _easi_g_est(Np, Dp, off, sw)
			_easi_g_parts(thm, al, yy, s, P, z, lnx, Q, zm, cix, ao, ystone, np,
				D0, ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, 1, Nm, Dm,
				off, Al, HS, yy)
			estm = _easi_g_est(Nm, Dm, off, sw)
			Gt[., j] = (estp - estm)' :/ (2 * h)
		}
		// the Jacobian in alpha (Phi, phi, the Mills ratios and y move)
		na = cols(al)
		Ga = J(nE, na, .)
		for (j = 1; j <= na; j++) {
			h = 1e-6 * max((1, abs(al[j])))
			alp = al; alp[j] = alp[j] + h
			alm = al; alm[j] = alm[j] - h
			_easi_g_parts(th, alp, yy, s, P, z, lnx, Q, zm, cix, ao, ystone, np,
				D0, ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, 1, Np, Dp,
				off, Al, HS, yy)
			estp = _easi_g_est(Np, Dp, off, sw)
			_easi_g_parts(th, alm, yy, s, P, z, lnx, Q, zm, cix, ao, ystone, np,
				D0, ki, ko, bidx, didx, R, py, zy, pz, ipz, lgQuant, 1, Nm, Dm,
				off, Al, HS, yy)
			estm = _easi_g_est(Nm, Dm, off, sw)
			Ga[., j] = (estp - estm)' :/ (2 * h)
		}
		IFe = Ts + IFt * Gt' + IFA * Ga'
		if (rows(st_numscalar("_easi_ifkeep"))) st_matrix("_easi_ifkeep_e" + sfx, IFe)
		Ve  = _easi_agg(IFe, vtype, clid, strid, fpc)
		sev = sqrt(diagonal(Ve))
		EIse  = sev[1::J]'
		EPRse = rowshape(sev[(rEP + 1)::(rEP + J * J)]', J)'
		EZse  = rowshape(sev[(rEZ + 1)::(rEZ + T * J)]', J)'
		if (doCQ) EPQse = rowshape(sev[(rEQ + 1)::nE]', J)'
	}

	st_matrix("_easi_ei"    + sfx, EI)
	st_matrix("_easi_eise"  + sfx, EIse)
	st_matrix("_easi_ep"    + sfx, EPRICE)
	st_matrix("_easi_epse"  + sfx, EPRse)
	st_matrix("_easi_ez"    + sfx, EZ)
	st_matrix("_easi_ezse"  + sfx, EZse)
	st_matrix("_easi_er"    + sfx, ER)
	st_matrix("_easi_semip" + sfx, EP)
	st_matrix("_easi_slut"  + sfx, EPS)
	st_matrix("_easi_cq"    + sfx, EPQ)
	st_matrix("_easi_epnc"   + sfx, EPRICE')
	st_matrix("_easi_epncse" + sfx, EPRse')
	st_matrix("_easi_epc"    + sfx, EPQ - I(J))
	st_matrix("_easi_epcse"  + sfx, EPQse)
	st_numscalar("_easi_chkE" + sfx, abs(ws * EI' - 1))
	st_numscalar("_easi_chkC" + sfx, max(abs(EPRICE * ws' + ws')))
}

// placeholders (missing) for the elasticity matrices, until step 3c
void _easi_g_elastna(real scalar J, real scalar T, string scalar sfx)
{
	string rowvector nm
	real scalar i
	nm = ("ei", "eise", "er", "ep", "epse", "semip", "slut", "cq",
		"epnc", "epncse", "epc", "epcse")
	for (i = 1; i <= cols(nm); i++) {
		st_matrix("_easi_" + nm[i] + sfx,
			(anyof(("ei", "eise", "er"), nm[i]) ? J(1, J, .) : J(J, J, .)))
	}
	st_matrix("_easi_ez" + sfx, J(T, J, .))
	st_matrix("_easi_ezse" + sfx, J(T, J, .))
	st_numscalar("_easi_chkE" + sfx, 0)
	st_numscalar("_easi_chkC" + sfx, 0)
}

// vce(bootstrap): the estimates of one replication, stacked -- e(b), the
// elasticities of the households (and of the market), the probit
// coefficients (with selection)
real rowvector _easi_bvec(real scalar mkt, real scalar sel)
{
	real rowvector v, a
	real matrix AL
	string rowvector sf
	real scalar i
	v  = st_matrix("e(b)")
	sf = (mkt ? ("", "_mkt") : "")
	for (i = 1; i <= cols(sf); i++) {
		v = v, st_matrix("e(elast_exp" + sf[i] + ")"),
			vec(st_matrix("e(elast_price_nc" + sf[i] + ")"))',
			vec(st_matrix("e(elast_price_c" + sf[i] + ")"))',
			vec(st_matrix("e(elast_demo" + sf[i] + ")"))'
	}
	if (sel) {
		AL = st_matrix("e(sel_alpha)")
		a  = vec(AL)'
		v  = v, select(a, a :< .)
	}
	return(v)
}

// the variances of the replications, under the layout of _easi_bvec
void _easi_bpost(real matrix BS, real scalar mkt, real scalar sel)
{
	real matrix VB, E, Ep, Ec, Ez
	real scalar K, J, T, pos, i, na
	string rowvector sf
	K = cols(st_matrix("e(b)"))
	J = cols(st_matrix("e(elast_exp)"))
	T = rows(st_matrix("e(elast_demo)"))
	VB = quadvariance(BS)
	st_matrix("_easi_bV", VB[| 1, 1 \ K, K |])
	pos = K
	sf = (mkt ? ("", "_mkt") : "")
	for (i = 1; i <= cols(sf); i++) {
		E  = sqrt(diagonal(VB[| pos + 1, pos + 1 \ pos + J, pos + J |]))'
		pos = pos + J
		Ep = rowshape(sqrt(diagonal(VB[| pos + 1, pos + 1 \ pos + J * J, pos + J * J |]))', J)'
		pos = pos + J * J
		Ec = rowshape(sqrt(diagonal(VB[| pos + 1, pos + 1 \ pos + J * J, pos + J * J |]))', J)'
		pos = pos + J * J
		Ez = rowshape(sqrt(diagonal(VB[| pos + 1, pos + 1 \ pos + T * J, pos + T * J |]))', J)'
		pos = pos + T * J
		st_matrix("_easi_bse_exp" + sf[i], E)
		st_matrix("_easi_bse_price_nc" + sf[i], Ep)
		st_matrix("_easi_bse_price_c" + sf[i], Ec)
		st_matrix("_easi_bse_demo" + sf[i], Ez)
	}
	if (sel) {
		na = cols(BS) - pos
		st_matrix("_easi_bVs", (VB[| 1, 1 \ K, K |], VB[| 1, pos + 1 \ K, pos + na |] \
			VB[| pos + 1, 1 \ pos + na, K |], VB[| pos + 1, pos + 1 \ pos + na, pos + na |]))
	}
}

// Point estimates, then standard errors.
//
// -legacy(elastse)- keeps the R package's heuristic: the median across
// households of the pointwise delta-method standard error of the SEMI-elasticity,
// divided by the mean budget share.  A bootstrap audit (tests/AUDIT.md, section
// I) shows that heuristic is wrong in both directions -- 40% too large on the
// expenditure elasticities, 14 times too large on the dropped good, and right on
// the price elasticities only because two errors cancel.
//
// By default we take the delta method on the quantity actually reported:
// se = sqrt(g'Vg) with g the gradient of the reported (weighted-mean) elasticity
// with respect to b.  That gradient is analytic (_easi_ejac): every elasticity
// is differentiated in completed coordinates and the adding-up completion is
// applied once, as a matrix.  It agrees with forward differences on
// _easi_epoint to the finite-difference noise floor on every specification
// tried, and it removes what was 86% of the running time.  The expenditure
// elasticities are then re-done as a ratio of two estimators (_easi_eise).
void _easi_elast(real colvector b, real matrix V, real matrix s, real matrix P,
	real matrix z, real colvector y, real colvector wt, real scalar R, real scalar py,
	real scalar zy, real scalar pz, real rowvector ipz,
	real scalar lgEZ, real scalar lgSE, real scalar lgB, real scalar lgESE,
	real scalar doESE, real scalar doCQ, real matrix IFb, real scalar vtype,
	real colvector clid, real colvector strid, real colvector fpc,
	string scalar sfx, | pointer(real matrix) rowvector DN, real matrix GI,
	real matrix DT)
{
	pointer(real matrix) rowvector A, Pk
	real matrix NN0, DD0, KDe
	real colvector swc, est0
	real rowvector Dbar
	real scalar jj
	real matrix B, Zc, bjk, bjr, gjt, hjt, DD, MAT
	real matrix EP, EPS, EPQ, EPRICE, EZ, EPse, EPRse, EZse, G, Gfd, Ve, EPQse
	real matrix Ts, Tb, Gy
	real colvector tot2, Ci, Di, Gi, Ui, Ei, Fi, Hv, cum, den1, Dn2
	real colvector e0, e1, bp, sev, est
	real rowvector ws, ER, EI, ERse, EIse, idx
	real scalar n, Jg, neq, T, nipz, k, onp, oynp, onpz, K, nE, h
	real scalar i, j, q, t, r, d, sIdx, den, mG, mCDG, Bsc, Hm, acc

	n = rows(s); Jg = cols(s); neq = Jg - 1; T = cols(z); nipz = cols(ipz)
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	K    = neq * k

	_easi_epoint(b, s, P, z, y, wt, R, py, zy, pz, ipz, lgEZ, lgB,
		EI, ER, EPRICE, EP, EPS, EPQ, EZ, ws, Dn2)

	if (lgESE) {
		// ---- standard errors (median over households of the delta method) ----
		ERse = J(1, neq + 1, 0); EIse = J(1, neq + 1, 0)
		for (i = 1; i <= neq; i++) {
			idx = J(1, 0, 0); MAT = J(n, 0, 0)
			for (r = 1; r <= R; r++) {
				idx = idx, ((i - 1) * k + 1 + r)
				MAT = MAT, (r :* y:^(r - 1))
			}
			if (zy) {
				for (t = 1; t <= T; t++) idx = idx, ((i - 1) * k + 1 + R + T + t)
				MAT = MAT, z
			}
			if (py) {
				for (j = 1; j <= neq; j++) idx = idx, ((i - 1) * k + oynp + j)
				MAT = MAT, P[, 1::neq]
			}
			DD = V[idx, idx]
			ERse[i] = _easi_dse(MAT, DD, wt)
			// the point estimate carries a 1/Dn2 factor, so its Jacobian does too
			if (lgB) EIse[i] = _easi_dse(MAT, DD, wt) / ws[i]
			else     EIse[i] = _easi_dse(MAT :/ Dn2, DD, wt) / ws[i]
		}
		// R bug 3: the last item's standard error is given a minus sign
		ERse[neq + 1] = (lgSE ? -1 : 1) * sqrt(sum(ERse[1::neq]:^2))
		EIse[neq + 1] = ERse[neq + 1] / ws[Jg]

		EPse  = J(neq + 1, neq + 1, 0)
		EPRse = J(neq + 1, neq + 1, 0)
		for (i = 1; i <= neq; i++) {
			for (j = 1; j <= neq; j++) {
				idx = ((i - 1) * k + onp + j)
				MAT = J(n, 1, 1)
				if (pz) {
					for (sIdx = 1; sIdx <= nipz; sIdx++) {
						idx = idx, ((i-1)*k + onpz + (sIdx-1)*neq + j)
					}
					MAT = MAT, z[, ipz]
				}
				if (py) {
					idx = idx, ((i - 1) * k + oynp + j)
					MAT = MAT, y
				}
				DD = V[idx, idx]
				EPse[i, j]  = _easi_dse(MAT, DD, wt)
				EPRse[i, j] = EPse[i, j] / ws[j]
			}
		}
		EPse[neq + 1, 1::neq] = sqrt(rowsum(EPse[1::neq, 1::neq]:^2))'
		EPse[., neq + 1]      = EPse[neq + 1, .]'
		// R bug 4: the squares are missing in these two closing cells
		EPse[neq + 1, neq + 1] = lgSE					///
			? sqrt(sum(EPse[neq + 1, 1::neq]))			///
			: sqrt(sum(EPse[neq + 1, 1::neq]:^2))
		for (i = 1; i <= neq + 1; i++) EPRse[i, neq + 1] = EPse[i, neq + 1] / ws[Jg]
		EPRse[neq + 1, 1::neq] = EPRse[1::neq, neq + 1]'
		EPRse[neq + 1, neq + 1] = lgSE				///
			? sqrt(sum(EPRse[1::neq, neq + 1]))			///
			: sqrt(sum(EPRse[1::neq, neq + 1]:^2))

		EZse = J(T, neq + 1, 0)
		for (i = 1; i <= neq; i++) {
			for (j = 1; j <= T; j++) {
				idx = ((i - 1) * k + 1 + R + j)
				MAT = J(n, 1, 1)
				if (zy) {
					idx = idx, ((i - 1) * k + 1 + R + T + j)
					MAT = MAT, y
				}
				if (pz & anyof(ipz, j)) {
					sIdx = selectindex(ipz :== j)[1]
					for (t = 1; t <= neq; t++) {
						idx = idx, ((i-1)*k + onpz + (sIdx-1)*neq + t)
					}
					MAT = MAT, P[, 1::neq]
				}
				DD = V[idx, idx]
				EZse[j, i] = _easi_dse(MAT, DD, wt)
			}
		}
		for (j = 1; j <= T; j++) EZse[j, neq + 1] = sqrt(sum(EZse[j, 1::neq]:^2))

	}
	else if (!doESE) {
		EIse = J(1, Jg, .); EPRse = J(Jg, Jg, .)
		EZse = J(T, Jg, .); ERse = J(1, Jg, .)
	}
	else {
		// ---- delta method on the reported elasticities -------------------
		e0 = _easi_estack(EI, EPRICE, EZ, EPQ, doCQ)
		nE = rows(e0)
		if (!lgB & !lgEZ) {
			// analytic, in completed coordinates times the completion map
			G = _easi_ejac(b, s, P, z, y, wt, Dn2, ws, R, py, zy, pz, ipz,
				doCQ)
		}
		// The legacy formulas keep the finite-difference Jacobian.  A test
		// can also request it beside the analytic one by defining the Stata
		// scalar _easi_jaccheck: both are then left in _easi_Gan / _easi_Gfd
		// (tests/test_step16.do).  Nothing else reads that scalar.
		if (lgB | lgEZ | rows(st_numscalar("_easi_jaccheck"))) {
			Gfd = J(nE, K, 0)
			for (j = 1; j <= K; j++) {
				h  = 1e-6 * max((1, abs(b[j])))
				bp = b
				bp[j] = bp[j] + h
				_easi_epoint(bp, s, P, z, y, wt, R, py, zy, pz, ipz,
					lgEZ, lgB, EI, ER, EPRICE, EP, EPS, EPQ, EZ, ws, Dn2)
				e1 = _easi_estack(EI, EPRICE, EZ, EPQ, doCQ)
				Gfd[, j] = (e1 - e0) :/ h
			}
			// restore the point estimates
			_easi_epoint(b, s, P, z, y, wt, R, py, zy, pz, ipz, lgEZ, lgB,
				EI, ER, EPRICE, EP, EPS, EPQ, EZ, ws, Dn2)
			if (lgB | lgEZ) G = Gfd
			else {
				st_matrix("_easi_Gan" + sfx, G)
				st_matrix("_easi_Gfd" + sfx, Gfd)
			}
		}

		if (lgB | lgEZ) {
			// the legacy formulas: delta method on beta alone, with
			// the expenditure elasticities as a ratio (easi 1.0.0)
			Ve  = G * V * G'
			sev = sqrt(diagonal(Ve))
			EIse  = sev[|1 \ Jg|]'
			_easi_eise(s, P[, 1::neq] :- P[, Jg], P[, Jg], z, y, wt,
				IFb, Dn2, ws, b, R, py, zy, vtype, clid, strid, fpc,
				lgB, EIse)
		}
		else {
			// Every elasticity is a ratio of weighted means over the
			// households: its influence function is the sampling term of
			// those means plus the estimation term IF(beta) G'.  easi
			// 1.0.0 did this for the expenditure elasticities only
			// (_easi_eise, which the stacked form reproduces) and took
			// G V G' for the others.  Under vce(conventional) the
			// coefficient part stays G V G', the conventional variance;
			// otherwise the whole influence function goes through the
			// design aggregator, whose coefficient part is exactly
			// G V G' (V is the aggregate of IF(beta)).
			Ts = _easi_esamp(b, s, P, z, y, wt, Dn2, R, py, zy, pz, ipz,
				doCQ, est)
			// pimpute(): the imputation of the prices moves the numerators of
			// the households whose price was filled
			if (args() > 24) {
				if (rows(DT)) {
					NN0 = .; DD0 = .; est0 = .
					(void) _easi_esamp(b, s, P, z, y, wt, Dn2, R, py, zy,
						pz, ipz, doCQ, est0, NN0, DD0)
					swc  = wt :/ quadsum(wt)
					Dbar = quadcolsum(swc :* DD0)
					Pk = J(1, cols(DN), NULL)
					for (jj = 1; jj <= cols(DN); jj++) {
						if (DN[jj] != NULL) Pk[jj] = &((swc :* (*DN[jj])) :/ Dbar)
					}
					KDe = _easi_impapply(Pk, rows(s), cols(Ts), GI, DT)
					Ts  = Ts + KDe
				}
			}
			// self-test: the numerators rebuilt for the influence
			// function give back the point estimates
			st_numscalar("_easi_esampchk" + sfx, max(abs(est - e0)))
			// the rest of the total derivative: y is a function of the
			// coefficients (the analytic Jacobian above holds it fixed)
			Gy = _easi_ejacy(b, s, P, z, y, wt, Dn2, R, py, zy, pz, ipz, doCQ)
			if (rows(st_numscalar("_easi_jaccheck"))) st_matrix("_easi_Gy" + sfx, Gy)
			G  = G + Gy
			Tb = IFb * G'
			if (rows(st_numscalar("_easi_ifkeep"))) {
				st_matrix("_easi_ifkeep_e" + sfx, Ts + Tb)
				st_matrix("_easi_ifkeep_e0" + sfx, e0')
			}
			if (vtype == 0) {
				Ve = quadcross(Ts, Ts) + quadcross(Ts, Tb)
				Ve = Ve + Ve' - quadcross(Ts, Ts) + G * V * G'
			}
			else Ve = _easi_agg(Ts + Tb, vtype, clid, strid, fpc)
			Ve  = (Ve + Ve') / 2
			sev = sqrt(diagonal(Ve))
			EIse = sev[|1 \ Jg|]'
		}
		EPRse = rowshape(sev[|Jg + 1 \ Jg + Jg * Jg|]', Jg)'
		EZse  = rowshape(sev[|Jg + Jg * Jg + 1 \ Jg + Jg * Jg + T * Jg|]', Jg)'
		ERse  = J(1, Jg, .)
		if (doCQ) {
			EPQse = rowshape(sev[|Jg + Jg * Jg + T * Jg + 1 \ nE|]', Jg)'
		}
	}
	if (!doCQ | doESE == 0 | lgESE) EPQse = J(Jg, Jg, .)

	st_matrix("_easi_ei"    + sfx, EI)
	st_matrix("_easi_eise"  + sfx, EIse)
	st_matrix("_easi_ep"    + sfx, EPRICE)
	st_matrix("_easi_epse"  + sfx, EPRse)
	st_matrix("_easi_ez"    + sfx, EZ)
	st_matrix("_easi_ezse"  + sfx, EZse)
	st_matrix("_easi_er"    + sfx, ER)
	st_matrix("_easi_semip" + sfx, EP)
	st_matrix("_easi_slut"  + sfx, EPS)
	st_matrix("_easi_cq"    + sfx, EPQ)

	// Orientation.  EPRICE is indexed [price, good] -- the second index is the
	// good, as the Cournot identity confirms to machine precision -- while EPQ
	// comes out [good, price].  The two stored matrices therefore read in
	// opposite directions, which is inherited from the R package and is a trap:
	// both axes carry the same good labels, so nothing tells the reader to
	// transpose.  The elast_price_nc / elast_price_c pair below is the one the
	// command displays and documents, always [good, price]: row j, column k is
	// the elasticity of good j with respect to the price of good k.  The legacy
	// names keep the legacy content.
	//
	// EPQ is eta^H + I, not eta^H: EPQ = diag(w)^-1 (EP + w'w) and the Slutsky
	// matrix in share form is EPS = EP + w'w - diag(w), so diag(w)^-1 EPS =
	// EPQ - I.  Subtracting the identity is not cosmetic -- it is the whole
	// diagonal, that is every own-price compensated elasticity.
	st_matrix("_easi_epnc"   + sfx, EPRICE')
	st_matrix("_easi_epncse" + sfx, EPRse')
	st_matrix("_easi_epc"    + sfx, EPQ - I(Jg))
	st_matrix("_easi_epcse"  + sfx, EPQse)

	// Aggregation identities, checked on the user's own data.  Engel:
	// sum_j w_j eta^x_j = 1.  Cournot: sum_j w_j eta^k_j = -w_k for every k.
	// Both hold exactly for the mean-based formulas, because sum_j aleph_j = 0
	// and sum_j a^jk = 0; measured residual is 1e-16 to 1e-10.  They are a
	// self-test of the elasticity code run on whatever data is at hand.
	st_numscalar("_easi_chkE" + sfx, abs(ws * EI' - 1))
	st_numscalar("_easi_chkC" + sfx, max(abs(EPRICE * ws' + ws')))
}



void _easi_run(string scalar svars, string scalar lpvars, string scalar lxvar,
	string scalar zvars, string scalar touse, string scalar wvar,
	real scalar neff, real scalar vtype, string scalar clvar,
	string scalar strvar, string scalar psuvar, string scalar fpcvar,
	real scalar R,
	real scalar py, real scalar pz, real scalar zy, string scalar sipz,
	real scalar lgInstr, real scalar lgQuant, real scalar lgEZ, real scalar lgSE,
	real scalar lgB, real scalar lgESE, real scalar doESE, real scalar doCQ,
	real scalar lgCondSE, real scalar maxsig,
	real scalar tol, real scalar maxit, real scalar noisy,
	string scalar bnm, string scalar Vnm, string scalar Snm,
	string scalar itnm, string scalar crnm, string scalar cvnm,
	real scalar doMkt, string scalar lp0vars, string scalar impgvars,
	real scalar doImp)
{
	real matrix s, p, z, np, X, Z, Y, Rr, V, Sigma, B, Aco, Bco, Pu, IFb
	real matrix LP0, GR, GI, DT, ZZiZX, Si, S0, KD, Pp, Pm, IFb0
	real colvector b1, rr
	pointer(real matrix) rowvector Dsc, DN
	real scalar jj, eps
	real colvector lnx, ystone, ytil, y, yold, yinst, b, bold, pAp, pBp
	real colvector wt, clid, strid, fpcv
	real rowvector ipz, ms
	real scalar n, J, neq, T, k, i, j, t, it, crit, conv
	real scalar onp, oynp, onpz, nipz, interact, ph, nphase

	s   = st_data(., svars,  touse)
	p   = st_data(., lpvars, touse)
	lnx = st_data(., lxvar,  touse)
	z   = st_data(., zvars,  touse)
	wt  = st_data(., wvar,   touse)
	clid  = (clvar  == "" ? J(0, 1, 0) : st_data(., clvar,  touse))
	strid = (strvar == "" ? J(0, 1, 0) : st_data(., strvar, touse))
	fpcv  = (fpcvar == "" ? J(0, 1, 0) : st_data(., fpcvar, touse))
	if (psuvar != "") clid = st_data(., psuvar, touse)
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))

	n = rows(s); J = cols(s); neq = J - 1; T = cols(z); nipz = cols(ipz)

	np     = p[, 1::neq] :- p[, J]
	ystone = lnx - rowsum(s :* p)
	ms     = quadcross(wt, s) :/ quadsum(wt)	// weighted mean shares
	ytil   = lnx - p * ms'

	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)

	Rr       = _easi_restrict(neq, k, onp, py, oynp, pz, onpz, nipz)
	Y        = s[, 1::neq]
	interact = (py | pz | zy)

	// Pendakur's scheme (EASI made Easier, appendix 8.2) runs two passes: the
	// first iterates y with the raw instrument y_tilda and produces p'Ap, p'Bp;
	// the instrument is then rebuilt ONCE as (y_tilda + p'Ap/2)/(1 - p'Bp/2),
	// powers and interactions included, and held fixed through the second pass.
	// The R package rebuilds it at the end of every first-pass iteration, never
	// refreshes the y_inst^r powers, and then resets its data frame, discarding
	// all of it -- so it estimates on the raw y_tilda.  -legacy(instr)-
	// reproduces that; otherwise we follow Pendakur.
	nphase = lgInstr ? 1 : 2
	yinst  = ytil
	pAp = J(n, 1, 0); pBp = J(n, 1, 0)

	if (noisy) printf("\n{txt}iterating on the implicit utility index y...\n")

	for (ph = 1; ph <= nphase; ph++) {
	y = ystone
	b = J(neq * k, 1, 0)
	crit = 1; conv = 0; it = 0
	Z = _easi_instr(yinst, z, np, R, py, zy, pz, ipz)
	if (noisy & nphase > 1) printf("{txt}  -- passe %1.0f sur %1.0f\n", ph, nphase)

	while (it < maxit) {
		it++
		bold = b
		yold = y

		X = _easi_design(y, z, np, R, py, zy, pz, ipz)
		_easi_3sls(X, Z, Y, Rr, wt, neff, 0, clid, strid, fpcv, maxsig, b, V, Sigma, Pu)

		_easi_quadf(b, np, z, R, py, pz, zy, ipz, lgQuant, pAp, pBp)

		y = (ystone :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)

		if (interact) {
			if (it > 1) crit = sum((b - bold):^2)
		}
		else crit = max(abs(y - yold))

		if (noisy) printf("{txt}  iteration %3.0f    criterion = %12.5e\n", it, crit)

		if (crit <= tol) {
			conv = 1
			break
		}
	}

	// between passes, rebuild the instrument once (b1, the coefficients
	// that build it, for the influence of the imputation of the prices)
	if (ph == 1) b1 = b
	if (ph < nphase) yinst = (ytil :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)
	}

	// final refit at the converged y
	X = _easi_design(y, z, np, R, py, zy, pz, ipz)
	_easi_3sls(X, Z, Y, Rr, wt, neff, vtype, clid, strid, fpcv, maxsig, b, V, Sigma, Pu)

	// y is a generated regressor: add its contribution to the Jacobian unless
	// the conditional covariance was asked for.  The elasticity standard errors
	// are delta-method on this V, so they inherit the correction.
	if (!lgCondSE) {
		_easi_gencorr(X, Z, Y, Rr, wt, np, z, y, pBp, Sigma, b, R,
			py, zy, pz, ipz, vtype, clid, strid, fpcv, V, Pu)
	}

	_easi_ifb(X, Z, Y, wt, Sigma, b, Pu, IFb)

	// pimpute(): the influence of the imputation of the prices.  The
	// map of the fills (GI, DT), then, good by good, the derivatives of the
	// scores and of the numerators of the elasticities in the filled prices
	// (central differences; each household's depend on its own prices only,
	// at fixed coefficients, Sigma, first stage and instrument coefficients)
	DN = J(1, J, NULL); GI = J(0, 0, .); DT = J(0, 4, .)
	if (doImp & vtype) {
		LP0 = st_data(., lp0vars, touse)
		GR  = st_data(., impgvars, touse)
		_easi_impmap(p, LP0, GR, wt, GI, DT)
	}
	if (rows(DT)) {
		ZZiZX = cholsolve(quadcross(Z, wt, Z), quadcross(Z, wt, X))
		if (hasmissing(ZZiZX)) ZZiZX = qrsolve(quadcross(Z, wt, Z), quadcross(Z, wt, X))
		Si = invsym(Sigma)
		S0 = _easi_impS(b, b1, nphase, Si, ZZiZX, s, p, lnx, z, wt, ms, R, py,
			zy, pz, ipz, lgQuant, k)
		st_numscalar("_easi_imp_schk", mreldif(S0 * Pu', IFb))
		Dsc = J(1, J, NULL)
		eps = 1e-5
		for (jj = 1; jj <= J; jj++) {
			rr = selectindex(GI[., jj] :> 0)
			if (rows(rr) == 0) continue
			Pp = p; Pp[rr, jj] = Pp[rr, jj] :+ eps
			Pm = p; Pm[rr, jj] = Pm[rr, jj] :- eps
			Dsc[jj] = &((_easi_impS(b, b1, nphase, Si, ZZiZX, s, Pp, lnx, z, wt,
				ms, R, py, zy, pz, ipz, lgQuant, k) - _easi_impS(b, b1, nphase,
				Si, ZZiZX, s, Pm, lnx, z, wt, ms, R, py, zy, pz, ipz, lgQuant,
				k)) :/ (2 * eps))
			DN[jj] = &((_easi_impNN(b, s, Pp, z, lnx, wt, R, py, zy, pz, ipz,
				lgQuant, lgEZ, lgB, doCQ) - _easi_impNN(b, s, Pm, z, lnx, wt, R,
				py, zy, pz, ipz, lgQuant, lgEZ, lgB, doCQ)) :/ (2 * eps))
		}
		KD   = _easi_impapply(Dsc, n, cols(S0), GI, DT)
		IFb0 = IFb
		IFb  = IFb + KD * Pu'
		V = V + _easi_agg(IFb, vtype, clid, strid, fpcv) -
			_easi_agg(IFb0, vtype, clid, strid, fpcv)
		V = (V + V') / 2
	}
	if (rows(st_numscalar("_easi_ifkeep"))) st_matrix("_easi_ifkeep_b", IFb)

	// the elasticities of the households: every mean weighted by the
	// weight of the estimation
	_easi_elast(b, V, s, p, z, y, wt, R, py, zy, pz, ipz, lgEZ, lgSE, lgB,
		lgESE, doESE, doCQ, IFb, vtype, clid, strid, fpcv, "", DN, GI, DT)
	// the elasticities of the market, from the same estimation: every mean
	// weighted by the weight times the total expenditure, which makes each
	// one a ratio of totals -- the elasticity of the total demand of the
	// sample (exp(ln x) up to a constant factor, which cancels; the largest
	// ln x is subtracted against overflow)
	if (doMkt) {
		_easi_elast(b, V, s, p, z, y, wt :* exp(lnx :- max(lnx)), R, py,
			zy, pz, ipz, lgEZ, lgSE, lgB, lgESE, doESE, doCQ, IFb, vtype,
			clid, strid, fpcv, "_mkt", DN, GI, DT)
	}

	st_matrix(bnm, b')
	st_matrix(Vnm, V)
	st_matrix(Snm, Sigma)
	st_numscalar(itnm, it)
	st_numscalar(crnm, crit)
	st_numscalar(cvnm, conv)
}


// Engel curves proper: the fitted budget shares over a grid of total
// expenditure, with demographics and prices held at their (weighted) means and
// the residual at zero.  This isolates the expenditure effect, which is what an
// Engel curve is; binning the sample by expenditure instead -- what the R
// package does -- also picks up however demographics happen to covary with it.
// At each grid point (w, y) solves the model's fixed point.
void _easi_engelgrid(string scalar svars, string scalar lpvars,
	string scalar lxvar, string scalar zvars, string scalar touse,
	string scalar wtvar, real scalar G, real scalar trim, real scalar R,
	real scalar py, real scalar pz, real scalar zy, string scalar sipz,
	real scalar lgQuant, string scalar outnm)
{
	real matrix s, p, z, zg, pg, npg, X, Bm, What, Wfull, SE, V, Vj, Vs, Aco, Bco
	real colvector lnx, wt, lnxg, y, ynew, ystone, pAp, pBp, b, q
	real rowvector ipz, zbar, pbar
	real scalar n, Jg, neq, T, k, i, j, t, m, it, nipz, onp, oynp, onpz

	s   = st_data(., svars,  touse)
	p   = st_data(., lpvars, touse)
	lnx = st_data(., lxvar,  touse)
	z   = st_data(., zvars,  touse)
	wt  = st_data(., wtvar,  touse)
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))
	b   = st_matrix("e(b)")'
	V   = st_matrix("e(V)")

	n = rows(s); Jg = cols(s); neq = Jg - 1; T = cols(z); nipz = cols(ipz)
	zbar = quadcross(wt, z) :/ quadsum(wt)
	pbar = quadcross(wt, p) :/ quadsum(wt)

	// the grid: equally spaced percentiles of observed total expenditure
	q    = (trim :+ (100 - 2 * trim) :* (((1::G) :- 0.5) :/ G)) :/ 100
	lnxg = J(G, 1, 0)
	for (i = 1; i <= G; i++) lnxg[i] = _easi_wquant(lnx, wt, q[i])

	zg  = J(G, 1, 1) * zbar
	pg  = J(G, 1, 1) * pbar
	npg = pg[, 1::neq] :- pg[, Jg]

	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	Bm   = rowshape(b', neq)'

	Aco = Bm[| onp + 1, 1 \ onp + neq, neq |]'
	pAp = J(G, 1, 0)
	for (i = 1; i <= neq; i++) {
		for (j = 1; j <= neq; j++) pAp = pAp + Aco[i, j] :* npg[, i] :* npg[, j]
	}
	if (pz) {
		for (t = 1; t <= nipz; t++) {
			Aco = Bm[| onpz + (t-1)*neq + 1, 1 \ onpz + t*neq, neq |]'
			for (i = 1; i <= neq; i++) {
				for (j = 1; j <= neq; j++) {
					pAp = pAp + Aco[i, j] :* npg[, i] :* npg[, j] :* zg[, ipz[t]]
				}
			}
		}
	}
	pBp = J(G, 1, 0)
	if (py) {
		Bco = Bm[| oynp + 1, 1 \ oynp + neq, neq |]'
		for (i = 1; i <= neq; i++) {
			for (j = 1; j <= neq; j++) pBp = pBp + Bco[i, j] :* npg[, i] :* npg[, j]
		}
	}
	if (lgQuant) {
		pAp = _easi_quant(pAp)
		pBp = _easi_quant(pBp)
	}

	y = lnxg
	for (it = 1; it <= 500; it++) {
		X     = _easi_design(y, zg, npg, R, py, zy, pz, ipz)
		What  = X * Bm
		Wfull = What, (1 :- rowsum(What))
		ystone = lnxg - rowsum(Wfull :* pg)
		ynew   = (ystone :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)
		if (max(abs(ynew - y)) < 1e-13) {
			y = ynew
			break
		}
		y = ynew
	}
	X     = _easi_design(y, zg, npg, R, py, zy, pz, ipz)
	What  = X * Bm
	Wfull = What, (1 :- rowsum(What))

	SE = J(G, Jg, 0)
	Vs = J(k, k, 0)
	for (i = 1; i <= neq; i++) {
		Vj = V[| (i-1)*k + 1, (i-1)*k + 1 \ i*k, i*k |]
		SE[, i] = sqrt(rowsum((X * Vj) :* X))
		for (m = 1; m <= neq; m++) {
			Vs = Vs + V[| (i-1)*k + 1, (m-1)*k + 1 \ i*k, m*k |]
		}
	}
	SE[, Jg] = sqrt(rowsum((X * Vs) :* X))

	st_matrix(outnm, (100 :* q, lnxg, y, Wfull, SE))
}

// weighted quantile
real scalar _easi_wquant(real colvector x, real colvector wt, real scalar q)
{
	real matrix XW
	real colvector cw
	real scalar n, i, tgt

	n   = rows(x)
	XW  = sort((x, wt), 1)
	cw  = runningsum(XW[, 2])
	tgt = q * cw[n]
	for (i = 1; i <= n; i++) {
		if (cw[i] >= tgt) return(XW[i, 1])
	}
	return(XW[n, 1])
}

// Fitted budget shares and their delta-method standard errors, for the Engel
// curves.  The share of the dropped good is 1 - sum of the others, so its
// variance is x'(sum_jk V_jk)x -- the full sum, covariances included, not the
// root-sum-of-squares the R package uses.
void _easi_engel(string scalar svars, string scalar lpvars, string scalar lxvar,
	string scalar zvars, string scalar touse, string scalar wvars,
	string scalar sevars, real scalar R, real scalar py, real scalar pz,
	real scalar zy, string scalar sipz, real scalar lgQuant)
{
	real matrix s, p, z, np, X, Bm, What, Wfull, SE, V, Vj, Vs, Aco, Bco
	real colvector lnx, ystone, y, pAp, pBp, b
	real rowvector ipz
	real scalar n, Jg, neq, T, k, i, j, t, m, nipz, onp, oynp, onpz

	s   = st_data(., svars,  touse)
	p   = st_data(., lpvars, touse)
	lnx = st_data(., lxvar,  touse)
	z   = st_data(., zvars,  touse)
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))
	b   = st_matrix("e(b)")'
	V   = st_matrix("e(V)")

	n = rows(s); Jg = cols(s); neq = Jg - 1; T = cols(z); nipz = cols(ipz)
	np     = p[, 1::neq] :- p[, Jg]
	ystone = lnx - rowsum(s :* p)

	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	Bm   = rowshape(b', neq)'

	Aco = Bm[| onp + 1, 1 \ onp + neq, neq |]'
	pAp = J(n, 1, 0)
	for (i = 1; i <= neq; i++) {
		for (j = 1; j <= neq; j++) pAp = pAp + Aco[i, j] :* np[, i] :* np[, j]
	}
	if (pz) {
		for (t = 1; t <= nipz; t++) {
			Aco = Bm[| onpz + (t-1)*neq + 1, 1 \ onpz + t*neq, neq |]'
			for (i = 1; i <= neq; i++) {
				for (j = 1; j <= neq; j++) {
					pAp = pAp + Aco[i, j] :* np[, i] :* np[, j] :* z[, ipz[t]]
				}
			}
		}
	}
	pBp = J(n, 1, 0)
	if (py) {
		Bco = Bm[| oynp + 1, 1 \ oynp + neq, neq |]'
		for (i = 1; i <= neq; i++) {
			for (j = 1; j <= neq; j++) pBp = pBp + Bco[i, j] :* np[, i] :* np[, j]
		}
	}
	if (lgQuant) {
		pAp = _easi_quant(pAp)
		pBp = _easi_quant(pBp)
	}
	y = (ystone :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)

	X     = _easi_design(y, z, np, R, py, zy, pz, ipz)
	What  = X * Bm
	Wfull = What, (1 :- rowsum(What))

	SE = J(n, Jg, 0)
	Vs = J(k, k, 0)
	for (i = 1; i <= neq; i++) {
		Vj = V[| (i-1)*k + 1, (i-1)*k + 1 \ i*k, i*k |]
		SE[, i] = sqrt(rowsum((X * Vj) :* X))
		for (m = 1; m <= neq; m++) {
			Vs = Vs + V[| (i-1)*k + 1, (m-1)*k + 1 \ i*k, m*k |]
		}
	}
	SE[, Jg] = sqrt(rowsum((X * Vs) :* X))

	st_store(., tokens(wvars),  touse, Wfull)
	st_store(., tokens(sevars), touse, SE)
}

// The layout of an estimation with selection, read from e(): the corrected
// goods (cix: 0 or their rank), the positions of b and delta in e(b), the
// probit coefficients stacked good by good (from e(sel_alpha)) and their
// offsets, the mask of the probit-only variables.
void _easi_g_layout(string scalar sgstr, string scalar zmstr, real scalar neq,
	real scalar k, real scalar nq, real rowvector cix, real rowvector ki,
	real rowvector ko, real rowvector bidx, real rowvector didx,
	real rowvector ao, real rowvector al, real matrix zm)
{
	real matrix AL
	real colvector col
	real scalar i, C, c

	AL  = st_matrix("e(sel_alpha)")
	cix = strtoreal(tokens(sgstr))
	C = 0
	for (i = 1; i <= neq; i++) {
		if (cix[i]) {
			C++
			cix[i] = C
		}
	}
	zm = (nq ? rowshape(strtoreal(tokens(zmstr)), neq) : J(neq, 0, .))
	ki = J(1, neq, k) + (cix :> 0)
	ko = J(1, neq, 0)
	for (i = 2; i <= neq; i++) ko[i] = ko[i - 1] + ki[i - 1]
	bidx = J(1, 0, .); didx = J(1, 0, .); ao = J(1, C, 0); al = J(1, 0, .)
	for (i = 1; i <= neq; i++) {
		bidx = bidx, (ko[i] + 1)..(ko[i] + k)
		if (!cix[i]) continue
		didx = didx, ko[i] + k + 1
		c = cix[i]
		col = AL[., i]
		col = select(col, col :< .)
		if (c < C) ao[c + 1] = ao[c] + rows(col)
		al = al, col'
	}
}

// the expected shares at the rows of X (the latent design) and S-rows of
// the probits, and their standard errors from the joint variance of
// (theta, alpha); y held fixed, as in the usual Engel curves
void _easi_g_ewse(real matrix X, real matrix P, real colvector lnx,
	real matrix z, real matrix Q, real matrix zm, real rowvector cix,
	real rowvector ki, real rowvector ko, real rowvector ao, real rowvector al,
	real colvector th, real rowvector bidx, real rowvector didx,
	real matrix Vj, real matrix E, real matrix SE)
{
	real matrix PH, PHI, M0, XB, F, G, Gs, S
	real rowvector dl
	real scalar n, neq, k, Kt, na, i, c, p

	n = rows(X); neq = cols(ki); k = cols(X); Kt = sum(ki); na = cols(al)
	dl = J(1, neq, 0)
	if (cols(didx)) dl[selectindex(cix :> 0)] = th[didx]'
	_easi_g_probs(al, P, lnx, z, Q, zm, cix, ao, PH, PHI, M0, XB)
	F  = X * rowshape(th[bidx]', neq)'
	E  = F
	SE = J(n, neq + 1, .)
	Gs = J(n, Kt + na, 0)
	for (i = 1; i <= neq; i++) {
		G = J(n, Kt + na, 0)
		if (cix[i]) {
			c = cix[i]
			E[., i] = PH[., c] :* F[., i] + dl[i] :* PHI[., c]
			G[| 1, ko[i] + 1 \ n, ko[i] + k |] = PH[., c] :* X
			G[., ko[i] + k + 1] = PHI[., c]
			S = _easi_seldesign(P, lnx, z, Q, (cols(Q) ? zm[i, .] : J(1, 0, .)))
			p = cols(S)
			G[| 1, Kt + ao[c] + 1 \ n, Kt + ao[c] + p |] =
				(PHI[., c] :* (F[., i] - dl[i] :* XB[., c])) :* S
		}
		else G[| 1, ko[i] + 1 \ n, ko[i] + k |] = X
		SE[., i] = sqrt(rowsum((G * Vj) :* G))
		Gs = Gs + G
	}
	E = E, (1 :- rowsum(E))
	SE[., neq + 1] = sqrt(rowsum((Gs * Vj) :* Gs))
}

// estat engel, at the means, with selection: the latent fixed point at each
// grid point (the residual at zero), then the expected shares with the
// probits at the means of prices, demographics and probit-only variables
void _easi_g_engelgrid(string scalar svars, string scalar lpvars,
	string scalar lxvar, string scalar zvars, string scalar qvars,
	string scalar touse, string scalar wtvar, real scalar G, real scalar trim,
	real scalar R, real scalar py, real scalar pz, real scalar zy,
	string scalar sipz, real scalar lgQuant, string scalar outnm,
	string scalar sgstr, string scalar zmstr)
{
	real matrix s, p, z, Q, zm, zg, pg, qg, npg, X, W, E, SE
	real colvector lnx, wt, lnxg, y, ynew, pAp, pBp, th, q
	real rowvector ipz, cix, ki, ko, bidx, didx, ao, al
	real scalar n, J, neq, T, k, i, it, nipz, onp, oynp, onpz, nq

	s   = st_data(., svars,  touse)
	p   = st_data(., lpvars, touse)
	lnx = st_data(., lxvar,  touse)
	z   = st_data(., zvars,  touse)
	Q   = (qvars == "" ? J(rows(s), 0, .) : st_data(., qvars, touse))
	wt  = st_data(., wtvar,  touse)
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))
	th  = st_matrix("e(b)")'
	n = rows(s); J = cols(s); neq = J - 1; T = cols(z); nipz = cols(ipz)
	nq = cols(Q)
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	_easi_g_layout(sgstr, zmstr, neq, k, nq, cix, ki, ko, bidx, didx, ao, al, zm)

	q    = (trim :+ (100 - 2 * trim) :* (((1::G) :- 0.5) :/ G)) :/ 100
	lnxg = J(G, 1, 0)
	for (i = 1; i <= G; i++) lnxg[i] = _easi_wquant(lnx, wt, q[i])
	zg  = J(G, 1, 1) * (quadcross(wt, z) :/ quadsum(wt))
	pg  = J(G, 1, 1) * (quadcross(wt, p) :/ quadsum(wt))
	qg  = (nq ? J(G, 1, 1) * (quadcross(wt, Q) :/ quadsum(wt)) : J(G, 0, .))
	npg = pg[, 1::neq] :- pg[, J]
	_easi_quadf(th[bidx], npg, zg, R, py, pz, zy, ipz, lgQuant, pAp, pBp)
	y = lnxg
	for (it = 1; it <= 500; it++) {
		X = _easi_design(y, zg, npg, R, py, zy, pz, ipz)
		W = X * rowshape(th[bidx]', neq)'
		W = W, (1 :- rowsum(W))
		ynew = (lnxg - rowsum(W :* pg) :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)
		if (max(abs(ynew - y)) < 1e-13) {
			y = ynew
			break
		}
		y = ynew
	}
	X = _easi_design(y, zg, npg, R, py, zy, pz, ipz)
	_easi_g_ewse(X, pg, lnxg, zg, qg, zm, cix, ki, ko, ao, al, th, bidx, didx,
		st_matrix("e(V_sel)"), E, SE)
	st_matrix(outnm, (100 :* q, lnxg, y, E, SE))
}

// estat engel, as observed, with selection: the expected shares of each
// household (as predict) and their standard errors
void _easi_g_engel(string scalar svars, string scalar lpvars, string scalar lxvar,
	string scalar zvars, string scalar qvars, string scalar touse,
	string scalar wvars, string scalar sevars, real scalar R, real scalar py,
	real scalar pz, real scalar zy, string scalar sipz, real scalar lgQuant,
	string scalar sgstr, string scalar zmstr)
{
	real matrix s, p, z, Q, zm, np, D0, PH, PHI, M0, X, E, SE
	real colvector lnx, ystone, y, th
	real rowvector ipz, cix, ki, ko, bidx, didx, ao, al, dl
	real scalar n, J, neq, T, k, i, nipz, onp, oynp, onpz, nq

	s   = st_data(., svars,  touse)
	p   = st_data(., lpvars, touse)
	lnx = st_data(., lxvar,  touse)
	z   = st_data(., zvars,  touse)
	Q   = (qvars == "" ? J(rows(s), 0, .) : st_data(., qvars, touse))
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))
	th  = st_matrix("e(b)")'
	n = rows(s); J = cols(s); neq = J - 1; T = cols(z); nipz = cols(ipz)
	nq = cols(Q)
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	_easi_g_layout(sgstr, zmstr, neq, k, nq, cix, ki, ko, bidx, didx, ao, al, zm)
	np     = p[, 1::neq] :- p[, J]
	ystone = lnx - rowsum(s :* p)
	D0 = J(n, neq, 0)
	for (i = 1; i <= neq; i++) if (cix[i]) D0[., i] = (s[., i] :<= 0)
	dl = J(1, neq, 0)
	dl[selectindex(cix :> 0)] = th[didx]'
	_easi_g_probs(al, p, lnx, z, Q, zm, cix, ao, PH, PHI, M0)
	y = _easi_g_y(ystone, ystone, np, z, D0, M0, cix, th[bidx], dl, R, py, zy,
		pz, ipz, lgQuant)
	X = _easi_design(y, z, np, R, py, zy, pz, ipz)
	_easi_g_ewse(X, p, lnx, z, Q, zm, cix, ki, ko, ao, al, th, bidx, didx,
		st_matrix("e(V_sel)"), E, SE)
	st_store(., tokens(wvars),  touse, E)
	st_store(., tokens(sevars), touse, SE)
}

// predict after an estimation with selection: the expected shares
// E[w_i] = Phi_i f_i + delta_i phi_i (the last good by adding-up), or the
// latent shares f, the completed implicit utility y, the residuals w - E[w].
// The probit coefficients are read from e(sel_alpha), the latent
// coefficients and the delta from e(b).
void _easi_g_pred(string scalar svars, string scalar lpvars, string scalar lxvar,
	string scalar zvars, string scalar qvars, string scalar touse,
	string scalar newvars, real scalar R, real scalar py, real scalar pz,
	real scalar zy, string scalar sipz, real scalar lgQuant, real scalar what,
	real scalar latent, string scalar sgstr, string scalar zmstr)
{
	real matrix s, p, z, np, Q, zm, AL, PH, PHI, M0, D0, F, E
	real colvector lnx, ystone, th, y, b, col
	real rowvector ipz, cix, ki, ko, bidx, didx, ao, al, dl
	real scalar n, J, neq, T, nipz, k, C, i, c, onp, oynp, onpz, nq

	s   = st_data(., svars,  touse)
	p   = st_data(., lpvars, touse)
	lnx = st_data(., lxvar,  touse)
	z   = st_data(., zvars,  touse)
	Q   = (qvars == "" ? J(rows(s), 0, .) : st_data(., qvars, touse))
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))
	th  = st_matrix("e(b)")'
	AL  = st_matrix("e(sel_alpha)")
	n = rows(s); J = cols(s); neq = J - 1; T = cols(z); nipz = cols(ipz)
	nq = cols(Q)
	np     = p[, 1::neq] :- p[, J]
	ystone = lnx - rowsum(s :* p)
	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)

	cix = strtoreal(tokens(sgstr))
	C = 0
	for (i = 1; i <= neq; i++) {
		if (cix[i]) {
			C++
			cix[i] = C
		}
	}
	zm = (nq ? rowshape(strtoreal(tokens(zmstr)), neq) : J(neq, 0, .))
	ki = J(1, neq, k) + (cix :> 0)
	ko = J(1, neq, 0)
	for (i = 2; i <= neq; i++) ko[i] = ko[i - 1] + ki[i - 1]
	bidx = J(1, 0, .); didx = J(1, 0, .)
	ao = J(1, C, 0); al = J(1, 0, .)
	D0 = J(n, neq, 0)
	for (i = 1; i <= neq; i++) {
		bidx = bidx, (ko[i] + 1)..(ko[i] + k)
		if (!cix[i]) continue
		didx = didx, ko[i] + k + 1
		c = cix[i]
		col = AL[., i]
		col = select(col, col :< .)
		if (c < C) ao[c + 1] = ao[c] + rows(col)
		al = al, col'
		D0[., i] = (s[., i] :<= 0)
	}
	b  = th[bidx]
	dl = J(1, neq, 0)
	dl[selectindex(cix :> 0)] = th[didx]'
	_easi_g_probs(al, p, lnx, z, Q, zm, cix, ao, PH, PHI, M0)
	y = _easi_g_y(ystone, ystone, np, z, D0, M0, cix, b, dl, R, py, zy, pz,
		ipz, lgQuant)
	if (what == 2) {
		st_store(., tokens(newvars), touse, y)
		return
	}
	F = _easi_design(y, z, np, R, py, zy, pz, ipz) * rowshape(b', neq)'
	E = F
	if (!latent) {
		for (i = 1; i <= neq; i++) {
			if (!cix[i]) continue
			E[., i] = PH[., cix[i]] :* F[., i] + dl[i] :* PHI[., cix[i]]
		}
	}
	E = E, (1 :- rowsum(E))
	if (what == 1) st_store(., tokens(newvars), touse, E)
	else           st_store(., tokens(newvars), touse, s - E)
}

void _easi_pred(string scalar svars, string scalar lpvars, string scalar lxvar,
	string scalar zvars, string scalar touse, string scalar newvars,
	real scalar R, real scalar py, real scalar pz, real scalar zy,
	string scalar sipz, real scalar lgQuant, real scalar what)
{
	real matrix s, p, z, np, X, B, Bfull, Shat, Aco, Bco
	real colvector lnx, ystone, y, pAp, pBp, b
	real rowvector ipz
	real scalar n, Jg, neq, T, k, i, j, t, nipz, onp, oynp, onpz

	s   = st_data(., svars,  touse)
	p   = st_data(., lpvars, touse)
	lnx = st_data(., lxvar,  touse)
	z   = st_data(., zvars,  touse)
	ipz = (sipz == "" ? J(1, 0, 0) : strtoreal(tokens(sipz)))
	b   = st_matrix("e(b)")'

	n = rows(s); Jg = cols(s); neq = Jg - 1; T = cols(z); nipz = cols(ipz)

	np     = p[, 1::neq] :- p[, Jg]
	ystone = lnx - rowsum(s :* p)

	onp  = 1 + R + T + (zy ? T : 0)
	oynp = onp + neq
	onpz = oynp + (py ? neq : 0)
	k    = onpz + (pz ? nipz * neq : 0)
	B    = rowshape(b', neq)'		// k x neq

	// the implicit utility index, straight from the coefficients
	Aco = B[| onp + 1, 1 \ onp + neq, neq |]'
	pAp = J(n, 1, 0)
	for (i = 1; i <= neq; i++) {
		for (j = 1; j <= neq; j++) pAp = pAp + Aco[i, j] :* np[, i] :* np[, j]
	}
	if (pz) {
		for (t = 1; t <= nipz; t++) {
			Aco = B[| onpz + (t-1)*neq + 1, 1 \ onpz + t*neq, neq |]'
			for (i = 1; i <= neq; i++) {
				for (j = 1; j <= neq; j++) {
					pAp = pAp + Aco[i, j] :* np[, i] :* np[, j] :* z[, ipz[t]]
				}
			}
		}
	}
	pBp = J(n, 1, 0)
	if (py) {
		Bco = B[| oynp + 1, 1 \ oynp + neq, neq |]'
		for (i = 1; i <= neq; i++) {
			for (j = 1; j <= neq; j++) pBp = pBp + Bco[i, j] :* np[, i] :* np[, j]
		}
	}
	if (lgQuant) {
		pAp = _easi_quant(pAp)
		pBp = _easi_quant(pBp)
	}
	y = (ystone :+ 0.5 :* pAp) :/ (1 :- 0.5 :* pBp)

	if (what == 2) {
		st_store(., tokens(newvars), touse, y)
		return
	}

	X     = _easi_design(y, z, np, R, py, zy, pz, ipz)
	Shat  = X * B
	Bfull = Shat, (1 :- rowsum(Shat))	// the last good by adding up

	if (what == 1) st_store(., tokens(newvars), touse, Bfull)
	else           st_store(., tokens(newvars), touse, s - Bfull)
}

end
