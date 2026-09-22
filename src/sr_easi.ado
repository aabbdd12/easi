*! version 2.0.0  21sep2026
*! sr_easi -- compatibility wrapper around -easi-
*!
*! The original sr_easi wrote an R script, shelled out to R to run the "easi"
*! R package, and read the results back.  R is no longer involved: everything is
*! computed by -easi- in Stata/Mata.  This wrapper exists so that do-files,
*! saved dialog projects (*.easi) and the WELCOM menu written for the old
*! command keep working unchanged.
*!
*! It calls -easi- in -compat- mode, which reproduces the R package bit for bit,
*! including its bugs.  New work should call -easi- directly; see -help easi-.

program sr_easi, eclass
	version 14.2

	#delimit ;
	syntax varlist(min=2 numeric) [if] [in] [aweight fweight pweight iweight] , [
		EXPenditure(varlist min=1 max=1 numeric)
		LNEXPenditure(varlist min=1 max=1 numeric)
		PRices(varlist numeric)
		LNPRices(varlist numeric)
		demographics(varlist numeric)
		RTool(string)
		power(int 5)
		INPY(int 0)
		INPZ(int 0)
		INZY(int 0)
		DEC(int 4)
		SNames(string)
		INISave(string)
		xfil(string)
		dislas(int 1)
		dregres(int 0)
		CORRected
		VCE(string)
	];
	#delimit cr

	if "`rtool'" != "" {
		di as txt "note: {bf:rtool()} is obsolete and ignored -- "	///
			  "R is no longer used"
	}

	* save the dialog inputs, exactly as before, if the helper is installed
	if "`inisave'" != "" {
		capture asdbsave_easi `0'
	}

	* -compat- unless the user explicitly asks for the corrected behaviour
	local mode compat
	if "`corrected'" != "" local mode

	local wgt
	if "`weight'" != "" local wgt [`weight'`exp']

	* pass an option through only when it has a value: -easi- rejects, say, an
	* empty expenditure()
	local opts
	foreach o in expenditure lnexpenditure prices lnprices demographics	///
		     snames vce {
		if `"``o''"' != "" local opts `opts' `o'(``o'')
	}

	#delimit ;
	easi `varlist' `if' `in' `wgt' , `opts'
		power(`power') inpy(`inpy') inpz(`inpz') inzy(`inzy')
		dec(`dec') dislas(`dislas') dregres(`dregres') `mode' ;
	#delimit cr

	ereturn local cmd "sr_easi"

	* Excel export, as before, if the WELCOM helper is installed
	if "`xfil'" != "" {
		tempname EI EP
		matrix `EI' = e(elast_income)
		matrix `EP' = e(elast_price)
		tokenize `varlist'
		capture noisily {
			mk_xtab_tr `1', matn(`EI') dec(`dec') xfil(`xfil')	///
				xshe(Table_02)					///
				xtit("Table 02: Expenditure elasticities")	///
				xlan(en) dste(0)
			mk_xtab_tr `1', matn(`EP') dec(`dec') xfil(`xfil')	///
				xshe(Table_03)					///
				xtit("Table 03: Price elasticities")		///
				xlan(en) dste(0)
		}
		if _rc {
			di as txt "note: Excel export skipped "		///
				  "({bf:mk_xtab_tr} not installed)"
		}
	}
end
