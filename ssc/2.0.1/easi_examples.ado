*! easi_examples 2.0.1  30sep2026  Abdelkrim Araar
*! The examples of help easi and help easidiag, run from their links.
*!   easi_examples #          run example # in the command window
*!   easi_examples #, db      open the dialog box of easi filled in for example #
*!   easi_examples #, do      open example # as a do-file in the Do-file Editor
*! The data in memory are never lost: the run keeps them (preserve) and gives
*! them back at the end, even after an error or a Break; the do-file does the
*! same; the dialog box, which needs the example data in memory, refuses to
*! replace data of the user that have unsaved changes (the example data it
*! loads are marked and can be replaced). Files written by the examples go to
*! Stata's temporary folder, c(tmpdir), never to the working folder.
*! The example data are ancillary files: they are read from the current
*! folder (where "ssc install easi, all" or "net get easi" copies them), else
*! from the SSC archive, else from GitHub (_easi_exload); nothing is written.
program define easi_examples
	version 14.2
	syntax anything(name=ex id="example number") [, DB DO NOEDIT]
	capture confirm integer number `ex'
	if _rc | !inrange(`ex', 1, 9) {
		di as err "easi_examples: the examples are numbered 1 to 9 (see help easi)"
		exit 198
	}
	if "`db'" != "" & "`do'" != "" {
		di as err "easi_examples: db and do cannot be combined"
		exit 198
	}
	local T = c(tmpdir)
	local T : subinstr local T "\" "/", all
	if substr("`T'", -1, 1) != "/" local T "`T'/"
	* the Canadian data of Lewbel and Pendakur (2009), prices and expenditure
	* in logarithms, and the Mexican cereals (ENIGH 2014) with their design,
	* both ancillary files of the package (see _easi_exload)
	local SH "sfoodh sfoodr srent soper sfurn scloth stranop srecr spers"
	local PR "pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers"
	local H "`SH', lnprices(`PR') lnexpenditure(log_y)"
	local HZ "age hsex carown time tran"
	local M "w1 w2 w3, lnprices(lp1 lp2 lp3) lnexpenditure(lx)"
	local data "hixdata"
	local n 0
	if `ex' == 1 {
		local title "The Canadian data of Lewbel and Pendakur (2009)"
		local c1 "easi `H' demographics(age hsex carown) power(3)"
		local c2 "easi, compensated checks stars"
		local n 2
	}
	else if `ex' == 2 {
		local title "The specification of Lewbel and Pendakur (2009): power 5, all interactions"
		local c1 "easi `H' demographics(`HZ') power(5) py zy pz"
		local n 1
	}
	else if `ex' == 3 {
		local data "mex_bench"
		local title "Survey design (Mexican cereals): the variance of the design"
		local c1 "svyset"
		local c2 "easi `M' demographics(z1 z2) power(3) py vce(svy) stars"
		local n 2
	}
	else if `ex' == 4 {
		local title "Reproduce the R package easi 0.21 (compat)"
		local c1 "easi `H' demographics(`HZ') power(5) py zy pz compat"
		local n 1
	}
	else if `ex' == 5 {
		local title "After estimation: fitted shares, the implicit utility, Engel curves"
		local c1 "easi `H' demographics(age hsex carown) power(3) notable"
		local c2 "predict double what*, shares"
		local c3 "predict double yhat, y"
		local c4 "summarize yhat what1-what3"
		local c5 `"estat engel, n(60) data("`T'easi_curves")"'
		local n 5
	}
	else if `ex' == 6 {
		local title "The tables in a file (Word, Excel, LaTeX, CSV, Markdown)"
		local c1 `"easi `H' demographics(age hsex carown) power(3) compensated notable saveres("`T'easi_tables.docx")"'
		local c2 `"easi, stars saveres("`T'easi_tables.xlsx") notable"'
		local n 2
	}
	else if `ex' == 7 {
		local data "mex_bench"
		local title "Diagnose a specification before estimating it"
		local c1 "easidiag w1 w2 w3, lnprices(lp1_raw lp2_raw lp3) lnexpenditure(lx_raw) demographics(z1 z2) power(3)"
		local c2 "easidiag `M' demographics(z1 z2) power(3)"
		local n 2
	}
	else if `ex' == 8 {
		local data "mex_bench"
		local title "The elasticities of the households, of the market and of the individuals"
		local c1 "easi `M' demographics(z1 z2) power(3) py vce(svy)"
		local c2 "easi, elasticities(market)"
		local c3 "easi `M' demographics(z1 z2) power(3) py vce(svy) hhsize(hhsize)"
		local n 3
	}
	else if `ex' == 9 {
		local data "mex_bench"
		local title "The households that do not buy: the selection of the buyers"
		local SV "selvars(age isMale)"
		local c1 "easidiag w1 w2 w3 [pw = sweight], lnprices(lp1 lp2 lp3) lnexpenditure(lx) demographics(z1 z2) power(3) `SV'"
		local c2 "easi `M' demographics(z1 z2) power(3) vce(svy) `SV'"
		local c3 "predict double Ew*, shares"
		local c4 "predict double f*, shares latent"
		local c5 "summarize w1 Ew1 f1 w2 Ew2 f2"
		local c6 "easi `M' demographics(z1 z2) power(3) `SV' vce(bootstrap, reps(50) seed(1) svy)"
		local n 6
	}

	* ---- as a do-file, in Stata's temporary folder ----
	if "`do'" != "" {
		local fn "`T'easi_example_`ex'.do"
		tempname fh
		file open `fh' using "`fn'", write text replace
		file write `fh' "* easi, example `ex': `title'" _n
		file write `fh' "* Written by easi_examples in Stata's temporary folder; save it elsewhere to keep it." _n
		file write `fh' "* preserve keeps the data in memory and gives them back when this do-file ends;" _n
		file write `fh' "* delete that line to keep working on the example data." _n
		file write `fh' "preserve" _n
		file write `fh' "* the example data: the current folder (where ssc install easi, all copies" _n
		file write `fh' "* them), else the SSC archive, else GitHub" _n
		local l = substr("`data'", 1, 1)
		file write `fh' `"capture use `data', clear"' _n
		file write `fh' `"if _rc capture use "http://fmwww.bc.edu/repec/bocode/`l'/`data'.dta", clear"' _n
		file write `fh' `"if _rc use "https://raw.githubusercontent.com/aabbdd12/easi/main/examples/`data'.dta", clear"' _n
		forvalues i = 1/`n' {
			file write `fh' `"`c`i''"' _n
		}
		file close `fh'
		if "`noedit'" == "" doedit "`fn'"
		di as txt "(example `ex' written to " as res `"`fn'"' as txt ")"
		exit
	}

	* ---- in the dialog box: needs the example data in memory ----
	* The dialog takes prices and expenditure in levels: they are built from
	* the logs of the example data (exp(), as examples/hixdata_for_dialog.do
	* does), which gives the same estimates.  Example 5 is predict and estat
	* engel after an estimation: commands, not a dialog.
	if "`db'" != "" {
		if `ex' == 5 {
			di as err "easi_examples: example 5 is predict and estat engel after an estimation; run it in the command window"
			exit 198
		}
		local isex : char _dta[easi_example]
		if c(changed) & "`isex'" != "1" {
			di as err "easi_examples, db: the data in memory have changes not saved;"
			di as err "save them (or clear) first: the dialog box needs the example data in memory"
			exit 4
		}
		_easi_exload `data'
		* the logs of the example, and their levels for the dialog
		if inlist(`ex', 1, 2, 4, 6) {
			local LP `PR'
			local LX log_y
			local SHD `SH'
		}
		else {
			local LP lp1 lp2 lp3
			local LX lx
			if `ex' == 7 {
				local LP lp1_raw lp2_raw lp3
				local LX lx_raw
			}
			local SHD w1 w2 w3
		}
		local EP ""
		foreach v of local LP {
			qui gen double e`v' = exp(`v')
			label variable e`v' "price (level), exp(`v')"
			local EP `EP' e`v'
		}
		qui gen double expend = exp(`LX')
		label variable expend "total expenditure (level), exp(`LX')"
		char _dta[easi_example] "1"
		di as txt "(example data loaded for the dialog box; prices and expenditure in levels:"
		di as txt " the exponentials of `LP' and `LX': `EP' expend)"
		db easi
		* Stata keeps the state of a dialog between two openings: every control
		* an example may set is first put back to its default
		.easi_dlg.main.cb_action.setvalue "est"
		.easi_dlg.main.name_snames.setvalue ""
		.easi_dlg.main.vl_inddemo.setvalue ""
		.easi_dlg.main.sp_pow.setvalue 5
		.easi_dlg.main.ck_inpy.setoff
		.easi_dlg.main.ck_inpz.setoff
		.easi_dlg.main.ck_inzy.setoff
		.easi_dlg.main.cb_vce.setvalue "robust"
		.easi_dlg.main.cb_elas.setvalue "default"
		.easi_dlg.main.vn_hhs.setvalue ""
		.easi_dlg.resop.ck_comp.setoff
		.easi_dlg.resop.ck_demoel.setoff
		.easi_dlg.resop.ck_checks.setoff
		.easi_dlg.resop.ck_noese.setoff
		.easi_dlg.resop.ck_stars.setoff
		.easi_dlg.resop.ck_notab.setoff
		.easi_dlg.resop.fi_save.setvalue ""
		.easi_dlg.resop.ck_compat.setoff
		.easi_dlg.main.sp_reps.setvalue 200
		.easi_dlg.main.ed_seed.setvalue ""
		.easi_dlg.main.ck_bsvy.setoff
		.easi_dlg.sel.vl_pimp.setvalue ""
		.easi_dlg.sel.ck_sel.setoff
		.easi_dlg.sel.ed_selg.setvalue ""
		.easi_dlg.sel.vl_all.setvalue ""
		.easi_dlg.sel.cb_ng.setvalue "0"
		* the example
		.easi_dlg.main.name_items.setvalue "`SHD'"
		.easi_dlg.main.name_prices.setvalue "`EP'"
		.easi_dlg.main.vn_hhexp.setvalue "expend"
		if `ex' == 1 {
			.easi_dlg.main.vl_inddemo.setvalue "age hsex carown"
			.easi_dlg.main.sp_pow.setvalue 3
			.easi_dlg.resop.ck_comp.seton
			.easi_dlg.resop.ck_checks.seton
			.easi_dlg.resop.ck_stars.seton
		}
		if inlist(`ex', 2, 4) {
			.easi_dlg.main.vl_inddemo.setvalue "`HZ'"
			.easi_dlg.main.sp_pow.setvalue 5
			.easi_dlg.main.ck_inpy.seton
			.easi_dlg.main.ck_inpz.seton
			.easi_dlg.main.ck_inzy.seton
		}
		if `ex' == 4 {
			.easi_dlg.resop.ck_compat.seton
		}
		if `ex' == 6 {
			.easi_dlg.main.vl_inddemo.setvalue "age hsex carown"
			.easi_dlg.main.sp_pow.setvalue 3
			.easi_dlg.resop.ck_comp.seton
			.easi_dlg.resop.ck_notab.seton
			.easi_dlg.resop.fi_save.setvalue "`T'easi_tables.docx"
		}
		if inlist(`ex', 3, 7, 8, 9) {
			.easi_dlg.main.vl_inddemo.setvalue "z1 z2"
			.easi_dlg.main.sp_pow.setvalue 3
		}
		if `ex' == 3 {
			.easi_dlg.main.ck_inpy.seton
			.easi_dlg.main.cb_vce.setvalue "svy"
			.easi_dlg.resop.ck_stars.seton
		}
		if `ex' == 7 {
			.easi_dlg.main.cb_action.setvalue "diag"
		}
		if `ex' == 8 {
			.easi_dlg.main.ck_inpy.seton
			.easi_dlg.main.cb_vce.setvalue "svy"
			.easi_dlg.main.cb_elas.setvalue "market"
		}
		if `ex' == 9 {
			.easi_dlg.main.cb_vce.setvalue "svy"
			.easi_dlg.sel.ck_sel.seton
			.easi_dlg.sel.vl_all.setvalue "age isMale"
		}
		exit
	}

	* ---- in the command window: the data in memory are kept ----
	preserve
	_easi_exload `data'
	local src "`r(source)'"
	di as txt _n "{hline 78}" _n "easi, example `ex': " as res "`title'" _n as txt "{hline 78}"
	di as txt `"(example data: `data'.dta, read from `src'; the data in memory come back at the end)"'
	forvalues i = 1/`n' {
		di as txt _n `". `c`i''"'
		`c`i''
	}
	if inlist(`ex', 5, 6) di as txt _n `"(files written to `T')"'
end

* ============================================================================
* the data of an example, an ancillary file of the package: from the current
* folder (where "ssc install easi, all" or "net get easi" copies it), else
* from the SSC archive, else from GitHub; nothing is written to disk
program define _easi_exload, rclass
	args f
	capture confirm file "`f'.dta"
	if !_rc {
		quietly use "`f'.dta", clear
		return local source "the current folder"
		exit
	}
	local l = substr("`f'", 1, 1)
	capture quietly use "http://fmwww.bc.edu/repec/bocode/`l'/`f'.dta", clear
	if !_rc {
		return local source "the SSC archive"
		exit
	}
	capture quietly use "https://raw.githubusercontent.com/aabbdd12/easi/main/examples/`f'.dta", clear
	if !_rc {
		return local source "GitHub"
		exit
	}
	di as err "easi_examples: `f'.dta is not in the current folder (`c(pwd)'),"
	di as err "  and neither the SSC archive nor GitHub could be reached."
	di as txt "  Copy the example data into the current folder with"
	di as txt `"  {stata "ssc install easi, all replace"} (from SSC), or"'
	di as txt `"  {stata "net get easi, from(https://raw.githubusercontent.com/aabbdd12/easi/main)"} (from GitHub)."'
	exit 601
end
