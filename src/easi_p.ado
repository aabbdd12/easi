*! version 0.2.0  21sep2026
*! predict after easi -- relay only
*!
*! The real work is -program Predict- inside easi.ado: Mata functions written in
*! an ado-file are private to that file, so the prediction code has to live
*! beside the routines it uses.  Keeping it there is what lets the module ship
*! as plain ado-files, with no Mata library to build.

program easi_p
	version 14.2
	if !inlist("`e(cmd)'", "easi", "sr_easi") error 301
	easi _predict `0'
end
