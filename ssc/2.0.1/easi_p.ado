*! version 1.0.0  23sep2026
*! License: GPL-3.0-or-later (see LICENSE); https://github.com/aabbdd12/easi
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
