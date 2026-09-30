*! version 1.0.0  23sep2026
*! License: GPL-3.0-or-later (see LICENSE); https://github.com/aabbdd12/easi
*! Pre-estimation diagnostic for the EASI demand system
*! Araar Abdelkrim

* Relay, exactly as easi_p.ado does for -predict-.  The Mata routines that do
* the work live inside easi.ado and are private to it, so the diagnostic runs
* as a hidden subcommand rather than duplicating them.
program easidiag, rclass
	version 14.2
	easi _diag `0'
	return add
end
