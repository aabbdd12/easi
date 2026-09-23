*! version 1.0.0  23sep2026
*! License: GPL-3.0-or-later (see LICENSE); https://github.com/aabbdd12/easi
*! estat after easi

program easi_estat, rclass
	version 14.2
	if !inlist("`e(cmd)'", "easi", "sr_easi") error 301

	gettoken sub rest : 0, parse(" ,")
	local l = max(3, length("`sub'"))

	if "`sub'" != "" & "`sub'" == substr("engel", 1, `l') {
		easi _engel `rest'
		return add
		exit
	}
	di as error "invalid subcommand `sub'; only {bf:engel} is available"
	exit 198
end
