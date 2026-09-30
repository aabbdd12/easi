# easi 2.0.1 -- the copy submitted to the SSC archive

The same code as easi 2.0.0 on GitHub; only the version stamps differ. It is
here to test the SSC package before the archive publishes it. As on SSC, the
package file lists every file with an `f` line: `net install` installs the
programs, the help and the dialog only; the example data (`hixdata.dta`,
`mex_bench.dta`), the tours and the technical note are ancillary files that
`net get` copies to the current folder, never to the system directories.

```stata
net install easi, from("https://raw.githubusercontent.com/aabbdd12/easi/main/ssc/2.0.1") replace
net get easi, from("https://raw.githubusercontent.com/aabbdd12/easi/main/ssc/2.0.1") replace
```

The examples of the help read the data from the current folder, else from the
SSC archive, else from GitHub. To go back to the GitHub version:

```stata
net install easi, from("https://raw.githubusercontent.com/aabbdd12/easi/v2.0.0") replace
```
