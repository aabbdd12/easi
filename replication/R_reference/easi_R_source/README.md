# easi 0.21 — R source code (recovered)

The CRAN binary build shipped in `easi/easi/` contains no source: the code lives in
the lazy-load database `R/easi.rdb` + `R/easi.rdx`. These files were recovered with

```r
e <- new.env(); lazyLoad("easi/easi/R/easi", envir = e)
for (n in ls(e)) writeLines(deparse(get(n, envir = e), width.cutoff = 200), paste0(n, ".R"))
```

so they are the *deparsed* (re-printed) bodies of the installed functions, not the
original `.R` files: formatting and comments are lost, semantics are not.

## The one patch applied

`easi.R` builds each equation's formula as

```r
formula(paste(paste("eqS", i, sep = ""), "<-", paste("s", i, sep = ""), "~", form6))
```

i.e. it hands `formula()` a string containing an **assignment**,
`"eqS1 <- s1 ~ + y1 + ..."`. `formula.character()` accepted that up to R 3.x; under
R >= 4.0 it fails with

```
Erreur dans class(ff) <- "formula" : tentative de changer un attribut en NULL
```

The line was replaced by

```r
as.formula(paste(paste("s", i, sep = ""), "~", form6))
```

The equation list is unnamed either way, so `systemfit` still labels the equations
`eq1`, `eq2`, ... — which is what the coefficient names (`eq1_np2`, ...) that the rest
of the package looks up depend on. **Estimates are unaffected**; this only makes the
code loadable on a modern R.

## Loading the package under R >= 4.0

`library(easi)` refuses to load the 2018 binary build ("package installed before
R 4.0.0"). Source these files instead, after loading the two real dependencies:

```r
library(systemfit); library(micEcon)
for (f in list.files("easi_R_source", full.names = TRUE)) source(f)
```

Verified working on R 4.3.0.
