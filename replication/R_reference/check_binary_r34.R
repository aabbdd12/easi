# ---------------------------------------------------------------------------
# Confirm the frozen reference with the ORIGINAL COMPILED PACKAGE.
#
# make_reference.R produces tests/out/R_*.csv by sourcing the package's own
# code under R 4.3.0, because the easi 0.21 binary does not load under R >= 4.
# This script runs the binary itself -- easi_0.21.zip, unmodified -- under
# R 3.4.4 (the R that sr_easi shipped with; installer in ../../easi/) and
# compares every table with the frozen reference.  It also shows two of the
# package's behaviours directly: interpz collapsed to its first element, and
# the negative standard error of the dropped good.
#
# Requirements: R 3.4.4 and Windows binaries of systemfit/micEcon/haven for
# R 3.4.  CRAN no longer serves them; the Posit Package Manager snapshot of
# 2018-06-01 still does, and that is where they are taken from, into a
# private library so that nothing else on the machine is touched.
#
# Run:  "C:\Program Files\R\R-3.4.4\bin\x64\Rscript.exe" check_binary_r34.R
#
# Result on 2026-09-22 (R 3.4.4, easi 0.21, systemfit 1.1.22):
#   coefficients      2.6e-09      elast_income     2.6e-11   _se  1.2e-13
#   elast_price       1.2e-09      _se  5.2e-12
#   elast_demo        6.1e-12      _se  2.8e-15
#   fit$interpz stored = 1 (1:5 requested); 40 coefficients per equation,
#   only z1*ps1..ps8 among the price x demographic terms.
#   income SE of the last good: -1.019 (the sign defect, from the binary).
# ---------------------------------------------------------------------------

cat(R.version.string, "\n")
.a   <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(.a)) dirname(normalizePath(sub("^--file=", "", .a[1]))) else getwd()
root <- normalizePath(file.path(here, "..", ".."))
out  <- file.path(root, "tests", "out")

lib <- file.path(Sys.getenv("LOCALAPPDATA"), "easi_r34_lib")
dir.create(lib, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
repo <- "https://packagemanager.posit.co/cran/2018-06-01"
need <- setdiff(c("systemfit", "micEcon", "haven"), rownames(installed.packages()))
if (length(need)) install.packages(need, lib = lib, repos = repo, type = "win.binary", quiet = TRUE)
if (!"easi" %in% rownames(installed.packages()))
  install.packages(file.path(root, "easi", "easi_0.21.zip"), lib = lib, repos = NULL,
                   type = "win.binary", quiet = TRUE)
suppressMessages(library(easi))
cat("easi", as.character(packageVersion("easi")), "| systemfit",
    as.character(packageVersion("systemfit")), "\n")

# ---- the reference specification, exactly as in make_reference.R ---------
data(hixdata); d <- hixdata
shares <- d[, 2:10]; log.price <- d[, 11:19]; log.exp <- d[, 20]; var.soc <- d[, 21:25]
t0 <- Sys.time()
est <- easi(shares = shares, log.price = log.price, var.soc = var.soc, log.exp = log.exp,
            y.power = 5, py.inter = TRUE, zy.inter = TRUE, pz.inter = TRUE, interpz = 1:5)
cat("time:", format(Sys.time() - t0), "\n")

# ---- E1: interpz collapsed to its first element ---------------------------
nm <- rownames(coef(est))
cat("interpz requested 1:5 ; est$interpz stored =", est$interpz, "\n")
cat("coefficients per equation:", sum(startsWith(nm, "eq1_")), "\n")
cat("eq1 price x demographic terms:",
    nm[startsWith(nm, "eq1_z") & grepl("*ps", nm, fixed = TRUE)], "\n")

# ---- every table against the frozen reference ----------------------------
rd  <- function(f) as.matrix(read.csv(file.path(out, f))[, -1])
cmp <- function(lbl, a, f) {
  a <- as.numeric(as.matrix(a)); r <- as.numeric(rd(f))
  cat(sprintf("  %-18s max |binary - frozen| = %s\n", lbl, format(max(abs(sort(a) - sort(r))), digits = 3)))
}
ei <- elastic(est, type = "income",       sd = TRUE)
ep <- elastic(est, type = "price",        sd = TRUE)
ez <- elastic(est, type = "demographics", sd = TRUE)
cmp("coef table",      coef(est),         "R_coef.csv")   # 320 x (est, se, t, p)
cmp("elast_income",    ei$ELASTINCOME,    "R_elast_income.csv")
cmp("elast_income_se", ei$ELASTINCOME_SE, "R_elast_income_se.csv")
cmp("elast_price",     ep$ELASTPRICE,     "R_elast_price.csv")
cmp("elast_price_se",  ep$ELASTPRICE_SE,  "R_elast_price_se.csv")
cmp("elast_demo",      ez$EZ,             "R_elast_demo.csv")
cmp("elast_demo_se",   ez$EZ_SE,          "R_elast_demo_se.csv")

# ---- L6: the negative standard error, straight from the binary ------------
cat("income SE of the last good (binary):",
    format(as.numeric(as.matrix(ei$ELASTINCOME_SE))[9], digits = 5), "\n")
