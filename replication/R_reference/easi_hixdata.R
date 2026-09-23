# ---------------------------------------------------------------------------
# The reference model on hixdata with the R package easi 0.21 -- nothing else.
#
# This is the model every comparison in the note is made against: the
# specification of the package vignette (easi/doc/easi.rnw, l. 853-872),
# which replicates Lewbel and Pendakur (2009):
#     9 goods, y.power = 5, 5 demographics, py.inter = zy.inter = pz.inter = TRUE,
#     interpz = 1:5
# In Stata:
#     easi sfoodh sfoodr srent soper sfurn scloth stranop srecr spers,
#          lnprices(pfoodh pfoodr prent poper pfurn pcloth ptranop precr ppers)
#          lnexpenditure(log_y) demographics(age hsex carown time tran)
#          power(5) py pz zy compat
#
# Two ways to have the package:
#   (a) the compiled package under R 3.4.4 (the R that sr_easi shipped with),
#       installed by check_binary_r34.R into %LOCALAPPDATA%\easi_r34_lib --
#       this script looks there first;
#   (b) any R >= 4: the package does not load, source its functions instead
#       (see make_reference.R and easi_R_source/).
#
# Run:  "C:\Program Files\R\R-3.4.4\bin\x64\Rscript.exe" easi_hixdata.R
# ---------------------------------------------------------------------------

lib <- file.path(Sys.getenv("LOCALAPPDATA"), "easi_r34_lib")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
if (!suppressWarnings(suppressMessages(require(easi, quietly = TRUE)))) {
  message("the compiled package is not installed: run check_binary_r34.R once,")
  message("or source easi_R_source/*.R as make_reference.R does")
  .a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(.a)) dirname(normalizePath(sub("^--file=", "", .a[1]))) else getwd()
  suppressMessages({library(systemfit); library(micEcon)})
  for (f in list.files(file.path(here, "easi_R_source"), pattern = "\\.R$", full.names = TRUE)) source(f)
  load(file.path(here, "..", "..", "easi", "easi", "data", "hixdata.rda"))
} else {
  data(hixdata)
}
cat(R.version.string, "| systemfit", as.character(packageVersion("systemfit")), "\n")

shares    <- hixdata[, 2:10]          # sfoodh ... spers
log.price <- hixdata[, 11:19]         # pfoodh ... ppers, already logs
log.exp   <- hixdata[, 20]            # log_y
var.soc   <- hixdata[, 21:25]         # age hsex carown time tran

t0  <- Sys.time()
est <- easi(shares = shares, log.price = log.price, var.soc = var.soc,
            log.exp = log.exp, y.power = 5,
            py.inter = TRUE, zy.inter = TRUE, pz.inter = TRUE, interpz = 1:5,
            labels.share = c("foodh", "foodr", "rent", "oper", "furn", "cloth", "tranop", "recr"),
            labels.soc   = c("age", "hsex", "carown", "time", "tran"))
cat("estimation time:", format(Sys.time() - t0), "\n")
cat("coefficients per equation:", sum(startsWith(rownames(coef(est)), "eq1_")),
    "(40 with interpz = 1:5: the package interacts prices with the first demographic only)\n")

ei <- elastic(est, type = "income",       sd = TRUE)
ep <- elastic(est, type = "price",        sd = TRUE)
ez <- elastic(est, type = "demographics", sd = TRUE)

cat("\nexpenditure elasticities\n");        print(round(ei$ELASTINCOME, 4))
cat("their standard errors\n");            print(round(ei$ELASTINCOME_SE, 4))
cat("\nuncompensated price elasticities, [price, good]\n"); print(round(ep$ELASTPRICE, 4))
cat("\ndemographic elasticities, [demographic, good]\n");    print(round(ez$EZ, 4))
