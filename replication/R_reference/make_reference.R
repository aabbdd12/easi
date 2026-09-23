# ---------------------------------------------------------------------------
# EASI - reference values produced by the R package easi 0.21
#
# Single reference example: hixdata, full specification of the package vignette
# (the one that replicates Lewbel & Pendakur, AER 2009):
#
#   9 goods, y.power = 5, 5 demographics,
#   py.inter = TRUE, zy.inter = TRUE, pz.inter = TRUE, interpz = 1:5
#
# Everything the Stata port has to match is written to ../out/ as plain CSV
# (full double precision, %.17g) plus hixdata.dta for Stata.
#
# Run:  Rscript make_reference.R
# ---------------------------------------------------------------------------

suppressMessages({library(systemfit); library(micEcon); library(haven)})

# script dir: --file= under Rscript, working dir otherwise
.a    <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here  <- if (length(.a)) dirname(normalizePath(sub("^--file=", "", .a[1]))) else getwd()
root  <- normalizePath(file.path(here, "..", ".."))
outd  <- file.path(root, "tests", "out")
datd  <- file.path(root, "tests", "data")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)
dir.create(datd, showWarnings = FALSE, recursive = TRUE)

for (f in list.files(file.path(here, "easi_R_source"), pattern = "\\.R$", full.names = TRUE))
  source(f)

load(file.path(root, "easi", "easi", "data", "hixdata.rda"))
d <- hixdata

# ---- the frozen estimation sample -----------------------------------------
# hixdata has no missing values and its shares sum to 1 (up to ~5e-8: the source
# data was stored in single precision), so the estimation sample IS the file:
# no listwise deletion, no subsampling, no RNG anywhere.
stopifnot(!anyNA(d), max(abs(rowSums(d[, 2:10]) - 1)) < 1e-6)
cat("max |sum(shares) - 1| =", sprintf("%.3e", max(abs(rowSums(d[, 2:10]) - 1))), "\n")

shares    <- d[, 2:10]                    # sfoodh ... spers
log.price <- d[, 11:19]                   # pfoodh ... ppers  (already logs)
log.exp   <- d[, 20]                      # log_y
var.soc   <- d[, 21:25]                   # age hsex carown time tran

labels.share <- c("foodh", "foodr", "rent", "oper", "furn",
                  "cloth", "tranop", "recr", "pers")
labels.soc   <- c("age", "hsex", "carown", "time", "tran")

# hixdata for Stata (values are bit-identical, both sides read the same numbers)
write_dta(d, file.path(datd, "hixdata.dta"), version = 13)

# ---- the reference specification ------------------------------------------
SPEC <- list(y.power = 5, py.inter = TRUE, zy.inter = TRUE, pz.inter = TRUE,
             interpz = 1:ncol(var.soc))
cat("spec: y.power=", SPEC$y.power, " py=", SPEC$py.inter,
    " zy=", SPEC$zy.inter, " pz=", SPEC$pz.inter,
    " interpz=1:", ncol(var.soc), "\n", sep = "")

t0  <- Sys.time()
est <- easi(shares = shares, log.price = log.price, var.soc = var.soc,
            log.exp = log.exp, y.power = SPEC$y.power,
            labels.share = labels.share, labels.soc = labels.soc,
            py.inter = SPEC$py.inter, zy.inter = SPEC$zy.inter,
            pz.inter = SPEC$pz.inter, interpz = SPEC$interpz)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("estimation elapsed:", round(secs, 1), "s\n")

ei <- elastic(est, type = "income",       sd = TRUE)
ep <- elastic(est, type = "price",        sd = TRUE)
ez <- elastic(est, type = "demographics", sd = TRUE)

# ---- dump ------------------------------------------------------------------
w17 <- function(x, file, rn = TRUE) {
  x <- as.matrix(x)
  write.table(format(x, digits = 17, scientific = TRUE, trim = TRUE),
              file.path(outd, file), sep = ",", quote = FALSE,
              row.names = rn, col.names = if (rn) NA else TRUE)
}

cf <- coef(est)                                   # Estimate / Std.Err / t / p
w17(cf,                     "R_coef.csv")
w17(vcov(est),              "R_vcov.csv")
w17(ei$ELASTINCOME,         "R_elast_income.csv")
w17(ei$ELASTINCOME_SE,      "R_elast_income_se.csv")
w17(ep$ELASTPRICE,          "R_elast_price.csv")
w17(ep$ELASTPRICE_SE,       "R_elast_price_se.csv")
w17(ez$EZ,                  "R_elast_demo.csv")
w17(ez$EZ_SE,               "R_elast_demo_se.csv")
w17(ei$ER,                  "R_semi_income.csv")     # semi-elasticities
w17(ep$EP,                  "R_semi_price.csv")
w17(ep$EPS,                 "R_slutsky.csv")
w17(ep$EPQ,                 "R_compensated_q.csv")
w17(cbind(y = est$y),       "R_y.csv",       rn = FALSE)   # implicit utility
w17(predict(est),           "R_fitted.csv",  rn = FALSE)   # fitted shares

meta <- data.frame(
  n          = length(log.exp),
  ngoods     = ncol(shares),
  neq        = est$neq,
  nsoc       = est$nsoc,
  y_power    = est$y.power,
  py_inter   = SPEC$py.inter,
  zy_inter   = SPEC$zy.inter,
  pz_inter   = SPEC$pz.inter,
  ncoef      = nrow(cf),
  k_per_eq   = nrow(cf) / est$neq,
  elapsed_s  = round(secs, 2),
  R_version  = paste(R.version$major, R.version$minor, sep = "."),
  systemfit  = as.character(packageVersion("systemfit")),
  stamp      = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write.csv(meta, file.path(outd, "R_meta.csv"), row.names = FALSE)

cat("\n--- reference written to", outd, "---\n")
print(t(meta))
cat("\nincome elasticities:\n"); print(round(ei$ELASTINCOME, 6))
