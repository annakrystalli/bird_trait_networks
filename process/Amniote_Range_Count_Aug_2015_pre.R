## ---- prN1-subset ----
data <- data[data$classx == "Aves",]

## ---- prN1-derive_repro.age.diff ----
add_keep.dat <- data.frame(repro.age.diff = data$count_female_maturity + data$count_male_maturity)
data[data == 0] <- NA
add_keep.vars <- "repro.age.diff"
