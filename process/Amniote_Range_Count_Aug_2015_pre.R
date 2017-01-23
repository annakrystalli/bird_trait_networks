## ---- prN1-subset ----
data <- data[data$class == "Aves",]

## ---- prN1-derive_repro.age.diff ----
data$repro.age.diff <- data$female_maturity + data.n$male_maturity
data[data == 0] <- NA

add_keep.vars <- "repro.age.diff"
