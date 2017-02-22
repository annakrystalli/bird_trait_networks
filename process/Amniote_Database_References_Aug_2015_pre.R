## ---- prR1-subset ----
data <- data[data$class == "Aves",]

## ---- prR1-derive_repro.age.diff ----
add_keep.dat <- data.frame(repro.age.diff = paste(data$female_maturity_d, 
                                                  data$male_maturity_d, sep = "; "))
add_keep.dat[is.na(data$female_maturity_d) | is.na(data$male_maturity_d),] <- NA
add_keep.vars <- "repro.age.diff"
