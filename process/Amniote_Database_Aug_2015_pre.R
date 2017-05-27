
### subset data by class to aves and process
data <- data[data$class == "Aves",]

# correct typo in "male_maturity_d"
data[paste(data$genus, "_", data$species, sep = "") == "Rhea_americana","male_maturity_d"] <- data[paste(data$genus, "_", data$species, sep = "") == "Rhea_americana","male_maturity_d"]/10

add_keep.dat <- data.frame(repro.age.diff = log(data$male_maturity_d/data$female_maturity_d))
add_keep.vars <- "repro.age.diff"











