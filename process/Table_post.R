data$value[data$var == "repro.age"] <- as.numeric(data$value[data$var == "repro.age"]) * 365
data$units[data$var == "repro.age"] <- "days"

# correct error in var entry for Branta sandvicensis

# correct data errors:
data <- data[!(data$species == "Phasianus_colchicus" & data$var == "fecundity"),]
data <- data[!(data$species == "Passer_domesticus" & data$var == "habitat"),]

# remove extinct species:
data <- data[!data$species == "Conuropsis_carolinensis",] #species extinct
data <- data[!data$species == "Podilymbus_gigas",] #species extinct
data <- data[!data$species == "Xenicus_longipes",] #species extinct
data <- data[!data$species == "Ectopistes_migratorius",] #species extinct
data <- data[!data$species == "Pezophaps_solitaria",] #species extinct
data <- data[!data$species == "Pinguinus_impennis",] #species extinct
data <- data[!data$species == "Porphyrio_albus",] #species extinct
