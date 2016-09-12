rm(list = ls())
wkf = "mgm"
param = "phylonetworker.R"
file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
source(file_setup_path)
source("~/Documents/workflows/bird_trait_networks/project_ui.R")



test.params <- expand.grid(vars = seq(0.2, 1, by = 0.2), spp = seq(0.2, 1, by = 0.2))

error <- testImpute(data, test.params, ms_vars, metadata, meta_types, mgm_types)

ggplot(error) + geom_point(aes(x = s.p, y = NRMSE, colour = factor(v.p))) 
ggplot(error) + geom_point(aes(x = v.p, y = NRMSE, colour = factor(s.p))) 


err.lm <- lm(NRMSE ~ s.p * v.p, data = data.frame(error))
summary(err.lm)
plot(err.lm)


load(file = paste(output.folder, "data/imputed_data/", an.ID, 
                  "-v", test.params$vars[i], 
                  "-s", test.params$spp[i],
                  ".Rdata", sep = ""))
