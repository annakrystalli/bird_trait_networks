rm(list = ls())
wkf = "mgm"
param = "var_impute_test.R"
file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
source(file_setup_path)
source(paste0(input.folder, "project_ui.R"))
source(paste0(input.folder, "test_impute_functions.R"))

par.pkgs <- c("foreach", "doParallel", "doMC")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(par.pkgs, character.only = T)

cores <- detectCores() - 1
registerDoMC(cores=cores)

data.t <- data[1:500, ms_vars[1:50]]

# names(data)[1:2]

t0 <- Sys.time()
error.var <- foreach(x = c("all", ms_vars[2:10]), .combine = rbind,
                             .inorder = F, .errorhandling = "remove") %dopar%{
                               impute_varerror(x, data.t)}

Sys.time() - t0


rownames(data) <- data$species

test.params <- expand.grid(vars = seq(0.2, 1, by = 0.2), spp = seq(0.2, 1, by = 0.2))

v.p <- 0.1
s.p <- 0.05

for(i in 1:nrow(test.params)){
  v.p <- test.params$vars[i]
  s.p <- test.params$spp[i]

error.s <- testImpute(data, v.p = v.p, s.p = s.p, ms_vars, metadata, meta_types, mgm_types)


}



impt <- g[c("imp_data","OOBerror")]
impt.out <- c(impt, list(vars = vars, spps = spps, v.p = v.p, s.p = s.p))
save(impt.out, file = paste(output.folder, "data/imputed_data/", an.ID, 
                            "-v", test.params$vars[i], 
                            "-s", test.params$spp[i],
                            ".Rdata", sep = ""))



error.s <- error[order(error$NRMSE),]
write.csv(error.s, paste(input.folder, "r data/imp_err_", an.ID, ".csv",
                         sep = ""))


ggplot(error) + geom_point(aes(x = s.p, y = NRMSE, colour = factor(v.p))) 
ggplot(error) + geom_point(aes(x = v.p, y = NRMSE, colour = factor(s.p))) 


err.lm <- lm(NRMSE ~ s.p * v.p, data = data.frame(error))
summary(err.lm)
plot(err.lm)


load(file = paste(output.folder, "data/imputed_data/", an.ID, 
                  "-v", test.params$vars[i], 
                  "-s", test.params$spp[i],
                  ".Rdata", sep = ""))



