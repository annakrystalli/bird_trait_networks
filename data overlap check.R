rm(list=ls())

source("~/Documents/worflows/bird_trait_networks/Setup.R")


#____________________________________________________________________________
### Check D1 against d0 data ###

# Load files

D0 <- read.csv(file = "D0.csv",fileEncoding = "mac")                    
D1 <- read.csv(file = "D1.csv", fileEncoding = "mac") 


# check species name overlap
sum(D0$species %in% D1$species)/length(D0$species)

D1 <- D1[D1$species %in% D0$species,]

D1[,!apply(D1, 2, FUN = function(x){all(is.na(x))})]

compare.ids <- which(vnames$D1 %in% names(D1))
bins <- 50

for(i in 1:length(compare.ids)){
  
  D0var <- vnames$code[compare.ids][i]
  D1var <- vnames$D1[compare.ids][i]
  
  type <- metadata$type[metadata$master.vname == D0var]
  
  D0dat <- numerise(na.omit(D0[D0$var == vnames$code[compare.ids][i], "value"]))
  D1dat <- numerise(na.omit(D1[,vnames$D1[compare.ids][i]]))
  
  if(type %in% c("Con", "Dis")){
    
    range <- range(c(D0dat, D1dat))
    breaks <- seq(range[1], range[2], length.out = bins)
    
    hist(D0dat, breaks = breaks, freq = F, col = makeTransparent("red", 20), 
         border = makeTransparent("red", 40), main = "",
         xlab = paste(metadata$descr[metadata$master.vname == D0var], 
                      metadata$unit[metadata$master.vname == D0var]))
    hist(D1dat, breaks = breaks, freq = F, add = T, col = makeTransparent("blue", 20), 
         border = makeTransparent("blue", 40))
  }
}
  
  
  
}
  
 library(knitr) 
knit2html(paste(output.folder, "Reports/D0_vs_D1_comparison.Rmd", sep = ""))
  


for(i in grep("maturity", names(D1))){
  x <- na.omit(D1[,i])
  if(i == 2){
    hist(x, freq = F, col = makeTransparent(i, 20), 
         border = makeTransparent(i, 40), main = "",
         xlab = names(D1)[i])}else{
           hist(x, freq = F, col = makeTransparent(i, 20), 
                border = makeTransparent(i, 40),
                add = T)
         }
    
}
legend("topright",fill = makeTransparent(grep("maturity", names(D1)), 40), 
       legend = grep("maturity", names(D1), value = T))

