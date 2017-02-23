# ---- load-sys.ref ----
if(exists("wkf")){
  if(wkf != "rmacro"){
    load(paste(input.folder, "metadata/fileEncodings.rda", sep = ""))
    
    metadata <- read.csv(paste(input.folder, "metadata/","metadata.csv", sep = ""), 
                         stringsAsFactors = F, fileEncoding = fileEncodings["metadata", "encoding"], 
                         na.strings=na.strings, strip.white = T, 
                         blank.lines.skip = T)
    
    data_log <- read.csv(file = paste(input.folder, "metadata/data_log.csv", sep = ""),
                         stringsAsFactors = F, fileEncoding = fileEncodings["data_log", "encoding"], 
                         na.strings=na.strings, strip.white = T, 
                         blank.lines.skip = T)
    
    vnames <- read.csv(paste(input.folder, "metadata/vnames.csv", sep = ""), 
                       stringsAsFactors = F, fileEncoding = fileEncodings["vnames", "encoding"],
                       na.strings=na.strings, strip.white = T, 
                       blank.lines.skip = T, header = T)
  }
}else{
  load(paste(input.folder, "metadata/fileEncodings.rda", sep = ""))
  
  metadata <- read.csv(paste(input.folder, "metadata/","metadata.csv", sep = ""), 
                       stringsAsFactors = F, fileEncoding = fileEncodings["metadata", "encoding"], 
                       na.strings=na.strings, strip.white = T, 
                       blank.lines.skip = T)
  
  data_log <- read.csv(file = paste(input.folder, "metadata/data_log.csv", sep = ""),
                       stringsAsFactors = F, fileEncoding = fileEncodings["data_log", "encoding"], 
                       na.strings=na.strings, strip.white = T, 
                       blank.lines.skip = T)
  
  vnames <- read.csv(paste(input.folder, "metadata/vnames.csv", sep = ""), 
                     stringsAsFactors = F, fileEncoding = fileEncodings["vnames", "encoding"],
                     na.strings=na.strings, strip.white = T, 
                     blank.lines.skip = T, header = T)
}
