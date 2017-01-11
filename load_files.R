# ---- load-sys.ref ----
metadata <- read.csv(paste(input.folder, "metadata/","metadata.csv", sep = ""), 
                     stringsAsFactors = F, fileEncoding = fileEncoding, 
                     na.strings=na.strings, strip.white = T, 
                     blank.lines.skip = T)

data_log <- read.csv(file = paste(input.folder, "metadata/data_log.csv", sep = ""),
                     stringsAsFactors = F, fileEncoding = fileEncoding, 
                     na.strings=na.strings, strip.white = T, 
                     blank.lines.skip = T)

vnames <- read.csv(paste(input.folder, "metadata/","vnames.csv", sep = ""), 
                   stringsAsFactors = F, fileEncoding = fileEncoding,
                   na.strings=na.strings, strip.white = T, 
                   blank.lines.skip = T, header = T)

