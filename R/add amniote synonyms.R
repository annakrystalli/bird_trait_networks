
# Load amniote csv and master synonyms file
dat <- read.csv(file = "/Users/Anna/Google Drive/bird trait networks/inputs/data/r data/match data/preAves_Synonyms.csv")
synonyms <- read.csv(file = "/Users/Anna/Google Drive/bird trait networks/inputs/data/r data/match data/pre_synonyms.csv")

# clean amniote synonyms
dat <- unique(dat[!dat[,1] == dat[,2],])

# add amniote synonyms to master
synonyms <- rbind(synonyms, data.frame(dat, dataset = "amniote1"), 
                  data.frame(dat[,2:1], dataset = "amniote2"))

# remove duplicates
synonyms <- synonyms[!duplicated(synonyms[,1:2]),] 

# save
write.csv(synonyms, file = "/Users/Anna/Google Drive/bird trait networks/inputs/data/r data/synonyms.csv")
