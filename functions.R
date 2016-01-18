
# Creates a folder in the data input folder for each qcref variable. These folders can then be 
# populated with qcref data to integrated into the long data format. 
setupInputFolder <- function(input.folder, qcmnames){

  lapply(qcmnames, FUN = function(x){dir.create(paste(input.folder, x, sep =""), 
                                               showWarnings = F)})
}

# Separates qcref variables from data. qcref variables defined in qcmnames. Columns are separated 
# if they are appended with the appropriate qcref variable label (eg ..._ref). 
separateQcRef <- function(dat, qcmnames = c("qc", "observer", "ref", "n", "notes")){
  
  qcref <- vector("list", length(qcmnames))
  names(qcref) <- qcmnames
  
  for(qc.cat in names(qcref)){
    dir.create(paste(qc.cat, "/", sep = ""), showWarnings = F) 
    
    if(any(names(dat) == qc.cat)){
      qcref[[qc.cat]] <- data.frame(species = dat$species, all = dat[, names(dat) == qc.cat])
      dat <- dat[, !names(dat) == qc.cat]
    }else{if(length(names(dat)[grep(paste("_", qc.cat, sep = ""), names(dat))]) > 0){
      qc.var <- names(dat)[grep(paste("_", qc.cat, sep = ""), names(dat))]
      qcref[[qc.cat]] <- dat[, c("species", qc.var)]
      names(qcref[[qc.cat]]) <- gsub(paste("_",qc.cat, sep = ""), "", names(qcref[[qc.cat]]))
      dat <- dat[, !names(dat) %in% qc.var]
    }
    }}
  return(list(data = dat, qcref = qcref))
}

# Extracts appropriate qcref value for data point
getQc <- function(qc.cat, qcref = qcref, spp = spp, var = var){
  if(is.null(qcref[[qc.cat]])){return(NA)}
  if(is.null(dim(qcref[[qc.cat]]))){return(qcref[[qc.cat]])}else{
    if("all" %in% names(qcref[[qc.cat]])){
      return(qcref[[qc.cat]]$all[match(spp, qcref[[qc.cat]]$species)])}else{
        return(qcref[[qc.cat]][cbind(match(spp, qcref[[qc.cat]]$species), 
                                     match(var, names(qcref[[qc.cat]])))])
      }
  }
  
}        

# Matches reference code to full reference. Looks for files with the same name
# as the data file, appended with _code for the file containing the codes and _source containing
# the full references. Source table can also be supplied. Takes filename without extension        
refCode2Source <- function(filename = "bird_ssd7", sources = NULL){
  
  code <- read.csv(paste("ref/", filename, "_code.csv", sep = ""), stringsAsFactors = F) 
  if(is.null(sources)){
    sources <- read.csv(paste("ref/", filename, "_source.csv", sep = ""), stringsAsFactors = F)
  }else{
    sources <- read.csv(paste("ref/", sources, "_source.csv", sep = ""), stringsAsFactors = F)
    
  }
  
  cols <- which(names(code) != "species")
  
  for(i in dim(sources)[1]:1){
    for(j in cols){
      code[,j]  <-  gsub(sources$code[i], sources$ref[i], code[,j])
      
    }}
  
  write.csv(code, paste("ref/", filename, ".csv", sep = ""), row.names = F)
}

# processes BirdFuncTxt file
processBirdFuncTxt <- function(){
  
  BirdFuncDat <- read.delim("raw data/BirdFuncDat.txt")
  code <- read.delim("raw data/BirdFuncDatSources.txt")
  
  
  # clean rows
  BirdFuncDat <- BirdFuncDat[-which(BirdFuncDat$Scientific == ""),]
  
  # rename variable columns
  names(BirdFuncDat)[names(BirdFuncDat) == 'Scientific'] <- "species"
  names(BirdFuncDat)[names(BirdFuncDat) == 'BodyMass.Value'] <- "unsexed.mass"
  names(BirdFuncDat)[names(BirdFuncDat) == 'BodyMass.Source'] <- "unsexed.mass_ref"
  names(BirdFuncDat)[names(BirdFuncDat) == "ForStrat.SpecLevel"] <- "ForStrat_qc"
  names(BirdFuncDat)[names(BirdFuncDat) == "BodyMass.SpecLevel"] <- "unsexed.mass_qc"
  
  names(BirdFuncDat) <- gsub(".Source", "_ref", names(BirdFuncDat))
  names(BirdFuncDat) <- gsub(".Certainty", "_qc", names(BirdFuncDat))
  names(BirdFuncDat) <- gsub(".EnteredBy", "_observer", names(BirdFuncDat))
  
  for(i in 1:length(code$Ref_ID)){
    for(var in grep("_ref", names(BirdFuncDat), value = T)){
      BirdFuncDat[,var] <- gsub(code$Ref_ID[i], code$Full.Reference[i], BirdFuncDat[,var])
    }
  }
  
  BirdFuncDat <- BirdFuncDat[,-c(1:7,9,39:40)]
  
  write.csv(BirdFuncDat, "standardised csv data/BirdFuncDat.csv", row.names = F)
}

# Checks whether qcref variable are appropriately allocated to data columns. Allows the
# assignment of qcref columns to more than one data variable columns. Produces
# appropriately named data.frame with qcref variables assigned to appropriate data variables. 
matchQcRef <- function(dat, file, qcref, var.omit, taxo.var, observer, qc, ref, n, notes, 
                       input.folder){
  

  
  for(qc.cat in names(qcref)){
    qc.cats <- qcref[[qc.cat]]
    
    # if qc.cat in qcref NULL, check whether there is a corresponding qc.cat file and assign
    # to qc.cats
    if(is.null(qc.cats)){
      if(!file %in% list.files(paste(input.folder, qc.cat, "/", sep =""))){}else{
        qc.cats <- read.csv(paste(input.folder, qc.cat, "/", file, sep = ""),
                            stringsAsFactors = F)
        
        # clean qc data
        qc.cats[qc.cats == ""] <- NA
        while(sum(na.omit(qc.cats == " ")) > 0){qc.cats[qc.cats == " "] <- NA}
      }}
    
    if(!is.null(qc.cats)){
      if("all" %in% names(qc.cats)){qcref[[qc.cat]] <- qc.cats[,c("species", "all")]}else{
        
        # vector of data vars to check for qc.caterences
        check.vars <- names(dat)[!(names(dat) %in% c(var.omit, taxo.var))]
        
        if(all(check.vars %in%  names(qc.cats))){qcref[[qc.cat]] <- qc.cats[c("species", check.vars)]}else{
          
          if(length(grep(paste("_", qc.cat, "_group", sep = ""), 
                         grep(gsub(".csv", "",file), list.files(paste(qc.cat, "/", sep ="")), 
                              value = T))) == 0){
            
            # create csv in which qcref group names can be assigned to individual data variables
            qc.cat.grp <- data.frame(var = check.vars, qc.cat.grp = "")
            qc.cat.grp$qc.cat.grp[qc.cat.grp$var %in% names(qc.cats)] <- qc.cat.grp$var[qc.cat.grp$var %in% names(qc.cats)]
            write.csv(qc.cat.grp, paste(paste(input.folder, qc.cat, "/", sep =""), 
                                        gsub(".csv", "",file), "_", qc.cat, "_group.csv", sep = ""),
                      row.names = F)
            
            stop(paste(qc.cat,".group file created, ","update _",qc.cat,".group file to proceed", sep = ""))
          }else{
            # load csv in which qcref group names are assigned to individual data variables
            qc.cat.grp  <- read.csv(paste(paste(input.folder, qc.cat, "/", sep =""), gsub(".csv", "",file), "_", 
                                          qc.cat, "_group.csv", sep = ""), 
                                    stringsAsFactors = F)
            
            # Check that all variables are assigned to valid qcref data columns
            if(!all(na.omit(qc.cat.grp$qc.cat.grp) %in% names(qc.cats))){
              stop(paste("assigned qcref group names does not match supplied qcref data names, update _",
                         qc.cat,".group file to proceed", sep = ""))}else{
                           # Make sure ALL variables have reference data         
                           if(qc.cat == "ref" & any(is.na(qc.cat.grp$qc.cat.grp))){
                             stop(paste("variables missing reference column. update ", 
                                        qc.cat,".group file to proceed", sep = ""))
                           }
                           
                           # Isolate variables to be assigned qcref data. Create new dataframe containing the 
                           # appropriate qcref column for each variable. 
                           # Name with data variables and update appropriate qcref slot           
                           check.vars <- qc.cat.grp$var[which(!is.na(qc.cat.grp$qc.cat.grp))]
                           qc.cat.table <- data.frame(species = dat$species, matrix(NA, nrow = dim(dat)[1], 
                                                                                    ncol = length(check.vars)))
                           names(qc.cat.table) <- c("species", check.vars)
                           qc.cat.table[,check.vars] <- qc.cats[match(qc.cat.table$species, qc.cats$species), 
                                                                qc.cat.grp$qc.cat.grp[match(check.vars, qc.cat.grp$var)]]
                           
                           print("vars matched successfully to _qc.cat_group")
                           
                           qcref[[qc.cat]] <- qc.cat.table
                         }}}
        
        
        
      }}
    
    # if value for qcref variable has been supplied through function, set as value for all
    # data points and variables
    if(is.null(qcref[[qc.cat]])){
      if(!is.null(get(qc.cat))){
        qcref[[qc.cat]] <- data.frame(species = dat$species, all = get(qc.cat))}else{
          if(qc.cat == "ref"){stop("Processing stopped: no reference information")}else{
            print(paste("Warning: NULL data for qc.cat:", qc.cat))
          }
        }
      
    }}
  return(qcref)
}


#process data set ready for matching. Calls separateQcRef to separate  qcref variables from data,
# matchQcRef to match qcref variables to data variables, prepares data for compiling and checks
# that metadata information has been supplied for all data variables. 
processDat <- function(file = "ASR_mortality_to_Anna_Gavin.csv",label = F,
                       taxo.dat, var.omit, input.folder = input.folder,
                       observer = NULL, qc = NULL, ref = NULL, n = NULL, notes = NULL,
                       master.vname = "master.vname"){
  
  
  dat <- read.csv(paste("csv/", file, sep = ""),  stringsAsFactors=FALSE)
  
  
  if(anyDuplicated(dat$species) > 0){stop("duplicate species name in match dat")}
  
  if(any(dat$species == "")){
    dat<- dat[-which(dat$species == ""),]}
  
  #Make sure there are no empty cells and replace any with NA cells
  dat[which(dat== "", arr.ind = T)] <- NA
  require(stringr)
  
  if(label){
    names(dat)[names(dat) != "species"] <- paste(names(dat), dat.ID, sep = "_")[
      names(dat) != "species"]
  }
  
  dat$species <- gsub(" ", "_", dat$species)
  
  # separate qcref data
  dat.l <- separateQcRef(dat)
  
  # create qcref object
  dat.l$qcref <- matchQcRef(dat = dat.l$data, file, qcref = dat.l$qcref, var.omit, taxo.var,
                            observer, qc, ref, n, notes, input.folder)
  
  if("parent.spp" %in% names(dat.l$data)){
    parent.spp <- dat.l$data$parent.spp
    dat.l$data <- dat.l$data[,names(dat.l$data) != "parent.spp"]
  }else{parent.spp <- NULL}
    
  # prepare dataset
  dat.l$data <- dataMatchPrep(dat.l$data)
  

  if(!is.null(parent.spp)){
    dat.l$data$parent.spp <- parent.spp
    dat.l$data$subspp[!is.na(dat.l$data$parent.spp)] <- T
  }
  
  # check metadata complete
  metadata <- read.csv(paste(input.folder, "metadata/metadata.csv", sep = ""),  stringsAsFactors=FALSE)
  
  metavar <- c("synonyms", "data.status", "species", "subspp", "parent.spp")
  
  if(!all(names(dat.l$data) %in% c(metadata[,master.vname], metavar))){
    print(names(dat.l$data)[!(names(dat.l$data) %in% c(metadata[,master.vname], metavar))])
    stop("metadata missing, metadata file needs updating")}
  
  return(dat.l)
}

# extracts taxonomic information for species. Matches to original taxonomy used on project so added  
# species are matched using parent.spp or syns information
spp2taxoMatch <- function(spp, parent.spp, taxo.table){
  
  if(is.null(taxo.table)){
  spp2taxo <- read.csv("r data/spp_to_taxo.csv", stringsAsFactors = F)}else{
  spp2taxo <- taxo.table
  }
  
  spp.id <- spp %in% spp2taxo$species
  pspp.id <- parent.spp %in% spp2taxo$species

  taxo.id <- spp.id == T | pspp.id == T
  
  if(!all(taxo.id)){
    stop(c("no spp2taxo data for species", unique(spp[!taxo.id])))}else{
      if(all(spp.id)){
        dat <- spp2taxo[match(spp, spp2taxo$species),]
      }else{
        sppp <- spp
        p <- parent.spp %in% spp2taxo$species & !spp %in% spp2taxo$species
        sppp[p] <- parent.spp[p]
        dat <- spp2taxo[match(sppp, spp2taxo$species),]
      
        if(is.null(taxo.table)){  
        # Update spp2taxo file
        add.dat <- cbind(species = spp[p],spp2taxo[match(parent.spp[p], spp2taxo$species),-1])
        add.dat$subspp <- TRUE
        add.dat$parent.spp <- parent.spp[p]
        write.csv(rbind(spp2taxo, add.dat), "r data/spp_to_taxo.csv", row.names = F)}
      }
    }
  
  dat <- data.frame(species = spp, dat[,c("order", "family")])  
  
  return(dat)
}

# Compiles dataset into format compatible with appending to master database. Takes 
# match object m.
matchMSToMaster <-  function(m, master, taxo.var = taxo.var, var.omit = var.omit, 
                             input.folder, output.folder, ignore.unmatched = F, synonyms,
                             taxo.table){
  
  data <- m$data
  sub <- m$sub
  set <- m$set
  qcref <- m$qcref
  spp.list <- m$spp.list
  
  # Check whether matching required and match
  unmatched <- get(sub)$species[!(get(sub)$species %in% get(set)$species)]
  if(length(unmatched) != 0){
    
    m <- dataSppMatch(m, unmatched = unmatched, ignore.unmatched = ignore.unmatched, synonyms = synonyms)
    data <- m$data
    
  }
  
  # Load taxo variable lookup table
  if(is.null(taxo.table)){
  spp2taxo <- read.csv(paste(input.folder, "r data/spp_to_taxo.csv", sep = ""), stringsAsFactors = F)}else{
    spp2taxo <- taxo.table
  }
  
  #make vector of data variables to be added
  match.vars <- names(data)[!names(data) %in% c(taxo.var, var.omit, c("synonyms", "data.status"))]
  match.dat <- data.frame(data[, match.vars])
  names(match.dat) <- match.vars
  
  #find non NA values in match data. Match arr.indices to spp and variable names (for QA)
  id <- which(!is.na(match.dat), arr.ind = T)
  spp <- as.character(data[,"species"][id[, "row"]])
  parent.spp <- as.character(data[,"parent.spp"][id[, "row"]])
  syns <- as.character(data[,"synonyms"][id[, "row"]])
  var <- as.character(match.vars[id[, "col"]])
  
  mdat <- try(cbind(spp2taxoMatch(spp, parent.spp, taxo.table), 
                subspp = data[id[,"row"], "subspp"],
                parent.spp = parent.spp,
                var = var, 
                value = match.dat[id], data = m$data.ID,
                synonyms = syns, 
                data.status = data[id[,"row"], "data.status"],
                qc = getQc("qc", qcref = qcref, spp = spp, var = var),
                observer = getQc("observer", qcref = qcref, spp = spp, var = var),
                ref = getQc("ref", qcref = qcref, spp = spp, var = var),
                n = getQc("n", qcref = qcref, spp = spp, var = var)))
  
  if(class(mdat) == "try-error"){return(m)}else{
  
    dir.create(paste(output.folder, "data/", sep = ""), showWarnings = F)
    dir.create(paste(output.folder, "data/match objects/", sep = ""), showWarnings = F)
    
  save(m, file = paste(output.folder, "data/match objects/", m$data.ID, " match object.RData", 
                       sep = ""))
  
  return(list(mdat = mdat, spp.list = m$spp.list))}
}

# Processes ITIS synonyms data into a species synonym dataset 
ITISlookUpData <- function(version=NULL){
  aves.names <- read.csv("r data/match data/Aves synonym data (ITIS).csv", stringsAsFactors=FALSE)
  aves.codes <- read.csv("r data/match data/spp code matches.csv", stringsAsFactors=FALSE)
  
  species <- aves.names$species[match(aves.codes$Main, aves.names$code)]
  synonyms <- aves.names$species[match(aves.codes$Synonym, aves.names$code)]
  
  itis.match <- data.frame(species, synonyms, stringsAsFactors = F)
  itis.match <- itis.match[complete.cases(itis.match),]
  
  if(version == 2){names(itis.match) <- c("synonyms", "species")}
  
  return(itis.match)}

# Outputs details of master species name match lookup
matchMetrics <- function(data, master, match.lengths){
  n <- min(nrow(master), nrow(data))
  print(paste("total data set:", length(data$species)))
  print(paste("total matched:", n - match.lengths["unmatched"]))            
  print(paste("total unmatched:", match.lengths["unmatched"]))                       
  print(paste("total ITIS matches:", 
              match.lengths["unmatched"] - match.lengths["unmatched1"]))
  print(paste("total Birdlife:", 
              match.lengths["unmatched1"] - match.lengths["unmatched2"]))
  print(paste("total master.match:", 
              match.lengths["unmatched2"] - match.lengths["unmatched3"]))
  print(paste("manual matches:", match.lengths["unmatched3"]))}




#Look up unmantched vector of species. If vector contains unmatched data species, 
#Prepares data for matching of species to master species name. Creates synonyms column to link back to original data and 
# and data status, used to indicate data has bee prepared but also whether state of data row is original or has been added.

dataMatchPrep <- function(data){
  if("data.status" %in% names(data)){
    print("Data already prep-ed")
    return(data)
  }else{
    dt <- data.frame(data[,names(data) != "species"], stringsAsFactors = FALSE)
    names(dt) <- names(data)[names(data) != "species"]
    data <- data.frame(species = data$species, synonyms = data$species, 
                       subspp = FALSE, parent.spp = NA,
                       data.status = "original", dt,
                       stringsAsFactors = FALSE)
    return(data)}
}


# Create match object
matchObj <- function(data.ID, spp.list, data = dl[[data.ID]], status = "unmatched", 
                     sub = data.match.params$sub[data.match.params$data.ID == data.ID],
                     qcref = qcref){
  
  
  if(sub == "spp.list"){
    set <- "data"
  }
  if(sub == "data"){
    set <- "spp.list"
  }   
  
  m <- list(data.ID = data.ID, spp.list = spp.list, data = data, sub = sub, 
            set = set, status = status, qcref = qcref)
  return(m)
}




# look up unmatched species in table (lookup.dat) of known match pairs. 
# If list is of unmatched data species (ie dataset species is a subset of spp.list), 
# match to known synonyms. If list is of unmatched spp.list species (ie spp.list a subset of 
# data$species), match to known species. 
sppMatch <- function(X, unmatched = unmatched, lookup.dat){
  
  data <- X$data
  spp.list <- X$spp.list
  sub <- X$sub
  set <- X$set
  
  if(sub == "spp.list"){
    lookup <- lookup.dat$species
    lookupin <- lookup.dat$synonyms
  }
  if(sub == "data"){
    lookup <- lookup.dat$synonyms
    lookupin <- lookup.dat$species
  }
  
  
  #Find synonyms of unmatched species 
  lookupin.match <- as.character(lookupin[lookup %in% unmatched]) #find syns of unmatched species with matches in lookup
  lookup.match <- as.character(lookup[lookup %in% unmatched]) # find unmatched spp with matches in lookup
  
  # Compile match pairs for species whose synonyms match species names in lookup.data
  spp <- lookup.match[lookupin.match %in% get(set)$species] #find unmacthed species with synonym matches in set
  syns <- lookupin.match[lookupin.match %in% get(set)$species] #trim synonym matches
  
  #Compile match pairs into dataframe 
  if(sub == "data"){match.data <- data.frame(synonyms = spp, species = syns, stringsAsFactors = F)}
  if(sub == "spp.list"){match.data <- data.frame(synonyms = syns, species = spp, stringsAsFactors = F)}
  
  # Check there are no duplicate species matches in match.data table. If so add duplicates to spp list
  # and keep datapoint under synonym species name. Mark as subspp 
  if(any(duplicated(match.data$species))){
    master.add <- match.data$synonyms[duplicated(match.data$species)]
    
    data$subspp[data$synonyms %in% master.add] <- TRUE
    data$parent.spp[match(master.add, data$synonyms)] <- match.data$species[
      match(master.add, match.data$synonyms)]
    
    spp.list <- rbind(spp.list, data.frame(species = master.add))
    
    #remove any species from match.data already added to master
    match.data <- match.data[-which(match.data$synonyms %in% master.add),]
  }
  
  #Check that data spp name change won't cause a single master species name to be associated with two
  # different data points. This would cause the original data point (which is a straught match to the dataset)
  # to be overwritten with most likely data for a subspecies. Identify such cases and add them to the master.
  # The data point can then form a straight match and is identified by the column subspp.
  
  if(any(match.data$species %in% data$species)){
    
    master.add <- match.data$synonyms[match.data$species %in% data$species]
    
    data$subspp[data$synonyms %in% master.add] <- TRUE
    data$parent.spp[match(master.add, data$synonyms)] <- match.data$species[
      match(master.add, match.data$synonyms)]
    
    spp.list <- rbind(spp.list, data.frame(species = master.add))
    
    #remove any species from match.data already added to master
    match.data <- match.data[-which(match.data$synonyms %in% master.add),]
  } 
  
  # Also check that any of the data species names about to be replaced do not already map to another species in the master.
  # Identify and duplicate any such datarows so that the orginal link to the master will remain. Species can be linked back
  # to data via the synonoys column in data. Additionally check against previously added species to avoid duplicate additions.
  if(any(match.data$synonyms %in% spp.list$species)){
    
    #find synonyms in master species
    data.add <- match.data[match.data$synonyms %in% spp.list$species,]
    
    if(nrow(data.add)>=1){
      #duplicate data row for species to be added to data
      add.dd <- data.frame(species = data.add$species, data[match(data.add$synonyms, 
                                                                  data$species),names(data) != "species"],
                           stringsAsFactors = F)
      add.dd$data.status <- "duplicate"
      
      #add to data
      data <- rbind(data, add.dd)
    }
    
    match.data <- match.data[-which(match.data$species %in% add.dd$species),]
  }   
  
  #update species names in data$species
  data$species[match(match.data$synonyms, data$species)] <- match.data$species
  
  
  m <- matchObj(data.ID=X$data.ID, spp.list = spp.list, data, status = X$status,
                sub = X$sub, qcref = X$qcref)
  
  return(m)}


# match data set to master species list using all available known match pair tables.
dataSppMatch <- function(m, unmatched, ignore.unmatched = ignore.unmatched, synonyms){
  
  
  sub <- m$sub
  set <- m$set

  # unmatched
  if(sub == "data"){
    rm <- synonyms$synonyms[synonyms$synonyms %in% unmatched & synonyms$species %in% c("Extinct", "New")]
    m$data <- m$data[!(m$data$species %in% rm),]
  }
  
    m <- sppMatch(m, unmatched = unmatched, lookup.dat = synonyms)
    
    #generate next unmatched species list
    unmatched <- m[[sub]]$species[!(m[[sub]]$species %in% m[[set]]$species)]
    
    
    # if no more species unmatched break loop
    if(length(unmatched) == 0){
      print(paste(data.ID, "match complete"))
     }else{
    
        print(paste("match incomplete,",length(unmatched), sub, "datapoints unmatched"))
        
        if(ignore.unmatched){if(sub == "data"){m$data <- m$data[!m$data$species %in% unmatched,]}}else{
      
        
    # if all match pair datasets checked and species remain unmatched write manual match spp list
      dir.create(paste(input.folder, "r data/", sep = ""), showWarnings = F)
      dir.create(paste(input.folder, "r data/match data/", sep = ""), showWarnings = F)
      
      write.csv(data.frame(synonyms = if(sub == "spp.list"){""}else{unmatched},
                           species = if(sub == "spp.list"){unmatched}else{""}),
                paste("r data/match data/",data.ID," mmatch.csv", sep = ""),
                row.names = F)
      print(paste("unmatched species list saved in file 'Data.IDmmatch.csv'"))
      stop("manually match and save as 'Data.IDmmatched.csv' to continue")}}
    
    if(sub == "spp.list"){m$data <- m$data[m$data$species %in% spp.list$species,]}  
    return(m)}


addVars <- function(data, master){
  
  vars <- names(data)[!(names(data) %in% c("species", "synonyms", "data.status"))]
  
  for(var in vars){
    if(is.character(data[,var])){data[,var][data[,var]==""] <- NA}
    add <- data.frame(spp = match(data$species, master$species), 
                      dat = data[,var], stringsAsFactors = F)
    add <- add[complete.cases(add),]
    #if(!any(is.na(as.numeric(add$dat)))){add$dat <- as.numeric(add$dat)}
    
    var.col <- data.frame(rep(NA, dim(master)[1]))
    var.col[add[,"spp"],] <- add[,"dat"]
    
    names(var.col)<- var
    
    master <- data.frame(master, var.col, stringsAsFactors = F)}
  
  return(master)}

# Matches sppecies names to identifiers in data. ids needs to be a 2 column dataframe. 
# str = c(index, species), file = path to data to be matched with column spp_no instead of species.
# writes processed data file to csv, appending "_SPP" to file name.

IDSppMatch <- function(file = "Display & resource_scores.csv", 
                       ids = read.csv("r data/id_to_spp.csv", stringsAsFactors = F)){
  dd <- read.csv(paste("standardised csv data/", file, sep = ""), stringsAsFactors = F)
  dd <- dd[!apply(dd[,names(dd) != "spp_no"], 1, FUN = function(x){all(is.na(x))}),]
  
  ids <- ids
  dd <- data.frame(species = ids$species[match(dd$spp_no, ids$index)], dd[,names(dd) != "spp_no"])
  
  write.csv(dd, paste("standardised csv data/", gsub(".csv", "",file), "_SPP.csv",sep = ""),
            row.names = F)
  
}


# Tests whether a proposed synonym/species has a match in the spp.list/data and updates the mmatched file for the data set.
# Takes a match object (x) all information needed is stored within the file
testSynonym <- function(spp, x){
  #identify next species being matched and print
  mmatch <- read.csv(paste("r data/match data/",x$data.ID," mmatched.csv", sep = ""),
                     stringsAsFactors=FALSE)  
  
  sub <- x$sub 
  
  
  if(sub == "spp.list"){
    set <- "data"
    lookup <- "species"
    lookupin <- "synonyms"
  }
  if(sub == "data"){
    set <- "spp.list"
    lookup <- "synonyms"
    lookupin <- "species"
  }   
  
  spp.m <- mmatch[[lookup]][min(which(mmatch[lookupin] == "" | is.na(mmatch[[lookupin]])))]
  print(paste("match", lookup, spp.m))
  
  #test potential synonym 
  if(spp %in% c("Extinct","New")){
    print(paste("match", lookupin, spp))
    mmatch[mmatch[lookup] == spp.m, lookupin] <- spp
    next.spp <-mmatch[[lookup]][min(which(mmatch[lookupin] == "" | is.na(mmatch[[lookupin]])))]
    
    print(paste("next", lookup, ":", next.spp))
    write.csv(mmatch, paste("r data/match data/",x$data.ID," mmatched.csv", sep = ""),
              row.names = F)
  }else{
    spp <- gsub(" ","_", spp)
    match <- any(spp %in% x[[set]]$species)  
    print(paste("match", lookupin, spp))
    print(match)
    
    if(match){
      
      mmatch[mmatch[lookup] == spp.m, lookupin] <- spp
      next.spp <-mmatch[[lookup]][min(which(mmatch[lookupin] == "" | is.na(mmatch[[lookupin]])))]
      
      print(paste("next", lookup, ":", next.spp))
      write.csv(mmatch, paste("r data/match data/",x$data.ID," mmatched.csv", sep = ""),
                row.names = F)
      
    }}
  
}
whichNext <- function(x = output){
  #identify next species being matched and print
  mmatch <- read.csv(paste("r data/match data/",x$data.ID," mmatched.csv", sep = ""),
                     stringsAsFactors=FALSE)  
  
  sub <- x$sub 
  
  
  if(sub == "spp.list"){
    set <- "data"
    lookup <- "species"
    lookupin <- "synonyms"
  }
  if(sub == "data"){
    set <- "spp.list"
    lookup <- "synonyms"
    lookupin <- "species"
  }   
  
  spp.m <- mmatch[[lookup]][min(which(mmatch[lookupin] == "" | is.na(mmatch[[lookupin]])))]
  print(paste("match", lookup, spp.m))
  
}


