
#' Setup meta folders in inputs folder.
#' 
#' Creates a folder in the data input folder for each meta variable. These folders can then be 
#' populated with meta data to integrated into the long data format. 
#' @param input.folder file path to input data folder.
#' @param meta.vars vector containing names of meta.vars
#' @keywords meta 
#' @export
#' setupInputFolder()

setupInputFolder <- function(input.folder, meta.vars){
  
  lapply(meta.vars, FUN = function(x){dir.create(paste(input.folder, x, sep =""), 
                                                 showWarnings = F)})
}

#' Create metadata list
#' 
#' Creates a named metadata list of the appropriate length to match supplied vector of meta.vars
#' @param meta.vars vector containing names of meta.vars
#' @keywords meta 
#' @export
#' createMeta()

createMeta <- function(meta.vars){
  
  meta <- vector("list", length(meta.vars))
  names(meta) <- meta.vars
  
  return(meta)
}

#############################################

#' Separate meta.vars
#' 
#' Separates meta variables from data into a `data` and a `meta` dataframe. Columns are separated if their name matches a meta.var 
#' or if they are appended with the appropriate meta variable label (eg ..._ref). Return a list
#' containing the data and separated metadata data.frames.
#' @param dat data.frame containg data in the form species (row) by variable (columns)
#' @param meta named metadata list to collect metadata
#' @keywords meta 
#' @export
#' separateDatMeta()

separateDatMeta <- function(dat, meta){
  
  for(meta.var in names(meta)){
    
    if(any(names(dat) == meta.var)){
      # if meta.var data for all variables is single column named `meta.var`
      meta[[meta.var]] <- data.frame(species = dat$species, all = dat[, names(dat) == meta.var])
      dat <- dat[, !names(dat) == meta.var]
    }else{# if meta.var data for individual variables is in columns named `_meta.var`
      if(length(names(dat)[grep(paste("_", meta.var, sep = ""), names(dat))]) > 0){
        vmeta.var <- names(dat)[grep(paste("_", meta.var, sep = ""), names(dat))]
        meta[[meta.var]] <- dat[, c("species", vmeta.var)]
        names(meta[[meta.var]]) <- gsub(paste("_",meta.var, sep = ""), "", names(meta[[meta.var]]))
        dat <- dat[, !names(dat) %in% vmeta.var]
      }
    }}
  return(list(data = dat, meta = meta))
}


#############################################

#' Get metadata values associated with observations
#' 
#' Extracts appropriate metadata values for specified species vs variable data points. 
#' Data points specified by vectors of species and variable names.   
#' @param input.folder file path to input data folder.
#' @param meta.var name of meta.var
#' @param meta metadata list
#' @param spp vector of species names
#' @param var vector of variable names
#' @keywords meta 
#' @export
#' getMeta()


getMeta <- function(meta.var, meta = meta, spp = spp, var = var){
  
  if(is.null(meta[[meta.var]])){return(NA)}
  if(is.null(dim(meta[[meta.var]]))){return(meta[[meta.var]])}else{
    if("all" %in% names(meta[[meta.var]])){
      return(meta[[meta.var]]$all[match(spp, meta[[meta.var]]$species)])}else{
        return(meta[[meta.var]][cbind(match(spp, meta[[meta.var]]$species), 
                                      match(var, names(meta[[meta.var]])))])
      }
  }
  
}        

#############################################

#' Substitute reference code with full reference
#' 
#' Sustitutes reference codes with full references across cells of a data.frame
#' @param ref.codes data.frame containing coded reference data. 
#' @param ref.table data.frame containing code to full reference look up table. Columns must 
#' be named `code` and `ref`.
#' @keywords meta 
#' @export
#' code2FullRef()

code2FullRef <- function(ref.codes, ref.table){
  
  cols <- which(names(ref.codes) != "species")
  
  for(i in dim(ref.table)[1]:1){
    for(j in cols){
      ref.codes[,j]  <-  gsub(ref.table$code[i], ref.table$ref[i], ref.codes[,j])
      
    }}

  return(ref.codes)
}


# Checks whether meta variable are appropriately allocated to data columns. Allows the
# assignment of meta columns to more than one data variable columns. Produces
# appropriately named data.frame with meta variables assigned to appropriate data variables. 
matchmeta <- function(dat, meta, meta.var, file, input.folder, metav.dd = NULL, 
                      ignore.vars = c(var.omit, taxo.var), write = T){
  
    metav.dd <- meta[[meta.var]]
    
    # if meta.var in meta NULL, check whether there is a corresponding meta.var file and assign
    # to metav.dd
    if(is.null(metav.dd)){
      if(is.null(file)){}else{
        if(!file %in% list.files(paste(input.folder, meta.var, "/", sep =""))){}else{
          metav.dd <- read.csv(file, stringsAsFactors = F)
        }}}
    
    if(!is.null(metav.dd)){
      if("all" %in% names(metav.dd)){meta[[meta.var]] <- metav.dd[,c("species", "all")]}else{
        
        # vector of data vars to check for meta.varerences
        check.vars <- names(dat)[!(names(dat) %in% ignore.vars)]
        
        # if all meta variable names match data var names, return data
        if(all(check.vars %in% names(metav.dd))){meta[[meta.var]] <- metav.dd[c("species", 
                                                                                 check.vars)]}else{
          # if
          if(length(grep(paste("_", meta.var, "_group", sep = ""), 
                         grep(gsub(".csv", "",file), 
                              list.files(paste(input.folder, meta.var, "/", sep ="")), 
                              value = T)
                         )
                    ) == 0){
            
            # create csv in which individual data variables can be assigned to meta group names
            meta.var.grp <- data.frame(var = check.vars, meta.var.grp = "")
            meta.var.grp$meta.var.grp[
              meta.var.grp$var %in% names(metav.dd)] <- meta.var.grp$var[
                meta.var.grp$var %in% names(metav.dd)]
            
            

            write.csv(meta.var.grp, paste(paste(input.folder, meta.var, "/", sep =""), 
                                          gsub(".csv", "",file), "_", meta.var, "_group.csv", 
                                          sep = ""),
                      row.names = F)
            
            stop(paste(meta.var,".group file created, ","update _",meta.var,
                       ".group file to proceed", sep = ""))
          }else{
            # load csv in which meta group names are assigned to individual data variables
            meta.var.grp  <- read.csv(paste(paste(input.folder, meta.var, "/", sep =""), gsub(".csv", "",file), "_", 
                                            meta.var, "_group.csv", sep = ""), 
                                      stringsAsFactors = F)
            
            # Check that all variables are assigned to valid meta data columns
            if(!all(na.omit(meta.var.grp$meta.var.grp) %in% names(metav.dd))){
              stop(paste("assigned meta group names does not match supplied meta data names, update _",
                         meta.var,".group file to proceed", sep = ""))}else{
                           # Make sure ALL variables have reference data         
                           if(meta.var == "ref" & any(is.na(meta.var.grp$meta.var.grp))){
                             stop(paste("variables missing reference column. update ", 
                                        meta.var,".group file to proceed", sep = ""))
                           }
                           
                           # Isolate variables to be assigned meta data. Create new dataframe containing the 
                           # appropriate meta column for each variable. 
                           # Name with data variables and update appropriate meta slot           
                           check.vars <- meta.var.grp$var[which(!is.na(meta.var.grp$meta.var.grp))]
                           meta.var.table <- data.frame(species = dat$species, matrix(NA, nrow = dim(dat)[1], 
                                                                                      ncol = length(check.vars)))
                           names(meta.var.table) <- c("species", check.vars)
                           meta.var.table[,check.vars] <- metav.dd[match(meta.var.table$species, metav.dd$species), 
                                                                    meta.var.grp$meta.var.grp[match(check.vars, meta.var.grp$var)]]
                           
                           print("vars matched successfully to _meta.var_group")
                           
                           meta[[meta.var]] <- meta.var.table
                         }}}
        
        
        
      }}
    
    # if value for meta variable has been supplied through function, set as value for all
    # data points and variables
    if(is.null(meta[[meta.var]])){
      if(!is.null(get(meta.var))){
        meta[[meta.var]] <- data.frame(species = dat$species, all = get(meta.var))}else{
          if(meta.var == "ref"){stop("Processing stopped: no reference information")}else{
            print(paste("Warning: NULL data for meta.var:", meta.var))
          }
        }
      
    }
    
    
  return(meta)
}


#process data set ready for matching. Calls separateDatMeta to separate  meta variables from data,
# matchmeta to match meta variables to data variables, prepares data for compiling and checks
# that metadata information has been supplied for all data variables. 
processDat <- function(file = "ASR_mortality_to_Anna_Gavin.csv", dat, label = F,
                       taxo.dat, var.omit, input.folder = input.folder, meta,
                       master.vname = "master.vname"){
  
  if(is.null(dat)){
    dat <- read.csv(paste(input.folder, "csv/", file, sep = ""),  stringsAsFactors=FALSE)}
  
  
  
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
  
  # separate meta data
  dat.l <- separateDatMeta(dat)
  
  for(meta.var in meta.vars){
    
    metav.dd <- read.csv(paste(input.folder, meta.var, "/", file, sep = ""),
             stringsAsFactors = F)
  
   # create meta object
    dat.l$meta <- matchmeta(dat = dat.l$data, meta.var = meta.var, meta = dat.l$meta, 
                            metav.dd = metav.dd)}
  
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
                             taxo.table, trim.dat = T, retain.dup = T){
  
  data <- m$data
  sub <- m$sub
  set <- m$set
  meta <- m$meta
  spp.list <- m$spp.list
  
  # Check whether matching required and match
  unmatched <- get(sub)$species[!(get(sub)$species %in% get(set)$species)]
  if(length(unmatched) != 0){
    
    m <- dataSppMatch(m, unmatched = unmatched, ignore.unmatched = ignore.unmatched, 
                      synonyms = synonyms, trim.dat = trim.dat, retain.dup = retain.dup)
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
                    qc = getMeta("qc", meta = meta, spp = spp, var = var),
                    observer = getMeta("observer", meta = meta, spp = spp, var = var),
                    ref = getMeta("ref", meta = meta, spp = spp, var = var),
                    n = getMeta("n", meta = meta, spp = spp, var = var)))
  
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
                     meta = meta){
  
  
  if(sub == "spp.list"){
    set <- "data"
  }
  if(sub == "data"){
    set <- "spp.list"
  }   
  
  m <- list(data.ID = data.ID, spp.list = spp.list, data = data, sub = sub, 
            set = set, status = status, meta = meta)
  return(m)
}




# look up unmatched species in table (lookup.dat) of known match pairs. 
# If list is of unmatched data species (ie dataset species is a subset of spp.list), 
# match to known synonyms. If list is of unmatched spp.list species (ie spp.list a subset of 
# data$species), match to known species. 
sppMatch <- function(X, unmatched = unmatched, lookup.dat, retain.dup = T){
  
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
  
  match.data <- unique(match.data)
  
  # Check there are no duplicate species matches in match.data table. If so, add duplicates to spp list
  # and keep datapoint under synonym species name. Mark as subspp ONLY IF retain.dup = T
  if(any(duplicated(match.data$species))){
    master.add <- match.data$synonyms[duplicated(match.data$species)]
    
    if(retain.dup){
      data$subspp[data$synonyms %in% master.add] <- TRUE
      data$parent.spp[match(master.add, data$synonyms)] <- match.data$species[
        match(master.add, match.data$synonyms)]
      
      spp.list <- rbind(spp.list, data.frame(species = master.add))}
    
    #remove any species from match.data already added to master
    match.data <- match.data[-which(match.data$synonyms %in% master.add),]
  }
  
  #Check that data spp name change won't cause a single master species name to be associated with two
  # different data points. This would cause the original data point (which is a straight match to the dataset)
  # to be overwritten with most likely data for a subspecies. Identify such cases and add them to the master.
  # The data point can then form a straight match and is identified by the column subspp.
  
  if(any(match.data$species %in% data$species)){
    
    master.add <- match.data$synonyms[match.data$species %in% data$species]
    
    if(retain.dup){
      data$subspp[data$synonyms %in% master.add] <- TRUE
      data$parent.spp[match(master.add, data$synonyms)] <- match.data$species[
        match(master.add, match.data$synonyms)]
      
      spp.list <- rbind(spp.list, data.frame(species = master.add))}
    
    #remove any species from match.data already added to master
    match.data <- match.data[-which(match.data$synonyms %in% master.add),]
  } 
  
  # Also check that any of the data species names about to be replaced do not already map 
  # to another species in the master. Identify and duplicate any such datarows so that the
  # orginal link to the master will remain. Species can be linked back
  # to data via the synonoys column in data. Additionally check against previously 
  # added species to avoid duplicate additions.
  
  if(any(duplicated(match.data$synonyms))){
    data.add <- match.data[duplicated(match.data$synonyms),]
    match.data <- match.data[!duplicated(match.data$synonyms),]
  }
  
  if(any(match.data$synonyms %in% spp.list$species)){
    
    #find synonyms in master species
    if(exists("data.add")){
      data.add <- rbind(data.add, match.data[match.data$synonyms %in% spp.list$species,])}else{
        data.add <- match.data[match.data$synonyms %in% spp.list$species,]}
  }
  
  
  if(exists("data.add")){
    #duplicate data row for species to be added to data
    add.dd <- data.frame(species = data.add$species, data[match(data.add$synonyms, 
                                                                data$species),names(data) != "species"],
                         stringsAsFactors = F)
    add.dd$data.status <- "duplicate"
    
    #add to data
    data <- rbind(data, add.dd)
    
    match.data <- match.data[-which(match.data$species %in% add.dd$species),]
  }   
  
  #update species names in data$species
  data$species[match(match.data$synonyms, data$species)] <- match.data$species
  
  
  m <- matchObj(data.ID=X$data.ID, spp.list = spp.list, data, status = X$status,
                sub = X$sub, meta = X$meta)
  
  return(m)}


# match data set to master species list using all available known match pair tables.
dataSppMatch <- function(m, unmatched, ignore.unmatched = ignore.unmatched, synonyms, 
                         trim.dat = T, retain.dup = T){
  
  
  sub <- m$sub
  set <- m$set
  
  # unmatched
  if(sub == "data"){
    rm <- synonyms$synonyms[synonyms$synonyms %in% unmatched & synonyms$species %in% c("Extinct", "New")]
    m$data <- m$data[!(m$data$species %in% rm),]
  }
  
  m <- sppMatch(m, unmatched = unmatched, lookup.dat = synonyms, retain.dup = retain.dup)
  
  #generate next unmatched species list
  unmatched <- m[[sub]]$species[!(m[[sub]]$species %in% m[[set]]$species)]
  
  
  # if no more species unmatched break loop
  if(length(unmatched) == 0){
    print(paste(m$data.ID, "match complete"))
    if(sub == "spp.list" & trim.dat){m$data <- m$data[m$data$species %in% m$spp.list$species,]}
  }else{
    
    print(paste("match incomplete,",length(unmatched), sub, "datapoints unmatched"))
    
    if(sub == "spp.list" & trim.dat){m$data <- m$data[m$data$species %in% m$spp.list$species,]}
    if(ignore.unmatched){if(sub == "data"){m$data <- m$data[!m$data$species %in% unmatched,]}}else{
      
      
      # if all match pair datasets checked and species remain unmatched write manual match spp list
      dir.create(paste(input.folder, "r data/", sep = ""), showWarnings = F)
      dir.create(paste(input.folder, "r data/match data/", sep = ""), showWarnings = F)
      
      write.csv(data.frame(synonyms = if(sub == "spp.list"){""}else{unmatched},
                           species = if(sub == "spp.list"){unmatched}else{""}),
                paste("r data/match data/",m$data.ID," mmatch.csv", sep = ""),
                row.names = F)
      print(paste("unmatched species list saved in file 'Data.IDmmatch.csv'"))
      stop("manually match and save as 'Data.IDmmatched.csv' to continue")}}
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


