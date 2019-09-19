# Function to generate 'probeFile' probe-level annotations from platform design (pd) 
# packages from BioConductor.  The generated 'probeFile' data is as a data.table in 
# an .rda file.

# This function is for use on ClariomD/XTA Affymetrix arrays.

# These 'probeFiles' are then converted to a custom 'probePackage' using AnnotationForge.
# The necessary functions for data parsing / annotation package creation can be found on Github:
# https://github.com/harrisgm/GCSscore-probeFile-functions


# Examples using platform design packages on Bioconductor, found at:
# https://www.bioconductor.org/packages/release/data/annotation/

# probeFile.hta20 <- ClariomDXTApFBuilder(chip.pd = "pd.hta.2.0")
# probeFile.mta10 <- ClariomDXTApFBuilder(chip.pd = "pd.mta.1.0")
# probeFile.rta10 <- ClariomDXTApFBuilder(chip.pd = "pd.rta.1.0")
# probeFile.clariomdhuman <- ClariomDXTApFBuilder(chip.pd = "pd.clariom.d.human")

# # load necessary libaries:
# library(data.table)
# # library(devtools)
#   # NOTE: devtools not needed if not installing a package from GitHub
# library(stringr)
# library(RSQLite)

# Function:
ClariomDXTApFBuilder <- function(chip.pd = NULL) {
  # Install necessary pd.* package if not already installed:
  if (!requireNamespace(chip.pd, quietly = TRUE)){
    BiocManager::install(chip.pd)
    print("installing necessary platform design (pd) package from Bioconductor")}
  
  # get chip name (without periods/dots) from chip.pd:
  chip <- str_remove_all(strsplit(chip.pd,"pd.")[[1]][2],"[.]")
  
  # Load .sqlite file from pd.mta.1.0 package:
  chip.sqlite <- paste(chip.pd,".sqlite",sep ="")
  chip.db <- system.file("extdata",chip.sqlite,package = chip.pd)
  ## connect to db
  con <- dbConnect(drv=RSQLite::SQLite(), dbname=chip.db)
  ## list all tables
  tables <- dbListTables(con)
  ## exclude sqlite_sequence (contains table information)
  tables <- tables[tables != "sqlite_sequence"]
  chip.pdinfo <- vector("list", length=length(tables))
  names(chip.pdinfo) <- tables
  ## create a data.frame for each table
  for (i in seq(along=tables)) {
    chip.pdinfo[[i]] <- dbGetQuery(conn=con, statement=paste("SELECT * FROM '", tables[[i]], "'", sep=""))
  }
  # Extract probe-level data:
  chip.pmfeature <- as.data.table(chip.pdinfo[["pmfeature"]])
  # NOTE: chip.probefile.BioC doesn't contain the named "features", it ONLY has the numbered ones
  # NOTE: chip.core.mps has TCid names in it!
  chip.core.mps <- as.data.table(chip.pdinfo[["core_mps"]])
  # Add PSR/JUC ids here as well:
  chip.featureSet <- as.data.table(chip.pdinfo[["featureSet"]])
  # Cut featureSet down to just 2 columns:
  chip.featureSet <- chip.featureSet[,c("fsetid","man_fsetid")]
  
  setkey(chip.pmfeature,fsetid)
  setkey(chip.core.mps,fsetid)
  setkey(chip.featureSet,fsetid)
  
  # Add in info from core.mps and featureSet:
  chip.pmfeature <- chip.core.mps[chip.pmfeature]
  chip.pmfeature <- chip.featureSet[chip.pmfeature]
  setkey(chip.pmfeature,fid)
  
  # Get PM sequence information from chip.pd: ----------------------------
  
  data("pmSequence", package = chip.pd, envir = environment())
  chip.pmSeq <- as.data.table(pmSequence)
  # Necessary for hta2.0: must remove duplicated 'fid' rows or merge will add duplicated rows!!
  # NOTE:  unique vs. !duplicated --> ALWAYS use (!duplicated)
  # chip.test.unique<- chip.pmSeq[unique(fid)]
  # chip.test.notdup <- chip.pmSeq[!duplicated(fid)]
  # These objects have the SAME nrow, BUT only the chip.test.notdup will merge without extra rows!
  # Also the chip.unqiue is LARGER in megabytes --> implies NAs must have been introduced with the unique() command
  chip.pmSeq <- chip.pmSeq[!duplicated(fid)]
  setkey(chip.pmSeq,fid)
  
  chip.pmfeature <- chip.pmSeq[chip.pmfeature]
  # add in GC.count info for each sequence in the pmfeature object:
  chip.pmfeature[,GC.count := str_count(sequence, "G|C")]
  
  # Load in pd.chip.1.0 “netaffxTranscript" data for locustype/category---
  
  load(system.file("extdata","netaffxTranscript.rda",package = chip.pd))
  # creates "netaffxTranscript" object
  load(system.file("extdata","netaffxProbeset.rda",package = chip.pd))
  # creates "netaffxTranscript" object
  
  affynet.TC <- as.data.table(netaffxTranscript@data)

  
  if (chip %like% "clarioms"){
    affynet.TC <- affynet.TC[,c("transcriptclusterid","category")]
    setkey(affynet.TC,transcriptclusterid)
    setkey(chip.pmfeature,transcript_cluster_id)
    probeFile <- affynet.TC[chip.pmfeature]
  }
  if (chip %like% "mta10"| chip %like% "hta20" | chip %like% "rta10" | chip %like% "clariomdhuman"){
    affynet.TC <- affynet.TC[,c("transcriptclusterid","category")]
    setkey(affynet.TC,transcriptclusterid)
    setkey(chip.pmfeature,transcript_cluster_id)
    chip.pmfeature <- affynet.TC[chip.pmfeature]
    
    affynet.PSR <- as.data.table(netaffxProbeset@data)
    affynet.PSR <- affynet.PSR[,c("probesetid","locustype","probesettype")]
    # table(table(affynet.PSR$locustype))
    # NOTE: for hta2.0, many rows affynet.PSR has multiple locustype entries per row:
    # Remove all but first 'locustype' entry:
    affynet.PSR[,locustype := sapply(strsplit(affynet.PSR$locustype, "///"), "[", 1 )]
    # Assign key of affynet.PSR to 'probesetid' for PSR/JUC 'locustype' and 'probesettype' merge:
    setkey(affynet.PSR,probesetid)
    # Change key of pmfeature to 'man_fsetid' (PSR/JUC) for 'locustype' and 'probesettype' merge:
    setkey(chip.pmfeature,man_fsetid)
    # Create master probeFile for using chip.pd:
    probeFile <- affynet.PSR[chip.pmfeature]
  }
  # key the probeFile to sort by fid (or other column type):
  setkey(probeFile,fid)
  
  # Remove unnecessary columns, if present:
  if (!is.na(match("atom", names(probeFile)))){
    probeFile[,c("atom") := NULL]
  }
  # If "row_names" column exists (such as in the HTA_2-0), remove it:
  if (!is.na(match("row_names", names(probeFile)))){
    probeFile[,c("row_names") := NULL]
  }
  
  # savepath <- system.file("data",package = "Sscore2")
    save(probeFile, file = paste(system.file("data/",package = "Sscore2"),
                               paste("Ss2.",chip,".probeFile.rda",sep=""),sep=""),
         compress = "gzip", compression_level = 1)
  
  print(paste("'probeFile' created and stored in: ",paste("Ss2.",chip,".probeFile.rda",sep=""),sep=""), sep="")
  # print("'probeFile' data.table returned to working environment")
  return(probeFile)	
}
