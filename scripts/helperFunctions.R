# Includes custom functions that do not require methylKit
# Operates on general objects that may be shared across different frameworks

#------------------------------------------------------------------------------------#

# Functions:

# getObject             Returns data object 
# getList               Returns list object
# getSampleFiles        Returns file names with full directory location 
# checkInputs           Validates input parameters
# checkConfifuration    Verifies project configuration file exists
# checkSamples          Verifies sample info file exists

#------------------------------------------------------------------------------------#

# Load Libraries

#------------------------------------------------------------------------------------#

############------------------------------------#
# FUNCTION # returns list object #
############------------------------------------#

getList <- function(temp) {
  return(lapply(temp, function(x) x))
}

############------------------------------------#
# FUNCTION # verifies sample files exist and returns list of files #
############------------------------------------#

getSampleFiles <- function(samples, inputDirectory) {
  files = list()
  for(i in samples) {
    temp <- paste0(inputDirectory,i)
    temp <- paste0(temp,"_val_1_bismark_bt2_pe.bismark.cov.gz")
    if(!file.exists(temp)) {
        print(paste0("Sample file not found: ", temp))
        quit()
    } else {
      files <- append(files,temp)
    }
  }
  return(files)
}

############------------------------------------#
# FUNCTION # validates input parameters #
############------------------------------------#

checkInputs <- function(args, executionConfiguration) {

  if(!file.exists(args[2])) {
	  print("output directory not found, creating one now")
	  dir.create(executionConfiguration$output)
  }

  if(executionConfiguration$plots == FALSE && executionConfiguration$tables == FALSE && executionConfiguration$all == FALSE) {
    print("No data output specified")
    quit()
  }

  temp <- list.files(args[1], "*.cov.gz", full=T)
  if(!length(temp) > 1) {
    print("No input files found")
    quit()
  }
}

############------------------------------------#
# FUNCTION # Verifies project configuration file exists #
############------------------------------------#

checkConfiguration <- function() {
  if(!file.exists("proj_config.yaml")) {
    print("Project configuration file not found")
    quit()
  } else {
      return(yaml::read_yaml("proj_config.yaml"))
  }
}

############------------------------------------#
# FUNCTION # Verifies sample info file exists #
############------------------------------------#

checkSamples <- function() {
  if(!file.exists("data/samples.info")) {
    print("Sample information file not found")
    quit()
  } else {
      return(suppressWarnings(read.csv("data/samples.info", header=TRUE)))
  }
}
