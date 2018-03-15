# These are functions to read a directory of Z2 files and return
# one big data.frame with the file contents, ID'd by the filename.
# The first one is to parse the Z2 file. The next one just calls
# this for everything in a directory, but you really oughta handle
# this differently. But whatever.

read_a_coulter_file <- function(path) {
  this_datar <- list()
  file_slurp <- readLines(path) # Now a character vector
  # 
  parameters <- file_slurp[ (grep("\\[Z2\\]",file_slurp)+1): # Find where this is
                (grep("\\[#Extra\\]",file_slurp)-1) # encoded in file
      ]
  #
  dilution_factor <- unique(as.numeric( 
    sub("DilF=\\s*","", # Make dilution factor numeric
      parameters[grep("^DilF=",parameters)] # Find dilution factor
      )
    ))
  if ( length(dilution_factor) > 1 ) { 
    return("Error: couldn't parse the dilution factor right")
    }
  #
  volume <- unique(as.numeric( 
    sub("Vol=\\s*","", # Make dilution factor numeric
      parameters[grep("^Vol=",parameters)] # Find dilution factor
      )
    ))
  if ( length(volume) > 1 ) { 
    return("Error: couldn't parse the volume right")
    }
  #
  bin_center <- as.numeric(
    file_slurp[ (grep("Bindiam",file_slurp)+1): # Find where this is
                (grep("Binunits",file_slurp)-1) # encoded in file
      ]
    )
  #
  bin_count <- as.numeric(
    file_slurp[ (grep("Binheight",file_slurp)+1): # Find where this is
                (grep("\\[end\\]",file_slurp)-1) # encoded in file
      ]
    )
  #
  return(
    data.frame(Path=path,BinCenter=bin_center,BinCount=bin_count,
      DilutionFactor=dilution_factor,Volume=volume)
    )
}

read_directory_of_Z2 <- function(data_dir,file_extension_pattern="*Z2") {
  return_df <- NULL
  paths_list <- list.files(path=data_dir,
    pattern=file_extension_pattern,full.names=T)
  #
  for (each_path in paths_list) {
  	tmp_df <- read_a_coulter_file(path=each_path)
    if (is.null(return_df)) {
      return_df <- tmp_df
    } else {
      return_df <- rbind(return_df,tmp_df)
    }
  }
  return(return_df)
}
