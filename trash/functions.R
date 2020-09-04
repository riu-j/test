LoadCsvFiles = function(folder_path){
  data_before <<- 
    read.csv(paste(folder_path, "data_before.csv", sep=""), stringsAsFactors=FALSE, na.strings =".")
  offsetT <<- read.csv(paste(folder_path, "offsetT.csv", sep=""), stringsAsFactors=FALSE)
  info <<- read.csv(paste(folder_path, "info.csv", sep=""), stringsAsFactors=FALSE)
  iteration <<- read.csv(paste(folder_path, "iteration.csv", sep=""), sep="", header=TRUE, skip=1, stringsAsFactors=FALSE)
}

SetGraphicsDevice =function(file_extension, file_name){
  break_point <- 0
  switch(file_extension,
         "svg" = svg(paste(file_name, ".svg", sep=""), width=12, height=8, bg="transparent"),
         "pdf" = pdf(paste(file_name, ".pdf", sep=""), width=12, height=8, bg="transparent"),
         break_point <- 1
  )
  if(break_point == 1)stop("Error: The file extension is not specified or supported.")
}