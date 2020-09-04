# Overview
# This code make control file for NONMEM and copy fortran source code.
# 
# Input
# initial_parameter.csv
# data.csv
# 
# Output
# control_NONMEM.txt
# PRED_SReFT.f90

#Library
library("dplyr")

print_list <- function(list){
  len <- length(list)
  output <- c()
  for(i in 1:len){
    output[i] <- list[[i]]
  }
  return(output)
}

proj <- "s139"
pred_type <- "old_ana_mean"
init <- read.csv(paste("../notes/", proj, "/", proj, "_initialprms.csv", sep=""), row.names=1, na=".") %>% signif(digits=4)
data <- read.csv(paste("../notes/", proj, "/", proj, "_data.csv", sep=""))
dv_name <- colnames(init) %>% gsub("__", " ", .) %>% gsub("_", "-", .)
cols <- colnames(data)
table_cols <- cols[-grep("mean", cols)] %>% .[-grep("count", .)] %>% .[-2]
ctl <- list()
numbm <- length(dv_name)

ctl["problem"] <- paste("$PROBLEM ", Sys.time(), "\n", paste(";", dv_name, "\n", sep="") %>% paste(collapse=""), "\n", sep="")
ctl["input"] <- paste("$INPUT", paste(cols, collapse=" ") %>% paste("\n\n", sep=""), sep=" ")
ctl["data"] <- paste("$DATA ", proj, "_data.csv ignore=@\n\n", sep="")
ctl["subroutine"] <- paste("$SUBROUTINE PRED=", proj, "_source.f90\n\n", sep="")
ctl["theta"] <-
  paste("$THETA\n" ,
        paste(init[1, 1], "fixed", init[1, 2:numbm] %>% paste(collapse=" ") %>% paste("\n", sep="")),
        init[2,] %>% paste(collapse=" ") %>% paste("\n", sep=""),
        init[3,] %>% paste(collapse=" ") %>% paste("\n\n", sep=""), sep="")
ctl["omega"] <- paste("$OMEGA\n0 fixed", "(0.01)x", numbm-1, "\n0.0001 (0 fixed)x", numbm-1, "\n0 fixed (0.0001)x", numbm-1, "\n\n")
ctl["sigma"] <- paste("$SIGMA (0.01)x", numbm, "\n\n", sep="")
ctl["estimation"] <- 
  paste("$ESTIMATION METHOD=1 MAXEVAL=999999 PRINT=1 NOABORT SIGDIGITS=2 FILE=", proj, "_iteration7.csv\n",
        "$ESTIMATION METHOD=1 MAXEVAL=999999 PRINT=1 NOABORT SIGDIGITS=3 FILE=", proj, "_iteration6.csv\n",
        "$ESTIMATION METHOD=1 MAXEVAL=999999 PRINT=1 NOABORT SIGDIGITS=4 FILE=", proj, "_iteration5.csv\n",
        "$ESTIMATION METHOD=1 MAXEVAL=999999 PRINT=1 NOABORT SIGDIGITS=5 FILE=", proj, "_iteration4.csv\n",
        "$ESTIMATION METHOD=1 MAXEVAL=999999 PRINT=1 NOABORT SIGDIGITS=4 FILE=", proj, "_iteration3.csv\n",
        "$ESTIMATION METHOD=1 MAXEVAL=999999 PRINT=1 NOABORT SIGDIGITS=3 FILE=", proj, "_iteration2.csv\n",
        "$ESTIMATION METHOD=1 MAXEVAL=999999 PRINT=1 NOABORT SIGDIGITS=2 FILE=", proj, "_iteration.csv FORMAT=,PE15.5 NOTITLE=1\n\n", sep="")
ctl["covariance"] <- "$COVARIANCE UNCONDITIONAL MATRIX=S\n\n"
ctl["table"] <- paste("$TABLE ", paste(table_cols, collapse=" "), "CIPRED CWRES\nNOPRINT FORMAT=,F10.5 FILE=s137_table.csv ONEHEADER NOTITLE NOAPPEND", sep="")

cat(print_list(ctl), file=paste("../notes/", proj, "/", proj, "_control.txt", sep=""), sep="")

file.copy(from=paste("../src_pred/", pred_type, ".f90", sep=""), to=paste("../notes/", proj, "/", sep=""))
file.rename(from=paste("../notes/", proj, "/", pred_type, ".f90", sep=""), to=paste("../notes/", proj, "/", proj, "_source.f90", sep=""))
