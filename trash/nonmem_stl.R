a <- c(1, 2, 3)
b <- c(4, 5, 6)
x <- sprintf("%d\n", a)
y <- print(a)
z <- print(b)
cat(a, file="test.txt")

library("tidyverse")

cols <- c("a", "b")
cols <- paste(cols, collapse=" ") %>% paste("\n\n", sep="")
proj <- "s137"
ctl <- list()

ctl["problem"] <- "$PROBLEM\n\n"
ctl["input"] <- paste("$INPUT", cols, sep=" ")
data <- paste("$DATA ", proj, "_data.csv ignore=@\n\n", sep="")
subroutine <- paste("$SUBROUTINE PRED=", proj, "_source.f90\n\n", sep="")
theta
omega
sigma <- paste("$SIGMA (0.01)x", numbm, sep="")
estimation
covariance <- "$COVARIANCE UNCONDITIONAL MATRIX=S\n\n"
table <- "$TABLE ID TIME BM DV CIPRED CWRES\nNOPRINT FORMAT=,F10.5 FILE=s137_table.csv ONEHEADER NOTITLE NOAPPEND"
cat(print_list(ctl), file="text.txt", sep="")

print_list <- function(list){
  len <- length(list)
  output <- c()
  for(i in 1:len){
    output[i] <- list[[i]]
  }
  return(output)
}

test <- function(list){
    list[[1]]
}
