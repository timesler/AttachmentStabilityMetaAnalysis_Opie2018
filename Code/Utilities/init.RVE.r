#### Load 3rd party packages ####

if (length(find.package("metafor", quiet = T)) == 0) {install.packages("metafor")}
if (length(find.package("robumeta", quiet = T)) == 0) {install.packages("robumeta")}
if (length(find.package("readxl", quiet = T)) == 0) {install.packages("readxl")}
if (length(find.package("dplyr", quiet = T)) == 0) {install.packages("dplyr")}
if (length(find.package("data.table", quiet = T)) == 0) {install.packages("data.table")}
if (length(find.package("nnet", quiet = T)) == 0) {install.packages("nnet")}
if (length(find.package("effects", quiet = T)) == 0) {install.packages("effects")}
if (length(find.package("meta", quiet = T)) == 0) {install.packages("meta")}
library(meta)
library(metafor)
library(robumeta)
library(readxl)
library(plyr)
library(dplyr)
library(data.table)
library(nnet)
library(effects)
library(meta)

#### Load utility functions ####
source("Code/Utilities/RVE.r")
source("Code/Utilities/forest.col.r")
source("Code/Utilities/trimfill.default.r")
source("Code/Utilities/robu.custom.r")
