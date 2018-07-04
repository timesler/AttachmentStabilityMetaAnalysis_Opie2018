# !diagnostics off
# April, 2017

source("Code/MetaAnalysis_4way_k.r")

k_RVEModel <- RVEModel

#### Read in data from EXCEL file ####

sheet <- "Processed - Stability"
classes <- c("B", "A", "C", "D")
identifier <- ".4"
fullRobu <- T

dataRaw <- read_excel("Data/20180121_Data_Extraction_Table.xlsx", sheet = sheet, col_names=TRUE, na = "NA")
dataRaw <- dataRaw[rowSums(is.na(dataRaw)) <= ncol(dataRaw)/2,]
names(dataRaw) <- make.names(names(dataRaw))
dataRaw$Study.ID <- match(dataRaw$Study.ID,unique(dataRaw$Study.ID))

dataRaw$Risk <- as.logical(dataRaw$Medical.Risk) | as.logical(dataRaw$Social.Risk)
dataRaw$Gender[dataRaw$Gender != "Female"] <- "Mixed"
dataRaw$IsAinsworthAinsworth <- dataRaw$Coding.Method == "Ainsworth/Ainsworth"
dataRaw$IsCassidyMarvinCassidyMarvin <- dataRaw$Coding.Method == "Cassidy-Marvin/Cassidy-Marvin"
dataRaw$IsAinsworthCassidyMarvin <- dataRaw$Coding.Method == "Ainsworth/Cassidy-Marvin"
dataRaw$IsAinsworthMainCassidy <- dataRaw$Coding.Method == "Ainsworth/Main-Cassidy"
dataRaw$YearBinary <- dataRaw$Year >= 2000
dataRaw$IRRBinary <- dataRaw$Interrater.Reliability...PA >= 0.85

moderators <- c(
  "Risk",
  "Year",
  "YearBinary",
  "Country",
  "Included.in.Prior.Meta.analysis.S.IS",
  "Social.Risk",
  "Medical.Risk",
  "Published",
  "Interrater.Reliability...PA",
  "IRRBinary",
  "Included.in.Prior.Meta.analysis.D.O",
  "Gender",
  "IsAinsworthAinsworth",
  "IsCassidyMarvinCassidyMarvin",
  "IsAinsworthCassidyMarvin",
  "IsAinsworthMainCassidy"
)

# Get relevant columns
data <- dataRaw[
  ,
  c(
    "Study.ID",
    "Subsample.ID",
    "Interval",
    "Subsample.Desciption",
    grep(identifier, names(dataRaw), value = T),
    moderators
  )
  ]
data[,grep(identifier, names(dataRaw), value = T)] <- lapply(
  data[,grep(identifier, names(dataRaw), value = T)],
  as.numeric
)
names(data)[grep(paste0("n", identifier), names(data))] <- "n"

# Filter to get just samples with contingency tables
data <- data[!is.na(data$n) & data$n > 0,]

# Get the unique intervals we have in this data (these are hard-coded as we want to add
# an empty section in the final plot for those intervals which aren't in the data)
intervals <- c(
  "Infancy - Toddlerhood",
  "Infancy - Preschool",
  "Infancy - School entry",
  "Toddlerhood - Preschool",
  "Toddlerhood - School entry",
  "Preschool - School entry",
  "overall"
)
data$Interval <- factor(data$Interval, levels <- intervals)

#### Build transitions tables ####
file <- "Output/CONT/4way_Results_Moderators.csv"
write.table(
  "4-way Moderator Results", 
  file = file, 
  row.names = F, 
  col.names = F, 
  sep = ","
)
mat <- matrix(
  0, 
  nrow = length(moderators),
  ncol = length(classes),
  dimnames = list(
    moderators,
    classes
  )
)
for (class in classes)
{
  data.int <- data
  data.weights <- k_RVEModel$dataList$overall %>%
    select(mdlWeights, StudyID, Sample, Subsample, Interval, n) %>%
    mutate(mdlWeights = sqrt(mdlWeights) / sum(sqrt(mdlWeights)) * length(unique(StudyID)))
  data.int <- data.int %>%
    left_join(
      data.weights,
      by = c(
        "Subsample.ID"="Sample",
        "Subsample.Desciption"="Subsample",
        "Interval"="Interval",
        "n"
      )
    ) %>%
    rename(weight = mdlWeights)
  
  class.n <- rowSums(data.int[,grep(paste0("^", class, "\\."), names(data.int))])
  data.int <- data.int[class.n > 0, ]
  
  col <- paste0(class, ".", class, identifier)
  cnt <- sum(data.int[[col]] * data.int$w.cnts)
  
  for (moderator in moderators)
  {
    if (fullRobu)
    {
      data.int$Meas <- data.int[[col]] / rowSums(
        data.int[,grep(paste0("^", class, "\\."), names(data.int), value = T)]
      )
      
      robuMdl <- robu.custom(
        formula      = as.formula(paste("Meas ~ 1 + ", moderator)),
        data         = data.int,
        studynum     = data.int$Study.ID,
        var.eff.size = 1/data.int$weight,
        rho          = 0.8,
        small        = T,
        modelweights = "CORR",
        userweights  = data.int$weight^2,
        userdfs      = k_RVEModel$REModels[[interval]]$reg_table$dfs
      )
      
      mat[moderator, class] <- round(robuMdl$reg_table$prob[2],3)
    }
  }
}

write.table(
  mat, 
  file = file, 
  row.names = T, 
  col.names = T, 
  sep = ",",
  append = T
)

