# April, 2017

source("Code/Utilities/init.RVE.r")

#### Read in data from EXCEL file ####

sheet <- "4-way"
zTrans <- T

dataRaw <- read_excel("Data/20180121_Data_Extraction_Table.xlsx", sheet = sheet, col_names=TRUE, na = "NA")
dataRaw <- dataRaw[rowSums(is.na(dataRaw)) != ncol(dataRaw),]

# Rename data columns
dataNames <- c(
  "StudyID",
  "Authors",
  "Sample",
  "Subsample",
  "Interval",
  "n",
  "Corr",
  "PA",
  "Year",
  "Country",
  "CodingMethod",
  "PriorInclusion",
  "SocialRisk",
  "MedicalRisk",
  "Published",
  "Interrater Reliability",
  "PriorInclusion",
  "Gender"
)
names(dataRaw)[1:length(dataNames)] <- dataNames
numericCols <- c(1, 3, 6, 7, 8, 9, 12, 13, 14, 15, 16, 17)
dataRaw[, numericCols] <- lapply(dataRaw[, numericCols], as.numeric)

# Fix some weird formatting in front and following names
dataRaw$Authors <- trimws(gsub("\r\n", "", dataRaw$Authors))

# Get the unique intervals we have in this data (these are hard-coded as we want to add
# an empty section in the final plot for those intervals which aren't in the data)
intervals <- c(
  "Infancy - Toddlerhood",
  "Infancy - Preschool",
  "Infancy - School entry",
  "Toddlerhood - Preschool",
  "Toddlerhood - School entry",
  "Preschool - School entry"
)
dataRaw$StudyID <- match(dataRaw$StudyID,unique(dataRaw$StudyID)) # convert to consequtive IDs
dataRaw$Interval <- factor(dataRaw$Interval, levels <- intervals)

# Add symbols indicating if studies were published or included prior
dataRaw$Authors[!as.logical(dataRaw$Published)] <- paste0(dataRaw$Authors[!as.logical(dataRaw$Published)],"^")

#### Build random-effects models ####

RVEModel <- MetaAnalysis.RVE(
  data = dataRaw,
  measureType = "k",
  measureCol = "Corr",
  n = "n",
  grping1 = "Interval",
  grping2 = "StudyID",
  published = "Published",
  PA = "PA",
  zTrans = zTrans
)

Summary.RVE(
  RVEModel,
  fileName = "Output/StatisticalResults_4x4.txt"
)

forest.RVE(
  RVEModel,
  studyNames = "Authors",
  sampleNames = "Subsample",
  effSizeName = "Cohen's k",
  fileName = "Figures/4x4_Forest.tiff"
) 

funnel.RVE(
  RVEModel,
  fileName = "Figures/4x4_Funnel.tiff",
  trimnfill = T,
  plotReg = T
)

if (zTrans) {
  transf <- tanh
} else {
  transf <- function(x) {x}
}
correlations <- lapply(RVEModel$REModels, function(x) transf(x$reg_table$b.r))
save(correlations, file = "Output/STAB/4-way/Correlations.RData")
