# Feb, 2018

source("Code/TransitionAnalysis_4way.r")
user_beta_4w <- user_beta

source("Code/TransitionAnalysis_BB.r")
user_beta_BB <- user_beta

source("Code/TransitionAnalysis_DD.r")
user_beta_DD <- user_beta

#### Read in data from EXCEL file ####

sheet <- "4-way"
zTrans <- F

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
  "PriorInclusionSIS",
  "SocialRisk",
  "MedicalRisk",
  "Published",
  "Interrater Reliability",
  "PriorInclusionODO",
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
  measureType = "p",
  measureCol = "AA Proportion",
  n = "AA n",
  grping1 = "Interval",
  grping2 = "StudyID",
  published = "Published",
  PA = "PA",
  zTrans = zTrans,
  user_beta = lapply(user_beta_4w, function(x) x$A)
)
RVEModel$data <- dataRaw[, c("StudyID", "Authors", "Sample", "Subsample", "Interval")] %>% left_join(RVEModel$data)
RVEModel$data$grp1 <- RVEModel$data$Interval
RVEModel$data$grp2 <- RVEModel$data$StudyID

forest.RVE(
  RVEModel,
  studyNames = "Authors",
  sampleNames = "Subsample",
  effSizeName = "Percentage",
  fileName = "Figures/AA_p_Forest.tiff"
)

RVEModel <- MetaAnalysis.RVE(
  data = dataRaw,
  measureType = "p",
  measureCol = "BB Proportion",
  n = "BB n",
  grping1 = "Interval",
  grping2 = "StudyID",
  published = "Published",
  PA = "PA",
  zTrans = zTrans,
  user_beta = lapply(user_beta_4w, function(x) x$B)
)
RVEModel$data <- dataRaw[, c("StudyID", "Authors", "Sample", "Subsample", "Interval")] %>% left_join(RVEModel$data)
RVEModel$data$grp1 <- RVEModel$data$Interval
RVEModel$data$grp2 <- RVEModel$data$StudyID

forest.RVE(
  RVEModel,
  studyNames = "Authors",
  sampleNames = "Subsample",
  effSizeName = "Percentage",
  fileName = "Figures/BB_p_Forest.tiff"
)

RVEModel <- MetaAnalysis.RVE(
  data = dataRaw,
  measureType = "p",
  measureCol = "CC Proportion",
  n = "CC n",
  grping1 = "Interval",
  grping2 = "StudyID",
  published = "Published",
  PA = "PA",
  zTrans = zTrans,
  user_beta = lapply(user_beta_4w, function(x) x$C)
)
RVEModel$data <- dataRaw[, c("StudyID", "Authors", "Sample", "Subsample", "Interval")] %>% left_join(RVEModel$data)
RVEModel$data$grp1 <- RVEModel$data$Interval
RVEModel$data$grp2 <- RVEModel$data$StudyID

forest.RVE(
  RVEModel,
  studyNames = "Authors",
  sampleNames = "Subsample",
  effSizeName = "Percentage",
  fileName = "Figures/CC_p_Forest.tiff"
)

RVEModel <- MetaAnalysis.RVE(
  data = dataRaw,
  measureType = "p",
  measureCol = "DD Proportion",
  n = "DD n",
  grping1 = "Interval",
  grping2 = "StudyID",
  published = "Published",
  PA = "PA",
  zTrans = zTrans,
  user_beta = lapply(user_beta_4w, function(x) x$D)
)
RVEModel$data <- dataRaw[, c("StudyID", "Authors", "Sample", "Subsample", "Interval")] %>% left_join(RVEModel$data)
RVEModel$data$grp1 <- RVEModel$data$Interval
RVEModel$data$grp2 <- RVEModel$data$StudyID

forest.RVE(
  RVEModel,
  studyNames = "Authors",
  sampleNames = "Subsample",
  effSizeName = "Percentage",
  fileName = "Figures/DD_p_Forest.tiff"
)

RVEModel <- MetaAnalysis.RVE(
  data = dataRaw,
  measureType = "p",
  measureCol = "ISIS Proportion",
  n = "ISIS n",
  grping1 = "Interval",
  grping2 = "StudyID",
  published = "Published",
  PA = "PA",
  zTrans = zTrans,
  user_beta = lapply(user_beta_BB, function(x) x$NotB)
)
RVEModel$data <- dataRaw[, c("StudyID", "Authors", "Sample", "Subsample", "Interval")] %>% left_join(RVEModel$data)
RVEModel$data$grp1 <- RVEModel$data$Interval
RVEModel$data$grp2 <- RVEModel$data$StudyID

forest.RVE(
  RVEModel,
  studyNames = "Authors",
  sampleNames = "Subsample",
  effSizeName = "Percentage",
  fileName = "Figures/ISIS_p_Forest.tiff"
)

RVEModel <- MetaAnalysis.RVE(
  data = dataRaw,
  measureType = "p",
  measureCol = "OO Proportion",
  n = "OO n",
  grping1 = "Interval",
  grping2 = "StudyID",
  published = "Published",
  PA = "PA",
  zTrans = zTrans,
  user_beta = lapply(user_beta_DD, function(x) x$NotD)
)
RVEModel$data <- dataRaw[, c("StudyID", "Authors", "Sample", "Subsample", "Interval")] %>% left_join(RVEModel$data)
RVEModel$data$grp1 <- RVEModel$data$Interval
RVEModel$data$grp2 <- RVEModel$data$StudyID

forest.RVE(
  RVEModel,
  studyNames = "Authors",
  sampleNames = "Subsample",
  effSizeName = "Percentage",
  fileName = "Figures/OO_p_Forest.tiff"
)

RVEModel <- MetaAnalysis.RVE(
  data = dataRaw,
  measureType = "k",
  measureCol = "Corr",
  n = "n",
  grping1 = "Interval",
  grping2 = "StudyID",
  published = "Published",
  PA = "PA",
  zTrans = T
)

funnel.RVE(
  RVEModel,
  fileName = "Figures/Four_way_k_Funnel.tiff",
  trimnfill = F,
  plotReg = T,
  effSizeName = expression(paste("Cohen's", ~italic(kappa)))
)
