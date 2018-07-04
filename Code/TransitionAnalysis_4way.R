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

# Get relevant columns
data <- dataRaw[
  ,
  c(
    "Study.ID",
    "Subsample.ID",
    "Interval",
    "Subsample.Desciption",
    grep(identifier, names(dataRaw), value = T)
  )
  ]
data[,grep(identifier, names(dataRaw), value = T)] <- lapply(
  data[,grep(identifier, names(dataRaw), value = T)],
  as.numeric
)
names(data)[ncol(data)] <- "n"

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
file <- "Output/CONT/4way_Results.csv"
fileFig <- "Output/CONT/4way_Results_figure.csv"
write.table(
  "4-way Transitions Tables", 
  file = file, 
  row.names = F, 
  col.names = F, 
  sep = ","
)
write.table(
  "4-way Transitions Figure", 
  file = fileFig, 
  row.names = F, 
  col.names = F, 
  sep = ","
)

matrices <- list()
matricesFig <- list()
user_beta <- list()
metrics <- c(".Count", ".Percentage", ".SE", ".CI_L", ".CI_U", ".ExpectedCount", ".ExpectedPercent",
             ".PercentDeviation", ".StandardizedResidual", ".AdjStandResidual")
metricsFig <- c(".Count", ".Percentage", ".CI_L", ".CI_U", ".ExpectedPercent", ".AdjStandResidual")
for (interval in intervals)
{
  user_beta[[interval]] <- list()
  
  mat <- matrix(
    0, 
    nrow = length(classes)*length(metrics) + 1,
    ncol = length(classes) + 1,
    dimnames = list(
      c(paste0(
        rep(classes, each = length(metrics)), 
        metrics), ".CountTotal"),
      c(classes, "")
    )
  )
  for (class1 in classes)
  {
    if (interval != "overall")
      data.int <- data[data$Interval == interval,]
    else
      data.int <- data
    data.weights <- k_RVEModel$dataList[[interval]] %>%
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
    
    class1.n <- rowSums(data.int[,grep(paste0("^", class1), names(data.int))])
    data.int <- data.int[class1.n > 0, ]
    
    data.int <- data.int %>%
      group_by(Study.ID) %>%
      arrange(Study.ID, -n) %>%
      mutate(w.cnts = as.integer(row_number() < 2))
    
    for (class2 in classes)
    {
      col <- paste0(class1, ".", class2, identifier)
      cnt <- sum(data.int[[col]] * data.int$w.cnts)
      
      mat[paste0(class1, ".Count"), class2] <- cnt
      mat[paste0(class1, ".Count"), ncol(mat)] <- cnt +
        mat[paste0(class1, ".Count"), ncol(mat)]
      mat[nrow(mat), class2] <- cnt +
        mat[nrow(mat), class2]
      mat[nrow(mat), ncol(mat)] <- cnt +
        mat[nrow(mat), ncol(mat)]
      
      if (fullRobu)
      {
        data.int$Meas <- data.int[[col]] / rowSums(
          data.int[,grep(paste0("^", class1), names(data.int), value = T)]
        )
        
        robuMdl <- robu.custom(
          formula      = Meas ~ 1,
          data         = data.int,
          studynum     = data.int$Study.ID,
          var.eff.size = 1/data.int$weight,
          rho          = 0.8,
          small        = T,
          modelweights = "CORR",
          userweights  = data.int$weight^2,
          userdfs      = k_RVEModel$REModels[[interval]]$reg_table$dfs
        )
        mat[paste0(class1, ".Percentage"), class2] <- robuMdl$reg_table$b.r*100
        mat[paste0(class1, ".SE"), class2] <- robuMdl$reg_table$SE*100
        mat[paste0(class1, ".CI_L"), class2] <- max(0, robuMdl$reg_table$CI.L*100)
        mat[paste0(class1, ".CI_U"), class2] <- min(100, robuMdl$reg_table$CI.U*100)
      }
    }
  }
  
  if (!fullRobu)
  {
    for (class1 in classes)
    {
      n <- mat[paste0(class1, ".Count"), ncol(mat)]
      props <- mat[paste0(class1, ".Count"), -ncol(mat)] / n
      SEs <- sqrt(props * (1 - props) / n)
      CI_L <- props - qnorm(.975) * SEs
      CI_U <- props + qnorm(.975) * SEs
      
      mat[paste0(class1, ".Percentage"), -ncol(mat)] <- props*100
      mat[paste0(class1, ".SE"), -ncol(mat)] <- SEs*100
      mat[paste0(class1, ".CI_L"), -ncol(mat)] <- CI_L*100
      mat[paste0(class1, ".CI_U"), -ncol(mat)] <- CI_U*100
    }
  }
  
  
  for (class1 in classes)
  {
    for (class2 in classes)
    {
      N <- mat[nrow(mat), ncol(mat)]
      cnt <- mat[paste0(class1, ".Count"), class2]
      row.cnt <- mat[paste0(class1, ".Count"), ncol(mat)]
      col.cnt <- mat[nrow(mat), class2]
      
      exp.cnt <- row.cnt * col.cnt / N
      exp.prop <- exp.cnt / row.cnt
      pc.dev <- (cnt - exp.cnt) / exp.cnt * 100
      std.res <- (cnt - exp.cnt) / sqrt(exp.cnt)
      adj.std.res <- std.res / sqrt((1 - exp.cnt / row.cnt) * (1 - exp.cnt / col.cnt))
      
      mat[paste0(class1, ".ExpectedCount"), class2] <- exp.cnt
      mat[paste0(class1, ".ExpectedPercent"), class2] <- exp.prop * 100
      mat[paste0(class1, ".PercentDeviation"), class2] <- pc.dev
      mat[paste0(class1, ".StandardizedResidual"), class2] <- std.res
      mat[paste0(class1, ".AdjStandResidual"), class2] <- adj.std.res
    }
  }
  matText <- matrix(
    "",
    nrow = length(classes)*length(metrics) + 1,
    ncol = length(classes) + 1,
    dimnames = list(
      c(paste0(
        rep(classes, each = length(metrics)),
        metrics), ".CountTotal"),
      c(classes, "")
    )
  )
  matText[grep("\\.Count", rownames(matText)), ] <- formatC(
    mat[grep("\\.Count", rownames(mat)), ], digits = 3, format = "f"
  )
  matText[grep("\\.Count", rownames(matText), invert = T), ] <- formatC(
    mat[grep("\\.Count", rownames(mat), invert = T), ], digits = 3, format = "f"
  )
  
  matTextFig <- matText[metrics %in% metricsFig, ]
  CI <- paste0("[", paste0(
    matTextFig[metricsFig == ".CI_L",-ncol(matTextFig)], ", ",
    matTextFig[metricsFig == ".CI_U",-ncol(matTextFig)]
  ), "]")
  matTextFig[metricsFig == ".CI_L", -ncol(matTextFig)] <- CI
  row.names(matTextFig)[metricsFig == ".CI_L"] <- paste0(row.names(matTextFig)[metricsFig == ".CI_L"], "U")
  matTextFig <- matTextFig[metricsFig != ".CI_U", ]
  
  matrices[[interval]] <- matText
  matricesFig[[interval]] <- matTextFig
  
  write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
  write.table(interval, file, row.names = F, col.names = F, sep = ",", append = T)
  write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
  
  write.table("", fileFig, row.names = F, col.names = F, sep = ",", append = T)
  write.table(interval, fileFig, row.names = F, col.names = F, sep = ",", append = T)
  write.table("", fileFig, row.names = F, col.names = F, sep = ",", append = T)
  
  write.table(matrix(c("", colnames(matText)), nrow = 1), 
              file, row.names = F, col.names = F, sep = ",", append = T)
  write.table(matText, file, row.names = T, col.names = F, sep = ",", append = T)
  
  write.table(matrix(c("", colnames(matTextFig)), nrow = 1), 
              fileFig, row.names = F, col.names = F, sep = ",", append = T)
  write.table(matTextFig, fileFig, row.names = T, col.names = F, sep = ",", append = T)
  
  for (class1 in classes)
  {
    if (interval != "overall")
      data.int <- data[data$Interval == interval,]
    else
      data.int <- data
    data.weights <- k_RVEModel$dataList[[interval]] %>%
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
    
    class1.n <- rowSums(data.int[,grep(paste0("^", class1), names(data.int))])
    data.int <- data.int[class1.n > 0, ]
    
    # data.int <- data.int %>%
    #   group_by(Study.ID) %>%
    #   mutate(weight = 1 / n())
    
    user_beta[[interval]][[class1]] <- list()
    
    user_beta[[interval]][[class1]]$b.r         <- mat[paste0(class1, ".Percentage"), class1]
    user_beta[[interval]][[class1]]$SE          <- mat[paste0(class1, ".SE"), class1]
    user_beta[[interval]][[class1]]$CI.L        <- mat[paste0(class1, ".CI_L"), class1]
    user_beta[[interval]][[class1]]$CI.U        <- mat[paste0(class1, ".CI_U"), class1]
    user_beta[[interval]][[class1]]$mdlWeights  <- data.int$weight/sum(data.int$weight)*100
  }
}

