MetaAnalysis.RVE <- function(
  data, 
  measureType = "r",
  measureCol = "Corr", 
  n = "n", 
  grping1, grping2,
  published = "Published",
  PA = NULL,
  zTrans = F,
  user_beta = NULL)
{
  notes <- c()
  
  # Remove empty samples
  data$n <- data[[n]]
  data <- subset(data, data$n > 0)
  
  # Define grouping, effect-size and published columns
  data$grp1 <- data[[grping1]]
  data$grp2 <- data[[grping2]]
  data$Meas <- data[[measureCol]]
  MeasInit <- data$Meas
  data$published <- as.logical(data[[published]])
  grping1_lvls <- levels(data$grp1)
  
  nName <- n
  
  # Get sample variances
  if (zTrans) {
    data <- escalc(
      ri        = Meas,
      ni        = n, 
      measure   = "ZCOR",
      data      = data,
      var.names = c("Meas", "Var")
    )
  } else if (measureType == "r") {
    data <- escalc(
      ri        = Meas,
      ni        = n, 
      measure   = "COR",
      data      = data,
      var.names = c("Meas", "Var")
    )
  } else if (measureType == "k") {
    pci <- (data[[PA]] - data$Meas) / (1 - data$Meas)
    data$Var <- data[[PA]] * (1 - data[[PA]]) / ((1 - pci) ^ 2 * data$n)
  } else if (measureType == "p") {
    data$Var <- data$Meas * (1 - data$Meas) / (data$n + 1e-4)
  } else {
    stop("Unknown 'measureType' specified.")
  }
  
  # Remove any undefined data (and take notes)
  if (any(is.na(data$Meas))) {
    notes <- c(notes, paste0(
      "\tThe following rows were removed because the effect-size was undefined:\n\n"))
    notes <- c(notes, paste0(
      "\t * Grp 1 - ", data$grp1[is.na(data$Meas)], ", Grp 2 - ", data$grp2[is.na(data$Meas)],
      " (Original effect-size: ", MeasInit[is.na(data$Meas)], ")\n"))
    notes <- c(notes, "\n")
    data <- data[!is.na(data$Meas),]
  }
 
  # Make sure none of variances are 0, otherwise analysis will fail
  data$Var[data$Var == 0] <- 1e-9
  
  # Build random-effects models using Robust Variance Estimation
  # Also run sensitivity analysis on rho and conduct Egger's regression test for bias
  REModels <- list()
  rhoSensitivity <- list()
  eggersRegTest <- list()
  dataList <- list()
  for (grping in grping1_lvls) {
    if (sum(data$grp1 == grping) > 1)
    {
      data.grp = data[data$grp1 == grping,]
      
      # Calculate random effects model
      argList <- list(
        formula      = Meas ~ 1,
        data         = data.grp,
        studynum     = data.grp$grp2,
        var.eff.size = data.grp$Var,
        rho          = 0.8,
        small        = T,
        modelweights = "CORR"
      )
      REModels[[grping]] <- do.call(robu, argList)
      dataList[[grping]] <- data.grp
      
      if (!is.null(user_beta))
      {
        REModels[[grping]]$reg_table$b.r  <- user_beta[[grping]]$b.r
        REModels[[grping]]$reg_table$SE   <- user_beta[[grping]]$SE
        REModels[[grping]]$reg_table$CI.L <- user_beta[[grping]]$CI.L
        REModels[[grping]]$reg_table$CI.U <- user_beta[[grping]]$CI.U
        REModels[[grping]]$data.full$r.weights <- user_beta[[grping]]$mdlWeights
      }
      
      dataList[[grping]]$mdlWeights <- REModels[[grping]]$data.full$r.weights
      dataList[[grping]]$mdlWeights2 <- REModels[[grping]]$data.full$weights
  
      # Perform sensitivity analysis
      tmp <- capture.output(rhoSensitivity[[grping]] <- sensitivity(REModels[[grping]]))
      
      if (sum(data.grp$published) > 2 && measureType != "p") {
        # Conduct regression test
        eggersRegTest[[grping]] <- robu(
          formula      = Meas ~ 1 + sqrt(Var),
          data         = data.grp[data.grp$published,],
          studynum     = grp2,
          var.eff.size = Var,
          rho          = 0.8,
          small        = T,
          modelweights = "CORR",
          userweights  = REModels[[grping]]$data.full$r.weights[data.grp$published]
        )
      } else {
        eggersRegTest[[grping]] <- list(reg_table = list(prob = c(NaN,NaN)))
      }
    } else {
      data.grp = data[data$grp1 == grping,]
      dataList[[grping]] <- data.grp
      if (nrow(data.grp) > 0) dataList[[grping]]$mdlWeights <- 1
      notes <- c(notes, paste0("\tSkipped '", grping, "' - ", sum(data$grp1 == grping), " studies.\n"))
    }
  }
  
  crossTab <- table(data$grp1)
  data.grp <- data[data$grp1 %in% names(crossTab)[crossTab > 0], ]
  argList <- list(
    formula      = Meas ~ 1,
    data         = data.grp,
    studynum     = data.grp$grp2,
    var.eff.size = data.grp$Var,
    rho          = 0.8,
    small        = T,
    modelweights = "CORR"
  )
  REModels$overall <- do.call(robu, argList)
  dataList$overall <- data.grp
  
  if (!is.null(user_beta))
  {
    REModels$overall$reg_table$b.r  <- user_beta$overall$b.r
    REModels$overall$reg_table$SE   <- user_beta$overall$SE
    REModels$overall$reg_table$CI.L <- user_beta$overall$CI.L
    REModels$overall$reg_table$CI.U <- user_beta$overall$CI.U
    REModels$overall$data.full$r.weights <- user_beta$overall$mdlWeights
  }
  
  dataList$overall$mdlWeights <- REModels$overall$data.full$r.weights

  tmp <- capture.output(rhoSensitivity$overall <- sensitivity(REModels$overall))
  
  if (measureType != "p")
    eggersRegTest$overall <- robu(
      formula      = Meas ~ 1 + sqrt(Var),
      data         = data.grp[data.grp$published,],
      studynum     = grp2,
      var.eff.size = Var,
      rho          = 0.8,
      small        = T,
      modelweights = "CORR",
      userweights  = REModels$overall$data.full$r.weights[data.grp$published]
    )
  
  return(
    list(
      REModels = REModels,
      rhoSensitivity = rhoSensitivity,
      eggersRegTest = eggersRegTest,
      data = data,
      dataList = dataList,
      notes = notes,
      measureType = measureType,
      zTrans = zTrans
    )
  )
}

Summary.RVE <- function(
  RVEModel,
  fileName,
  studyName = "Authors"
)
{
  output <- c()
  
  # Unpack
  REModels <- RVEModel$REModels
  rhoSensitivity <- RVEModel$rhoSensitivity
  eggersRegTest <- RVEModel$eggersRegTest
  notes <- RVEModel$notes
  zTrans <- RVEModel$zTrans
  measureType <- RVEModel$measureType
  dataList <- RVEModel$dataList
  
  output <- c(output, "Notes\n-----\n", paste(notes, collapse = ""), "\n")
  
  # Get first level groupings
  grping1_lvls <- names(REModels)
  
  for (grping in grping1_lvls)
  {
    # Print group heading
    output <- c(output, "\n\n========", grping, "========\n\n")
    
    # Print model statistics and tests
    Q <- t(REModels[[grping]]$data.full$e) %*% 
      diag(REModels[[grping]]$data.full$weights) %*% 
      REModels[[grping]]$data.full$e
    
    if (zTrans) {
      output <- c(output, "\tEffect-size (z):          \t", REModels[[grping]]$reg_table$b.r, "\n")
      output <- c(output, "\tEffect-size (eff-size):   \t", tanh(REModels[[grping]]$reg_table$b.r), "\n")
      output <- c(output, "\tStandard error (z):       \t", REModels[[grping]]$reg_table$SE, "\n")
    } else if (measureType == "p") {
      output <- c(output, "\tEffect-size (log-odds):   \t", REModels[[grping]]$reg_table$b.r, "\n")
      output <- c(output, "\tEffect-size (eff-size):   \t", 1/(1+exp(-REModels[[grping]]$reg_table$b.r)), "\n")
      output <- c(output, "\tStandard error (log-odds):\t", REModels[[grping]]$reg_table$SE, "\n")
    } else {
      output <- c(output, "\tEffect-size (eff-size):   \t", REModels[[grping]]$reg_table$b.r, "\n")
      output <- c(output, "\tStandard error (eff-size):\t", REModels[[grping]]$reg_table$SE, "\n")
    }
    output <- c(output, "\tk (ind. samples):              \t", length(unique(dataList[[grping]]$grp2)), "\n")
    output <- c(output, "\tn (participants):         \t", sum(aggregate(n ~ grp2, data = dataList[[grping]], max)$n), "\n")
    output <- c(output, "\tp-value:                  \t", REModels[[grping]]$reg_table$prob, "\n")
    output <- c(output, "\tdf:                       \t", REModels[[grping]]$reg_table$dfs, "\n")
    output <- c(output, "\tI.sq:                     \t", REModels[[grping]]$mod_info$I.2, "\n")
    output <- c(output, "\tQ-test:                   \t", Q, "\n")
    output <- c(output, "\tTau.sq:                   \t", REModels[[grping]]$mod_info$tau.sq, "\n\n")
    
    output <- c(output, "\tEgger's Regression Test for Bias\n\n")
    
    if (is.nan(eggersRegTest[[grping]]$reg_table$prob[2])) {
      output <- c(output, "\t\tEgger's test not defined (probably not enough data)\n\n")
    } else {
      output <- c(output, "\t\tp-value:         \t", eggersRegTest[[grping]]$reg_table$prob[2], 
                  if (eggersRegTest[[grping]]$reg_table$prob[2] < 0.05) "*" else "", "\n")
      output <- c(output, "\t\tdf:              \t", eggersRegTest[[grping]]$reg_table$dfs[2], "\n\n")
    }
    
    # Print sensitivity analysis
    output <- c(output, "\tIntra-study clustering sensitivity analysis:\n")
    output <- c(output, paste0("\t\t", capture.output(rhoSensitivity[[grping]]), sep = "\n"))
  }
  
  cat(output)
  capture.output(cat(output), file = fileName)
  
  flush.console()
  invisible()
}


forest.RVE <- function(
  RVEModel,
  studyNames,
  sampleNames,
  effSizeName,
  fileName,
  extraCol = NULL,
  extraColName = NULL,
  fullData = NULL
) 
{
  # Unpack
  REModels <- RVEModel$REModels
  data <- RVEModel$data
  zTrans <-RVEModel$zTrans
  data$study <- as.character(data[[studyNames]])
  data$sample <- data[[sampleNames]]
  grping1_lvls <- levels(data$grp1)
  
  measScale <- 1
  if (RVEModel$measureType == "p")
  {
    measScale <- 100
  }
  data$Meas <- data$Meas * measScale
  
  if (!is.null(extraCol)) {
    data$extra <- data[[extraCol]]
  }
  
  if (zTrans) {
    transf <- tanh
  } else {
    transf <- function(x) {x}
  }
  
  # Get weights; summary effects and standard errors; and n sums for plotting
  effSize <- c()
  standErr <- c()
  ci.lb <- c()
  ci.ub <- c()
  nSum <- c()
  data$weight <- NA
  for (grping in grping1_lvls)
  {
    indx <- data$grp1 == grping & !is.na(data$n)
    if (sum(indx) > 1) {
      data$weight[indx] <- REModels[[grping]]$data.full$r.weights / 
        sum(REModels[[grping]]$data.full$r.weights) * 100
      
      effSize[grping] <- REModels[[grping]]$reg_table$b.r
      standErr[grping] <- REModels[[grping]]$reg_table$SE
      ci.lb[grping] <- REModels[[grping]]$reg_table$CI.L
      ci.ub[grping] <- REModels[[grping]]$reg_table$CI.U
      nSum[grping] <- sum(aggregate(n ~ grp2, data = data[indx,], max)$n)
    } else {
      data$weight[indx] <- 100
      effSize[grping] <- NA
      standErr[grping] <- NA
      ci.lb[grping] <- NA
      ci.ub[grping] <- NA
      nSum[grping] <- NA
    }
  }
  
  # Sort data by level 1 and level 2 grouping and add an id
  data <- data %>% arrange(as.integer(grp1), study)
  data$ID <- 1:nrow(data)
  
  # Define headings for studies with multiple endpoints
  for (grping1 in grping1_lvls)
  {
    data.grp <- data[data$grp1 == grping1,]
    for (grping2 in unique(data.grp$grp2))
    {
      data.grp2 <- data.grp[data.grp$grp2 == grping2,]
      if (nrow(data.grp2) > 1) {
        minID <- min(data.grp2$ID)
        study <- data.grp2$study[1]
        indx <- data$grp1 == grping1 & data$grp2 == grping2
        data$study[indx] <- paste0("    ", data$sample[indx])
        data <- rbindlist(
          list(data, 
               data.frame(grp1 = grping1, grp2 = grping2, study = study, ID = minID - 0.5, 
                          n = max(data.grp2$n), weight = sum(data.grp2$weight), stringsAsFactors = F)),
          fill = T)
      }
    }
  }
  # Re-sort data
  data <- data %>% arrange(ID)
  
  # Define forest plot row formatting
  headingLocations <- c()
  summaryLocations <- c()
  sampleLocations <- c()
  
  for (grping in grping1_lvls)
  {
    if (length(headingLocations) == 0)
    {
      headingLocations <- 1
    }
    else
    {
      headingLocations <- c(headingLocations, max(c(sampleLocations, summaryLocations)) + 2)
    }
    sampleLocations <- c(
      sampleLocations, 
      max(headingLocations) + seq(1,sum(data$grp1 == grping),length = sum(data$grp1 == grping))
    )
    # if (sum(data$grp1 == grping) != 1)
      summaryLocations <- c(summaryLocations, max(headingLocations) + sum(data$grp1 == grping) + 1)
  }
  maxLocations <- max(summaryLocations, sampleLocations, headingLocations)
  headingLocations <- maxLocations - headingLocations + 1
  sampleLocations <- maxLocations - sampleLocations + 1
  summaryLocations <- maxLocations - summaryLocations + 1
  
  # Add indent to study names for formatting
  data$study <- paste0("    ", data$study)

  # Split out study headings and samples
  studyHeadingLocations <- sampleLocations[is.na(data$Meas)]
  sampleLocations <- sampleLocations[!is.na(data$Meas)]
  studyHeadings <- data$study[is.na(data$Meas)]
  studyHeadingsN <- data$n[is.na(data$Meas)]
  studyHeadingsWeight <- round(data$weight[is.na(data$Meas)], 2)
  data <- data[!is.na(data$Meas),]
  
  # Reset zero variances to zero and unity correlations to unity
  if (!zTrans) {
    data$Meas[data$Meas == 1] <- 1-1e-8
    data$Var[data$Var < 1e-8] <- 0
  }
  
  # Initialise forest plot with sample effect-sizes and CIs
  cex = 1.5
  op <- par(mar=c(4,4,0,2))
  ylim <- c(-1, maxLocations+3)
  efac <- 0.4
  digits <- if (RVEModel$measureType == "p") 1 else 2
  slab.col <- rep("black", nrow(data))
  slab.col[grep("        ", data$study)] <- "dimgray"
  anno.pos <- 2.5 * measScale
  if (is.null(extraCol)) {
    ilab = cbind(data$n, format(round(data$weight,2), trim = F))
    ilab.xpos <- c(-1.9, -1.2) * measScale + if (RVEModel$measureType == "p") measScale else 0
    ilab.pos <- 2
    xlim <- c(if (RVEModel$measureType == "p") -4 else -5, 2.5) * measScale
    colLabels <- c("N", "Weight", paste(effSizeName, "[95% CI]"))
  } else {
    ilab = cbind(data$n, format(round(data$weight,2), trim = F), data$extra)
    ilab.xpos <- c(c(-1.9, -1.2) * measScale + if (RVEModel$measureType == "p") measScale else 0, 3 * measScale)
    ilab.pos <- 2
    xlim <- c(if (RVEModel$measureType == "p") -4 else -5, 3) * measScale
    colLabels <- c("N", "Weight", extraColName, paste(effSizeName, "[95% CI]"))
  }
  
  
  tiff(filename = fileName,
       width=18, height = 3+25*maxLocations/75, units = "cm", res = 600,
       compression = "lzw", pointsize = 6)
  forest.col(
    x = data$Meas,
    vi = data$Var * measScale^2,
    xlim = xlim,
    alim = c(if (RVEModel$measureType == "p") 0 else -1, 1) * measScale,
    ylim = ylim,
    digits = digits,
    xlab = effSizeName,
    slab = data$study,
    ilab = ilab,
    ilab.xpos = ilab.xpos,
    ilab.pos = ilab.pos,
    transf = transf,
    rows = sampleLocations,
    efac = efac,
    cex = cex,
    slab.col = slab.col,
    anno.pos = anno.pos,
    refline = if (RVEModel$measureType == "p") c(0, 0.25, 0.5, 0.75, 1)*measScale else 0
  )
  
  # Add study headings (level 2 groupings)
  if (length(studyHeadings) > 0) {
    text(
      x = xlim[1],
      y = studyHeadingLocations,
      labels = studyHeadings,
      cex = cex,
      pos = 4
    )
    text(
      x = ilab.xpos[1],
      y = studyHeadingLocations,
      labels = studyHeadingsN,
      cex = cex,
      pos = 2
    )
    text(
      x = ilab.xpos[2],
      y = studyHeadingLocations,
      labels = format(studyHeadingsWeight),
      cex = cex,
      pos = 2
    )
  }
  
  # Add level 1 grouping headings
  text(
    x = xlim[1],
    y = headingLocations,
    labels = grping1_lvls,
    cex = cex,
    pos = 4,
    font = 4
  )
  
  # Add level 1 group summaries and overall summary
  grpIncl <- grping1_lvls %in% names(REModels)
  text(
    x = c(rep(ilab.xpos[1], sum(grpIncl)), rep(ilab.xpos[2], sum(grpIncl))),
    y = rep(summaryLocations[grpIncl], 2),
    labels = c(nSum[grpIncl], rep("100.00", sum(grpIncl))),
    cex = cex,
    pos = ilab.pos,
    font = 2
  )
  text(
    x = c(ilab.xpos[1], ilab.xpos[2]),
    y = ylim[1],
    labels = c(sum(aggregate(n ~ grp2, data = data, max)$n), "100.00"),
    cex = cex,
    pos = ilab.pos,
    font = 2
  )
  
  if(sum(!grpIncl) > 0)
    text(
      x = c(rep(ilab.xpos[1], sum(!grpIncl)), rep(ilab.xpos[2], sum(!grpIncl))),
      y = rep(summaryLocations[!grpIncl], 2),
      labels =  rep("N/A", sum(!grpIncl)*2),
      cex = cex,
      pos = ilab.pos,
      font = 2
    )
  
  addpoly.col(
    x = effSize[grpIncl],
    ci.lb = ci.lb[grpIncl],
    ci.ub = ci.ub[grpIncl],
    rows = summaryLocations[grpIncl],
    cex = cex,
    mlab = rep("    RE Model", sum(grpIncl)),
    transf = transf,
    efac = efac, 
    digits = digits,
    font = 2,
    anno.pos = anno.pos
  )
  addpoly.col(
    x = REModels$overall$reg_table$b.r,
    sei = REModels$overall$reg_table$SE,
    rows = ylim[1],
    cex = cex,
    mlab = "Overall RE Model",
    transf = transf,
    efac = efac, 
    digits = digits,
    font = 2,
    anno.pos = anno.pos
  )
  
  # Add column labels
  text(
    x = c(ilab.xpos, anno.pos),
    y = maxLocations + 2,
    labels = colLabels,
    cex = cex,
    pos = ilab.pos,
    font = 2
  )
  
  dev.off()
  
  invisible()
}

funnel.RVE <- function(
  RVEModel,
  fileName,
  trimnfill = T,
  plotReg = T,
  effSizeName = "Cohen's kappa"
)
{
  # Unpack
  data <- RVEModel$data
  zTrans <- RVEModel$zTrans
  REModels <- RVEModel$REModels
  eggersRegTest <- RVEModel$eggersRegTest
  grping1_lvls <- names(REModels)
  grpIncl <- grping1_lvls %in% levels(data$grp1)
  overallB <- REModels$overall$reg_table$b.r
  data$standErr <- sqrt(data$Var)
  maxStandErr <- max(data$standErr, na.rm = T)
  
  if (plotReg)
  {
    # Get Egger's regression line
    xEgg <- seq(-1,1,0.001)
    mdl <- REModels$overall
    egg <- eggersRegTest$overall
    intercept <- -egg$reg_table$b.r[1]/egg$reg_table$b.r[2]
    gradient <- 1/egg$reg_table$b.r[2]
    regLine <- xEgg * gradient + intercept
    print(paste0("Egger's regression p-value: ", egg$reg_table$prob[2]))
    
    if (zTrans)
      regLine <- tanh(regLine)
  }
  
  if (trimnfill)
  {
    tnf <- meta::trimfill.default(
      data$Meas[data$published], 
      data$standErr[data$published],
      left = T,
      sm = if (zTrans) "ZCOR" else "COR",
      ma.fixed = F,
      type = "L"
    )
    
    tnf_Meas <- tnf$TE[tnf$trimfill]
    tnf_standErr <- tnf$seTE[tnf$trimfill]
    tnf_study <- tnf$studlab[tnf$trimfill]
    
    if (zTrans)
    {
      tnf_standErr <- (tanh(tnf_Meas + tnf_standErr) - tanh(tnf_Meas - tnf_standErr)) / 2
      tnf_Meas <- tanh(tnf_Meas)
    }
  }
  
  if (zTrans)
  {
    data$standErr <- (tanh(data$Meas + sqrt(data$Var)) - tanh(data$Meas - sqrt(data$Var))) / 2
    data$Meas <- tanh(data$Meas)
    maxStandErr <- max(data$standErr)
    overallB <- tanh(overallB)
  }
  
  op <- par(mar = c(4,4,1,2))
  tiff(filename = fileName,
       width = 15, height = 12, units = "cm", res = 600,
       compression = "lzw", pointsize = 8)
  
  funnel(
    x = c(-99, -99),
    vi = c(99, 99),
    yaxis = "sei",
    xlim = c(-1, 1),
    ylim = c(0, maxStandErr),
    xlab = effSizeName,
    back = "white",
    shade = F,
    hlines = "white",
    pch = 1,
    cex = 1,
    refline = overallB
  )
  box()
  
  points(
    x = data$Meas[data$published],
    y = data$standErr[data$published],
    # pch = as.integer(data$grp1[data$published]),
    pch = 1,
    cex = 1
  )
  points(
    x = data$Meas[!data$published],
    y = data$standErr[!data$published],
    pch = 16,
    cex = 1
  )

  if (plotReg)
  {
    lines(xEgg[regLine >= 0], regLine[regLine >= 0], lty = 2)
  }
  
  if (trimnfill && length(tnf_Meas)>0)
  {
    points(
      x = tnf_Meas,
      y = tnf_standErr,
      pch = 17,
      cex = 1
    )
    
    crossTab <- table(data$grp1)
    data.grp <- data[data$grp1 %in% names(crossTab)[crossTab > 0], ]
    tnf.data.grp <- data.frame(
      grp2 = tnf_study,
      Meas = tnf_Meas,
      n = data$n[as.numeric(gsub("Filled: ", "", tnf_study))]
    )
    data.grp <- rbind.fill(data.grp, tnf.data.grp)
    
    data.grp <- escalc(
      ri        = Meas,
      ni        = n, 
      measure   = if (zTrans) "ZCOR" else "COR",
      data      = data.grp,
      var.names = c("Meas", "Var")
    )
    argList <- list(
      formula      = Meas ~ 1,
      data         = data.grp,
      studynum     = data.grp$grp2,
      var.eff.size = data.grp$Var,
      rho          = 0.8,
      small        = T,
      modelweights = "CORR"
    )
    adjMdl <- do.call(robu, argList)
    lines(c(adjMdl$reg_table$b.r, adjMdl$reg_table$b.r), c(0, 10), lty = 2)
  }
  
  legend(
    x = "topleft",
    legend = eval(parse(text = 
      paste0("expression('Published', 'Unpublished', 'Original estimate','Eggers regression line ('*italic('p')*' value = '*", as.character(round(egg$reg_table$prob[2], 3)), "*')')"))),
    pch = c(1, 16, NA, NA),
    lty = c(0, 0, 1, 2),
    col = "black",
    pt.cex = 1,
    bty = "n"
  )
  
  dev.off()
  
  invisible()
}

SensitivityAnalysis.RVE <- function(
  data,
  measureType = "r",
  measureCol = "Corr",
  sensitivityCol,
  n = "n",
  grping1, grping2,
  PA = NULL,
  zTrans = F,
  grpd = T)
{
  notes <- c()
  
  # Remove empty samples
  data$n <- data[[n]]
  data <- data[data$n > 0,]
  
  data$grp1 <- data[[grping1]]
  data$grp2 <- data[[grping2]]
  data$Meas <- data[[measureCol]]
  data$Sens <- data[[sensitivityCol]]
  IsBinary <- length(na.omit(unique(data$Sens))) <= 2
  MeasInit <- data$Meas
  SensInit <- data$Sens
  grping1_lvls <- levels(data$grp1)
  
  nName <- n
  
  # Get sample variances
  if (zTrans) {
    data <- escalc(
      ri        = Meas,
      ni        = n, 
      measure   = "ZCOR",
      data      = data,
      var.names = c("Meas", "Var")
    )
  } else if (measureType == "r") {
    data <- escalc(
      ri        = Meas,
      ni        = n, 
      measure   = "COR",
      data      = data,
      var.names = c("Meas", "Var")
    )
    data$Meas[data$Meas == 1] <- 1 - 1e-6
    data$Meas[data$Meas == 0] <- 1e-6
  } else if (measureType == "k") {
    pci <- (data[[PA]] - data$Meas) / (1 - data$Meas)
    data$Var <- data[[PA]] * (1 - data[[PA]]) / ((1 - pci) ^ 2 * data$n)
  } else if (measureType == "p") {
    data <- escalc(
      xi = Meas,
      ni = n,
      measure = "PLO",
      data = data,
      var.names = c("Meas", "Var")
    )
  } else {
    stop("Unknown 'measureType' specified.")
  }
  
  # Remove any undefined data (and take notes)
  if (any(is.na(data$Meas))) {
    notes <- c(notes, paste0(
      "\tThe following rows were removed because the effect-size was undefined:\n\n"))
    notes <- c(notes, paste0(
      "\t * ", data$grp1[is.na(data$Meas)], " - ", data$grp2[is.na(data$Meas)],
      " (", data$Meas[is.na(data$Meas)], ")\n"))
    notes <- c(notes, "\n")
    data <- data[!is.na(data$Meas),]
  }
  if (any(is.na(data$Sens))) {
    notes <- c(notes, paste0(
      "\tThe following rows were removed because the sensitivity variable was undefined:\n\n"))
    notes <- c(notes, paste0(
      "\t * ", data$grp1[is.na(data$Sens)], " - ", data$grp2[is.na(data$Sens)],
      " (", data$Sens[is.na(data$Sens)], ")\n"))
    notes <- c(notes, "\n")
    data <- data[!is.na(data$Sens),]
  }
  
  # Make sure none of variances are 0, otherwise analysis will fail
  data$Var[data$Var == 0] <- 1e-6
  
  # Remove first level groups with only 1 study
  for (grping in grping1_lvls)
    if (sum(data$grp1 == grping) < 1)
    {
      notes <- c(notes, paste0("\tRemoved '", grping, "' - ", sum(data$grp1 == grping), " studies.\n"))
      data <- data[data$grp1 != grping,]
      grping1_lvls <- grping1_lvls[grping1_lvls != grping]
    }
  
  # Build random-effects models using Robust Variance Estimation
  REModels <- list()
  dataList <- list()
  for (grping in grping1_lvls) {
    data.grp = data[data$grp1 == grping,]
    
    if ((!IsBinary || all(table(data.grp$Sens) >= 2)) && length(unique(data.grp$Sens)) > 1 && nrow(data.grp) > 2) {
      # Calculate random effects model
      argList <- list(
        formula      = Meas ~ 1 + Sens,
        data         = data.grp,
        studynum     = data.grp$grp2,
        var.eff.size = data.grp$Var,
        rho          = 0.8,
        small        = T,
        modelweights = "CORR"
      )
      REModels[[grping]] <- do.call(robu, argList)
      dataList[[grping]] <- data.grp
      dataList[[grping]]$mdlWeights <- REModels[[grping]]$data.full$r.weights
    } else {
      notes <- c(notes, paste0("\tSkipped '", grping, "' due to lack of data (we require at least 1 per group and 2 overall)\n"))
    }
  }

  if (grpd)
  {
    formula <- as.formula(Meas ~ 1 + Sens + grp1)
  }
  else
  {
    formula <- as.formula(Meas ~ 1 + Sens)
  }
  res <- try(
    {
      dataOA <- data
      if ((!IsBinary || all(table(dataOA$Sens) >= 2)) && length(unique(dataOA$Sens)) > 1 && nrow(dataOA) >= 2) {
        argList <- list(
          formula      = formula,
          data         = dataOA,
          studynum     = dataOA$grp2,
          var.eff.size = dataOA$Var,
          rho          = 0.8,
          small        = T,
          modelweights = "CORR"
        )
        REModels$overall <- do.call(robu, argList)
        dataList$overall <- dataOA
        dataList$overall$mdlWeights <- REModels$overall$data.full$r.weights
      }
    }
  )
  if (inherits(res, "try-error") || any(is.nan(REModels$overall$reg_table$prob))) {
    crossTab <- table(data$grp1)
    dataOA <- data[data$grp1 %in% names(crossTab)[crossTab > 1],]
    if ((!IsBinary || all(table(dataOA$Sens) >= 2)) && length(unique(dataOA$Sens)) > 1 && nrow(dataOA) >= 2 &&
        length(REModels) > 0) {
      argList <- list(
        formula      = formula,
        data         = dataOA,
        studynum     = dataOA$grp2,
        var.eff.size = dataOA$Var,
        rho          = 0.8,
        small        = T,
        modelweights = "CORR"
      )
      REModels$overall <- do.call(robu, argList)
      dataList$overall <- dataOA
      dataList$overall$mdlWeights <- REModels$overall$data.full$r.weights
    }
  }
  
  
  return(
    list(
      REModels = REModels,
      data = data,
      dataList = dataList,
      notes = notes,
      measureType = measureType,
      sensitivityCol = sensitivityCol,
      zTrans = zTrans
    )
  )
}

SensSummary.RVE <- function(
  RVEModel,
  fileName,
  isBinary = T,
  overallOnly = F,
  studyName = "Authors"
)
{
  output <- c()
  
  # Unpack
  REModels <- RVEModel$REModels
  notes <- RVEModel$notes
  zTrans <- RVEModel$zTrans
  measureType <- RVEModel$measureType
  sensCol <- RVEModel$sensitivityCol
  data <- RVEModel$data
  dataList <- RVEModel$dataList
  
  if (length(REModels) == 0) {
    output <- c(output, "Notes\n-----\n", paste(notes, collapse = ""), "\n\n")
    output <- c(output, "No models built - see notes\n")
    
    cat(output)
    capture.output(cat(output), file = fileName)
    
    return(invisible())
  } 
  
  output <- c(output, "\n\n========", sensCol, "========")
  output <- c(output, "\nNotes\n-----\n", paste(notes, collapse = ""), "\n")
  
  # Determine relvant levels
  if (isBinary) {
    factors <- c(0, 1)
  } else {
    factors <- c(
      min(data$Sens), 
      (min(data$Sens) + max(data$Sens))/2, 
      max(data$Sens)
    )
    output <- c(output, "\tLine of best fit estimations\n")
  }
  
  # Get first level groupings
  grping1_lvls <- names(REModels)
  if (overallOnly) {
    grping1_lvls <- "overall"
  }
  
  probVec <- c()
  for (grping in grping1_lvls)
  {
    # Print group heading
    output <- c(output, "\n", grping, "\n-------\n\n")
    
    regTable <- REModels[[grping]]$reg_table
    data.grp <- dataList[[grping]]
    if (is.character(data.grp$Sens)) {
      data.grp$Sens <- as.integer(as.factor(data.grp$Sens))
    }
    if (isBinary) {
      posGrp <- data.grp$Sens == 1
    } else {
      posGrp <- rep(T, nrow(data.grp))
    }
    
    output <- c(output, "\tk_total:      \t", length(unique(data.grp$grp2)), "\n")
    output <- c(output, "\tk_posGrp:     \t", length(unique(data.grp$grp2[posGrp])), "\n")
    output <- c(output, "\tn_total:      \t", sum(aggregate(n ~ grp2, data = data.grp, max)$n), "\n")
    output <- c(output, "\tn_posGrp:     \t", sum(aggregate(n ~ grp2, data = data.grp[posGrp,], max)$n), "\n")
    output <- c(output, "\tb_0:          \t", regTable$b.r[1], "\n")
    output <- c(output, "\tb_1:          \t", regTable$b.r[2], "\n")
    output <- c(output, "\tCI(b_0):      \t(", regTable$CI.L[1], ", ", regTable$CI.U[1], ")\n")
    output <- c(output, "\tCI(b_1):      \t(", regTable$CI.L[2], ", ", regTable$CI.U[2], ")\n")
    output <- c(output, "\tSE(b_0):      \t", regTable$SE[1], "\n")
    output <- c(output, "\tSE(b_1):      \t", regTable$SE[2], "\n")
    output <- c(output, "\tdf_0:         \t", regTable$dfs[1], "\n")
    output <- c(output, "\tdf_1:         \t", regTable$dfs[2], "\n")
    output <- c(output, "\tt_0:          \t", regTable$t[1], "\n")
    output <- c(output, "\tt_1:          \t", regTable$t[2], "\n")
    output <- c(output, "\tp_0:          \t", regTable$prob[1], "\n")
    output <- c(output, "\tp_1:          \t", regTable$prob[2], "\n")
    output <- c(output, "\tTauSq:        \t", REModels[[grping]]$mod_info$tau.sq, "\n")
    output <- c(output, "\tISq:          \t", REModels[[grping]]$mod_info$I.2, "\n\n")
    
    roundDig <- 3
    if (is.na(regTable$prob[1]))
    {
      probString0 <- probString1 <- "NA"
    }
    else
    {
      probString0 <- paste0(round(regTable$prob[1],roundDig), if (regTable$prob[1] < 0.01) "***" else if (regTable$prob[1] < 0.05) "**" else if (regTable$prob[1] < 0.1) "*" else "")
      probString1 <- paste0(round(regTable$prob[2],roundDig), if (regTable$prob[2] < 0.01) "***" else if (regTable$prob[2] < 0.05) "**" else if (regTable$prob[2] < 0.1) "*" else "")
    }
    if (grping == "overall")
    {
      desc <- c(
        sensCol,
        length(unique(data.grp$grp2)),
        length(unique(data.grp$grp2[posGrp])),
        sum(aggregate(n ~ grp2, data = data.grp, max)$n),
        sum(aggregate(n ~ grp2, data = data.grp[posGrp,], max)$n),
        paste0(round(regTable$b.r[1],roundDig), " [", round(regTable$CI.L[1],roundDig), ", ", round(regTable$CI.U[1],roundDig), "]"),
        paste0(round(regTable$b.r[2],roundDig), " [", round(regTable$CI.L[2],roundDig), ", ", round(regTable$CI.U[2],roundDig), "]"),
        # round(regTable$SE[1],roundDig),
        # round(regTable$SE[2],roundDig),
        round(regTable$dfs[1],roundDig),
        round(regTable$dfs[2],roundDig),
        round(regTable$t[1],roundDig),
        round(regTable$t[2],roundDig),
        probString0,
        probString1,
        round(REModels[[grping]]$mod_info$tau.sq,roundDig),
        round(REModels[[grping]]$mod_info$I.2,roundDig)
      )
    }
    probVec[[grping]] <- probString1
  }
    
  
  cat(output)
  capture.output(cat(output), file = fileName)
  
  flush.console()
  return(list(text = output, matrix = probVec, desc = desc))
}

logit <- function(p)
{
  log(p / (1 - p))
}

sigmoid <- function(x)
{
  1 / (1 + exp(-x))
}

# Wilson score method
#
# Wilson, E. B. ‘Probable inference, the law of succession, and statistical inference’,
# Journal of the AmericanStatistical Association
# 
p_CI <- function(n, p, z)
{
  list(
    (2 * n * p + z^2 + z * sqrt(z^2 + 4 * n * p * (1 - p))) / (2 * (n + z^2)),
    (2 * n * p + z^2 - z * sqrt(z^2 + 4 * n * p * (1 - p))) / (2 * (n + z^2))
  )
}