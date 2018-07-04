forest.col <- function (x, vi, sei, ci.lb, ci.ub, annotate = TRUE, showweights = FALSE, 
    xlim, alim, clim, ylim, at, steps = 5, level = 95, refline = 0, 
    digits = 2L, width, xlab, slab, ilab, ilab.xpos, ilab.pos, 
    subset, transf, atransf, targs, rows, efac = 1, pch = 15, 
    psize, col, slab.col, lty, cex, cex.lab, cex.axis, anno.pos, ...) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(transf)) 
        transf <- FALSE
    if (missing(atransf)) 
        atransf <- FALSE
    transf.char <- deparse(substitute(transf))
    atransf.char <- deparse(substitute(atransf))
    if (is.function(transf) && is.function(atransf)) 
        stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")
    if (missing(targs)) 
        targs <- NULL
    if (missing(at)) 
        at <- NULL
    if (missing(ilab)) 
        ilab <- NULL
    if (missing(ilab.xpos)) 
        ilab.xpos <- NULL
    if (missing(ilab.pos)) 
        ilab.pos <- NULL
    if (missing(subset)) 
        subset <- NULL
    if (missing(psize)) 
        psize <- NULL
    if (missing(col)) 
        col <- NULL
    if (missing(slab.col)) 
        slab.col <- NULL
    if (missing(cex)) 
        cex <- NULL
    if (missing(cex.lab)) 
        cex.lab <- NULL
    if (missing(cex.axis)) 
        cex.axis <- NULL
    if (missing(lty)) {
        lty <- c("solid", "solid")
    }
    else {
        if (length(lty) == 1L) 
            lty <- c(lty, "solid")
    }
    if (length(efac) == 1L) 
        efac <- rep(efac, 2)
    if (length(digits) == 1L) 
        digits <- c(digits, digits)
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    yi <- x
    if (is.null(attr(yi, "measure"))) {
        measure <- "GEN"
    }
    else {
        measure <- attr(yi, "measure")
    }
    if (hasArg(ci.lb) && hasArg(ci.ub)) {
        if (length(ci.lb) != length(ci.ub)) 
            stop("Length of 'ci.lb' and 'ci.ub' do not match.")
        if (missing(vi) && missing(sei)) {
            vi <- ((ci.ub - ci.lb)/(2 * qnorm(alpha/2, lower.tail = FALSE)))^2
        }
        else {
            if (missing(vi)) 
                vi <- sei^2
        }
        if (length(ci.lb) != length(vi)) 
            stop("Length of 'vi' (or 'sei') does not match length of ('ci.lb', 'ci.ub') pairs.")
    }
    else {
        if (missing(vi)) {
            if (missing(sei)) {
                stop("Must specify either 'vi', 'sei', or ('ci.lb', 'ci.ub') pairs.")
            }
            else {
                vi <- sei^2
                ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
                  sei
                ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
                  sei
            }
        }
        else {
            ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
                sqrt(vi)
            ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
                sqrt(vi)
        }
    }
    if (length(yi) != length(vi)) 
        stop("Length of 'yi' does not match the length of 'vi', 'sei', or the ('ci.lb', 'ci.ub') pairs.")
    k <- length(yi)
    if (missing(slab)) {
        if (!is.null(attr(yi, "slab"))) {
            slab <- attr(yi, "slab")
        }
        else {
            slab <- paste("Study", seq_len(k))
        }
    }
    else {
        if (length(slab) == 1 && is.na(slab)) 
            slab <- rep("", k)
    }
    if (length(yi) != length(slab)) 
        stop("Number of outcomes does not correspond to the length of the 'slab' argument.")
    if (is.null(dim(ilab))) 
        ilab <- cbind(ilab)
    if (length(pch) == 1L) 
        pch <- rep(pch, k)
    if (length(pch) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the 'pch' argument.")
    if (!is.null(psize)) {
        if (length(psize) == 1L) 
            psize <- rep(psize, k)
        if (length(psize) != length(yi)) 
            stop("Number of outcomes does not correspond to the length of the 'psize' argument.")
    }
    if (!is.null(col)) {
        if (length(col) == 1L) 
            col <- rep(col, k)
        if (length(col) != length(yi)) 
            stop("Number of outcomes does not correspond to the length of the 'col' argument.")
    }
    else {
        col <- rep("black", k)
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        ci.lb <- ci.lb[subset]
        ci.ub <- ci.ub[subset]
        slab <- slab[subset]
		slab.col <- slab.col[subset]
        ilab <- ilab[subset, , drop = FALSE]
        pch <- pch[subset]
        psize <- psize[subset]
        col <- col[subset]
    }
    k <- length(yi)
    if (missing(rows)) {
        rows <- k:1
    }
    else {
        if (length(rows) == 1L) 
            rows <- rows:(rows - k + 1)
    }
    if (length(rows) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the 'rows' argument.")
    yi <- yi[k:1]
    vi <- vi[k:1]
    ci.lb <- ci.lb[k:1]
    ci.ub <- ci.ub[k:1]
    slab <- slab[k:1]
	slab.col <- slab.col[k:1]
    ilab <- ilab[k:1, , drop = FALSE]
    pch <- pch[k:1]
    psize <- psize[k:1]
    col <- col[k:1]
    rows <- rows[k:1]
    yivi.na <- is.na(yi) | is.na(vi)
    if (any(yivi.na)) {
        not.na <- !yivi.na
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            ci.lb <- ci.lb[not.na]
            ci.ub <- ci.ub[not.na]
            slab <- slab[not.na]
			slab.col <- slab.col[not.na]
            ilab <- ilab[not.na, , drop = FALSE]
            pch <- pch[not.na]
            psize <- psize[not.na]
            col <- col[not.na]
            rows.new <- rows
            rows.na <- rows[!not.na]
            for (j in seq_len(length(rows.na))) {
                rows.new[rows >= rows.na[j]] <- rows.new[rows >= 
                  rows.na[j]] - 1
            }
            rows <- rows.new[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    if (is.function(transf)) {
        if (is.null(targs)) {
            yi <- sapply(yi, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
        }
        else {
            yi <- sapply(yi, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
        }
    }
    tmp <- .psort(ci.lb, ci.ub)
    ci.lb <- tmp[, 1]
    ci.ub <- tmp[, 2]
    if (!missing(clim)) {
        clim <- sort(clim)
        if (length(clim) != 2L) 
            stop("Argument 'clim' must be of length 2.")
        ci.lb[ci.lb < clim[1]] <- clim[1]
        ci.ub[ci.ub > clim[2]] <- clim[2]
    }
    if (showweights) {
        weights <- 1/vi
        weights <- 100 * weights/sum(weights, na.rm = TRUE)
    }
    if (is.null(psize)) {
        if (any(vi <= 0, na.rm = TRUE)) {
            psize <- rep(1, k)
        }
        else {
            wi <- 1/sqrt(vi)
            psize <- wi/sum(wi, na.rm = TRUE)
            psize <- (psize - min(psize, na.rm = TRUE))/(max(psize, 
                na.rm = TRUE) - min(psize, na.rm = TRUE))
            psize <- (psize * 1) + 0.5
            if (all(is.na(psize))) 
                psize <- rep(1, k)
        }
    }
    rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
    if (annotate) {
        if (showweights) {
            plot.multp.l <- 2
            plot.multp.r <- 2
        }
        else {
            plot.multp.l <- 1.2
            plot.multp.r <- 1.2
        }
    }
    else {
        plot.multp.l <- 1.2
        plot.multp.r <- 0.4
    }
    if (missing(xlim)) {
        xlim <- c(min(ci.lb, na.rm = TRUE) - rng * plot.multp.l, 
            max(ci.ub, na.rm = TRUE) + rng * plot.multp.r)
        xlim <- round(xlim, digits[2])
    }
    alim.spec <- TRUE
    if (missing(alim)) {
        if (is.null(at)) {
            alim <- range(pretty(x = c(min(ci.lb, na.rm = TRUE), 
                max(ci.ub, na.rm = TRUE)), n = steps - 1))
            alim.spec <- FALSE
        }
        else {
            alim <- range(at)
        }
    }
    alim <- sort(alim)
    xlim <- sort(xlim)
    if (xlim[1] > min(yi, na.rm = TRUE)) {
        xlim[1] <- min(yi, na.rm = TRUE)
    }
    if (xlim[2] < max(yi, na.rm = TRUE)) {
        xlim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] < xlim[1]) {
        xlim[1] <- alim[1]
    }
    if (alim[2] > xlim[2]) {
        xlim[2] <- alim[2]
    }
    if (missing(ylim)) {
        ylim <- c(0.5, k + 3)
    }
    else {
        ylim <- sort(ylim)
    }
    if (is.null(at)) {
        if (alim.spec) {
            at <- seq(from = alim[1], to = alim[2], length.out = steps)
        }
        else {
            at <- pretty(x = c(min(ci.lb, na.rm = TRUE), max(ci.ub, 
                na.rm = TRUE)), n = steps - 1)
        }
    }
    else {
        at[at < alim[1]] <- alim[1]
        at[at > alim[2]] <- alim[2]
        at <- unique(at)
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits[2], 
                format = "f", drop0trailing = ifelse(class(digits) == 
                  "integer", TRUE, FALSE))
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits[2], format = "f", drop0trailing = ifelse(class(digits) == 
                  "integer", TRUE, FALSE))
        }
    }
    else {
        at.lab <- formatC(at.lab, digits = digits[2], format = "f", 
            drop0trailing = ifelse(class(digits) == "integer", 
                TRUE, FALSE))
    }
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, 3, 1, 1)
    par.mar.adj[par.mar.adj < 0] <- 0
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", 
        yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", col = "black", 
        ...)
    abline(h = ylim[2] - 2, lty = lty[2], col = "black", ...)
    if (is.numeric(refline))
    {
      segments(refline, ylim[1] - 5, refline, ylim[2] - 2, 
               lty = "dotted", col = "black", ...)
      segments(refline[ceiling(length(refline)/2)], ylim[1] - 5,
               refline[ceiling(length(refline)/2)], ylim[2] - 2, 
               col = "black", ...)
    }
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    if (is.null(cex)) {
        lheight <- strheight("O")
        cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * 
            k * lheight), 1)
    }
    if (is.null(cex)) {
        cex <- par("cex") * cex.adj
    }
    else {
        if (is.null(cex.lab)) 
            cex.lab <- cex
        if (is.null(cex.axis)) 
            cex.axis <- cex
    }
    if (is.null(cex.lab)) 
        cex.lab <- par("cex") * cex.adj
    if (is.null(cex.axis)) 
        cex.axis <- par("cex") * cex.adj
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis, 
        col = "black", ...)
    if (missing(xlab)) 
        xlab <- .setlab(measure, transf.char, atransf.char, gentype = 1)
    mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
        line = par("mgp")[1] - 0.5, cex = cex.lab, col = "black", 
        ...)
    for (i in seq_len(k)) {
        if (is.na(yi[i]) || is.na(ci.lb[i]) || is.na(ci.ub[i])) 
            next
        if (ci.lb[i] >= alim[2]) {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac[2], rows[i] - 
                  (height/150) * cex * efac[2], rows[i]), col = col[i], 
                border = col[i], ...)
            next
        }
        if (ci.ub[i] <= alim[1]) {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac[2], rows[i] - 
                  (height/150) * cex * efac[2], rows[i]), col = col[i], 
                border = col[i], ...)
            next
        }
        segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], 
            alim[2]), rows[i], lty = lty[1], col = col[i], ...)
        if (ci.lb[i] >= alim[1]) {
            segments(ci.lb[i], rows[i] - (height/150) * cex * 
                efac[1], ci.lb[i], rows[i] + (height/150) * cex * 
                efac[1], col = col[i], ...)
        }
        else {
          segments(alim[1], rows[i] - (height/150) * cex * 
                     efac[1], alim[1], rows[i] + (height/150) * cex * 
                     efac[1], col = col[i], ...)
            # polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
            #     (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
            #     (xlim[2] - xlim[1]), alim[1]), y = c(rows[i], 
            #     rows[i] + (height/150) * cex * efac[2], rows[i] - 
            #       (height/150) * cex * efac[2], rows[i]), col = col[i], 
            #     border = col[i], ...)
        }
        if (ci.ub[i] <= alim[2]) {
            segments(ci.ub[i], rows[i] - (height/150) * cex * 
                efac[1], ci.ub[i], rows[i] + (height/150) * cex * 
                efac[1], col = col[i], ...)
        }
        else {
          segments(alim[2], rows[i] - (height/150) * cex * 
                     efac[1], alim[2], rows[i] + (height/150) * cex * 
                     efac[1], col = col[i], ...)
            # polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
            #     (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
            #     (xlim[2] - xlim[1]), alim[2]), y = c(rows[i], 
            #     rows[i] + (height/150) * cex * efac[2], rows[i] - 
            #       (height/150) * cex * efac[2], rows[i]), col = col[i], 
            #     border = col[i], ...)
        }
    }
    text(xlim[1], rows, slab, pos = 4, cex = cex, col = slab.col, 
        ...)
    if (!is.null(ilab)) {
        if (is.null(ilab.xpos)) 
            stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'.")
        if (length(ilab.xpos) != ncol(ilab)) 
            stop(paste0("Number of 'ilab' columns (", ncol(ilab), 
                ") does not match length of 'ilab.xpos' argument (", 
                length(ilab.xpos), ")."))
        if (!is.null(ilab.pos) && length(ilab.pos) == 1) 
            ilab.pos <- rep(ilab.pos, ncol(ilab))
        for (l in seq_len(ncol(ilab))) {
          if (l <= 2)
            text(ilab.xpos[l], rows, ilab[, l], pos = ilab.pos[l], 
                cex = cex, col = slab.col, ...)
          else
            text(ilab.xpos[l], rows, ilab[, l], pos = ilab.pos[l], 
                 cex = cex, ...)
        }
    }
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                annotext <- cbind(sapply(yi, atransf), sapply(ci.lb, 
                  atransf), sapply(ci.ub, atransf))
            }
            else {
                annotext <- cbind(sapply(yi, atransf, targs), 
                  sapply(ci.lb, atransf, targs), sapply(ci.ub, 
                    atransf, targs))
            }
            tmp <- .psort(annotext[, 2:3])
            annotext[, 2:3] <- tmp
        }
        else {
            annotext <- cbind(yi, ci.lb, ci.ub)
        }
        if (showweights) 
            annotext <- cbind(weights, annotext)
        annotext <- formatC(annotext, format = "f", digits = digits[1])
        if (missing(width)) {
            width <- apply(annotext, 2, function(x) max(nchar(x)))
        }
        else {
            if (length(width) == 1L) 
                width <- rep(width, ncol(annotext))
        }
        for (j in 1:ncol(annotext)) {
            annotext[, j] <- formatC(annotext[, j], width = width[j])
        }
        if (showweights) {
            annotext <- cbind(annotext[, 1], "%   ", annotext[, 
                2], " [", annotext[, 3], ", ", annotext[, 4], 
                "]")
        }
        else {
            annotext <- cbind(annotext[, 1], " [", annotext[, 
                2], ", ", annotext[, 3], "]")
        }
        annotext <- apply(annotext, 1, paste, collapse = "")
        text(x = anno.pos, rows, labels = annotext, pos = 2, cex = cex, 
            col = col, ...)
    }
    for (i in seq_len(k)) {
        if (is.na(yi[i])) 
            next
        if (yi[i] >= alim[1] && yi[i] <= alim[2]) 
            points(yi[i], rows[i], pch = pch[i], cex = cex * 
                psize[i], col = col[i], ...)
    }
    res <- list(xlim = par("usr")[1:2], alim = alim, at = at, 
        ylim = ylim, rows = rows, cex = cex, cex.lab = cex.lab, 
        cex.axis = cex.axis)
    invisible(res)
}


.psort <- function(x,y) {
  
  ### t(apply(xy, 1, sort)) would be okay, but problematic if there are NAs;
  ### either they are removed completely (na.last=NA) or they are always put
  ### first/last (na.last=FALSE/TRUE); but we just want to leave the NAs in
  ### their position!
  
  if (is.null(x) || length(x) == 0) ### need to catch this
    return(NULL)
  
  if (missing(y)) {
    if (is.matrix(x)) {
      xy <- x
    } else {
      xy <- rbind(x) ### in case x is just a vector
    }
  } else {
    xy <- cbind(x,y)
  }
  
  n <- nrow(xy)
  
  for (i in 1:n) {
    if (anyNA(xy[i,]))
      next
    xy[i,] <- sort(xy[i,])
  }
  
  colnames(xy) <- NULL
  
  return(xy)
  
}

addpoly.col <- function (x, vi, sei, ci.lb, ci.ub, rows = -1, level = 95, annotate = TRUE, 
                         digits = 2, width, mlab, transf, atransf, targs, efac = 1, 
                         col, border, cex, anno.pos, ...) 
{
  na.act <- getOption("na.action")
  if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
                            "na.pass"))) 
    stop("Unknown 'na.action' specified under options().")
  if (missing(transf)) 
    transf <- FALSE
  if (missing(atransf)) 
    atransf <- FALSE
  if (is.function(transf) && is.function(atransf)) 
    stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")
  if (missing(targs)) 
    targs <- NULL
  if (missing(mlab)) 
    mlab <- NULL
  if (missing(cex)) 
    cex <- NULL
  if (missing(col)) 
    col <- "black"
  if (missing(border)) 
    border <- "black"
  alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
  yi <- x
  if (hasArg(ci.lb) && hasArg(ci.ub)) {
    if (length(ci.lb) != length(ci.ub)) 
      stop("Length of 'ci.lb' and 'ci.ub' do not match.")
    if (missing(vi) && missing(sei)) {
      vi <- ((ci.ub - ci.lb)/(2 * qnorm(alpha/2, lower.tail = FALSE)))^2
    }
    else {
      if (missing(vi)) 
        vi <- sei^2
    }
    if (length(ci.lb) != length(vi)) 
      stop("Length of 'vi' (or 'sei') does not match length of ('ci.lb', 'ci.ub') pairs.")
  }
  else {
    if (missing(vi)) {
      if (missing(sei)) {
        stop("Must specify either 'vi', 'sei', or ('ci.lb', 'ci.ub') pairs.")
      }
      else {
        vi <- sei^2
        ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
          sei
        ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
          sei
      }
    }
    else {
      ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
        sqrt(vi)
      ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
        sqrt(vi)
    }
  }
  k <- length(yi)
  if (is.null(rows)) {
    rows <- -1:(-k)
  }
  else {
    if (length(rows) == 1L) 
      rows <- rows:(rows - k + 1)
  }
  if (length(rows) != length(yi)) 
    stop("Number of outcomes does not correspond to the length of the 'rows' argument.")
  yivi.na <- is.na(yi) | is.na(vi)
  if (any(yivi.na)) {
    not.na <- !yivi.na
    if (na.act == "na.omit") {
      yi <- yi[not.na]
      vi <- vi[not.na]
      ci.lb <- ci.lb[not.na]
      ci.ub <- ci.ub[not.na]
      mlab <- mlab[not.na]
      rows.new <- rows
      rows.na <- rows[!not.na]
      for (j in seq_len(length(rows.na))) {
        rows.new[rows <= rows.na[j]] <- rows.new[rows <= 
                                                   rows.na[j]] + 1
      }
      rows <- rows.new[not.na]
    }
    if (na.act == "na.fail") 
      stop("Missing values in results.")
  }
  k <- length(yi)
  if (is.function(transf)) {
    if (is.null(targs)) {
      yi <- sapply(yi, transf)
      ci.lb <- sapply(ci.lb, transf)
      ci.ub <- sapply(ci.ub, transf)
    }
    else {
      yi <- sapply(yi, transf, targs)
      ci.lb <- sapply(ci.lb, transf, targs)
      ci.ub <- sapply(ci.ub, transf, targs)
    }
  }
  tmp <- .psort(ci.lb, ci.ub)
  ci.lb <- tmp[, 1]
  ci.ub <- tmp[, 2]
  par.usr <- par("usr")
  height <- par.usr[4] - par.usr[3]
  cex.adj <- min(1, 20/height)
  xlim <- par.usr[1:2]
  if (is.null(cex)) 
    cex <- par("cex") * cex.adj
  if (annotate) {
    if (is.function(atransf)) {
      if (is.null(targs)) {
        annotext <- cbind(sapply(yi, atransf), sapply(ci.lb, 
                                                      atransf), sapply(ci.ub, atransf))
      }
      else {
        annotext <- cbind(sapply(yi, atransf, targs), 
                          sapply(ci.lb, atransf, targs), sapply(ci.ub, 
                                                                atransf, targs))
      }
      tmp <- .psort(annotext[, 2:3])
      annotext[, 2:3] <- tmp
    }
    else {
      annotext <- cbind(yi, ci.lb, ci.ub)
    }
    annotext <- formatC(annotext, format = "f", digits = digits)
    if (missing(width) || is.null(width)) {
      width <- apply(annotext, 2, function(x) max(nchar(x)))
    }
    else {
      if (length(width) == 1L) 
        width <- rep(width, ncol(annotext))
    }
    for (j in 1:ncol(annotext)) {
      annotext[, j] <- formatC(annotext[, j], width = width[j])
    }
    annotext <- cbind(annotext[, 1], " [", annotext[, 2], 
                      ", ", annotext[, 3], "]")
    annotext <- apply(annotext, 1, paste, collapse = "")
    text(x = anno.pos, rows, labels = annotext, pos = 2, cex = cex, 
         ...)
  }
  for (i in seq_len(k)) {
    polygon(x = c(ci.lb[i], yi[i], ci.ub[i], yi[i]), y = c(rows[i], 
                                                           rows[i] + (height/100) * cex * efac, rows[i], rows[i] - 
                                                             (height/100) * cex * efac), col = col, border = border, 
            ...)
    if (!is.null(mlab)) 
      text(xlim[1], rows[i], mlab[i], pos = 4, cex = cex, 
           ...)
  }
}