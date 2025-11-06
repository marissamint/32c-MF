


f.roughGate <- function(x, channel, thrd = 2.5, tinypeak.removal = 1 / 25, twin.factor = 0.98, alpha = 1 / 10) {
  
  if(class(x) == "flowFrame" | class(x) == "cytoframe"  ) {
    cellcount <- nrow(x)
    cat("perform roughGate ",identifier(x), "on", channel, "using", thrd, "as thrd",  "\n"  )
  }else if (class(x) == "CellPopulation") {
    cellcount <- x@cell.count
    cat("perform roughGate on","cell pop", "on", channel, "using", thrd, "as thrd",  "\n"  )
    
  }
  
  if (cellcount <= 20) {  
    cat("0 cell, return -Inf\n")
    cut = -Inf} else {
      pks <- getPeaks(x, channel, tinypeak.removal = tinypeak.removal)
      cat("Find peaks at:", pks$Peaks, "on", channel,  "\n")
      # print(pks$Peaks)
      if (length(pks$Peaks) == 1) { 
        if (pks$Peaks >= thrd) {
          cut <- deGate(x, channel, use.upper = F, upper = F, tinypeak.removal = tinypeak.removal, alpha = alpha)
          cat("Only one POS peak identified, find cut at", cut, "based on infelction point of POS Population\n")
          
        } else {
          cut <- deGate(x, channel, use.upper = T, upper = T, tinypeak.removal = tinypeak.removal, alpha = alpha) # ??
          cat("Only one NEG peak identified, find cut at", cut, "based on infelction point of Neg Population\n")
        }
      } else {
        if (max(pks$Peaks) < thrd) {
          cut <- deGate(x, channel, use.upper = T, upper = T, tinypeak.removal = tinypeak.removal, alpha = alpha)
          cat("No POS peaks, find cut at", cut, "based on infelction point of Neg Population\n")
        } else if (min(pks$Peaks) >= thrd) {
          cut <- deGate(x, channel, use.upper = T, upper = F, tinypeak.removal = tinypeak.removal, alpha = alpha)
          # if(cut > min(pks$Peaks)) {
          #   cut <- deGate(x, channel, use.upper = T, upper = F, tinypeak.removal = tinypeak.removal, alpha = alpha)
          # }
          
          cat("No NEG peaks, find cut at", cut, "based on infelction point of POS Population\n")
          
        } else {
          
          cut <- deGate(x, channel, tinypeak.removal = tinypeak.removal, twin.factor = twin.factor, all.cuts = T)
          cat("Find both POS and NEG peaks, find cuts at", cut,  ", choosing one cloest to threshold point\n")
          cut <- cut[which.min(abs(cut - thrd))]
        }
      }
    }
  
  return(cut)
}





f.remove_pop_outliers <- function(obj, channels, axis.x = T, axis.y = T, quantiles = c(0.0001, 0.9999)) {
  if (axis.x == T) {
    min.x <- deGate(obj, channels[1], use.percentile = T, percentile = quantiles[1])
    max.x <- deGate(obj, channels[1], use.percentile = T, percentile = quantiles[2])
    #obj <- flowDensity(obj, channels, position = c(T, NA), gates = c(min.x, NA))
    #obj <- flowDensity(obj, channels, position = c(F, NA), gates = c(max.x, NA))
  }
  
  if (axis.y == T) {
    min.y <- deGate(obj, channels[2], use.percentile = T, percentile = quantiles[1])
    max.y <- deGate(obj, channels[2], use.percentile = T, percentile = quantiles[2])
    #obj <- flowDensity(obj, channels, position = c(NA, T), gates = c(NA, min.y))
    #obj <- flowDensity(obj, channels, position = c(NA, F), gates = c(NA, max.y))
  }
  
  obj <- flowDensity(obj, channels, position = c(T, T), gates = c(min.x, min.y))
  obj <- flowDensity(obj, channels, position = c(F, F), gates = c(max.x, max.y))
  
  
  
}









# x <-fa
# ur =T
# lr =T

f.GatePrettify <- function(x, ll = F, lr = F, ul = F, ur = F) {
  
  if(max(x) == -Inf | is.na(max(x))) {
    x
  } else {
    g <- x[1:nrow(x)-1, ] # remove last row
    #g<- g[order(g[,1], decreasing = T), ]
    
    #max.x.index <- which.max(g[, 1])
    #reorder to make the biggest x on the top
    g <- rbind(g[which.max(g[, 1]): nrow(g),], g[1:which.max(g[, 1])-1,]) 
    
    
    max.x.index <- 1
    min.x.index <- which.min(g[, 1])
    min.y.index <- which.min(g[, 2])
    max.y.index <- which.max(g[, 2])
    
    g <- rbind(g,g[1,])
    if(max.y.index ==1) {max.y.index <- nrow(g)} # fixed the bug
    
    g.lr <- g[1:min.y.index,, drop= F] 
    g.ll <- g[min.y.index : min.x.index,, drop= F] 
    g.ul <- g[min.x.index : max.y.index,,drop= F]   # upper left
    g.ur <- g[max.y.index : nrow(g),, drop= F]  # upper right
    
    
    # plot(g.lr)
    # plot(g.ll)
    # plot( g.ul)
    # plot(g.ur)
    
    if(lr ==T ) { # lower right corner
      g.lr <- rbind(g.lr[1, ], 
                    c(g[max.x.index,1], g.lr[nrow(g.lr), 2]), 
                    g.lr[nrow(g.lr), ])
    }
    
    if(ll ==T ) {
      g.ll <- rbind(g.ll[1, ], 
                    c(g[min.x.index,1], g.ll[1, 2]), 
                    g.ll[nrow(g.ll), ])
    }
    
    if(ul == T) {
      g.ul <- rbind(g.ul[1, ], 
                    c(g[min.x.index,1], g.ul[nrow(g.ul), 2]), 
                    g.ul[nrow(g.ul), ])
    }
    
    if(ur == T) {
      g.ur <- rbind(g.ur[1, ], 
                    c(g[max.x.index,1], g.ur[1, 2]), 
                    g.ur[nrow(g.ur), ])
    }
    
    g2 <- rbind(g.lr, g.ll, g.ul, g.ur)
    g2 <- g2[!duplicated(g2),]
    g2 <- rbind(g2,g2[1,])
    
    return(g2)
    
  }
  
  
}


# this is an uility function for the general gating
#' @param pops a vector of population names
#' @param obj gatingset
f.extract_cut_from_gs <- function(obj, pops, channel) {
  lapply(obj, function(x) {
    gates <- unlist(lapply(pops, function(g) {
      g <- gh_pop_get_gate(x, g)
      g@min[channel]
    }))
    
    min(gates[!is.infinite(gates)])
  })
}


f.percent_annotation <- function(numerator, denominator, gh,  decimal = 1) {
  pct <- NULL
  for (n in 1:length(numerator)) {
    if( !(is.na(numerator[n]) | numerator[n] == "NA")  ) {
      v <- paste0(round(gh_pop_get_count(gh, numerator[n], xml = F) / gh_pop_get_count(gh, denominator, xml = F) * 100, decimal))
    } else {v <- NA}
    pct <- c(pct, v)
  }
  
  positions <- c("topleft", "topright", "bottomleft", "bottomright")
  adj <- list(c(0.8, -0.2), c(-0.1, -0.2), c(0.8, 1.1), c(-0.1, 1.1))
  ann.l <- list()
  for (p in 1:length(pct)) {
    ann.l <- c(ann.l, list(
      c(pct[p], 
        positions[p],
        adj[p]
        )
    )
    )
  }
  
  
  return(ann.l)
}


f.add_pct_to_plot <- function(percentage_list, background_col = "white", background_alpha = 0) {
  
  bgcolor <- alpha(background_col, background_alpha)

  lapply(percentage_list, function(x) {
    if(!is.na(x[[1]])) {
      
    legend(x[[2]],#pct_antn[[i]][2], 
           x[[1]], 
           #col = "red", 
           box.col = bgcolor, 
           bg = bgcolor, 
           #xjust = 0, yjust = 0, 
           adj = x[[3]],
           text.font = 1, cex = 2, 
           text.col = alpha("black", 0.7)
    )
    }
    
  })
  
  
}







