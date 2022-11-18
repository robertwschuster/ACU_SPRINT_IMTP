#-----------------------------------------------------------------------------------------
#
# Functions for IMTP analysis Shinny app
# Robert Schuster (ACU SPRINT)
# August 2022
#
#-----------------------------------------------------------------------------------------

# Basic functions ------------------------------------------------------------------------
# Find peaks
# https://github.com/stas-g/findPeaks
# a 'peak' is defined as a local maxima with m points either side of it being smaller than it. 
# hence, the bigger the parameter m, the more stringent is the peak finding procedure
find_peaks <- function(x, m = 3) {
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

# moving standard deviation
movstd <- function(vec, width)      
  return(c(rep(NA,width-1), sapply(seq_along(vec[width:length(vec)]),      
                                   function(i) sd(vec[i:(i+width-1)]))))

# find rep around a local maximum
find_rep <- function(fz, fmaxi, width = 500) {
  cutoff <- movstd(fz[1:fmaxi], width)
  coi <- which(cutoff == min(cutoff, na.rm = T)) + 1
  cutoff <- min(cutoff, na.rm = T)
  
  baseline <- mean(fz[(coi - width):(coi - 1)])
  
  out <- list(start = max(which(fz[1:fmaxi] <= (baseline + cutoff*5))),
              end = min(which(fz[fmaxi:length(fz)] <= (baseline + cutoff*5))) + fmaxi)
  return(out)
}


# Load data ------------------------------------------------------------------------------
importTrial <- function(filepath,filename) {
  cols <- c("Time","Z.Left","Z.Right","Left","Right")
  fn <- basename(filename)
  # determine length of header
  skip <- read.csv(filepath, nrows = 20, header = F, blank.lines.skip = F, 
                   col.names = paste0("V",seq_len(max(count.fields(filepath, sep = ',')))), fill = T)
  if (any(skip == "Time", na.rm = T)) {
    skip <- which(skip == 'Time')-1
    
    # load data (without header) and keep only Time, Left and Right columns
    df <- read.csv(filepath, skip = skip)
    if (any(names(df) %in% cols, na.rm = T)) {
      df <- df[,which(names(df) %in% cols)] # keep only left, right and time columns
      colnames(df) <- gsub('Z.','',colnames(df)) # remove 'Z." from column names
      
      df$Total <- df$Left + df$Right
      
      # load header
      header <- read.csv(filepath, nrows = skip-1, header = F)
      # find weight for force normalisation
      if (any(header == "Weight", na.rm = T)) {
        bodymass <- as.numeric(na.omit(header[header == 'Weight',2])) * 9.81
      } else {
        bodymass <- median(df$Total[df$Total >= 300])
      }
      # find frequency
      freq <- as.numeric(na.omit(header[header == 'Frequency',2]))
      if (freq < 1000) {
        message(paste(f,": Sampling frequency =",freq,"\nFrequencies > 1000 Hz are recommended"))
      }
      
      data <- list("df" = df, "bodymass" = bodymass, "freq" = freq, "fn" = fn)
      
      return(data)
      
    } else {
      # print warning but continue execution with other files
      warning(paste("The file ", fn, " does not contain the required data and was not processed"))
    }
  } else {
    # print warning but continue execution with other files
    warning(paste("The file", fn, "does not contain the required data and was not processed"))
  }
}


# Determine number of repetitions --------------------------------------------------------
nReps <- function(data) {
  # find peaks (at least 2s apart from each other)
  d <- data$freq*2 # frequency * 2s
  pks <- find_peaks(data$df$Total, m = d)
  # determine a cutoff under which all other peaks cannot be considered IMTP trials
  cutoff <- (max(data$df$Total[pks])-data$bodymass)/2 + data$bodymass # half of BM normalised max force
  # cutoff <- max(data$df$Total[pks])-250 # within 250 N of max value
  # only keep peaks above cutoff and more than 5s apart
  pks <- pks[which(data$df$Total[pks] > cutoff)]
  d <- data$freq*5 # frequency * 5s
  if (any(diff(pks) < d)) {
    pks <- pks[-(which(diff(pks) < d)+1)]
  }
  
  # If trial contains more than one rep, cut trial into separate reps
  # find valleys left and right of peaks to determine start and end of rep
  for (p in 1:length(pks)) {
    # define arbitrary period to look for start and end of rep
    # start of period
    if (p > 1) { # if not first rep
      ds <- min(data$freq*7.5, (pks[p]-pks[p-1])/2) # start of period is 5S before max
    } else {
      ds <- min(data$freq*7.5, pks[p]-1) # else start = either start of trial or 5s before max
    }
    # end of period
    if (p < length(pks)) { # if not last rep
      de <- min(data$freq*5, (pks[p+1]-pks[p])/2) # end of period = halfway between two adjacent maxes
    } else {
      de <- min(data$freq*5, length(data$df$Total)-pks[p]) # else end of period = either end of trial or 5s after max
    }
    w = (pks[p]-ds):(pks[p]+de)
    
    # start and end of each rep are valleys closest to max
    r <- find_rep(data$df$Total[w], (pks[p] - w[1]) + 1, 500)
    rs <- r$start
    re <- r$end
    
    # redefine rep
    # start of rep
    ps <- min(data$freq*1.5, rs-1)
    # end of rep
    pe <- min(data$freq*1.5, length(w)-re) # end of rep = either end of trial or 1.5s after end of pull
    w <- w[rs-ps]:w[re+pe]
    
    if (length(pks) > 1) {
      n <- paste0(data$fn,'_',p)
      data[[n]] <- data$df[w,]
      
      data$fmaxi[p] <- which(data[[n]]$Total == max(data[[n]]$Total))
      if (p == length(pks)) {
        data$df <- NULL
      }
    } else {
      data[[data$fn]] <- data$df[w,]
      data$fmaxi <- which(data[[data$fn]]$Total == max(data[[data$fn]]$Total))
      
      data$df <- NULL
    }
  }
  return(data)
}


# Determine start and end of pull --------------------------------------------------------
pull <- function(data,wp) {
  reps <- names(data)[which(grepl(data$fn,names(data),fixed = T))]
  pull <- matrix(0,length(reps),2)
  colnames(pull) <- c('start','end')
  
  for (r in 1:length(reps)) {
    rn <- reps[r]
    # noise = 5 SD of 1 s weighing period before pull
    # start and end of pull = last and first instance when force < start of pull threshold
    p <- find_rep(data[[rn]]$Total, data$fmaxi[r], data$freq*wp)
    pull[r,1] <- p$start
    pull[r,2] <- p$end
  }
  data$pull <- pull
  return(data)
}


# Check for indicators of poor trial -----------------------------------------------------
qualityCheck <- function(data) {
  reps <- names(data)[which(grepl(data$fn,names(data),fixed = T))]
  data$warn <- list()
  for (r in 1:length(reps)) {
    rn <- reps[r]
    sp <- data$pull[r,1]
    ep <- data$pull[r,2]
    bm <- mean(data[[rn]]$Total[1:(data$freq*1)])
    i <- 1
    
    # unstable weighing period before pull (change in force > 50 N)
    if (min(movstd(data[[rn]]$Total[1:sp],1000), na.rm = T) > 10) {
      data$warn[[rn]][[i]] <- "This rep does not have a stable weighing period"
      i <- i+1
    }
    
    # comparable force before and after pull
    pt <- mean(data[[rn]]$Total[1:sp]) * 0.1 # threshold = 10% of force before pull
    if (mean(data[[rn]]$Total[1:sp]) > (mean(data[[rn]]$Total[ep:length(data[[rn]]$Total)]) + pt) ||
        mean(data[[rn]]$Total[1:sp]) < (mean(data[[rn]]$Total[ep:length(data[[rn]]$Total)]) - pt)) {
      data$warn[[rn]][[i]] <- "This rep does not have equal forces before and after the pull"
      i <- i+1
    }
    
    # peak force during the last 5% of pull
    if (data$fmaxi[r] >= (sp + (ep-sp)*0.95)) {
      data$warn[[rn]][[i]] <- "In this rep peak force occurs at the end of the pull"
      i <- i+1
    }
    
    # countermovement before pull
    # CM threshold = mean force of 1st second - 5 SD of force during 1st second
    cmt <- mean(data[[rn]]$Total[1:(data$freq*1)]) - (sd(data[[rn]]$Total[1:(data$freq*1)])*5)
    if (any(data[[rn]]$Total[(sp - data$freq*0.5):sp] < cmt)) {
      data$warn[[rn]][[i]] <- "This rep contains a countermovement prior to the pull"
      i <- i+1
    }
  }
  return(data)
}


# Performance metrics --------------------------------------------------------------------
perfMetrics <- function(data) {
  reps <- names(data)[which(grepl(data$fn,names(data),fixed = T))]
  pm <- as.data.frame(matrix(0, length(reps),13))
  colnames(pm) <- c('peak force', 'relative peak force', 'time to peak force', paste('RFD',seq(50,250,50)), 
                    paste('J-',seq(50,250,50)))
  rownames(pm) <- reps
  
  for (r in 1:length(reps)) {
    rn <- reps[r]
    sp <- data$pull[r,1]
    ep <- data$pull[r,2]
    # Peak force
    pm[r,1] <- data[[rn]]$Total[data$fmaxi[r]]
    # Relative peak force
    pm[r,2] <- pm[r,1]/(data$bodymass/9.81)
    # Time to peak force
    pm[r,3] <- (data$fmaxi[r] - sp) / data$freq
    # RFD 50, 100, 150, 200 and 250 ms (RFD = change in force / change in time)
    pm[r,4] <- (data[[rn]]$Total[sp+(data$freq*0.05)] - data[[rn]]$Total[sp]) / 0.05
    pm[r,5] <- (data[[rn]]$Total[sp+(data$freq*0.10)] - data[[rn]]$Total[sp]) / 0.10
    pm[r,6] <- (data[[rn]]$Total[sp+(data$freq*0.15)] - data[[rn]]$Total[sp]) / 0.15
    pm[r,7] <- (data[[rn]]$Total[sp+(data$freq*0.20)] - data[[rn]]$Total[sp]) / 0.20
    pm[r,8] <- (data[[rn]]$Total[sp+(data$freq*0.25)] - data[[rn]]$Total[sp]) / 0.25
    # Impulse after 50, 100, 150, 200 and 250 ms
    pm[r,9] <- sum(diff(data[[rn]]$Time[sp:(sp+data$freq*0.05)]) * (head(data[[rn]]$Total[sp:(sp+data$freq*0.05)],-1) + tail(data[[rn]]$Total[sp:(sp+data$freq*0.05)],-1)))/2
    pm[r,10] <- sum(diff(data[[rn]]$Time[sp:(sp+data$freq*0.10)]) * (head(data[[rn]]$Total[sp:(sp+data$freq*0.10)],-1) + tail(data[[rn]]$Total[sp:(sp+data$freq*0.10)],-1)))/2
    pm[r,11] <- sum(diff(data[[rn]]$Time[sp:(sp+data$freq*0.15)]) * (head(data[[rn]]$Total[sp:(sp+data$freq*0.15)],-1) + tail(data[[rn]]$Total[sp:(sp+data$freq*0.15)],-1)))/2
    pm[r,12] <- sum(diff(data[[rn]]$Time[sp:(sp+data$freq*0.20)]) * (head(data[[rn]]$Total[sp:(sp+data$freq*0.20)],-1) + tail(data[[rn]]$Total[sp:(sp+data$freq*0.10)],-1)))/2
    pm[r,13] <- sum(diff(data[[rn]]$Time[sp:(sp+data$freq*0.25)]) * (head(data[[rn]]$Total[sp:(sp+data$freq*0.25)],-1) + tail(data[[rn]]$Total[sp:(sp+data$freq*0.25)],-1)))/2
  }
  # data$pm <- data.frame(pm)
  data$pm <- round(pm[,1:13],2)
  return(data)
}

