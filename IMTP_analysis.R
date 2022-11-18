#-----------------------------------------------------------------------------------------
#
# Load and analyze IMTP data (*.csv files)
# Robert Schuster (ACU SPRINT)
# July 2022
#
#-----------------------------------------------------------------------------------------


## clear environment
rm(list = ls())


# Functions ------------------------------------------------------------------------------
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
# select files to analyze
flist <- choose.files(caption = "Choose IMTP files")

# load data and header separately
skip <- list()
data <- list()
cols <- c("Time","Z.Left","Z.Right","Left","Right")
header <- list()
for (f in flist) {
  fn <- basename(f)
  # determine length of header
  skip[[fn]] <- read.csv(f, nrows = 20, header = F, blank.lines.skip = F,
                   col.names = paste0("V",seq_len(max(count.fields(f, sep = ',')))), fill = T)
  if (any(skip[[fn]] == "Time", na.rm = T)) {
    skip[[fn]] <- which(skip[[fn]] == 'Time')-1
    
    # load data (without header) and keep only Time, Left and Right columns
    data[[fn]] <- read.csv(f, skip = skip[[fn]])
    if (any(names(data[[fn]]) %in% cols, na.rm = T)) {
      data[[fn]] <- data[[fn]][,which(names(data[[fn]]) %in% cols)] # keep only left, right and time columns
      colnames(data[[fn]]) <- gsub('Z.','',colnames(data[[fn]])) # remove 'Z." from column names
      
      data[[fn]]$Total <- data[[fn]]$Left + data[[fn]]$Right
      
      # load header
      header[[fn]] <- read.csv(f, nrows = skip[[fn]]-1, header = F)

    } else {
      # print warning but continue execution with other files
      warning(paste("The file ", fn, " does not contain the required data and was not processed"))
      flist <- flist[-which(flist == f)] # remove the file from file list
    }
  } else {
    # print warning but continue execution with other files
    warning(paste("The file", fn, "does not contain the required data and was not processed"))
    flist <- flist[-which(flist == f)] # remove the file from file list
  }
}
rm(skip,cols,f,fn)

fpath <- unique(dirname(flist))
fnames <- basename(flist)

bodymass <- list()
freq <- list()
for (f in fnames) {
  # find weight for force normalisation
  if (any(header[[f]] == "Weight", na.rm = T)) {
    bodymass[[f]] <- as.numeric(na.omit(header[[f]][header[[f]] == 'Weight',2])) * 9.81
  } else {
    bodymass[[f]] <- median(data[[f]]$Total[data[[f]]$Total >= 300])
  }
  # find frequency
  freq[[f]] <- as.numeric(na.omit(header[[f]][header[[f]] == 'Frequency',2]))
  if (freq[[f]] < 1000) {
    message(paste(f,": Sampling frequency =",freq[[f]],"\nFrequencies > 1000 Hz are recommended"))
  }
}
rm(f,header)


# Determine number of repetition per file ------------------------------------------------
fmaxi <- list()
for (f in fnames) {
  # find peaks (at least 2s apart from each other)
  d <- freq[[f]]*2 # frequency * 2s
  pks <- find_peaks(data[[f]]$Total, m = d)
  plot(x = data[[f]]$Time, y = data[[f]]$Total, type = "l", lwd = 2, 
       xlab = "Time [s]", ylab = "Force [N]", main = f)
  points(x = data[[f]]$Time[c(1:nrow(data[[f]]))[pks]], y = data[[f]]$Total[pks], col = "red")
  # determine a cutoff under which all other peaks cannot be considered IMTP trials
  cutoff <- (max(data[[f]]$Total[pks])-bodymass[[f]])/2 + bodymass[[f]] # half of BM normalised max force
  # cutoff <- max(data[[fn]]$Total[pks])-250 # within 250 N of max value
  abline(h = cutoff, col = "green", lty = "dashed")
  # only keep peaks above cutoff and more than 5s apart
  pks <- pks[which(data[[f]]$Total[pks] > cutoff)]
  d <- freq[[f]]*5 # frequency * 5s
  if (any(diff(pks) < d)) {
    pks <- pks[-(which(diff(pks) < d)+1)]
  }
  points(x = data[[f]]$Time[c(1:nrow(data[[f]]))[pks]], y = data[[f]]$Total[pks], col = "red", pch = 16)

  
  # find valleys left and right of peaks to determine start and end of rep
  for (p in 1:length(pks)) {
    # define arbitrary period to look for start and end of rep
    # start of period
    if (p > 1) { # if not first rep
      ds <- min(freq[[f]]*7.5, (pks[p]-pks[p-1])/2) # start of period is 5S before max
    } else {
      ds <- min(freq[[f]]*7.5, pks[p]-1) # else start = either start of trial or 5s before max
    }
    
    # end of period
    if (p < length(pks)) { # if not last rep
      de <- min(freq[[f]]*5, (pks[p+1]-pks[p])/2) # end of period = halfway between two adjacent maxes
    } else {
      de <- min(freq[[f]]*5, length(data[[f]]$Total)-pks[p]) # else end of period = either end of trial or 5s after max
    }
    w = (pks[p]-ds):(pks[p]+de)
    
    # start and end of each rep
    r <- find_rep(data[[f]]$Total[w], (pks[p] - w[1]) + 1, 500)
    rs <- r$start
    re <- r$end
    
    # redefine rep
    # start of period
    ps <- min(freq[[f]]*1.5, rs-1)
    # end of period
    pe <- min(freq[[f]]*1.5, length(w)-re) # end of rep = either end of trial or 1.5s after end of pull
    w <- w[rs-ps]:w[re+pe]
    
    if (length(pks) > 1) {
      n <- paste0(f,'_',p)
      data[[n]] <- data[[f]][w,]
      
      bodymass[[n]] <- bodymass[[f]]
      freq[[n]] <- freq[[f]]
      
      fmaxi[[n]] <- which(data[[n]]$Total == max(data[[n]]$Total))
      if (p == length(pks)) {
        data[[f]] <- NULL
        fmaxi[[f]] <- NULL
        freq[[f]] <- NULL
        bodymass[[f]] <- NULL
      }
    } else {
      data[[f]] <- data[[f]][w,]
      fmaxi[[f]] <- which(data[[f]]$Total == max(data[[f]]$Total))
    }
  }
}
fmaxi <- fmaxi[match(names(data), names(fmaxi))]
rm(f,d,pks,cutoff,p,ds,de,w,rs,re,ps,pe,n)


# Determine start and end of pull --------------------------------------------------------
pull <- matrix(0,length(data),2)
colnames(pull) <- c('start','end')
for (f in 1:length(data)) {
  fn <- names(data)[f]
  # noise = 5 SD of 1 s weighing period before pull
  # start and end of pull = last and first instance when force < start of pull threshold
  p <- find_rep(data[[fn]]$Total, fmaxi[[fn]], freq[[fn]]*1)
  pull[f,1] <- p$start
  pull[f,2] <- p$end
}
rm(f,fn,p)


# Check for indicators of poor trial -----------------------------------------------------

ptq <- matrix('_',length(data),4)
colnames(ptq) <- c('Unstable weighing','Unequal forces','Peak force at end of pull','Countermovement')
rownames(ptq) <- names(data)
for (f in 1:length(data)) {
  sp <- pull[f,1]
  ep <- pull[f,2]
  bm <- mean(data[[f]]$Total[1:(freq[[f]]*1)])
  
  plot(x = data[[f]]$Time, y = data[[f]]$Total, type = "l", lwd = 2, 
       xlab = "Time [s]", ylab = "Force [N]", main = names(data)[f])
  abline(h = bm, col = "green", lwd = 2, lty = "dashed")
  abline(v = data[[f]]$Time[sp], col = "red", lwd = 2, lty = "dotted")
  abline(v = data[[f]]$Time[ep], col = "red", lwd = 2, lty = "dotted")
  points(x = data[[f]]$Time[fmaxi[[f]]], y = data[[f]]$Total[fmaxi[[f]]], col = "blue", pch = 8, lwd = 2)
  
  # unstable weighing period before pull (change in force > 50 N)
  if (min(movstd(data[[f]]$Total[1:sp],1000), na.rm = T) > 10) {
    ptq[f,1] <- 'X'
    # message(paste("Warning:", names(data)[f], "does not have a stable weighing period"))
    # # ask whether to continue processing or not
    # resp <- readline("Would you like to continue processing anyway? (Y/N) ")
    # if (grepl(resp, 'n', ignore.case = TRUE)) {
    #   stop("Processing stopped")
    # }
  }
  
  # comparable force before and after pull
  pt <- mean(data[[f]]$Total[1:sp]) * 0.1 # threshold = 10% of force before pull
  if (mean(data[[f]]$Total[1:sp]) > (mean(data[[f]]$Total[ep:length(data[[f]]$Total)]) + pt) ||
      mean(data[[f]]$Total[1:sp]) < (mean(data[[f]]$Total[ep:length(data[[f]]$Total)]) - pt)) {
    ptq[f,2] <- 'X'
    # message(paste("Warning:", names(data)[f], "does not have equal forces before and after the pull"))
    # # ask whether to continue processing or not
    # resp <- readline("Would you like to continue processing anyway? (Y/N) ")
    # if (grepl(resp, 'n', ignore.case = TRUE)) {
    #   stop("Processing stopped")
    # }
  }
  
  # peak force during the last 5% of pull
  if (fmaxi[[f]] >= (sp + (ep-sp)*0.95)) {
    ptq[f,3] <- 'X'
    # message(paste("Warning: In", names(data)[f], "peak force occurs at the end of the pull"))
    # # ask whether to continue processing or not
    # resp <- readline("Would you like to continue processing anyway? (Y/N) ")
    # if (grepl(resp, 'n', ignore.case = TRUE)) {
    #   stop("Processing stopped")
    # }
  }
  
  # countermovement before pull
  # CM threshold = mean force of 1st second - 5 SD of force during 1st second
  cmt <- mean(data[[f]]$Total[1:(freq[[f]]*1)]) - (sd(data[[f]]$Total[1:(freq[[f]]*1)])*5)
  if (any(data[[f]]$Total[(sp - freq[[f]]*0.5):sp] < cmt)) {
    ptq[f,4] <- 'X'
    # message(paste("Warning:", names(data)[f], "contains a countermovement prior to the pull"))
    # # ask whether to continue processing or not
    # resp <- readline("Would you like to continue processing anyway? (Y/N) ")
    # if (grepl(resp, 'n', ignore.case = TRUE)) {
    #   stop("Processing stopped")
    # }
  }
}
rm(f,sp,ep,bm,pt,cmt,resp)


# Extract performance metrics ------------------------------------------------------------
# impulse
# relative peak force
RFD <-  matrix(0,length(data),5)
colnames(RFD) <- paste('RFD',seq(50,250,50))
IMP <-  matrix(0,length(data),5)
colnames(IMP) <- paste('J-',seq(50,250,50))
pF <-  matrix(0,length(data),3)
colnames(pF) <- c('peak force','relative peak force','time to peak force')
for (f in 1:length(data)) {
  sp <- pull[f,1]
  ep <- pull[f,2]
  # RFD 50, 100, 150, 200 and 250 ms (RFD = change in force / change in time)
  RFD[f,1] <- (data[[f]]$Total[sp+(freq[[f]]*0.05)] - data[[f]]$Total[sp]) / 0.05
  RFD[f,2] <- (data[[f]]$Total[sp+(freq[[f]]*0.10)] - data[[f]]$Total[sp]) / 0.10
  RFD[f,3] <- (data[[f]]$Total[sp+(freq[[f]]*0.15)] - data[[f]]$Total[sp]) / 0.15
  RFD[f,4] <- (data[[f]]$Total[sp+(freq[[f]]*0.20)] - data[[f]]$Total[sp]) / 0.20
  RFD[f,5] <- (data[[f]]$Total[sp+(freq[[f]]*0.25)] - data[[f]]$Total[sp]) / 0.25
  # Impulse from 50, 100, 150, 200 and 250 ms
  IMP[f,1] <- sum(diff(data[[f]]$Time[sp:(sp+freq[[f]]*0.05)]) * (head(data[[f]]$Total[sp:(sp+freq[[f]]*0.05)],-1) + tail(data[[f]]$Total[sp:(sp+freq[[f]]*0.05)],-1)))/2
  IMP[f,2] <- sum(diff(data[[f]]$Time[sp:(sp+freq[[f]]*0.10)]) * (head(data[[f]]$Total[sp:(sp+freq[[f]]*0.10)],-1) + tail(data[[f]]$Total[sp:(sp+freq[[f]]*0.10)],-1)))/2
  IMP[f,3] <- sum(diff(data[[f]]$Time[sp:(sp+freq[[f]]*0.15)]) * (head(data[[f]]$Total[sp:(sp+freq[[f]]*0.15)],-1) + tail(data[[f]]$Total[sp:(sp+freq[[f]]*0.15)],-1)))/2
  IMP[f,4] <- sum(diff(data[[f]]$Time[sp:(sp+freq[[f]]*0.20)]) * (head(data[[f]]$Total[sp:(sp+freq[[f]]*0.20)],-1) + tail(data[[f]]$Total[sp:(sp+freq[[f]]*0.20)],-1)))/2
  IMP[f,5] <- sum(diff(data[[f]]$Time[sp:(sp+freq[[f]]*0.25)]) * (head(data[[f]]$Total[sp:(sp+freq[[f]]*0.25)],-1) + tail(data[[f]]$Total[sp:(sp+freq[[f]]*0.25)],-1)))/2
  # Peak force
  pF[f,1] <- data[[f]]$Total[fmaxi[[f]]]
  # Relative peak force
  pF[f,2] <- pF[f,1]/(bodymass[[f]]/9.81)
  # Time to peak force
  pF[f,3] <- (fmaxi[[f]] - sp) / freq[[f]]
}
rm(f,sp,ep)


# Export results -------------------------------------------------------------------------
ex <- readline("Do you want to export the results? (Y/N) ")
if (grepl(ex,'y',ignore.case = T)) {
  # tbl <- cbind(pF,RFD,IMP)
  tbl <- cbind(pF,RFD,IMP,ptq)
  rownames(tbl) <- names(data)
  # file name
  fn <- readline("Enter the name of the file you want to save the results to: ")
  # file path
  fp <- dirname(flist[1])
  
  # write.xlsx(tbl, file = paste0(fp,'/',fn,'.xlsx'))
  write.csv(tbl, file = paste0(fp,'/',fn,'.csv'))
  rm(tbl)
}

