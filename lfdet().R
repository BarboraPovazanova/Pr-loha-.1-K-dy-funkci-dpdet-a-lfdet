lfdet <- function(daily_discharge, distribution, sri_scale, first_month=1){
  
  library(hydroTSM)
  library(xts)
  library(dplyr)
  library(readr)
  library(lubridate)
  library(SCI)
  library(evd)
  library(ggplot2)
  
  # oznacenia suborov
  SRI_lab <- paste0('SRI',sri_scale,'month')
  LF_lab <- paste0('LF',sri_scale)
  Graph_lab <- paste0('SRI',sri_scale,'_graph')
  
  winDialog(type = 'ok', message = 'Input file needs to be a .txt/.csv with 2 columns date & Qd and tab separator, date must be in dd.mm.yy format')
  daily_discharge <- choose.files(caption = 'Choose a file with daily flow time series')
  
  Qd <- read_table2(daily_discharge )
  colnames(Qd) <- c('datum', 'Qd')
  
  Qd$Qd <- as.numeric(Qd$Qd)
  Qd$datum <- dmy(Qd$datum)
  
  Qd.xts <- xts(x = Qd$Qd,order.by = Qd$datum)
  Qm <- daily2monthly(Qd.xts, FUN = sum)
  Qm <- data.frame(datum=index(Qm), Qm=coredata(Qm))
  
  for (i in 1:length(Qm$Qm)) {
    Qm[i,'Qm'] <- ifelse(Qm[i,'Qm']==0,0.000001,Qm[i,'Qm'])
  }
  
  sri.para <- fitSCI(Qm$Qm,first.mon=first_month,time.scale=sri_scale,distr=distribution,p0=FALSE, 
                     start.fun = dist.start, start.fun.fix = TRUE)
  
  sri <- transformSCI(Qm$Qm,first.mon=first_month,obj=sri.para) 
  
  Qm$SRI <- sri
  
  if (sri_scale != 1) {
    slice_vector <- seq(1,(sri_scale-1),1)  
    Qm <-
      Qm |>
      slice(-slice_vector)
  }
  
  SRInmonth <- Qm
  
  ## Identifikacia LF
  
  SRInmonth$ifLF0 <- ifelse(SRInmonth$SRI<0, 0, 1)
  
  ## Dlzka suchych obdobi 
  
  # funkcia na kumulativny sucet nulovych hodnot
  cumul_zeros <- function(x)  {
    x <- !x
    rl <- rle(x)
    len <- rl$lengths
    v <- rl$values
    cumLen <- cumsum(len)
    z <- x
    iDrops <- c(0, diff(v)) < 0
    z[ cumLen[ iDrops ] ] <- -len[ c(iDrops[-1],FALSE) ]
    x*cumsum(z)
  }
  
  SRInmonth$LF_nmonth <- cumul_zeros(SRInmonth$ifLF0)
  vektor <- numeric(nrow(SRInmonth))
  
  
  for (i in 1:nrow(SRInmonth)) {
    vektor[i] <- ifelse(SRInmonth[i+1, 'LF_nmonth'] == 0, SRInmonth[i, 'LF_nmonth'], 0)
    
  }
  
  SRInmonth$LF_length <- vektor
  
  # pridanie posledneho riadku
  SRInmonth[nrow(SRInmonth), 'LF_length'] <- SRInmonth[nrow(SRInmonth), 'LF_nmonth']
  
  ### Pomenovanie suchych obdobi 
  
  onlyLF <- SRInmonth[SRInmonth$LF_length>0, ] # vyber iba suchych obdobi
  
  LF_name <- numeric(nrow(onlyLF)) # prazdny vektor
  numbers <- 1:nrow(onlyLF) # vektor s cislami 
  
  # priradenie oznacenia kazdemu suchemu obdobiu
  
  for (i in 1:nrow(onlyLF)) {
    LF_name[i] <- paste0('LF',numbers[i])
  }
  
  onlyLF$LF_name <- LF_name # zapisanie do df
  
  SRInmonth <- merge(SRInmonth, onlyLF[, c('datum', 'LF_name')], all = T)
  
  ## Dlzka suchych obdobi 
  
  SRInmonth$datum <- as.Date(SRInmonth$datum)
  onlyLF$datum <- as.Date(onlyLF$datum)
  
  LF_start <- Date(length(onlyLF$datum))
  
  for (i in 1:length(onlyLF$datum)) {
    LF_start[i] <- onlyLF[i,'datum'] %m-% months(onlyLF[i,'LF_length']-1)
  }
  
  onlyLF$LF_start <- LF_start
  
  onlyLF <- 
    onlyLF |> 
    rename(LF_end=datum)
  
  ## Severity 
  
  LF_severity <- numeric(nrow(onlyLF))
  
  for (i in 1:nrow(onlyLF)) {
    LF_severity[i] <- 
      sum(SRInmonth[SRInmonth$datum >= onlyLF[i,'LF_start']&SRInmonth$datum <= onlyLF[i,'LF_end'], 'SRI'])
  }
  onlyLF$LF_severity <- LF_severity
  
  # zapis suboru len s potrebnymi stlpcami
  LF <- 
    onlyLF |> 
    dplyr::select(LF_name,LF_start, LF_end, LF_length, LF_severity )
  
  # graf
  
  LF_ave_date <- lubridate::Date(nrow(LF))
  
  for (i in 1:nrow(LF)) {
    LF_ave_date[i] <- mean.Date(c(LF[i,'LF_start'], LF[i, 'LF_end']))
  }
  LF$LF_ave_date <- LF_ave_date
  
  
  
  SRIgraf <-
    SRInmonth |>
    mutate(pos = SRI >= 0) |>
    ggplot() +
    geom_col(aes(x = datum, y = SRI, fill = pos)) +
    scale_fill_manual(values = c('firebrick2', 'dodgerblue'), guide = 'none') +
    labs(title = paste0('SRI index for scale = ',sri_scale, ' month(s)'), caption = 'calculated with monthly runoff data') +
    xlab('Year') + ylab(SRI_lab) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9)) +
    scale_x_date(name = 'Year', date_breaks = '1 years') +
    scale_y_continuous(breaks = seq(round(min(SRInmonth$SRI, na.rm = T),1), round(max(SRInmonth$SRI, na.rm = T),1), by = 0.5)) +
    geom_hline(yintercept = 0, color = 'red')+
    geom_text(data = LF, mapping = aes(x = LF_ave_date, y = 0.2, label = LF_name), angle = 90, size = 2.5)
  
  SRInmonth <- 
    SRInmonth |> 
    rename(DATE=datum)
  
  colnames(SRInmonth)[colnames(SRInmonth) == "SRI"] =SRI_lab
  LF <- LF |> 
    dplyr::select(-'LF_ave_date')
  
  LFSRIlist <- list(LF, SRInmonth, SRIgraf)
  names(LFSRIlist) <- c(LF_lab, SRI_lab, Graph_lab)
  return(LFSRIlist)
  
}

