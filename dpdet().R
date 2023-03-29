dpdet <- function(precip, spi_scale){
  
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggtext)
  library(SPEI)
  library(hydroTSM)
  library(lubridate)
  library(xts)
  
  # oznacenia suborov
  SPI_lab <- paste0('SPI',spi_scale,'month')
  DP_lab <- paste0('DP',spi_scale)
  Graph_lab <- paste0('SPI',spi_scale,'_graph')
  
  winDialog(type = 'ok', message = 'Input file needs to be a .txt/.csv with 2 columns: date & RR and tab separator, date must be in dd.mm.yy format')
  precip <- choose.files(caption = 'Choose a file with daily precipitation time series') 
  
  RR <- read_table2(precip )
  
  colnames(RR) <- c('datum', 'RR')
  RR$datum <- dmy(RR$datum)
  RR$RR <- as.numeric(RR$RR)
  
  RR.xts <- xts(x = RR$RR,order.by = as.POSIXct(RR$datum)) 
  RR.month.xts <- daily2monthly(RR.xts, FUN = sum) 
  RR.month <- data.frame(date=index(RR.month.xts), RR=coredata(RR.month.xts))
  
  for (i in 1:length(RR.month$RR)) {
    RR.month[i,'RR'] <- ifelse(RR.month[i,'RR']==0,0.00001,RR.month[i,'RR'])
  }
  
  ## SPIn vypocet
  
  SPInmonth <- spi(data = RR.month$RR, scale = spi_scale, na.rm = T)
  SPInmonth <- SPInmonth[['fitted']]
  
  RR.month$SPInmonth <- as.vector(SPInmonth)
  
  RR.month <- 
    RR.month |> 
    rename(datum=date)
  
  if (spi_scale != 1) {
    slice_vector <- seq(1,(spi_scale-1),1)  
    RR.month <-
      RR.month |>
      slice(-slice_vector)
  }
  
  SPInmonth <- RR.month
  
  ## Identifikacia DP
  
  SPInmonth$ifDP0 <- ifelse(SPInmonth$SPInmonth<0, 0, 1)
  
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
  
  SPInmonth$DP_nmonth <- cumul_zeros(SPInmonth$ifDP0)
  pokus <- numeric(nrow(SPInmonth))
  
  for (i in 1:nrow(SPInmonth)) {
    pokus[i] <- ifelse(SPInmonth[i+1, 'DP_nmonth'] == 0, SPInmonth[i, 'DP_nmonth'], 0)
  }
  
  SPInmonth$DP_length <- pokus
  
  # pridanie posledneho riadku
  SPInmonth[nrow(SPInmonth), 'DP_length'] <- SPInmonth[nrow(SPInmonth), 'DP_nmonth']
  
  ### Pomenovanie suchych obdobi 
  
  onlyDP <- SPInmonth[SPInmonth$DP_length>0, ] # vyber iba suchych obdobi
  
  DP_name <- numeric(nrow(onlyDP)) 
  numbers <- 1:nrow(onlyDP)  
  
  # priradenie oznacenia kazdemu suchemu obdobiu
  
  for (i in 1:nrow(onlyDP)) {
    DP_name[i] <- paste0('DP',numbers[i])
  }
  
  onlyDP$DP_name <- DP_name 
  
  SPInmonth <- merge(SPInmonth, onlyDP[, c('datum', 'DP_name')], all = T)
  
  
  ## Dlzka suchych obdobi 
  
  SPInmonth$datum <- as.Date(SPInmonth$datum)
  onlyDP$datum <- as.Date(onlyDP$datum)
  
  DP_start <- Date(length(onlyDP$datum))
  
  for (i in 1:length(onlyDP$datum)) {
    DP_start[i] <- onlyDP[i,'datum'] %m-% months(onlyDP[i,'DP_length']-1)
  }
  
  onlyDP$DP_start <- DP_start
  
  onlyDP <- 
    onlyDP |> 
    rename(DP_end=datum)
  
  ## Výdatnosť (severity) 
  
  DP_severity <- numeric(nrow(onlyDP))
  
  for (i in 1:nrow(onlyDP)) {
    DP_severity[i] <- 
      sum(SPInmonth[SPInmonth$datum >= onlyDP[i,'DP_start']&SPInmonth$datum <= onlyDP[i,'DP_end'], 'SPInmonth'])
  }
  onlyDP$DP_severity <- DP_severity
  
  # zapis suboru len s potrebnymi stlpcami
  DP <- 
    onlyDP |> 
    dplyr::select(DP_name,DP_start, DP_end, DP_length, DP_severity )
  
  # graf
  
  DP_ave_date <- lubridate::Date(nrow(DP))
  
  for (i in 1:nrow(DP)) {
    DP_ave_date[i] <- mean.Date(c(DP[i,'DP_start'], DP[i, 'DP_end']))
  }
  DP$DP_ave_date <- DP_ave_date
  
  SPIgraf <-
    SPInmonth |>
    mutate(pos = SPInmonth >= 0) |>
    ggplot() +
    geom_col(aes(x = datum, y = SPInmonth, fill = pos)) +
    scale_fill_manual(values = c('firebrick2', 'dodgerblue'), guide = 'none') +
    labs(title = paste0('SPI index for scale = ',spi_scale, ' month(s)'), caption = 'calculated with monthly precipitation data') +
    xlab('Year') + ylab(SPI_lab) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9)) +
    scale_x_date(name = 'Year', date_breaks = '1 years') +
    scale_y_continuous(breaks = seq(round(min(SPInmonth$SPInmonth, na.rm = T),1), round(max(SPInmonth$SPInmonth, na.rm = T),1), by = 0.5)) +
    geom_hline(yintercept = 0, color = 'red')+
    geom_text(data = DP, mapping = aes(x = DP_ave_date, y = 0.2, label = DP_name), angle = 90, size = 2.5)
  
  SPInmonth <- 
    SPInmonth |> 
    rename(DATE=datum)
  
  colnames(SPInmonth)[colnames(SPInmonth) == "SPInmonth"] =SPI_lab
  
  DP <- DP |> 
    dplyr::select(-'DP_ave_date')
  
  DPSPIlist <- list(DP, SPInmonth, SPIgraf)
  names(DPSPIlist) <- c(DP_lab, SPI_lab, Graph_lab)
  return(DPSPIlist) 
}

