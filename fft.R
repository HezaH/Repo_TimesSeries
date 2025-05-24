source("header.R")
library(daltoolbox)
library(harbinger)
library(dalevents)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
options(scipen=999)

fft_harmonics <- function(x, n = NULL) {
  minx <- min(x)
  x <- x - minx
  
  dff = fft(x)
  if (is.null(n)) {
    n <- length(dff)/2 - 1    
  }
  
  t = seq(from = 1, to = length(x))
  ndff = array(data = 0, dim = c(length(t), 1L))
  ndff[1] = dff[1] #Always, it's the DC component
  if(n != 0){
    ndff[2:(n+1)] = dff[2:(n+1)] #The positive frequencies always come first
    #The negative ones are trickier
    ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)]
  }
  indff = fft(ndff/length(x), inverse = TRUE)
  
  return(Mod(indff) + minx)
}

gen_data <- function() {
  require(lubridate)
  require(dplyr)
  require(daltoolbox)
  
  cases_motifs <- list()
  
  { # mitdb
    carrega_data <- function(name) {
      load(sprintf("data-gen/%s.rdata", name))
      load(sprintf("data-gen/%s-anot.rdata", name))
      colnames(data) <- c("v1", "v2")
      data$i <- 1:nrow(data)
      data <- merge(x = data, y = anot, by.x = "i", by.y = "anot", all.x = TRUE)
      data$symbol[is.na(data$symbol)] <- "N"
      data$event <- (data$symbol != "N") & (data$symbol != "/")
      return(data)
    }
    
    compute_range <- function(data, time_str) {
      total <- period_to_seconds(hms("0:30:05"))
      sec   <- period_to_seconds(hms("0:00:10"))
      
      time_val <- period_to_seconds(hms(time_str))
      n <- nrow(data)
      
      # Calcula o índice alvo proporcionalmente ao total de linhas
      idx_target <- time_val * n / total
      sec_idx <- sec * n / total
      
      data_subset <- data[(data$i > idx_target - sec_idx) & (data$i < idx_target + sec_idx), ]
      return(data_subset)
    }
    
    data <- carrega_data("100")
    motif_data <- compute_range(data, "00:25:18")
    motif_data <- motif_data %>% select(serie = v1, event = event, symbol = symbol)
    cases_motifs$mitdb100 <- motif_data
    
    data <- carrega_data("102")
    motif_data <- compute_range(data, "00:04:56")
    motif_data <- motif_data %>% select(serie = v2, event = event, symbol = symbol)
    cases_motifs$mitdb102 <- motif_data
  }
  return(cases_motifs)
}

# Carrega os dados de exemplo
cases_motifs <- gen_data()


# Para cada caso na lista, calcula e plota o espectro da FFT
for (i in 1:length(cases_motifs)) {
  # Obtém o nome da série (se definido na lista)
  nome_serie <- names(cases_motifs)[i]
  # Extrai os dados do caso atual
  data <- cases_motifs[[i]]
  
  # Extraí a série numérica dos dados (coluna "serie")
  y <- data$serie
  
  yts <- ts(y )
  xts <- time(yts)
  
  # Determine the number of harmonics to include based on the significant frequency components
  periodogram <- spec.pgram(serie, plot=TRUE)
  harmonics <- length(periodogram$freq)
  yhat <- fft_harmonics(y, harmonics)
  print(sum(abs(y-yhat)))
  
  grf <- autoplot(ts(yhat, ))
  grf <- grf + theme_bw(base_size = 10)
  grf <- grf + theme(plot.title = element_blank())
  grf <- grf + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
  grf <- grf + ylab("value")
  grf <- grf + xlab("time")
  grf <- grf + geom_point(aes(y=yts),size = 0.5, col="black") 
  grf <- grf + labs(caption = sprintf("(a): %d harmonics", harmonics)) 
  grf <- grf + theme(plot.caption = element_text(hjust = 0.5))
  grf <- grf  + font
  grfa <- grf
  plot(grfa)
  
  df <- data.frame(x = periodogram$freq, y = periodogram$spec)
  harmonics <- as.integer(sqrt(length(periodogram$freq)))
  spec <- df$y[harmonics]
  
  for (i in harmonics:nrow(df)) {
    spec <- df$y[i]
    significant_freq <- which(df$y > spec)
    if(i >= max(significant_freq)) {
      harmonics <- i
      break
    }
  }
  
  # periodogram 
  grf <- ggplot(df, aes(x = x, y = y)) + geom_line() + geom_point(size=0.5) + scale_y_log10()
  grf <- grf + theme_bw(base_size = 10)
  grf <- grf + theme(plot.title = element_blank())
  grf <- grf + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
  grf <- grf + ylab("frequency")
  grf <- grf + xlab("spectrum")  
  grf <- grf + geom_hline(yintercept = spec, col="darkgrey", size = 0.5, linetype="dotted")
  grf <- grf + labs(caption = "(b)") 
  grf <- grf + theme(plot.caption = element_text(hjust = 0.5))
  grf <- grf  + font
  grfb <- grf
  plot(grfb)
  
  yhat <- fft_harmonics(y, harmonics)
  print(c(harmonics, sum(abs(y-yhat))))
  tolerance <- ceiling(0.03*length(yts))
  
  grf <- autoplot(ts(yts ))
  grf <- grf + theme_bw(base_size = 10)
  grf <- grf + theme(plot.title = element_blank())
  grf <- grf + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
  grf <- grf + ylab("value")
  grf <- grf + xlab("time")
  grf <- grf + geom_point(aes(y=yts),size = 0.5, col="black") 
  grf <- grf + geom_line(aes(y=yhat), linetype = "dashed", col="darkblue") 
  grf <- grf + geom_vline(xintercept = xts[tolerance], col="darkgrey", size = 0.5, linetype="dotted")
  grf <- grf + geom_vline(xintercept = xts[length(xts)-tolerance], col="darkgrey", size = 0.5, linetype="dotted")
  grf <- grf + labs(caption = sprintf("(c): %d harmonics", harmonics)) 
  grf <- grf + theme(plot.caption = element_text(hjust = 0.5))
  grf <- grf  + font
  grfc <- grf
  plot(grf)
  
  grf <- autoplot(ts(y - yhat ))
  grf <- grf + theme_bw(base_size = 10)
  grf <- grf + theme(plot.title = element_blank())
  grf <- grf + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
  grf <- grf + ylab("residual")
  grf <- grf + xlab("time")
  grf <- grf + geom_point(size = 0.5, col="black") 
  grf <- grf + geom_vline(xintercept = xts[tolerance], col="darkgrey", size = 0.5, linetype="dotted")
  grf <- grf + geom_vline(xintercept = xts[length(xts)-tolerance], col="darkgrey", size = 0.5, linetype="dotted")
  grf <- grf + labs(caption = sprintf("(d)")) 
  grf <- grf + theme(plot.caption = element_text(hjust = 0.5))
  grf <- grf  + font
  grfd <- grf
  plot(grfd)
  
  save_file = paste0("figures/fft_", nome_serie, ".png")
  mypng(file = save_file, width = 1280, height=1440)
  gridExtra::grid.arrange(grfa, grfb, grfc, grfd, 
                          layout_matrix = matrix(c(1,2,3,4), byrow = TRUE, ncol = 1))
  dev.off() 
  
}
