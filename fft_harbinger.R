source("header.R")
library(daltoolbox)
library(harbinger)
library(dalevents)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
options(scipen=999)


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
    yts <- ts(y, start = c(1850, 1))
    xts <- time(yts)

    compute_cut_index <- function(freqs) {
    cutindex <- which.max(freqs)
    if (min(freqs) != max(freqs)) {
        threshold <- mean(freqs) + 2.698 * sd(freqs)
        freqs[freqs < threshold] <- min(freqs) + max(freqs)
        cutindex <- which.min(freqs)
    }
    return(cutindex)
    }

    fft_signal <- stats::fft(yts)

    spectrum <- base::Mod(fft_signal) ^ 2
    half_spectrum <- spectrum[1:(length(yts) / 2 + 1)]

    cutindex <- compute_cut_index(half_spectrum)
    print(cutindex)
    n <- length(fft_signal)

    fft_signal[1:cutindex] <- 0
    fft_signal[(n - cutindex):n] <- 0


    residual <- - base::Re(stats::fft(fft_signal, inverse = TRUE) / n)

    yhat <- yts - residual


    grf <- autoplot(ts(yts, start = c(1850, 1)))
    grf <- grf + theme_bw(base_size = 10)
    grf <- grf + theme(plot.title = element_blank())
    grf <- grf + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
    grf <- grf + ylab("value")
    grf <- grf + xlab("time")
    grf <- grf + geom_point(aes(y=yts),size = 0.5, col="black") 
    grf <- grf + geom_line(aes(y=yhat), linetype = "dashed", col="darkblue") 
    grf <- grf + theme(plot.caption = element_text(hjust = 0.5))
    grf <- grf  + font
    grfc <- grf
    plot(grf)

    grf <- autoplot(ts(residual, start = c(1850, 1)))
    grf <- grf + theme_bw(base_size = 10)
    grf <- grf + theme(plot.title = element_blank())
    grf <- grf + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
    grf <- grf + ylab("residual")
    grf <- grf + xlab("time")
    grf <- grf + geom_point(size = 0.5, col="black") 
    grf <- grf + labs(caption = sprintf("(d)")) 
    grf <- grf + theme(plot.caption = element_text(hjust = 0.5))
    grf <- grf  + font
    grfd <- grf
    plot(grfd)

    save_file = paste0("figures/fft_harbinger_", nome_serie, ".png")
    mypng(file = save_file, width = 1280, height=1440)
    gridExtra::grid.arrange(grfc, grfd, 
                            layout_matrix = matrix(c(1,2), byrow = TRUE, ncol = 1))
    dev.off() 
}