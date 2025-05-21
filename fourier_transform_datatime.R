source("header.R")
library(daltoolbox)
library(harbinger)
library(dalevents)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)

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
  serie <- data$serie
  
  # Calcula a FFT da série
  fft_values <- fft(serie)
  
  # Calcula as amplitudes (módulos dos coeficientes complexos)
  amplitudes <- Mod(fft_values)
  
  # Número total de pontos
  n <- length(serie)
  
  # Se você souber a taxa de amostragem (fs, em amostras por segundo), defina-a.
  # Exemplo: fs <- 100
  # Caso contrário, trabalharemos com o índice da frequência.
  
  # Cria o vetor de frequências.
  # Se fs estiver definido, converta índices para frequências reais:
  # freq <- (0:(n-1)) * fs / n
  # Caso contrário, usamos:
  freq <- 0:(n - 1)
  
  # Plot do espectro completo de amplitudes
  plot(freq, amplitudes, type = "h", 
       xlab = "Frequência (ou índice)", 
       ylab = "Amplitude", 
       main = paste("Espectro de Amplitude (FFT) -", nome_serie))
  
  # Limita a visualização até a frequência de Nyquist
  # Se fs estiver definido, a frequência de Nyquist é fs/2; caso contrário, usamos n/2.
  nyquist <- if (exists("fs")) fs / 2 else floor(n / 2)
  plot(freq[1:nyquist], amplitudes[1:nyquist], type = "h",
       xlab = "Frequência (Hz ou índice)", 
       ylab = "Amplitude", 
       main = paste("Espectro de Amplitude (até Nyquist) -", nome_serie))
}
