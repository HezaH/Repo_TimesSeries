source("header.R")
library(daltoolbox)
library(harbinger)
library(dalevents)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)

plot_motif_html <- function(motif_data, plot_title = "Motif", save_file = NULL) {
  # Carrega/garante que os pacotes necessários estão carregados
  library(dplyr)
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
  
  # Adiciona um índice para representar a posição ou o tempo
  motif_data <- motif_data %>% mutate(Index = row_number())
  
  # Cria o gráfico utilizando ggplot2
  p <- ggplot(motif_data, aes(x = Index, y = serie)) +
    geom_line(color = "blue") +
    geom_point(data = filter(motif_data, !is.na(diff_index)), 
               aes(x = Index, y = serie), 
               color = "red", size = 3) +
    labs(title = plot_title, x = "Índice", y = "Valor da Série") +
    theme_minimal(base_size = 16)
  
  # Converte o gráfico para o formato interativo plotly
  interactive_plot <- ggplotly(p)
  
  # Exibe o gráfico interativo na sessão atual
  print(interactive_plot)
  
  # Se save_file não for NULL, salva o gráfico interativo em um arquivo HTML
  if (!is.null(save_file)) {
    saveWidget(interactive_plot, file = save_file, selfcontained = TRUE)
    message("Gráfico salvo em: ", save_file)
  }
  
  return(interactive_plot)
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
      
      data_subset <- data[ (data$i > idx_target - sec_idx) & (data$i < idx_target + sec_idx), ]
      return(data_subset)
    }
    
    # Exemplo para 2 casos
    data <- carrega_data("100")
    motif_data <- compute_range(data, "00:25:18")
    motif_data <- motif_data |> select(serie = v1, event = event, symbol = symbol)
    cases_motifs$mitdb100 <- motif_data
    
    data <- carrega_data("102")
    motif_data <- compute_range(data, "00:04:56")
    motif_data <- motif_data |> select(serie = v2, event = event, symbol = symbol)
    cases_motifs$mitdb102 <- motif_data
  }
  return(cases_motifs)
}

# Carrega os dados de exemplo
cases_motifs <- gen_data()

# Para cada caso, agrupa os valores, extrai o grupo de maior valor 
# e plota os dados mantendo o eixo x completo (respeitando os índices ausentes)
for (i in 1:length(cases_motifs)) {
  nome_serie <- names(cases_motifs)[i]
  # Extraímos os dados do caso atual
  data <- cases_motifs[[i]]
  
  # Cria a coluna "tempo" com o índice original
  data$tempo <- 1:nrow(data)
  
  # Define a coluna de eventos (aqui, apenas FALSE para demonstração)
  data$event <- FALSE
  
  # Definindo o intervalo dos valores da série conforme informado (-2.7 a 1.3)
  min_val <- min(data$serie, na.rm = TRUE)
  max_val <- max(data$serie, na.rm = TRUE)
  
  # Cria os pontos de quebra (breaks)
  # Para 5 grupos, precisamos de 6 pontos, por isso usamos length.out = 6
  breaks <- seq(min_val, max_val, length.out = 6)
  
  # A função cut categoriza os valores da coluna "serie" conforme esses intervalos
  data$group <- cut(data$serie, breaks = breaks, include.lowest = TRUE, right = TRUE)
  
  # Conta a frequência de ocorrências em cada grupo
  freq_groups <- table(data$group)
  
  # Exibe o nome da série e os resultados 
  cat("Frequência dos grupos para", nome_serie, ":\n") 
  print(freq_groups) 
  cat("\n")
  
  # Identifica o grupo com valores mais altos (o último nível)  
  highest_group <- tail(levels(data$group), n = 1)
  
  # Cria uma nova coluna que mantém os valores se pertencem ao grupo mais alto, ou NA caso contrário
  data$serie_high <- ifelse(data$group == highest_group, data$serie, NA)
  
  # Inicializa a variável auxiliar e a coluna diff_index
  last_time <- 0
  data$diff_index <- NA  # Cria/limpa a coluna diff_index
  
  # Vamos percorrer o dataframe usando um loop while
  j <- 1
  n <- nrow(data)
  while(j <= n) {
    if (!is.na(data$serie_high[j])) {
      # Se o valor não for NA, é o começo de um bloco de valores consecutivos
      block_start <- j
      # Avança enquanto os valores em serie_high forem não NA
      while(j <= n && !is.na(data$serie_high[j])) {
        j <- j + 1
      }
      block_end <- j - 1  # Último índice do bloco
      # Vetor com os índices do bloco
      block_indices <- block_start:block_end
      # Extrai os valores desse bloco
      block_values <- data$serie_high[block_indices]
      # Identifica, dentro do bloco, o índice relativo onde o valor absoluto é o maior
      max_val_index <- which.max(abs(block_values))
      # Converte para o índice global do dataframe
      global_row <- block_indices[max_val_index]
      # Calcula a diferença: índice atual menos o último índice registrado
      data$diff_index[global_row] <- data$tempo[global_row] - last_time
      # Atualiza last_time para o índice atual
      last_time <- data$tempo[global_row]
    } else {
      # Se for NA, avança para a próxima linha
      j <- j + 1
    }
  }
  
  print(data)
  
  # Se desejar, defina o argumento 'save_file' com o nome (e caminho) do arquivo HTML a ser salvo.
  plot_html <- plot_motif_html(data,
                               plot_title = paste("Grupo com valores mais altos para", nome_serie),
                               save_file = paste0("figures/max_values_", nome_serie, ".html"))
  
  print(plot_html)
  print(data %>% filter(!is.na(diff_index)))
  
  media_diff <- data %>% 
    filter(!is.na(diff_index)) %>% 
    summarise(mean_diff = mean(diff_index, na.rm = TRUE))
  
  print(media_diff)
}
