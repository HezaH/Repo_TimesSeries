source("header.R")
library(daltoolbox)
library(harbinger)
library(dalevents)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)

# Função de plotagem que exibe e opcionalmente salva o gráfico
plot_motif <- function(motif_data, plot_title = "Motif", save_file = NULL) {
  # Adiciona um índice (proxy para o tempo/posição)
  motif_data <- motif_data %>% mutate(Index = row_number())
  
  p <- ggplot(motif_data, aes(x = Index, y = serie)) +
    geom_line(color = "blue") +
    geom_point(data = filter(motif_data, event), 
               aes(x = Index, y = serie), 
               color = "red", size = 3) +
    labs(title = plot_title, x = "Índice", y = "Valor da Série") +
    theme_minimal(base_size = 16)
  
  print(p)
  
  # Se o parâmetro save_file não for NULL, salva o gráfico em arquivo
  if (!is.null(save_file)) {
    ggsave(filename = save_file, plot = p, width = 10, height = 6, dpi = 300)
    message("Gráfico salvo em: ", save_file)
  }
  
  return(p)
}

gen_data <- function() {
  require(lubridate)
  require(dplyr)
  require(daltoolbox)
  
  cases_motifs <- list()
  
  { #mitdb
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

# Função para gerar o plot da série com os pontos de detecção quando presentes
plot_detection_points <- function(ts_data, detection, title = "Detection Plot") {
  # Carrega o pacote ggplot2, se ainda não estiver carregado
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 é necessário para esta função. Instale-o usando install.packages('ggplot2').")
  
  # Cria um data frame com o índice (proxy para tempo) e os valores da série
  df <- data.frame(
    tempo = seq_along(ts_data),
    valor = ts_data
  )
  
  # Cria o gráfico base: linha azul para a série
  p <- ggplot2::ggplot(df, ggplot2::aes(x = tempo, y = valor)) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::labs(title = title, x = "Tempo", y = "Valor") +
    ggplot2::theme_minimal()
  
  # Se o data frame 'detection' contiver a coluna 'idx' e se houver eventos marcados, sobrepõe os pontos em vermelho.
  if("idx" %in% names(detection) && nrow(detection[detection$event, ]) > 0) {
    idx_events <- detection$idx[detection$event]
    df_det <- data.frame(
      tempo = df$tempo[idx_events],
      valor = ts_data[idx_events]
    )
    p <- p + ggplot2::geom_point(data = df_det, ggplot2::aes(x = tempo, y = valor),
                                 color = "red", size = 3)
  }
  
  return(p)
}

# Carrega os dados de exemplo
cases_motifs <- gen_data()

# Plot e salvamento dos gráficos:
p100 <- plot_motif(cases_motifs$mitdb100, "MITDB Record 100", save_file = "plot_mitdb100.png")
p102 <- plot_motif(cases_motifs$mitdb102, "MITDB Record 102", save_file = "plot_mitdb102.png")
print(p100)
print(p102)

# Exemplo completo: para cada caso em cases_motifs, para cada tratamento e modelo execute:
for (i in 1:length(cases_motifs)) {
  # Extraímos os dados do caso atual:
  data <- cases_motifs[[i]]
  # Cria a coluna de tempo como índice
  data$tempo <- 1:nrow(data)
  # Define a coluna de eventos (aqui apenas FALSE para demonstrar; ajuste conforme necessário)
  data$event <- FALSE
  
  # Dados originais: vetor numérico da série
  ts_numeric <- data$serie
  
  # Tratamento 0: Normalização
  preproc <- minmax()   # ou ts_norm_gminmax()
  preproc <- fit(preproc, ts_numeric)
  ts_data_norm <- transform(preproc, ts_numeric)
  
  # Tratamento 1: Suavização (média móvel)
  window_size <- 5
  ts_ma <- stats::filter(ts_numeric, rep(1/window_size, window_size), sides = 2)
  ts_ma[is.na(ts_ma)] <- ts_numeric[is.na(ts_ma)]
  preproc_ma <- minmax()
  preproc_ma <- fit(preproc_ma, ts_ma)
  ts_data_ma <- transform(preproc_ma, ts_ma)
  
  # Tratamento 2: Diferenciação
  ts_diff <- diff(ts_numeric, differences = 1)
  # Para manter o alinhamento, remova o primeiro elemento das outras variáveis (aqui, usamos data$tempo)
  # Supondo que 'data_filter' seja o próprio 'data'
  data_filter <- data
  # data_filter_diff: removendo o primeiro elemento
  data_filter_diff <- data_filter[-1, ]
  preproc_diff <- minmax()
  preproc_diff <- fit(preproc_diff, ts_diff)
  ts_data_diff <- transform(preproc_diff, ts_diff)
  
  # Tratamento 3: Ajuste de Outliers via IQR
  Q1 <- quantile(ts_numeric, 0.25, na.rm = TRUE)
  Q3 <- quantile(ts_numeric, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  upper_bound <- Q3 + 1.5 * IQR_val
  ts_outlier <- ts_numeric
  ts_outlier[ts_numeric < lower_bound] <- lower_bound
  ts_outlier[ts_numeric > upper_bound] <- upper_bound
  preproc_out <- minmax()
  preproc_out <- fit(preproc_out, ts_outlier)
  ts_data_out <- transform(preproc_out, ts_outlier)
  
  # Cria a lista de tratamentos: adiciona também a série original ("orig")
  treatments <- list(
    orig = ts_numeric,       # Série original sem modificação
    norm = ts_data_norm,     # Tratamento 0: Normalização
    suav = ts_data_ma,       # Tratamento 1: Suavização (média móvel)
    diff = ts_data_diff,     # Tratamento 2: Diferenciação
    out  = ts_data_out       # Tratamento 3: Ajuste de Outliers via IQR
  )
  
  modelos <- list(
    harb     = list(def = function() harbinger(), name = "harb"),
    hdis_sax = list(def = function() hdis_sax(13, 10), name = "hdis_sax"),
    hdis_mp  = list(def = function() hdis_mp(mode = "stamp", w = 10, qtd = 5), name = "hdis_mp"),
    hmo_mp   = list(def = function() hmo_mp(mode = "stamp", w = 10, qtd = 5), name = "hmo_mp"),
    hmo_sax  = list(def = function() hmo_sax(10, 13), name = "hmo_sax")
  )
  
  results <- list()
  
  for (trt_name in names(treatments)) {
    ts_data <- treatments[[trt_name]]
    results[[trt_name]] <- list()
    
    for (mdl in names(modelos)) {
      model_info <- modelos[[mdl]]
      
      # Instancia o modelo usando a função definida
      model_obj <- model_info$def()
      
      # Ajusta o modelo à série tratada
      model_obj <- fit(model_obj, ts_data)
      
      # Executa a detecção dos eventos na série tratada
      detection <- detect(model_obj, ts_data)
      
      # Armazena os resultados (modelo e detecção)
      results[[trt_name]][[mdl]] <- list(
        model = model_obj,
        detection = detection
      )
      
      # Exibe um resumo no console
      cat("Tratamento:", trt_name, " | Modelo:", mdl, "\n")
      print(detection[detection$event, ])
      
      # Gera o plot nativo da detecção utilizando a função de plot criada
      p_obj <- plot_detection_points(ts_data, detection,
                                     title = paste0(mdl, " - ", trt_name))
      print(p_obj)
    }
  }
}
cat("\nAll models processed. Plots should be saved in your working directory.\n")
