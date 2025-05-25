source("header.R")
library(daltoolbox)
library(harbinger)
library(dalevents)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
options(scipen=999)
library(ggpmisc)

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
    geom_point(data = filter(motif_data, event), 
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
plot_detection_points <- function(ts_data, detection, title = "Detection Plot", 
                                  width = 10, height = 6, dpi = 300) {
  # Carrega o pacote ggplot2, se não estiver carregado
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 é necessário para esta função. Instale-o com install.packages('ggplot2').")
  
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
  
  # Se houver detecções, sobrepõe os pontos em vermelho
  if("idx" %in% names(detection) && nrow(detection[detection$event, ]) > 0) {
    idx_events <- detection$idx[detection$event]
    df_det <- data.frame(
      tempo = df$tempo[idx_events],
      valor = ts_data[idx_events]
    )
    p <- p + ggplot2::geom_point(data = df_det, ggplot2::aes(x = tempo, y = valor),
                                 color = "red", size = 3)
  }
  
  # Se save_file não for NULL, salva o gráfico
  save_file = paste0("figures/", title, ".png")
  ggplot2::ggsave(filename = save_file, plot = p, width = width, height = height, dpi = dpi)
  message("Gráfico salvo em: ", save_file)
  
  return(p)
}

# Carrega os dados de exemplo
cases_motifs <- gen_data()

# Plot e salvamento dos gráficos:
p100 <- plot_motif(cases_motifs$mitdb100, "MITDB Record 100", save_file = "figures/plot_mitdb100.png")
p100_html <- plot_motif_html(cases_motifs$mitdb100, "MITDB Record 100", save_file = "figures/plot_mitdb100.html")

p102 <- plot_motif(cases_motifs$mitdb102, "MITDB Record 102", save_file = "figures/plot_mitdb102.png")
p102_html <- plot_motif_html(cases_motifs$mitdb102, "MITDB Record 102", save_file = "figures/plot_mitdb102.html")

# Exibe os gráficos
print(p100)
print(p100_html)

print(p102_html)
print(p102)

# Exemplo completo: para cada caso em cases_motifs, para cada tratamento e modelo execute:
for (i in 1:length(cases_motifs)) {
  nome_serie <- names(cases_motifs)[i]
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
  
  # Tratamento 1: Global Norma #! Erro
  ngminmax <- ts_norm_gminmax(remove_outliers = FALSE)
  ngminmax <- fit(ngminmax, ts_numeric)
  dgminmax <- transform(ngminmax, ts_numeric)
  
  ### Tratamento 2: Diff com Scaller
  ndiff <- ts_norm_diff(remove_outliers = FALSE)
  ndiff <- fit(ndiff, ts_numeric)
  ddiff <- transform(ndiff, ts_numeric)
  
  ### Tratamento 3: Slidiing Window MinMax
  nswminmax <- ts_norm_swminmax(remove_outliers = FALSE)
  nswminmax <- fit(nswminmax, ts_numeric)
  dswminmax <- transform(nswminmax, ts_numeric)
  
  ### Tratamento 4: AdaptineNormalization
  n_an <- ts_norm_an()
  n_an <- fit(n_an, ts_numeric)
  d_an <- transform(n_an, ts_numeric)
  
  # Tratamento 5: Zscore
  preproc <- zscore()
  preproc <- fit(preproc, ts_numeric)
  ts_data_zscore <- transform(preproc, ts_numeric)
  
  # Tratamento 6: Suavização (média móvel)
  window_size <- 5
  ts_ma <- stats::filter(ts_numeric, rep(1/window_size, window_size), sides = 2)
  ts_ma[is.na(ts_ma)] <- ts_numeric[is.na(ts_ma)]
  preproc_ma <- minmax()
  preproc_ma <- fit(preproc_ma, ts_ma)
  ts_data_ma <- transform(preproc_ma, ts_ma)
  
  # Tratamento 7: Diferenciação
  ts_diff <- diff(ts_numeric, differences = 1)
  # Para manter o alinhamento, remova o primeiro elemento das outras variáveis (aqui, usamos data$tempo)
  # Supondo que 'data_filter' seja o próprio 'data'
  data_filter <- data
  # data_filter_diff: removendo o primeiro elemento
  data_filter_diff <- data_filter[-1, ]
  preproc_diff <- minmax()
  preproc_diff <- fit(preproc_diff, ts_diff)
  ts_data_diff <- transform(preproc_diff, ts_diff)
  
  # Tratamento 8: Ajuste de Outliers via IQR
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
    orig = data$serie,       # Série original sem modificação
    norm = ts_data_norm,     # Normalização
    ngminmax = dgminmax,     # Global Norma #! Erro
    ddiff = ddiff,           # Diff com Scaller
    swminmax = dswminmax,    # Sliding Window MinMax
    an = d_an,               # Adaptive Normalization
    zscore = ts_data_zscore, # Z-score
    ma = ts_data_ma,         # Suavização (média móvel)
    diff = ts_data_diff,     # Diferenciação
    out = ts_data_out        # Ajuste de Outliers via IQR
  )
  
  modelos <- list(
    harb     = list(def = function() harbinger(), name = "harb"),
    hdis_sax = list(def = function() hdis_sax(26, 25), name = "hdis_sax"),
    hdis_mp_stamp  = list(def = function() hdis_mp(mode = "stamp", w = 25, qtd = 1), name = "hdis_mp_stamp"),
    hmo_mp_stamp   = list(def = function() hmo_mp(mode = "stamp", w = 25, qtd = 3), name = "hmo_mp_stamp"),
    hdis_mp_stomp  = list(def = function() hdis_mp(mode = "stomp", w = 25, qtd = 1), name = "hdis_mp_stomp"),
    hmo_mp_stomp   = list(def = function() hmo_mp(mode = "stomp", w = 25, qtd = 3), name = "hmo_mp_stomp"),
    hdis_mp_scrimp  = list(def = function() hdis_mp(mode = "scrimp", w = 25, qtd = 1), name = "hdis_mp_scrimp"),
    hmo_mp_scrimp   = list(def = function() hmo_mp(mode = "scrimp", w = 25, qtd = 3), name = "hmo_mp_scrimp"),
    hmo_sax  = list(def = function() hmo_sax(26, 25), name = "hmo_sax"),
    hmo_xsax  = list(def = function() hmo_xsax(37,3,3), name = "hmo_xsax")
  )
  
  results <- list()
  
  for (trt_name in names(treatments)) {
    ts_data <- treatments[[trt_name]]
    results[[trt_name]] <- list()
    
    for (mdl in names(modelos)) {
      model_info <- modelos[[mdl]]
      
      # Ajusta o modelo à série tratada
      model_obj <- fit(model_info$def(), ts_data)
      
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
      print(detection |> dplyr::filter(event==TRUE))
      if (grepl("hdis_", mdl)){      
        evaluation <- evaluate(model_obj, detection$event, dataset$event)
        print(evaluation$confMatrix)
      }
      
      
      # Gera o plot nativo da detecção utilizando a função de plot criada
      p_obj <- plot_detection_points(ts_data, detection,
                                     title = paste0(nome_serie, mdl, " - ", trt_name))
      print(p_obj)
    }
  }
}
cat("\nAll models processed. Plots should be saved in your working directory.\n")
