source("header.R")
# install.packages("devtools")
library(devtools)
library(daltoolbox) 
library(harbinger) 
#library(dalevents)
library(htmlwidgets)
library(tspredit)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
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

# Função de plot atualizada para receber true_event e desenhar em vermelho
plot_detection_points_html <- function(ts_data, detection, true_event,
                                       title = "Detection Plot", 
                                       width = 10, height = 6, dpi = 300, 
                                       save_file = NULL) {
  if (!requireNamespace("ggplot2", quietly=TRUE)) stop("Instale ggplot2")
  if (!requireNamespace("plotly", quietly=TRUE))  stop("Instale plotly")
  if (!requireNamespace("htmlwidgets", quietly=TRUE)) stop("Instale htmlwidgets")
  
  vec <- as.numeric(ts_data)
  df  <- data.frame(
    tempo = seq_along(vec),
    valor = vec,
    true_event = as.logical(true_event),
    detect_event = detection$event
  )
  
  p <- ggplot2::ggplot(df, aes(x=tempo, y=valor)) +
    ggplot2::geom_line(color="blue") +
    ggplot2::labs(title=title, x="Tempo", y="Valor") +
    ggplot2::theme_minimal()
  
  # pontos preditos (roxo)
  df_pred <- df[df$detect_event, , drop=FALSE]
  if (nrow(df_pred) > 0) {
    p <- p + ggplot2::geom_point(
      data = df_pred, aes(x=tempo, y=valor),
      color = "purple", size = 3, inherit.aes = FALSE
    )
  }
  
  # pontos reais (vermelho)
  df_true <- df[df$true_event, , drop=FALSE]
  if (nrow(df_true) > 0) {
    p <- p + ggplot2::geom_point(
      data = df_true, aes(x=tempo, y=valor),
      color = "red", size = 2, shape = 4, inherit.aes = FALSE
    )
  }
  
  interactive_plot <- plotly::ggplotly(p)
  print(interactive_plot)
  
  if (!is.null(save_file)) {
    htmlwidgets::saveWidget(interactive_plot, file=save_file, selfcontained=TRUE)
    message("Gráfico salvo em HTML: ", save_file)
  }
  
  return(interactive_plot)
}

gen_data <- function() {
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
start_time <- Sys.time()
# Carrega os dados de exemplo
cases_motifs <- gen_data()

# Exemplo completo: para cada caso em cases_motifs, para cada tratamento e modelo execute:
for (i in 1:length(cases_motifs)) {
  nome_serie <- names(cases_motifs)[i]
  
  # Extraímos os dados do caso atual:
  data <- cases_motifs[[i]]
  event_info <- data$event
  
  # Cria a coluna "tempo" com o índice original
  data$tempo <- 1:nrow(data)
  
  num_true <- sum(data$event, na.rm = TRUE)
  if (num_true < 3){
    num_true = 3
  }
  
  # Definindo o intervalo dos valores da série conforme informado (-2.7 a 1.3)
  min_val <- min(data$serie, na.rm = TRUE)
  max_val <- max(data$serie, na.rm = TRUE)
  
  
  # Cria os pontos de quebra (breaks)
  # Para 5 grupos, precisamos de 6 pontos, por isso usamos length.out = 6
  breaks <- seq(min_val, max_val, length.out = 7)
  
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
  
  # Plot e salvamento dos gráficos:
  p1 <- plot_motif(data, paste0("MITDB Record ", nome_serie), 
                   save_file = paste0("figures/", nome_serie, ".png"))
  p1_html <- plot_motif_html(data, paste0("MITDB Record ", nome_serie), 
                             save_file = paste0("figures/", nome_serie, ".html"))
  
  # Se desejar, defina o argumento 'save_file' com o nome (e caminho) do arquivo HTML a ser salvo.
  plot_html <- plot_motif_html(data,
                               plot_title = paste("Grupo com valores mais altos para", nome_serie),
                               save_file = paste0("figures/max_values_", nome_serie, ".html"))
  
  # Exibe os gráficos
  print(p1)
  print(p1_html)
  print(plot_html)
  print(data %>% filter(!is.na(diff_index)))
  
  media_diff <- data %>% 
    filter(!is.na(diff_index)) %>% 
    summarise(mean_diff = mean(diff_index, na.rm = TRUE))
  
  media_diff <- as.integer(media_diff[[1]])
  print(media_diff)
  
  # Dados originais: vetor numérico da série
  ts_numeric <- data$serie
  data_size <- length(ts_numeric)
  ts_numeric_complete <- ts_data(data$serie, data_size)
  
  # Tratamento 0: Normalização
  preproc <- minmax()   # ou ts_norm_gminmax()
  preproc <- fit(preproc, ts_numeric)
  ts_data_norm <- transform(preproc, ts_numeric)
  
  # Tratamento 1: Global Norma 
  ngminmax <- ts_norm_gminmax()
  ngminmax <- fit(ngminmax, ts_numeric)
  dgminmax <- transform(ngminmax, ts_numeric)
  
  ### Tratamento 2: Diff com Scaller
  ndiff <- ts_norm_diff()
  ndiff <- fit(ndiff, ts_numeric_complete)
  ddiff <- transform(ndiff, ts_numeric_complete)
  
  ### Tratamento 3: Slidiing Window MinMax
  nswminmax <- ts_norm_swminmax()
  nswminmax <- fit(nswminmax, ts_numeric_complete)
  dswminmax <- transform(nswminmax, ts_numeric_complete)
  
  ### Tratamento 4: AdaptineNormalization
  n_an <- ts_norm_an()
  n_an <- fit(n_an, ts_numeric_complete)
  d_an <- transform(n_an, ts_numeric_complete)
  
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
    ngminmax = dgminmax,     # Global Norma 
    ddiff = ddiff,           # Diff com Scaller
    swminmax = dswminmax,    # Sliding Window MinMax
    an = d_an,               # Adaptive Normalization
    zscore = ts_data_zscore, # Z-score
    ma = ts_data_ma,         # Suavização (média móvel)
    diff = ts_data_diff,     # Diferenciação
    out = ts_data_out        # Ajuste de Outliers via IQR
  )
  
  if (nome_serie == "mitdb100") { #discord 100
    modelos <- list(
      harb     = list(def = function() harbinger(), name = "harb"),
      hdis_sax = list(def = function() hdis_sax(26, media_diff), name = "hdis_sax"),
      hdis_mp_stamp  = list(def = function() hdis_mp(mode = "stamp", w = media_diff, qtd = num_true), name = "hdis_mp_stamp"),
      hdis_mp_stomp  = list(def = function() hdis_mp(mode = "stomp", w = media_diff, qtd = num_true), name = "hdis_mp_stomp"),
      hdis_mp_scrimp  = list(def = function() hdis_mp(mode = "scrimp", w = media_diff, qtd = num_true), name = "hdis_mp_scrimp")
    )
  } else if (nome_serie == "mitdb102") { #motif 102
    modelos <- list(
      harb     = list(def = function() harbinger(), name = "harb"),
      hdis_sax = list(def = function() hdis_sax(26, media_diff), name = "hdis_sax"),
      hdis_mp_stamp  = list(def = function() hdis_mp(mode = "stamp", w = media_diff, qtd = num_true), name = "hdis_mp_stamp"),
      hmo_mp_stamp   = list(def = function() hmo_mp(mode = "stamp", w = media_diff, qtd = num_true), name = "hmo_mp_stamp"),
      hdis_mp_stomp  = list(def = function() hdis_mp(mode = "stomp", w = media_diff, qtd = num_true), name = "hdis_mp_stomp"),
      hmo_mp_stomp   = list(def = function() hmo_mp(mode = "stomp", w = media_diff, qtd = num_true), name = "hmo_mp_stomp"),
      hdis_mp_scrimp  = list(def = function() hdis_mp(mode = "scrimp", w = media_diff, qtd = num_true), name = "hdis_mp_scrimp"),
      hmo_mp_scrimp   = list(def = function() hmo_mp(mode = "scrimp", w = media_diff, qtd = num_true), name = "hmo_mp_scrimp"),
      hmo_sax  = list(def = function() hmo_sax(26, media_diff), name = "hmo_sax"),
      hmo_xsax  = list(def = function() hmo_xsax(37,3,num_true), name = "hmo_xsax")
    )
  }
  
  results <- list()
  for (trt_name in names(treatments)) {
    # extrai o vetor puro e o vetor true_event do dataset
    raw_vec   <- as.numeric(treatments[[trt_name]])
    true_evt  <- data$event     # vetor lógico original do dataset
    
    # se o tratamento for diff/subtração (ddiff) e reduzir 1 ponto, alinha o true_event:
    if (trt_name == "ddiff" || trt_name == "diff") {
      # diff de ordem 1 reduz comprimento em 1
      if (length(true_evt) == length(raw_vec) + 1) {
        true_evt <- true_evt[-1]  # remove primeiro para alinhar com diff
      }
    }
    
    results[[trt_name]] <- list()
    
    for (mdl in names(modelos)) {
      cat(sprintf("Processing treatment: %s with model: %s\n", trt_name, mdl))
      model_info <- modelos[[mdl]]
      
      # Ajusta o modelo
      model_obj <- fit(model_info$def(), raw_vec)
      
      # Detect com alinhamento de comprimento
      detection <- detect(model_obj, raw_vec)
      if (length(detection$event) == length(raw_vec) - 1) {
        detection$event <- c(FALSE, detection$event)
      }
      if (length(detection$event) != length(raw_vec)) {
        stop(sprintf("Comprimento inválido: detect retornou %d flags, mas vetor tem %d pontos",
                     length(detection$event), length(raw_vec)))
      }
      
      # Armazena resultados
      results[[trt_name]][[mdl]] <- list(
        model     = model_obj,
        detection = detection,
        true_evt  = true_evt
      )
      
      # Plot predições vs reais
      p_obj <- plot_detection_points_html(
        ts_data     = raw_vec,
        detection   = detection,
        true_event  = true_evt,
        title       = paste0(nome_serie, " - ", mdl, " (", trt_name, ")"),
        save_file   = paste0("figures/", nome_serie, "_", mdl, "_", trt_name, ".html")
      )
      print(p_obj)
    }
  }
  
}
cat("\nAll models processed. Plots should be saved in your working directory.\n")
# Marca fim
end_time <- Sys.time()

# Calcula diferença
total_time <- end_time - start_time
print(total_time)
