# Instalar pacotes se necessário
# install.packages("ggplot2")
library(ggplot2)

# 1. Gerar série temporal
set.seed(42)
n <- 60
tempo <- 1:n
serie <- sin(0.2 * tempo) + rnorm(n, mean = 0, sd = 0.2)

# 2. Normalização z-score
serie_norm <- scale(serie)

# 3. PAA (bloco de 5)
bloco <- 5
blocos <- seq(1, n, by = bloco)
paa <- sapply(blocos, function(i) mean(serie_norm[i:(i + bloco - 1)]))
paa_rep <- rep(paa, each = bloco)

# 4. SAX com 3 símbolos: breakpoints para distribuição normal: [-0.43, 0.43]
breakpoints <- c(-Inf, -0.43, 0.43, Inf)
symbols <- c("a", "b", "c")
sax <- cut(paa, breaks = breakpoints, labels = symbols, include.lowest = TRUE)
sax_labels <- rep(as.character(sax), each = bloco)

# 5. Criar data frame para plot
df <- data.frame(
  Tempo = tempo,
  Serie = as.numeric(serie_norm),
  PAA = paa_rep,
  SAX = sax_labels
)

# 6. Plot
ggplot(df, aes(x = Tempo)) +
  geom_line(aes(y = Serie), color = "gray70", size = 1.2) +
  geom_line(aes(y = PAA), color = "orange", linetype = "dashed", size = 1) +
  geom_text(aes(y = 2, label = SAX), size = 5) +
  labs(title = "Representação PAA e SAX em Série Temporal",
       y = "Valor (normalizado)") +
  ylim(-2.5, 2.5) +
  theme_minimal()