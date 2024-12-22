#### PACOTE --------------------------------------------------------------------

suppressWarnings({
  library(data.table)
  library(tidyverse)
  library(ggfortify)
  library(forecast)
  library(tseries)
  library(readxl)
  library(lmtest)
  library(urca)
  library(zoo)
})

#### FUNÇÕES -------------------------------------------------------------------

# Leitura e formatação dos dados da EPE
leitura_dados <- function(regiao,arquivo, sheet){
  consumo <- read_xls(paste0("../Dados/",arquivo),
                      sheet = sheet)
  
  consumo <- consumo[-c(1:3),]
  
  names(consumo) <- paste0("V_",seq(1:dim(consumo)[2]))
  
  consumo <- consumo %>% 
    mutate(
      is_ano = ifelse(nchar(str_remove(V_2,"\\*")) == 4,T,F),
      is_mes = ifelse(nchar(V_2) == 3,T,F)) %>% 
    filter(is_ano == T | is_mes == T | V_1 == regiao) %>% 
    mutate(base = ifelse(is.na(V_1),V_2,V_1))
  
  if (all(!regiao %in% c("Sudeste","Centro-Oeste"))) {
    consumo <- consumo %>% 
      group_by(base) %>%
      mutate(row_n = row_number(),
             row_n = ifelse(row_n > 2, 
                            ifelse(row_n%%2 > 0,1,2),
                            row_n),
             V_2 = str_remove(V_2,"\\*")) %>% 
      filter(row_n == 1) %>% 
      pmap_dfr(., ~ na.locf(c(...)) %>%
                 as.list %>%
                 as_tibble) %>% 
      select(-c(V_1, V_14,is_ano,is_mes,base,row_n))
  } else {
    consumo <- consumo %>%
      mutate(V_2 = str_remove(V_2,"\\*")) %>% 
      pmap_dfr(., ~ na.locf(c(...)) %>%
                 as.list %>%
                 as_tibble) %>% 
      select(-c(V_1, V_14,is_ano,is_mes,base))
  }
  
  names(consumo) <- consumo[2,]
  
  consumo <- consumo %>% 
    filter(JAN != "JAN") %>% 
    mutate(is_ano = ifelse(nchar(JAN) == 4,T,F)) %>% 
    gather("key","value",-is_ano) %>% 
    mutate(
      is_ano = ifelse(nchar(value) == 4,paste0(str_to_title(key),"/",value,"/01"),NA),
      is_ano = na.locf(is_ano)) %>% 
    filter(nchar(value) > 4) %>% 
    mutate(
      Data = as.Date(is_ano,"%b/%Y/%d"),
      value = as.numeric(value),
      Regiao = regiao) %>% 
    select(Regiao,Data,value) %>% 
    arrange(Data)
  
  return(consumo)
}

# Cria gráficos da serie e as previsões individuais
plot_series <- function(nordeste_df,previsao,nome){
  # Convertendo o objeto de previsão para um data frame com datas
  previsao_df <- data.frame(
    time = as.Date(as.yearmon(time(previsao$mean))), # Converte para Date
    value = as.numeric(previsao$mean),
    lower_80 = as.numeric(previsao$lower[, 1]),
    upper_80 = as.numeric(previsao$upper[, 1]),
    lower_95 = as.numeric(previsao$lower[, 2]),
    upper_95 = as.numeric(previsao$upper[, 2])
  )
  
  cols <- c("Observado" = "black","Previsto" = "#f04546")
  
  ggplot() +
    geom_line(data = nordeste_df, aes(x = time, y = value, color = "Observado"), linetype = "solid") +
    geom_line(data = previsao_df, aes(x = time, y = value, color = "Previsto"), size = .7) +
    geom_vline(xintercept = min(previsao_df$time), color = "blue", linetype = "dashed") +
    geom_ribbon(data = previsao_df, aes(x = time, ymin = lower_80, ymax = upper_80), fill = "gray", alpha = 0.3) +
    geom_ribbon(data = previsao_df, aes(x = time, ymin = lower_95, ymax = upper_95), fill = "gray", alpha = 0.2) +
    labs(
      x = "Tempo",
      y = "Valores",
      title = paste0("Série temporal observada e previsão ",nome),
      color = ""
    ) +
    xlim(c(as.Date("2020-01-01"),as.Date("2024-10-01")))+
    scale_colour_manual(values=cols)+
    theme_light()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position=c(.12,.85))
}

#### LEITURA DOS DADOS ---------------------------------------------------------

nordeste <- leitura_dados("Nordeste",
                          "CONSUMO MENSAL DE ENERGIA ELÉTRICA POR CLASSE.xls",
                          "RESIDENCIAL")

nordeste$value <- (nordeste$value)/1000 # Transformando MWh em GWh 

nordeste_all.ts <- ts(nordeste$value, start = c(2004, 1), frequency = 12)

# Amostra teste (2004-01-01 até 2023-12-01)
nordeste.ts <- subset(nordeste_all.ts, end = 240) 

# Amostra de validação ( > 2023-12-01)
nordeste.ts.prev <- subset(nordeste_all.ts, start = 241)

summary(nordeste.ts) # Estatísticas sumárias da série
sd(nordeste.ts) # Desvio padrão da série

#### ESTACIONARIEDADE E PLOTS DA SÉRIE -----------------------------------------

# Gráfico da série com os correlogramas
ggtsdisplay(nordeste.ts, 
            theme = theme_light()+theme(panel.grid.minor = element_blank(),
                                        panel.grid.major.x = element_blank(),))

# Teste ADF: Rejeitamos a estacionariedade pois  0.3866 > -2,88 
ur.fit = ur.df(nordeste.ts, type="drift", lags=12, selectlags="AIC")
summary(ur.fit) 

# Gráfico da série com os correlogramas após aplicar a primeira diferença
nordeste_diff <- diff(nordeste.ts)
ggtsdisplay(nordeste_diff, 
            theme = theme_light()+theme(panel.grid.minor = element_blank(),
                                        panel.grid.major.x = element_blank(),))

# Testando estacionariedade com a primeira diferença
# Não rejeitamos a estacionariedade pois -4.3404 < -2,88 
ur.fit2 = ur.df(nordeste_diff, type="drift", lags=12, selectlags="AIC")
summary(ur.fit2) 

#### TESTES PARA ENCONTRAR 'd' e 'D' -------------------------------------------

# Critérios para a escolha dos parâmetros:
# d: valor suficiente para tornar a série estacionária
# p: lags ultrapassam os limites no PACF 
# q: lags ultrapassam os limites no ACF

# D: valor suficiente para retirar a sazonalidade da série
# P: Identificar lags sazonais significativos no PACF 
# Q: Identificar lags sazonais significativos no ACF

# Testes que escolhem o valor d
data.frame(
  teste = c("ADF", "KPSS", "PP"),
  valor.d = c(ndiffs(nordeste.ts, test = "adf"),
              ndiffs(nordeste.ts, test = "kpss"),
              ndiffs(nordeste.ts, test = "pp"))
)

# Teste que escolhem o valor de D
data.frame(
  teste = c("OCSB", "SEAS", "CH"),
  valor.D = c(nsdiffs(nordeste.ts, test = "ocsb"),
              nsdiffs(nordeste.ts, test = "seas"),
              nsdiffs(nordeste.ts, test = "ch"))
)

#### AJUSTANDO MODELOS SARIMA --------------------------------------------------

# SARIMA(1,1,0)(2,0,0) with drift
ajuste1 <- auto.arima(nordeste.ts, d=1, D=0)

# Airline(0,1,1)(0,1,1)
ajuste2 <- auto.arima(nordeste.ts, d=1, D=1)

# SARIMA(1,0,1)(0,1,2) with drift
ajuste3 <- auto.arima(nordeste.ts)

# SARIMA(1,1,1)(1,1,1)
ajuste4 <- arima(nordeste.ts, order = c(1,1,1), seasonal = c(1,1,1))

#### DIAGNÓSTICO DOS MODELOS SARIMA --------------------------------------------

diag1 <- ggtsdiag(ajuste1, gof.lag=15) +
  theme_light()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

diag2 <- ggtsdiag(ajuste2, gof.lag=15) +
  theme_light()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

diag3 <- ggtsdiag(ajuste3, gof.lag=15) +
  theme_light()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

# SARIMA(1,1,1)(1,1,1) passa no teste de Ljung-Box até o nono lag
diag4 <- ggtsdiag(ajuste4, gof.lag=15) +
  theme_light()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

diag2[[1]] <- diag2[[1]] + ggtitle("Residuos Padronizados")
diag2[[2]] <- diag2[[2]] + ggtitle("ACF dos Resíduos")
diag2[[3]] <- diag2[[3]] + ylab("p-valor") + 
  ggtitle("p-valores para a estatística Ljung-Box")

diag4[[1]] <- diag4[[1]] + ggtitle("Residuos Padronizados")
diag4[[2]] <- diag4[[2]] + ggtitle("ACF dos Resíduos")
diag4[[3]] <- diag4[[3]] + ylab("p-valor") + 
  ggtitle("p-valores para a estatística Ljung-Box")

diag2
diag4

# Teste de Breusch-Pagan
bptest(residuals(ajuste2) ~ fitted(ajuste2), studentize = F)
# Teste de Koenker
bptest(residuals(ajuste2) ~ fitted(ajuste2), studentize = T)

# Teste de Breusch-Pagan
bptest(residuals(ajuste4) ~ fitted(ajuste2), studentize = F)
# Teste de Koenker
bptest(residuals(ajuste4) ~ fitted(ajuste2), studentize = T)

#### ALISAMENTO EXPONENCIAL DE HOLT WINTERS ------------------------------------

alisamento_add <- HoltWinters(nordeste.ts, seasonal = "additive")

alisamento_mul <- HoltWinters(nordeste.ts, seasonal = "multiplicative",
                              optim.start = c(alpha = 0.3047762, 
                                              beta = 0.006364901, 
                                              gamma = 0.4088794))

#### PREDIÇÃO ------------------------------------------------------------------

pred_add <- forecast(alisamento_add, h = 10) # HW aditivo 
pred_mul <- forecast(alisamento_mul, h = 10) # HW multiplicativo

pred_AIRLINE <- forecast(ajuste2, h = 10) # Airline(0,1,1)(0,1,1)
pred_SARIMA <- forecast(ajuste4, h = 10) # SARIMA(1,1,1)(1,1,1)

# Acurácia do modelo
acuracia <- bind_rows(
  accuracy(forecast(pred_add, h = 10), nordeste.ts.prev) %>% 
    as.data.frame() %>% mutate(Modelo = "Holt Winters Aditivo"),
  accuracy(forecast(pred_mul, h = 10), nordeste.ts.prev) %>% 
    as.data.frame() %>% mutate(Modelo = "Holt Winters Multiplicativo"),
  accuracy(forecast(pred_AIRLINE, h = 10), nordeste.ts.prev) %>% 
    as.data.frame() %>% mutate(Modelo = "Airline(0,1,1)(0,1,1)"),
  accuracy(forecast(pred_SARIMA, h = 10), nordeste.ts.prev) %>% 
    as.data.frame() %>% mutate(Modelo = "SARIMA(2,1,1)(0,1,2)")) %>% 
  mutate(Conjunto = rep(c("Treinamento","Validação"),4)) %>% 
  select(Conjunto,Modelo,RMSE,MAPE)

acuracia %>% 
  filter(Conjunto == "Validação")

#### GRÁFICOS DE PREDIÇÃO ------------------------------------------------------

# Dataframe da série
nordeste_df <- data.frame(
  time = as.Date(as.yearmon(time(nordeste_all.ts))), # Converte para Date
  value = as.numeric(nordeste_all.ts)
)

# Plots
p1 <- plot_series(nordeste_df,pred_mul,"Holt Winters Multiplicativo")+
  xlim(c(as.Date("2022-12-01"),as.Date("2024-10-01")))+
  ylim(c(2200,3500))

p2 <- plot_series(nordeste_df,pred_AIRLINE,"Airline(0,1,1)(0,1,1)")+
  xlim(c(as.Date("2022-12-01"),as.Date("2024-10-01")))+
  ylim(c(2200,3500))

p3 <- plot_series(nordeste_df,pred_SARIMA,"SARIMA(1,1,1)(1,1,1)")+
  xlim(c(as.Date("2022-12-01"),as.Date("2024-10-01")))+
  ylim(c(2200,3500))

gridExtra::grid.arrange(p1,p2,p3,ncol = 2)

# Dataset para comparação dos modelos
previsao_df <- data.frame(
  time = as.Date(as.yearmon(time(pred_mul$mean))), # Converte para Date
  `HW Mult.` = as.numeric(pred_mul$mean),
  `Airline` = as.numeric(pred_AIRLINE$mean),
  `SARIMA` = as.numeric(pred_SARIMA$mean)
) %>% 
  gather("key","value",-time) %>% 
  mutate(key = ifelse(key == "HW.Mult.","HW Mult.",key))

# Plot para comparação dos modelos
ggplot() +
  geom_line(data = nordeste_df, aes(x = time, y = value), size = 1, linetype = "dashed") +
  geom_line(data = previsao_df, aes(x = time, y = value, color = key),size = .6) +
  geom_vline(xintercept = min(previsao_df$time), color = "blue", linetype = "dashed") +
  labs(
    x = "Tempo",
    y = "Valores",
    color = ""
  ) +
  xlim(c(as.Date("2023-12-01"),as.Date("2024-10-01")))+
  ylim(c(2200,3700))+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position=c(.9,.80))
