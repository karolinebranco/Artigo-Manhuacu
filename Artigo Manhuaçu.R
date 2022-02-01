rm(list=ls()) #Limpa a memória do R

#Carregando os pacotes necessários
library(evd)
library(R2OpenBUGS)
library(readxl)
library(Kendall) 
library(coda)

#Carregando o conjunto de dados
dados <- read_excel("~/unifal/Mestrado/Análise de Extremos/dados_trabalho_final.xlsx")
View(dados)

#Verificando as estatísticas descritivas 
summary(dados$max_anual)

#Selecionando a coluna com os máximos anuais para a realização dos testes de independência e tendência
precipitacao <- dados$max_anual

#======================================================================
##    Verificando se as observações são independentes
##                   Testes de Ljung-Box 
#======================================================================

Box.test(precipitacao,type=("Ljung-Box")) #Se for maior que 5% os dados são independentes

#======================================================================
# Verificando Tendência Teste Mannkendall #######
##======================================================================

MannKendall(precipitacao) #Se for maior que 5% não há evidência de tendência nos dados

#======================================================================
# Distribuição GEV sem Tendência
#======================================================================
#Priori Não informativa
sink("Gev_1_NI.txt")
cat("
    model {
    
    # Verossimilhança
    for (i in 1:n) {
    y[i] ~ dgev(mu, sigma, eta)
    }
   
    # Prior
    mu ~ dnorm(0,0.0001) # coloca a média e a precisão, lembrando que a precisão é 1/variância
    sigma ~ dlnorm(-4.60517,0.217142) #precisa fazer ln dos valores de média e de precisão
    eta ~ dnorm(0,0.01)
    
  # Para Fazer preditiva
 
nr2<- mu + ((sigma/eta)*(pow(-log(1-1/2),-eta)-1))
nr4<- mu + ((sigma/eta)*(pow(-log(1-1/4),-eta)-1))
nr6<- mu + ((sigma/eta)*(pow(-log(1-1/6),-eta)-1))
nr10<- mu + ((sigma/eta)*(pow(-log(1-1/10),-eta)-1))
 
}
    ",fill=TRUE)
sink()

y <- dados$max_anual[1:41] #Dados de 1970 a 2010 para estimar os parâmetros da GEV e os níveis de retorno
dados_bug<- list(y=y,n=length(y))
dados_bug
fgev(y) #Estimativas de verossimilhança 

inits <- function(){ list(mu=66.08, sigma=15.15, eta=0.0701)} #Criando uma veriável com as estimativas de verossimilhança para serem empregadas na inferência

params <- c("mu","sigma", "eta")
nc = 1 #Numero de cadeias
ni = 27000 #Tamanho da cadeia quase final
nb = 7000 #Numero de simulação que serão descartadas
nt =  4#Thinning rate

#======================================================================
# Iniciar o Amostrador para a geração das cadeias posteriori dos parâmetros da GEV

gev.bayes_L1 = bugs(data = dados_bug, inits = inits, 
                    parameters =c(params,"nr2","nr4","nr6","nr10"),
                    model = "Gev_1_NI.txt",
                    n.thin = nt, n.chains = nc,
                    n.burnin = nb, n.iter = ni, codaPkg=FALSE, debug=T)

print(gev.bayes_L1, dig = 4) #Mostra as estimativas médias para os parâmetros e níveis de retorno

post_gb_L1<-as.mcmc(gev.bayes_L1$sims.matrix[,]) # salva a saída como cadeia mcmc

#Intervalo de credibilidade HPD
HPDinterval(post_gb_L1)

#Critérios de Convergência das cadeias posteriori dos parâmetros e níveis de retorno
geweke.diag(post_gb_L1) # Geweke
raftery.diag(post_gb_L1) # Raftery & Lewis
heidel.diag(post_gb_L1) # Heidelberger & Welch  

##### Calculando EMP - Erro Médio de Predição
resumo2=print(gev.bayes_L1,dig=4) #Arredondando os resultados das estimativas para 4 casas decimais 
ypred<-c(resumo2$mean$nr2,resumo2$mean$nr4,resumo2$mean$nr6,resumo2$mean$nr10) #Salva o valor médio predito para cada nível de retorno

VP = ypred;VP #Salvando ypred na nova variável VP
Vobs = c(97.9, 97.9,99.1,105.1) #Valores máximos observados na série de dados

EpGev= abs((Vobs-VP)/Vobs) #Calculando o Erro Médio de Predição (EMP)

round(mean(EpGev)*100,2) #Colocando o EMP em % e em duas casas decimais

################################Priori informativa##############################
################################ Priori Lavras #################################
#Médias: mu=68,9, sigma=2,90 e eta=-0,05
#Variâncias: mu=10,48, sigma= 5,64 e eta=0,02
#Vamos considerar a variância original e a flexibilização da variância 2 e 4 vezes
#No modelo criado, fornece as informações de média e precisão para cada parâmetro
#A precisão é o inverso da variância


sink("Gev_Lavras.txt")
cat("
    model {
    
    # Verossimilhança
    for (i in 1:n) {
    y[i] ~ dgev(mu, sigma, eta)
    }
   
    # Prior - usa a precosao ao invés da variância - a precisao é o inverso da variancia
    mu ~ dnorm(68.9,  0.09541985) 
    sigma ~ dlnorm( 1.0647, 0.4127062) #quando usa dlnorm eu preciso aplicar log no valor dos parâmetros informativos, no caso, da precisão, por exemplo é 1/log(informacao)
    eta ~ dnorm(-0.05,50) 
    
  # Para Fazer preditiva
 
nr2<- mu + ((sigma/eta)*(pow(-log(1-1/2),-eta)-1))
nr4<- mu + ((sigma/eta)*(pow(-log(1-1/4),-eta)-1))
nr6<- mu + ((sigma/eta)*(pow(-log(1-1/6),-eta)-1))
nr10<- mu + ((sigma/eta)*(pow(-log(1-1/10),-eta)-1))
 
}
    ",fill=TRUE)
sink()

fgev(y) #escolher chutes iniciais

inits <- function(){ list(mu=66.08, sigma=15.15, eta=0.07601)} #Chutes iniciais

params <- c("mu","sigma", "eta")
nc = 1 #Numero de cadeias
ni = 38000 #Tamanho da cadeia quase final
nb = 7500 #Numero de simulação que serão descartadas
nt =  5#Thinning rate

#======================================================================
# Inicie o Amostrador

gev.bayes_L1 = bugs(data = dados_bug, inits = inits, 
                    parameters =c(params,"nr2","nr4","nr6","nr10"),
                    model = "Gev_Lavras.txt",
                    n.thin = nt, n.chains = nc,
                    n.burnin = nb, n.iter = ni, codaPkg=FALSE, debug=T)

print(gev.bayes_L1, dig = 4)

post_gb_L1<-as.mcmc(gev.bayes_L1$sims.matrix[,]) # salva a saída como cadeia mcmc

HPDinterval(post_gb_L1) 
geweke.diag(post_gb_L1)
raftery.diag(post_gb_L1)
heidel.diag(post_gb_L1)

##### Calculando EMP - Erro Médio de Predição
resumo2=print(gev.bayes_L1,dig=4)
ypred<-c(resumo2$mean$nr2,resumo2$mean$nr4,resumo2$mean$nr6,resumo2$mean$nr10) 

VP = ypred;VP
Vobs = c(97.9, 97.9,99.1,105.1)

EpGev= abs((Vobs-VP)/Vobs)

round(mean(EpGev)*100,2)

################################Priori informativa##############################
################################ Priori Juiz de Fora #################################
#Médias: mu=76,8169, sigma=19,2028 e eta=-0,0088
#Vamos considerar a variância original e a flexibilização da variância 2 e 4 vezes
#Variâncias: mu= 9,8804, sigma= 5,220759 e eta= 0,01309623


sink("Gev_JF.txt")
cat("
    model {
    
    # Verossimilhança
    for (i in 1:n) {
    y[i] ~ dgev(mu, sigma, eta)
    }
   
    # Prior - usa a precosao ao invés da variância - a precisao é o inverso da variancia
    mu ~ dnorm(76.8179,   0.1012105) 
    sigma ~ dlnorm(2.9551,  0.6050914) 
    eta ~ dnorm(-0.0088, 76.35785) 
    
  # Para Fazer preditiva
 
nr2<- mu + ((sigma/eta)*(pow(-log(1-1/2),-eta)-1))
nr4<- mu + ((sigma/eta)*(pow(-log(1-1/4),-eta)-1))
nr6<- mu + ((sigma/eta)*(pow(-log(1-1/6),-eta)-1))
nr10<- mu + ((sigma/eta)*(pow(-log(1-1/10),-eta)-1))
 
}
    ",fill=TRUE)
sink()

fgev(y) 
inits <- function(){ list(mu=66.08, sigma=15.15, eta=0.0701)} 

params <- c("mu","sigma", "eta")
nc = 1 #Numero de cadeias
ni = 55000 #Tamanho da cadeia quase final
nb = 10000 #Numero de simulação que serão descartadas
nt =  10#Thinning rate

#======================================================================
# Inicie o Amostrador

gev.bayes_L1 = bugs(data = dados_bug, inits = inits, 
                    parameters =c(params,"nr2","nr4","nr6","nr10"),
                    model = "Gev_JF.txt",
                    n.thin = nt, n.chains = nc,
                    n.burnin = nb, n.iter = ni, codaPkg=FALSE, debug=T)

print(gev.bayes_L1, dig = 4)

library(coda)

post_gb_L1<-as.mcmc(gev.bayes_L1$sims.matrix[,]) 

HPDinterval(post_gb_L1) 
geweke.diag(post_gb_L1)
raftery.diag(post_gb_L1)
heidel.diag(post_gb_L1)

##### Calculando EMP - Erro Médio de Predição
resumo2=print(gev.bayes_L1,dig=4)
ypred<-c(resumo2$mean$nr2,resumo2$mean$nr4,resumo2$mean$nr6,resumo2$mean$nr10) 

VP = ypred;VP
Vobs = c(97.9, 97.9,99.1,105.1)

EpGev= abs((Vobs-VP)/Vobs)

round(mean(EpGev)*100,2)

######## Usando o modelo escolhido para a realização das previsões, 
######## Aplicando na série histórica completa (1970-2020)

#Médias: mu=76,8169, sigma=19,2028 e eta=-0,0088
#Variâncias: mu= 9,8804, sigma= 5,220759 e eta= 0,01309623


sink("Gev_modelo_escolhido.txt")
cat("
    model {
    
    # Verossimilhança
    for (i in 1:n) {
    y[i] ~ dgev(mu, sigma, eta)
    }
   
    # Prior - usa a precosao ao invés da variância - a precisao é o inverso da variancia
    mu ~ dnorm(76.8179, 0.05060524) 
    sigma ~ dlnorm(2.9551, 0.4262956) #quando usa dlnorm eu preciso aplicar ln no valor dos parâmetros informativos, no caso, da precisão, por exemplo é 1/log(informacao)
    eta ~ dnorm(-0.0088, 38.17893) 
    
  # Para Fazer preditiva
 
nr2<- mu + ((sigma/eta)*(pow(-log(1-1/2),-eta)-1))
nr5<- mu + ((sigma/eta)*(pow(-log(1-1/5),-eta)-1))
nr10<- mu + ((sigma/eta)*(pow(-log(1-1/10),-eta)-1))
nr25<- mu + ((sigma/eta)*(pow(-log(1-1/25),-eta)-1))
nr50<- mu + ((sigma/eta)*(pow(-log(1-1/50),-eta)-1))
nr100<- mu + ((sigma/eta)*(pow(-log(1-1/100),-eta)-1))
nr200<- mu + ((sigma/eta)*(pow(-log(1-1/200),-eta)-1))
 
}
    ",fill=TRUE)
sink()
y <- dados$max_anual[1:51]
dados_bug<- list(y=y,n=length(y))
dados_bug
fgev(y) #escolher chutes iniciais

inits <- function(){ list(mu=67.37, sigma=15.90467, eta=0.04911)} #Chutes iniciais

params <- c("mu","sigma", "eta")
nc = 1 #Numero de cadeias
ni = 55000 #Tamanho da cadeia quase final
nb = 8000 #Numero de simulação que serão descartadas
nt =  8#Thinning rate

#======================================================================
# Inicie o Amostrador

gev.bayes_L1 = bugs(data = dados_bug, inits = inits, 
                    parameters =c(params,"nr2","nr5","nr10","nr25","nr50","nr100","nr200"),
                    model = "Gev_modelo_escolhido.txt",
                    n.thin = nt, n.chains = nc,
                    n.burnin = nb, n.iter = ni, codaPkg=FALSE, debug=T)

print(gev.bayes_L1, dig = 4)

library(coda)
post_gb_L1<-as.mcmc(gev.bayes_L1$sims.matrix[,]) # salva a saída como cadeia mcmc
HPDinterval(post_gb_L1) 
geweke.diag(post_gb_L1)
raftery.diag(post_gb_L1)
heidel.diag(post_gb_L1)











