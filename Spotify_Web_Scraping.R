################################################
# Início : 03/04/2019                          #
# Última modificação : 12/01/20                #
# Autor : Pedro Henrique Sodré Puntel          #
# Email : pedro.puntel@gmail.com               #
# Tema : Praticando Web Scraping com o Spotify #
# Script Econding : UTF-8                      #
################################################

################
## Referências :
################
# > https://cran.r-project.org/web/packages/Rspotify/Rspotify.pdf

##############
## Descrição :
##############
# Este script aborda processos para coleta de dados do Spotify similar àquele feito com
# o site da Billboard. Neste, por meio da utilização do pacote pacote "Rspotify",
# desenvolvido por um ex-aluno da minha mesma instituição  de ensino, Tiago Mendes Dantas,
# será criada uma rotina que extraí músicas da chart "Global Top 50" e "Top 50 By Country",
# no intuito de investigar o quão "único" é o gosto musical de cada país, bem como o seus
# níveis de afinidade tratando-se de gosto musical.
#
# Adicionalmente, a rotina será também responsável pela criação de uma sócio-matriz a qual
# servirá de base para a aplicação de futuros modelos estatísticos.

##########
## Setup :
##########
# Pacotes
library("rio")
library("dplyr")
library("Rspotify")
library("compiler")

# Autenticação com API do Spotify
token_data = import(file = file.choose(), format = "csv") %>% as.data.frame()
my_auth = spotifyOAuth(token_data$App_Name, token_data$Public_Key, token_data$Private_Key)

#################################
## Countries Playlists Analysis :
#################################
# Data Frame com as "Unique IDs" das charts que serão scrapeadas
global_top50 = import(file = file.choose(), format = "csv") %>% as.data.frame()
country_top_50 = import(file = file.choose(), format = "csv") %>% as.data.frame()

# Rotina que implementa o processo de coleta
spotify_countries_scrap = function() {
  
  # Obtenção atualizada da playlist "Global Top 50"
  global_df = matrix(NA, nrow=50, ncol=1) %>% as.data.frame()
  global_df = getPlaylistSongs("spotifycharts", global_top50[1,1], token=my_auth)$tracks[1:50]
  
  # Obtenção atualizada da playlist "Top 50 by Country"
  top50_df = matrix(NA, nrow=50, ncol=50) %>% as.data.frame()
  for(i in 1:ncol(top50_df)) {
    playlist_songs = getPlaylistSongs("spotifycharts", country_top_50$Country_Top50_Track_ID[i], token=my_auth)
    top50_df[1,i] = country_top_50$Country[i]
    top50_df[,i] = playlist_songs$tracks[1:50]
  }
  colnames(top50_df) = country_top_50$Country
  
  # Montagem da sócio-matriz associada à chart "Top 50 by Country"
  top50_social_matrix = matrix(NA, nrow=50, ncol=50)
  for (i in 1:50) {
    for (j in 1:50) {
      counter = length(intersect(top50_df[,i], top50_df[,j]))
      top50_social_matrix[i,j] = counter
    }
  }
  colnames(top50_social_matrix) = colnames(top50_df)
  rownames(top50_social_matrix) = colnames(top50_df)
  
  # Número de músicas em comum que cada país compartilha com a chart "Global Top 50"
  #    >> "Quão similar ao gosto musical mundial é o gosto de cada país?" <<
  unique_index = matrix(NA, nrow=50, ncol=2) %>% as.data.frame()
  for(i in 1:50) {
    
    # Inverso da proporção de músicas da chart global que pertencem à chart do país
    unique_index[i,2] = 1/((length(intersect(global_df, top50_df[, i])))/50)
  }
  colnames(unique_index) = c("Country","Score")
  unique_index$Country = as.vector(colnames(top50_social_matrix))
  unique_index$Score = as.numeric(as.vector(unique_index$Score))
  
  # Normaliza e ordena de forma decrescente do índice para melhor entendimento
  unique_index$Score = (unique_index$Score-min(unique_index$Score))/(max(unique_index$Score)-min(unique_index$Score))
  unique_index = unique_index[order(unique_index$Score, decreasing=T),]
  
  # Retorna as informações
  return(list("Global_Top_50" = global_df,
              "Top50_Country" = top50_df,
              "Top50_Social_Matrix" = top50_social_matrix,
              "Country_Unique_Index" = unique_index))
  
}
spotify_countries_scrap = cmpfun(spotify_countries_scrap)

# Extraindo os dados
spotify_data = spotify_countries_scrap()

