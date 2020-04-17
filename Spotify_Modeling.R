######################################################################
# Início : 15/01/2020                                                #
# Última modificação : 15/01/2020                                    #
# Autor : Pedro Henrique Sodré Puntel                                #
# Email : pedro.puntel@gmail.com                                     # 
# Tema : Modelagem Estatística de Redes Sociais com dados do Spotify #
# Script Encoding : UTF-8                                            #
######################################################################

##############
## Descrição :
##############
# Neste scrpit, faremos a modelagem estatística em cima das sócio-matrizes
# extraídas do Spotify. Como de costume, será feita primeiramente uma análise
# exploratória dos dados, no intuito de obter insights iniciais sobre a estrutura
# de relacionamentos implícita nos gostos musicais dos países.
#
# Em um segundo momento, serão aplicados aos dados os modelos de Escalonamento
# Multidimensional e Modelo de Distâncias Latentes, sendo este último utilizando 
# os métodos de estimação por Máxima Verossimilhança e estimação Bayesiana pelo
# método de MCMC. Demais detalhes sobre os modelos podem ser encontrados nas suas
# respectivas sub-seções.

##########
## Setup :
##########
# Pacotes utilizados
library("rio")
library("dplyr")
library("plotly")
library("igraph")
library("MCMCpack")
library("latentnet")

# Importação da sócio-matrix
top50_social_matrix = import(file.choose(), format="csv") %>% as.matrix()
unique_index = import(file.choose(), format="csv") %>% as.data.frame()
unique_index$Score = as.vector(unique_index$Score) %>% as.numeric()

#########################
## Análise Exploratória :
#########################
# Heatmap
plot_ly(z=top50_social_matrix, x=colnames(top50_social_matrix), y=colnames(top50_social_matrix),
        type="heatmap", zauto=F, zmin=0, zmax=50, colors=heat.colors(n=256)) %>%
  layout(title="Heatmap - Músicas compartilhadas no Top 50 entre os países",
         xaxis=list(categoryorder="array", categoryarray=colnames(top50_social_matrix)),
         yaxis=list(categoryorder="array", categoryarray=colnames(top50_social_matrix))) %>%
  colorbar(title="Músicas compartihadas")

# Barchart dos países quanto ao seu número de músicas compartilhadas
aux_vec = sort(colSums(top50_social_matrix), decreasing=T)-50 # desconsidera a diagonal da sociomatriz
plot_ly(as.data.frame(aux_vec), x=names(aux_vec), y=aux_vec, type="bar", marker=list(color="blue")) %>%
  layout(title="Países que mais compartilham músicas no seu Top 50 com os demais",
         xaxis=list(categoryorder="array", categoryarray=names(aux_vec)),
         yaxis=list(title="# Músicas compartilhadas", categoryorder="array", categoryarray=aux_vec)) %>%
  rangeslider()

# Barchart do Índice de autenticidade do gosto musical dos países
plot_ly(data=unique_index, x=unique_index$Country, y=unique_index$Score, type="bar", marker=list(color="green")) %>%
  layout(title="Índice de autenticidade do gosto musical dos países",
         xaxis=list(categoryorder="array", categoryarray=unique_index$Country),
         yaxis=list(title="Score de Autenticidade", categoryorder="array", categoryarray=unique_index$Score)) %>%
  rangeslider()

# Limpando memória
rm(aux_vec)

#############################################
## Modelo de Escalonamento Multidimensional :
#############################################

# O Escalonamento Multidimensional Clássico supõem existência de uma matriz de distâncias
# D[n x n], cujas entradas d_ij quantificam a distância entre o elemento i e o elemento j
# sob o contexto dos dados originais.
#
# O objetivo do método é então encontrar uma configuração de pontos no espaço p-dimensional,
# de tal forma que as coordenadas dos pontos, ao longo da dimensão p, produzam uma matriz de
# distâcias Euclidianas, cujos elementos estejam mais próximos possível dos de D.
#
# A justificativa para a implementação deste método é imediata :
# 
#       Gostaríamos de verificar a similiarede/dissimilarides dos países da rede
#          quanto ao seu padrão de conexões ou seja, quais seriam os países que
#                        mais assemelham-se em termos de perfil ? 
#
# Tão imediata quanto, é a sua interpretação :
#
#   Um grupo de dois ou mais países podem ser considerados semelhantes se a distância 
#      entre estes no espaço euclidano p-dimensianal (p <= n) é pequena e dissimilares
#                                  caso contrário.
#
# Em nossa aplicação do modelo, chegamos a conclusão de que o Escalonamento Multidimensional,
# quando aplicado para com uma matriz de adjacência valorada, acaba introduzindo um viés no
# sentido que indivíduos de uma rede cujo padrão de conexões é similar, acabam ficando muito
# próximos uns dos outros por mais que estes mesmos, entre si, não estabelecam nenhum
# tipo de relação. No nosso contexto, por exemplo, pode-se verificar este efeito com
# os países "não-populares" que ficam todos próximos uns dos outros (aglomerado),
# mesmo que entre eles, não existam conexões.
#
# A saída alternativa adotada foi realizar o Escalonamento Multidimensional em cima
# de uma "Sócio-Matriz de Correlações" cujas entradas i e j quantificam o Coeficiente
# de Correlação de Pearson entre os países i e j.

# Função que implementa o Escalonamento Multidimensional
mds_fit = function(socio_mat, use_risk_mat = F, use_cor_mat = T) {
  
  # Escalonamento Multidimensional Clássico
  if (isFALSE(use_risk_mat)) {
    
    if(isTRUE(use_cor_mat)) {
      
      # Computa o Escalonamento Multidimensional com a matriz de correlações
      mds_fit = dist(cor(socio_mat), method = "euclidean", diag = T, upper = T) %>% cmdscale(k = 2) %>% as.matrix()
      
    } else {
      
      # Computa o Escalonamento Multidimensional com a matriz original
      mds_fit = dist(socio_mat, method = "euclidean", diag = T, upper = T) %>% cmdscale(k = 2) %>% as.matrix()
      
    }
    
    # Nomeia as colunas do Data Table resultante para melhor manipulação
    colnames(mds_fit) = c("x.coord","y.coord")
    
    # Simplifica as arestas da socio-matriz para plotagem
    edges = graph.adjacency(as.matrix(round(socio_mat/10,0)), mode = "undirected", diag = F)
    
    # Grafo interativo
    tkplot(edges, vertex.size = 30, edge.width = 1, vertex.label = colnames(socio_mat),
           layout = mds_fit, vertex.color = "yellow", vertex.label.cex = 1, vertex.label.color = "black",
           edge.color = "grey")
    
    # Retorna os parâmetros
    list("par" = mds_fit) %>% return()
    
  }
  
  # Escalonamento Multidimensional com a Matriz de Riscos
  else {
    
    # Suponha que uma epidemia atingiu 4 cidades distintas e gostaríamos de quantificar
    # qual destas foi a mais tingidda/contaminda. É natural que de prontidão, a nossa
    # solução seja contar o número de contágios da doença que ocorreram em cada uma
    # das cidades e assim, verificar qual delas foi a mais impactada.
    #
    # Porém, existe um problema com essa abordagem: estamos desconsiderando a população
    # de cada cidade para efeito de comparação. Isto pode ser potencialmente problemático
    # uma vez que 50 casos de contágio em um população de 100 pessoas é muito mais preocupante
    # do que 500 casos em uma população de 100000 pessoas. 
    #
    # É justamente neste sentido que convêm realizarmos o Escalonamento Multidimensional não
    # diretamente ao nosso conjunto de dados mas sim com uma "Matriz de Riscos". Seu procediemnto
    # de montagem é descrito a seguir :
    #
    #  Seja A(nxn) a sociomatriz de interesse. De forma a entender proporcionalmente o quão
    #  forte ou fraca é a conexão entre dados dois atores da rede, montemos uma nova matriz
    #  P(nxn) onde cada entrada p_ij quantifica o quão acima/abaixo da média (ou do esperado)
    #  se dá a relação entre tais atores i e j. 
    #    
    #  Para tal, definemos uma taxa de normalização alpha, dada pela razão entre a soma das
    #  entradas da sociomatriz A, pela a soma das quantidades as os tais atores poderiam,
    #  supostamente, estabelecer uns com os outros.
    #    
    #  Assim, preenchemos por fim a matriz P onde cada entrada p_ij seja então a razão entre
    #  as entradas da sociomatriz (valor real) pelos valores esperados para tais conexões
    #  (valor teórico).
    #
    # Interpretação : "Quão acima da média se dá a relação entre os atores i e j da rede ?"
    
    # Sub-rotina responsável pela construção Matriz de Mínimos
    Build.Minimum_Matrix = function(socio_mat) {
      
      # Define a estrutura da Matriz de Mínimos associada
      min_mat = matrix(NA, ncol = ncol(socio_mat), nrow = nrow(socio_mat))
      
      # Diagonal da Sociomatriz
      m_diag = diag(socio_mat)
      
      # Constrói a Matriz de Mínimos
      for (i in 1:ncol(socio_mat)) {
        for (j in 1:nrow(socio_mat)) {
          min_mat[i,j] <- min(m_diag[i], m_diag[j])
        }
      }
      
      # Nomeia as linhas e as Colunas da Matriz de Mínimos
      colnames(min_mat) = colnames(socio_mat)
      rownames(min_mat) = rownames(socio_mat)
      
      # Retorna a Matriz de Mínimos
      return(min_mat)
      
    }
    
    # Computa a soma das entradas do triângulo inferior da sociomatriz
    sum.conex = socio_mat[lower.tri(socio_mat, diag = F)] %>% sum()
    
    # Constrói a Matriz de  Mínimos
    min_mat = Build.Minimum_Matrix(socio_mat)
    
    # Computa a soma das entradas do triângulo inferior da matriz de mínimos
    sum.min_conex = min_mat[lower.tri(min_mat, diag = F)] %>% sum()
    
    # Calculo da Taxa de Normalização
    # . Valor entre 0 e 1
    # . mais próximo de 1 indicam que estabeleceram-se mais conexões do que o esperado
    # . mais afastado de 1 indicam que estabeleceram-se menos conexões do que o esperado
    taxa.norm = sum.conex / sum.min_conex
    
    # Montagem da Matriz de Riscos
    risk.mat = taxa.norm * min_mat
    colnames(risk.mat) = colnames(socio_mat)
    rownames(risk.mat) = rownames(socio_mat)
    
    # Computa o Escalonamento Multidimensional
    mds_fit = dist(risk.mat, method = "euclidean", diag = F, upper = T) %>% cmdscale(k = 2) %>% as.matrix()
    
    # Nomeia as colunas do Data Table resultate para melhor manipulação
    colnames(mds_fit) = c("x.coord","y.coord")
    
    # Simplifica as arestas da socio-matriz para plotagem
    edges = graph.adjacency(as.matrix(round(socio_mat/10,0)), mode = "undirected", diag = F)
    
    # Grafo interativo
    tkplot(edges, vertex.size = 30, edge.width = 1, vertex.label = colnames(socio_mat), layout = mds_fit,
           vertex.color = "yellow", vertex.label.cex = 1, vertex.label.color = "black", edge.color = "grey")
    
    # Retorna os parâmetros
    list("par" = mds_fit) %>% return()
    
  }
  
}
mds_fit = compiler::cmpfun(mds_fit)

# Resultados
spotify_mds = mds_fit(top50_social_matrix, use_risk_mat=F, use_cor_mat=T)
spotify_mds = mds_fit(top50_social_matrix, use_risk_mat=F, use_cor_mat=F)

#########################################
## Estimação por Máxima Verossimilhança :
#########################################
# Os modelos de Espaços Latentes propõem uma classe de modelos nos quais a probabilidade
# de uma relação entre os atores de uma rede é uma função não-linear das suas respectivas
# posições em um espaço intitulado "espaço social". Dentre as possíves interpretações para
# o referido espaço, existe aquela de um espaço não-observacional que retrata o padrão de
# conectividade dos atores da rede como base nas suas posições. Intuitivamente, atores mais
# próximos são mais próprios de estabelecerem uma relação do que aqueles mais afastados.
# 
# Estatisticamente, o modelo caracteriza-se por assumir que cada vértice possui uma posição
# desconhecida no espaço social, e que uma relação entre dois vértices são condicionalmente
# independentes dadas estas suas respectivas posições. Assim, a probabilidade de uma conexão
# entre estes indivíduos é modelada no seguinte formato: 
#
#                                  logit(pij) ~ Q - |Zi-Zj|
#
#               * Q = nível de atividade social de cada ator (covariável) 
#               * |Zi-Zj| = distância entre os vertices i e j no espaço social
#
# 'Q' pode ser pensado também como o grau de extroversão do i-ésimo ator da rede, de forma
# que a distância entre um ator aos demais seja "ponderada" também pelo quão propício tal ator
# é de estabelecer uma conexão.O
#
# Cada 'Zi' e 'Zj' são vetores de dimensão 2 (coordenadas no plano bi-dimensional).
#
# As distâncias |Zi-Zj| entre os atores são invariantes quanto a reflexão, rotação ou deslo-
# camento dos eixos coordenados. Como consequência, é comum fixar 3 dos parâmetros para fins de
# estimação. Valores iniciais candidatos às posições dos atores no espaço latente são aqueles provenientes
# do escalonamento multidimensional.

# Rotina que implementa o Modelo de Distâncias Latentes
ldm_fit = function(socio_mat, Is.direct = F, Is.Weighted = T, Has.Loops = F, to.fit = c("ergmm","bayesian")) {
  
  # Sub-rotina que implementa o modelo com estimação dos parâmetros por Máxima Verossimilhança
  ergmm_mle = function(socio_mat, Is.direct, Is.Weighted, Has.Loops) {
    
    # Objeto tipo network à ser manipulado pela função ergmm()
    net_aux = network::network(socio_mat, directed = Is.direct, loops = Has.Loops, matrix.type = "adjacency")
    
    # Atribui pesos as arestas da rede (caso esta seja de fato valorada)
    if (isTRUE(Is.Weighted)) {
      network::set.edge.value(net_aux, attrname = "weight", value = socio_mat)
    }
    
    # Valores iniciais provenientes do Escalonamento Multidimensional
    par_vec = dist(socio_mat, method = "euclidean", diag = F, upper = T) %>%  cmdscale(k = 2) %>% as.vector()
    par_vec[(length(par_vec)+1)] = 0
    
    # Parâmetros estimados
    ergmm_coords = ergmm(net_aux ~ euclidean(d = 2), verbose = 1, family = "Poisson", user.start = list(par_vec), tofit = "mle")
    
    # Simplifica as arestas da socio-matriz para plotagem, mantendo somente significativas (> 10 meses no Top 20)
    edges = graph.adjacency(as.matrix(round(socio_mat/10,0)), mode = "undirected", diag = F)
    
    # Grafo interativo
    tkplot(edges, vertex.size = 30, edge.width = 1, vertex.label = colnames(socio_mat), layout = ergmm_coords$mle$Z,
           vertex.color = "Green", vertex.label.cex = 1, vertex.label.color = "black", edge.color = "grey", axes = T)
    
    # Retorna as informações
    return(list("output" = ergmm_coords, "par" = ergmm_coords$mle$Z))
    
  }
  
  # Sub-rotina que implementa o modelo com estimação dos parâmetros por simulação estocástica bayesiana
  bayesian_estimation = function(socio_mat, Is.direct, Is.Weighted, Has.Loops) {
    
    # Método de simulação estocástica cujo objetivo é obter amostras da distribuição
    # "a posteriori" quando as mesmas não possuem um forma fechada, ou seja, quando não
    # pertencem à nenhuma distribuição conhecida. De forma simples, o método busca estimar
    # a função de distribuição "a posteriori" possibilitando assim inferência acerca dos
    # parâmetros de interesse.
    #
    # O objetivo do método de MCMC é construir uma cadeia de Markov cuja distribuição limite
    # seja igual a distribuição de interesse. O algortimo realiza um número finito de 
    # sucessivas de amostagens desta cadeia, de forma a esperar-se que a sua distribuição
    # convirja para a tal distribuição "a posteriori". 
    #
    # Dentre os métodos mais utlizados para a construção desta cadeia de Markov, destacam-se
    # o "Amostrador de Gibbs" e o método de "Metropolis-Hastings" :
    #
    # O método Amostrador de Gibbs é um método iterativo de amostragem de uma cadeia de Markov,
    # cujas probabilidades de transição são definidas de acordo com a distribuição condicional
    # completa dos parâmetros | conjunto de dados observados.
    #
    # O método Metropolis-Hastings é um método que extrai amostras de qualquer distribuição de
    # probabilidade P(x), dado o valor de uma função f(x) que seja proporcional é densidade P.
    # O algortimo funciona de forma iterativa, amostrando da função proporcional f(x), de forma
    # que o próximo valor amostrado seja exclusivamente dependente do seu anterior, formando 
    # assim uma sequência de amostras que concomitam em uma cadeia de Markov. A diferença,
    # porém, reside em uma espécie de "filtro" que avalia a plausabilidade dos valores
    # amostrados, de forma que onde em cada iteração, o algortimo aceita/rejeita determinada
    # amostra com certa probabilidade. Caso aceita, seu valor será utilizado na próxima
    # iteração e caso contrário, descarta-se este mantendo aquele da iteração anterior.
    #
    # Assim, os processos de estimação descritos acima serão aplicados no R com o auxílio dos
    # pacotes 'network', 'latentnet' e 'MCMCpack' como segue :
    
    # Objeto tipo 'network' a ser manipulado pela função ergmm()
    net_aux = network::network(socio_mat, directed = Is.direct, loops = Has.Loops, matrix.type = "adjacency")
    
    # Atribui pesos as arestas da rede (caso esta seja de fato valorada)
    if (isTRUE(Is.Weighted)) {
      network::set.edge.value(net_aux, attrname = "weight", value = socio_mat)
    }
    
    # Exponential Random Graph Mixed Models (ERGMM's)
    # 
    # . pmode = Posterior mode (MAP)
    # . mcmc = Monte Carlo Markov Chain
    # . mkl = Kullback-Leiber weighted likelihood
    # . mle = Maximum Likelihood Estimation
    # . procrustes = Procrustean Analysis
    #
    # Dentre os ERGMM's aplicados, nos interessa somente as estimações referente é moda
    # a posteriori, Máxima Verossimilhança e Procrustes. Uma breve explicação destes
    # modelos é feita a seguir (exceto sobre a Máxima Verossimilhança que já sabemos) :
    #
    # Uma estimativa de Máxima Densidade a Posteriori (MAP) é uma estimativa de uma
    # quantidade desconhecida, que é igual a moda da distribuição posterior. O MAP pode ser
    # usado para obter uma estimativa pontual de uma quantidade não observada com base em
    # dados empíricos. Esta é relacionada ao método de estimação por máxima verossimilhança,
    # mas emprega um objetivo de otimização que incorpora uma distribuição a priori (esta
    # fornece informações adicionais através de um evento prévio) sobre a quantidade que se
    # deseja estimar.
    
    # Procrustes Analysis
    #
    # Distances between a set of points in Euclidean space are invariant under rotation,
    # reflection, and translation. Therefore, for each matrix of latent positions Z, there
    # is an infinite number of other positions giving the same log-likelihood.
    #
    # A confidence region that includes two equivalent positions Z1 and Z2 is, in a sense,
    # overestimating the variability in the unknown latent positions, because these are
    # identical for Z1 and Z2.
    #
    # Fortunately, this problem can be resolved by perfoming inference on a set that we call 
    # "equivalence classes of latent positions" : Let [Z] be the set of positions equivalent
    # to matrix of latent positions Z under the previous mentioned operations. For each [Z],
    # there is one set of distances between the nodes. We oftne refet to this class of positions
    # as a "configuration model".
    # 
    # Therefore, We make inference on these configurations by performing inference on particular
    # elements that are comparable across configurations. For example, given an configuration set
    # [Z], we select for inference the configuration Z* that minimizes the sum of squared positional
    # difference. It turns out that Z* ends up beign a "Procrustean" transformation of Z, therefore
    # being the element of the set [Z] which is closest to the the original "fixed" set.
    #
    # We typically take for this orginal "fixed" set the MLE estimates of the latent positions.
    
    # Valores iniciais provenientes do Escalonamento Multidimensional
    par_vec = dist(socio_mat, method = "euclidean", diag = F, upper = T) %>% cmdscale(k = 2) %>% as.vector()
    par_vec[(length(par_vec)+1)] = 0
    
    # Fits the Model
    model_fit = ergmm(net_aux ~ euclidean(d = 2), family = "Poisson", user.start = list(par_vec), tofit = c("pmode"))
    
    # Simplifica as arestas da socio-matriz para plotagem, mantendo somente significativas (> 10 meses no Top 20)
    edges = graph.adjacency(as.matrix(round(socio_mat/10,0)), mode = "undirected", diag = F)
    
    # Grafo interativo
    tkplot(edges, vertex.size = 30, edge.width = 1, vertex.label = colnames(socio_mat), layout = model_fit$Results$pmode$Z,
           vertex.color = "red", vertex.label.cex = 1, vertex.label.color = "black", edge.color = "grey")
    
    # Retorna as informações
    return(list("info" = model_fit, "par" = model_fit$pmode$Z))
    
  }
  
  if (to.fit == "ergmm") {
    return(list("Ergmm.Mle"=ergmm_mle(socio_mat = socio_mat,
                                     Is.direct = Is.direct,
                                     Is.Weighted = Is.Weighted,
                                     Has.Loops = Has.Loops)))
    
  } else {
    return(list("Bayesian.Fit"=bayesian_estimation(socio_mat = socio_mat,
                                                  Is.direct = Is.direct,
                                                  Is.Weighted = Is.Weighted,
                                                  Has.Loops = Has.Loops)))
    
  }
  
}
ldm_fit = compiler::cmpfun(ldm_fit)

# Resultados
spotify_ldm_mle = ldm_fit(top50_social_matrix,F,T,F,"ergmm")
spotify_ldm_bay = ldm_fit(top50_social_matrix,F,T,F,"bayesian")
