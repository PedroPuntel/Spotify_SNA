#########################################################################
# Data de início :  21/01/2020                                          #
# Última atualização : 23/01/2020                                       #
# Autor : Pedro Henrique Sodré Puntel                                   #
# Email : pedro.puntel@gmail.com                                        #
# Tema : Spotify - Estatísticas Descritivas (Análise de Redes Sociais)  #
# Script Encoding : UTF-8                                               #
#########################################################################

################
## Referências :
################
# Network Science
# > https://www.sci.unich.it/~francesc/teaching/network/
#
# Introduction to Social Network Methods 
# > http://faculty.ucr.edu/~hanneman/nettext/
#
# Social Network Analysis : A methodological introduction 
# > http://courses.washington.edu/ir2010/readings/butts.pdf
#
# Medidas de Centralidade em Grafos 
# > http://objdig.ufrj.br/60/teses/coppe_m/LeandroQuintanilhaDeFreitas.pdf
# 
# Comunity Analysis Algorithms in R using the igraph package
# > https://stackoverflow.com/questions/9471906/what-are-the-differences-between-community-detection-algorithms-in-igraph#
#
# R Tutorial: How to identify communities of items in networks :
# > https://psych-networks.com/r-tutorial-identify-communities-items-networks/
# 
# Introduction to R and Network Analysis :
# > http://kateto.net/ruworkshop
#
# Polnet Network Visualization with R :
# > http://kateto.net/network-visualization

##############
## Descrição :
##############
# Neste script, serão aplicados alguns conceitos estudados sobre a modelagem 
# matemática/estatística de redes sociais.
#
# De forma geral, as análises a seguir serão conduzidas estudando-se propriedades
# espectrais da sócio-matriz que comporta os países da chart Top 50 by Country do
# Spotify. É com base nestas propriedades que derivam-se diversas medidas capazes
# de fornecer insights valiosos sobre os relacionamentos entre os vértices de uma
# rede.

##########
## Setup :
##########
# Pacotes utlizados
library("dplyr")
library("igraph")

# Opções adicionais (desabilita a notação científica padrão do R e configura o número de casas decimais para 4)
options(scipen = 999, digits = 4)

# Importa a sócio-matriz
social_matrix = rio::import(file=file.choose(), format="csv") %>% as.matrix()
rownames(social_matrix) = colnames(social_matrix)

# Grafo que modelará a rede
spotify_graph = graph_from_adjacency_matrix(social_matrix, mode = "undirected", diag = F, weighted = T)

# Função auxiliar que mapeia os valores de um vetor em um determinado intervalo
mapToRange = function(x, from, to) {return( (x-min(x)) / max(x-min(x)) * (to-from) + from)} 

# Tamanho dos vértices da rede proporcional ao Score de Autenticidade
unique_index = rio::import(file=file.choose(), format="csv") %>% as.data.frame()
unique_index = unique_index[order(unique_index$Country),]
V(spotify_graph)$size = mapToRange(unique_index$Score,15,35)

###############################
## Graph Level Indices (GLIs) :
###############################
# Rotina que computa características globais do Grafo
graph_glis = function(g, use.loops, has.directions, has.weights) {
  
  # Sub-rotina que retorna características básicas do grafo
  graph_stats = function(g) {
    
    # Número de arestas
    edge_count = ecount(g)
    
    # Número de vértices
    vertex_count = vcount(g)
    
    # Densidade do Grafo 
    g_density = graph.density(g, loops = use.loops)
    
    # Caminho Geodésico médio
    g_mean_geodesic = mean_distance(g, directed = has.directions, unconnected = isFALSE(is.connected(g, mode = "strong")))
    
    # Maior caminho geodésico
    g_big_geodesic = farthest_vertices(g, directed = has.directions, unconnected = isFALSE(is.connected(g, mode = "strong")))
    
    # Diâmetro do Grafo
    g_diameter = diameter(g, directed = has.directions, unconnected = is.connected(g, mode = "strong"))
    
    # Transitividade do Grafo
    # . Quão interligadas são as conexões entre os vértices do grafo. Na literatura, para grafos não-direcionados,
    # a transitividade se dá pela razão entre  o número de loops de tamanho 3 (x -- > y --> z -- > x) e o total
    # de arestas do grafo.
    if (isTRUE(has.directions)) {
      g_trans = transitivity(g, type = "directed", isolates = NaN)
    } else {
      g_trans = transitivity(g, type = "undirected", isolates = NaN)
    }
    
    # Pareamento Assortativo/Dissortativo
    #
    #  Um fenômeno comum em redes sociais é a tendência natural dos atores em associar-se
    # com aqueles semelhantes a si mesmos de alguma forma. Tal tendência é chamada de 
    # homofilia ou mistura combinatória. É possível ainda, porém de forma mais rara, também
    # encontramos uma mistura desassortativa, ou seja, a tendência dos se associarem com
    # outros diferentes de si.
    #
    #  O pareamento assortativo de atores em uma rede pode ser identificado por meio de 
    # características enumerativas ou escalares. Uma característica enumerativa tem um
    # conjunto finito de valores possíveis. Exemplos são gênero, raça e nacionalidade.
    # Dada uma característica enumerativa, cada nó da rede pode ser atribuído a um tipo de
    # acordo com o valor da característica para o nó. Assim, dizemos que uma rede é então
    # assortativa se uma fração significativa das arestas ligam entre vértices do mesmo tipo.
    #
    #  Matematicamente, uma medidade de assortatividade é definida pelo número de arestas
    # que ocorrem entre vértices do mesmo tipo diminúidas daquelas as quais esperaríamos 
    # encontrar em uma disposição aleatória dos vértices na rede. Indo além, pode-se provar 
    # que a assortatividade de um vértice é dada pela sua covariância com os demais pares e
    # que ao ser 'normalizada' coincide com o coeficiente de correalção de Pearson.
    g_assortativity = assortativity_degree(g, directed = has.directions)
    
    # Tabela com as informações
    out_df = c(edge_count, vertex_count, g_density, g_mean_geodesic, g_big_geodesic, g_diameter, g_trans, g_assortativity) %>%
      matrix(nrow=1, ncol=8, byrow = F) 
    colnames(out_df) = c("Edge Count","Vertex Count","Density","Avg. Geodesic","Biggest Geodesic","Diameter","Transitivity","Assortativity")
    
    # Retorna as informações
    return(out_df)
    
  }
  
  ## Retorna os resultados
  list("Graph Stats" = graph_stats(g)) %>% return()
  
}
graph_glis = compiler::cmpfun(graph_glis)

# Análise das informações
g_glis = graph_glis(spotify_graph,F,F,T)
g_glis$`Graph Stats`

##############################
## Node Level Indices (NLIs) :
##############################
# Rotina que computa diversas medidas associadas a cada vértice do grafo
graph_nlis = function(g, use.loops, has.directions, has.weights) {
  
  # Centralidade de Grau dos vértices do grafo
  # . Número de vizinhos/vértices adjacentes de cada vértice do grafo
  node_degree = degree(g, loops = use.loops, mode = "all", normalized = T)
  
  # Centralidade do autovetor dos vértices do grafo
  #
  # . Nem todos os vértices de uma rede são equivalentes, alguns são mais importantes
  # do que outros. Sendo assim, um vértice que estabelece algum tipo de relação com 
  # um outro vértice "importante" da rede merece atenção se comparada com as demais.
  #
  # . "A node receiving many links does not necessarily have a high eigenvector centrality
  # (it might be that all linkers have low or null eigenvector centrality). Moreover, a
  # node with high eigenvector centrality is not necessarily highly linked (the node might
  # have few but important linkers)."
  node_egvc = eigen_centrality(g, directed = has.directions, scale = T)
  
  # Centralidade de Page Rank dos vértices do grafo
  #
  # . Um problema com a medida de centralidade de autovetor é que se um vértice com alta 
  # centralidade referencia muitos outros, então estes últimos por sua vez também terão
  # um alta centralidade. Porém, é razoável pensar que vértices referenciados por outros
  # altamente centrais são de menor importância: o ganho de centralidade do referenciado
  # deve ser menosprezado se quem o referencia, referencia também vários outros. Assim,
  # a cálculo da centralidade Page Rank leva em consieração o número de arestas indicentes
  # ao vértice, a propensão daquele que referencia, e por fim o valor da sua centralidade.
  #
  # . Variados estudos sugerem o valor 0.85 para o parâmetro 'damping' do Page Rank.
  node_pgrk = page_rank(g, algo = "power", directed = has.directions, damping = 0.85)
  
  # Centralidade de Kleingberg dos vértices do grafo
  #
  # . É razoável assumir que um vértice também é importante por referenciar outros
  # vértices importantes da rede. Em redes de co-autorias, por exemplo, artigos podem
  # referenciar vários outros os quais julgam sere 'autoridades' no assunto e por isso,
  # podem ser vistos como importantes. Seguindo essa lógica, é natural dividirmos os atores
  # de uma rede em duas categorias: os 'chefes', que são referenciados por muitos outros,
  # e aqueles 'subordinados', que apontam para os 'chefes'.
  node_klgbr = authority.score(g, scale = T)
  
  # Centralidade de Proximidade dos vértices do grafo
  #
  # . Computa a distância média que um vértice possui dos demais, retornando valores menores
  # para os vértices que estão separados por uma curta distância geodésica. Vértices com
  # essa característica possuem um melhor acesso a informação dissipada na rede. Em redes
  # sociais, atores com menores índices de centralidade de proximidade são mais propícios
  # a influenciar os demais vértices.
  node_clo = closeness(g, mode = "all", normalized = T)
  
  # Centralidade de Intermediação dos vértices do grafo
  #
  # . Quantifica o quão provável um dado vértice é de intermediar a ligação entre
  # outros dois vértices. Nesse sentido, vértices com alto índice de entrelaçamento
  # influenciam diretamente no fluxo de informação pela rede. Analagomente, estes são
  # aqueles cuja remoção da rede irá pertubar a comunição inter-vértices.
  node_btw = betweenness(g, directed = has.directions, normalized = T)
  
  # Tabela com as medidas de centralidade de cada vértice da rede
  Degree = node_degree %>% as.list() %>% unlist() %>% as.numeric()
  EigenVec = node_egvc$vector %>% as.list() %>% unlist() %>% as.numeric()
  PageRank = node_pgrk$vector %>% as.list() %>% unlist() %>% as.numeric()
  Kleinberg = node_klgbr$vector %>% as.list() %>% unlist() %>% as.numeric()
  Closeness = node_clo %>% as.list() %>% unlist() %>% as.numeric()
  Betweenness = node_btw %>% as.list() %>% unlist() %>% as.numeric()
  Centrality_Table = cbind(Degree, EigenVec, PageRank, Kleinberg, Closeness, Betweenness) %>% as.matrix()
  rownames(Centrality_Table) = V(g)$name
  
  ## Retorna a tabela
  list("Node_Centrality_Table" = Centrality_Table) %>% return()
  
}
graph_nlis = cmpfun(graph_nlis)

# Resultados
spotify_nlis = graph_nlis(spotify_graph, T, F, T)
spotify_centr_table = spotify_nlis$Node_Centrality_Table

################################
## Comunity Detection Analysis :
################################

# Rotina que implementa diversos algoritmos para detecção de comunidades
graph_cda = function(g, use.loops, has.directions, has.weights) {
  
  # Maximização de Modularidade
  # 
  # . A identificação de comunidades de atores em uma rede é feita com base na
  # maximização de uma medida determinada "modularidade". Tal medida computa a
  # diferença entre o número de arestas que ocorrem entre vértices de uma mesma
  # comunidade diminuída e o número esperado de arestas entre os mesmos no caso
  # de uma disposição aleatória dos vértices pela rede. Nesse sentido, tal medida
  # assume valores positivos se existem mais arestas conectando vértices "de um
  # mesmo tipo" do que o esperado e valores negativos caso contrário.
  
  # Método Fast-Greedy
  # 
  # . Método de abordagem hierárquica, aglomerativo, que atua diretamente na otimização
  # da função modularidade associada a rede. Por ser aglomerativo, o algoritmo inicia
  # com cada nó em uma comunidade separada, e mescla iterativamente as mesmas. Cada fusão
  # é sempre localmente ótima e o algoritmo termina quando não é mais possível aumentar a
  # modularidade.
  #
  # . O método tem como vantagem a sua velocidade sendo geralmente tentado como uma primeira
  # aproximação por não possuir parâmetros ajustáveis (como no K-Means, por exemplo). No entanto,
  # sabe-se o método sofre de um limite de resolução, isto é, comunidades abaixo de um determinado
  # limite (threshold) de tamanho serão sempre mescladas com comunidades vizinhas.
  fast_greedy = cluster_fast_greedy(g, merges = T, membership = T)
  V(g)$fast_greedy = fast_greedy$membership
  
  # Método do Autovetor Líder
  #
  # . Método de abordagem hierárquica, divisivo, que otimiza a função de modularidade com base
  # na avaliação do autovetor associado ao maior autovalor da chamada matriz de modularidade.
  # Por envolver cálculos de autovetored, tal método pode não funcionar em gráficos degenerados.
  #
  # . A vantagem deste algortimo é que para grafos não degenerados, este atinge um score de
  # modularidade geralmente maior do que outros métodos (como o Fast-Greedy, pode exemplo),
  # porém ao custo de ser mais lento. O algortimo agrupa os vértices da rede baseando-se nos
  # sinais das entradas do autovetor calculado, de forma que caso todas as entradas do autovetor
  # sejam de mesmo sinal, a rede possui estruturação única.
  egv_comunnity = cluster_leading_eigen(g, steps = 10, callback = NULL)
  V(g)$egv_comunnity = egv_comunnity$membership
  
  # Método da Caminhanda Aleatória
  #
  # .Tenta encontrar comunidades em um grafo com base no sorteio aleatório de camminhos
  # entre os vértices. A premissa do método é que caminhadas curtas tendem a permanecer
  # na mesma comunidade.
  rdm_walk = cluster_walktrap(g, steps = 10, modularity = T, merges = T)
  V(g)$rdm_walk = rdm_walk$membership
  
  # Retorna os grafos associados a cada algoritmo de clusterização
  g_layout = layout_randomly(g, dim = 2)
  
  n_clusters = as.factor(V(g)$fast_greedy) %>% levels() %>% as.numeric() %>% max()
  g_colors = RColorBrewer::brewer.pal(n_clusters, name = "Set1")
  plot.igraph(g, main = "Comunity Analysis - Fast Greedy", vertex.shape = "circle", 
              vertex.color = g_colors[V(g)$fast_greedy], vertex.size = 20,
              vertex.label = V(g)$name, vertex.label.font = 2, vertex.label.cex = 0.6,
              vertex.label.color = "Black", edge.color = "Grey", edge.width = 1.25, edge.lty = 1,
              layout = g_layout)
  
  n_clusters = as.factor(V(g)$egv_comunnity) %>% levels() %>% as.numeric() %>% max()
  g_colors = RColorBrewer::brewer.pal(n_clusters, name = "Set1")
  plot.igraph(g, main = "Comunity Analysis - Leading Eigen Vector", vertex.shape = "circle", 
              vertex.color = g_colors[V(g)$egv_comunnity], vertex.size = 20,
              vertex.label = V(g)$name, vertex.label.font = 2, vertex.label.cex = 0.6,
              vertex.label.color = "Black", edge.color = "Grey", edge.width = 1.25, edge.lty = 1,
              layout = g_layout)
  
  n_clusters = as.factor(V(g)$rdm_walk) %>% levels() %>% as.numeric() %>% max()
  g_colors = RColorBrewer::brewer.pal(n_clusters, name = "Set1")
  plot.igraph(g, main = "Comunity Analysis - Random Walk", vertex.shape = "circle", 
              vertex.color = g_colors[V(g)$rdm_walk], vertex.size = 20,
              vertex.label = V(g)$name, vertex.label.font = 2, vertex.label.cex = 0.6,
              vertex.label.color = "Black", edge.color = "Grey", edge.width = 1.25, edge.lty = 1,
              layout = g_layout)
  
  # Tabela com as informações dos algoritmos
  Membership_Table = cbind(V(g)$name, fast_greedy$membership, egv_comunnity$membership, rdm_walk$membership) %>%
    as.data.frame()
  colnames(Membership_Table) = c("Vertex Label", "Fast-Greedy", "Leading Eigen Vector", "Random Walk")
  
  # Retorna as informações
  list("Fast_Greedy" = fast_greedy,
       "Lead_Eigen" = egv_comunnity,
       "Random_Walk" = rdm_walk,
       "Membership_Table" = Membership_Table) %>% return()
  
}
graph_cda = compiler::cmpfun(graph_cda)

# Resultados
spotify_cda = graph_cda(spotify_graph, F, F, T)

##################
## Visualizações :
##################

# Alguns algoritmos do pacote igraph para disposição dos vértices na rede 
#
# . Davidson-Harel Layout
# > https://igraph.org/r/doc/layout_with_dh.html
# 
# . Fruchterman-Reingold
# > https://igraph.org/r/doc/layout_with_fr.html
#
# . Graph-Opt Layout
# > https://igraph.org/r/doc/layout_with_graphopt.html
# 
# . Kamada-Kawai Layout
# > https://igraph.org/r/doc/layout_with_kk.html
#
# . Large-Graph Layout
# > https://igraph.org/r/doc/layout_with_lgl.html
#
# . Multidimensional Scaling Layout
# > https://igraph.org/r/doc/layout_with_mds.html

# Grafos estáticos
par(mfcol=c(1,2))
dh_layout = layout_with_dh(spotify_graph)
plot.igraph(spotify_graph, main = "Davidson-Harrel Layout", vertex.size = V(spotify_graph)$size, vertex.color = "dark orange",
            vertex.frame.color = "black", vertex.shape = "circle", vertex.label = NA,
            vertex.label.font = 2, vertex.label.cex = 0.6, vertex.label.dist = 0,  vertex.label.color = "black",
            edge.color = "darkgrey", edge.width = 0.5, edge.lty = 1, edge.label = NA,  layout = dh_layout)
plot.igraph(spotify_graph, main = "Davidson-Harrel Layout", vertex.shape = "none",
            vertex.label = V(spotify_graph)$name, vertex.label.font = 2, vertex.label.cex = 0.6,
            vertex.label.color = "black", edge.color = "white", layout = dh_layout)

par(mfcol=c(1,2))
mds_layout = layout_with_mds(spotify_graph)
plot.igraph(spotify_graph, main = "Multidimensional Scaling Layout", vertex.size = V(spotify_graph)$size, vertex.color = "yellow",
            vertex.frame.color = "black", vertex.shape = "circle", vertex.label = NA,
            vertex.label.font = 2, vertex.label.cex = 0.6, vertex.label.dist = 0,  vertex.label.color = "black",
            edge.color = "darkgrey", edge.width = 0.5, edge.lty = 1, edge.label = NA,  layout = mds_layout)
plot.igraph(spotify_graph, main = "Multidimensional Scaling Layout", vertex.shape = "none",
            vertex.label.font = 2, vertex.label.cex = 0.6, vertex.label = V(spotify_graph)$name,
            vertex.label.dist = 0, vertex.label.color = "black", edge.color = "white", layout = mds_layout)

par(mfcol=c(1,2))
fr_layout = layout_with_fr(spotify_graph, dim = 2)
plot.igraph(spotify_graph, main = "Fruchterman-Reingold Layout", vertex.size = V(spotify_graph)$size, vertex.color = "green",
            vertex.frame.color = "black", vertex.shape = "circle", vertex.label = NA,
            vertex.label.font = 2, vertex.label.cex = 0.6, vertex.label.dist = 0,  vertex.label.color = "black",
            edge.color = "darkgrey", edge.width = 0.5, edge.lty = 1, edge.label = NA,  layout = fr_layout)
plot.igraph(spotify_graph, main = "Fruchterman-Reingold Layout", vertex.shape = "none",
            vertex.label = V(spotify_graph)$name, vertex.label.font = 2, vertex.label.cex = 0.6,
            vertex.label.dist = 0,  vertex.label.color = "black", edge.color = "white", layout = fr_layout)

par(mfcol=c(1,2))
kk_layout = layout_with_kk(spotify_graph, dim = 2)
plot.igraph(spotify_graph, main = "Kamada-Kawai Layout", vertex.size = V(spotify_graph)$size, vertex.color = "red",
            vertex.frame.color = "black", vertex.shape = "circle", vertex.label = NA,
            vertex.label.font = 2, vertex.label.cex = 0.6, vertex.label.dist = 0,  vertex.label.color = "black",
            edge.color = "darkgrey", edge.width = 0.5, edge.lty = 1, edge.label = NA,  layout = kk_layout)
plot.igraph(spotify_graph, main = "Kamada-Kawai Layout", vertex.shape = "none",
            vertex.label = V(spotify_graph)$name, vertex.label.font = 2, vertex.label.cex = 0.6,
            vertex.label.dist = 0,  vertex.label.color = "black", edge.color = "white", layout = kk_layout)

par(mfcol=c(1,2))
lgl_layout = layout_with_lgl(spotify_graph)
plot.igraph(spotify_graph, main = "Large-Graph Layout", vertex.size = V(spotify_graph)$size, vertex.color = "blue",
            vertex.frame.color = "black", vertex.shape = "circle", vertex.label = NA,
            vertex.label.font = 2, vertex.label.cex = 0.6, vertex.label.dist = 0,  vertex.label.color = "black",
            edge.color = "darkgrey", edge.width = 0.5, edge.lty = 1, edge.label = NA,  layout = lgl_layout)
plot.igraph(spotify_graph, main = "Large-Graph Layout", vertex.shape = "none",
            vertex.label = V(spotify_graph)$name, vertex.label.font = 2, vertex.label.cex = 0.6,
            vertex.label.dist = 0,  vertex.label.color = "black", edge.color = "white", layout = lgl_layout)

par(mfcol=c(1,2))
opt_layout = layout_with_graphopt(spotify_graph)
plot.igraph(spotify_graph, main = "Graph-Opt Layout", vertex.size = V(spotify_graph)$size, vertex.color = "purple",
            vertex.frame.color = "black", vertex.shape = "circle", vertex.label = NA,
            vertex.label.font = 2, vertex.label.cex = 0.6, vertex.label.dist = 0,  vertex.label.color = "black",
            edge.color = "darkgrey", edge.width = 0.5, edge.lty = 1, edge.label = NA,  layout = opt_layout)
plot.igraph(spotify_graph, main = "Graph-Opt Layout", vertex.shape = "none",
            vertex.label = V(spotify_graph)$name, vertex.label.font = 2, vertex.label.cex = 0.6,
            vertex.label.dist = 0,  vertex.label.color = "black", edge.color = "white", layout = opt_layout)

