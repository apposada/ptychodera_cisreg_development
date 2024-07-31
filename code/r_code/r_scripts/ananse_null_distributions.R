milestones = sort(sapply(seq(1:9, by = .5),function(x){x*c(10,100,1000,10000)}))
gene_pool_1 <- V(g_LG)$name[V(g_LG)$name %in% pfla_attributes_list[[2]]$id]
gene_pool_2 <- V(g_LG)$name[!(V(g_LG)$name %in% pfla_attributes_list[[2]]$id)]
lg_30 <- list()
j = 0
set.seed(1234)
for(i in 1:100000){
  vids <- c(sample(gene_pool_1,15),sample(gene_pool_2,15))
  lg_30[[i]] <- induced_subgraph(
    g_LG,
    vids = V(g_LG)[V(g_LG)$name %in% vids],
    impl = "auto"
  )
  if(i - j == 100){
    message(i)
    j = j+100
  }
}
ananse_null_distr_30genes <-
  sapply(
    lg_30,
    function(x){
      length(which(igraph::degree(x, mode = "all")==0)) / vcount(x)
    }
  )
ananse_null_distr_30genes <- ananse_null_distr_30genes[complete.cases(ananse_null_distr_30genes)]

plot(density(ananse_null_distr_30genes))
length(ananse_null_distr_30genes[ananse_null_distr_30genes==0])/length(ananse_null_distr_30genes)

lg_50 <- list()
j = 1
set.seed(1234)
for(i in 1:100000){
  vids <- c(sample(gene_pool_1,25),sample(gene_pool_2,25))
  lg_50[[i]] <- induced_subgraph(
    g_LG,
    vids = V(g_LG)[V(g_LG)$name %in% vids],
    impl = "auto"
  )
  if(i - j == 100){
    message(i)
    j = j+100
  }
}
ananse_null_distr_50genes <-
  sapply(
    lg_50,
    function(x){
      length(which(igraph::degree(x)==0)) / vcount(x)
    }
  )
ananse_null_distr_50genes <- ananse_null_distr_50genes[complete.cases(ananse_null_distr_50genes)]

plot(density(ananse_null_distr_50genes))
length(ananse_null_distr_50genes[ananse_null_distr_50genes==0])/length(ananse_null_distr_50genes)

save(ananse_null_distr_50genes, ananse_null_distr_30genes, file = "outputs/rda/ananse_null_distributions.rda")