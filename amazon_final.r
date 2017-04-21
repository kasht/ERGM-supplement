#setwd("C:/Users/Sergii/Dropbox/MTSSC/new/code/livejournal")
setwd("/home/svb1/local_TERGM/livejournal_clean")
library(hergm)
library(parallel)

set.seed(123456)

#number of cores, set equal to 1 for Windows
#Run with mc.cores = 45 for the paper
mc.cores = 45;

#Preprocessing

# Read in the data
x <- scan("com-amazon.top5000.cmty.txt", what="", sep="\n")
# Separate elements by one or more whitepace
y <- strsplit(x, "[[:space:]]+")

group_sizes <- unlist(lapply(y, length))


#Sampling nodes
set.seed(123456)
sampling_order <- 1:5000
iterator <- 0
node_indeces <- c()
group_indeces <- c()
group_iterator <- 0
n_nodes <- 0
while (group_iterator < 500){
  iterator <- iterator + 1
  target_group <- y[[sampling_order[iterator]]]
  temp_length <- length(target_group)
  if ((temp_length < 10) || (temp_length > 80)) next
  if (sum(target_group %in% node_indeces) == 0) {
    group_iterator <- group_iterator + 1
    node_indeces <- c(node_indeces,target_group)
    n_nodes <- n_nodes + temp_length
    print(c(iterator, n_nodes))
    group_indeces <- c(group_indeces, rep(group_iterator, temp_length))
  }
}

t.mtx <- readRDS("amazon_graph.rds")

output = t.mtx[(t.mtx[,1] %in% node_indeces + t.mtx[,2] %in% node_indeces)==2,]
rm(t.mtx)

save(output, file = "amazon_edges_500.RData")
nodes <- list(node_indeces = node_indeces, group_indeces = group_indeces)
save(nodes,file = "amazon_nodes_500.RData")

#Estimation
load("amazon_edges_500.RData")
load("amazon_nodes_500.RData")

node_indeces <- nodes$node_indeces
group_indeces <- nodes$group_indeces
n_nodes <- length(node_indeces)

storage.mode(output) <- "character"
storage.mode(group_indeces) <- "character"
storage.mode(node_indeces) <- "character"


network <- as.network.matrix(output, directed = FALSE)

names <- network::get.vertex.attribute(network, "vertex.names")
new_order <- numeric(n_nodes)
for (i in 1:n_nodes){
  new_order[i] <- which(node_indeces == names[i])  
}

network::set.vertex.attribute(network, "true_partition", 
                              value = nodes$group_indeces[new_order])
true_partition <- nodes$group_indeces[new_order]


######################################################################################
######################################################################################
##############
##############        Model estimation
##############
######################################################################################
######################################################################################

# SBM model estimation

formula <- network ~ edges


result_sbm <- hergm(formula,
                    method = "ml",
                    parameterization = "size",
                    estimate_parameters = TRUE,
                    max_number = 500,
                    parallel = mc.cores,
                    verbose = 0,
                    sample_size = 1e+5,
                    sample_size_multiplier_block = NULL,
                    n_em_step_max = 500,
                    max_iter = 100)

params_sbm <- c(result_sbm$results$parameters,rep(0,4),result_sbm$results$between_parameter)
st.error_sbm <- c(result_sbm$results$st.error,rep(0,4),result_sbm$results$st.error.between)

indicator <- as.numeric(result_sbm$results$partition)

# Full model estimation

formula <- network ~ edges + gwesp(fixed = FALSE, cutoff = 12) + gwdegree(fixed = FALSE, cutoff = 20)

#Starting values of parameters were obtained by experiments with marginal models having
#only one geometrically weighted term unconstrained while the other one had a fixed decay
init_estimate <- c(-1.403, 0.291, 1.161, 1.086, 0.760)

result_full <- hergm(formula,
                     method = "ml",
                     parameterization = "size",
                     estimate_parameters = TRUE,
                     indicator = indicator,
                     max_number = 500,
                     parallel = mc.cores,
                     verbose = 0,
                     sample_size = 1e+5,
                     sample_size_multiplier_block = NULL,
                     max_iter = 100,
                     initial_estimate = init_estimate)

params <- result_full$results$parameters
st.error_ergm <- c(result_full$results$st.error, result_full$results$st.error.between)

params_ergm <- c(params,result_full$results$between_parameter)



######################################################################################
######################################################################################
##############
##############        Goodness-of-fit
##############
######################################################################################
######################################################################################

gof_new <- function(network, indicator, eta, formula, n_trials, threshold = 20, 
                    parameterization, mc.cores) {
  
  num_clust <- length(table(indicator))
  
  net_list <- rep(list(NULL), num_clust)
  for (k in 1:num_clust) {
    print(k)
    if (sum(indicator == k) > 1) {
      v_id <- which(indicator == k)
      net_list[[k]] <- get.inducedSubgraph(network, v = v_id)
    }
  }
  
  stats <- list()
  
  stats$dist <- matrix(0,n_trials, threshold)
  stats$esp <- matrix(0,n_trials, threshold)
  stats$degree <- matrix(0,n_trials, threshold)
  stats$dsp <- matrix(0,n_trials, threshold)
  stats$tr_edges <- matrix(0,n_trials, 1)
  stats$edges <- matrix(0,n_trials, 1)
  
  stats$obs.dist <- numeric(threshold)
  stats$obs.esp <- numeric(threshold)
  stats$obs.degree <- numeric(threshold)
  stats$obs.dsp <- numeric(threshold)
  stats$obs.tr_edges <- numeric(1)
  stats$obs.edges <- numeric(1)
  
  
  if (parameterization == "size"){
    model <- ergm.getmodel(formula, network)
    no_scaling <- c()
    if (length(model$etamap$curved) > 0) {
      for (i in 1:length(model$etamap$curved))
      {
        no_scaling <- c(no_scaling, model$etamap$curved[[i]]$from[2])
      }
    }
  } else {
    no_scaling <- NULL
  }
  
  
  rhs <- paste(deparse(formula[[3]]), collapse = "")
  
  
  test <- mclapply(1:num_clust,gof_one_clust,net_list, no_scaling, parameterization, 
                   eta, rhs, threshold, n_trials, mc.cores = mc.cores, mc.preschedule = F)
  
  for (i in 1:num_clust){
    #test[[1]] <- 
    stats$dist <- stats$dist + test[[i]]$sim.dist[,1:threshold]
    stats$esp <- stats$esp + test[[i]]$sim.espart[,1:threshold]
    stats$degree <- stats$degree + test[[i]]$sim.deg[,1:threshold]
    stats$dsp <- stats$dsp + test[[i]]$sim.dspart[,1:threshold]
    stats$tr_edges <- stats$tr_edges + rowSums(test[[i]]$sim.espart[,-1])
    stats$edges <- stats$edges + rowSums(test[[i]]$sim.espart)
    
    stats$obs.dist <- stats$obs.dist + test[[i]]$obs.dist[1:threshold]
    stats$obs.esp <- stats$obs.esp + test[[i]]$obs.esp[1:threshold]
    stats$obs.degree <- stats$obs.degree + test[[i]]$obs.deg[1:threshold]
    stats$obs.dsp <- stats$obs.dsp + test[[i]]$obs.dsp[1:threshold]
    stats$obs.tr_edges <- stats$obs.tr_edges + sum(test[[i]]$obs.esp[-1])
    stats$obs.edges <- stats$obs.edges + sum(test[[i]]$obs.esp)
  }
  
  return(stats)
}

gof_one_clust <-function(k, net_list, no_scaling, parameterization, eta, rhs, threshold, n_trials) {
  cat(paste("\n\n Simulating neighborhood: ", k))
  cur_net <- net_list[[k]]
  rm(net_list)
  cur_eta <- eta
  
  if (parameterization == "size"){
    if (!is.null(no_scaling)){
      cur_eta[setdiff(1:length(cur_eta),no_scaling)] <- eta[setdiff(1:length(cur_eta),no_scaling)]*
        log(cur_net$gal$n)
    }else {
      cur_eta[1:length(eta)] <- cur_eta[1:length(eta)] * log(cur_net$gal$n)
    }
  }
  
  
  formula <- as.formula(paste("cur_net ~ ",rhs,sep=""))
  
  test <- gof(formula, coef = cur_eta, GOF = ~ degree + espartners + distance + dspartners,
              control = control.gof.formula(nsim=n_trials,MCMC.burnin=5e+4))
  
  if (length(test$sim.dist[1,]) < threshold + 1){
    n_dim <- threshold + 1 - length(test$sim.dist[1,])
    test$sim.dist <- cbind(test$sim.dist[,1:(length(test$sim.dist[1,])-1)], matrix(0,n_trials,n_dim))
    test$obs.dist <- c(test$obs.dist[1:(length(test$obs.dist)-1)], numeric(n_dim))
  }
  
  if (length(test$sim.espart[1,]) < threshold + 1){
    n_dim <- threshold + 1 - length(test$sim.espart[1,])
    test$sim.espart <- cbind(test$sim.espart[,1:(length(test$sim.espart[1,])-1)], matrix(0,n_trials,n_dim))
    test$obs.esp <- c(test$obs.esp[1:(length(test$obs.esp)-1)], numeric(n_dim))
  }
  
  if (length(test$sim.deg[1,]) < threshold + 1){
    n_dim <- threshold + 1 - length(test$sim.deg[1,])
    test$sim.deg <- cbind(test$sim.deg[,1:(length(test$sim.deg[1,])-1)], matrix(0,n_trials,n_dim))
    test$obs.deg <- c(test$obs.deg[1:(length(test$obs.deg)-1)], numeric(n_dim))
  }
  
  if (length(test$sim.dspart[1,]) < threshold + 1){
    n_dim <- threshold + 1 - length(test$sim.dspart[1,])
    test$sim.dspart <- cbind(test$sim.dspart[,1:(length(test$sim.dspart[1,])-1)], matrix(0,n_trials,n_dim))
    test$obs.dsp <- c(test$obs.dsp[1:(length(test$obs.dsp)-1)], numeric(n_dim))
  }
  test
}

set.seed(123456)

full_new_stats <- gof_new(network = network, 
                          indicator = indicator,
                          eta = result_full$results$parameters,
                          formula = result_full$formula,
                          n_trials = 1000,
                          threshold = 20,
                          parameterization = "size",
                          mc.cores = mc.cores)



sbm_stats <- gof_new(network = network, 
                     indicator = indicator,
                     eta = result_sbm$results$parameters,
                     formula = result_sbm$formula,
                     n_trials = 1000,
                     threshold = 20,
                     parameterization = "size",
                     mc.cores = mc.cores)








# Plot for paper


par(mfrow=c(1, 1)) 

pdf("gof_ergm.pdf", width=16, height=4)
range_esp <- c(1:13)
range_dist <- c(1:10)
range_dsp <- c(1:13)
par(mfrow=c(1, 4)) 
opar <- par() 
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.8,3.6,1.1,1.1))
x_axis_size <- 1.9
#color <- "lightgray"
color <- "black"

plot(as.numeric(full_new_stats$obs.dist[range_dist]),type='l',ylim = c(0,90000), ylab = "", 
     xlab = "Geodesic distance",  
     #    xaxt="n", 
     cex.lab=x_axis_size)
boxplot(full_new_stats$dist[,range_dist],col=color, add = TRUE,  xaxt="n",outcex=0.3,pch=16,
        medcol = color, whiskcol = color, staplecol = color, 
        outcol = color, boxcol = color,
        pars = list(boxwex = 0.6))
lines(range_dist,as.numeric(full_new_stats$obs.dist[range_dist]),lwd=2,col='red')
#axis(1, at=1:20*1,labels=1:20*1, col.axis="black", las=1)


plot(as.numeric(full_new_stats$obs.dsp[range_dsp]),type='l',ylim = c(0,55000), ylab = "", 
     xlab = "Number of dyadwise shared partners",  
     xaxt="n", 
     cex.lab=x_axis_size)
boxplot(full_new_stats$dsp[,range_dsp],col=color, add = TRUE,  xaxt="n",outcex=0.3,pch=16,
        medcol = color, whiskcol = color, staplecol = color, 
        outcol = "white", boxcol = color,
        pars = list(boxwex = 0.6))
lines(range_dsp,as.numeric(full_new_stats$obs.dsp[range_dsp]),lwd=2,col='red')
axis(1, at=seq(1,13,2),labels=seq(1,13,2)-1, col.axis="black", las=1)



plot(as.numeric(full_new_stats$obs.esp[range_esp]),type='l',ylim = c(0,12000), ylab = "", 
     xlab = "Number of edgewise shared partners",  
     xaxt="n", 
     cex.lab=x_axis_size)
boxplot(full_new_stats$esp[,range_esp],col=color, add = TRUE,  xaxt="n",outcex=0.3,pch=16,
        medcol = color, whiskcol = color, staplecol = color, 
        outcol = "white", boxcol = color,
        pars = list(boxwex = 0.6))
lines(range_esp,as.numeric(full_new_stats$obs.esp[range_esp]),lwd=2,col='red')
axis(1, at=seq(1,13,2),labels=seq(1,13,2)-1, col.axis="black", las=1)

hist(full_new_stats$tr_edges, xlim=c(22000,32000), xlab = "Transitive edges",  ylab= "",
     cex.lab=x_axis_size, main = "")
abline(v = full_new_stats$obs.tr_edges, col = 'red')
box(col = 'black')


dev.off()





pdf("gof_sbm.pdf", width=16, height=4)
range_esp <- c(1:13)
range_dist <- c(1:10)
range_dsp <- c(1:13)
par(mfrow=c(1, 4)) 
opar <- par() 
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.8,3.6,1.1,1.1))

plot(as.numeric(full_new_stats$obs.dist[range_dist]),type='l',ylim = c(0,90000), ylab = "", 
     xlab = "Geodesic distance",  
     #    xaxt="n", 
     cex.lab=x_axis_size)
boxplot(sbm_stats$dist[,range_dist],col='black', add = TRUE,  xaxt="n",outcex=0.3,pch=16,
        medcol = "black", whiskcol = "black",  staplecol = "black",outcol = "black",
        pars = list(boxwex = 0.6))
lines(range_dist,as.numeric(full_new_stats$obs.dist[range_dist]),lwd=2,col='red')
#axis(1, at=1:20*1,labels=1:20*1, col.axis="black", las=1)

plot(as.numeric(full_new_stats$obs.dsp[range_dsp]),type='l',ylim =c(0,55000), ylab = "", 
     xlab = "Number of dyadwise shared partners",  
     xaxt="n", 
     cex.lab=x_axis_size)
boxplot(sbm_stats$dsp[,range_dsp],col='black', add = TRUE,  xaxt="n",outcex=0.3,pch=16,
        medcol = "black", whiskcol = "black",  staplecol = "black",outcol = "white",
        pars = list(boxwex = 0.6))
lines(range_dsp,as.numeric(full_new_stats$obs.dsp[range_dsp]),lwd=2,col='red')
axis(1, at=seq(1,13,2),labels=seq(1,13,2)-1, col.axis="black", las=1)

plot(as.numeric(full_new_stats$obs.esp[range_esp]),type='l',ylim = c(0,12000), ylab = "", 
     xlab = "Number of edgewise shared partners",  
     xaxt="n", 
     cex.lab=x_axis_size)
boxplot(sbm_stats$esp[,range_esp],col='black', add = TRUE,  xaxt="n",outcex=0.3,pch=16,
        medcol = "black", whiskcol = "black",  staplecol = "black",outcol = "white",
        pars = list(boxwex = 0.6))
lines(range_esp,as.numeric(full_new_stats$obs.esp[range_esp]),lwd=2,col='red')
axis(1, at=seq(1,13,2),labels=seq(1,13,2)-1, col.axis="black", las=1)


hist(sbm_stats$tr_edges, xlim=c(22000,32000), xlab = "Transitive edges", ylab= "",
     cex.lab=x_axis_size, main = "")
abline(v = sbm_stats$obs.tr_edges, col = 'red')
box(col = 'black')

dev.off()


#Table from Application section
final_table <- matrix(0,6,4)
final_table[,1] <- params_sbm
final_table[,2] <- st.error_sbm
final_table[,3] <- params_ergm
final_table[,4] <- st.error_ergm

print(final_table)
