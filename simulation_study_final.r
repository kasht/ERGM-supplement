library(hergm)
library(parallel)

set.seed(123456)

#number of cores, set equal to 1 for Windows
#Run with mc.cores = 45 for the paper
mc.cores = 1;

#Set path for the root directory

path <- "/home/svb1/local_TERGM"

#The code assumes the following folder structure

# path/networks/small/same
# path/networks/small/diff
# path/networks/large/same
# path/networks/large/diff
# path/results/small/bayesian
# path/results/small/30_nodes
# path/results/large
# path/results_dif_sizes/small/bayesian
# path/results_dif_sizes/small/30_nodes
# path/results_dif_sizes/large2

generate_network <- function(index,formula, coef_wi, coef_bw, true_memb){
  set.seed(index)
  network <- simulate_hergm(formula, coef_wi, coef_bw, true_memb, parameterization =  "size")
  name = paste("network",index,sep=".")
  save(network, file=name)
}


run_simulation_fast <- function(index,path, 
                                n_em_step_max, parameterization, n_clusters, 
                                max_iter, parallel,
                                sample_size = 1000, burnin = 2e+4, interval = 500){
  load(paste(path,paste("network",index,sep="."),sep = "/"))
  set.seed(index)
  
  a <- Sys.time()
  
  result <- hergm(network ~ edges + transitiveties, 
                  method = "ml", 
                  parameterization = parameterization, 
                  estimate_parameters = TRUE, 
                  max_number = n_clusters, 
                  parallel = parallel, 
                  verbose = 2,
                  n_em_step_max = n_em_step_max,
                  max_iter = max_iter,
                  sample_size = sample_size,
                  burnin = burnin,
                  interval = interval
  )
  b <- Sys.time()
  
  time_diff <- difftime(b,a,units = "secs")
  ###################################################
  
  #g <- asIgraph(as.network(network))
  
  #b <- cluster_walktrap(g, steps = 4)
  
  #c <- cluster_spinglass(g,spins = n_clusters)
  
  
  d<- spec_clust(network, n_clusters)
  
  n <- network$gal$n
  name <- paste("output",n, "fast", parameterization, index, sep = ".")
  output <- list(network = network, result = result, 
                 #memb_walktrap =  b$membership,
                 time_diff = time_diff,
                 #memb_spinglass = c$membership, 
                 memb_spec_clust = d
  )
  
  save(output, file=name)
}

run_simulation_bayesian <- function(index, path, parameterization, max_number){
  load(paste(path,paste("network",index,sep="."),sep = "/"))
  
  set.seed(index)
  
  a <- Sys.time()
  
  result <- hergm(network ~ edges_ij + transitiveties_ijk, 
                  verbose = -1,max_number = max_number, initialize = TRUE,
                  relabel = 1, number_runs = 3)
  
  b <- Sys.time()
  
  time_diff <- difftime(b,a,units = "secs")
  
  n <- network$gal$n
  name <- paste("bayesian.output",n, parameterization, index, sep = ".")
  
  output <- list(result = result, time_diff = time_diff)
  save(output, file=name)
}

#GENERATING NETWORKS
n_iterations <- 500


setwd(paste(path, "/networks/small/same", sep = ""))

formula = as.formula("network ~ edges + transitiveties")
coef_wi = c(-1, 0.5)/log(10)
coef_bw = -3/log(30)
true_memb <- c(rep(1,10),rep(2,10),rep(3,10))

mclapply(X = 1:n_iterations,FUN = generate_network, formula = formula, coef_wi = coef_wi,
         coef_bw = coef_bw, true_memb = true_memb, mc.cores = 50)

setwd(paste(path, "/networks/small/diff", sep = ""))

formula = as.formula("network ~ edges + transitiveties")
coef_wi = c(-1, 0.5)/log(10)
coef_bw = -3/log(30)
true_memb <- c(rep(1,5),rep(2,10),rep(3,15))

mclapply(X = 1:n_iterations,FUN = generate_network, formula = formula, coef_wi = coef_wi,
         coef_bw = coef_bw, true_memb = true_memb, mc.cores = 50)

setwd(paste(path, "/networks/large/same", sep = ""))

formula = as.formula("network ~ edges + transitiveties")
coef_wi = c(-2, 1)/log(25)
coef_bw = -4/log(2500)
true_memb <- rep(1:100,each = 25)

mclapply(X = 1:n_iterations,FUN = generate_network, formula = formula,
         coef_wi = coef_wi,
         coef_bw = coef_bw, true_memb = true_memb, mc.cores = 50)

setwd(paste(path, "/networks/large/diff", sep = ""))

formula = as.formula("network ~ edges + transitiveties")
coef_wi = c(-2, 1)/log(25)
coef_bw = -4/log(2500)
true_memb <- c(rep(1:20,rep(15,20)), rep(21:40,rep(20,20)), rep(41:60,rep(25,20)),
               rep(61:80,rep(30,20)), rep(81:100,rep(35,20)))

mclapply(X = 1:n_iterations,FUN = generate_network, formula = formula,
         coef_wi = coef_wi,
         coef_bw = coef_bw, true_memb = true_memb, mc.cores = 50)



#ESTIMATING

setwd(paste(path, "/results/small/bayesian", sep = ""))

mclapply(X = 1:n_iterations,FUN = run_simulation_bayesian,
         path = paste(path, "/networks/small/same", sep = ""),
         parameterization = "standard", max_number = 3, mc.cores = 36)

setwd(paste(path, "/results_dif_sizes/small/bayesian", sep = ""))

mclapply(X = 1:n_iterations,FUN = run_simulation_bayesian,
         path = paste(path, "/networks/small/diff", sep = ""),
         parameterization = "size", max_number = 3, mc.cores = 36)


setwd(paste(path, "/results/small/30_nodes", sep = ""))

mclapply(X = 1:n_iterations,FUN = run_simulation_fast,
         path = paste(path, "/networks/small/same", sep = ""), n_em_step_max = 200,
         parameterization = "standard", n_clusters = 3, max_iter = 3, parallel = 1,
         sample_size = 1000, burnin = 2e+4, interval = 500,
         mc.cores = 50)

setwd(paste(path, "/results_dif_sizes/small/30_nodes", sep = ""))

mclapply(X = 1:n_iterations,FUN = run_simulation_fast,
         path = paste(path, "/networks/small/diff", sep = ""), n_em_step_max = 200,
         parameterization = "size", n_clusters = 3, max_iter = 3, parallel = 1,
         sample_size = 1000, burnin = 2e+4, interval = 500,
         mc.cores = 50)

setwd(paste(path, "/results/large", sep = ""))

mclapply(X = 1:n_iterations,FUN = run_simulation_fast,
         path = paste(path, "/networks/large/same", sep = ""), n_em_step_max = 800,
         parameterization = "standard", n_clusters = 100, max_iter = 2,
         sample_size = 500, burnin = 5e+4, interval = 1000, parallel = 1,
         mc.cores = 50)

setwd(paste(path, "/results_dif_sizes/large2", sep = ""))

mclapply(X = 1:n_iterations,FUN = run_simulation_fast,
         path = paste(path, "/networks/large/diff", sep = ""), n_em_step_max = 2400,
         parameterization = "size", n_clusters = 100, max_iter = 2,
         sample_size = 250, burnin = 5e+4, interval = 1000, parallel = 1,
         mc.cores = 50)


#PROCESSING
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#Find the names of all files
process <- function(path, true_memb, pattern = "^o.+?\\d$") {
  setwd(path)
  files <- list.files(pattern = pattern)
  n_runs <- length(files)
  print(n_runs)
  coefs <- matrix(0,n_runs,2)
  accuracy_em <- matrix(0,n_runs,1)
  accuracy_spec_clust <- matrix(0,n_runs,1)
  accuracy_walk_trap <- matrix(0,n_runs,1)
  partition <- list()
  times <- numeric(n_runs)
  #accuracy_sim_ann <- matrix(0,n_runs,1)
  counter <- 0
  for (i in files){
    counter <- counter + 1
    print(counter)
    load(i)
    coefs[counter,] <- as.numeric(output$result$results$parameters)
    accuracy_em[counter,] <- as.numeric(find_phi_coefficient(true_memb, 
                                                             output$result$results$partition))
    #accuracy_sim_ann[counter,] <- as.numeric(find_phi_coefficient(output$true_memb, output$memb_spinglass))
    accuracy_walk_trap[counter,] <- as.numeric(find_phi_coefficient(true_memb, output$memb_walktrap))
    accuracy_spec_clust[counter,] <- as.numeric(find_phi_coefficient(true_memb, output$memb_spec_clust))
    partition[[counter]] <- output$result$results$partition
    times[counter] <- output$time_diff
  }
  obj <- list(coefs,accuracy_em,accuracy_spec_clust,accuracy_walk_trap, partition, times)
  save(obj,file = "output")
  print("Done")
}

process_bayesian <- function(path, true_memb, pattern = "^o.+?\\d$") {
  setwd(path)
  files <- list.files(pattern = pattern)
  n_runs <- length(files)
  print(n_runs)
  coefs <- matrix(0,n_runs,6)
  accuracy_bayes <- matrix(0,n_runs,1)
  partition <- list()
  #accuracy_sim_ann <- matrix(0,n_runs,1)
  counter <- 0
  n_true_clusters <- length(table(true_memb))
  times <- numeric(n_runs)
  for (i in files){
    counter <- counter + 1
    print(counter)
    load(i)
    partition[[counter]] <- apply(output$result$relabeled.indicator,2,Mode)
    cluster_sizes <- table(factor(partition[[counter]],levels = 1:n_true_clusters))
    coefs_temp <- apply(output$result$relabeled.hergm_theta,2,median)
    coefs_temp[which(cluster_sizes==0 | cluster_sizes==1)] <- NA
    coefs_temp[which(cluster_sizes==0| cluster_sizes==1)] <- NA
    coefs[counter,] <- coefs_temp[c(1:3,5:7)]
    
    accuracy_bayes[counter,] <- as.numeric(find_phi_coefficient(true_memb, 
                                                                partition[[counter]]))
    times[counter] <- output$time_diff
    
  }
  obj <- list(coefs,accuracy_bayes, partition, times)
  save(obj,file = "output")
  print("Done")
}

process_bayesian(paste(path, "/results/small/bayesian", sep = ""), pattern = "\\.[0-9]+",
                 true_memb = c(rep(1,10),rep(2,10),rep(3,10)))

process(paste(path, "/results/small/30_nodes", sep = ""),
        true_memb = c(rep(1,10),rep(2,10),rep(3,10)))

process_bayesian(paste(path, "/results_dif_sizes/small/bayesian", sep = ""), pattern = "\\.[0-9]+",
                 true_memb = c(rep(1,5),rep(2,10),rep(3,15)))

process(paste(path, "/results_dif_sizes/small/30_nodes", sep = ""),
        true_memb = c(rep(1,5),rep(2,10),rep(3,15)))



process(paste(path, "/results/large", sep = ""),
        true_memb = rep(1:100,each = 25))


process(paste(path, "/results_dif_sizes/large2", sep = ""),
        true_memb = c(rep(1:20,rep(15,20)), rep(21:40,rep(20,20)), rep(41:60,rep(25,20)),
                      rep(61:80,rep(30,20)), rep(81:100,rep(35,20))))

#Plotting
#setwd("C:/Users/Sergii/Dropbox/MTSSC/new/code/variational/sim_study_results/")
setwd(path)
library(car)

load(paste(path, "/results/small/30_nodes/output", sep = ""))
small_same <- obj
load(paste(path, "/results_dif_sizes/small/30_nodes/output", sep = ""))
small_diff <- obj
load(paste(path, "/results/small/bayesian/output", sep = ""))
bayes_same <- obj
load(paste(path, "/results_dif_sizes/small/bayesian/output", sep = ""))
bayes_diff <- obj

#Figure 1 
pdf("small_30_boxplot.pdf",width=12,height=7) 
#Small accuracy plot


layout(1)

par(mgp=c(3,1.6,0))
par(mar=c(5.1,4.1,4.1,2.1))
bar1 <- boxplot(small_same[[2]],plot=FALSE)
bar2 <- boxplot(bayes_same[[2]],plot=FALSE)
bar3 <- boxplot(small_same[[3]],plot=FALSE)
bar4 <- boxplot(small_diff[[2]],plot=FALSE)
bar5 <- boxplot(bayes_diff[[2]],plot=FALSE)
bar6 <- boxplot(small_diff[[3]],plot=FALSE)
boxplot(list(small_same[[2]], bayes_same[[2]], small_same[[3]], small_diff[[2]], 
             bayes_diff[[2]],  small_diff[[3]]),
        names = c("Two-step approach:\nbalanced case", "Bayesian approach:\nbalanced case",
                  "Spectral clustering:\nbalanced case", 
                  "Two-step approach:\nunbalanced case", "Bayesian approach:\nunbalanced case",
                  "Spectral clustering:\nunbalanced case"),
        cex.main = 1,
        outline=FALSE,ylim=c(0,1))
set.seed(1)
points(jitter(rep(1, length(bar1$out)),factor = 4, amount = 0.3), bar1$out)
points(jitter(rep(2, length(bar2$out)),factor = 4, amount = 0.3), bar2$out)
points(jitter(rep(3, length(bar3$out)),factor = 4), bar3$out)
points(jitter(rep(4, length(bar4$out)),factor = 4), bar4$out)
points(jitter(rep(5, length(bar5$out)),factor = 4), bar5$out)
points(jitter(rep(6, length(bar6$out)),factor = 4, amount = 0.3), bar6$out)
dev.off()

#Computing times

comp_times <- matrix(0,2,2)

comp_times[1,1] <- mean(small_same[[6]])
comp_times[2,1] <- mean(small_diff[[6]])
comp_times[1,2] <- mean(bayes_same[[4]])
comp_times[2,2] <- mean(bayes_diff[[4]])

print(comp_times)

##############################################################################################
##############################################################################################
#############################
############################# Large study
#############################
##############################################################################################
##############################################################################################

load(paste(path, "/results/large/output", sep = ""))
large_same <- obj
load(paste(path, "/results_dif_sizes/large2/output", sep = ""))
large_diff <- obj


#Large accuracy plot
#Figure 2
pdf("large_accuracy.pdf",width=10,height=7)


layout(1)

par(mgp=c(3,1.6,0))
par(mar=c(5.1,4.1,4.1,2.1))
boxplot(list(large_same[[2]],large_same[[3]], large_diff[[2]], large_diff[[3]]), 
        names = c("Two-step approach:\nbalanced case", "Spectral clustering:\nbalanced case", 
                  "Two-step approach:\nunbalanced case", "Spectral clustering:\nunbalanced case"),
        cex.main = 1, ylim = c(0,1))
dev.off()

#Figure 3

pdf("parameters.pdf",width=12,height=8)
layout(matrix(c(1,2,3,4),2,2))
par(mar=c(3.5,4,2,2))
range <- 0.65
dataEllipse(small_same[[1]][,1]/log(10), small_same[[1]][,2]/log(10), cex = 0.75, robust = T, 
            cex.main = 1.3,pch=19,
            levels=c(0.95),add=F,
            ylim = c(0.5/log(10)-range,0.5/log(10)+range),
            xlim = c(-1/log(10)-range,-1/log(10)+range),
            col= c("black","grey"), center.cex = 0, grid = F,
            xlab=expression(theta[1]),ylab=expression(theta[2]), cex.lab = 2,
            mgp = c(2.3, 0.5, 0))
title("n = 30, K = 3, balanced case", cex.main = 1.6)
points(t(as.matrix(c(-1,0.5)/log(10))),col="red",pch=19,cex=2)


dataEllipse(small_diff[[1]][,1], small_diff[[1]][,2], robust = T, 
            cex = 0.75,cex.main = 1.3,pch=19,
            levels=c(0.95),add=F,
            ylim = c(0.5/log(10)-range,0.5/log(10)+range),
            xlim = c(-1/log(10)-range,-1/log(10)+range),
            col= c("black","grey"), center.cex = 0, grid = F,
            xlab=expression(theta[1]),ylab=expression(theta[2]), cex.lab = 2,
            mgp = c(2.3, 0.5, 0))
title("n = 30, K = 3, unbalanced case", cex.main = 1.6)
points(t(as.matrix(c(-1,0.5)/log(10))),col="red",pch=19,cex=2)

dataEllipse(large_same[[1]][,1]/log(25), large_same[[1]][,2]/log(25), cex = 0.75, robust = T, cex.main = 1.5,
            levels=c(0.95),add=F,
            ylim = c(0.5/log(10)-range,0.5/log(10)+range),
            xlim = c(-1/log(10)-range,-1/log(10)+range),
            pch=19,
            col= c("black","grey"), center.cex = 0, grid = F,
            xlab=expression(theta[1]),ylab=expression(theta[2]), cex.lab = 2,mgp = c(2.3, 0.5, 0))
title("n = 2,500, K = 100, balanced case", cex.main = 1.6)
points(t(as.matrix(c(-2,1)/log(25))),col="red",pch=19,cex=2)

dataEllipse(large_diff[[1]][,1], large_diff[[1]][,2], robust = T, cex = 0.75,cex.main = 1.5,
            levels=c(0.95),add=F,
            ylim = c(0.5/log(10)-range,0.5/log(10)+range),
            xlim = c(-1/log(10)-range,-1/log(10)+range),
            pch=19,
            col= c("black","grey"), center.cex = 0, grid = F,
            xlab=expression(theta[1]),ylab=expression(theta[2]), cex.lab = 2,mgp = c(2.3, 0.5, 0))
title("n = 2,500, K = 100, unbalanced case",cex.main = 1.6)
points(t(as.matrix(c(-2,1)/log(25))),col="red",pch=19,cex=2)
dev.off()

