Massive-scale estimation of exponential-family random graph models with additional structure

Contents:
I.   R package hergm.
II.  R scripts.
III. Data.

Please note: 

1. To execute the attached R scripts, please make sure that the following R packages are installed: hergm (included), parallel, car.

2. The simulations and data analyses were executed on a 48-core machine. 
The 48-core machine finished the simulation study and data analyses within three weeks. 
To facilitate reproducibility on a single-core machine, we set the number of cores in the attached R scripts to 1: mc.cores <- 1. 
However, please note that reproducing the simulation and application results on a single-core machine may take up to 24 months.

3. The simulations and data analyses require ~ 8GB RAM for working with large matrices and ~ 12 GB disk space to store the simulated networks used in the simulation study.

I. R package hergm:

hergm_3.2-0.tar.gz

II. R scripts:

- simulation_study_final.r: Simulation study with output:

	* small_30_boxplot.pdf: Agreement of estimated and data-generating neighborhood structure in terms of Yule's coefficient for small networks, corresponding to Figure 1.

	* final_table: Computing time in seconds: two-step likelihood-based approach versus Bayesian approach, corresponding to Table 2.

	* large_accuracy.pdf: Agreement of estimated and data-generating neighborhood structure in terms of Yule's coefficient for large networks, corresponding to Figure 2.

	* large_accuracy.pdf: Estimates of within-neighborhood edge and transitive edge parameters, corresponding to Figure 3.
	
- amazon_final.r: Data analyses with output:

	* gof_ergm.pdf: plot of goodness-of-fit of curved exponential-family random graph model, corresponding to Figure 4.

	* gof_sbm.pdf: plot of goodness-of-fit of stochastic block model, corresponding to Figure 5.

	* final_table: estimates of parameters and standard errors of both curved exponential-family random graph model and stochastic block model, corresponding to Table 3.

III. Data: Amazon product network:

- amazon_graph.rds: Amazon product network collected by Yang and Leskovic (2015) and downloaded from

http://snap.stanford.edu/data/com-Amazon.html

- com-amazon.top5000.cmty.txt: top 5000 communities according to Yang and Leskovic (2015).

- Amazon product network corresponding to the top 500 networks with 10 to 80 products, where the ranking of communities is based on Yang and Leskovic (2015).
The resulting network can be found in the following two files:
 
  * amazon_nodes_500.RData: node indices and ground-truth neighborhood labels.

  * amazon_edges_500.RData: edgelist.

