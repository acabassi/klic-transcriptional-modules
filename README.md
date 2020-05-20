# Transcriptional module discovery

The code contained in this folder can be used to reproduce the transcriptional module discovery presented in Cabassi and Kirk (2020).

## Folder structure

### R files

0) **Marina data**. This is the cluster analysis of the expression dataset of Granovskaia et al. (2010) that contains measurements related to 551 genes whose expression profiles have been measured at 41 different time points of the cell cycle.
1) **Harbison data**. This is the cluster analysis of the ChIP-chip dataset of Harbison et al. (2004) which provides binding information for 117 transcriptional regulators for the same genes considered above. The data were discretised as in Kirk et al. (2012). The clustering is done via (a) Bayesian hierarchical clustering, (b) partitioning around medoids with Gower's distance, (c) a greedy Bayesian non-parametric clustering algorithm.
2) **COCA**. Integrative clustering is performed using Cluster-Of-Clusters Analysis (Cabassi and Kirk, 2020) and the clusters found in the previous steps.
3) **KLIC**. Integrative clustering is performed using Kernel Learning Integrative Clustering (Cabassi and Kirk, 2020).

### Subfolders

- **GBNP**: MATLAB code needed to run the greedy Bayesian non-parametric clustering algorithm.
- **Cophenetic correlation coefficients**: cophenetic correlation coefficients of all kernel matrices produced by KLIC.
- **Figures**: all figures produced for the data analysis.
- **Functions**: some auxiliary functions.
- **GOTO scores**: MATLAB code needed to compute GOTO scores.

## References

- Cabassi, A. and Kirk, P.D., 2020. Multiple kernel learning for integrative consensus clustering of 'omic datasets. arXiv preprint arXiv:1904.07701.

- Granovskaia, M. V. et al. (2010). High-resolution transcription atlas of the mitotic cell cycle in budding yeast. Genome biology, 11(3):R24.

- Harbison, C. T. et al. (2004). Transcriptional regulatory code of a eucaryotic genome. Nature, 431(7004):99–104. 

- Kirk, P. D. W. et al. (2012). Bayesian correlated clustering to integrate multiple datasets. Bioinformatics, 28(24):3290–3297.

## Contact

Please feel free to contact me (Alessandra) should you need any help using this code. You can find my up-to-date email address in my GitHub profile.
