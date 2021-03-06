# CI-for-fraction-who-benefit
This code corresponds to the paper "Constructing a Confidence Interval for the Fraction who Benefit From Treatment, Using Randomized Trial Data" by Emily Huang, Ethan Fang, Daniel Hanley, and Michael Rosenblum.

The folder entitled functions

This folder gives the code for our method in the case of no restrictions and assignment probability = 1/2. We demonstrate how to apply the code with an example data set. 

The implementation of our method is presented in Appendix J of the Supplementary Materials of our paper. We provide both R and MATLAB code implementing our method. The MATLAB code is in the script entitled “CI_noRes.m”. This code outputs the confidence interval constructed using our method \hat{CI}_n, in the case of no restrictions and assignment probability theta = 1/2. We demonstrate how to apply the code with an example data set taken from our simulations. We also provide R code in the file entitled “CI_noRes.R”. This R code performs the same tasks as the MATLAB code. The results from the R and MATLAB code can differ slightly since they use different quadratic program solvers.

The folder entitled CLEARIII

There are three subfolders:
1) “constructDataset”: We format the CLEAR III data set before applying our method and the m-out-of-n bootstrap. The data set is not included because it is proprietary. 
2) “applyMethod”: We apply our method \hat{CI}_n to the CLEAR III data set. The MATLAB code is provided separately for each outcome.
3) “applyMoutOfN”: We apply the m-out-of-n bootstrap, using different choices for m including m = n, m = 0.9n, m = 0.75n, m = 0.5n, m = 0.25n. A function for applying the m-out-of-n bootstrap for these choices of m is in the R script “mOfn_code.R”. The code for applying the m-out-of-n bootstrap to the CLEAR III data set is in the R script “computingCIs.R”.
4) “applyHorowitzManski”: We apply the method by Horowitz and Manski (2000). R code is provided.


The folder entitled limitDist_underNull

In the subfolder “binary_noRes_5050”, we take 100,000 draws from the asymptotic distribution of the test statistic under Simulation Setting A when psi = 0.5. We use the code “Fig2.R” to generate Figure 2 of the main paper (“Fig2.pdf”), which is a histogram of these 100,000 draws. 

The folder entitled simulations

There is a subfolder for each setting of the simulation studies. Setting A corresponds to “binary_noRes_5050”; Setting B corresponds to “binary_noHarm_5050”; Setting C corresponds to “binary_noRes_7550”; Setting D corresponds to “MISTIE_RICV5”. Within each subfolder, there are three folders:
1) “testDatasets”, in which the 5000 (or 1000 in the case of Setting D) trials are simulated
2) “ourMethod”, which applies our method \hat{CI}_n to each of the simulated trials
3) “mOfn”, which applies the competitor method using m-out-of-n bootstrap to each of the simulated trials
4) “manski_horowitz”, which applies the competitor method proposed by Horowitz and Manski (2000) to each of the simulated trials
5) “chernozhukov” (only in the “binary_noRes_5050” subfolder), which applies the competitor method proposed by Chernozhukov et al (2013) to the simulated trials

The code “plotResults.R” plots the results. It generates Figures 3-4 from the main paper, and Figures 2-16 and Tables 1-4 from the Supplementary Materials.

All data sets and results are included in the folders for Settings A, B, and C. The results are presented in matrices; each matrix row corresponds to a specific simulated trial. The data sets for Setting D are not included since the MISTIE II data is proprietary.

In the SettingD folder, there is an R script entitled “SuppMat_Fig1_code.R” that generates a plot of the distribution of the outcome RICV5 observed in the MISTIE II trial. This plot is Figure 1 of the Supplementary Materials.

The folder entitled power

This folder corresponds to the simulation study on power in the Supplementary Materials of our paper. There are two subfolders called “n10000” and “plots”. In the “n10000” subfolder, we provide the code used to conduct the simulation study. In the code, we simulate many trial data sets and perform our hypothesis test using each data set at certain values of psi (see Supplementary Materials for details). In the “plots” subfolder, we visualize the results of the simulation study. 
