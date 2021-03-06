Last updated: April 29, 2020

This folder contains tests of the G-PCCA algorithm on several models. Read this file to learn how to run the algorithm.

The original algorithm that I got from Bernhard Reuter is run from the file main_interactive.m. It is completely interactive and requires the user to input the count matrix and select several options while it is running. A slightly modified version is in main.m which takes the file name of the count matrix and kmin and kmax as arguments. 

I wanted to automate the process so I split up the algorith into 2 parts:

1. step_1.m: Perform the Schur decomposition of the count matrix and compute minChi for the user-selected cluster numbers. The (most important) outputs of this part are the Schur vectors, eigenvalues, and minChi values.

	The input to step 1 is the name of the count matrix and the minimum and maximum number of clusters for which to compute minChi values. The matrix should be in a tsv file with .txt extension. If the count matrix file is called countmatrix.txt, then:
	
	The output of step 1 is the Schur vectors, in a file countmatrix-etc-X.txt. The eigenvalues are in the file countmatrix-etc-RR.txt. The full transition probability matrix is in the file countmatrix-etc-P.txt. The minChi values are in countmatrix-etc-minChi-n=kmin-kmax.txt and are plotted in countmatrix-etc-minChi-n=kmin-kmax.png. Other outputs have similar names and are detailed at the beginning of the matlab scripts.
	
	example command for running from the command line if I am running from this folder and have put the count matrix in a folder called Results/countmatrix (which is what the shell script that I use to run G-PCCA does):
	
	matlab -r "step_1 "Results/countmatrix/countmatrix" kmin kmax; exit"  

2. step_2.m and step_2_klist.m: Perform G-PCCA for the user-selected cluster number. This uses the Schur vectors and eigenvalues computed in step_1. The output are the optimized reduced transition probability matrices foe the selected cluster numbers, if those numbers do not split up 2x2 blocks of eigenvalues in the Schur decomposition (pairs of complex eigenvalues). step_2 takes as input a minimum cluster number kmin and a maximum cluster number kmax and optimized for [kmin, kmax]. step_2_klist takes as input a list of cluster values which need not be continuous.

	The input to step 2 is the name of the count matrix, and either a kmin and kmax of the cluster numbers for which to optimize or a list of cluster numbers. The script will retrieve the eigenvalues and Schur vectors that it needs using the count matrix name supplied.
	
	The output is the results of the first and second optimization steps for the k values that were supplied and which were not disallowed because they break up 2x2 clusters of complex eigenvalues in the Schur form. The first optimization is performed using Nelder-Mead and the resulting reduced transition probability matrices  are in countmatrix-etc-nelder-mead-etc-n=k-Pc.txt and plotted in the .pdf with the same name. The second optimization is performed using Gauss-Newton and the resulting Pc are in countmatrix-etc-gauss-newton-n=k-Pc.txt and plotted in the .pdf with the same name. The crispness of the model is in files ending in n=k-crispness.txt and, if there are more than one value of k that the optimization was done for, crispness is plotted vs cluster number in the file ending in -crispness-figure-n=kmin-kmax.png. Again, there are many other outputs to discover upon inspection.
	
	example command corresponding to the one above:
	
	matlab -r "step_2 "Results/countmatrix/countmatrix" kmin kmax; exit"
	or
	matlab -r "step_2_klist "Results/countmatrix/countmatrix" [ k1 k2 k3 k4 ]; exit"
	
	note that the list of k values must have spaces between the opening bracket and first entry and the the last entry and closing bracket when running from the command line or else the list will not be parsed properly by the program! 
	good: [ k1 k2 k3 k4 ]
	bad:  [k1 k2 k3] or [ k1 k2 k3] etc
	
In addition, I've written some batch file and bash scripts to facilitate running the code on midway. The workflow for running G-PCCA on midway that I have set up is as follows:

- Put the count matrices to be analyzed in the folder ~/G-PCCA/Count_Matrices. 	COUNT MATRICES MUST HAVE .TXT EXTENSION
- Symlink to this folder (~/G-PCCA) in the scratch folder (~/scratch/G-PCCA)
	- NOTE ON THIS: creating the results files in a symlinked folder has weird behavior. It may be better to just copy the code and count matrices to wherever it needs to be run (scratch) and then copy the results to a long term storage location. 
- Modify run_both_steps_sbatch.sh with the appropriate parameters for the name of the count matrix and k values for each step of the analysis
- Running run_both_steps_sbatch.sh will write sbatch files for each set of parameters and submit jobs for step 1 and step 2. It makes sure that step 1 has successfully completed before running step 2. Right now step 2 will still run if step 1 fails, but it will also fail promptly. A to-do item is to modify this so that the step 2 job gets cancelled if step 1 fails.

