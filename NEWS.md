# AFTCoop News

# AFTCoop 0.2.4 (Development Version) : 

- Improved functions marginal_ranking.R and variable_screening.R

# AFTCoop 0.2.3 (Development Version) : Branched from 0.2.1 

- To ensure compatibility with R 4.5.2 + Apple Silicon M4.

we removed:  mclapply(...) 

and replaced with: cl <- makeCluster(...);   parLapply(...) ; stopCluster(cl)

- Now, it uses the same parallel interface without need of swithching on windows. 

- The parallel interface is still based on the use of the library "parallel".

# AFTCoop 0.2.2 (Development Version)

- To ensure compatibility with R 4.5.2 + Apple Silicon M4.

We used tryCatch() to handle error message in the parallelism over rho-values and k-fold. 
Then we force sequential mode in case of errors--> but it becomes slower.

# AFTCoop 0.2.1 (Development Version)

- Modified the logical parameter "parallel" into two logical parameters "parallel_rho" and "parallel_cv", to better handle parallelism over the vector of rho values and the cross-validation.

# AFTCoop 0.2.0 (Development Version)

- Modified the call to the evaluation of the maximun eigenvalue of the Hessian in the main_functions using RSpectra.

- Corrected some typos in the documentation 

# AFTCoop 0.1.1 (Development Version)

- This is the first online version. 

### 🏛 Project Funding

This work is supported by the PRIN 2022 PNRR P2022BLN38 project, *Computational approaches for the integration of multi-omics data* funded by European Union - Next Generation EU, CUP **B53D23027810001**.

