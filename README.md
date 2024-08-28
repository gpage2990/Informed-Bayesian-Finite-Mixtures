# Informed-Bayesian-Finite-Mixtures
This folder contains all code created to carry out simulations, galaxy data analysis, and functional data analysis that are detailed in paper entitled "Informed Bayesian Finite Mixture Models via Asymmetric Dirichlet Priors".  More details associated with each object in the repository follow


All R scripts contained in this folder require the "micsPack" R-package that is available at "https://github.com/gpage2990/miscPack".  This package can be installed after loading the "devtools" library and using the following in an R session

install_github("gpage2990/miscPack")

Once the package is installed, the "dataAnalysis_Galaxy.R" script is self-contained and should run without any user input.  

Running the univariate simulation study will require user input.  In particularly the R-script "simulationStudy_asym2.R" is set up to be run using "gnu parallel" (https://www.gnu.org/software/parallel/) by executing the "execute_simulation.sh" bash script.  This can be done using the following terminal command

./execute_simulation.sh | parallel &

To facilitate reproducibility, we provide the results from the simulation study in the folder "ss_results".  The figures contained in the paper can be reproduced using code that appears after the simulation study in R-script "simulationStudy_asym2.R".


Running the functional data simulation by way of the R-script "simulationStudy_asym_multi.R" and the biomechanics analysis via "dataAnalysis_Biomechanic_rev.R" is not possible at this time as the data have not been made publicly available.  That said, I include them for completeness.  The simulation results for the functional type data are included in the folder "ss_multi_results" and the code after that which runs executes the simulation in the file ""simulationStudy_asym_multi.R" can be used to read in results from the folder and produce figures in the manuscript.



