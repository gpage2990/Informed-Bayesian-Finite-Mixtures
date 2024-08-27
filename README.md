# Informed-Bayesian-Finite-Mixtures
This folder contains all code required to carry out simulation and galaxy data analysis detailed in paper entitled "Informed Bayesian Finite Mixture Models via Asymmetric Dirichlet Priors".  
More details associated with each object in the fold follow


All three scripts contained in this folder require the "mixcPack" R-package that is available at "https://github.com/gpage2990/miscPack".  Once the package is installed, the "dataAnalysis_Galaxy.R" script is self-contained and should run without any user input.  

Running the simulation study will require user input.  In particularly the R-script "simulationStudy_asym2.R" is set up to be run using "gnu parallel" (https://www.gnu.org/software/parallel/) by executing the "execute_simulation.sh" bash script.  This can be done using the following terminal command

./execute_simulation.sh | parallel &

To facilitate reproducibility, we provide the results from the simulation study in the zipped folder "ss_results14.zip"


Running the "dataAnalysis_Biomechanic.R" script at this time is not possible as I have not yet been granted permission to make the data public by my collaborators.  But, I include the the R-script nonetheless. 