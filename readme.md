This document describes the download and usage of the code in this repository, both for the calculation of homology as well as for the S test.

Preliminaries:
The cloning of a github depository is a problem with a well-known solution (See e.g. https://help.github.com/articles/fork-a-repo ).

After you have cloned (or downloaded) the depository to a folder of your liking, you will need to run it. You can run it in R on its own or with RStudio (I recommend RStudio).

R:

For the homology calculation, source the three source-files which are in the folder R_Code, namely boundary_calculations.R, found_functions.R and homology_calculation.R, starting with homology_calculation.R. Before the first time, I recommend installing potentially missing packages with the following command:
install.packages(c("Matrix", "MASS", "numbers", "compiler"))
After sourcing the three files, you can then calculate the three homologies by running the homology and degenerate_homology functions from the command line as follows:
n-th Rack homology group of R_k: homology(k,n,FALSE)
n-th Quandle homology group of R_k: homology(k,n,TRUE)
n-th Degenerate homology group of R_k: degenerate_homology(k,n)

Should you wish to calculate the homology group of a different quandle/biquandle, please change the up_action and down_action functions in the boundary_calculations.R file and re-source afterwards. As an example, the up_action for the group quandle over S_3 has been included as an action matrix which is commented out. In this case, one would have to change call the homology and degenerate_homology functions similarly as above, the only difference being that k is then the order of the underlying group.

For the S-test, source the file S.R (which can be found in the folder S_test, which is itself in the folder R_code). Change the up_action and down_action to the actions you would like to check, then source the file again. After this, run 
S_test(k), where k is the order of the underlying group. Currently, the file is only set up for checking the dihedral quandle, but it is easy to adapt.

RStudio:

There are two ways of using this with RStudio. The first is exactly as with R above. The second one involves opening the appropiate project, namely R code.Rproj in the folder R_code for the homology calculations and S test.Rproj in the folder S test. After loading, source the files as explained above, possibly installing the missing packages before usage, as required. Then run the calculations as detailed above.


Question/Problems?
Send an email to a.wenzel(at)sussex.ac.uk
