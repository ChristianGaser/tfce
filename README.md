# TFCE: Threshold-Free Cluster Enhancement
This toolbox is a an extension to [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (Wellcome Department of Cognitive Neurology) to provide non-parametric statistics based on threshold-free cluster enhancement (TFCE). The idea of the TFCE approach is to combine focal effects with large voxel height as well as broad effects. In contrast to common approaches that consider cluster-based thresholding no initial (and arbitrary) cluster-forming threshold is necessary. Furthermore, TFCE is fairly robust to nonstationarity in the data, which can be often seen in VBM data.

The toolbox can be applied to (almost) any existing statistical (parametric) SPM design for 3D volume or surface data.

It is developed by Christian Gaser (Jena University Hospital, Departments of Psychiatry and Neurology) and free but copyright software, distributed under the terms of the [GNU General Public Licence](http://www.gnu.org/licenses/gpl-2.0.html) as published by the Free Software Foundation; either version 2 of the Licence, or (at your option) any later version.

## Download
[TFCE toolbox](http://141.35.69.218/tfce/tfce_latest.zip)

Older version can be obtained [here](http://141.35.69.218/tfce/).

## Quick Start Guide
### Estimate TFCE

Use any existing 2nd-level Model from a previous (parametric) SPM12 analysis of volume or surface data.

### Select SPM.mat
Select the SPM.mat file that contains the design specification from a previous (parametric) estimation, where all required contrasts are already specified.

### Select an additional mask
Select an additional mask image or surface to restrict your analysis. As default the mask in the analysis folder is used. Here you can select a mask to additionally restrict the analysis to regions of interest (i.e. small volume/surface correction).

### Results title
Heading on results page - determined automatically if left empty

### Contrast index
Index(es) of contrast according to the contrast manager: Each contrast in SPM is indicated by a sequential number that is displayed in the first column of the contrast manager. You can enter one or more contrasts. If only one number is entered, and this number is "Inf", you can select one or more contrasts interactively using the contrast manager. Do not define here the contrast itself. This should be done in the contrast manager, that is automatically called if "Inf" is kept as entry.

### Number of permutations
In order to obtain reliable estimates you need about 5000-10000 permutations. There is also an option to interrrupt the permutation process and to save the results at this step to take a first look on your results. If the number of maximal possible permutations is smaller, then this number is used resulting in an exact permutation test.

### Permutation method to deal with nuisance variables
A number of methods are available to obtain parameter estimates and construct a reference distribution in the presence of nuisance variables. Smith permutation method is used if any nuisance variables exist and is selected by default. If no nuisance variables were found in the model then Draper-Stoneman method is automatically used. Freedman-Lane is another permutation method to deal with nuisance parameters. However, behaviour of that method was found to be strange under some circumstances and you have to apply this method very carefully. It's only necessary to change the permutation method if a large discrepancy between parametric and non-parametric statistic was found, which is indicated at the Matlab command line. See also [Winkler et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.01.060) for more information.

### TBSS data
Use 2D optimization (e.g. for TBSS DTI data) with internal TFCE parameters H=2, E=1.

### Weighting of cluster size
The idea of the TFCE approach is to combine focal effects with large voxel height as well as broad effects. The weighting of these effects is defined using the parameters E (extent) and H (height). [Smith and Nichols](https://doi.org/10.1016/j.neuroimage.2008.03.061) (Neuroimage 2009) empirically estimated E=0.5 and H=2 for volume data to provide good statistical power. However, the empirically derived values found to be very sensitive for local effects, but not for broader effects. Thus, you can try to change the weighting if you expect more broader effects or even very broad effects by changing the weighting parameter E. Please note that for surfaces and TBSS data the weightings are different from that of 3D volume data with E=1 and H=2 and will be fixed and not changed by this setting.
