# Environment and Installation

1. Make sure to use the right CMS_SW version: `CMSSW_10_6_30`

2. Clone the git project.


# Create Templates for the single saturated strip clusters

The correction is based on cross-talk inversion. To do so, templates are computed layer by layer to adjust the correction depending on the cluster shape. These templates are established thanks to samples of two generated muons back-to-back (with a $p_T$ from 50 to 200 GeV and $\eta$ from -3 to 3), with PU under Run2 UL conditions. The cluster collection of interest is an intermediate collection stored up to the AOD level which provides information on the simulated deposited charge on the strips associated to the muons, with no electronic saturation applied yet, without noise, without PU contamination, without zero suppression.


Use the script `Script_CreateTemplate.sh` to create the templates. It will run the code `CreateTemplate.C` on samples of two muons back-to-back and store those templates as `.txt` files in the folder `Template_correction`.


It also creates `Check_Template.root` if you want to view the coefficient of correction used according to the layer and the shape.


# Test the correction on the samples itself

A little `.C` file exists to test and view the quality of the correction itself on single and double consecutive saturated strip clusters (on the output file `TestOnSigdigi.root`).


Use the script `Script_TestNewCorr.sh` to do so.
