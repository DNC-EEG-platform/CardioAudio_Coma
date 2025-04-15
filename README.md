# CardioAudio_Coma
REFERENCE

Repository containing the Matlab code for the scientific article :

Cardiac signals inform auditory regularity processing in the absence of consciousness.

Andria Pelentritou, Christian Pfeiffer, Manuela Iten, Matthias Haenggi, Frédéric Zubler, Sophie Schwartz, Marzia De Lucia

DEPENDENCIES

Matlab (version 2019b or later; The MathWorks, Natick, MA)

Fieldtrip (version used 20201205; https://www.fieldtriptoolbox.org/download/)

EEGLAB (version 13.4.4b; https://sccn.ucsd.edu/eeglab/download.php)

CONTENTS

Run_ControlAnalyses.m = performs quality control analyses related to experimental design from trigger information and ECG data

Run_ClusterAnalyses.m = performs condition contrast using cluster permutation statistical analyses on preprocessed ERP data

run_clustperm.m = function called in RunClusterAnalyses.m to run the cluster permutation statistical analysis step in fieldtrip

Run_LinearSVMClassifier.m = performs multivariate decoding linear SVM classification of preprocessed ERP data

roc2.m = function called in Run_LinearSVMClassifier.m to compute the ROC curve

rocarea2.m = function called in Run_LinearSVMClassifier.m to compute the AUC

Extract_CardiacResponse.m = detects omission trial related RR intervals from ECG data

INPUT DATA

Input data paths have to be set within each script.

Data are not publicly available due to ethical constraints but can be made available upon reasonable request.

AUTHORS

Author: Andria Pelentritou

PI: Marzia De Lucia

Brain-Body and Consciousness Laboratory (LNCC)

Department of Clinical Neuroscience (DNC)

Lausanne University Hospital (CHUV) and University of Lausanne (UNIL)

Mont-Paisible 16, CH-1011 Lausanne, Switzerland

Email: andria.pelentritou@gmail.com

Last updated: April 2024

