# CereBRLSIM
This repo corresponds to the manuscript: Physics-Informed Digital Twin Can Predict Cerebral Blood Flow and Cerebral Vascular Regulation in Neurocritical Care Patients by Briggs et al. Here, we construct and validate a digital twin of cerebral hemodynamics that can be used to estimate cerebral vascular regulation and cerebral blood flow. 

# CereBRLSIM.m
This script contains the CereBRLSIM ODE model. In our studies, this script was simulated using Runge Kutta (2,2). 

## InVivo_Experiments
This folder contains code for in vivo experiments corresponding to figures 2-3 in the manuscript. There are three dynamic experiments (corresponding to fig. 2) titled "dCA.m, FMD.m, NVC.m" each of these solve to a slightly different CereBRLSIM script depending on the predominant driver in the experiment. The only difference between the CereBRLSIM scripts are the way the drivers are handled. 

Besides handling of the drivers that are unique to the in vivo experiments, there is a slight difference between how normalization factors Smyo are treated in this folder and the primary CereBRLSIM given in the parent directory and used for modeling neurocritical care patient data. A fundamental challenge in the development of this digital twin is that the values that the CVR parameters would take in patient data. Because the raw values were not the main focus in the in vivo experiments, we set them to 0 or 1, representing "off" or "on". All of these experiments were set with Smyo = 1/r0^4. However, after the neurocritical care experiments were conducted, it became clear that with Smyo = 1/r0^4, the value of Cmyo was an order of magnitude smaller than the Cendo and Cmeta. Therefore, we changed Smyo = 10/r0^4 to make the values comparable (as discussed in the manuscript methods section). To remain true to the manuscript description, we left the in vivo tests as is. To replicate the in vivo tests using Smyo = 10/r0^4, Cmyo will simply need to be set to 10 instead of 1. Because the in vivo experiments were conducted to validate the CereBRLSIM model and not estimate or focus on CVR parameter values, this makes fact has negligable impact on the results. 

Real data are avaliable via reasonable request to Phillip Ainslie (co-author on paper). Fake data are provided in folder "FakeData" to provide idea of structure.
