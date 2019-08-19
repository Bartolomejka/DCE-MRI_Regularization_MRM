# Spatially regularized estimation in DCE-MRI
This repository contains implementation of the proximal Newton algorithm to estimate perfusion parameters in DCE-MRI with spatial regularization in the form of total variation. The algorithm is described in the paper:
*M. Bartoš, P. Rajmic, M. Šorel, M. Mangová, O. Keunen and R. Jiřík. Spatially regularized estimation of the tissue homogeneity model parameters in DCE-MRI using proximal minimization. Magnetic Resonance in Medicine, 2019; 00: 1– 16. https://doi.org/10.1002/mrm.27874.*
Preprint available at: *http://www.utko.feec.vutbr.cz/~rajmic/papers/Bartos_etal_RegularizedDCEMRI_web.pdf.*

## Files and directories
`/Data` - folder with input DCE-MRI data  
`/Viewer` - viewer of the results  
`start_DCEMRI_estimation.m` - **starting script** defining input and additional parameters, running the main algorithm and visualizing/saving results  
`proximalLM_ChP.m` - main algorithm - proximal Levenberg-Marquardt combined with Primal-Dual Chambolle-Pock algorithm to solve the denoising sub-problem  
`PD_denoising_H_TV.m` - Primal-Dual Chambolle-Pock algorithm to solve the denoising sub-problem  
`y_minus_Cfit_L2_2.m` - Least-Mean-Squares function evaluating model and computing Hessian using Gauss-Newton approximation  
`C_fit_FT.m` - Pharmacokinetic model computing the convolution of AIF and IRF in the Fourier domain  
`TH_Sourbron2_FT.m` - Tissue Homogeneity model - IRF  
`aif_nonblind_real.m` - Arterial Input Function (AIF) model. It is no analytical model - just handling of sampled AIF.  
`fgrad_1.m` - computes image gradients on a mask. Non-mask version comes originally from Kristian Bredies (TGV JPEG reconstruction).  
`bdiv_1.m` - computes divergence (adjoint to gradients) on a mask. Non-mask version comes originally from Kristian Bredies (TGV JPEG reconstruction).  
`function_stdEst2D.m` - noise standard deviation estimator based on the Median of Absolute Deviation (MAD). It comes from Alessandro Foi.  


## Run estimation
The estimation of the perfusion parameters is done by the script `start_DCEMRI_estimation.m`. The script inicializes the necessary constants and variables.

### Initialization
At the beginnig of the starting script, a dataset is chosen by un/commenting the line defining the path to the file. Additionally, an initial starting point `x0_irf` is defined here. E.g. the last element in the vector is the bolus arrival time, which should be updated in case of new dataset.

### Regularization
The regularization can be influenced by scalar parameter `gamma_general` at the beginnig of the script `proximalLM_ChP.m`. Larger value means stronger total variation regularization. The regularization can be switched off by setting `TV_regularization=false`, leading to standard Levenberg-Marquardt algorithm.

### View results
The saved results of the computation can be opened by a provided interactive viewer located in  folder. It can be run by `VISUALIZATION.M` from `/Viewer/gui_mri/visualization`. Just load the file and then press *Show* button.

Alternativelly, the results can be in a standard way loaded into Matlab and then the perfusion parameter maps can be created by e.g. `map = form_image(cell2mat(tissue(:,5)),tissue(:,2));` to create the map of the perfusion parameter in the fifth column. The columns are described in `info.tissue.descr`.

### Other options
Other parameters influencing the performance of the algorithm are described in the code.

---
*Terms of Use*

This code can be freely used for research purposes by academic organizations.
If you use it in your research, please cite the paper above.
