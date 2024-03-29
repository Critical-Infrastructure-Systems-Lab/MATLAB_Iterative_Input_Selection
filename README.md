MATLAB_Iterative_Input_Selection
==============================

"Update August, 2015: the toolbox can now implement the IIS algorithm for classification problems. The MATLAB_ExtraTrees toolbox has been accordingly updated as well". 

The MATLAB_Iterative_Input_Selection toolbox is a MatLab implementation of the Iterative Input Selection (IIS) algorithm proposed by Galelli and Castelletti (2013a). The underlying regression method adopted by the IIS algorithm is an ensemble of Extra-Trees (Geurts et al., 2006; Galelli and Castelletti, 2013b). The user is referred to the original publication for details regarding the IIS algorithm.  

The MATLAB_Iterative_Input_Selection toolbox requires the MATLAB_ExtraTrees toolbox, which can be found at https://github.com/rtaormina/MATLAB_ExtraTrees.


!!!!! ATTENTION !!!!!
A FASTER VERSION OF THE IIS ALGORITHM THAT EMPLOYS Extra-Trees WRITTEN IN C IS NOW AVAILABLE AT https://github.com/Critical-Infrastructure-Systems-Lab/Iterative_Input_Selection.



Contents of MATLAB_Iterative_Input_Selection :
* `script_example.m`: show how to use the available functions on a sample dataset (Friedman_dataset.txt).
* `crossvalidation_extra_tree_ensemble.m`: run a k-fold cross-validation for an ensemble of Extra-Trees.
* `input_ranking.m`: rank the input variables.
* `iterative_input_selection.m`: run the IIS algorithm.
* `visualize_inputSel.m`: visualize the results obtained with multiple runs of the IIS algorithm.
* `shuffle_data.m`: shuffle the observations of the sample dataset.
* `Rt2_fit.m`: compute the coefficient of determination R2.
* `Friedman_dataset.txt`: sample dataset, with 10 candidate inputs (first 10 columns) and 1 output (last column). The observations, arranged by rows, are 250.

Based on work from the following papers:

- Galelli, S., Humphrey, G.B., Maier, H.R., Castelletti, A., Dandy, G.C., Gibbs, M.S. (2014) An evaluation framework for input variable selection algorithms for environmental data-driven models (2014). Environmental Modelling & Software, 62, 33-51 ([Link to Paper](https://www.sciencedirect.com/science/article/abs/pii/S1364815214002394)).
- Galelli, S., and A. Castelletti (2013a), Tree-based iterative input variable selection for hydrological modeling, Water Resour. Res., 49(7), 4295-4310 ([Link to Paper](http://onlinelibrary.wiley.com/doi/10.1002/wrcr.20339/abstract)).
- Galelli, S., and A. Castelletti (2013b), Assessing the predictive capability of randomized tree-based ensembles in streamflow modelling, Hydrol. Earth Syst. Sci., 17, 2669-2684 ([Link to Paper](http://www.hydrol-earth-syst-sci.net/17/2669/2013/hess-17-2669-2013.html)).
- Geurts, P., D. Ernst, and L. Wehenkel (2006), Extremely randomized trees, Mach. Learn., 63(1), 3-42 ([Link to Paper](http://link.springer.com/article/10.1007/s10994-006-6226-1)).

Acknowledgements: to Dr. Matteo Giuliani (Politecnico di Milano), Riccardo Taormina (TU Delft), and Ahmad Alsahaf 
(Politecnico di Milano).

Copyright 2014 Stefano Galelli.

This file is part of MATLAB_Iterative_Input_Selection.

MATLAB_Iterative_Input_Selection is free software: you can redistribute
it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MATLAB_Iterative_Input_Selection. If not, see <http://www.gnu.org/licenses/>.

