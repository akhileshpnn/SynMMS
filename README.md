# SynMMS

This library is developed to produce the results for the paper, 

Accepted in Nature Communication : Link will be added soon

previous version,

[Konstantin Gavriljuk, Farid Ghasemalizadeh, Bruno Scocozza, Akhilesh P. Nandan, Hans Seidel, Malte Schmick, Aneta Koseska, Philippe I. H. Bastiaens. "A synthetic morphogenic perceptory cell](https://doi.org/10.1101/481887).


Python Code
===================

The code has been tested successfully with versions of Python (at least 2.7) on Windows and macOS.

-------------------------
Using the framework
-------------------------

The subfolders contains codes to generate the results obtained by numerical simulations in the main and supplementary figures.

Reaction Diffusion 1D : Kymographs in Fig. 1g, 8b and c, Supplementary Fig. 1h, i, 9d-g\
MT Monte Carlo               : Supplementary Fig. 1g \
Aurora Reaction Diffusion 2D_2 variable + Recurrence Quantifications : Fig. 3i, Supplementary Fig. 4h. The spatial recurrence quantification code is used for the quantification shown in Fig. 3f \
Aurora Reaction Diffusion 2D_3 variable                   : Supplementary Fig. 4g

-----------------------------------------------------------

Each folder has a file with name ending '_main.py'. This file combines different objects (functions) defined in linked files in the same folder. For example, in the folder SynMMS/Reaction Diffusion 1D,
the 'reaction_diffusion1D_main.py' loads the model equations and parameters from 'model.py' file. 
Sample code :
```python
pbc = PeriodicBoundaryConditions()
model = SynMMS()
light = TimeVaryingCue()
rd = ReactionDiffusion1D(model, pbc, light)
rd.simulate(save=True, animate=None)
```
One can either animate/save or do both simultaneosly. In the case when 'save=True', a subfolder needs to be created to save the files in '.npy' format.
'plot_kymograph.py' can be used to generate the kymographs from 'saved data' folder.

The subfolder 'Aurora Reaction Diffusion 2D_2 variable + Recurrence Quantifications' has two main files. 'reaction_diffusion2D_main.py' to genrate the 2D pattern.
'spatial_recurrence_plot_main.py' to caculate Information Entropy value corresponding to an already saved spatial pattern. 

Excecuting codes via command prompt (eg. Anaconda command prompt) is recommended.


Matlab Code
===================

The folder "CellularAutomata_RDSim" contains Matlab R2018a scripts and c++ code for compilation in Matlab to simulate the reaction diffusion system of the tubulin-stathmin association/dissociation cycle and the stathmin phosphorylation cycle in 1D and 2D (Fig 4f, Supplementary Fig. 6i-k and 8i). This code has been successfully tested with Matlab R2018a - R2020a under Mac OS X and Windows 10. It requires a c++ compiler (xCode for Mac OS X, Visual Studio 2019 for Windows) to be configured to compile and execute c++ code from within Matlab.

Please refer to the respective Matlab scripts on how to use the simulations to generate figure panels.
