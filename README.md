# SynMMS

This library is developed to produce the results for the paper, 

Nat comm link

previous version,

[Konstantin Gavriljuk, Farid Ghasemalizadeh, Bruno Scocozza, Akhilesh P. Nandan, Hans Seidel, Malte Schmick, Aneta Koseska, Philippe I. H. Bastiaens. "A synthetic morphogenic perceptory cell](https://doi.org/10.1101/481887).

-------------------------
Requirements
-------------------------

The code has been tested successfully with versions of Python (at least 2.7) on Windows and macOS.

Using the framework
===================

The subfolders contains codes to generate the results obtained by numerical simulations in the main and supplementary figures.

SynMMs Reaction Diffusion 1D : Codes to generate kymographs in Fig. 1g, Fig. 8b and c, Supplementary Fig. 1h and i, Supplementary Fig. 9d-f\
MT Monte Carlo               : Codes to generate Supplementary Fig. 1g \
Aurora Reaction Diffusion 2D_2 variable + Recurrence Quantifications : Codes to generate Fig. 3i, Supplementary Fig. 4h\
Aurora Reaction Diffusion 2D_3variable                   : Codes to generate Fig. 4g

-----------------------------------------------------------

Each folder has a file having name ending with '_main.py'. It combines different objects (functions) defined in other files in the same folder. For example, in the folder SynMMS/SynMMs Reaction Diffusion 1D,
the 'reaction_diffusion1D_main.py' load the model equations and parameters from 'model.py' file. 
Sample code :
```python
pbc = PeriodicBoundaryConditions()
model = SynMMS()
light = TimeVaryingCue()
rd = ReactionDiffusion1D(model, pbc, light)
rd.simulate(save=True, animate=None)
```
One can either animate/save or do both simultaneosly. In the case when 'save=True', a subfolder needs to be created to save the files in '.npy' format.
'plot_kymograph.py' can be used to generate the kymographs from saved data.

The subfolder 'Aurora Reaction Diffusion 2D_2 variable + Recurrence Quantifications' has two main files. 'reaction_diffusion2D_main.py' to genrate the 2D pattern.
'spatial_recurrence_plot_main.py' to caculate Information Entropy value corresponding to an already saved spatial pattern. 

Excecuting codes via command prompt (eg. Anaconda command prompt) is recommended.

