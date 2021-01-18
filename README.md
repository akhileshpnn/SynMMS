# SynMMS

This library is developed to produce the results for the paper 

Nat comm link

previous version on

[Konstantin Gavriljuk, Farid Ghasemalizadeh, Bruno Scocozza, Akhilesh P. Nandan, Hans Seidel, Malte Schmick, Aneta Koseska, Philippe I. H. Bastiaens. "A synthetic morphogenic perceptory cell](https://doi.org/10.1101/481887).

-------------------------
Requirements
-------------------------

The code has been tested successfully with versions of Python (at least 2.7) on Windows and macOS.

Using the framework
===================

The subfolders contains codes that generates theoretical results in the main and supplimentary figures. 

SynMMs ReactionDiffusion1D : Codes to generate and plot kymographs in Fig. 1, Fig. 8b and c, Supplimentary Fig. 1h and i, Supplimentary Fig. 9d-f
MT Monte Carlo             : Codes to generate Supplementary Fig. 1g 
Aurora ReactionDiffusion2D + RecurrenceQuantifications : Codes to generate Fig. 3i, Supplementary Fig. 4h
Aurora ReactionDiffusion2D_3variable                   : Codes to generate Fig. 4g

-----------------------------------------------------------

Each folder has a file with filename ending with '_main.py'. It combines different objects defined in other files. For example, in the folder SynMMS/SynMMs ReactionDiffusion1D,
the 'reaction_diffusion1D_main.py' load the model equations and parameters from 'model.py' file. 
Sample code :
```python
pbc = PeriodicBoundaryConditions()
model = SynMMS()
pheromone = TimeVaryingPheromone()
rd = ReactionDiffusion1D(model, pbc, pheromone)
rd.simulate(save=True, animate=None)
```
One can either animate/save or do both simultaneosly. In the case when 'save=True', a subfolder needs to be created to save the files in '.npy' format.
'plot_kymo.py' can be used to generate the kymographs from saved data.

The subfolder 'Aurora ReactionDiffusion2D + RecurrenceQuantifications' has two main files. 'reaction_diffusion2D_main.py' to genrate the 2D pattern.
'spatial_recurrence_plot_quantifications.py' to caculat Information Entropy value corresponding to an already saved pattern. 

Excecuting codes via command prompt (eg. Anaconda cmmand Prompt) is recommended.
