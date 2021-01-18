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
One can either animate/save or do both simultaneosly. In the case when 'save=True', a subfolder needs to be defined to save the files in '.npy' format.
'plot_kymo.py' can be used to generate the kymographs from saved data.

Different types of experiments can be simulated, such as streams of pulses for example:

```matlab
experiments.stream_pulse_experiment(n_pulses, lambda, t_total_min)
```

or different models can be used:

```matlab
model = models.model_2comp();
```

To run a stochastic single-molecule simulation simply define an animation and plot the execution:

```matlab
absa = abs_animation('ics_hl',1,'g1',4.95); % ics_hl=1 for high Ra initial conditions; also set g1 parameter
absa.plot_all();
```

Each run simulates 2s - add more 'n_parts' if you wish to simulate for longer, rather than increasing the duration in the agent_based_simulation class:

```matlab
absa = abs_animation('ics_hl',1,'g1',4.95,'rep',2,'n_parts',2); % different repetition (rep=2); simulate 4s (2x2s)
absa.plot_all();
```

Different repetitions with the same settings can be simulated as shown above. If the repetition exists already, it will be loaded instead.
Simulations from the paper can only be loaded with versions of Matlab ≥ R2018b.
