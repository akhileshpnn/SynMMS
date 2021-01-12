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

SynMMs ReactionDiffusion1D : 


-----------------------------------------------------------

In general, for running simulations of a customized model, it needs to be defined first, similar to the way the other models are defined (labels, df_model function, etc).
Same holds for the experiment type (pulsed, sustained, step-wise increasing/decreasing, etc.).
Then the two objects are combined in a model_simulation object that can run deterministic/stochastic simulations and plot results.
Sample code:

```matlab
model = models.simple_dnf_model;  % create model object
model.par.g1 = 2.957;  % set regime of operation by setting the parameter(s); in this case the criticality regime
mpe = experiments.multi_pulse_experiment(2);  % create experiment object - multiple pulses at regular intervals
ms = model_simulation(model, mpe);  % combine the two in a model_simulation object 
ms.simulate();  % run deterministic simulations
ms.plot_fraction_active();  % plot the results
```

One can also animate a stochastic simulation, where the state space is updated, resulting in temporal changes:

```matlab
ms.stochastic_simulation();
ms.animate();
```

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
