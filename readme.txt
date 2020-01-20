README

-Package Description

	-DCS-LTP package for running simulations of detailed biophysical neurons with the NEURON simulation environment as a python module.  The package is intended for simulating synaptic inputs onto pyramidal neurons during application of extracellular electric fields

-General Workflow
	
	-run_control.py is used to specify parameters and procdures for specific experiments
		-see run_control.ExperimentsParallel for running parallel simulations

	-neurons.py contains methods for creating cells.  These are python objects that contain a collection hoc sections and their inserted mechanisms

	-run.py contains methods for running simulations and saving data

	-param.py contains parameter dictionaries and methods for updating parameters for specific simulations

	-stims.py contains methods for creating stimulation objects (e.g. extracellular field, synaptic stimulation) based on parmaters specified paramter dictionaries

	-analysis.py contains functions for post-simulation analysis

	-analysis_exp.py procedures for analyzing data from specific experiments (those specified by run_control)

