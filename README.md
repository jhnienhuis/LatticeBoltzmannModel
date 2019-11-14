# lattice_boltzmann_model
Fluid flow model using the lattice boltzmann concepts
see, Nienhuis et al., JGR 2014, https://doi.org/10.1002/2014JF003158

Central function is LBM.m

LBM.m starts a model run using settings described in InitializeStruct.m and InitializeLBM.m
LBM.m then uses these settings to run the model in RunLBM.m

The user can initialize multiple model runs using LBM_ParameterRun.m
LBM_ParameterRun.m uses a list of settings such as bed geometries or boundary conditions to start a LBM.m model run for each setting.
