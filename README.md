# Short-term plasticity (STP) models for Neural Simulation Technology (NEST) platform. 
## Background
STP is one of synapses’ fundamental properties. Since Tsodyks and Markram proposed mathematical models (TM models) to account for observations, they have been considered the ‘standard model’, but more recent experimental studies suggested that more processes can be involved in regulating STP; see Hennig 2013 for the details. 
To characterize STP in the mouse primary visual cortex (V1) and study their functions in visual perception, we implemented generalized STP models (Hennig 2013) for NEST. This module (allenmodule) can be installed as an independent expansion without modifying the NEST source codes. 
## How to use the allenmodule
We have summarized the instructions provided by NEST developers. For more details, please visit [their github page] (https://nest.github.io/nest-simulator/extension_modules). directly. 
### Install
1.	Copy the folder (AllenModule) into a directory.
2.	Make a build directory (e.g., mmb).
3.	cd mmb
4.	cmake -Dwith-nest=${NEST_INSTALL_DIR}/bin/nest-config ../AllenModule
5.	make
6.	make install
### Using it in the NEST
For SLI, (allenmodule) Install.
For pynest, nest.Install("allenmodule"). 
## The summary of generalized STP model
This STP model implemented 5 gating processes (or dynamics), depression, facilitation, use-dependent replenishment, desensitization and slow-modulation of release probability. By setting the corresponding parameters to 0, each gating process can be turned on. The mathematical descriptions of the gating processes are provided in our CNS 2018 poster (available in the ‘examples’ folder). 
## Examples
The example folder has a script to replicate STPs which were presented at the CNS 2018 meeting. 
Usage: python syn_test.py pvalb_pvalb_1. 
Any 11 synapse types in the CNS_param.py can be chosen for now. The parameters have been obtained by our pipeline project analysis; see (Seeman et al., 2018) and CNS poster. Please note that the parameters are provided for illustration purposes, not as a part of our official data release. 
## Additional Notes
1.	During the compile, several warnings may appear. We are in the process of addressing them.
2.	The current module is compatible with NEST-2.16.0, and another version to be compatible with NEST-2.14.0 will be added to the repository. 
## Reference
Hennig, M.H. (2013). Theoretical models of synaptic short term plasticity. Front. Comput. Neurosci. 7, 1–10.

Seeman, S.C., Campagnola, L., Davoudian, P.A., Hoggarth, A., Hage, T.A., Bosma-Moody, A., Baker, C.A., Lee, J.H., Mihalas, S., Teeter, C., et al. (2018b). Sparse recurrent excitatory connectivity in the microcircuit of the adult mouse and human cortex. Elife 7, 1–27.

Tsodyks, M., Uziel, A., and Markram, H. (2000). Synchrony Generation in Recurrent Networks with Frequency-Dependent Synapses. 20, 1–5.

