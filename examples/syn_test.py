import nest
import numpy as np
nest.Install("allenmodule")

import pylab
from numpy import exp
import json
import sys
from CNS_param import *

target=str(sys.argv[1])
params=List[target]

'''
Synpase_internal_to_NEST_reserved={'tau_f_': 'tau_fac', 'tau_r_':'tau_rec','tau_fdr_':'tau_1', 'tau_r0_':'tau_2', 'tau_i_':'tau_eta',
  'p0_':'x', 'n_':'n', 'p_':'p', 'a_fdr':'y_0', 'a_D_':'alpha_1', 'a_i_':'alpha_2', 
  'S_':'S', 'p0bar_':'alpha', 'tau_d_':'z'}

This mapping is necessary to use synapse model without modifying NEST source codes (there are predefined variable names because NEST is written in C++. 
'''



  
def get_synaptic_w(STimes):

	h=0.1
	nest.ResetKernel()
	nest.SetKernelStatus({"resolution": h})

	spks=nest.Create('spike_generator',1,{'spike_times':STimes})
	pre=nest.Create('parrot_neuron',1)
	post=nest.Create('iaf_psc_exp',1)
	detector=nest.Create('weight_recorder')

       

	synapse_params={'weight':1000.,'tau_fac':params['Tau_f'], 'weight_recorder':detector[0],
                'tau_rec':params['Tau_r'],'tau_1':params['Tau_FDR'],'tau_2':params['Tau_r0'],
                'tau_eta':params['Tau_i'],'x':params['p0'],'n':1.0,
                'y_0':params['a_FDR'],'alpha_1':params['a_D'],'alpha_2':params['a_i'],
                'S':1.0,'alpha':params['p0bar'],'z':params['Tau_D'],'p':params['p0bar']}  

	nest.CopyModel("aibs_synapse", "syn",synapse_params)
	nest.Connect(spks,pre)
	nest.Connect(pre,post,syn_spec="syn")
	#nest.Connect(detector,pre)

	nest.Simulate(5000.0)

	weights=nest.GetStatus(detector,"events")[0]["weights"]

	weights=np.array(weights)
	weights=weights/weights[0]

	return weights

def get_times(freq,delay=250.):
	dt=1000./freq
	times=[]
        times.append(1.0) # currently ai_synapse disregards a 0-ms spike 
        for xin in xrange(7):
		times.append(times[-1]+dt)
	times.append(times[-1]+delay)
	for xin in xrange(3):
		times.append(times[-1]+dt)

	return times
	
       
Freqs=[10.,20.,50.,100.,200.]
Delays=[250] #.,500.,1000.,2000.,4000.]
colors={'10.0':'--mo','20.0':'--ro','50.0':'--go','100.0':'--bo','200.0':'--ko'}
for fin in Freqs:
	if fin==50.0:
		del_loop=Delays
	else:
		del_loop=[Delays[0]]
	for din in del_loop:
		ST=get_times(fin,din)

		weights=get_synaptic_w(ST)
		weights=np.array(weights)
		weights=weights/weights[0]
		pylab.semilogx(ST,weights,colors[str(fin)],label=str(fin)+'Hz stim')
pylab.xlabel('Time (ms)')
pylab.ylabel('Normalized PSP')
pylab.legend()


pylab.show()

