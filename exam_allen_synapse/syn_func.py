import nest
import json
import numpy as np


from syn_param import *

nest.Install("allenmodule")


def main():
	h=0.1
	nest.ResetKernel()
	nest.SetKernelStatus({"resolution": h})

	spks=nest.Create('spike_generator',1,{'spike_times':STimes})
	
	pyr=nest.Create('iaf_psc_exp',80)
        pv=nest.Create('iaf_psc_exp',20)

	detector=nest.Create('weight_recorder')



if __name__=='__main__':
    main()


