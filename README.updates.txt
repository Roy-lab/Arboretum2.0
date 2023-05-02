original
        SpeciesClusterManager.C:	L304
		if(sample<10)

_4neuron
	SpeciesClusterManager.C:	L304                           
		if(sample<3)    // neuron MODIFICATION
_4neuron_sparse
        SpeciesClusterManager.C:	L304
                if(sample<1)    // neuron MODIFICATION


Fixation of "remove the whole OG" message (Feb 27 2021)
#######################################################
[Yesterday 4:56 PM] SUSHMITA ROY
    hi Junha Shin I was able to get over the message of "remove the whole OG". That happens when the probability values become 0 and there is nothing to compute for the probabilities. I put a check so if things are really small, just set to a small value like 1e-300.. 
    the code is in /mnt/dv/wid/projects5/Roy-singlecell/sr_work/arboretum_debug/arboretum-repop The change is GammaManger.C and Expert.C
    GammaManager.C only estimateNonLeafAlpha change is to look at

