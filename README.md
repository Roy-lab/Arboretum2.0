# Arboretum (version 2.4)

### How to compile the executable
---
type "make" in the arboretum directory.
```
  cd arboretum/
  make
```
The executable will be named arboretum.
<br>
<br>

### How to get the usage description
---
You can obtain the usage for this version with
```
  ./arboretum --help
```
Detail description for the input arguments is available at [README.txt](https://github.com/Roy-lab/Arboretum2.0/blob/Version2.4/README.txt)
<br>
<br>

### Differences / Key debugging log
---
- Addition of the usage of **fixed covariance** values
  + commented as "FIXED COVARIANCE" in the code.

* Repopulating the ancestral information at the **allspecies_clusterassign_lca_brk.txt** file
  + commented as "repop" in the several parts of the code.

- Allowance of very small cluster (previously the program ignored clusters < 10 genes)
  + commented as "neuron MODIFICATION" in SpeciesClusterManager.C

* DEBUG: gene tree misusing issue (found by Jon Ide)
  + commented as "JI mod" in GeneTreeManager.C

- DEBUG: fixation of "remove the whole OG" message (Feb 27 2021)
  + code changes GammaManager.C and Expert.C

* DEBUG: fixation of disappearing genes in the final result (Apr 4 2021)
  + code changes GammaManager.C

- DEBUG: fixation for mergeData issue (Aug 19 2021)
  + commented as "source spc only, 210819" in Framework.C
