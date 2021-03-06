The .zip archive dowloadable from this repository includes one directory, named Arboretum, with the C++ code for the updated version of the Arboretum program.

To compile the executable, type "make" in the arboretum directory. The executable will be named arboretum. The code should take not more than 1 minute to compile in my experience.

You can obtain the usage for this version with ./arboretum --help.

The main difference in this update of arboretum is that 1) the command arguments are now formally managed in the code by the getopt package, meaning each argument in the command will be proceeded with a "-x X" flag. These arguments are as follows:

Required options:
-s              File listing species in the orthology.
-e              File listing the orthology relationships of genes across species.
-k              Number of clusters in each species.
-t              Species tree file of species which are represented by the input data in the analysis
-c              The input cluster assignment and expression data for each species in the analysis.
This file should be of the form: species <tab> cluster_assign_file <tab> expression_data_file
-r              Option to randomize input cluster assignments rseed|none|yes.
-o              Output directory.
-m              Defines the mode in which the algorithm is to be used; learn|generate|visualize|crossvalidate are the options.
-b              A well annotated species to which gene names of other species will be mapped in the *_clusterassign.txt output.
-i              Initialization method for transition probabilities for cluster membership across species, uniform|branchlength.
-p              Initial transition probability values, either a single value or a file defining values every specie tree branch.
The init--type setting is for specifying how the diagonal cluster membership transition probabilities across species will be initialized. If inittype is uniform then the -p argument should be the default initial value for diagonal transition probabilities on all branches of the species tree. If the branchlength option is used, then non-uniform transition probabilities will be taken from a file set by the -p option.


The following are non-required options.
-g              Directory containing gene trees, which can be used to directly define the GeneTreeManager gene mappings.
-d,--conv-thresh                This is the change in the likelihood score between iterations that defines the convergence of the algorithm.
-n              The maximum number of iterations if convergence condition is not met.
-v              This option turns on a cross validation test of the clustering; and defines the number of partitions to use.
-w              A true|false option to overwrite the output directory if it already exists, the default is false.
-l		A true|false option to run the last EM optimization step in SpeciesClusterManager::dumpAllInferredClusterAssignments(). Default is True.
-f		A true|false option to run the merged clustering of the data to generate input clusterassignments. Default is False.
 When this option is used the input configuration (-c) file should be in the form: species <tab> expression_data_file <endl>
-u		A true|false option for updating the cluster means from the source species. Note this should only be attempted if the species all have the same number of data points. Otherwise it is not stringently requred that all species do have the same number of data points per gene, but all genes represented within each species must of course have the same number of data points. 
-h,--help		Reports this usage information.
A note on prediction mode: this is most commonly used when you have an arboretum clustering result to which you have applied a reordering, and hence need to re-infer assignments for ancestral assignments, or convesely when you want to include additional ortho-groups in your analysis and infer the cluster assignments for new genes based on the existing model.

The short tutorial on usisng this updated code is as follows:

The initial (required) options are those that you are already used to, the formats of the inputs and the fundamentals of the Arboretum algorithm are not changed.

The -o option is the output directory that you want Arboretum to write out results to. This directory will be generated if it doesn't exist, and if it exists, you have the option to overwrite an existing directory with the new -w true|false option. This is intended to protect anyone from overwritting results accidentally.

The -g argument option will point to the gene tree directory that you want to use as the input. This is new, but it's one of the most useful options that have been added to the command line.

The "-f true" argument option is the new option that will allow you to tell the program to merged the data across species, and to generate the initial species clusterassignments from GMM clustering of the meregd data. Three sets of outputs will be produced:

        1. In the Arboretum output directory a mergedData.txt file will appear, which is the data merged across all species. There will be several companion files that begin with "mergedData" that format this information as input for the GMM clustering process.

        2. In the output directory there will be the mergedClustering sub-directory, which contains the GMM clustering results fo the merged data.

        3. In the Arboretum output directory the initial clusterassignments for each species will appear with the name <species>_initial_clusterassign.txt.

Note that if you use "-f true", the -c input file should be of the form "species\texpression data file\n", instead of "species\tinput clusterassignment file\texpression data file\n" with "-1 false" (default). In short there will only two columns with with the species names and expression data files, none for cluster assignment file. A typical usage in a case of running the mereged clustering within Arboretum is as follows:

./arboretum -s ${SPECIESORDER} -e ${OGIDS} -k ${k} -t ${TREE} -c ${CONF} -r none -o ${DIR} -m learn -b ${BSPECIES} -i uniform -p 0.80 -g ${GENETREES} -n 100 -f true -w true

To use the example files from the Arboretum Website at http://pages.discovery.wisc.edu/~sroy/arboretum/example_inputs.tgz, an example command could be:

./arboretum -s specorder_allclade.txt -e OGid_members.txt -k 5 -t species_prob_heat8spec.txt -c cluster_conf.txt - r none -o result_dir -m learn -b Scer -i uniform -p 0.8 -g data/TREES -n 100 -f false -w true

This will create the output directory result_dir is it already exists, or otherwise write over it with the -w true option. Details of the file formats can also be studied from these examples.

==============================================================================================================================

A note about GSL installation and compilation of Arboretum:

the GSL library is a dependency of Arboretum. In some UNIX systems it will be readilly available for compilation in the lib and include directory of this repository for a working version of this code. However in general, that may not automatically work on all platfoms and you can install the GSL library in a given system (tested in OS X) with the following steps:

wget http://mirror.clarkson.edu/gnu/gsl/gsl-latest.tar.gz
tar xvfz gsl-latest.tar.gz 
cd gsl-2.6/
./configure 
make
make install

This will place an installation of the GSL library at the /usr/local/lib and /use/local/include, respectively (again, tested in OS X). This means that if you're not using the provided GSL lib and include files, you may need to update LIBPATH and INCLPATH1 in the Makefile for the compilation to work. 

