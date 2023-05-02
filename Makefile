LFLAG = -lgsl -lgslcblas 
SRC = Framework.C  SpeciesDistManager.C	SpeciesClusterManager.C	Expert.C Matrix.C GeneExpManager.C MappedOrthogroupReader.C     MappedOrthogroup.C GeneMap.C GeneTreeManager.C GammaManager.C Gamma.C GeneTree.C GeneNameMapper.C NewickReader.C BioNetwork.C  BioNode.C Gene.C  GeneManager.C  Graph.C  Interaction.C  InteractionManager.C  Node.C  Path.C  Protein.C  ProteinManager.C BFGSWrapper.C ClusterManager.C Kmeans.C  ExpertL.C MotifManager.C MotifRegressor.C Randomizer.C TestData.C  Distance.C HyperGeomPval.C common/Variable.C common/VariableManager.C common/Evidence.C common/EvidenceManager.C common/Error.C

LIBPATH = lib
INCLPATH1 =include
INCLPATH2 = common

CC=g++
CFLAGS = -g

condSpecLearner: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH1) -I $(INCLPATH2)  -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o arboretum

