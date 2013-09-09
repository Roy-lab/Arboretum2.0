LFLAG = -lgsl -lgslcblas 
SRC = Framework.C  SpeciesDistManager.C	SpeciesClusterManager.C	Expert.C Matrix.C GeneExpManager.C MappedOrthogroupReader.C     MappedOrthogroup.C GeneMap.C GeneTreeManager.C GammaManager.C Gamma.C GeneTree.C GeneNameMapper.C NewickReader.C
LIBPATH = lib
INCLPATH1 =include
INCLPATH2 = common

CC=g++
CFLAGS = -g

condSpecLearner: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH1) -I $(INCLPATH2)  -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o incAncClust

