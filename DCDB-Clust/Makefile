#*************************************************************************
#*                                                                       *
#*     Makefile for clustering software with DC-Clust and DCDB-methods   *
#*     (last modified 24.02.2016).                                       *
#*                                                                       *
#*************************************************************************
# 
FF = gfortran -O2
OBJ = parameters.o initclust.o functionmod.o clusteringmod.o dcclust_method.o dcfun.o subpro.o dcdb.o dcclustering.o

all:	dcclustering


parameters.o: parameters.f95
	$(FF) -c parameters.f95

initclust.o: r_precision.mod initclust.f95
	$(FF) -c initclust.f95

functionmod.o: r_precision.mod param.mod initclust.mod functionmod.f95
	$(FF) -c functionmod.f95

clusteringmod.o: r_precision.mod param.mod initclust.mod clusteringmod.f95
	$(FF) -c clusteringmod.f95

dcclust_method.o: r_precision.mod param.mod initclust.mod initdcclust.mod functionmod.mod dcclust_method.f95
	$(FF) -c dcclust_method.f95

dcfun.o: r_precision.mod functionmod.mod dcfun.f95
	$(FF) -c dcfun.f95

subpro.o: r_precision.mod param.mod subpro.f95
	$(FF) -c subpro.f95

dcdb.o: r_precision.mod param.mod initdcdb.mod exe_time.mod dc_fun.mod subpro.mod dcdb.f95
	$(FF) -c dcdb.f95

dcclustering.o: r_precision.mod param.mod initclust.mod functionmod.mod clusteringmod.mod dcclust_method.mod initdcclust.mod initdcdb.mod exe_time.mod subpro.mod dcdb_mod.mod dcclustering.f95
	$(FF) -c dcclustering.f95

dcclustering: $(OBJ)
	$(FF) -o dcclustering $(OBJ)

clean:	
	rm dcclust $(OBJ) *.mod *~