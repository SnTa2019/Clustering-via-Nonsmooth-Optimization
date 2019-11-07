#*************************************************************************
#*                                                                       *
#*     Makefile for nonsmooth clustering software with LMBM,             *
#*     (last modified 31.05.2016).                                       *
#*                                                                       *
#*************************************************************************
# 
#FF = gfortran -O2 -ffast-math
FF = gfortran -O2
OBJ = parameters.o initlmbmclust.o functionmod.o lmbmclustmod.o objfun.o subpro.o lmbm.o lmbmclust.o

all:	lmbmclust


parameters.o: parameters.f95
	$(FF) -c parameters.f95

initlmbmclust.o: r_precision.mod initlmbmclust.f95
	$(FF) -c initlmbmclust.f95

functionmod.o: r_precision.mod param.mod initclust.mod functionmod.f95
	$(FF) -c functionmod.f95

lmbmclustmod.o: r_precision.mod param.mod initclust.mod lmbmclustmod.f95
	$(FF) -c lmbmclustmod.f95

objfun.o: r_precision.mod functionmod.mod objfun.f95
	$(FF) -c objfun.f95

subpro.o: r_precision.mod param.mod subpro.f95
	$(FF) -c subpro.f95

lmbm.o: r_precision.mod param.mod initclust.mod initlmbm.mod exe_time.mod obj_fun.mod subpro.mod lmbm.f95
	$(FF) -c lmbm.f95

lmbmclust.o: r_precision.mod param.mod exe_time.mod initclust.mod initlmbm.mod clusteringmod.mod lmbm_mod.mod lmbmclust.f95
	$(FF) -c lmbmclust.f95

lmbmclust: $(OBJ)
	$(FF) -o lmbmclust $(OBJ)

clean:	
	rm lmbmclust $(OBJ) *.mod
	echo Clean done	