!*************************************************************************
!*                                                                       *
!*     LMBM-Clust - Nonsmooth Optimization based Incremental             *
!*                  Clustering Software using LMBM (version 2.1)         *
!*                                                                       *
!*     by Napsu Karmitsa 2016 (last modified 11.11.2017).                *
!*     The code is partly based on clustering algorithms by              *
!*     Adil Bagirov 2015.                                                *
!*                                                                       *
!*     The software is free for academic teaching and research           *
!*     purposes but I ask you to refer the appropriate references        *
!*     given below, if you use it.                                       *
!*                                                                       *
!*************************************************************************
!*
!*
!*     Codes included:
!*
!*     lmbmclust.f95         - Mainprogram for clustering software
!*                             (this file).
!*     parameters.f95        - Parameters. Inludes modules:
!*                               - r_precision - Precision for reals,
!*                               - param - Parameters,
!*                               - exe_time - Execution time.
!*     initlmbmclust.f95     - initialization of clustering parameters and LMBM.
!*                             Includes modules:
!*                               - initclust - Initialization of parameters for 
!*                                             clustering.
!*                               - initlmbm  - Initialization of LMBM.
!*     lmbmclustmod.f95      - Subroutines for clustering software.
!*     functionmod.f95       - Computation of clustering function and gradient values.
!*     lmbm.f95              - LMBM - limited memory bundle method.
!*     objfun.f95            - Computation of the function and subgradients values.
!*     subpro.f95            - subprograms for LMBM.
!*
!*     Makefile              - makefile.
!*
!*     iris.txt              - sample data set
!*
!*     To use the software modify initlmbmclust.f95 as needed.
!*
!*
!*     References:
!*
!*     for LMBM-Clust:
!*
!*       N. Karmitsa, A. Bagirov and S. Taheri, "Clustering in Large Data Sets with 
!*       the Limited Memory Bundle Method", Pattern Recognition, Vol. 83, pp. 245-259, 
!*       2018.
!*
!*       N. Karmitsa, A. Bagirov and S. Taheri, "MSSC Clustering of Large Data using
!*       the Limited Memory Bundle Method" TUCS Technical Report, No. 1164, 
!*       Turku Centre for Computer Science, Turku, 2016.
!*
!*     for LMBM:
!*       N. Haarala, K. Miettinen, M.M. Mäkelä, "Globally Convergent Limited Memory 
!*       Bundle Method for Large-Scale Nonsmooth Optimization", Mathematical Programming, 
!*       Vol. 109, No. 1,pp. 181-205, 2007. DOI 10.1007/s10107-006-0728-2.
!*
!*       M. Haarala, K. Miettinen, M.M. Mäkelä, "New Limited Memory Bundle Method for 
!*       Large-Scale Nonsmooth Optimization", Optimization Methods and Software, 
!*       Vol. 19, No. 6, pp. 673-692, 2004. DOI 10.1080/10556780410001689225.
!*
!*       N. Karmitsa, "Numerical Methods for Large-Scale Nonsmooth Optimization" in 
!*       Big Data Optimization. A. Emrouznejad (eds.), Springer, 2016.
!*
!*     for NSO clustering:
!*       A. Bagirov, N. Karmitsa, M.M. Mäkelä, "Introduction to nonsmooth optimization:
!*       theory, practice and software", Springer, 2014.
!*
!*
!*     Acknowledgements:
!*
!*     The work was financially supported by the Academy of Finland (Project No. 289500).
!*
!*
