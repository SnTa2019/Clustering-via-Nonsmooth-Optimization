!*************************************************************************
!*                                                                       *
!*     DC Optimization based incremental clustering using L2 norm        *
!*                                                                       *
!*     Original code with DC-Clust algorithm by Adil Bagirov 2015.       *                     
!*                                                                       *
!*     Fortran 95 version and optimization with DCD-Bundle method        *
!*     by Napsu Karmitsa 2016 (last modified 24.02.2016).                *
!*                                                                       *
!*     The software includes both DCD-Bundle and DC-Clust algorithms.    * 
!*     You can switch between the algorithms by changing parameter       *
!*     “optmet” in initclust.f95. The software is free for academic      *
!*     teaching and research purposes but I ask you to refer the         *
!*     the references given below, if you use it.                        *
!*                                                                       *
!*************************************************************************
!*
!*
!*     Codes included:
!*
!*     dcclustering.f95      - Mainprogram for clustering software.
!*     parameters.f95        - Parameters. Inludes modules:
!*                               - r_precision - Precision for reals,
!*                               - param - Parameters,
!*                               - exe_time - Execution time.
!*     initclust.f95         - Initialization of clustering parameters and
!*                             optimization methods. Includes modules:
!*                               - initclust - Initialization of parameters for
!*                                 clustering.
!*                               - initdcclust - Initialization of DC-Clust -solver.
!*                               - initdcdb - Initialization of DCDB -solver.
!*     clusteringmod.f95     - Subroutines for clustering software.
!*     functionmod.f95       - Computation of function and (sub)gradients values for
!*                             clustering software.
!*     dcclust_method.f95    - DC-Clust method.
!*     dcdb.f95              - DCD-Bundle method.
!*     dcfun.f95             - Computation of the function and (sub)gradients
!*                             values for DCD-Bundle.
!*     subpro.f95            - subprograms for DCD-Bundle.
!*
!*     iris.txt              - sample data set
!*
!*     Makefile              - makefile.
!*
!*
!*     To use the software modify initclust.f95 (and dcclustering.f95) as needed
!*     and give your data set similar to iris.txt here.
!*
!*
!*     References:
!*
!*     "New diagonal bundle method for clustering problems in very large data sets",
!*     Napsu Karmitsa, Adil Bagirov and Sona Taheri, European Journal of Operational 
!*     Research, Vol. 263, No. 2, pp. 367-379, 2017.
!*
!*     "Diagonal Bundle Method for Solving the Minimum Sum-of-Squares Clustering 
!*     Problems”, Napsu Karmitsa, Adil Bagirov, and Sona Taheri, TUCS Technical 
!*     Report, No. 1156, Turku Centre for Computer Science, Turku, 2016.
!*
!*     "Nonsmooth DC programming approach to the minimum sum-of-squares 
!*     clustering problems", Adil Bagirov, Sona Taheri and Julien Ugon, 
!*     Pattern Recognition, Vol. 53, pp. 12–24, 2016.
!*
!*
!*     Acknowledgements: 
!*
!*     The work was financially supported by the Academy of Finland (Project No. 289500)
!*     and Australian Research Counsil's Discovery Projects funding scheme (Project No.
!*     DP140103213).
!*
!*************************************************************************

