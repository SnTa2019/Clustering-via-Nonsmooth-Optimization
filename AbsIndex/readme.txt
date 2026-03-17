!****************************************************************************************
!*     AbsIndex: Absolute indices                                                       *
!*                                                                                      *
!*     The code is written by Prof. Adil Bagirov (March 2026)                           * 
!*     It is for determining compactness, separability and number of clusters           *                                         
!*                                                                                      *
!*     The software is free for academic teaching and research purposes,                *
!*     but please refer to the references given below, if you use it.                   *
!****************************************************************************************
!*     The user provides:  
!*     Data set  : (open(78,file='A1(k=20).txt',status='old',form='formatted'))
!*     nft       : Number of features, without label/output 
!*     nclust    : Maximum number of clusters
!*     eps1      : Value of epsilon (between 0.1 and 0.00001) 
!*
!*     Output files: 
!*     open(50,file='Compactness_clusters-A1(k=20).txt')
!*     open(51,file='Margins_between_clusters-A1(k=20).txt')
!*     open(49,file='Allres-A1(k=20).txt')
!****************************************************************************************
!*     Reference:
!*    "Bagirov, M., Aliguliyev, R., Sultanova, N., Taheri, S., (2025) 
!*     Absolute indices for determining compactness, separability and number of clusters" 
!*          
!*     Acknowledgements: 
!*     The work was financially supported by the Australian Government through the
!*     Australian Research Council's Discovery Projects funding scheme (Project No.
!*     DP190100580). 
!****************************************************************************************
