--------------------
------ vault storage
--------------------


-- profondeur
!DEEP=PLTUB*NRANG
!NDEP=NRANG*NCELL+1

-- largeur
!WIDE=PTTUB*NLIGN
!NWID=NLIGN*NCELT+1

-- hauteur
!HIGHT=PVTUB*NZTUB
!NHIG=NZTUB+1

-- jeu amont 1
!JAMU=-1*GAPGRE
!NAMU=NELGRE+1

-- jeu aval 1
!JAVU=DEEP+GAPREG
!NAVU=NELREG+1

-- grille amont
!AMON=JAMU-ABAMON
!NAMO=5

-- grille aval 
!AVAL=JAVU+ABAMON
!NAVA=5

-- epaisseur des cheminees
!SCHEM=EPSCHEM
!NCHEM=6

-- jeu amont 2
!JAMD=AMON-GAPCHEG
!NJEU2=NELCHEG+1

-- jeu aval 2
!JAVD=AVAL+GAPGECH
!NJEU3=NELGECH+1

-- mur amont 
!MAMO=JAMD-SCHEM

-- mur aval
!MAVA=JAVD+SCHEM

-- hauteur cheminee d alimentation
!HALIM=HCHALIM
!NALIM=11

-- cote cheminee
!YMCH=WIDE*(1.-1./RCONVER)*0.5
!YPCH=WIDE*(1.+1./RCONVER)*0.5

-- hauteur convergent
!HCONV=HCONVER
!NCONV=11

-- hauteur cheminee d evacuation
!HEVAC=HCHEVAC
!NEVAC=13

-- reference point
!RPOI=11

-- reference surface
!RIN=5
!ROUT=1
!RBOT=10
!RTOP=11 
!RSYM=12
!RMUR=13

-- references cellules
!RALI=0
!REVA=0
!RCAL=3
!RCEV=6
!RAMO=2
!RAVA=4
!RSTO=8
!RJM=7
!RJV=9

POINTS
     2   CART    0                      $ IMPRE TYPCOO ICONST
$ NUM  REF      X       Y       Z
  1    RPOI     0       0       0
  2    RPOI     WIDE    0       0
  3    RPOI     WIDE    DEEP    0      
  4    RPOI     0       DEEP    0
  5    RPOI     WIDE    AMON    0
  6    RPOI     0       AMON    0
  7    RPOI     WIDE    AVAL    0
  8    RPOI     0       AVAL    0
  11   RPOI     0       0       HIGHT  
  12   RPOI     WIDE    0       HIGHT
  13   RPOI     WIDE    DEEP    HIGHT  
  14   RPOI     0       DEEP    HIGHT
  15   RPOI     WIDE    AMON    HIGHT
  16   RPOI     0       AMON    HIGHT
  17   RPOI     WIDE    AVAL    HIGHT
  18   RPOI     0       AVAL    HIGHT
  21   RPOI     WIDE    MAMO    HIGHT
  22   RPOI     0       MAMO    HIGHT
  23   RPOI     WIDE    MAVA    HIGHT
  24   RPOI     0       MAVA    HIGHT
  25   RPOI     WIDE    MAMO    0
  26   RPOI     0       MAMO    0
  27   RPOI     WIDE    MAVA    0
  28   RPOI     0       MAVA    0
  
  44   RPOI     YMCH       MAVA    HEVAC
  43   RPOI     YPCH    MAVA    HEVAC
  37   RPOI     YPCH    JAVD    HEVAC
  38   RPOI     YMCH       JAVD    HEVAC
  36   RPOI     YMCH       JAMD    HALIM
  35   RPOI     YPCH    JAMD    HALIM
  41   RPOI     YPCH    MAMO    HALIM
  42   RPOI     YMCH       MAMO    HALIM
  
  84   RPOI     YMCH       MAVA    HCONV
  83   RPOI     YPCH    MAVA    HCONV
  77   RPOI     YPCH    JAVD    HCONV
  78   RPOI     YMCH       JAVD    HCONV
  76   RPOI     YMCH       JAMD    HCONV
  75   RPOI     YPCH    JAMD    HCONV
  81   RPOI     YPCH    MAMO    HCONV
  82   RPOI     YMCH       MAMO    HCONV
  
  51   RPOI     WIDE    JAMD    HIGHT
  52   RPOI     0       JAMD    HIGHT
  53   RPOI     WIDE    JAVD    HIGHT
  54   RPOI     0       JAVD    HIGHT
  
  55   RPOI     WIDE    JAMD    0
  56   RPOI     0       JAMD    0
  57   RPOI     WIDE    JAVD    0
  58   RPOI     0       JAVD    0
  
  59   RPOI     WIDE    JAMU    HIGHT
  60   RPOI     0       JAMU    HIGHT
  61   RPOI     WIDE    JAVU    HIGHT
  62   RPOI     0       JAVU    HIGHT
  
  63   RPOI     WIDE    JAMU    0
  64   RPOI     0       JAMU    0
  65   RPOI     WIDE    JAVU    0
  66   RPOI     0       JAVU    0  
 0                                      $ FIN DES POINTS

LIGNES
     2     0                                $ IMPRE ICONST
$ NO NOEUDS NP1  NP2 NREF NFFRON RAISON
  1   NWID  1    2    RBOT    0   1
  2   NDEP  2    3    RSYM    0   1
  3   NWID  3    4    RBOT    0   1
  4   NDEP  4    1    RSYM    0   1
  
  11  NWID  11   12   RTOP    0   1
  12  NDEP  12   13   RSYM    0   1
  13  NWID  13   14   RTOP    0   1
  14  NDEP  14   11   RSYM    0   1

  5   NAMO  63   5    RSYM    0   1
  6   NAMO  64   6    RSYM    0   1 
  15  NAMO  59   15   RSYM    0   1
  16  NAMO  60   16   RSYM    0   1 
  
  7   NAVA  65   7    RSYM    0   1
  8   NAVA  66   8    RSYM    0   1 
  17  NAVA  61   17   RSYM    0   1
  18  NAVA  62   18   RSYM    0   1 

  21  NHIG  1    11   RSYM    0   1 
  22  NHIG  2    12   RSYM    0   1 
  23  NHIG  3    13   RSYM    0   1 
  24  NHIG  4    14   RSYM    0   1 
  25  NHIG  5    15   RIN    0   1 
  26  NHIG  6    16   RIN    0   1 
  27  NHIG  7    17   ROUT    0   1 
  28  NHIG  8    18   ROUT    0   1 
  
  9   NWID  5    6    RMUR    0   1
  10  NWID  15  16    RMUR    0   1
  19  NWID  7    8    RMUR    0   1
  20  NWID  17  18    RMUR    0   1

  31  NCHEM  25  55    RSYM     0   1
  32  NJEU2  51  15    RSYM     0   1
  51  NCHEM  21  51    RSYM     0   1
  33  NCHEM  26  56    RSYM     0   1
  34  NJEU2  52  16    RSYM     0   1
  52  NCHEM  22  52    RSYM     0   1
  35  NCHEM  27  57    RSYM     0   1
  36  NJEU3  53  17    RSYM     0   1
  53  NCHEM  23  53    RSYM     0   1
  37  NCHEM  28  58    RSYM     0   1
  38  NJEU3  54  18    RSYM     0   1    
  54  NCHEM  24  54    RSYM     0   1    

  39  NWID  25  26    RMUR    0   1
  40  NWID  27  28    RMUR    0   1
  29  NWID  21  22    RMUR    0   1
  30  NWID  23  24    RMUR    0   1

  41  NWID  44  43    RMUR    0   1
  42  NCHEM  43  37    RSYM   0   1  
  43  NWID  37  38    RMUR    0   1
  44  NCHEM  38  44    RSYM   0   1  

  81  NWID  84  83    RMUR    0   1
  82  NCHEM  83  77    RSYM   0   1  
  83  NWID  77  78    RMUR    0   1
  84  NCHEM  78  84    RSYM   0   1  

  45  NWID  42  41    RMUR    0   1
  46  NCHEM  41  35    RSYM   0   1  
  47  NWID  35  36    RMUR    0   1
  48  NCHEM  36  42    RSYM   0   1 

  85  NWID  82  81    RMUR    0   1
  86  NCHEM  81  75    RSYM   0   1  
  87  NWID  75  76    RMUR    0   1
  88  NCHEM  76  82    RSYM   0   1 
  
  55  NWID  52  51    RMUR    0   1
  56  NWID  53  54    RMUR    0   1  
  57  NWID  55  56    RMUR    0   1
  58  NWID  57  58    RMUR    0   1 

  59  NJEU2   55  5    RSYM     0   1  
  60  NJEU2   56  6    RSYM     0   1  
  61  NJEU3   57  7    RSYM     0   1  
  62  NJEU3   58  8    RSYM     0   1  
   
  63  NAMU   12  59    RSYM     0   1  
  64  NAMU   11  60    RSYM     0   1  
  65  NAVU   61  13    RSYM     0   1  
  66  NAVU   62  14    RSYM     0   1  
   
  67  NAMU    2  63    RSYM     0   1  
  68  NAMU    1  64    RSYM     0   1  
  69  NAVU    3  65    RSYM     0   1  
  70  NAVU    4  66    RSYM     0   1  

  71  NWID  60  59    RMUR    0   1
  72  NWID  61  62    RMUR    0   1  
  73  NWID  64  63    RMUR    0   1
  74  NWID  65  66    RMUR    0   1

0                                        $ FIN DES LIGNES

  
  



