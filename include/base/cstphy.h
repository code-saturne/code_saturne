!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

!                              cstphy.h
!===============================================================================

! KELVIN

!       TKELVI       --> =  273,15
!       TKELVN       --> = -273,15

double precision   tkelvi            , tkelvn
parameter (        tkelvi =  273.15d0, tkelvn = -273.15d0)

! CALORIES

!       1 cal = XCAL2J J

double precision   xcal2j
parameter (        xcal2j = 4.1855d0)

! STEPHAN BOLTZMANN

double precision   stephn
parameter (        stephn = 5.6703d-8)

! GRAVITE

double precision   gx,gy,gz
common / rgravi /  gx,gy,gz

! CONSTANTES PHYSIQUES DU FLUIDE
!   IXYZP0 : INDICATEUR DE REMPLISSAGE DE XYZP0
!   RO0    : MASSE VOLUMIQUE    DE REFERENCE
!   VISCL0 : VISCOSITE          DE REFERENCE
!   P0     : PRESSION TOTALE    DE REFERENCE
!   PRED0  : PRESSION REDUITE   DE REFERENCE
!   XYZP0  : POSITION PRESSION  DE REFERENCE
!   T0     : TEMPERATURE        DE REFERENCE
!   CP0    : CHALEUR SPECIFIQUE DE REFERENCE

integer           ixyzp0(nphsmx)
common / icstfl / ixyzp0
double precision  ro0(nphsmx)    , viscl0(nphsmx),                &
                  p0 (nphsmx)    , pred0 (nphsmx),                &
                  xyzp0(3,nphsmx), t0    (nphsmx),                &
                  cp0(nphsmx)
common / rcstfl / ro0        , viscl0        ,                    &
                  p0         , pred0         ,                    &
                  xyzp0      , t0            ,                    &
                  cp0

! TURBULENCE
!   IVISLS = 0 : VISCOSITE LAMINAIRE CONSTANTE = VISLS0
!   XKAPPA : CST DE KARMAN (~0.42)
!   CSTLOG : CST DE LA LOI LOG: 1/XKAPPA*LOG(YPLUS) + CSTLOG (~5.2)
!   YPLULI : YPLUS LIMITE 1./XKAPPA OU 10.88 SI IDEUCH=2
!   *POW   : COEFF WERNER AND WENGLE
!   CMU025 = CMU**0.25
!   CE1, CE2, SIGMAK, SIGMAE :
!            CONSTANTES DU K-EPSILON
!   C*RIJ* : CONSTANTES DU Rij-EPSILON STANDARD (LRR)
!   CSSG*  : CONSTANTES SPECIFIQUES DU RIJ-EPSILON SSG
!   CV2F*  : CONSTANTES SPECIFIQUES DU V2F PHI-MODEL
!   CKW*   : CONSTANTES SPECIFIQUES DU K-OMEGA SST
!            (SK=SIGMA_K, SW=SIGMA_W, BT=BETA, GM=GAMMA)
!   ALMAX  : ECHELLE DE LONGUEUR TURBULENTE
!   UREF   : VITESSE DE REFERENCE
!   XLOMLG : LONGUEUR POUR LONGUEUR DE MELANGE
!   XLESFL, ALES, BLES
!       DELTA = XLESFL * (ALES*VOLUME)^BLES (largeur du filtre utilise
!       en fonction du volume de la cellule)
!   CSMAGO
!       La constante de Smagorinsky theorique vaut 0.18
!       pour un canal plan, on prendra cependant plutot 0.065
!   XLESFD
!       Dans le cas d un modele dynamique, XLESFD est le rapport entre la
!       largeur du filtre explicite et celle du filtre implicite
!   SMAGMX
!       Constante de Smagorinsky maximale souhaitee (on peut prendre 10*CSMAGO)
!   IDRIES
!       Amortissement Van Driest active (=1) ou non (=0)
!   CDRIES
!       Constante de Van Driest dans (1-exp(-y+/CDRIES))
!   CE4    : Coefficient du terme interfacial dans k-eps
!            (Ce coefficient sert en Lagrangien)
!   VOLMIN : VOLUME DE CONTROLE MINIMAL
!   VOLMAX : VOLUME DE CONTROLE MAXIMAL
!   VOLTOT : VOLUME TOTAL DU DOMAINE

double precision  xkappa , cstlog , ypluli(nphsmx)  ,             &
                  apow   , bpow   , cpow   , dpow   ,             &
                  cmu    , cmu025 , ce1    , ce2    , ce4    ,    &
                  sigmak , sigmae ,                               &
                  crij1  , crij2  , crij3  , crijep , csrij  ,    &
                  crijp1 , crijp2 ,                               &
                  cssge2 , cssgs1 , cssgs2 ,                      &
                  cssgr1 , cssgr2 , cssgr3 , cssgr4 , cssgr5 ,    &
                  cv2fa1 , cv2fe2 , cv2fmu , cv2fc1 , cv2fc2 ,    &
                  cv2fct , cv2fcl , cv2fet ,                      &
                  ckwsk1 , ckwsk2 , ckwsw1 , ckwsw2 , ckwbt1 ,    &
                  ckwbt2 , ckwgm1 , ckwgm2 , ckwa1  , ckwc1  ,    &
                  volmin , volmax , voltot ,                      &
                  almax (nphsmx)  , uref  (nphsmx),               &
                  xlomlg(nphsmx)  ,                               &
                  xlesfl(nphsmx)  , ales  (nphsmx), bles(nphsmx), &
                  csmago(nphsmx)  , cdries(nphsmx),               &
                  xlesfd(nphsmx)  , smagmx(nphsmx)
common / rturbu / xkappa , cstlog , ypluli ,                      &
                  apow   , bpow   , cpow   , dpow   ,             &
                  cmu    , cmu025 , ce1    , ce2    , ce4    ,    &
                  sigmak , sigmae ,                               &
                  crij1  , crij2  , crij3  , crijep , csrij  ,    &
                  crijp1 , crijp2 ,                               &
                  cssge2 , cssgs1 , cssgs2 ,                      &
                  cssgr1 , cssgr2 , cssgr3 , cssgr4 , cssgr5 ,    &
                  cv2fa1 , cv2fe2 , cv2fmu , cv2fc1 , cv2fc2 ,    &
                  cv2fct , cv2fcl , cv2fet ,                      &
                  ckwsk1 , ckwsk2 , ckwsw1 , ckwsw2 , ckwbt1 ,    &
                  ckwbt2 , ckwgm1 , ckwgm2 , ckwa1  , ckwc1  ,    &
                  volmin , volmax , voltot ,                      &
                  almax           , uref            ,             &
                  xlomlg          ,                               &
                  xlesfl          , ales            , bles      , &
                  csmago          , cdries          ,             &
                  xlesfd          , smagmx

! CONSTANTES POUR LES SCALAIRES

! ISCSTH :
!   -1 : DE TYPE TEMPERATURE EN C (      CP POUR LA LOI DE PAROI)
!    0 : SCALAIRE PASSIF      (IE PAS DE CP POUR LA LOI DE PAROI)
!    1 : DE TYPE TEMPERATURE EN K (      CP POUR LA LOI DE PAROI)
!    2 : ENTHALPIE            (IE PAS DE CP POUR LA LOI DE PAROI)
!      LA DISTINCTION C/K SERT EN RAYONNEMENT
! IVISLS : SI POSITIF STRICTEMENT, INDIQUE QUE LA VISCOSITE ASSOCIEE
!            AU SCALAIRE EST VARIABLE, ET LA VALEUR EST LE NUMERO
!            D'ORDRE DE LA VISCOSITE DANS LE TABLEAU DES VISCOSITES
!            VARIABLES
! IVISSA : COMME IVISLS SAUF QUE SERT AU STOCKAGE DE LA VISCOSITE AU
!          PAS DE TEMPS PRECEDENT
! ICLVFL : 0 : CLIPPING DES VARIANCES A ZERO
!          1 : CLIPPING DES VARIANCES A ZERO ET A f(1-f)
!          2 : CLIPPING DES VARIANCES A MAX(ZERO,SCAMIN) ET SCAMAX
! ISCAVR : NUMERO DU SCALAIRE ASSOCIE A LA VARIANCE OU ZERO
!          SI LE SCALAIRE N'EST PAS UNE VARIANCE
! IPHSCA : NUMERO DE LA PHASE PORTEUSE
! SCAMIN, SCAMAX : MIN ET MAX POUR CLIPPING DES SCALAIRES
!                  ON NE CLIPPE QUE SI SCAMIN < SCAMAX
! VISLS0 : VISCOSITE DES SCALAIRES SI CONSTANTE
! SIGMAS : PRANDTL DES SCALAIRES
! RVARFL : COEFF DE DISSIPATION DES VARIANCES

integer           iscsth(nscamx),ivisls(nscamx),ivissa(nscamx),   &
                  iclvfl(nscamx),                                 &
                  iscavr(nscamx),iphsca(nscamx)
double precision  scamin(nscamx),scamax(nscamx),                  &
                  visls0(nscamx),sigmas(nscamx),                  &
                  rvarfl(nscamx)
common / iscala / iscsth        ,ivisls        ,ivissa        ,   &
                  iclvfl        ,                                 &
                  iscavr        ,iphsca
common / rscala / scamin        ,scamax        ,                  &
                  visls0        ,sigmas        ,                  &
                  rvarfl

! FIN
