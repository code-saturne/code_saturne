c@a
c@versb
C-----------------------------------------------------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2008 EDF S.A., France
C
C     contact: saturne-support@edf.fr
C
C     The Code_Saturne Kernel is free software; you can redistribute it
C     and/or modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2 of
C     the License, or (at your option) any later version.
C
C     The Code_Saturne Kernel is distributed in the hope that it will be
C     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with the Code_Saturne Kernel; if not, write to the
C     Free Software Foundation, Inc.,
C     51 Franklin St, Fifth Floor,
C     Boston, MA  02110-1301  USA
C
C-----------------------------------------------------------------------
c@verse
C                              cstphy.h
C***********************************************************************
C
C KELVIN
C
C       TKELVI       --> =  273,15
C       TKELVN       --> = -273,15
C
      DOUBLE PRECISION   TKELVI            , TKELVN
      PARAMETER (        TKELVI =  273.15D0, TKELVN = -273.15D0)
C
C CALORIES
C
C       1 cal = XCAL2J J
C
      DOUBLE PRECISION   XCAL2J
      PARAMETER (        XCAL2J = 4.1855D0)
C
C STEPHAN BOLTZMANN
C
      DOUBLE PRECISION   STEPHN
      PARAMETER (        STEPHN = 5.6703D-8)
C
C GRAVITE
C
      DOUBLE PRECISION   GX,GY,GZ
      COMMON / RGRAVI /  GX,GY,GZ
C
C CONSTANTES PHYSIQUES DU FLUIDE
C   IXYZP0 : INDICATEUR DE REMPLISSAGE DE XYZP0
C   RO0    : MASSE VOLUMIQUE    DE REFERENCE
C   VISCL0 : VISCOSITE          DE REFERENCE
C   P0     : PRESSION TOTALE    DE REFERENCE
C   PRED0  : PRESSION REDUITE   DE REFERENCE
C   XYZP0  : POSITION PRESSION  DE REFERENCE
C   T0     : TEMPERATURE        DE REFERENCE
C   CP0    : CHALEUR SPECIFIQUE DE REFERENCE
C
      INTEGER           IXYZP0(NPHSMX)
      COMMON / ICSTFL / IXYZP0
      DOUBLE PRECISION  RO0(NPHSMX)    , VISCL0(NPHSMX),
     &                  P0 (NPHSMX)    , PRED0 (NPHSMX),
     &                  XYZP0(3,NPHSMX), T0    (NPHSMX),
     &                  CP0(NPHSMX)
      COMMON / RCSTFL / RO0        , VISCL0        ,
     &                  P0         , PRED0         ,
     &                  XYZP0      , T0            ,
     &                  CP0
C
C TURBULENCE
C   IVISLS = 0 : VISCOSITE LAMINAIRE CONSTANTE = VISLS0
C   XKAPPA : CST DE KARMAN (~0.42)
C   CSTLOG : CST DE LA LOI LOG: 1/XKAPPA*LOG(YPLUS) + CSTLOG (~5.2)
C   YPLULI : YPLUS LIMITE 1./XKAPPA OU 10.88 SI IDEUCH=2
C   *POW   : COEFF WERNER AND WENGLE
C   CMU025 = CMU**0.25
C   CE1, CE2, SIGMAK, SIGMAE :
C            CONSTANTES DU K-EPSILON
C   C*RIJ* : CONSTANTES DU Rij-EPSILON STANDARD (LRR)
C   CSSG*  : CONSTANTES SPECIFIQUES DU RIJ-EPSILON SSG
C   CV2F*  : CONSTANTES SPECIFIQUES DU V2F PHI-MODEL
C   CKW*   : CONSTANTES SPECIFIQUES DU K-OMEGA SST
C            (SK=SIGMA_K, SW=SIGMA_W, BT=BETA, GM=GAMMA)
C   ALMAX  : ECHELLE DE LONGUEUR TURBULENTE
C   UREF   : VITESSE DE REFERENCE
C   XLOMLG : LONGUEUR POUR LONGUEUR DE MELANGE
C   XLESFL, ALES, BLES
C       DELTA = XLESFL * (ALES*VOLUME)^BLES (largeur du filtre utilise
C       en fonction du volume de la cellule)
C   CSMAGO
C       La constante de Smagorinsky theorique vaut 0.18
C       pour un canal plan, on prendra cependant plutot 0.065
C   XLESFD
C       Dans le cas d un modele dynamique, XLESFD est le rapport entre la
C       largeur du filtre explicite et celle du filtre implicite
C   SMAGMX
C       Constante de Smagorinsky maximale souhaitee (on peut prendre 10*CSMAGO)
C   IDRIES
C       Amortissement Van Driest active (=1) ou non (=0)
C   CDRIES
C       Constante de Van Driest dans (1-exp(-y+/CDRIES))
C   CE4    : Coefficient du terme interfacial dans k-eps
C            (Ce coefficient sert en Lagrangien)
C   VOLTOT : VOLUME TOTAL DU DOMAINE
C
      DOUBLE PRECISION  XKAPPA , CSTLOG , YPLULI(NPHSMX)  ,
     &                  APOW   , BPOW   , CPOW   , DPOW   ,
     &                  CMU    , CMU025 , CE1    , CE2    , CE4    ,
     &                  SIGMAK , SIGMAE ,
     &                  CRIJ1  , CRIJ2  , CRIJ3  , CRIJEP , CSRIJ  ,
     &                  CRIJP1 , CRIJP2 ,
     &                  CSSGE2 , CSSGS1 , CSSGS2 ,
     &                  CSSGR1 , CSSGR2 , CSSGR3 , CSSGR4 , CSSGR5 ,
     &                  CV2FA1 , CV2FE2 , CV2FMU , CV2FC1 , CV2FC2 ,
     &                  CV2FCT , CV2FCL , CV2FET ,
     &                  CKWSK1 , CKWSK2 , CKWSW1 , CKWSW2 , CKWBT1 ,
     &                  CKWBT2 , CKWGM1 , CKWGM2 , CKWA1  , CKWC1  ,
     &                  VOLTOT ,
     &                  ALMAX (NPHSMX)  , UREF  (NPHSMX),
     &                  XLOMLG(NPHSMX)  ,
     &                  XLESFL(NPHSMX)  , ALES  (NPHSMX), BLES(NPHSMX),
     &                  CSMAGO(NPHSMX)  , CDRIES(NPHSMX),
     &                  XLESFD(NPHSMX)  , SMAGMX(NPHSMX)
      COMMON / RTURBU / XKAPPA , CSTLOG , YPLULI ,
     &                  APOW   , BPOW   , CPOW   , DPOW   ,
     &                  CMU    , CMU025 , CE1    , CE2    , CE4    ,
     &                  SIGMAK , SIGMAE ,
     &                  CRIJ1  , CRIJ2  , CRIJ3  , CRIJEP , CSRIJ  ,
     &                  CRIJP1 , CRIJP2 ,
     &                  CSSGE2 , CSSGS1 , CSSGS2 ,
     &                  CSSGR1 , CSSGR2 , CSSGR3 , CSSGR4 , CSSGR5 ,
     &                  CV2FA1 , CV2FE2 , CV2FMU , CV2FC1 , CV2FC2 ,
     &                  CV2FCT , CV2FCL , CV2FET ,
     &                  CKWSK1 , CKWSK2 , CKWSW1 , CKWSW2 , CKWBT1 ,
     &                  CKWBT2 , CKWGM1 , CKWGM2 , CKWA1  , CKWC1  ,
     &                  VOLTOT ,
     &                  ALMAX           , UREF            ,
     &                  XLOMLG          ,
     &                  XLESFL          , ALES            , BLES      ,
     &                  CSMAGO          , CDRIES          ,
     &                  XLESFD          , SMAGMX
C
C CONSTANTES POUR LES SCALAIRES
C
C ISCSTH :
C   -1 : DE TYPE TEMPERATURE EN C (      CP POUR LA LOI DE PAROI)
C    0 : SCALAIRE PASSIF      (IE PAS DE CP POUR LA LOI DE PAROI)
C    1 : DE TYPE TEMPERATURE EN K (      CP POUR LA LOI DE PAROI)
C    2 : ENTHALPIE            (IE PAS DE CP POUR LA LOI DE PAROI)
C      LA DISTINCTION C/K SERT EN RAYONNEMENT
C IVISLS : SI POSITIF STRICTEMENT, INDIQUE QUE LA VISCOSITE ASSOCIEE
C            AU SCALAIRE EST VARIABLE, ET LA VALEUR EST LE NUMERO
C            D'ORDRE DE LA VISCOSITE DANS LE TABLEAU DES VISCOSITES
C            VARIABLES
C IVISSA : COMME IVISLS SAUF QUE SERT AU STOCKAGE DE LA VISCOSITE AU
C          PAS DE TEMPS PRECEDENT
C ICLVFL : 0 : CLIPPING DES VARIANCES A ZERO
C          1 : CLIPPING DES VARIANCES A ZERO ET A f(1-f)
C          2 : CLIPPING DES VARIANCES A MAX(ZERO,SCAMIN) ET SCAMAX
C ISCAVR : NUMERO DU SCALAIRE ASSOCIE A LA VARIANCE OU ZERO
C          SI LE SCALAIRE N'EST PAS UNE VARIANCE
C IPHSCA : NUMERO DE LA PHASE PORTEUSE
C SCAMIN, SCAMAX : MIN ET MAX POUR CLIPPING DES SCALAIRES
C                  ON NE CLIPPE QUE SI SCAMIN < SCAMAX
C VISLS0 : VISCOSITE DES SCALAIRES SI CONSTANTE
C SIGMAS : PRANDTL DES SCALAIRES
C RVARFL : COEFF DE DISSIPATION DES VARIANCES
C
      INTEGER           ISCSTH(NSCAMX),IVISLS(NSCAMX),IVISSA(NSCAMX),
     &                  ICLVFL(NSCAMX),
     &                  ISCAVR(NSCAMX),IPHSCA(NSCAMX)
      DOUBLE PRECISION  SCAMIN(NSCAMX),SCAMAX(NSCAMX),
     &                  VISLS0(NSCAMX),SIGMAS(NSCAMX),
     &                  RVARFL(NSCAMX)
      COMMON / ISCALA / ISCSTH        ,IVISLS        ,IVISSA        ,
     &                  ICLVFL        ,
     &                  ISCAVR        ,IPHSCA
      COMMON / RSCALA / SCAMIN        ,SCAMAX        ,
     &                  VISLS0        ,SIGMAS        ,
     &                  RVARFL
C
C FIN
c@z
