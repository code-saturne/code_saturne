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
                        SUBROUTINE IMPINI
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  IMPRESSION DES PARAMETRES DE CALCUL APRES INTERVENTION
CFONC    UTILISATEUR DANS usini1
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "cstnum.h"
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "albase.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
      INCLUDE "radiat.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagdim.h"
      INCLUDE "lagran.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
C
C VARIABLES LOCALES
C
      CHARACTER        NAME*300, CHAINE*80
      INTEGER          IOK20 , IOK21 , IOK30 , IOK31 , IOK50 , IOK60
      INTEGER          II    , JJ    , IVAR  , IPHAS , IIESCA, IEST
      INTEGER          IPP   , IWAR  , IMOM
      INTEGER          NBCCOU
C
C***********************************************************************
C
C
C=======================================================================
C 1. INTRODUCTION
C=======================================================================
C
      WRITE(NFECRA,1000)
C
      IF (IPPMOD(ICOD3P).NE.-1) THEN
        WRITE(NFECRA,1010)
        WRITE(NFECRA,1020) IPPMOD(ICOD3P)
        WRITE(NFECRA,1060) INDJON
      ELSE IF (IPPMOD(ICOEBU).NE.-1) THEN
        WRITE(NFECRA,1010)
        WRITE(NFECRA,1030) IPPMOD(ICOEBU), CEBU
        WRITE(NFECRA,1060) INDJON
      ELSE IF (IPPMOD(ICP3PL).NE.-1) THEN
        WRITE(NFECRA,1010)
        WRITE(NFECRA,1040) IPPMOD(ICP3PL)
      ELSE IF (IPPMOD(ICFUEL).NE.-1) THEN
        WRITE(NFECRA,1010)
        WRITE(NFECRA,1050) IPPMOD(ICFUEL)
      ENDIF

C
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(
     &                                                                /,
     &' ===========================================================' ,/,
     &                                                                /,
     &'               RESUME DES PARAMETRES DE CALCUL'               ,/,
     &'               ==============================='               ,/,
     &                                                                /,
     &' -----------------------------------------------------------' ,/)
C
 9900 FORMAT(
     &                                                                /,
     &' -----------------------------------------------------------' ,/)
 1010 FORMAT(
     &                                                                /,
     &' ** PHYSIQUE PARTICULIERE :'                                  ,/,
     &'    ---------------------'                                    ,/)
 1020 FORMAT(
     &' --- Flamme de diffusion : Chimie 3 points'                   ,/,
     &'       OPTION = ',4X,I10                                       /)
 1030 FORMAT(
     &' --- Flamme premelangee : Modele EBU'                         ,/,
     &'       OPTION = ',4X,I10                                      ,/,
     &'       CEBU   = ',E14.5                                        /)
 1040 FORMAT(
     &' --- Charbon pulverise : Modele Combustible moyen local'      ,/,
     &'       OPTION = ',4X,I10                                       /)
 1050 FORMAT(
     &' --- Fuel              : Modele Combustible moyen local'       ,/,
     &'       OPTION = ',4X,I10                                       /)
 1060 FORMAT(
     &' --- Janaf ou non (dans ce cas tabulation utilisateur)'       ,/,
     &'       INDJON = ',4X,I10,    ' (1: Janaf, 0: utilisateur)'    ,/)
C
#else
C
 1000 FORMAT(
     &                                                                /,
     &' ===========================================================' ,/,
     &                                                                /,
     &'               CALCULATION PARAMETERS SUMMARY'                ,/,
     &'               =============================='                ,/,
     &                                                                /,
     &' -----------------------------------------------------------' ,/)
C
 9900 FORMAT(
     &                                                                /,
     &' -----------------------------------------------------------' ,/)
 1010 FORMAT(
     &                                                                /,
     &' ** SPECIFIC PHYSICS:'                                        ,/,
     &'    ----------------'                                         ,/)
 1020 FORMAT(
     &' --- Diffusion Flame: 3 Point Chemistry'                      ,/,
     &'       OPTION = ',4X,I10                                       /)
 1030 FORMAT(
     &' --- Premixed Flame: EBU Model'                               ,/,
     &'       OPTION = ',4X,I10                                      ,/,
     &'       CEBU   = ',E14.5                                        /)
 1040 FORMAT(
     &' --- Pulverized Coal: Local Mean Combustible Model'           ,/,
     &'       OPTION = ',4X,I10                                       /)
 1050 FORMAT(
     &' --- Fuel:            Local Mean Combustible Model'           ,/,
     &'       OPTION = ',4X,I10                                       /)
 1060 FORMAT(
     &' --- Janaf or not (user tabulation required in this case)'    ,/,
     &'       INDJON = ',4X,I10,    ' (1: Janaf, 0: user)'           ,/)
C
#endif
C
C=======================================================================
C 2. DEFINITION GENERALE DU CAS
C=======================================================================
C
C --- Dimensions
C
      WRITE(NFECRA,1500)
      WRITE(NFECRA,1510) NPRFML,NFML
      WRITE(NFECRA,1520) NPHAS,NVAR,NSCAL,NSCAUS,NSCAPP,
     &                   NPROCE,NPROFA,NPROFB
C
      WRITE(NFECRA,9900)
C
C
#if defined(_CS_LANG_FR)
C
 1500 FORMAT(
     &'                                                             ',/,
     &' ** DIMENSIONS                                               ',/,
     &'    ----------                                               ',/)
 1510 FORMAT(
     &' --- Geometrie'                                               ,/,
     &'       NPRFML = ',4X,I10,    ' (Nb de proprietes de famille )',/,
     &'       NFML   = ',4X,I10,    ' (Nb de familles              )',/)
 1520 FORMAT(
     &' --- Physique'                                                ,/,
     &'       NPHAS  = ',4X,I10,    ' (Nb de phases                )',/,
     &'       NVAR   = ',4X,I10,    ' (Nb de variables             )',/,
     &'       NSCAL  = ',4X,I10,    ' (Nb de scalaires             )',/,
     &'       NSCAUS = ',4X,I10,    ' (Nb de scalaires utilisateur )',/,
     &'       NSCAPP = ',4X,I10,    ' (Nb de scalaires phys. part. )',/,
     &'       NPROCE = ',4X,I10,    ' (Nb de proprietes (cellules) )',/,
     &'       NPROFA = ',4X,I10,    ' (Nb de proprietes (faces int))',/,
     &'       NPROFB = ',4X,I10,    ' (Nb de proprietes (faces brd))',/)
C
#else
C
 1500 FORMAT(
     &                                                                /,
     &' ** DIMENSIONS'                                               ,/,
     &'    ----------'                                               ,/,
 1510 FORMAT(
     &' --- Geometry'                                                ,/,
     &'       NCEL   = ',4X,I10,    ' (Nb local active cells       )',/,
     &'       NFAC   = ',4X,I10,    ' (Nb local internal faces     )',/,
     &'       NFABOR = ',4X,I10,    ' (Nb local boundary faces     )',/,
     &'       NNOD   = ',4X,I10,    ' (Nb local vertices           )',/,
     &'       LNDFAC = ',4X,I10,    ' (Loc internal face->vtx size )',/,
     &'       LNDFBR = ',4X,I10,    ' (Loc boundary face->vtx size )',/,
     &'       NPRFML = ',4X,I10,    ' (Nb max. family properties   )',/,
     &'       NFML   = ',4X,I10,    ' (Nb families                 )',/)
 1520 FORMAT(
     &' --- Physics'                                                 ,/,
     &'       NPHAS  = ',4X,I10,    ' (Nb phases                   )',/,
     &'       NVAR   = ',4X,I10,    ' (Nb variables                )',/,
     &'       NSCAL  = ',4X,I10,    ' (Nb scalars                  )',/,
     &'       NSCAUS = ',4X,I10,    ' (Nb user scalars             )',/,
     &'       NSCAPP = ',4X,I10,    ' (Nb specific physics scalars )',/,
     &'       NPROCE = ',4X,I10,    ' (Nb cell properties          )',/,
     &'       NPROFA = ',4X,I10,    ' (Nb internal face properties )',/,
     &'       NPROFB = ',4X,I10,    ' (Nb boundary face properties )',/)
C
#endif
C
C=======================================================================
C 3. MODELISATION PHYSIQUE
C=======================================================================
C
C --- Proprietes physiques
C
      WRITE(NFECRA,2000)
      WRITE(NFECRA,2010) GX,GY,GZ
C
      DO IPHAS = 1, NPHAS
          WRITE(NFECRA,2020) IPHAS,
     &                       RO0   (IPHAS),VISCL0(IPHAS),
     &                       CP0   (IPHAS),ICP   (IPHAS),
     &                       P0    (IPHAS),PRED0 (IPHAS),
     &                       T0    (IPHAS),
     &                       IROVAR(IPHAS),IVIVAR(IPHAS),
     &                       (XYZP0(II,IPHAS),II=1,3)
      ENDDO
C
      IF (NSCAPP.GT.0) WRITE(NFECRA,2030) DIFTL0
C
C
      WRITE(NFECRA,9900)
C
C
#if defined(_CS_LANG_FR)
C
 2000 FORMAT(
     &                                                                /,
     &' ** PROPRIETES PHYSIQUES'                                     ,/,
     &'    --------------------'                                     ,/)
 2010 FORMAT(
     &'       GX     = ', E14.5,    ' (Composante x de la gravite  )',/,
     &'       GY     = ', E14.5,    ' (Composante y de la gravite  )',/,
     &'       GZ     = ', E14.5,    ' (Composante z de la gravite  )',/)
 2020 FORMAT(
     &'  -- Phase : ',I10                                            ,/,
     &                                                                /,
     &'       RO0    = ', E14.5,    ' (Masse volumique     de ref. )',/,
     &'       VISCL0 = ', E14.5,    ' (Visc. molec. dynam. de ref. )',/,
     &'       CP0    = ', E14.5,    ' (Chal. Spec.     de reference)',/,
     &'       ICP    = ',4X,I10,    ' (> 0 : CP variable   (usphyv))',/,
     &'       P0     = ', E14.5,    ' (Pression totale de reference)',/,
     &'       PRED0  = ', E14.5,    ' (Press. reduite  de reference)',/,
     &'       T0     = ', E14.5,    ' (Temperature     de reference)',/,
     &                                                                /,
     &'       IROVAR = ',4X,I10,    ' (Masse vol.  cst (0) ou non(1)',/,
     &'       IVIVAR = ',4X,I10,    ' (Visc molec. cst (0) ou non(1)',/,
     &/,
     &'       Point de reference initial pour la pression           ',/,
     &'       XYZP0  = ', E14.5, E14.5, E14.5                          )
 2030 FORMAT(
     &'       DIFTL0 = ', E14.5,    ' (Diff. dynam.    de reference)',/)
C
#else
C
 2000 FORMAT(
     &                                                                /,
     &' ** PHYSICAL PROPERTIES'                                      ,/,
     &'    -------------------'                                      ,/)
 2010 FORMAT(
     &'       GX     = ', E14.5,    ' (Gravity x component         )',/,
     &'       GY     = ', E14.5,    ' (Gravity y component         )',/,
     &'       GZ     = ', E14.5,    ' (Gravity z component         )',/)
 2020 FORMAT(
     &'  -- Phase: ',I10                                             ,/,
     &                                                                /,
     &'       RO0    = ', E14.5,    ' (Reference density           )',/,
     &'       VISCL0 = ', E14.5,    ' (Ref. molecular dyn. visc.   )',/,
     &'       CP0    = ', E14.5,    ' (Ref. specific heat          )',/,
     &'       ICP    = ',4X,I10,    ' (> 0: variable CP (usphyv)   )',/,
     &'       P0     = ', E14.5,    ' (Ref. total pressure         )',/,
     &'       PRED0  = ', E14.5,    ' (Ref. reduced pressure       )',/,
     &'       T0     = ', E14.5,    ' (Ref. temperature            )',/,
     &                                                                /,
     &'       IROVAR = ',4X,I10,    ' (Density constant(0) or not(1)',/,
     &'       IVIVAR = ',4X,I10,    ' (Molec. visc cst.(0) or not(1)',/,
     &/,
     &'       Initial reference point for pressure'                  ,/,
     &'       XYZP0  = ', E14.5, E14.5, E14.5                          )
 2030 FORMAT(
     &'       DIFTL0 = ', E14.5,    ' (Ref. dynamic diffusivity    )',/)
C
#endif
C
C --- Turbulence
C
      WRITE(NFECRA,2510)
C
C   - Modeles
C
      DO IPHAS = 1, NPHAS
        WRITE(NFECRA,2515)IPHAS,
     &    ITURB(IPHAS),IDEUCH(IPHAS),YPLULI(IPHAS),ILOGPO(IPHAS),
     &    IGRHOK(IPHAS),ISCALT(IPHAS)
        IF(ITURB(IPHAS).EQ.10) THEN
          WRITE(NFECRA,2516)
     &      XLOMLG(IPHAS)
        ELSEIF(ITURB(IPHAS).EQ.20) THEN
          WRITE(NFECRA,2517)
     &      ALMAX(IPHAS) ,UREF(IPHAS)  ,
     &      ICLKEP(IPHAS),IKECOU(IPHAS),IGRAKE(IPHAS)
          IF (IKECOU(IPHAS).EQ.0 .AND. IDTVAR.GE.0) THEN
            WRITE(NFECRA,2527) RELAXV(IK(IPHAS)),RELAXV(IEP(IPHAS))
          ELSE
            WRITE(NFECRA,2529)
          ENDIF
        ELSEIF(ITURB(IPHAS).EQ.21) THEN
          WRITE(NFECRA,2518)
     &      ALMAX(IPHAS) ,UREF(IPHAS)  ,
     &      ICLKEP(IPHAS),IKECOU(IPHAS),IGRAKE(IPHAS)
          IF (IKECOU(IPHAS).EQ.0.AND. IDTVAR.GE.0) THEN
            WRITE(NFECRA,2527) RELAXV(IK(IPHAS)),RELAXV(IEP(IPHAS))
          ELSE
            WRITE(NFECRA,2529)
          ENDIF
        ELSEIF(ITURB(IPHAS).EQ.30) THEN
          WRITE(NFECRA,2519)
     &      ALMAX(IPHAS) ,UREF(IPHAS)  ,
     &      IRIJNU(IPHAS),IRIJRB(IPHAS),IRIJEC(IPHAS),
     &      IDIFRE(IPHAS),IGRARI(IPHAS),ICLSYR(IPHAS),ICLPTR(IPHAS)
        ELSEIF(ITURB(IPHAS).EQ.31) THEN
          WRITE(NFECRA,2520)
     &      ALMAX(IPHAS) ,UREF(IPHAS)  ,
     &      IRIJNU(IPHAS),IRIJRB(IPHAS),
     &      IGRARI(IPHAS),ICLSYR(IPHAS),ICLPTR(IPHAS)
        ELSEIF(ITYTUR(IPHAS).EQ.4) THEN
          WRITE(NFECRA,2521)
     &      CSMAGO(IPHAS),XLESFL(IPHAS),ALES(IPHAS),
     &      BLES(IPHAS),IDRIES(IPHAS),CDRIES(IPHAS), XLESFD(IPHAS),
     &      SMAGMX(IPHAS), IVRTEX
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
          WRITE(NFECRA,2522)
     &      ALMAX(IPHAS) ,UREF(IPHAS)  ,
     &      ICLKEP(IPHAS),IKECOU(IPHAS),IGRAKE(IPHAS)
          IF (IKECOU(IPHAS).EQ.0 .AND. IDTVAR.GE.0) THEN
            WRITE(NFECRA,2527) RELAXV(IK(IPHAS)),RELAXV(IEP(IPHAS))
          ELSE
            WRITE(NFECRA,2529)
          ENDIF
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
          WRITE(NFECRA,2523)
     &      ALMAX(IPHAS) ,UREF(IPHAS)  ,
     &      IKECOU(IPHAS),IGRAKE(IPHAS)
          IF (IKECOU(IPHAS).EQ.0 .AND. IDTVAR.GE.0) THEN
            WRITE(NFECRA,2528) RELAXV(IK(IPHAS)),RELAXV(IOMG(IPHAS))
          ELSE
            WRITE(NFECRA,2529)
          ENDIF
        ENDIF
      ENDDO
C
C   - Constantes
C
      WRITE(NFECRA,2530)XKAPPA,CSTLOG,APOW,BPOW
C
      IOK20 = 0
      IOK21 = 0
      IOK30 = 0
      IOK31 = 0
      IOK50 = 0
      IOK60 = 0
      DO IPHAS = 1, NPHAS
        IF(ITURB(IPHAS).EQ.20) THEN
          IOK20 = 20
        ENDIF
        IF(ITURB(IPHAS).EQ.21) THEN
          IOK21 = 21
        ENDIF
        IF(ITURB(IPHAS).EQ.30) THEN
          IOK30 = 30
        ENDIF
        IF(ITURB(IPHAS).EQ.31) THEN
          IOK31 = 31
        ENDIF
        IF(ITURB(IPHAS).EQ.50) THEN
          IOK50 = 50
        ENDIF
        IF(ITURB(IPHAS).EQ.60) THEN
          IOK60 = 60
        ENDIF
      ENDDO
      IF(IOK20.GT.0) THEN
        WRITE(NFECRA,2531)CE1,CE2,SIGMAK,SIGMAE,CMU
      ENDIF
      IF(IOK21.GT.0) THEN
        WRITE(NFECRA,2532)CE1,CE2,SIGMAK,SIGMAE,CMU
      ENDIF
      IF(IOK30.GT.0) THEN
        WRITE(NFECRA,2533)CE1,CE2,CRIJ1,CRIJ2,CRIJ3,CRIJEP,CSRIJ,
     &                    CRIJP1,CRIJP2,CMU
      ENDIF
      IF(IOK31.GT.0) THEN
        WRITE(NFECRA,2534)CSSGS1,CSSGS2,CSSGR1,CSSGR2,CSSGR3,CSSGR4,
     &       CSSGR5,CSRIJ,CRIJ3,CE1,CSSGE2,SIGMAE,CMU
      ENDIF
      IF(IOK50.GT.0) THEN
        WRITE(NFECRA,2535) CV2FA1,CV2FE2,SIGMAK,SIGMAE,CV2FMU,CV2FCT,
     &       CV2FCL,CV2FET,CV2FC1,CV2FC2
      ENDIF
      IF(IOK60.GT.0) THEN
        WRITE(NFECRA,2536) CKWSK1,CKWSK2,CKWSW1,CKWSW2,CKWBT1,CKWBT2,
     &       CKWGM1,CKWGM2,CKWA1,CKWC1
      ENDIF
C
      WRITE(NFECRA,9900)
C
C
#if defined(_CS_LANG_FR)
C
 2510 FORMAT(
     &                                                                /,
     &' ** TURBULENCE                                               ',/,
     &'    ----------                                               ',/)
 2515 FORMAT(
     &' --- Phase : ',I10                                            ,/,
     &                                                                /,
     &'   - Communs                                                 ',/,
     &'       ITURB  = ',4X,I10,    ' (Modele de turbulence        )',/,
     &'       IDEUCH = ',4X,I10,    ' (0: modele a une echelle     )',/,
     &'                               (1: modele a deux echelles   )',/,
     &'                               (2: loi de paroi invariante  )',/,
     &'       YPLULI = ', E14.5,    ' (Y plus limite               )',/,
     &'       ILOGPO = ',4X,I10,    ' (0: loi puissance (interdite  ',/,
     &'                                              en k-epsilon) )',/,
     &'                               (1: loi log une echelle      )',/,
     &'       IGRHOK = ',4X,I10,    ' (1: Grad (rho k ) calcule    )',/,
     &'       ISCALT = ',4X,I10,    ' (Numero du scalaire temp     )',/)
 2516 FORMAT(
     &'   - Longueur de melange (ITURB = 10)                        ',/,
     &'       XLOMLG = ', E14.5,    ' (Longueur caracteristique    )',/)
 2517 FORMAT(
     &'   - k-epsilon           (ITURB = 20)                        ',/,
     &'       ALMAX  = ', E14.5,    ' (Longueur caracteristique    )',/,
     &'       UREF   = ', E14.5,    ' (Vitesse  caracteristique    )',/,
     &'       ICLKEP = ',4X,I10,    ' (Mode de clipping k-epsilon  )',/,
     &'       IKECOU = ',4X,I10,    ' (Mode de couplage k-epsilon  )',/,
     &'       IGRAKE = ',4X,I10,    ' (Prise en compte de gravite  )')
 2518 FORMAT(
     &'   - k-epsilon production lineaire (ITURB = 21)              ',/,
     &'       ALMAX  = ', E14.5,    ' (Longueur caracteristique    )',/,
     &'       UREF   = ', E14.5,    ' (Vitesse  caracteristique    )',/,
     &'       ICLKEP = ',4X,I10,    ' (Mode de clipping k-epsilon  )',/,
     &'       IKECOU = ',4X,I10,    ' (Mode de couplage k-epsilon  )',/,
     &'       IGRAKE = ',4X,I10,    ' (Prise en compte de gravite  )')
 2519 FORMAT(
     &'   - Rij-epsilon         (ITURB = 30)                        ',/,
     &'       ALMAX  = ', E14.5,    ' (Longueur caracteristique    )',/,
     &'       UREF   = ', E14.5,    ' (Vitesse  caracteristique    )',/,
     &'       IRIJNU = ',4X,I10,    ' (Stabilisation matricielle   )',/,
     &'       IRIJRB = ',4X,I10,    ' (Reconstruction aux bords    )',/,
     &'       IRIJEC = ',4X,I10,    ' (Termes d echo de paroi      )',/,
     &'       IDIFRE = ',4X,I10,    ' (Traitmnt du tenseur de diff.)',/,
     &'       IGRARI = ',4X,I10,    ' (Prise en compte de gravite  )',/,
     &'       ICLSYR = ',4X,I10,    ' (Implicitation en symetrie   )',/,
     &'       ICLPTR = ',4X,I10,    ' (Implicitation en paroi      )',/)
 2520 FORMAT(
     &'   - Rij-epsilon SSG     (ITURB = 31)                        ',/,
     &'       ALMAX  = ', E14.5,    ' (Longueur caracteristique    )',/,
     &'       UREF   = ', E14.5,    ' (Vitesse  caracteristique    )',/,
     &'       IRIJNU = ',4X,I10,    ' (Stabilisation matricielle   )',/,
     &'       IRIJRB = ',4X,I10,    ' (Reconstruction aux bords    )',/,
     &'       IGRARI = ',4X,I10,    ' (Prise en compte de gravite  )',/,
     &'       ICLSYR = ',4X,I10,    ' (Implicitation en symetrie   )',/,
     &'       ICLPTR = ',4X,I10,    ' (Implicitation en paroi      )',/)
 2521 FORMAT(
     &'   - LES                 (ITURB = 40, 41)                    ',/,
     &'                               (Modele de sous-maille       )',/,
     &'                               (40 Modele de Smagorinsky    )',/,
     &'                               (41 Modele dynamique         )',/,
     &'       CSMAGO = ', E14.5,    ' (Constante de Smagorinski    )',/,
     &'       XLESFL = ', E14.5,    ' (La largeur du filtre en une )',/,
     &'       ALES   = ', E14.5,    ' (cellule s''ecrit            )',/,
     &'       BLES   = ', E14.5,    ' (XLESFL*(ALES*VOLUME)**(BLES))',/,
     &'       IDRIES = ',4X,I10,    ' (=1 Amortissement Van Driest )',/,
     &'       CDRIES = ', E14.5,    ' (Constante de Van Driest     )',/,
     &'       XLESFD = ', E14.5,    ' (Rapport entre le filtre     )',/,
     &'                               (explicite et le filtre LES  )',/,
     &'                               (valeur conseillee 1.5       )',/,
     &'       SMAGMX = ', E14.5,    ' (Smagorinsky max dans le cas )',/,
     &'                               (du modele dynamique         )',/,
     &'       IVRTEX = ',4X,I10,    ' (Utilisation de la methode   )',/,
     &'                               (des vortex                  )')
 2522 FORMAT(
     &'   - v2f phi-model       (ITURB = 50)                        ',/,
     &'       ALMAX  = ', E14.5,    ' (Longueur caracteristique    )',/,
     &'       UREF   = ', E14.5,    ' (Vitesse  caracteristique    )',/,
     &'       ICLKEP = ',4X,I10,    ' (Mode de clipping k-epsilon  )',/,
     &'       IKECOU = ',4X,I10,    ' (Mode de couplage k-epsilon  )',/,
     &'       IGRAKE = ',4X,I10,    ' (Prise en compte de gravite  )')
 2523 FORMAT(
     &'   - k-omega SST         (ITURB = 60)                        ',/,
     &'       ALMAX  = ', E14.5,    ' (Longueur caracteristique    )',/,
     &'       UREF   = ', E14.5,    ' (Vitesse  caracteristique    )',/,
     &'       IKECOU = ',4X,I10,    ' (Mode de couplage k-omega    )',/,
     &'       IGRAKE = ',4X,I10,    ' (Prise en compte de gravite  )')
 2527 FORMAT(
     &'       RELAXV = ', E14.5,    ' pour k       (Relaxation)     ',/,
     &'       RELAXV = ', E14.5,    ' pour epsilon (Relaxation)     ',/)
 2528 FORMAT(
     &'       RELAXV = ', E14.5,    ' pour k     (Relaxation)       ',/,
     &'       RELAXV = ', E14.5,    ' pour omega (Relaxation)       ',/)
 2529 FORMAT(/)
C
 2530 FORMAT(
     &' --- Constantes                                              ',/,
     &                                                                /,
     &'   - Communs                                                 ',/,
     &'       XKAPPA = ', E14.5,    ' (Constante de Von Karman     )',/,
     &'       CSTLOG = ', E14.5,    ' (U+=Log(y+)/kappa +CSTLOG    )',/,
     &'       APOW   = ', E14.5,    ' (U+=APOW (y+)**BPOW (W&W law))',/,
     &'       BPOW   = ', E14.5,    ' (U+=APOW (y+)**BPOW (W&W law))',/)
 2531 FORMAT(
     &'   - k-epsilon           (ITURB = 20)'                        ,/,
     &'       CE1    = ', E14.5,    ' (Cepsilon 1 : coef de Prod.  )',/,
     &'       CE2    = ', E14.5,    ' (Cepsilon 2 : coef de Diss.  )',/,
     &'       SIGMAK = ', E14.5,    ' (Prandtl relatif a k         )',/,
     &'       SIGMAE = ', E14.5,    ' (Prandtl relatif a epsilon   )',/,
     &'       CMU    = ', E14.5,    ' (Constante Cmu               )',/)
 2532 FORMAT(
     &'   - k-epsilon production lineaire (ITURB = 21)'              ,/,
     &'       CE1    = ', E14.5,    ' (Cepsilon 1 : coef de Prod.  )',/,
     &'       CE2    = ', E14.5,    ' (Cepsilon 2 : coef de Diss.  )',/,
     &'       SIGMAK = ', E14.5,    ' (Prandtl relatif a k         )',/,
     &'       SIGMAE = ', E14.5,    ' (Prandtl relatif a epsilon   )',/,
     &'       CMU    = ', E14.5,    ' (Constante Cmu               )',/)
 2533 FORMAT(
     &'   - Rij-epsilon std     (ITURB = 30)'                        ,/,
     &'       CE1    = ', E14.5,    ' (Cepsilon 1 : coef de Prod.  )',/,
     &'       CE2    = ', E14.5,    ' (Cepsilon 2 : coef de Diss.  )',/,
     &'       CRIJ1  = ', E14.5,    ' (Coef terme lent             )',/,
     &'       CRIJ2  = ', E14.5,    ' (Coef terme rapide           )',/,
     &'       CRIJ3  = ', E14.5,    ' (Coef terme de gravite       )',/,
     &'       CRIJEP = ', E14.5,    ' (Coef diffusion epsilon      )',/,
     &'       CSRIJ  = ', E14.5,    ' (Coef diffusion Rij          )',/,
     &'       CRIJP1 = ', E14.5,    ' (Coef lent pour echo de paroi)',/,
     &'       CRIJP2 = ', E14.5,    ' (Coef rapide    echo de paroi)',/,
     &'       CMU    = ', E14.5,    ' (Constante Cmu               )',/)
 2534 FORMAT(
     &'   - Rij-epsilon SSG     (ITURB = 31)'                        ,/,
     &'       CSSGS1 = ', E14.5,    ' (Coef Cs1                    )',/,
     &'       CSSGS2 = ', E14.5,    ' (Coef Cs2                    )',/,
     &'       CSSGR1 = ', E14.5,    ' (Coef Cr1                    )',/,
     &'       CSSGR2 = ', E14.5,    ' (Coef Cr2                    )',/,
     &'       CSSGR3 = ', E14.5,    ' (Coef Cr3                    )',/,
     &'       CSSGR4 = ', E14.5,    ' (Coef Cr4                    )',/,
     &'       CSSGR5 = ', E14.5,    ' (Coef Cr5                    )',/,
     &'       CRIJS  = ', E14.5,    ' (Coef Cs diffusion de Rij    )',/,
     &'       CRIJ3  = ', E14.5,    ' (Coef terme de gravite       )',/,
     &'       CE1    = ', E14.5,    ' (Coef Ceps1                  )',/,
     &'       CSSGE2 = ', E14.5,    ' (Coef Ceps2                  )',/,
     &'       SIGMAE = ', E14.5,    ' (Coef sigma_eps              )',/,
     &'       CMU    = ', E14.5,    ' (Constante Cmu               )',/)
 2535 FORMAT(
     &'   - v2f phi-model       (ITURB = 50)'                        ,/,
     &'       CV2FA1 = ', E14.5,    ' (a1 pour calculer Cepsilon1  )',/,
     &'       CV2FE2 = ', E14.5,    ' (Cepsilon 2 : coef de Diss.  )',/,
     &'       SIGMAK = ', E14.5,    ' (Prandtl relatif a k         )',/,
     &'       SIGMAE = ', E14.5,    ' (Prandtl relatif a epsilon   )',/,
     &'       CV2FMU = ', E14.5,    ' (Constante Cmu               )',/,
     &'       CV2FCT = ', E14.5,    ' (Constante CT                )',/,
     &'       CV2FCL = ', E14.5,    ' (Constante CL                )',/,
     &'       CV2FET = ', E14.5,    ' (Constante C_eta             )',/,
     &'       CV2FC1 = ', E14.5,    ' (Constante C1                )',/,
     &'       CV2FC2 = ', E14.5,    ' (Constante C2                )',/)
 2536 FORMAT(
     &'   - k-omega SST         (ITURB = 60)'                        ,/,
     &'       CKWSK1 = ', E14.5,    ' (Constante sigma_k1          )',/,
     &'       CKWSK2 = ', E14.5,    ' (Constante sigma_k2          )',/,
     &'       CKWSW1 = ', E14.5,    ' (Constante sigma_omega1      )',/,
     &'       CKWSW2 = ', E14.5,    ' (Constante sigma_omega2      )',/,
     &'       CKWBT1 = ', E14.5,    ' (Constante beta1             )',/,
     &'       CKWBT2 = ', E14.5,    ' (Constante beta2             )',/,
     &'       CKWGM1 = ', E14.5,    ' (Constante gamma1            )',/,
     &'       CKWGM2 = ', E14.5,    ' (Constante gamma2            )',/,
     &'       CKWA1  = ', E14.5,    ' (Cste a1 pour calculer mu_t  )',/,
     &'       CKWC1  = ', E14.5,    ' (Cste c1 pour limiteur prod  )',/,
     &'       CMU    = ', E14.5,    ' (Cste Cmu (ou Beta*) pour    )',/,
     &'                                    conversion omega/epsilon)',/)
C
#else
C
 2510 FORMAT(
     &                                                                /,
     &' ** TURBULENCE'                                               ,/,
     &'    ----------'                                               ,/)
 2515 FORMAT(
     &' --- Phase: ',I10                                             ,/,
     &                                                                /,
     &'   - Commons                                                 ',/,
     &'       ITURB  = ',4X,I10,    ' (Turbulence model            )',/,
     &'       IDEUCH = ',4X,I10,    ' (0: one-scale model          )',/,
     &'                               (1: two-scale model          )',/,
     &'                               (2: invariant wall function  )',/,
     &'       YPLULI = ', E14.5,    ' (Limit Y+                    )',/,
     &'       ILOGPO = ',4X,I10,    ' (0: power law (forbidden for'  ,/,
     &'                                              k-epsilon)    )',/,
     &'                               (1: one-scale log law        )',/,
     &'       IGRHOK = ',4X,I10,    ' (1: computed Grad(rho k)     )',/,
     &'       ISCALT = ',4X,I10,    ' (Temperature salar number    )',/)
 2516 FORMAT(
     &'   - Mixing length       (ITURB = 10)'                        ,/,
     &'       XLOMLG = ', E14.5,    ' (Characteristic length       )',/)
 2517 FORMAT(
     &'   - k-epsilon           (ITURB = 20)'                        ,/,
     &'       ALMAX  = ', E14.5,    ' (Characteristic length       )',/,
     &'       UREF   = ', E14.5,    ' (Characteristic velocity     )',/,
     &'       ICLKEP = ',4X,I10,    ' (k-epsilon clipping model    )',/,
     &'       IKECOU = ',4X,I10,    ' (k-epsilon coupling mode     )',/,
     &'       IGRAKE = ',4X,I10,    ' (Account for gravity         )')
 2518 FORMAT(
     &'   - Linear production k-epsilon (ITURB = 21)'                ,/,
     &'       ALMAX  = ', E14.5,    ' (Characteristic length       )',/,
     &'       UREF   = ', E14.5,    ' (Characteristic velocity     )',/,
     &'       ICLKEP = ',4X,I10,    ' (k-epsilon clipping model    )',/,
     &'       IKECOU = ',4X,I10,    ' (k-epsilon coupling mode     )',/,
     &'       IGRAKE = ',4X,I10,    ' (Account for gravity         )')
 2519 FORMAT(
     &'   - Rij-epsilon         (ITURB = 30)'                        ,/,
     &'       ALMAX  = ', E14.5,    ' (Characteristic length       )',/,
     &'       UREF   = ', E14.5,    ' (Characteristic velocity     )',/,
     &'       IRIJNU = ',4X,I10,    ' (Matrix stabilization        )',/,
     &'       IRIJRB = ',4X,I10,    ' (Reconstruct at boundaries   )',/,
     &'       IRIJEC = ',4X,I10,    ' (Wall echo terms             )',/,
     &'       IDIFRE = ',4X,I10,    ' (Handle diffusion tensor     )',/,
     &'       IGRARI = ',4X,I10,    ' (Prise en compte de gravite  )',/,
     &'       ICLSYR = ',4X,I10,    ' (Symmetry implicitation      )',/,
     &'       ICLPTR = ',4X,I10,    ' (Wall implicitation          )',/)
 2520 FORMAT(
     &'   - SSG Rij-epsilon     (ITURB = 31)'                        ,/,
     &'       ALMAX  = ', E14.5,    ' (Characteristic length       )',/,
     &'       UREF   = ', E14.5,    ' (Characteristic velocity     )',/,
     &'       IRIJNU = ',4X,I10,    ' (Matrix stabilization        )',/,
     &'       IRIJRB = ',4X,I10,    ' (Reconstruct at boundaries   )',/,
     &'       IGRARI = ',4X,I10,    ' (Account for gravity         )',/,
     &'       ICLSYR = ',4X,I10,    ' (Symmetry implicitation      )',/,
     &'       ICLPTR = ',4X,I10,    ' (Wall implicitation          )',/)
 2521 FORMAT(
     &'   - LES                 (ITURB = 40, 41)'                    ,/,
     &'                               (Sub-grid scale model        )',/,
     &'                               (40 Smagorinsky model        )',/,
     &'                               (41 Dynamic model            )',/,
     &'       CSMAGO = ', E14.5,    ' (Smagorinsky constant        )',/,
     &'       XLESFL = ', E14.5,    ' (Filter with in a cell is    )',/,
     &'       ALES   = ', E14.5,    ' (written as                  )',/,
     &'       BLES   = ', E14.5,    ' (XLESFL*(ALES*VOLUME)**(BLES))',/,
     &'       IDRIES = ',4X,I10,    ' (=1 Van Driest damping       )',/,
     &'       CDRIES = ', E14.5,    ' (Van Driest constant         )',/,
     &'       XLESFD = ', E14.5,    ' (Ratio between the explicit  )',/,
     &'                               (filter and LES filter       )',/,
     &'                               (recommended value: 1.5      )',/,
     &'       SMAGMX = ', E14.5,    ' (Max Smagonsky in the        )',/,
     &'                               (dynamic model case          )',/,
     &'       IVRTEX = ',4X,I10,    ' (Use the vortex method       )')
 2522 FORMAT(
     &'   - v2f phi-model       (ITURB = 50)                        ',/,
     &'       ALMAX  = ', E14.5,    ' (Characteristic length       )',/,
     &'       UREF   = ', E14.5,    ' (Characteristic velocity     )',/,
     &'       ICLKEP = ',4X,I10,    ' (k-epsilon clipping model    )',/,
     &'       IKECOU = ',4X,I10,    ' (k-epsilon coupling mode     )',/,
     &'       IGRAKE = ',4X,I10,    ' (Account for gravity         )')
 2523 FORMAT(
     &'   - k-omega SST         (ITURB = 60)                        ',/,
     &'       ALMAX  = ', E14.5,    ' (Characteristic length       )',/,
     &'       UREF   = ', E14.5,    ' (Characteristic velocity     )',/,
     &'       IKECOU = ',4X,I10,    ' (k-epsilon coupling mode     )',/,
     &'       IGRAKE = ',4X,I10,    ' (Account for gravity         )')
 2527 FORMAT(
     &'       RELAXV = ', E14.5,    ' for k        (Relaxation)     ',/,
     &'       RELAXV = ', E14.5,    ' for epsilon  (Relaxation)     ',/)
 2528 FORMAT(
     &'       RELAXV = ', E14.5,    ' for k      (Relaxation)       ',/,
     &'       RELAXV = ', E14.5,    ' for omega  (Relaxation)       ',/)
 2529 FORMAT(/)
C
 2530 FORMAT(
     &' --- Constants'                                               ,/,
     &                                                                /,
     &'   - Commons                                                 ',/,
     &'       XKAPPA = ', E14.5,    ' (Von Karman constant         )',/,
     &'       CSTLOG = ', E14.5,    ' (U+=Log(y+)/kappa +CSTLOG    )',/,
     &'       APOW   = ', E14.5,    ' (U+=APOW (y+)**BPOW (W&W law))',/,
     &'       BPOW   = ', E14.5,    ' (U+=APOW (y+)**BPOW (W&W law))',/)
 2531 FORMAT(
     &'   - k-epsilon           (ITURB = 20)'                        ,/,
     &'       CE1    = ', E14.5,    ' (Cepsilon 1: production coef.)',/,
     &'       CE2    = ', E14.5,    ' (Cepsilon 2: dissipat.  coef.)',/,
     &'       SIGMAK = ', E14.5,    ' (Prandtl relative to k       )',/,
     &'       SIGMAE = ', E14.5,    ' (Prandtl relative to epsilon )',/,
     &'       CMU    = ', E14.5,    ' (Cmu constant                )',/)
 2532 FORMAT(
     &'   - Linear production k-epsilon (ITURB = 21)'                ,/,
     &'       CE1    = ', E14.5,    ' (Cepsilon 1: production coef.)',/,
     &'       CE2    = ', E14.5,    ' (Cepsilon 2: dissipat.  coef.)',/,
     &'       SIGMAK = ', E14.5,    ' (Prandtl relative to k       )',/,
     &'       SIGMAE = ', E14.5,    ' (Prandtl relative to epsilon )',/,
     &'       CMU    = ', E14.5,    ' (Cmu constant                )',/)
 2533 FORMAT(
     &'   - Rij-epsilon         (ITURB = 30)'                        ,/,
     &'       CE1    = ', E14.5,    ' (Cepsilon 1: production coef.)',/,
     &'       CE2    = ', E14.5,    ' (Cepsilon 2: dissipat.  coef.)',/,
     &'       CRIJ1  = ', E14.5,    ' (Slow term coefficient       )',/,
     &'       CRIJ2  = ', E14.5,    ' (Fast term coefficient       )',/,
     &'       CRIJ3  = ', E14.5,    ' (Gravity term coefficient    )',/,
     &'       CRIJEP = ', E14.5,    ' (Epsilon diffusion coeff.    )',/,
     &'       CSRIJ  = ', E14.5,    ' (Rij diffusion coeff.        )',/,
     &'       CRIJP1 = ', E14.5,    ' (Slow coeff. for wall echo   )',/,
     &'       CRIJP2 = ', E14.5,    ' (Fast coeff. for wall echo   )',/,
     &'       CMU    = ', E14.5,    ' (Cmu constant                )',/)
 2534 FORMAT(
     &'   - SSG Rij-epsilon     (ITURB = 31)'                        ,/,
     &'       CSSGS1 = ', E14.5,    ' (Cs1 coeff.                  )',/,
     &'       CSSGS2 = ', E14.5,    ' (Cs2 coeff.                  )',/,
     &'       CSSGR1 = ', E14.5,    ' (Cr1 coeff.                  )',/,
     &'       CSSGR2 = ', E14.5,    ' (Cr2 coeff.                  )',/,
     &'       CSSGR3 = ', E14.5,    ' (Cr3 coeff.                  )',/,
     &'       CSSGR4 = ', E14.5,    ' (Cr4 coeff.                  )',/,
     &'       CSSGR5 = ', E14.5,    ' (Cr5 coeff.                  )',/,
     &'       CRIJS  = ', E14.5,    ' (Rij Cs diffusion coeff.     )',/,
     &'       CRIJ3  = ', E14.5,    ' (Gravity term coeff.         )',/,
     &'       CE1    = ', E14.5,    ' (Ceps1 coeff.                )',/,
     &'       CSSGE2 = ', E14.5,    ' (Ceps2 coeff.                )',/,
     &'       SIGMAE = ', E14.5,    ' (sigma_eps coeff.            )',/,
     &'       CMU    = ', E14.5,    ' (Cmu constant                )',/)
 2535 FORMAT(
     &'   - v2f phi-model       (ITURB = 50)                        ',/,
     &'       CV2FA1 = ', E14.5,    ' (a1 to calculate Cepsilon1   )',/,
     &'       CV2FE2 = ', E14.5,    ' (Cepsilon 2: dissip. coeff.  )'/,
     &'       SIGMAK = ', E14.5,    ' (Prandtl relative to k       )',/,
     &'       SIGMAE = ', E14.5,    ' (Prandtl relative to epsilon )',/,
     &'       CV2FMU = ', E14.5,    ' (Cmu constant                )',/,
     &'       CV2FCT = ', E14.5,    ' (CT constant                 )',/,
     &'       CV2FCL = ', E14.5,    ' (CL constant                 )',/,
     &'       CV2FET = ', E14.5,    ' (C_eta constant              )',/,
     &'       CV2FC1 = ', E14.5,    ' (C1 constant                 )',/,
     &'       CV2FC2 = ', E14.5,    ' (C2 constant                 )',/)
 2536 FORMAT(
     &'   - k-omega SST         (ITURB = 60)                        ',/,
     &'       CKWSK1 = ', E14.5,    ' (sigma_k1 constant           )',/,
     &'       CKWSK2 = ', E14.5,    ' (sigma_k2 constant           )',/,
     &'       CKWSW1 = ', E14.5,    ' (sigma_omega1 constant       )',/,
     &'       CKWSW2 = ', E14.5,    ' (sigma_omega2 constant       )',/,
     &'       CKWBT1 = ', E14.5,    ' (beta1 constant              )',/,
     &'       CKWBT2 = ', E14.5,    ' (beta2 constant              )',/,
     &'       CKWGM1 = ', E14.5,    ' (gamma1 constant             )',/,
     &'       CKWGM2 = ', E14.5,    ' (gamma2 constant             )',/,
     &'       CKWA1  = ', E14.5,    ' (a1 constant to compute mu_t )',/,
     &'       CKWC1  = ', E14.5,    ' (c1 const. for prod. limiter )',/,
     &'       CMU    = ', E14.5,    ' (Cmu (or Beta*) constant for )',/,
     &'                                    omega/epsilon conversion)',/)
C
#endif
C
C --- Viscosite secondaire
C
      WRITE(NFECRA,2610)
      DO IPHAS = 1, NPHAS
        WRITE(NFECRA,2620) IPHAS, IVISSE(IPHAS)
      ENDDO
C
      WRITE(NFECRA,9900)
C
C
#if defined(_CS_LANG_FR)
C
 2610 FORMAT(
     &                                                                /,
     &' ** VISCOSITE SECONDAIRE                                     ',/,
     &'    --------------------                                     ',/)
 2620 FORMAT(
     &' --- Phase : ',I10                                            ,/,
     &'       IVISSE = ',4X,I10,    ' (1 : pris en compte          )',/)
C
#else
C
 2610 FORMAT(
     &                                                                /,
     &' ** SECONDARY VISCOSITY'                                      ,/,
     &'    -------------------'                                      ,/)
 2620 FORMAT(
     &' --- Phase: ',I10                                             ,/,
     &'       IVISSE = ',4X,I10,    ' (1: accounted for            )',/)
C
#endif
C
C --- Rayonnement thermique
C
      IF (IIRAYO.EQ.1) THEN
        WRITE(NFECRA,2630)
        DO IPHAS = 1, NPHAS
          WRITE(NFECRA,2640) IPHAS, IRAYON(IPHAS), ISCALT(IPHAS),
     &                       ISCSTH(ISCALT(IPHAS))
        ENDDO
        WRITE(NFECRA,2650) ISUIRD, NFREQR, NDIREC,
     &                     IDIVER, IMODAK, IIMPAR, IIMLUM
        WRITE(NFECRA,2660)
        DO IPHAS = 1, NPHAS
          DO II = 1,NBRAYP
            IF(IRAYVP(II,IPHAS).EQ.1) WRITE(NFECRA,2661)NBRVAP(II,IPHAS)
          ENDDO
          DO II = 1,NBRAYF
            IF(IRAYVF(II,IPHAS).EQ.1) WRITE(NFECRA,2662)NBRVAF(II,IPHAS)
          ENDDO
        ENDDO
C
        WRITE(NFECRA,9900)
C
      ENDIF
C
#if defined(_CS_LANG_FR)
C
 2630 FORMAT(
     &                                                                /,
     &' ** TRANSFERTS THERMIQUES RADIATIFS                          ',/,
     &'    -------------------------------                          ',/)
 2640 FORMAT(
     &' --- Phase : ',I10                                            ,/,
     &'       IRAYON = ',4X,I10,    ' (0 : non ; 1 : DOM ; 2 : P-1 )',/,
     &'       ICSALT = ',4X,I10,    ' (Num du sca thermique associe)',/,
     &'       ISCSTH = ',4X,I10,    ' (-1 : T(C) ; 1 : T(K) ; 2 : H)',/)
 2650 FORMAT(
     &' --- Options :                                               ',/,
     &'       ISUIRD = ',4X,I10,    ' (0 : pas de suite ; 1 : suite)',/,
     &'       NFREQR = ',4X,I10,    ' (Frequence pass. rayonnement )',/,
     &'       NDIREC = ',4X,I10,    ' (32 ou 128 directions(si DOM))',/,
     &'       IDIVER = ',4X,I10,    ' (0 1 ou 2: calcul TS radiatif)',/,
     &'       IMODAK = ',4X,I10,    ' (1: modak coef absor; 0 sinon)',/,
     &'       IIMPAR = ',4X,I10,    ' (0 1 ou 2: impr Tempera paroi)',/,
     &'       IIMLUM = ',4X,I10,    ' (0 1 ou 2: impr infos solveur)',/)
 2660 FORMAT(
     &' --- Sorties graphiques :                                    '  )
 2661 FORMAT(
     &'       NBRVAP = ',4X,A40                                        )
 2662 FORMAT(
     &'       NBRVAF = ',4X,A40                                        )
C
#else
C
 2630 FORMAT(
     &                                                                /,
     &' ** RADIATIVE THERMAL TRANSFER'                               ,/,
     &'    --------------------------'                               ,/)
 2640 FORMAT(
     &' --- Phase: ',I10                                             ,/,
     &'       IRAYON = ',4X,I10,    ' (0: no; 1: DOM; 2: P-1       )',/,
     &'       ICSALT = ',4X,I10,    ' (Assoc. thermal scalar num.  )',/,
     &'       ISCSTH = ',4X,I10,    ' (-1: T(C); 1: T(K); 2: H     )',/)
 2650 FORMAT(
     &' --- Options:'                                                ,/,
     &'       ISUIRD = ',4X,I10,    ' (0: no restart; 1: restart   )',/,
     &'       NFREQR = ',4X,I10,    ' (Radiation pass frequency    )',/,
     &'       NDIREC = ',4X,I10,    ' (32 or 128 directions(if DOM))',/,
     &'       IDIVER = ',4X,I10,    ' (0 1 or 2: compute radiat. ST)',/,
     &'       IMODAK = ',4X,I10,    ' (1: modak absor coef; 0 else )',/,
     &'       IIMPAR = ',4X,I10,    ' (0 1 or 2: print wall temp.  )',/,
     &'       IIMLUM = ',4X,I10,    ' (0 1 or 2: print solver info )',/)
 2660 FORMAT(
     &' --- Graphical output:'                                         )
 2661 FORMAT(
     &'       NBRVAP = ',4X,A40                                        )
 2662 FORMAT(
     &'       NBRVAF = ',4X,A40                                        )
C
#endif
C
C --- Compressible
C
      IF (IPPMOD(ICOMPF).GE.0) THEN
        WRITE(NFECRA,2700)
        DO IPHAS = 1, NPHAS
          WRITE(NFECRA,2710) IPHAS,
     &         ICV(IPHAS),
     &         IVISCV(IPHAS),VISCV0(IPHAS),
     &         ICFGRP(IPHAS)
        ENDDO
C
        WRITE(NFECRA,9900)
C
      ENDIF
C
#if defined(_CS_LANG_FR)
C
 2700 FORMAT(
     &                                                                /,
     &' ** COMPRESSIBLE : donnees complementaires                   ',/,
     &'    ------------                                             ',/)
 2710 FORMAT(
     &' --- Phase : ',I10                                            ,/,
     &'       ICV    = ',4X,I10,    ' (0 : Cv cst ; 1 : variable   )',/,
     &'       IVISCV = ',4X,I10,    ' (0 : kappa cst ; 1 : variable ',/,
     &'                                kappa : viscosite en volume  ',/,
     &'                                en kg/(m s)                 )',/,
     &'       VISCV0 = ',E14.5 ,    ' (Valeur de kappa si cst      )',/,
     &'       ICFGRP = ',4X,I10,    ' (1 : C.L. pression avec effet ',/,
     &'                                hydrostatique dominant      )',/)
C
#else
C
 2700 FORMAT(
     &                                                                /,
     &' ** COMPRESSIBLE: additional data'                            ,/,
     &'    ------------                                             ',/)
 2710 FORMAT(
     &' --- Phase : ',I10                                            ,/,
     &'       ICV    = ',4X,I10,    ' (0: Cv cst; 1: variable      )',/,
     &'       IVISCV = ',4X,I10,    ' (0: kappa cst; 1: variable'    ,/,
     &'                                kappa: volume viscosity'      ,/,
     &'                                in kg/(m.s)                 )',/,
     &'       VISCV0 = ',E14.5 ,    ' (kappa value if constant     )',/,
     &'       ICFGRP = ',4X,I10,    ' (1: pressure BC with dominant' ,/,
     &'                                hydrostatic effect          )',/)
C
#endif
C
C=======================================================================
C 4. DISCRETISATION DES EQUATIONS
C=======================================================================
C
C --- Marche en temps
C
      WRITE(NFECRA,3000)
C
C     Stationnaire
      IF (IDTVAR.LT.0) THEN

C   - Parametres du pas de temps
C
        WRITE(NFECRA,3010) IDTVAR, RELXST
C
C   - Champ de vitesse fige
C
        WRITE(NFECRA,3030) ICCVFG
C
C   - Coefficient de relaxation
C
        WRITE(NFECRA,3011)
        DO IPP = 2, NVPPMX
          II = ITRSVR(IPP)
          IF(II.GE.1) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,3012) CHAINE(1:8),RELAXV(II)
          ENDIF
        ENDDO
        WRITE(NFECRA,3013)
C
C     Instationnaire
      ELSE
C
C   - Parametres du pas de temps
C
        WRITE(NFECRA,3020) IDTVAR,IPTLRO,COUMAX,FOUMAX,
     &       VARRDT,DTMIN,DTMAX,DTREF
C
C   - Champ de vitesse fige
C
        WRITE(NFECRA,3030) ICCVFG
C
C   - Coef multiplicatif du pas de temps
C
        WRITE(NFECRA,3040)
        DO IPP = 2, NVPPMX
          II = ITRSVR(IPP)
          IF(II.GE.1) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,3041) CHAINE(1:8),ISTAT(II),CDTVAR(II)
          ENDIF
        ENDDO
        WRITE(NFECRA,3042)
C
C
C   - Coefficient de relaxation de la masse volumique
C
        IF (NSCAPP.GT.0) WRITE(NFECRA,3050) SRROM
C
C   - Ordre du schema en temps
C
        WRITE(NFECRA,3060)
        DO IPHAS = 1, NPHAS
          WRITE(NFECRA,3061) IPHAS, ISCHTP(IPHAS)
        ENDDO
        WRITE(NFECRA,3062)
C
      ENDIF
C
      WRITE(NFECRA,9900)
C
#if defined(_CS_LANG_FR)
C
 3000 FORMAT(
     &                                                                /,
     &' ** MARCHE EN TEMPS'                                          ,/,
     &'    ---------------'                                          ,/)
 3010 FORMAT(
     &'    ALGORITHME STATIONNAIRE'                                  ,/,
     &                                                                /,
     &' --- Parametres globaux'                                      ,/,
     &                                                                /,
     &'       IDTVAR = ',4X,I10,    ' (-1: algorithme stationnaire )',/,
     &'       RELXST = ', E14.5,    ' (Coef relaxation de reference)',/,
     &'                                                             ',/)
 3011 FORMAT(
     &' --- Coefficient de relaxation par variable'                  ,/,
     &                                                                /,
     &'---------------------'                                        ,/,
     &' Variable  RELAXV'                                            ,/,
     &'---------------------'                                          )
 3012 FORMAT(
     & 1X,    A8,      E12.4                                           )
 3013 FORMAT(
     &'----------------------------'                                 ,/,
     &                                                                /,
     &'       RELAXV =  [0.,1.]       (coeff de relaxation         )',/)
 3020 FORMAT(
     &'    ALGORITHME INSTATIONNAIRE'                                ,/,
     &                                                                /,
     &' --- Parametres du pas de temps'                              ,/,
     &                                                                /,
     &'       IDTVAR = ',4X,I10,    ' (0 cst;1,2 var(tps,tps-espace)',/,
     &'       IPTLRO = ',4X,I10,    ' (1 : clipping de DT lie a rho)',/,
     &'       COUMAX = ', E14.5,    ' (Courant maximum cible       )',/,
     &'       FOUMAX = ', E14.5,    ' (Fourier maximum cible       )',/,
     &'       VARRDT = ', E14.5,    ' (En DT var, accroissement max)',/,
     &'       DTMIN  = ', E14.5,    ' (Pas de temps min            )',/,
     &'       DTMAX  = ', E14.5,    ' (Pas de temps max            )',/,
     &'       DTREF  = ', E14.5,    ' (Pas de temps de reference   )',/,
     &                                                                /,
     &'       En pas de temps non constant (IDTVAR = 1 ou 2),'       ,/,
     &'         lorsque la valeur de COUMAX ou FOUMAX est negative'  ,/,
     &'         ou nulle, la limitation du pas de temps associee (au',/,
     &'         nombre de Courant et de Fourier, respectivement)    ',/,
     &'         n ''est pas prise en compte.'                        ,/)
 3030 FORMAT(
     &' --- Champ de vitesse fige'                                   ,/,
     &                                                                /,
     &'       ICCVFG = ',4X,I10,    ' (1 : champ de vitesse fige   )',/)
 3040 FORMAT(
     &' --- Proprietes par variable'                                 ,/,
     &                                                                /,
     &'----------------------------'                                 ,/,
     &' Variable  ISTAT      CDTVAR'                                 ,/,
     &'----------------------------'                                   )
 3041 FORMAT(
     & 1X,    A8,    I7,      E12.4                                    )
 3042 FORMAT(
     &'----------------------------'                                 ,/,
     &                                                                /,
     &'       ISTAT  =  0 ou  1       (1 pour instationnaire       )',/,
     &'       CDTVAR >  0             (coeff mult. du pas de temps )',/)
C
 3050 FORMAT(
     &'--- Coefficient de relaxation'                                ,/,
     &'    RHO(n+1)=SRROM*RHO(n)+(1-SRROM)*RHO(n+1)'                 ,/,
     &'       SRROM  = ',E14.5                                       ,/)
C
 3060 FORMAT(
     &' --- Ordre du schema en temps de base'                          )
 3061 FORMAT(
     &'     Phase : ',I10                                            ,/,
     &'       ISCHTP = ',4X,I10,    ' (1 : ordre 1 ; 2 : ordre 2   )'  )
 3062 FORMAT(
     &'                                                             '  )
C
#else
C
 3000 FORMAT(
     &                                                                /,
     &' ** TIME STEPPING'                                            ,/,
     &'    -------------'                                            ,/)
 3010 FORMAT(
     &'    STEADY ALGORITHM'                                         ,/,
     &                                                                /,
     &' --- Global parameters'                                       ,/,
     &                                                                /,
     &'       IDTVAR = ',4X,I10,    ' (-1: steady algorithm        )',/,
     &'       RELXST = ', E14.5,    ' (Reference relaxation coeff. )',/,
     &'                                                             ',/)
 3011 FORMAT(
     &' --- Per variable relaxation coefficient'                     ,/,
     &                                                                /,
     &'---------------------'                                        ,/,
     &' RELAXV variable'                                             ,/,
     &'---------------------'                                          )
 3012 FORMAT(
     & 1X,    A8,      E12.4                                           )
 3013 FORMAT(
     &'----------------------------'                                 ,/,
     &                                                                /,
     &'       RELAXV =  [0.,1.]       (relaxation coefficient      )',/)
 3020 FORMAT(
     &'    UNSTEADY ALGORITHM'                                       ,/,
     &                                                                /,
     &' --- Time step parameters'                                    ,/,
     &                                                                /,
     &'       IDTVAR = ',4X,I10,    ' (0 cst; 1,2 var (t, t-space  )',/,
     &'       IPTLRO = ',4X,I10,    ' (1: rho-related DT clipping  )',/,
     &'       COUMAX = ', E14.5,    ' (Maximum target CFL          )',/,
     &'       FOUMAX = ', E14.5,    ' (Maximum target Fourier      )',/,
     &'       VARRDT = ', E14.5,    ' (For var. DT, max. increase  )',/,
     &'       DTMIN  = ', E14.5,    ' (Minimum time step           )',/,
     &'       DTMAX  = ', E14.5,    ' (Maximum time step           )',/,
     &'       DTREF  = ', E14.5,    ' (Reference time step         )',/,
     &                                                                /,
     &'       With a non-constant time step (IDTVAR = 1 or 2),'      ,/,
     &'         when the value of COUMAX or FOUMAX is negative'      ,/,
     &'         or zero, the associated time step limitation (for'   ,/,
     &'         CFL and Fourier respectively) is ignored.'           ,/)
 3030 FORMAT(
     &' --- Frozen velocity field'                                   ,/,
     &                                                                /,
     &'       ICCVFG = ',4X,I10,    ' (1: frozen velocity field    )',/)
 3040 FORMAT(
     &' --- Per-variable properties'                                 ,/,
     &                                                                /,
     &'----------------------------'                                 ,/,
     &' Variable  ISTAT      CDTVAR'                                 ,/,
     &'----------------------------'                                   )
 3041 FORMAT(
     & 1X,    A8,    I7,      E12.4                                    )
 3042 FORMAT(
     &'----------------------------'                                 ,/,
     &                                                                /,
     &'       ISTAT  =  0 ou  1       (1 for unsteady              )',/,
     &'       CDTVAR >  0             (time step multiplier        )',/)
C
 3050 FORMAT(
     &'--- Relaxation coefficient'                                   ,/,
     &'    RHO(n+1)=SRROM*RHO(n)+(1-SRROM)*RHO(n+1)'                 ,/,
     &'       SRROM  = ',E14.5                                       ,/)
C
 3060 FORMAT(
     &' --- Order of base time stepping scheme'                        )
 3061 FORMAT(
     &'     Phase: ',I10                                             ,/,
     &'       ISCHTP = ',4X,I10,    ' (1: order 1; 2: order 2      )'  )
 3062 FORMAT(
     &'                                                             '  )
C
#endif
C
C --- Convection Diffusion
C
      WRITE(NFECRA,4000)
C
      WRITE(NFECRA,4010)
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,4020) CHAINE(1:8),
     &                       ICONV(II),IDIFF(II),IDIFFT(II),
     &                       ISCHCV(II),ISSTPC(II),
     &                       BLENCV(II),THETAV(II)
        ENDIF
      ENDDO
      WRITE(NFECRA,4030)
C
      WRITE(NFECRA,9900)
C
C --- Stokes
C
      WRITE(NFECRA,4110) IPHYDR,ICALHY,IPRCO,IPUCOU,NTERUP
      DO IPHAS = 1, NPHAS
        WRITE(NFECRA,4111)IPHAS,
     &                    IREVMC(IPHAS)
        IF (IDTVAR.GE.0) THEN
          WRITE(NFECRA,4112) RELAXV(IPR(IPHAS)),ARAK(IPHAS)
        ELSE
          WRITE(NFECRA,4113) ARAK(IPHAS)*RELAXV(IU(IPHAS))
        ENDIF
        WRITE(NFECRA,4114)ISTMPF(IPHAS),THETFL(IPHAS),
     &                    IROEXT(IPHAS),THETRO(IPHAS),
     &                    IVIEXT(IPHAS),THETVI(IPHAS),
     &                    ICPEXT(IPHAS),THETCP(IPHAS),
     &                    THETSN(IPHAS),THETST(IPHAS),EPSUP(IPHAS)
      ENDDO
C
      WRITE(NFECRA,9900)
C
C
#if defined(_CS_LANG_FR)
C
 4000 FORMAT(
     &                                                                /,
     &' ** CONVECTION - DIFFUSION'                                   ,/,
     &'    ----------------------'                                   ,/)
 4010 FORMAT(
     &'-------------------------------------------------------------',/,
     &' Variable  ICONV  IDIFF IDIFFT ISCHCV ISSTPC   BLENCV  THETAV',/,
     &'-------------------------------------------------------------'  )
 4020 FORMAT(
     & 1X,    A8,    I7,    I7,    I7,    I7,    I7,    E9.2,    E9.2  )
 4030 FORMAT(
     &'-------------------------------------------------------------',/,
     &                                                                /,
     &'       ICONV  =  0 ou  1       (1 pour convection branchee  )',/,
     &'       IDIFF  =  0 ou  1       (1 pour diff. tot branchee   )',/,
     &'       IDIFFT =  0 ou  1       (1 pour diff. turb. branchee )',/,
     &'       ISCHCV =  0 ou  1       (SOLU ou CD                  )',/,
     &'       ISSTPC =  0 ou  1       (1 : sans test de pente      )',/,
     &'       BLENCV =  [0.;1.]       (1-proportion d upwind       )',/,
     &'       THETAV =  [0.;1.]       (0.5 Crank-Nicolson/AB       )',/,
     &'                               (theta pour les termes de    )',/,
     &'                               (convection diffusion utilise)',/,
     &'                               ((1-theta)ancien+theta nouveau',/)
C
 4110 FORMAT(
     &                                                                /,
     &' ** STOKES'                                                   ,/,
     &'    ------'                                                   ,/,
     &'       IPHYDR = ',4X,I10,  ' (1 : prise en compte explicite  ',/,
     &'                ',14X,     '      de l''equilibre entre grad ',/,
     &'                ',14X,     '      de pression et termes      ',/,
     &'                ',14X,     '      sources de gravite et de   ',/,
     &'                ',14X,     '      pertes de charge          )',/,
     &'       ICALHY = ',4X,I10,  ' (1 : calcul de la pression      ',/,
     &'                ',14X,     '      hydrostatique pour les     ',/,
     &'                ',14X,     '      conditions de Dirichlet en ',/,
     &'                ',14X,     '      sortie sur la pression    )',/,
     &'       IPRCO  = ',4X,I10,  ' (1 : avec pression-continuite  )',/,
     &'       IPUCOU = ',4X,I10,  ' (1 : avec couplage U-P renforce)',/,
     &'       NTERUP = ',4X,I10,  ' (n : avec n sweep sur navsto    ',/,
     &'                ',14X,     '      pour couplage vites/pressio',/)
 4111 FORMAT(
     &'  -- Phase : ',I10                                            ,/,
     &                                                                /,
     &'       IREVMC = ',4X,I10,    ' (Mode de reconstruction vites)',/)
 4112 FORMAT(
     &'       RELAXV = ', E14.5,    ' pour la pression (relaxation) ',/,
     &'       ARAK   = ', E14.5,    ' (Facteur d Arakawa           )',/)
 4113 FORMAT(
     &'       ARAK   = ', E14.5,    ' (Facteur d Arakawa           )',/)
 4114 FORMAT(
     &'       ISTMPF = ',4X,I10,    ' (schema en temps pour le flux ',/,
     &'                ',14X,       ' (0 : explicite (THETFL = 0   )',/,
     &'                ',14X,       ' (1 : schema std (Saturne 1.0 )',/,
     &'                ',14X,       ' (2 : ordre 2   (THETFL = 0.5 )',/,
     &'       THETFL = ', E14.5,    ' (theta pour flux de masse    )',/,
     &'       IROEXT = ',4X,I10,    ' (extrap. masse volumique      ',/,
     &'                ',14X,       ' (0 : explicite                ',/,
     &'                ',14X,       ' (1 : n+thetro avec thetro=1/2 ',/,
     &'                ',14X,       ' (2 : n+thetro avec thetro=1   ',/,
     &'       THETRO = ', E14.5,    ' (theta pour masse volumique   ',/,
     &'                               ((1+theta)nouveau-theta ancien',/,
     &'       IVIEXT = ',4X,I10,    ' (extrap. viscosite totale     ',/,
     &'                ',14X,       ' (0 : explicite                ',/,
     &'                ',14X,       ' (1 : n+thetvi avec thetro=1/2 ',/,
     &'                ',14X,       ' (2 : n+thetvi avec thetro=1   ',/,
     &'       THETVI = ', E14.5,    ' (theta pour viscosite totale  ',/,
     &'                               ((1+theta)nouveau-theta ancien',/,
     &'       ICPEXT = ',4X,I10,    ' (extrap. chaleur specifique   ',/,
     &'                ',14X,       ' (0 : explicite                ',/,
     &'                ',14X,       ' (1 : n+thetcp avec thetro=1/2 ',/,
     &'                ',14X,       ' (2 : n+thetcp avec thetro=1   ',/,
     &'       THETCP = ', E14.5,    ' (theta schema chaleur spec    ',/,
     &'                               ((1+theta)nouveau-theta ancien',/,
     &'       THETSN = ', E14.5,    ' (theta schema T.S. Nav-Stokes)',/,
     &'                               ((1+theta)nouveau-theta ancien',/,
     &'       THETST = ', E14.5,    ' (theta schema T.S. Turbulence)',/,
     &'                               ((1+theta)nouveau-theta ancien',/,
     &'       EPSUP  = ', E14.5,    ' (Test d''arret du couplage    ',/,
     &'                ',14X,       '  vitesse/pression            )',/)

C
#else
C
 4000 FORMAT(
     &                                                                /,
     &' ** CONVECTION - DIFFUSION'                                   ,/,
     &'    ----------------------'                                   ,/)
 4010 FORMAT(
     &'-------------------------------------------------------------',/,
     &' Variable  ICONV  IDIFF IDIFFT ISCHCV ISSTPC   BLENCV  THETAV',/,
     &'-------------------------------------------------------------'  )
 4020 FORMAT(
     & 1X,    A8,    I7,    I7,    I7,    I7,    I7,    E9.2,    E9.2  )
 4030 FORMAT(
     &'-------------------------------------------------------------',/,
     &                                                                /,
     &'       ICONV  =  0 ou  1       (1 for convection active     )',/,
     &'       IDIFF  =  0 ou  1       (1 for total diffusion active)',/,
     &'       IDIFFT =  0 ou  1       (1 for turbulent diff. active)',/,
     &'       ISCHCV =  0 ou  1       (SOLU or CD                  )',/,
     &'       ISSTPC =  0 ou  1       (1: no slope test            )',/,
     &'       BLENCV =  [0.;1.]       (1-upwind proportion         )',/,
     &'       THETAV =  [0.;1.]       (0.5 Crank-Nicolson/AB       )',/,
     &'                               (theta for convection-       )',/,
     &'                               (diffusion terms uses        )',/,
     &'                               ((1-theta).old+theta.new     )',/)
C
 4110 FORMAT(
     &                                                                /,
     &' ** STOKES'                                                   ,/,
     &'    ------'                                                   ,/,
     &'       IPHYDR = ',4X,I10,  ' (1: account for explicit'        ,/,
     &'                ',14X,     '     balance between pressure'    ,/,
     &'                ',14X,     '     gradient, gravity source'    ,/,
     &'                ',14X,     '     terms, and head losses     )',/,
     &'       ICALHY = ',4X,I10,  ' (1: compute hydrastatic'        ,/,
     &'                ',14X,     '     pressure for Dirichlet'      ,/,
     &'                ',14X,     '     conditions for pressure'     ,/,
     &'                ',14X,     '     on outlet                  )',/,
     &'       IPRCO  = ',4X,I10,  ' (1: pressure-continuity        )',/,
     &'       IPUCOU = ',4X,I10,  ' (1: reinforced U-P coupling    )',/,
     &'       NTERUP = ',4X,I10,  ' (n: n sweeps on navsto for'      ,/,
     &'                ',14X,     '     velocity/pressure coupling )',/)
 4111 FORMAT(
     &'  -- Phase: ',I10                                             ,/,
     &                                                                /,
     &'       IREVMC = ',4X,I10,    ' (Velocity reconstruction mode)',/)
 4112 FORMAT(
     &'       RELAXV = ', E14.5,    ' for pressure (relaxation)'     ,/,
     &'       ARAK   = ', E14.5,    ' (Arakawa factor              )',/)
 4113 FORMAT(
     &'       ARAK   = ', E14.5,    ' (Arakawa factor              )',/)
 4114 FORMAT(
     &'       ISTMPF = ',4X,I10,    ' (time scheme for flow'         ,/,
     &'                ',14X,       ' (0: explicit (THETFL = 0     )',/,
     &'                ',14X,       ' (1: std scheme (Saturne 1.0  )',/,
     &'                ',14X,       ' (2: 2nd-order (THETFL = 0.5  )',/,
     &'       THETFL = ', E14.5,    ' (theta for mass flow         )',/,
     &'       IROEXT = ',4X,I10,    ' (density extrapolation'        ,/,
     &'                ',14X,       ' (0: explicit'                  ,/,
     &'                ',14X,       ' (1: n+thetro with thetro=1/2'  ,/,
     &'                ',14X,       ' (2: n+thetro with thetro=1'    ,/,
     &'       THETRO = ', E14.5,    ' (theta for density'            ,/,
     &'                               ((1+theta).new-theta.old'      ,/,
     &'       IVIEXT = ',4X,I10,    ' (total viscosity extrapolation',/,
     &'                ',14X,       ' (0: explicit'                  ,/,
     &'                ',14X,       ' (1: n+thetvi with thetro=1/2'  ,/,
     &'                ',14X,       ' (2: n+thetvi with thetro=1'    ,/,
     &'       THETVI = ', E14.5,    ' (theta for total viscosity'    ,/,
     &'                               ((1+theta).new-theta.old'      ,/,
     &'       ICPEXT = ',4X,I10,    ' (specific heat extrapolation'  ,/,
     &'                ',14X,       ' (0: explicit'                  ,/,
     &'                ',14X,       ' (1: n+thetcp with thetro=1/2'  ,/,
     &'                ',14X,       ' (2: n+thetcp with thetro=1'    ,/,
     &'       THETCP = ', E14.5,    ' (specific heat theta-scheme'   ,/,
     &'                               ((1+theta).new-theta.old'      ,/,
     &'       THETSN = ', E14.5,    ' (Nav-Stokes S.T. theta scheme)',/,
     &'                               ((1+theta).new-theta.old'      ,/,
     &'       THETST = ', E14.5,    ' (Turbulence S.T. theta-scheme)',/,
     &'                               ((1+theta).new-theta.old'      ,/,
     &'       EPSUP  = ', E14.5,    ' (Velocity/pressure coupling'   ,/,
     &'                ',14X,       '  stop test                   )',/)
C
#endif
C
C --- Calcul des gradients
C
      WRITE(NFECRA,4500)
C
      WRITE(NFECRA,4510) IMRGRA, ANOMAX
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,4520) CHAINE(1:8),
     &      NSWRGR(II),NSWRSM(II),EPSRGR(II),EPSRSM(II),EXTRAG(II)
        ENDIF
      ENDDO
      WRITE(NFECRA,4511)
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,4521) CHAINE(1:8),
     &      IRCFLU(II),IMLIGR(II),CLIMGR(II)
        ENDIF
      ENDDO
      WRITE(NFECRA,4530)
C
      WRITE(NFECRA,9900)
C
C --- Interpolation face des viscosites
C
      WRITE(NFECRA,4810) IMVISF
C
      WRITE(NFECRA,9900)
C
C --- Estimateurs d'erreurs pour Navier-Stokes
C
      IIESCA = 0
      DO IPHAS = 1, NPHAS
        DO IEST = 1, NESTMX
          IF(IESCAL(IEST,IPHAS).GT.0) THEN
            IIESCA = 1
          ENDIF
        ENDDO
      ENDDO
C
      IF(IIESCA.GT.0) THEN
        WRITE(NFECRA,4820)
        DO IPHAS = 1, NPHAS
          WRITE(NFECRA,4821)IPHAS
          DO IEST = 1, NESTMX
            WRITE(NFECRA,4822)IEST, IESCAL(IEST,IPHAS)
          ENDDO
          WRITE(NFECRA,4823)
        ENDDO
        WRITE(NFECRA,4824)IESPRE,IESDER,IESCOR,IESTOT
        WRITE(NFECRA,9900)
      ENDIF
C
C --- Calcul de moyennes temporelles
C
      WRITE(NFECRA,4900) NBMOMT
      IF(NBMOMT.GT.0) THEN
        WRITE(NFECRA,4901)
        DO IMOM = 1, NBMOMT
          WRITE(NFECRA,4920)IMOM,IMOOLD(IMOM),
     &         NTDMOM(IMOM),(IDFMOM(JJ,IMOM),JJ=1,NDGMOX)
        ENDDO
        WRITE(NFECRA,4930)
      ENDIF
      WRITE(NFECRA,9900)
C
C
C --- Calcul de la distance a la paroi
C
      IF(INEEDY.EQ.1) THEN
C
        WRITE(NFECRA,4950) ICDPAR
        IF(ABS(ICDPAR).EQ.1) THEN
          WRITE(NFECRA,4951)
     &        NITMAY, NSWRSY, NSWRGY, IMLIGY, IRCFLY, ISCHCY,
     &        ISSTPY, IMGRPY, IWARNY, NTCMXY,
     &        BLENCY, EPSILY, EPSRSY, EPSRGY, CLIMGY, EXTRAY, COUMXY,
     &        EPSCVY, YPLMXY
        ENDIF
        WRITE(NFECRA,9900)
C
      ENDIF
C
C
#if defined(_CS_LANG_FR)
C
 4500 FORMAT(
     &                                                                /,
     &' ** CALCUL DES GRADIENTS'                                     ,/,
     &'    --------------------'                                     ,/)
 4510 FORMAT(
     &'       IMRGRA = ',4X,I10,    ' (Mode de reconstruction      )',/,
     &'       ANOMAX = ',E14.5 ,    ' (Angle de non ortho. limite  )',/,
     &'                               (pour moindres carres etendu )',/,
     &'                                                             ',/,
     &'-----------------------------------------------------------'  ,/,
     &' Variable NSWRGR NSWRSM      EPSRGR      EPSRSM      EXTRAG'  ,/,
     &'-----------------------------------------------------------'    )
 4520 FORMAT(
     & 1X,    A8,    I7,    I7,      E12.4,      E12.4,      E12.4     )
 4511 FORMAT(
     &'-----------------------------------------------'              ,/,
     &                                                                /,
     &'-----------------------------------'                          ,/,
     &' Variable IRCFLU IMLIGR      CLIMGR'                          ,/,
     &'-----------------------------------'                            )
 4521 FORMAT(
     & 1X,    A8,    I7,    I7,      E12.4                             )
 4530 FORMAT(
     &'-----------------------------------'                          ,/,
     &                                                                /,
     &'       NSWRGR =                (nb sweep reconstruction grad)',/,
     &'       NSWRSM =                (nb sweep reconstruction smb )',/,
     &'       EPSRGR =                (precision reconstruction gra)',/,
     &'       EPSRSM =                (precision reconstruction smb)',/,
     &'       EXTRAG =  [0.;1.]       (extrapolation des gradients )',/,
     &'       IRCFLU =  0 ou  1       (reconstruction des flux     )',/,
     &'       IMLIGR =  < 0, 0 ou 1   (methode de limit. des grad  )',/,
     &'       CLIMGR =  > 1 ou 1      (coef de limitation des grad )',/)
C
 4810 FORMAT(
     &                                                                /,
     &' ** INTERPOLATION FACE'                                       ,/,
     &'    ------------------'                                       ,/,
     &'       IMVISF = ',4X,I10,    ' (0 arithmetique              )',/)
C
 4820 FORMAT(
     &                                                                /,
     &' ** ESTIMATEURS D''ERREUR POUR NAVIER-STOKES'                 ,/,
     &'    ----------------------------------------'                 ,/)
 4821 FORMAT(
     &' --- Phase : ',I10                                            ,/,
     &                                                                /,
     &'----------------------------------------'                     ,/,
     &' Estimateur      IESCAL (mode de calcul)'                     ,/,
     &'----------------------------------------'                       )
 4822 FORMAT(
     & 1X,     I10,2X,    I10                                          )
 4823 FORMAT(
     &'----------------------------------------'                       )
 4824 FORMAT(
     &                                                                /,
     &' Estimateurs possibles :'                                     ,/,
     &' ',I2,' =IESPRE : prediction'                                 ,/,
     &'            L''estimateur est base sur la grandeur'           ,/,
     &'            I = rho_n (u*-u_n)/dt + rho_n u_n grad u*'        ,/,
     &'              - rho_n div (mu+mu_t)_n grad u* + grad P_n'     ,/,
     &'              - reste du smb(u_n, P_n, autres variables_n)'   ,/,
c     &'            Idealement nul quand les methodes de reconstruction
c     &'                sont parfaites et le systeme est resolu exactement
     &' ',I2,' =IESDER : derive'                                     ,/,
     &'            L''estimateur est base sur la grandeur'           ,/,
     &'            I = div (flux de masse corrige apres etape'       ,/,
     &'                                               de pression)'  ,/,
     &'            Idealement nul quand l''equation de Poisson est'  ,/,
     &'              resolue exactement'                             ,/,
     &' ',I2,' =IESCOR : correction'                                 ,/,
     &'            L''estimateur est base sur la grandeur'           ,/,
     &'            I = div (rho_n u_(n+1))'                          ,/,
     &'            Idealement nul quand l''equation de Poisson est'  ,/,
     &'              resolue exactement et que le passage des flux'  ,/,
     &'              de masse aux faces vers les vitesses au centre' ,/,
     &'              se fait dans un espace de fonctions'            ,/,
     &'              a divergence nulle'                             ,/,
     &' ',I2,' =IESTOT : total'                                      ,/,
     &'            L''estimateur est base sur la grandeur'           ,/,
     &'            I = rho_n (u_(n+1)-u_n)/dt'                       ,/,
     &'                                 + rho_n u_(n+1) grad u_(n+1)',/,
     &'              - rho_n div (mu+mu_t)_n grad u_(n+1)'           ,/,
     &'                                               + gradP_(n+1)' ,/,
     &'              - reste du smb(u_(n+1), P_(n+1),'               ,/,
     &'                                          autres variables_n)',/,
     &'             Le flux du terme convectif est calcule a partir' ,/,
     &'               de u_(n+1) pris au centre des cellules (et'    ,/,
     &'               non pas a partir du flux de masse aux faces'   ,/,
     &'               actualise)'                                    ,/,
     &                                                                /,
     &' On evalue l''estimateur selon les valeurs de IESCAL :'       ,/,
     &'   IESCAL = 0 : l''estimateur n''est pas calcule'             ,/,
     &'   IESCAL = 1 : l''estimateur    est     calcule,'            ,/,
     &'               sans contribution du volume  (on prend abs(I))',/,
     &'   IESCAL = 2 : l''estimateur    est     calcule,'            ,/,
     &'               avec contribution du volume ("norme L2")'      ,/,
     &'               soit abs(I)*SQRT(Volume_cellule),'             ,/,
     &'               sauf pour IESCOR : on calcule'                 ,/,
     &'                 abs(I)*Volume_cellule pour mesurer'          ,/,
     &'                 un ecart en kg/s'                            ,/)
C
 4900 FORMAT(
     &                                                                /,
     &' ** CALCUL DES MOYENNES TEMPORELLES (MOMENTS)'                ,/,
     &'    -----------------------------------------'                ,/,
     &                                                                /,
     &'       NBMOMT = ',4X,I10,    ' (Nombre de moments           )'  )
 4901 FORMAT(
     &                                                                /,
     &'------------------------------------------------------'       ,/,
     &' IMOM IMOOLD NDTMOM IDFMOM'                                   ,/,
     &'------------------------------------------------------'         )
 4920 FORMAT(
     & 1X,I4,    I7,    I7,5(I7)                                       )
 4930 FORMAT(
     &'------------------------------------------------------'       ,/,
     &                                                                /,
     &'       IMOM   = 0 ou > 0       (numero du moment            )',/,
     &'       IMOOLD =-1 ou > 0       (ancien moment correspondant )',/,
     &'                               (  en suite de calcul ou     )',/,
     &'                               (-1 si le moment est         )',/,
     &'                               (  reinitialise              )',/,
     &'       NDTMOM = 0 ou > 0       (numero du pas de temps de   )',/,
     &'                               (debut de calcul du moment   )',/,
     &'       IDFMOM = 0 ou > 0       (numero des variables        )',/,
     &'                               (composant le moment         )',/)
C
 4950 FORMAT(
     &                                                                /,
     &' ** CALCUL DE LA DISTANCE A LA PAROI'                         ,/,
     &'    --------------------------------'                         ,/,
     &                                                                /,
     &'       ICDPAR = ',4X,I10,    ' ( 1: std et relu      si suite',/,
     &'                               (-1: std et recalcule si suite',/,
     &'                               ( 2: old et relu      si suite',/,
     &'                               (-2: old et recalcule si suite',/)
4951  FORMAT(
     &                                                                /,
     &'       NITMAY = ',4X,I10,    ' (Nb iter pour resolution iter.',/,
     &'       NSWRSY = ',4X,I10,    ' (Nb iter pour reconstr. smb. )',/,
     &'       NSWRGY = ',4X,I10,    ' (Nb iter pour reconstr. grd. )',/,
     &'       IMLIGY = ',4X,I10,    ' (Methode de limitation grd.  )',/,
     &'       IRCFLY = ',4X,I10,    ' (Reconst. flux conv. diff.   )',/,
     &'       ISCHCY = ',4X,I10,    ' (Schema convectif            )',/,
     &'       ISSTPY = ',4X,I10,    ' (Utilisation test de pente   )',/,
     &'       IMGRPY = ',4X,I10,    ' (Algorithme multigrille      )',/,
     &'       IWARNY = ',4X,I10,    ' (Niveau d''impression        )',/,
     &'       NTCMXY = ',4X,I10,    ' (Nb iter pour convection stat.',/,
     &                                                                /,
     &'       BLENCY = ',E14.5 ,    ' (Prop. ordre 2 schema convect.',/,
     &'       EPSILY = ',E14.5 ,    ' (Precision solveur iteratif  )',/,
     &'       EPSRSY = ',E14.5 ,    ' (Precision reconst. smb.     )',/,
     &'       EPSRGY = ',E14.5 ,    ' (Precision reconst. grd.     )',/,
     &'       CLIMGY = ',E14.5 ,    ' (Coeff. pour limitation grd. )',/,
     &'       EXTRAY = ',E14.5 ,    ' (Coeff. pour extrapolation grd',/,
     &'       COUMXY = ',E14.5 ,    ' (Courant max pour convection )',/,
     &'       EPSCVY = ',E14.5 ,    ' (Precision pour convect. stat.',/,
     &'       YPLMXY = ',E14.5 ,    ' (y+ max avec influence amort.)',/)
C
#else
C
 4500 FORMAT(
     &                                                                /,
     &' ** GRADIENTS CALCULATION'                                    ,/,
     &'    ---------------------'                                    ,/)
 4510 FORMAT(
     &'       IMRGRA = ',4X,I10,    ' (Reconstruction mode         )',/,
     &'       ANOMAX = ',E14.5 ,    ' (Non-ortho angle: limit for  )',/,
     &'                               (least squares ext. neighbors)',/,
     &                                                                /,
     &'-----------------------------------------------------------'  ,/,
     &' Variable NSWRGR NSWRSM      EPSRGR      EPSRSM      EXTRAG'  ,/,
     &'-----------------------------------------------------------'    )
 4520 FORMAT(
     & 1X,    A8,    I7,    I7,      E12.4,      E12.4,      E12.4     )
 4511 FORMAT(
     &'-----------------------------------------------'              ,/,
     &                                                                /,
     &'-----------------------------------'                          ,/,
     &' Variable IRCFLU IMLIGR      CLIMGR'                          ,/,
     &'-----------------------------------'                            )
 4521 FORMAT(
     & 1X,    A8,    I7,    I7,      E12.4                             )
 4530 FORMAT(
     &'-----------------------------------'                          ,/,
     &                                                                /,
     &'       NSWRGR =                (nb sweep gradient reconstr. )',/,
     &'       NSWRSM =                (nb sweep rhs reconstrcution )',/,
     &'       EPSRGR =                (grad. reconstruction prec.  )',/,
     &'       EPSRSM =                (rhs   reconstruction prec.  )',/,
     &'       EXTRAG =  [0.;1.]       (gradients extrapolation     )',/,
     &'       IRCFLU =  0 ou  1       (flow reconstruction         )',/,
     &'       IMLIGR =  < 0, 0 ou 1   (gradient limitation method  )',/,
     &'       CLIMGR =  > 1 ou 1      (gradient limitation coeff.  )',/)
C
 4810 FORMAT(
     &                                                                /,
     &' ** FACE INTERPOLATION'                                       ,/,
     &'    ------------------'                                       ,/,
     &'       IMVISF = ',4X,I10,    ' (0 arithmetic                )',/)
C
 4820 FORMAT(
     &                                                                /,
     &' ** ERROR ESTIMATORS FOR NAVIER-STOKES'                       ,/,
     &'    ----------------------------------'                       ,/)
 4821 FORMAT(
     &' --- Phase: ',I10                                             ,/,
     &                                                                /,
     &'------------------------------------------'                   ,/,
     &' Estimateur      IESCAL (calculation mode)'                   ,/,
     &'------------------------------------------'                     )
 4822 FORMAT(
     & 1X,     I10,2X,    I10                                          )
 4823 FORMAT(
     &'----------------------------------------'                       )
 4824 FORMAT(
     &                                                                /,
     &' Possible estimators:'                                        ,/,
     &' ',I2,' =IESPRE: prediction'                                  ,/,
     &'            The estimatore is based on the quantity'          ,/,
     &'            I = rho_n (u*-u_n)/dt + rho_n u_n grad u*'        ,/,
     &'              - rho_n div (mu+mu_t)_n grad u* + grad P_n'     ,/,
     &'              - remainder of rhs(u_n, P_n, other variables_n)'  ,/,
     &' ',I2,' =IESDER: drift'                                       ,/,
     &'            The estimator is based on quantity'               ,/,
     &'            I = div (mass flow corrected after pressure step)',/,
     &'            Ideally zero when Poisson''s equation is'         ,/,
     &'              resolved exactly'                               ,/,
     &' ',I2,' =IESCOR: correction'                                  ,/,
     &'            The estimator is based on quantity'               ,/,
     &'            I = div (rho_n u_(n+1))'                          ,/,
     &'            Ideally zero when Poisson''s equation is'         ,/,
     &'              resolved exactly and the passage from mass flow',/,
     &'              at faces to velocity at cell centers is done'   ,/,
     &'              in a function space with zero divergence'       ,/,
     &' ',I2,' =IESTOT: total'                                       ,/,
     &'            Estimator is based on the quantity'               ,/,
     &'            I = rho_n (u_(n+1)-u_n)/dt'                       ,/,
     &'                                 + rho_n u_(n+1) grad u_(n+1)',/,
     &'              - rho_n div (mu+mu_t)_n grad u_(n+1)'           ,/,
     &'                                               + gradP_(n+1)' ,/,
     &'              - rmainder of rhs(u_(n+1), P_(n+1),'            ,/,
     &'                                          other variables_n)',/,
     &'             The convective term flow is calculated from'     ,/,
     &'               u_(n+1) taken at the cell centers (and not'    ,/,
     &'               from the updated mass flow at faces)'          ,/,
     &                                                                /,
     &' We evaluate the estimator based on values of IESCAL:'        ,/,
     &'   IESCAL = 0: the estimator is not calculated'               ,/,
     &'   IESCAL = 1: the estimator is calculated, with no'          ,/,
     &'               contribution from the volume (we take abs(I))' ,/,
     &'   IESCAL = 2: the estimator is calculated,'                  ,/,
     &'               with contribution from the volume ("L2 norm")' ,/,
     &'               that is abs(I)*SQRT(Cell_volume),'             ,/,
     &'               except for IESCOR: we calculate'               ,/,
     &'                 abs(I)*Cell_volume to measure'               ,/,
     &'                 a difference in kg/s'                        ,/)
C
 4900 FORMAT(
     &                                                                /,
     &' ** CALCULATION OF TEMPORAL MEANS (MOMENTS)'                  ,/,
     &'    ---------------------------------------'                  ,/,
     &                                                                /,
     &'       NBMOMT = ',4X,I10,    ' (Number of moments           )'  )
 4901 FORMAT(
     &                                                                /,
     &'------------------------------------------------------'       ,/,
     &' IMOM IMOOLD NDTMOM IDFMOM'                                   ,/,
     &'------------------------------------------------------'         )
 4920 FORMAT(
     & 1X,I4,    I7,    I7,5(I7)                                       )
 4930 FORMAT(
     &'------------------------------------------------------'       ,/,
     &                                                                /,
     &'       IMOM   = 0 ou > 0       (moment number               )',/,
     &'       IMOOLD =-1 ou > 0       (old moment corresponding    )',/,
     &'                               (  to calculation restart or )',/,
     &'                               (-1 if the moment is         )',/,
     &'                               (  reinitialized             )',/,
     &'       NDTMOM = 0 ou > 0       (moment calculation starting )',/,
     &'                               (time step number            )',/,
     &'       IDFMOM = 0 ou > 0       (number of variables of      )',/,
     &'                               (which the moment is composed)',/)
C
 4950 FORMAT(
     &                                                                /,
     &' ** WALL DISTANCE COMPUTATION'                                ,/,
     &'    -------------------------'                                ,/,
     &                                                                /,
     &'       ICDPAR = ',4X,I10,    ' ( 1: std, reread if restart'   ,/,
     &'                               (-1: std, recomputed if restrt',/,
     &'                               ( 2: old, reread if restart'   ,/,
     &'                               (-2: old, recomputed if restrt',/)
4951  FORMAT(
     &                                                                /,
     &'       NITMAY = ',4X,I10,    ' (Nb iter for iter resolution )',/,
     &'       NSWRSY = ',4X,I10,    ' (Nb iter for rhs reconstr.   )',/,
     &'       NSWRGY = ',4X,I10,    ' (Nb iter for grad. reconstr. )',/,
     &'       IMLIGY = ',4X,I10,    ' (Gradient limitation method  )',/,
     &'       IRCFLY = ',4X,I10,    ' (Conv. Diff. flow reconstr.  )',/,
     &'       ISCHCY = ',4X,I10,    ' (Convective scheme           )',/,
     &'       ISSTPY = ',4X,I10,    ' (Slope tet use               )',/,
     &'       IMGRPY = ',4X,I10,    ' (Multigrid algorithm         )',/,
     &'       IWARNY = ',4X,I10,    ' (Verbosity level             )',/,
     &'       NTCMXY = ',4X,I10,    ' (Nb iter for steady convect. )',/,
     &                                                                /,
     &'       BLENCY = ',E14.5 ,    ' (2nd order conv. scheme prop.)',/,
     &'       EPSILY = ',E14.5 ,    ' (Iterative solver precision  )',/,
     &'       EPSRGY = ',E14.5 ,    ' (Gradient reconstr. precision)',/,
     &'       CLIMGY = ',E14.5 ,    ' (Coeff. for grad. limitation )',/,
     &'       EXTRAY = ',E14.5 ,    ' (Coeff. for grad. extrapolat.)',/,
     &'       COUMXY = ',E14.5 ,    ' (Max CFL for convection      )',/,
     &'       EPSCVY = ',E14.5 ,    ' (Precision for steady conv.  )',/,
     &'       YPLMXY = ',E14.5 ,    ' (y+ max w. damping influence )',/)
C
#endif
C
C
C=======================================================================
C 5. SOLVEURS
C=======================================================================
C
C --- Solveurs iteratifs de base
C
      WRITE(NFECRA,5010)
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,5020) CHAINE(1:8),IRESOL(II),
     &                     NITMAX(II),EPSILO(II),IDIRCL(II)
        ENDIF
      ENDDO
      WRITE(NFECRA,5030)
C
      WRITE(NFECRA,9900)
C
C
C --- Multigrille
C
      WRITE(NFECRA,5510)NCEGRM, NGRMAX
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,5520) CHAINE(1:8),
     &      IMGR(II),NCYMAX(II),NITMGF(II)
        ENDIF
      ENDDO
      WRITE(NFECRA,5530)
C
      WRITE(NFECRA,9900)
C
C
#if defined(_CS_LANG_FR)
C
 5010 FORMAT(
     &                                                                /,
     &' ** SOLVEURS ITERATIFS DE BASE'                               ,/,
     &'    --------------------------'                               ,/,
     &                                                                /,
     &'------------------------------------------'                   ,/,
     &' Variable IRESOL NITMAX      EPSILO IDIRCL'                   ,/,
     &'------------------------------------------'                     )
 5020 FORMAT(
     & 1X,    A8,    I7,    I7,      E12.4,    I7                      )
 5030 FORMAT(
     &'-----------------------------------'                          ,/,
     &                                                                /,
     &'       IRESOL =            -1  (choix automatique du solveur)',/,
     &'                IPOL*1000 + 0  (gradient conjugue           )',/,
     &'                            1  (jacobi                      )',/,
     &'                IPOL*1000 + 2  (bigradient conjugue         )',/,
     &'                  avec IPOL    (degre du preconditionnement )',/,
     &'       NITMAX =                (nb d iterations max         )',/,
     &'       EPSILO =                (precision resolution        )',/,
     &'       IDIRCL = 0 ou 1         (decalage de la diagonale si'  ,/,
     &'                                 ISTAT=0 et pas de Dirichlet)',/)
C
 5510 FORMAT(
     &                                                                /,
     &' ** MULTIGRILLE'                                              ,/,
     &'    -----------'                                              ,/,
     &                                                                /,
     &'       NCEGRM = ',4X,I10,    ' (Nb cell max mail grossier   )',/,
     &'       NGRMAX = ',4X,I10,    ' (Nb max de niveaux de mail   )',/,
     &'------------------------------                               ',/,
     &' Variable   IMGR NCYMAX NITMGF                               ',/,
     &'------------------------------                               '  )
 5520 FORMAT(
     & 1X,    A8,    I7,    I7,    I7                                  )
 5530 FORMAT(
     &'------------------------------'                               ,/,
     &                                                                /,
     &'       IMGR   =  0 ou 1        (1 : activation du mltgrd    )',/,
     &'       NCYMAX =                (Nb max de cycles            )',/,
     &'       NITMGF =                (Nb max d iter sur mail fin  )',/)
C
#else
C
 5010 FORMAT(
     &                                                                /,
     &' ** BASE ITERATIVE SOLVERS'                                   ,/,
     &'    ----------------------'                                   ,/,
     &                                                                /,
     &'------------------------------------------'                   ,/,
     &' Variable IRESOL NITMAX      EPSILO IDIRCL'                   ,/,
     &'------------------------------------------'                     )
 5020 FORMAT(
     & 1X,    A8,    I7,    I7,      E12.4,    I7                      )
 5030 FORMAT(
     &'-----------------------------------'                          ,/,
     &                                                                /,
     &'       IRESOL =            -1  (automatic solver choice     )',/,
     &'                IPOL*1000 + 0  (p conjuguate gradient       )',/,
     &'                            1  (Jacobi                      )',/,
     &'                IPOL*1000 + 2  (bicgstab                    )',/,
     &'                  avec IPOL    (preconditioning degree      )',/,
     &'       NITMAX =                (max number of iterations    )',/,
     &'       EPSILO =                (resolution precision        )',/,
     &'       IDIRCL = 0 ou 1         (shift diagonal if          '  ,/,
     &'                                 ISTAT=0 and no Dirichlet   )',/)
C
 5510 FORMAT(
     &                                                                /,
     &' ** MULTIGRID'                                                ,/,
     &'    ---------'                                                ,/,
     &                                                                /,
     &'       NCEGRM = ',4X,I10,    ' (Max nb cells coarsest grid  )',/,
     &'       NGRMAX = ',4X,I10,    ' (Max numeber of levels       )',/,
     &'------------------------------                               ',/,
     &' Variable   IMGR NCYMAX NITMGF                               ',/,
     &'------------------------------                               '  )
 5520 FORMAT(
     & 1X,    A8,    I7,    I7,    I7                                  )
 5530 FORMAT(
     &'------------------------------'                               ,/,
     &                                                                /,
     &'       IMGR   =  0 ou 1        (1: multigrid activated      )',/,
     &'       NCYMAX =                (Max number  of cycles       )',/,
     &'       NITMGF =                (Max nb iter on coarsest grid)',/)
C
#endif
C
C=======================================================================
C 6. SCALAIRES
C=======================================================================
C
C --- Scalaires
C
      IF(NSCAL.GE.1) THEN
        WRITE(NFECRA,6000)
        WRITE(NFECRA,6010)ITBRRB
        WRITE(NFECRA,6011)
        DO II = 1, NSCAL
          CHAINE=NOMVAR(IPPRTP(ISCA(II)))
          WRITE(NFECRA,6021) CHAINE(1:8),II,IPHSCA(II),ISCSTH(II),
     &                       IVISLS(II),VISLS0(II),SIGMAS(II)
        ENDDO
        WRITE(NFECRA,6031)
        WRITE(NFECRA,6012)
        DO II = 1, NSCAL
          CHAINE=NOMVAR(IPPRTP(ISCA(II)))
          WRITE(NFECRA,6022) CHAINE(1:8),II,IPHSCA(II),ISCAVR(II),
     &                       RVARFL(II)
        ENDDO
        WRITE(NFECRA,6032)
        WRITE(NFECRA,6013)
        DO II = 1, NSCAL
          CHAINE=NOMVAR(IPPRTP(ISCA(II)))
          WRITE(NFECRA,6023) CHAINE(1:8),II,IPHSCA(II),ICLVFL(II),
     &                       SCAMIN(II),SCAMAX(II)
        ENDDO
        WRITE(NFECRA,6033)
        WRITE(NFECRA,6030)
        WRITE(NFECRA,6040)
        DO II = 1, NSCAL
          WRITE(NFECRA,6041) II,THETSS(II),IVSEXT(II),THETVS(II)
        ENDDO
        WRITE(NFECRA,6042)
C
        WRITE(NFECRA,9900)
C
      ENDIF
C
C
#if defined(_CS_LANG_FR)
C
 6000 FORMAT(
     &                                                                /,
     &' ** SCALAIRES'                                                ,/,
     &'    ---------'                                                ,/)
 6010 FORMAT(
     &'       ITBRRB = ',4X,I10,    ' (Reconstruction T ou H au brd)',/)
 6011 FORMAT(
     &'-------------------------------------------------------------',/,
     &' Variable Numero IPHSCA ISCSTH IVISLS      VISLS0      SIGMAS',/,
     &'-------------------------------------------------------------'  )
 6021 FORMAT(
     & 1X,    A8,    I7,    I7,    I7,    I7,      E12.4,      E12.4   )
 6031 FORMAT(
     &'-------------------------------------------------------------',/)
 6012 FORMAT(
     &'------------------------------------------'                   ,/,
     &' Variable Numero IPHSCA ISCAVR      RVARFL'                   ,/,
     &'------------------------------------------'                     )
 6022 FORMAT(
     & 1X,    A8,    I7,    I7,    I7,      E12.4                      )
 6032 FORMAT(
     &'------------------------------------------'                   ,/)
 6013 FORMAT(
     &'------------------------------------------------------'       ,/,
     &' Variable Numero IPHSCA ICLVFL      SCAMIN      SCAMAX'       ,/,
     &'------------------------------------------------------'         )
 6023 FORMAT(
     & 1X,    A8,    I7,    I7,    I7,      E12.4,      E12.4          )
 6033 FORMAT(
     &'------------------------------------------------------'       ,/)
 6030 FORMAT(
     &'-------------------------------------------------------------',/,
     &                                                                /,
     &'       Le numero indique pour chaque scalaire le rang'        ,/,
     &'         dans la liste de tous les scalaires. Les scalaires'  ,/,
     &'         utilisateurs sont places en tete, de 1 a NSCAUS. Les',/,
     &'         scalaires physique particuliere sont a la fin, de'   ,/,
     &'         NSCAUS+1 a NSCAPP+NSCAUS=NSCAL.'                     ,/,
     &                                                                /,
     &'       IPHSCA =                (Phase porteuse              )',/,
     &'       ISCSTH = -1,0, 1 ou 2   (T (C), Passif, T (K) ou H   )',/,
     &'       IVISLS = 0 ou >0        (Viscosite constante ou non  )',/,
     &'       VISLS0 = >0             (Viscosite de reference      )',/,
     &'       SIGMAS = >0             (Schmidt                     )',/,
     &'       ISCAVR = 0 ou >0        (Scalaire associe si variance)',/,
     &'       RVARFL = >0             (Rf, cf dissipation variance )',/,
     &'       ICLVFL = 0, 1 ou 2      (Mode de clipping variance   )',/,
     &'       SCAMIN =                (Valeur min autorisee        )',/,
     &'       SCAMAX =                (Valeur max autorisee        )',/,
     &'        Pour les variances, SCAMIN est ignore et SCAMAX n est',/,
     &'          pris en compte que si ICLVFL = 2'                   ,/)
 6040 FORMAT(
     &'------------------------------------------------------'       ,/,
     &'   Scalaire      THETSS    IVSEXT      THETVS'                ,/,
     &'------------------------------------------------------'         )
 6041 FORMAT(
     & 1X,     I10,      E12.4,      I10,      E12.4                   )
 6042 FORMAT(
     &'------------------------------------------------------'       ,/,
     &                                                                /,
     &'       THETSS =                (theta pour termes sources   )',/,
     &'                               ((1+theta)nouveau-theta ancien',/,
     &'       IVSEXT =                (extrap. viscosite totale    )',/,
     &'                               (0 : explicite               )',/,
     &'                               (1 : n+thetvs avec thetvs=1/2' ,/,
     &'                               (2 : n+thetvs avec thetvs=1  )',/,
     &'       THETVS =                (theta pour diffusiv. scalaire',/,
     &'                               ((1+theta)nouveau-theta ancien',/)
C
#else
C
 6000 FORMAT(
     &                                                                /,
     &' ** SCALARS'                                                  ,/,
     &'    -------'                                                  ,/)
 6010 FORMAT(
     &'       ITBRRB = ',4X,I10,    ' (T or H reconstruction at bdy)',/)
 6011 FORMAT(
     &'-------------------------------------------------------------',/,
     &' Variable Number IPHSCA ISCSTH IVISLS      VISLS0      SIGMAS',/,
     &'-------------------------------------------------------------'  )
 6021 FORMAT(
     & 1X,    A8,    I7,    I7,    I7,    I7,      E12.4,      E12.4   )
 6031 FORMAT(
     &'-------------------------------------------------------------',/)
 6012 FORMAT(
     &'------------------------------------------'                   ,/,
     &' Variable Number IPHSCA ISCAVR      RVARFL'                   ,/,
     &'------------------------------------------'                     )
 6022 FORMAT(
     & 1X,    A8,    I7,    I7,    I7,      E12.4                      )
 6032 FORMAT(
     &'------------------------------------------'                   ,/)
 6013 FORMAT(
     &'------------------------------------------------------'       ,/,
     &' Variable Number IPHSCA ICLVFL      SCAMIN      SCAMAX'       ,/,
     &'------------------------------------------------------'         )
 6023 FORMAT(
     & 1X,    A8,    I7,    I7,    I7,      E12.4,      E12.4          )
 6033 FORMAT(
     &'------------------------------------------------------'       ,/)
 6030 FORMAT(
     &'-------------------------------------------------------------',/,
     &                                                                /,
     &'       For each scalar, the number indicates it''s rank'      ,/,
     &'         in the list of all scalars. User scalars are placed' ,/,
     &'         first, from 1 to NSCAUS. Specific physics scalars'   ,/,
     &'         are placed at the end, from'                         ,/,
     &'         NSCAUS+1 to NSCAPP+NSCAUS=NSCAL.'                    ,/,
     &                                                                /,
     &'       IPHSCA =                (Carrier phase               )',/,
     &'       ISCSTH = -1,0, 1 ou 2   (T (C), Passive, T (K) or H  )',/,
     &'       IVISLS = 0 ou >0        (Viscosity: constant or not  )',/,
     &'       VISLS0 = >0             (Reference viscosity         )',/,
     &'       SIGMAS = >0             (Schmidt                     )',/,
     &'       ISCAVR = 0 ou >0        (Associat. scalar if variance)',/,
     &'       RVARFL = >0             (Rf, cf variance dissipation )',/,
     &'       ICLVFL = 0, 1 ou 2      (Variance clipping mode      )',/,
     &'       SCAMIN =                (Min authorized value        )',/,
     &'       SCAMAX =                (Max authorized value        )',/,
     &'        For variances, SCAMIN is ignored and SCAMAX is used'  ,/,
     &'          only if ICLVFL = 2'                                 ,/)
 6040 FORMAT(
     &'------------------------------------------------------'       ,/,
     &'   Scalar        THETSS    IVSEXT      THETVS'                ,/,
     &'------------------------------------------------------'         )
 6041 FORMAT(
     & 1X,     I10,      E12.4,      I10,      E12.4                   )
 6042 FORMAT(
     &'------------------------------------------------------'       ,/,
     &                                                                /,
     &'       THETSS =                (theta for source terms      )',/,
     &'                               ((1+theta).new-theta.old     )',/,
     &'       IVSEXT =                (extrap. total viscosity     )',/,
     &'                               (0: explicit                 )',/,
     &'                               (1: n+thetvs with thetvs=1/2 )',/,
     &'                               (2: n+thetvs with thetvs=1   )',/,
     &'       THETVS =                (theta pour diffusiv. scalaire',/,
     &'                               ((1+theta)nouveau-theta ancien',/)
C
#endif
C
C=======================================================================
C 7. GESTION DU CALCUL
C=======================================================================
C
C --- Gestion du calcul
C
      WRITE(NFECRA,7000)
C
C   - Suite de calcul
C
      WRITE(NFECRA,7010) ISUITE, ILEAUX, IECAUX
      IF(ISUITE.EQ.1.AND.NSCAL.GE.1) THEN
        WRITE(NFECRA,7020)
        DO II = 1, NSCAL
          IVAR = ISCA(II)
          CHAINE=NOMVAR(IPPRTP(IVAR))
          WRITE(NFECRA,7030) CHAINE(1:8),II,ISCOLD(II)
        ENDDO
        WRITE(NFECRA,7040)
      ENDIF
C
C   - Duree du calcul
C
      WRITE(NFECRA,7110) INPDT0,NTMABS
C
C   - Marge en temps CPU
C
      WRITE(NFECRA,7210) TMARUS
C
      WRITE(NFECRA,9900)
C
C
#if defined(_CS_LANG_FR)
C
 7000 FORMAT(
     &                                                                /,
     &' ** GESTION DU CALCUL'                                        ,/,
     &'    -----------------'                                        ,/)
 7010 FORMAT(
     &' --- Suite de calcul'                                         ,/,
     &'       ISUITE = ',4X,I10,    ' (1 : suite de calcul         )',/,
     &'       ILEAUX = ',4X,I10,    ' (1 : lecture  de suiamx aussi)',/,
     &'       IECAUX = ',4X,I10,    ' (1 : ecriture de suiavx aussi)',/,
     &                                                                /,
     &'       suiamx et suiavx sont les fichiers suite auxiliaires.)',/)
 7020 FORMAT(
     &'       ISCOLD(I) : Dans le calcul precedent, numero du'       ,/,
     &'                   scalaire correspondant au scalaire I du'   ,/,
     &'                   calcul courant :'                          ,/,
     &                                                                /,
     &'---------------------------------------------------'          ,/,
     &' Scalaire Numero   <-   Numero de l''ancien scalaire'         ,/,
     &'---------------------------------------------------'            )
 7030 FORMAT(
     & 1X,    A8,    I7,    7X,    I7                                  )
 7040 FORMAT(
     &'---------------------------------------------------'          ,/,
     &                                                                /,
     &'   La table precedente (ISCOLD) donne la correspondance des'  ,/,
     &'     scalaires du calcul courant avec ceux du calcul'         ,/,
     &'     precedent apres intervention eventuelle de l utilisateur',/,
     &'   Il s''agit de numeros de 1 a NSCAL qui reperent le'        ,/,
     &'     scalaire dans la liste de tous les scalaires'            ,/,
     &'     utilisateur+physique particuliere.'                      ,/,
     &'     .-999 est la valeur par defaut si l utilisateur n est'   ,/,
     &'           intervenu. Le numero des correspondant sera'       ,/,
     &'           complete a la lecture du fichier suite, selon'     ,/,
     &'           le nombre de scalaires disponibles, en utilisant'  ,/,
     &'           la loi (nouveau scalaire ii <- ancien scalaire ii)',/,
     &'     .   0 si l utilisateur souhaite que le scalaire du'      ,/,
     &'           calcul courant n ait pas de correspondant (i.e.'   ,/,
     &'           soit un nouveau scalaire).'                        ,/,
     &'     .   n > 0 si l utilisateur souhaite que le scalaire du'  ,/,
     &'           calcul courant ait pour correspondant le scalaire' ,/,
     &'           n du calcul precedent.'                            ,/)
 7110 FORMAT(
     &' --- Duree du calcul'                                         ,/,
     &'     La numerotation des pas de temps et la mesure du temps'  ,/,
     &'       physique simule sont des valeurs absolues'             ,/,
     &'       et non pas des valeurs relatives au calcul en cours.'  ,/,
     &                                                                /,
     &'       INPDT0 = ',4X,I10,    ' (1 : calcul a zero pas de tps)',/,
     &'       NTMABS = ',4X,I10,    ' (Pas de tps final demande    )',/)
 7210 FORMAT(
     &' --- Marge en temps CPU'                                      ,/,
     &'       TMARUS = ', E14.5,    ' (Marge CPU avant arret       )',/)
C
#else
C
 7000 FORMAT(
     &                                                                /,
     &' ** CALCULATION MANAGEMENT'                                   ,/,
     &'    ----------------------'                                   ,/)
 7010 FORMAT(
     &' --- Restarted calculation'                                   ,/,
     &'       ISUITE = ',4X,I10,    ' (1: restarted calculuation   )',/,
     &'       ILEAUX = ',4X,I10,    ' (1: also read  suiamx        )',/,
     &'       IECAUX = ',4X,I10,    ' (1: also write suiavx        )',/,
     &                                                                /,
     &'       suiamx and suiavx are the auxiliary restart files.'    ,/)
 7020 FORMAT(
     &'       ISCOLD(I): In the previous calculation, number of'     ,/,
     &'                  the scalar corresponding to scalar I in'    ,/,
     &'                  the current calculation:'                   ,/,
     &                                                                /,
     &'---------------------------------------------------'          ,/,
     &' Scalar   Number   <-   Old scalar number'                    ,/,
     &'---------------------------------------------------'            )
 7030 FORMAT(
     & 1X,    A8,    I7,    7X,    I7                                  )
 7040 FORMAT(
     &'---------------------------------------------------'          ,/,
     &                                                                /,
     &'   The preceding array (ISCOLD) defines the correspondance'   ,/,
     &'     of scalars in the current calculation with those of the' ,/,
     &'     previous calculation after possible user intervention.'  ,/,
     &'   It consists of numbers 1 to NSCAL which locate the scalar' ,/,
     &'     in the list of all user+specific physics scalars.'       ,/,
     &'     .-999 is the default value if the user has not'          ,/,
     &'           intervened. The correspondant''s numbers will'     ,/,
     &'           be completed when reading the restart file, based' ,/,
     &'           on the number of available scalars, using the'     ,/,
     &'           rule (new scalar ii <- old scalar ii)'             ,/,
     &'     .   0 if the user does not want the scalar of the'       ,/,
     &'           current calculation to have a correspondant (i.e.' ,/,
     &'           it is a new scalar).'                              ,/,
     &'     .   n > 0 if the user wants the scalar of the current'   ,/,
     &'           calculation to correspond to scalar n in the'      ,/,
     &'           previous calculation.'                             ,/)
 7110 FORMAT(
     &' --- Calculation time'                                        ,/,
     &'     The numbering of time steps and the measure of simulated',/,
     &'       physical time are absolute values, and not values'     ,/,
     &'       relative to the current calculation.'                  ,/,
     &                                                                /,
     &'       INPDT0 = ',4X,I10,    ' (1: 0 time step calcuation   )',/,
     &'       NTMABS = ',4X,I10,    ' (Final time step required    )',/)
 7210 FORMAT(
     &' --- CPU time margin'                                         ,/,
     &'       TMARUS = ', E14.5,    ' (CPU time margin before stop )',/)
C
#endif
C
C=======================================================================
C 8. ENTREES SORTIES
C=======================================================================
C
      WRITE(NFECRA,7500)
C
C   - Fichier suite
C
      WRITE(NFECRA,7510) NTSUIT
C
C   - Fichiers Ensight
C
      WRITE(NFECRA,7520) NTCHR
      DO II = 2, NVPPMX
        IF(ICHRVR(II).EQ.1) THEN
          NAME = NOMVAR(II)
          WRITE(NFECRA,7521) II,NAME  (1:16)
        ENDIF
      ENDDO
      WRITE(NFECRA,7522)
C
C   - Fichiers historiques
C
      WRITE(NFECRA,7530) NTHIST,NCAPT,NTHSAV
      DO II = 2, NVPPMX
        IF(IHISVR(II,1).NE.0) THEN
          NAME = NOMVAR(II)
          WRITE(NFECRA,7531) II,NAME  (1:16),IHISVR(II,1)
        ENDIF
      ENDDO
      WRITE(NFECRA,7532)
C
C   - Fichiers listing
C
      WRITE(NFECRA,7540) NTLIST
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IWAR = IWARNI(II)
        ELSE
          IWAR = -999
        ENDIF
        IF(ILISVR(IPP).EQ.1) THEN
          NAME = NOMVAR(IPP)
          WRITE(NFECRA,7531) IPP,NAME  (1:16),IWAR
        ENDIF
      ENDDO
      WRITE(NFECRA,7532)
C
C   - Post-traitement automatique (bord)
C
      WRITE(NFECRA,7550)   'IPSTDV',IPSTDV,
     &                     'IPSTYP',IPSTYP,
     &                     'IPSTCL',IPSTCL,
     &                     'IPSTFT',IPSTFT,
     &                     'IPSTFO',IPSTFO,
     &                     'IPSTDV'
C
      WRITE(NFECRA,9900)
C
C
#if defined(_CS_LANG_FR)
C
 7500 FORMAT(
     &                                                                /,
     &' ** ENTREES SORTIES'                                          ,/,
     &'    ---------------'                                          ,/)
 7510 FORMAT(
     &' --- Fichier suite'                                           ,/,
     &'       NTSUIT = ',4X,I10,    ' (Periode de sauvegarde)'       ,/)
 7520 FORMAT(
     &' --- Fichiers EnSight, MED, ou CGNS'                          ,/,
     &'       NTCHR  = ',4X,I10,    ' (Periode de sortie    )'       ,/,
     &                                                                /,
     &'       Numero Nom'                                              )
 7521 FORMAT(
     &'   ',     I10,1X,          A16                                  )
 7522 FORMAT(
     &'         --           --'                                     ,/)
 7530 FORMAT(
     &' --- Fichiers historiques'                                    ,/,
     &'       NTHIST = ',4X,I10,    ' (Periode de sortie    )'       ,/,
     &'       NCAPT  = ',4X,I10,    ' (Nombre de capteurs   )'       ,/,
     &'       NTHSAV = ',4X,I10,    ' (Periode de sauvegarde)'       ,/,
     &                                                                /,
     &'       Numero Nom                   Nb. sondes (-1 : toutes)'   )
 7531 FORMAT(
     &'   ',     I10,1X,          A16,6X,         I10                  )
 7532 FORMAT(
     &'         --           --                --'                   ,/)
 7540 FORMAT(
     &' --- Fichiers listing'                                        ,/,
     &'       NTLIST = ',4X,I10,    ' (Periode de sortie    )'       ,/,
     &                                                                /,
     &'       Numero Nom                 Niveau d''impression IWARNI',/,
     &'                                      (-999 : non applicable)',/)
 7550 FORMAT(
     &' --- Variables supplementaires en post-traitement'            ,/,
     &'       ',A6,' = ',4X,I10,    ' (Produit des valeurs suivantes',/,
     &'                                selon activation ou non'      ,/,
     &'       ',A6,' = ',4X,I10,    ' (Yplus          au bord'       ,/,
     &'       ',A6,' = ',4X,I10,    ' (Variables      au bord'       ,/,
     &'       ',A6,' = ',4X,I10,    ' (Flux thermique au bord'       ,/,
     &'       ',A6,' = ',4X,I10,    ' (Force exercee  au bord'       ,/,
     &'  et   ',A6,' =              1 (Pas de sortie supplementaire' ,/)
C
#else
C
 7500 FORMAT(
     &                                                                /,
     &' ** INPUT-OUTPUT'                                             ,/,
     &'    ------------'                                             ,/)
 7510 FORMAT(
     &' --- Restart file'                                            ,/,
     &'       NTSUIT = ',4X,I10,    ' (Checkpoint frequency )'       ,/)
 7520 FORMAT(
     &' --- EnSight, MED, or CGNS output'                            ,/,
     &'       NTCHR  = ',4X,I10,    ' (Output frequency     )'       ,/,
     &                                                                /,
     &'       Number Name'                                             )
 7521 FORMAT(
     &'   ',     I10,1X,          A16                                  )
 7522 FORMAT(
     &'         --           --'                                     ,/)
 7530 FORMAT(
     &' --- Probe history files'                                     ,/,
     &'       NTHIST = ',4X,I10,    ' (Output frequency     )'       ,/,
     &'       NCAPT  = ',4X,I10,    ' (Number of probes     )'       ,/,
     &'       NTHSAV = ',4X,I10,    ' (Checkpoint frequency )'       ,/,
     &                                                                /,
     &'       Number Name                  Nb. probes (-1: all)'       )
 7531 FORMAT(
     &'   ',     I10,1X,          A16,6X,         I10                  )
 7532 FORMAT(
     &'         --           --                --'                   ,/)
 7540 FORMAT(
     &' --- Log files'                                               ,/,
     &'       NTLIST = ',4X,I10,    ' (Output frequency     )'       ,/,
     &                                                                /,
     &'       Number Name                IWARNI verbosity level'     ,/,
     &'                                      (-999: not applicable)' ,/)
 7550 FORMAT(
     &' --- Additional post-processing variables'                    ,/,
     &'       ',A6,' = ',4X,I10,    ' (Product of the following     ',/,
     &'                                values based on activation  )',/,
     &'       ',A6,' = ',4X,I10,    ' (Yplus          on boundary  )',/,
     &'       ',A6,' = ',4X,I10,    ' (Variables      on boundary  )',/,
     &'       ',A6,' = ',4X,I10,    ' (Thermal flow   on boundary  )',/,
     &'       ',A6,' = ',4X,I10,    ' (Force exerted  on boundary  )',/,
     &'  and  ',A6,' =              1 (No additional output        )',/)
C
#endif
C
C=======================================================================
C 9. FICHIERS
C=======================================================================
C
C --- Fichiers
C
      WRITE(NFECRA,7600)
      WRITE(NFECRA,7610)FICGEO,FICAMO,FICAMX,FICSTP,
     &                  IMPGEO,IMPSTP
      WRITE(NFECRA,7611)FICAVA,FICAVX
      WRITE(NFECRA,7612)FICFPP,FICAMR,FICJNF,FICAVR,
     &                  IMPFPP,       IMPJNF
      WRITE(NFECRA,7613)FICAML, FICMLS, FICAVL, FICVLS, FICLAL,
     &                  IMPLAL, IMPLI1, IMPLI2
C
C Autres fichiers lagrangien
      WRITE(NFECRA,7615) IMPLA1,IMPLA2,IMPLA3,IMPLA4,IMPLA5(1),
     &                   IMPLA5(2),IMPLA5(3),IMPLA5(4),
     &                   IMPLA5(5),IMPLA5(6),IMPLA5(7),
     &                   IMPLA5(8),IMPLA5(9),IMPLA5(10),
     &                   IMPLA5(11),IMPLA5(12),IMPLA5(13),
     &                   IMPLA5(14),IMPLA5(15)
C
C Fichiers hist
      WRITE(NFECRA,7620)(IMPHIS(II),II=1,2     )
C
      WRITE(NFECRA,7630) NUSHMX
      WRITE(NFECRA,7631)(FICUSH(II),II=1,NUSHMX)
      WRITE(NFECRA,7632)(IMPUSH(II),II=1,NUSHMX)
C
C
C Fichiers utilisateurs
      WRITE(NFECRA,7640) NUSRMX
      WRITE(NFECRA,7631)(FICUSR(II),II=1,NUSRMX)
      WRITE(NFECRA,7632)(IMPUSR(II),II=1,NUSRMX)
C
      WRITE(NFECRA,7650) IFOENV
C
      WRITE(NFECRA,9900)
C
C
#if defined(_CS_LANG_FR)
C
 7600 FORMAT(
     &                                                                /,
     &' ** FICHIERS'                                                 ,/,
     &'    --------'                                                 ,/)
 7610 FORMAT(
     &' --- Fichiers standard'                                       ,/,
     &                                                                /,
     &'         GEOMETRIE AMONT PPAL  AMONT AUX       STOP'          ,/,
     &'Nom    ',     4(5X,A6)                                        ,/,
     &'Unite  ',     I11,        11X,       11X,       I11           ,/)
 7611 FORMAT(
     &'         AVAL PPAL   AVAL AUX'                                ,/,
     &'Nom    ',     2(5X,A6)                                        ,/)
 7612 FORMAT(
     &' --- Fichiers rayonnement'                                    ,/,
     &                                                                /,
     &'         DONNEES    AMONT    JANAF     AVAL                  ',/,
     &'Nom    ',     4(3X,A6)                                        ,/,
     &'Unite  ',     I9,   9X,     I9                                ,/)
 7613 FORMAT(
     &' --- Fichiers Lagrangiens'                                    ,/,
     &                                                                /,
     &'         ---- AMONT ----    ---- AVAL -----    ---- POST ----',/,
     &'         CALCUL    STAT.    CALCUL    STAT.    LISTING  HISTO',/,
     &'Nom   ',     2(3X,A6),  1X,      2(3X,A6),   4X,       A6     ,/,
     &'Unite ',         18X ,  1X,           18X,   6X,     3 I4     ,/)
 7615 FORMAT(
     &' --- Autres fichiers pour le module Lagrangien'               ,/,
     &                                                                /,
     &'Unite  ',     5I9                                             ,/,
     &'       ',     5I9                                             ,/,
     &'       ',     5I9                                             ,/,
     &'       ',     4I9                                             ,/)
C
 7620 FORMAT(
     &' --- Fichiers developpeurs pour historiques'                  ,/,
     &'                IMPHIS(1)         IMPHIS(2)'                  ,/,
     &'Unite  ',     2(9X,I9)                                        ,/)
 7630 FORMAT(
     &' --- Fichiers utilisateurs pour historiques'                  ,/,
     &'                 ',I10   ,' fichiers'                         ,/,
     &' Nom et unite'                                                  )
 7631 FORMAT(
     &               6(3X,A6)                                          )
 7632 FORMAT(
     &               6I9                                               )
 7640 FORMAT(
     &                                                                /,
     &' --- Fichiers utilisateurs libres'                            ,/,
     &'                 ',I10   ,' fichiers'                         ,/,
     &' Nom et unite'                                                  )
 7650 FORMAT(
     &                                                                /,
     &' --- Fichiers preprocesseur'                                  ,/,
     &'   Formats : 1=Lecture dus fichiers preprocessor_output'      ,/,
     &'             0=Solveur autonome'                              ,/,
     &                                                                /,
     &'Format ',I10                                                  ,/)
C
#else
C
 7600 FORMAT(
     &                                                                /,
     &' ** FILES'                                                    ,/,
     &'    -----'                                                    ,/)
 7610 FORMAT(
     &' --- Standard files'                                          ,/,
     &                                                                /,
     &'         GEOMETRY  MAIN RESTART  AUX RESTART   STOP'          ,/,
     &'Name   ',     4(5X,A6)                                        ,/,
     &'Unit   ',     4I11                                            ,/,
     &'Format ',     I11,      7X,A4,     7X,A4,       I11           ,/)
 7611 FORMAT(
     &'         MAIN CHECKPOINT  AUX CHECKPOINT'                     ,/,
     &'Name   ',     2(5X,A6)                                        ,/,
     &'Unit   ',     2I11                                            ,/,
     &'Format ',     2I11                                            ,/)
 7612 FORMAT(
     &' --- Radiative files'                                         ,/,
     &                                                                /,
     &'         DATA       RESTART  JANAF                 CHECKPOINT',/,
     &'Name   ',     3(3X,A6),     9X, 9X, (3X,A6)                   ,/,
     &'Unit   ',     3    I9 ,     9X, 9X,     I9                    ,/,
     &'Format ',     I9,3X,A4,2X,      I9,            9X, 9X,     I9 ,/)
 7613 FORMAT(
     &' --- Lagrangian files'                                        ,/,
     &                                                                /,
     &'         --- RESTART ---    -- CHECKPOINT -    ---- POST ----',/,
     &'         CALCUL.   STAT.    CALCUL.   STAT.    LOG      HIST.',/,
     &'Name  ',     2(3X,A6),  1X,      2(3X,A6),   4X,       A6     ,/,
     &'Unit  ',     2    I9 ,  1X,       2    I9,   6X,     3 I4     ,/,
     &'Format',5X,A4,5X,A4,2X, 4X,      I4,5X,I4,   6X,     3 I4     ,/)
 7615 FORMAT(
     &' --- Other files for Lagrangian module'                       ,/,
     &                                                                /,
     &'Unit   ',     5I9                                             ,/,
     &'       ',     5I9                                             ,/,
     &'       ',     5I9                                             ,/,
     &'       ',     4I9                                             ,/)
C
 7620 FORMAT(
     &' --- Developper files for probe history'                      ,/,
     &'                IMPHIS(1)         IMPHIS(2)'                  ,/,
     &'Unit   ',     2(9X,I9)                                        ,/)
 7630 FORMAT(
     &' --- User files for probe history'                            ,/,
     &'                 ',I10   ,' files'                            ,/,
     &' Name and unit'                                                 )
 7631 FORMAT(
     &               6(3X,A6)                                          )
 7632 FORMAT(
     &               6I9                                               )
 7640 FORMAT(
     &                                                                /,
     &' --- Free user files'                                         ,/,
     &'                 ',I10   ,' files'                            ,/,
     &' Name and unit'                                                 )
 7650 FORMAT(
     &                                                                /,
     &' --- Preprocessor files'                                      ,/,
     &'   Formats: 1=Read preprocessor_output files'                 ,/,
     &'            0=Standalone solver'                              ,/,
     &                                                                /,
     &'Format ',I10                                                  ,/)
C
#endif
C
C
C=======================================================================
C 10. COUPLAGES
C=======================================================================
C
C
C --- Couplage SYRTHES
C
C     RECUPERATION DU NOMBRE DE CAS DE COUPLAGE
C
      CALL NBCSYR (NBCCOU)
C     ===========
C
      IF (NBCCOU .GE. 1) THEN
C
        WRITE (NFECRA, 8000)
        WRITE (NFECRA, 8010) NBCCOU
C
        WRITE(NFECRA,8030)
        DO II = 1, NSCAL
          CHAINE=NOMVAR(IPPRTP(ISCA(II)))
          WRITE(NFECRA,8031) CHAINE(1:8),II,ICPSYR(II)
        ENDDO
        WRITE(NFECRA,8032)
C
        WRITE(NFECRA,9900)
C
      ENDIF
C
C
#if defined(_CS_LANG_FR)
C
 8000 FORMAT(
     &                                                                /,
     &' ** COUPLAGE SYRTHES'                                         ,/,
     &'    ----------------'                                         ,/)
 8010 FORMAT(
     &'       NBCCOU = ',4X,I10,    ' (Nombre de couplages         )',/)
 8030 FORMAT(
     &                                                                /,
     &'  -- Scalaires couples'                                       ,/,
     &'-----------------------'                                      ,/,
     &' Scalaire Numero ICPSYR'                                      ,/,
     &'-----------------------'                                        )
 8031 FORMAT(
     & 1X,    A8,    I7,    I7                                         )
 8032 FORMAT(
     &'-----------------------'                                      ,/,
     &                                                                /,
     &'       ICPSYR = 0 ou 1         (1 : scalaire couple SYRTHES )',/)
C
#else
C
 8000 FORMAT(
     &                                                                /,
     &' ** SYRTHES COUPLING'                                         ,/,
     &'    ----------------'                                         ,/)
 8010 FORMAT(
     &'       NBCCOU = ',4X,I10,    ' (Number of couplings         )',/)
 8030 FORMAT(
     &                                                                /,
     &'  -- Coupled scalars'                                         ,/,
     &'-----------------------'                                      ,/,
     &' Scalar   Number ICPSYR'                                      ,/,
     &'-----------------------'                                        )
 8031 FORMAT(
     & 1X,    A8,    I7,    I7                                         )
 8032 FORMAT(
     &'-----------------------'                                      ,/,
     &                                                                /,
     &'       ICPSYR = 0 or 1         (1: scalar coupled to SYRTHES)',/)
C
#endif
C
C=======================================================================
C 11. Lagrangien
C=======================================================================
C
C --- Lagrangien
C
      IF (IILAGR.NE.0) THEN
        WRITE(NFECRA,8100) ILPHAS, IILAGR, ISUILA, ISUIST, IPHYLA
C
        IF (IPHYLA.EQ.1) THEN
          WRITE(NFECRA,8105) IDPVAR, ITPVAR, IMPVAR
        ENDIF
C
        WRITE(NFECRA,8106) NBPMAX, NVLS, ISTTIO, INJCON, IROULE
C
        IF (IPHYLA.EQ.2) THEN
          WRITE(NFECRA,8111) IENCRA
          DO II = 1 , NCHARB
            WRITE(NFECRA,8112) II, TPRENC(II),  II
          ENDDO
          DO II = 1 , NCHARB
            WRITE(NFECRA,8113) II, VISREF(II), II
          ENDDO
        ENDIF
C
        IF (IILAGR.EQ.2) THEN
          WRITE(NFECRA,8120) NSTITS, LTSDYN, LTSMAS, LTSTHE
        ENDIF
C
        WRITE(NFECRA,8130) ISTALA
        IF (ISTALA.EQ.1) THEN
          WRITE(NFECRA,8135) SEUIL, IDSTNT, NSTIST, NVLSTS
        ENDIF
C
        WRITE(NFECRA,8140) IDISTU, IDIFFL, MODCPL
C
        IF (MODCPL.GT.0) WRITE(NFECRA,8141) IDIRLA
C
        WRITE(NFECRA,8142) NORDRE, ILAPOI
C
        WRITE(NFECRA,8150) IENSI1, IENSI2
        IF (IENSI1.EQ.1 .OR. IENSI2.EQ.1) THEN
          WRITE(NFECRA,8155) NBVIS, NVISLA, IVISV1, IVISV2,
     &                       IVISTP, IVISDM, IVISTE, IVISMP
          IF (IPHYLA.EQ.2) THEN
            WRITE(NFECRA,8156) IVISHP, IVISDK, IVISCH, IVISCK
          ENDIF
        ENDIF
C
        WRITE(NFECRA,8160) IENSI3
        IF (IENSI3.EQ.1) THEN
          WRITE(NFECRA,8165) SEUILF, NSTBOR,
     &          INBRBD, IFLMBD, IANGBD, IVITBD, IENCBD, NUSBOR
        ENDIF
C
        WRITE(NFECRA,9900)
C
      ENDIF
C
#if defined(_CS_LANG_FR)
C
 8100 FORMAT(
     &                                                                /,
     &' ** ECOULEMENT DIPHASIQUE LAGRANGIEN'                         ,/,
     &'    --------------------------------'                         ,/,
     &' --- Numero de la phase continue : ',I10                      ,/,
     &'       IILAGR = ',4X,I10,    ' (0 : Lagrangien desactive     ',/,
     &'                ',14X,       '  1 : one way coupling         ',/,
     &'                ',14X,       '  2 : two way coupling         ',/,
     &'                ',14X,       '  3 : sur champs figes        )',/,
     &'       ISUILA = ',4X,I10,    ' (0 : pas de suite ; 1 : suite)',/,
     &'       ISUIST = ',4X,I10,    ' (1 : suite de calcul stats et' ,/,
     &'                ',14X,       '      TS de couplage retour   )',/,
     &' --- Physique particuliere associee aux particules :'         ,/,
     &'       IPHYLA = ',4X,I10,    ' (0 : pas d eqn supplementaires',/,
     &'                ',14X,       '  1 : eqns sur Dp Tp Mp        ',/,
     &'                ',14X,       '  2 : particules de charbon   )'  )
 8105 FORMAT(
     &'       IDPVAR = ',4X,I10,    ' (1 eqn diametre Dp,   0 sinon)',/,
     &'       ITPVAR = ',4X,I10,    ' (1 eqn temperature Tp,0 sinon)',/,
     &'       IMPVAR = ',4X,I10,    ' (1 eqn masse Mp,      0 sinon)'  )
 8106 FORMAT(
     &' --- Parametres Globaux :                                    ',/,
     &'       NBPMAX = ',4X,I10,    ' (nb max de part par iteration)',/,
     &'       NVLS   = ',4X,I10,    ' (nb var particulaires suppl. )',/,
     &'       ISTTIO = ',4X,I10,    ' (1 phase porteuse stationnair)',/,
     &'       INJCON = ',4X,I10,    ' (1 injection continue,0 sinon)',/,
     &'       IROULE = ',4X,I10,    ' (2 clonage/fusion avec calc Y+',/,
     &'                                1 clonage/fusion sans calc Y+',/,
     &'                                0 sinon                     )'  )
C
 8111 FORMAT(
     &' --- Options Charbon :'                                       ,/,
     &'       IENCRA = ',4X,I10,    ' (1 : encrassement si charbon )'  )
C
 8112 FORMAT(
     &'       TPRENC(',I1,') = ', E11.5,
     &                              ' (temp seuil pour encrassement' ,/,
     &'                ',14X,       '  charbon ',I1,8X,'           )'  )
C
 8113 FORMAT(
     &'       VISREF(',I1,') = ', E11.5,
     &                         ' (viscosite critique charbon ',I1,')'  )
C
 8120 FORMAT(
     &' --- Options Couplage Retour :'                               ,/,
     &'       NSTITS = ',4X,I10,    ' (iter de debut moy. en temps )',/,
     &'       LTSDYN = ',4X,I10,    ' (1 couplage retour dynamique )',/,
     &'       LTSMAS = ',4X,I10,    ' (1 couplage retour massique  )',/,
     &'       LTSTHE = ',4X,I10,    ' (1 couplage retour thermique )'  )
C
 8130 FORMAT(
     &' --- Options Statistiques :'                                  ,/,
     &'       ISTALA = ',4X,I10,    ' (1 : calcul de statistiques  )'  )

 8135 FORMAT(
     &'       SEUIL  = ', E14.5,    ' (val min de prise en compte  )',/,
     &'       IDSTNT = ',4X,I10,    ' (iter de debut de calcul stat)',/,
     &'       NSTIST = ',4X,I10,    ' (iter de debut moy. en temps )',/,
     &'       NVLSTS = ',4X,I10,    ' (nb var statistiques suppl.  )'  )

 8140 FORMAT(
     &' --- Options Dispersion Turbulente :                         ',/,
     &'       IDISTU = ',4X,I10,    ' (1 : prise en compte; 0 sinon)',/,
     &'       IDIFFL = ',4X,I10,    ' (1 dispersion =diffusion turb)',/,
     &'       MODCPL = ',4X,I10,    ' (iter lag debut model complet)'  )
 8141 FORMAT(
     &'       IDIRLA = ',4X,I10,    ' (1 2 ou 3 : dir principal ect)'  )
 8142 FORMAT(
     &' --- Options Numeriques :                                    ',/,
     &'       NORDRE = ',4X,I10,    ' (1 ou 2 ordre schema en temps)',/,
     &'       ILAPOI = ',4X,I10,    ' (1 corr. vit instantannees   )'  )
C
 8150 FORMAT(
     &' --- Options Postprocessing Trajectoires/Deplacement :'       ,/,
     &'       IENSI1 = ',4X,I10,    ' (1 : post mode trajectoires  )',/,
     &'       IENSI2 = ',4X,I10,    ' (1 : post mode deplacements  )'  )
C
 8155 FORMAT(
     &'       NBVIS  = ',4X,I10,    ' (nb part max visualisables   )',/,
     &'       NVISLA = ',4X,I10,    ' (periode d acquisition, 0 non)',/,
     &'       IVISV1 = ',4X,I10,    ' (1 : vitesse fluide vu, 0 non)',/,
     &'       IVISV2 = ',4X,I10,    ' (1 : vitesse particule, 0 non)',/,
     &'       IVISTP = ',4X,I10,    ' (1 : temps de sejour  , 0 non)',/,
     &'       IVISDM = ',4X,I10,    ' (1 : diametre part.   , 0 non)',/,
     &'       IVISTE = ',4X,I10,    ' (1 : temperature part., 0 non)',/,
     &'       IVISMP = ',4X,I10,    ' (1 : masse particule  , 0 non)'  )
C
 8156 FORMAT(
     &'       IVISHP = ',4X,I10,    ' (1 : temp/enthal pour charbon)',/,
     &'       IVISDK = ',4X,I10,    ' (1 : diam coeur retrecissant )',/,
     &'       IVISCH = ',4X,I10,    ' (1 : masse de charbon actif  )',/,
     &'       IVISCK = ',4X,I10,    ' (1 : masse de coke           )'  )
C
 8160 FORMAT(
     &' --- Options Stat des Interactions Particules/Frontieres :'   ,/,
     &'       IENSI3 = ',4X,I10,    ' (1 calcul stat parietales    )'  )
C
 8165 FORMAT(
     &'       SEUILF = ', E14.5,    ' (val min de prise en compte  )',/,
     &'       NSTBOR = ',4X,I10,    ' (iter de debut moy. en temps )',/,
     &'       INBRBD = ',4X,I10,    ' (1 : enr. nb d interactions  )',/,
     &'       IFLMBD = ',4X,I10,    ' (1 : enr. flux de masse part.)',/,
     &'       IANGBD = ',4X,I10,    ' (1 : enr. angle d interaction)',/,
     &'       IVITBD = ',4X,I10,    ' (1 : enr. vitesse interaction)',/,
     &'       IENCBD = ',4X,I10,    ' (1 : masse de charbon encrass)',/,
     &'       NUSBOR = ',4X,I10,    ' (1 : enr. infos user suppl.  )'  )
C
#else
C
 8100 FORMAT(
     &                                                                /,
     &' ** TWO-PHASE LANGRANGIEN FLOW'                               ,/,
     &'    --------------------------'                               ,/,
     &' --- Number of the continuous phase: ',I10                    ,/,
     &'       IILAGR = ',4X,I10,    ' (0: Lagrangian deactivated    ',/,
     &'                ',14X,       '  1: one way coupling          ',/,
     &'                ',14X,       '  2: two way coupling          ',/,
     &'                ',14X,       '  3: on frozen fields         )',/,
     &'       ISUILA = ',4X,I10,    ' (0: no restart; 1: restart   )',/,
     &'       ISUIST = ',4X,I10,    ' (1: restart stats and return'  ,/,
     &'                ',14X,       '     coupling ST              )',/,
     &' --- Specific physics associated with particles:'             ,/,
     &'       IPHYLA = ',4X,I10,    ' (0: no additional equations'   ,/,
     &'                ',14X,       '  1: equations on Dp Tp Mp'     ,/,
     &'                ',14X,       '  2: coal particles           )'  )
 8105 FORMAT(
     &'       IDPVAR = ',4X,I10,    ' (1 eqn diameter Dp,      or 0)',/,
     &'       ITPVAR = ',4X,I10,    ' (1 eqn temperature Tp,   or 0)',/,
     &'       IMPVAR = ',4X,I10,    ' (1 eqn mass Mp,          or 0)'  )
 8106 FORMAT(
     &' --- Global parameters:                                      ',/,
     &'       NBPMAX = ',4X,I10,    ' (nb max parts per iteration  )',/,
     &'       NVLS   = ',4X,I10,    ' (nb add. suppl. variables    )',/,
     &'       ISTTIO = ',4X,I10,    ' (1 steady carrier phase      )',/,
     &'       INJCON = ',4X,I10,    ' (1 continuous injection, or 0)',/,
     &'       IROULE = ',4X,I10,    ' (2 clone/merge with Y+ calc'   ,/,
     &'                                1 clone/merge without Y+ calc',/,
     &'                                0 otherwise                 )'  )
C
 8111 FORMAT(
     &' --- Coal options:'''''''''                                   ,/,
     &'       IENCRA = ',4X,I10,    ' (1: fouling if coal          )'  )
C
 8112 FORMAT(
     &'       TPRENC(',I1,') = ', E11.5,
     &                              ' (threshold temp. for coal'     ,/,
     &'                ',14X,       '  fouling ',I1,8X,'           )'  )
C
 8113 FORMAT(
     &'       VISREF(',I1,') = ', E11.5,
     &                         ' (critical coal viscosity    ',I1,')'  )
C
 8120 FORMAT(
     &' --- Return coupling options:'                                ,/,
     &'       NSTITS = ',4X,I10,    ' (start iter for time average )',/,
     &'       LTSDYN = ',4X,I10,    ' (1 dynamic return coupling   )',/,
     &'       LTSMAS = ',4X,I10,    ' (1 mass return coupling      )',/,
     &'       LTSTHE = ',4X,I10,    ' (1 thermal return coupling   )'  )
C
 8130 FORMAT(
     &' --- Statistics options:'                                     ,/,
     &'       ISTALA = ',4X,I10,    ' (1: compute statistics       )'  )

 8135 FORMAT(
     &'       SEUIL  = ', E14.5,    ' (minimum value for handling  )',/,
     &'       IDSTNT = ',4X,I10,    ' (start iter for stat calc    )',/,
     &'       NSTIST = ',4X,I10,    ' (start iter for time avergage)',/,
     &'       NVLSTS = ',4X,I10,    ' (nb add statistics variables )'  )

 8140 FORMAT(
     &' --- Turbulent dispersion options:                           ',/,
     &'       IDISTU = ',4X,I10,    ' (1: accounted for; 0 otherws.)',/,
     &'       IDIFFL = ',4X,I10,    ' (1 dispersion = turb diffus. )',/,
     &'       MODCPL = ',4X,I10,    ' (complete model lag start it.)'  )
 8141 FORMAT(
     &'       IDIRLA = ',4X,I10,    ' (1 2 ou 3: main ect dir)'        )
 8142 FORMAT(
     &' --- Numerical options:                                      ',/,
     &'       NORDRE = ',4X,I10,    ' (1 or 2 time scheme order    )',/,
     &'       ILAPOI = ',4X,I10,    ' (1 inst. velcity corr.       )'  )
C
 8150 FORMAT(
     &' --- Trajectory/displacement postprocessing options:'         ,/,
     &'       IENSI1 = ',4X,I10,    ' (1: post trajectories mode   )',/,
     &'       IENSI2 = ',4X,I10,    ' (1: post displcaments mode   )'  )
C
 8155 FORMAT(
     &'       NBVIS  = ',4X,I10,    ' (max nb visualizatlbe parts  )',/,
     &'       NVISLA = ',4X,I10,    ' (acquisition period, 0 none  )',/,
     &'       IVISV1 = ',4X,I10,    ' (1: fluid velocity vu, 0 none)' ,/,
     &'       IVISV2 = ',4X,I10,    ' (1: particle velocity, 0 none)',/,
     &'       IVISTP = ',4X,I10,    ' (1: resident time,     0 none)',/,
     &'       IVISDM = ',4X,I10,    ' (1: particle diameter, 0 none)',/,
     &'       IVISTE = ',4X,I10,    ' (1: part. temperature, 0 none)',/,
     &'       IVISMP = ',4X,I10,    ' (1: particle mass,     0 none)'  )
C
 8156 FORMAT(
     &'       IVISHP = ',4X,I10,    ' (1: temp/enthalpy for coal   )',/,
     &'       IVISDK = ',4X,I10,    ' (1: shrinking core diameter  )',/,
     &'       IVISCH = ',4X,I10,    ' (1: active coal mass         )',/,
     &'       IVISCK = ',4X,I10,    ' (1: coke mass                )'  )
C
 8160 FORMAT(
     &' --- Statistics options for particles/boundary interaction:'  ,/,
     &'       IENSI3 = ',4X,I10,    ' (1 calculate wall stats      )'  )
C
 8165 FORMAT(
     &'       SEUILF = ', E14.5,    ' (minimul value for handlin   )',/,
     &'       NSTBOR = ',4X,I10,    ' (start iter for time average )',/,
     &'       INBRBD = ',4X,I10,    ' (1: nb interactions rec.     )',/,
     &'       IFLMBD = ',4X,I10,    ' (1: particle mass flow rec.  )',/,
     &'       IANGBD = ',4X,I10,    ' (1: interaction angle rec.   )',/,
     &'       IVITBD = ',4X,I10,    ' (1: interaction velocity rec.)',/,
     &'       IENCBD = ',4X,I10,    ' (1: fouling coal mass        )',/,
     &'       NUSBOR = ',4X,I10,    ' (1: additional user info rec.)'  )
C
#endif
C
C=======================================================================
C 12. METHODE ALE
C=======================================================================
C --- Activation de la methode ALE
C
      WRITE(NFECRA,8210)
      WRITE(NFECRA,8220) IALE, NALINF
C
      WRITE(NFECRA,9900)
C
C
#if defined(_CS_LANG_FR)
C
 8210 FORMAT(
     &                                                                /,
     &' ** METHODE ALE (MAILLAGE MOBILE)'                            ,/,
     &'    -----------                                              ',/)
 8220 FORMAT(
     &'       IALE   = ',4X,I10,    ' (1 : activee                 )',/
     &'       NALINF = ',4X,I10,    ' (Iterations d''initialisation ',/
     &'                                                   du fluide)',/)
C
#else
C
 8210 FORMAT(
     &                                                                /,
     &' ** ALE METHOD (MOVING MESH)'                                 ,/,
     &'    -----------                                              ',/)
 8220 FORMAT(
     &'       IALE   = ',4X,I10,    ' (1: activated                )',/
     &'       NALINF = ',4X,I10,    ' (Fluid initialization'         ,/
     &'                                                  iterations)',/)
C
#endif
C
C=======================================================================
C 13. FIN
C=======================================================================
C
C
      RETURN
      END
c@z
