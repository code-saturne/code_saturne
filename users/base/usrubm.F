c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS                  Code_Saturne version 1.3
C                      ------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2007 EDF S.A., France
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
                        SUBROUTINE USRUBM
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , STATIS ,
     &   RDEVEL , RTUSER , RA     )
C     ---------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC POUR LA SORTIE RUBENS
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VARIABLES STAT. LAGRANGIEN !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET !    !     !                                      !
CARGU ! SURFAC       ! TR !  -> ! VECTEUR SURFACE DES FACES INTERNES   !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! SURFBO       ! TR !  -> ! VECTEUR SURFACE DES FACES DE BORD    !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! CDGFAC       ! TR !  -> ! CENTRE DE GRAVITE DES FACES INTERNES !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS (OPTIONNEL)    !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! STATIS       ! TR !  -> ! MOYENNES STATISTIQUES                !
CARGU !(NCELET,NVLSTA!    !     !    (LAGRANGIEN)                      !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
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
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , NVLSTA
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IEL, NUMERO, NUMPDT, IMPRUB, IPHAS
      INTEGER          IDEBIA, IDEBRA
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C
C INITIALISATION
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IMPRUB = IMPUSR(1)
C
C=======================================================================
C 1. OUVERTURE DU FICHIER
C=======================================================================
C
C
      OPEN(FILE='rubmail',UNIT=IMPRUB,FORM='formatted')
C
C
C=======================================================================
C 2. ENTETE : TITRE, NOM DES VARIABLES
C=======================================================================
C
C TITRE
      WRITE(IMPRUB,1001)
      WRITE(IMPRUB,2001)
      WRITE(IMPRUB,1002)
C
 1001 FORMAT('TITRE_DEBUT')
 2001 FORMAT('VORTEX DE TAYLOR GREEN')
 1002 FORMAT('TITRE_FIN')
C
C NOM DES VARIABLES
      WRITE(IMPRUB,1011)
      WRITE(IMPRUB,2011)
      WRITE(IMPRUB,2012)
      WRITE(IMPRUB,2013)
      WRITE(IMPRUB,2014)
      WRITE(IMPRUB,2015)
      WRITE(IMPRUB,2016)
      WRITE(IMPRUB,2017)
      WRITE(IMPRUB,2018)
      WRITE(IMPRUB,1012)
C
 1011 FORMAT('DEFINITION_VARIABLE_DEBUT')
 2011 FORMAT('p               pa              ')
 2012 FORMAT('vitesse_u       m/s             ')
 2013 FORMAT('vitesse_v       m/s             ')
 2014 FORMAT('vitesse_w       m/s             ')
 2015 FORMAT('scalaire_1      co              ')
 2016 FORMAT('scalaire_2      co              ')
 2017 FORMAT('scalaire_3      co              ')
 2018 FORMAT('scalaire_4      co              ')
 1012 FORMAT('DEFINITION_VARIABLE_FIN')
C
C
C=======================================================================
C 3. COORDONNEES, CONTOUR
C=======================================================================
C
C COORDONNEES
      WRITE(IMPRUB,1021)
      DO IEL = 1, NCEL
           IF(ABS(XYZCEN(1,IEL)).LE.1.D-99)XYZCEN(1,IEL)=0.D0
           IF(ABS(XYZCEN(2,IEL)).LE.1.D-99)XYZCEN(2,IEL)=0.D0
           WRITE(IMPRUB,1100) XYZCEN(1,IEL),XYZCEN(2,IEL)
      ENDDO
      WRITE(IMPRUB,1022)
C
 1021 FORMAT('COORDONNEES_DEBUT')
 1022 FORMAT('COORDONNEES_FIN')
 1100 FORMAT(2E12.4)
C
CMOC GABARIT
CMO      WRITE(IMPRUB,1031)
CMO      WRITE(IMPRUB,   *)' 0.0000E-00  0.0000E-00'
CMO      WRITE(IMPRUB,   *)' 550.00E-03  0.0000E-00'
CMO      WRITE(IMPRUB,   *)' 550.00E-03 60.0000E-03'
CMO      WRITE(IMPRUB,   *)'  50.00E-03 60.0000E-03'
CMO      WRITE(IMPRUB,   *)'  50.00E-03 10.0000E-03'
CMO      WRITE(IMPRUB,   *)' 0.0000E-00 10.0000E-00'
CMO      WRITE(IMPRUB,1032)
CMOC
CMO 1031 FORMAT('GABARIT_DEBUT')
CMO 1032 FORMAT('GABARIT_FIN')
C
      NUMPDT = 0
C
C=======================================================================
C 4. DEBUT DU PAS DE TEMPS ET NUMERO DU PAS DE TEMPS
C=======================================================================
C
      WRITE(IMPRUB,3001)
C
      NUMPDT = NUMPDT + 1
      WRITE(IMPRUB,   *) NUMPDT
C
 3001 FORMAT('PAS_DE_TEMPS_DEBUT')
C
C
C=======================================================================
C 5. VALEUR DES VARIABLES AU CHOIX MAIS COHERENT AVEC LE NBRE DE NOMS
C=======================================================================
C
C VALEURS DES VARIABLES
C
      IPHAS = 1
C
C ------- COMPTEUR DU NOMBRE DE VARIABLES
C
      NUMERO = 0
C
C-----------------------------------------------------------------------
C ------ p1
C
      WRITE(IMPRUB,4001)
C
      DO IEL = 1, NCEL
        IF(ABS(RTP(IEL,IPR(IPHAS))).LE.1.D-99)RTP(IEL,IPR(IPHAS))=0.D0
        WRITE(IMPRUB,4100) RTP(IEL,IPR(IPHAS))
      ENDDO
C
      WRITE(IMPRUB,4002)
C
 4001 FORMAT('VARIABLE_DEBUT_1')
 4002 FORMAT('VARIABLE_FIN_1')
C
C-----------------------------------------------------------------------
C ------ vitesse u
C
      WRITE(IMPRUB,4003)
C
      DO IEL = 1, NCEL
        IF(ABS(RTP(IEL,IU(IPHAS))).LE.1.D-99)RTP(IEL,IU(IPHAS))=0.D0
        WRITE(IMPRUB,4100) RTP(IEL,IU(IPHAS))
      ENDDO
C
      WRITE(IMPRUB,4004)
C
 4003 FORMAT('VARIABLE_DEBUT_2')
 4004 FORMAT('VARIABLE_FIN_2')
C
C-----------------------------------------------------------------------
C ------ vitesse v
C
      WRITE(IMPRUB,4005)
C
      DO IEL = 1, NCEL
        IF(ABS(RTP(IEL,IV(IPHAS))  ).LE.1.D-99)RTP(IEL,IV(IPHAS))  =0.D0
        WRITE(IMPRUB,4100) RTP(IEL,IV(IPHAS))
      ENDDO
C
      WRITE(IMPRUB,4006)
C
 4005 FORMAT('VARIABLE_DEBUT_3')
 4006 FORMAT('VARIABLE_FIN_3')
C
C-----------------------------------------------------------------------
C ------ vitesse w
C
      WRITE(IMPRUB,4007)
C
      DO IEL = 1, NCEL
        IF(ABS(RTP(IEL,IW(IPHAS))  ).LE.1.D-99)RTP(IEL,IW(IPHAS))  =0.D0
        WRITE(IMPRUB,4100) RTP(IEL,IW(IPHAS))
      ENDDO
C
      WRITE(IMPRUB,4008)
C
 4007 FORMAT('VARIABLE_DEBUT_4')
 4008 FORMAT('VARIABLE_FIN_4')
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C ------ scalaire 1
C
      WRITE(IMPRUB,4009)
C
      DO IEL = 1, NCEL
        IF(ABS(RTP(IEL,ISCA(1))  ).LE.1.D-99)RTP(IEL,ISCA(1))  =0.D0
        WRITE(IMPRUB,4100) RTP(IEL,ISCA(1))
      ENDDO
C
      WRITE(IMPRUB,4010)
C
 4009 FORMAT('VARIABLE_DEBUT_5')
 4010 FORMAT('VARIABLE_FIN_5')
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C ------ scalaire 2
C
      WRITE(IMPRUB,4011)
C
      DO IEL = 1, NCEL
        IF(ABS(RTP(IEL,ISCA(2))  ).LE.1.D-99)RTP(IEL,ISCA(2))  =0.D0
        WRITE(IMPRUB,4100) RTP(IEL,ISCA(2))
      ENDDO
C
      WRITE(IMPRUB,4012)
C
 4011 FORMAT('VARIABLE_DEBUT_6')
 4012 FORMAT('VARIABLE_FIN_6')
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C ------ scalaire 3
C
      WRITE(IMPRUB,4013)
C
      DO IEL = 1, NCEL
        IF(ABS(RTP(IEL,ISCA(3))  ).LE.1.D-99)RTP(IEL,ISCA(3))  =0.D0
        WRITE(IMPRUB,4100) RTP(IEL,ISCA(3))
      ENDDO
C
      WRITE(IMPRUB,4014)
C
 4013 FORMAT('VARIABLE_DEBUT_7')
 4014 FORMAT('VARIABLE_FIN_7')
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C ------ scalaire 4
C
      WRITE(IMPRUB,4015)
C
      DO IEL = 1, NCEL
        IF(ABS(RTP(IEL,ISCA(4))  ).LE.1.D-99)RTP(IEL,ISCA(4))  =0.D0
        WRITE(IMPRUB,4100) RTP(IEL,ISCA(4))
      ENDDO
C
      WRITE(IMPRUB,4016)
C
 4015 FORMAT('VARIABLE_DEBUT_8')
 4016 FORMAT('VARIABLE_FIN_8')
C
C-----------------------------------------------------------------------
C=======================================================================
C 6. FIN   DU PAS DE TEMPS
C=======================================================================
C
      WRITE(IMPRUB,3002)
 3002 FORMAT('PAS_DE_TEMPS_FIN')
C
 4100 FORMAT(E12.4)
C
C
C=======================================================================
C 7. FERMETURE ET SORTIE
C=======================================================================
C
      CLOSE(IMPRUB)
C
      RETURN
      END
c@z
