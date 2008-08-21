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
                        SUBROUTINE COUPBI
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NFABOR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ICODCL ,
     &   IDEVEL , ITUSER , IA     ,
     &   RCODCL ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC LECTURE DE DONNEES RELATIVES A UN COUPLAGE AVEC SYRTHES
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ICODCL       ! TE ! <-  ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! RCODCL       ! TR ! <-  ! VALEUR DES CONDITIONS AUX LIMITES    !
CARGU !  (NFABOR,NVAR!    !     !  AUX FACES DE BORD                   !
CARGU !              !    !     ! RCODCL(1) = VALEUR DU DIRICHLET      !
CARGU !              !    !     ! RCODCL(2) = VALEUR DU COEF. D'ECHANGE!
CARGU !              !    !     !  EXT. (INFINIE SI PAS D'ECHANGE)     !
CARGU !              !    !     ! RCODCL(3) = VALEUR DE LA DENSITE DE  !
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2 OU      !
CARGU !              !    !     !  HAUTEUR DE RUGOSITE (m) si ICODCL=6 !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0, IDBRA0
      INTEGER          NFABOR
      INTEGER          NVAR , NSCAL , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          LL, NBCCOU, INBCOU, INBCOO, NBFCOU
      INTEGER          IFAC, ILOC, ISCAL , IPHAS
      INTEGER          IDEBIA, IDEBRA, IFINIA, IFINRA, IPFCOU, ITHPAR
      INTEGER          ICLDEF
      INTEGER          MODE
      DOUBLE PRECISION TEMPER, ENTHAL

C
C***********************************************************************
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C     COUPLAGE SYRTHES : RECUPERATION DE LA TEMPERATURE DE PAROI
C=======================================================================
C
C     RECUPERATION DU NOMBRE DE CAS DE COUPLAGE
C
      CALL NBCSYR (NBCCOU)
C     ===========
C
C---> BOUCLE SUR LES COUPLAGES : RECEPTION DES TABLEAUX TPAROI POUR
C                                CHAQUE COUPLAGE ET APPLICATION DE
C                                LA CONDITION LIMITE CORRESPONDANTE
C
      DO INBCOU = 1, NBCCOU
C
C        NOMBRE DE FACES DE BORD PAR CAS DE COUPLAGE
        INBCOO = INBCOU
        CALL NBFSYR (INBCOO, NBFCOU)
C       ===========
C
C        GESTION MEMOIRE POUR RECEVOIR LE TABLEAU
C
        IPFCOU = IDEBIA
        IFINIA = IPFCOU + NBFCOU
C
        ITHPAR = IDEBRA
        IFINRA = ITHPAR + NBFCOU
C
        CALL IASIZE('COUPBI',IFINIA)
C       ===========
        CALL RASIZE('COUPBI',IFINRA)
C       ===========
C
C        LECTURE DU MESSAGE (TEMPERATURE PAROI) ET
C        INTERPOLATION DES VALEURS AUX SOMMETS AUX FACES
        CALL VARSYI (INBCOU, NBFCOU, RA(ITHPAR))
C       ===========
C
C        ON IMPOSE LA TEMPERATURE A LA PAROI
C
        INBCOO = INBCOU
        CALL LFASYR(INBCOO, NBFCOU, IA(IPFCOU))
C       ===========
C
C        TYPE DE CONDITION PAR DEFAUT
        ICLDEF = 5
C
C
        DO ISCAL = 1, NSCAL
C
          IF(ICPSYR(ISCAL).EQ.1) THEN
C
C --- Pour les scalaires qui sont couples a SYRTHES
C       on impose une condition de Dirichlet aux faces couplees
C     Pour le moment, on ne passe ici qu'une seule fois,
C       etant entendu que un seul scalaire est couple a SYRTHES
C     Pour le module compressible, on resout en energie, mais on
C       stocke la temperature a part, pour que ce soit plus clair
C       dans les conditions aux limites

C
            LL = ISCA(ISCAL)
            IF(IPPMOD(ICOMPF).GE.0) THEN
              IPHAS = IPHSCA(ISCAL)
              IF(ISCAL.EQ.IENERG(IPHAS)) THEN
                LL = ISCA(ITEMPK(IPHAS))
              ELSE
                WRITE(NFECRA,9000)IENERG(IPHAS),ISCAL
                CALL CSEXIT (1)
              ENDIF
            ENDIF

C
            DO ILOC = 1, NBFCOU
C
              IFAC = IA(IPFCOU+ILOC-1)
C
              IF ((ICODCL(IFAC,LL) .NE. 1) .AND.
     &            (ICODCL(IFAC,LL) .NE. 5) .AND.
     &            (ICODCL(IFAC,LL) .NE. 6)) ICODCL(IFAC,LL) = ICLDEF
C
              RCODCL(IFAC,LL,1) = RA(ITHPAR+ILOC-1)
              RCODCL(IFAC,LL,2) = RINFIN
              RCODCL(IFAC,LL,3) = 0.D0
C
            ENDDO
C
C           Conversion eventuelle temperature -> enthalpie
C
            IF(ISCSTH(LL).EQ.2) THEN
C
              DO ILOC = 1, NBFCOU
                IFAC = IA(IPFCOU+ILOC-1)
C
                TEMPER = RCODCL(IFAC,LL,1)
                MODE   = -1
                CALL USTHHT(MODE,ENTHAL,TEMPER)
C               ===========
                RCODCL(IFAC,LL,1) = ENTHAL
C
              ENDDO
C
            ENDIF
C
          ENDIF
C
        ENDDO
C
      ENDDO
C
C======================================================================
C     FIN DES COUPLAGES DE BORD
C======================================================================
C
      RETURN
C
C FORMATS
C
#if defined(_CS_LANG_FR)
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU COUPLAGE SYRTHES              ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Avec le module compressible, seul le scalaire ',I10        ,/,
     &'@    peut etre couple a SYRTHES. Ici, on cherche a coupler   ',/,
     &'@    le scalaire ',I10                                        ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE SYRTHES COUPLING                  ',/,
     &'@    ********                                                ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  With the compressible module, only the scalar ',I10        ,/,
     &'@    may be coupled with SYRTHES. Here, one tries to couple  ',/,
     &'@    with the scalar ',I10                                    ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#endif
C
      END
c@z
