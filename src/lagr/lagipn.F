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
                        SUBROUTINE LAGIPN
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCELET , NCEL   ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NPAR1  , NPAR2  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   RTP    ,
     &   ETTP   , TEPA   , VAGAUS ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC      Initialisation de la vitesse fluide vu pour les nouvelles
CFONC      particules.
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
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NPAR1 ,NPAR2 ! E  !  -> ! BORNE MIN ET MAX DES PARTICULES      !
CARGU !              !    !     !    A INITIALISER                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! VAGAUS       ! TR ! <-  ! VARIABLES ALEATOIRES GAUSSIENNES     !
CARGU !(NBPMAX,NVGAUS!    !     !                                      !
CARGU ! W1...W3(NCEL)! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
C
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NCELET , NCEL
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NPAR1 , NPAR2
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ITEPA(NBPMAX,NIVEP)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION RTP(NCELET,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION VAGAUS(NBPMAX,*)
      DOUBLE PRECISION W1(NCELET) ,  W2(NCELET) ,  W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL , NPT , NOMB , IPHAS
      DOUBLE PRECISION TU , D2S3
C
C***********************************************************************
C
C=======================================================================
C 0.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IPHAS = ILPHAS
      D2S3 = 2.D0 / 3.D0
C
C=======================================================================
C 2. SIMULATION DES VITESSES TURBULENTES FLUIDES INSTANTANNEES VUES
C    PAR LES PARTICULES SOLIDES LE LONG DE LEUR TRAJECTOIRE.
C=======================================================================
C
      IF (IDISTU.EQ.1) THEN
C
        IF (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &       .OR. ITURB(IPHAS).EQ.60) THEN
          DO IEL = 1,NCEL
            W1(IEL) = RTP(IEL,IK(IPHAS))
          ENDDO
        ELSE IF (ITYTUR(IPHAS).EQ.3) THEN
          DO IEL = 1,NCEL
            W1(IEL) = 0.5D0 * ( RTP(IEL,IR11(IPHAS))
     &                        + RTP(IEL,IR22(IPHAS))
     &                        + RTP(IEL,IR33(IPHAS)) )
          ENDDO
        ELSE
          WRITE(NFECRA,9000) IILAGR, IDISTU, IPHAS, ITURB(IPHAS)
          CALL CSEXIT (1)
C         ===========
        ENDIF
      ELSE
        DO IEL = 1,NCEL
          W1(IEL) = 0.D0
        ENDDO
      ENDIF
C
C
C---> CALCUL DES TIRAGES ALEATOIRES
C     CALCUL DU TEMPS CARACTERISTIQUE DES PARTICULES
C     remarque : NORMALEN est dans le fichier ZUFALL.F
C     ^^^^^^^^
C
      IF (IDISTU.EQ.1) THEN
        NOMB = NPAR2-NPAR1+1
        CALL NORMALEN (NOMB,VAGAUS(NPAR1,1))
        CALL NORMALEN (NOMB,VAGAUS(NPAR1,2))
        CALL NORMALEN (NOMB,VAGAUS(NPAR1,3))
      ELSE
        DO NPT = NPAR1,NPAR2
          VAGAUS(NPT,1) = 0.D0
          VAGAUS(NPT,2) = 0.D0
          VAGAUS(NPT,3) = 0.D0
        ENDDO
      ENDIF
C
      DO NPT = NPAR1,NPAR2
C
        IEL = ITEPA(NPT,JISOR)
C
        TU = SQRT( D2S3*W1(IEL) )
C
        ETTP(NPT,JUF) = RTP(IEL,IU(IPHAS)) + VAGAUS(NPT,1)*TU
        ETTP(NPT,JVF) = RTP(IEL,IV(IPHAS)) + VAGAUS(NPT,2)*TU
        ETTP(NPT,JWF) = RTP(IEL,IW(IPHAS)) + VAGAUS(NPT,3)*TU
C
      ENDDO
C
C--------
C FORMATS
C--------
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGIPN)                                    ',/,
     &'@                                                            ',/,
     &'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,
     &'@    DE TURBULENCE SELECTIONNE.                              ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@   Le module Lagrangien a ete active avec IILAGR = ',I10     ,/,
     &'@     et la dispersion turbulente est prise en compte        ',/,
     &'@                                     avec IDISTU = ',I10     ,/,
     &'@   Le modele de turbulence active pour la phase ',I6         ,/,
     &'@     correspond a ITURB  = ',I10                             ,/,
     &'@   Or, les seuls traitements de la turbulence compatibles   ',/,
     &'@     avec le module Lagrangien et la dispersion turbulente  ',/,
     &'@     sont k-epsilon et Rij-epsilon, v2f et k-omega.         ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,
     &'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,
     &'@  USINI1.                                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      END
c@z
