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
                        SUBROUTINE LAGFCH
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  , IBORD  , IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , ETTPA  , TEPA   , STATIS , STATIV ,
     &   TAUP   , TLAG   , PIIL   ,
     &   TSUF   , TSUP   , BX     , TSFEXT ,
     &   VAGAUS , GRADPR , GRADVF ,
     &   ROMP   , FEXTLA ,
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
CFONC     CALCUL DES FORCES DLVO
CFONC
CFONC        - FORCES DE VAN DER WAALS
CFONC        - FORCES ELECTROSTATIQUES
CFONC
CFONC     ELLES DOIVENT ETRE CONNUES EN CHAQUE CELLULE
CFONC       ET ETRE HOMOGENES A LA GRAVITE (M/S2)
CFONC
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC           !
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR           !
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! IBORD        ! TE ! <-  ! SI NORDRE=2, CONTIENT LE NUMERO DE LA!
CARGU !   (NBPMAX)   !    !     !   FACE D'INTERACTION PART/FRONTIERE  !
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
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS                !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! STATIS       ! TR !  -> ! CUMUL POUR LES MOYENNES DES          !
CARGU !(NCELET,NVLSTA!    !     !    STATISTIQUES VOLUMIQUES           !
CARGU ! STATIV       ! TR !  -> ! CUMUL POUR LES VARIANCES DES         !
CARGU !(NCELET,      !    !     !    STATISTIQUES VOLUMIQUES           !
CARGU !   NVLSTA-1)  !    !     !                                      !
CARGU ! TAUP(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE DYNAMIQUE      !
CARGU ! TLAG(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE FLUIDE         !
CARGU ! PIIL(NBPMAX,3! TR !  -> ! TERME DANS L'INTEGRATION DES EDS Up  !
CARGU ! TSUP(NBPMAX,3! TR !  -> ! PREDICTION 1ER SOUS-PAS POUR         !
CARGU !              !    !     !   LA VITESSE DES PARTICULES          !
CARGU ! TSUF(NBPMAX,3! TR !  -> ! PREDICTION 1ER SOUS-PAS POUR         !
CARGU !              !    !     !   LA VITESSE DU FLUIDE VU            !
CARGU ! BX(NBPMAX,3,2! TR !  -> ! CARACTERISTIQUES DE LA TURBULENCE    !
CARGU ! TSFEXT(NBPMAX! TR !  -> ! INFOS POUR LE COUPLAGE RETOUR        !
CARGU ! VAGAUS       ! TR !  -> ! VARIABLES ALEATOIRES GAUSSIENNES     !
CARGU !(NBPMAX,NVGAUS!    !     !                                      !
CARGU ! GRADPR(NCEL,3! TR !  -> ! GRADIENT DE PRESSION                 !
CARGU ! GRADVF(NCEL,3! TR !  -> ! GRADIENT DE LA VITESSE DU FLUIDE     !
CARGU ! ROMP         ! TR !  -  ! MASSE VOLUMIQUE DES PARTICULES       !
CARGU ! FEXTLA       ! TR ! <-  ! CHAMP DE FORCES EXTERIEUR            !
CARGU !(NCELET,3)    !    !     !    UTILISATEUR (M/S2)                !
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
      INCLUDE "numvar.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cpincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ITEPA(NBPMAX,NIVEP) , IBORD(NBPMAX)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION STATIS(NCELET,*),STATIV(NCELET,*)
      DOUBLE PRECISION TAUP(NBPMAX) , TLAG(NBPMAX,3)
      DOUBLE PRECISION PIIL(NBPMAX,3) , BX(NBPMAX,3,2)
      DOUBLE PRECISION TSUF(NBPMAX,3) , TSUP(NBPMAX,3)
      DOUBLE PRECISION TSFEXT(NBPMAX)
      DOUBLE PRECISION VAGAUS(NBPMAX,*)
      DOUBLE PRECISION GRADPR(NCELET,3) , GRADVF(NCELET,9)
      DOUBLE PRECISION ROMP(NBPMAX)
      DOUBLE PRECISION FEXTLA(NBPMAX,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA,IFINIA , IFINRA
      INTEGER          IDPPAR, INXPAR, INYPAR, INZPAR
      INTEGER          IP , IEL , IPHAS , MODE
C
      DOUBLE PRECISION VAL , TEMPF , DNORM
      DOUBLE PRECISION DEBYE, AA
C
C***********************************************************************
C
C=======================================================================
C 0.  GESTION MEMOIRE ET INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPHAS  = 1
C
C=======================================================================
C 1. CALCUL DE LA DISTANCE A LA PAROI + NORMAL A LA PAROI
C=======================================================================
C
      IFINIA = IDEBIA
      IDPPAR = IDEBRA
      INXPAR = IDPPAR + NBPART
      INYPAR = INXPAR + NBPART
      INZPAR = INYPAR + NBPART
      IFINRA = INZPAR + NBPART
      CALL RASIZE('LAGFCH',IFINRA)
C     ===========
C
      DO IP = 1,NBPART
        RA(IDPPAR+IP-1) = 0.D0
        RA(INXPAR+IP-1) = 0.D0
        RA(INYPAR+IP-1) = 0.D0
        RA(INZPAR+IP-1) = 0.D0
      ENDDO
C
      CALL USLADP
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , ETTPA  , TEPA   , STATIS ,
     &   TAUP   , TLAG   , PIIL   ,
     &   VAGAUS , GRADPR , GRADVF ,
     &   ROMP   ,
     &   RA(IDPPAR) ,RA(INXPAR) , RA(INYPAR) , RA(INZPAR) ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 2. FORCES DE VAN DER WAALS
C    Pour etre homogene a des m/s2 on divise par la masse
C=======================================================================
C
      DO IP = 1,NBPART
C
C Force = -A/6 dp/2 /D**2
C
        IF ( RA(IDPPAR+IP-1) .GT. DPARMN ) THEN
C
          VAL = (CSTHAM*ETTP(IP,JDP)/2.D0)
     &         /(6.D0*RA(IDPPAR+IP-1)*RA(IDPPAR+IP-1))
C
          DNORM = SQRT( RA(INXPAR+IP-1)*RA(INXPAR+IP-1)
     &                 +RA(INYPAR+IP-1)*RA(INYPAR+IP-1)
     &                 +RA(INZPAR+IP-1)*RA(INZPAR+IP-1) )
C
C Attention la normale est oriente du fluide vers l'exterieur
C
          AA = DNORM*ETTP(IP,JMP)
C
          FEXTLA(IP,1) = FEXTLA(IP,1) + VAL*RA(INXPAR+IP-1) /AA
          FEXTLA(IP,2) = FEXTLA(IP,2) + VAL*RA(INYPAR+IP-1) /AA
          FEXTLA(IP,3) = FEXTLA(IP,3) + VAL*RA(INZPAR+IP-1) /AA
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 3. FORCES ELECTROSTATIQUES
C    Pour etre homogene a des m/s2 on divise par la masse
C=======================================================================
C
      DO IP = 1,NBPART
C
        IEL = ITEPA(IP,JISOR)
C
C Calcul de la temperature du fluide en fonction du type
C d'ecoulement
C
        IF ( RA(IDPPAR+IP-1) .GT. DPARMN ) THEN
C
          IF ( IPPMOD(ICP3PL).GE.0 .OR.
     &         IPPMOD(ICPL3C).GE.0      ) THEN
C
            TEMPF = PROPCE(IEL,IPPROC(ITEMP1)) - TKELVI
C
          ELSE IF ( IPPMOD(ICOD3P).GE.0 .OR.
     &              IPPMOD(ICOEBU).GE.0 .OR.
     &              IPPMOD(IELARC).GE.0 .OR.
     &              IPPMOD(IELJOU).GE.0      ) THEN
C
            TEMPF = PROPCE(IEL,IPPROC(ITEMP)) - TKELVI
C
          ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.-1 ) THEN
            TEMPF = RTP(IEL,ISCA(ISCALT(IPHAS)))
C
          ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.1 ) THEN
            TEMPF = RTP(IEL,ISCA(ISCALT(IPHAS))) - TKELVI
C
          ELSE IF ( ISCSTH(ISCALT(IPHAS)).EQ.2 ) THEN
            MODE = 1
            CALL USTHHT (MODE, RTP(IEL,ISCA(ISCALT(IPHAS))), TEMPF)
C           ===========
          ELSE
            TEMPF = T0(IPHAS)
          ENDIF
C
C FORCE :
C
C Longueur de Debye
C
          IF (FION .NE. 0 .AND. CSTFAR .GT. 0.D0) THEN
            DEBYE  = SQRT( (EPSEAU*EPSVID*RR*TEMPF)
     &              /(2000.D0*CSTFAR*CSTFAR*FION) )
          ELSE
            WRITE(NFECRA,9001) FION,CSTFAR
            CALL CSEXIT(1)
          ENDIF
C
          IF ( DEBYE .GT. 0.D0 ) THEN
            DEBYE = SQRT(DEBYE)
          ELSE
            WRITE(NFECRA,9002) IP,DEBYE,TEMPF,FION,EPSEAU,EPSVID
            CALL CSEXIT(1)
          ENDIF
C
          VAL = -4.D0*PI*EPSEAU*EPSVID*PHI1*PHI2*(ETTP(IP,JDP)/2.D0)
     &            *EXP(-RA(IDPPAR+IP-1)/DEBYE)
     &            /DEBYE
C
          DNORM = SQRT( RA(INXPAR+IP-1)*RA(INXPAR+IP-1)
     &                 +RA(INYPAR+IP-1)*RA(INYPAR+IP-1)
     &                 +RA(INZPAR+IP-1)*RA(INZPAR+IP-1) )
C
C Attention la normale est oriente du fluide vers l'exterieur
C
          FEXTLA(IP,1)= FEXTLA(IP,1)+VAL*RA(INXPAR+IP-1)
     &                 /(DNORM*ETTP(IP,JMP))
          FEXTLA(IP,2)= FEXTLA(IP,2)+VAL*RA(INYPAR+IP-1)
     &                 /(DNORM*ETTP(IP,JMP))
          FEXTLA(IP,3)= FEXTLA(IP,3)+VAL*RA(INZPAR+IP-1)
     &                 /(DNORM*ETTP(IP,JMP))
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 4. FORCES D'ADHESION
C=======================================================================
C
      IF ( DCOUP .GT. 0.D0 ) THEN
        GAMASV = CSTHAM/(24.D0*PI*DCOUP*DCOUP)
      ELSE
        WRITE(NFECRA,9010) DCOUP
        CALL CSEXIT(1)
      ENDIF
C
      DO IP = 1,NBPART
C
C Force = 3*PI*(Dp/2)*Gamma_SV + SIG2*PI*(Dp/2)/Eps0
C
        IF ( RA(IDPPAR+IP-1) .LE. DPARMN ) THEN
C
          VAL = 3.D0*PI*(ETTP(IP,JDP)/2.D0)*GAMASV
     &        + SIGCH*SIGCH*PI*(ETTP(IP,JDP)/2.D0)/EPSVID
C
          DNORM = SQRT( RA(INXPAR+IP-1)*RA(INXPAR+IP-1)
     &                 +RA(INYPAR+IP-1)*RA(INYPAR+IP-1)
     &                 +RA(INZPAR+IP-1)*RA(INZPAR+IP-1) )
C
C Attention la normale est oriente du fluide vers l'exterieur
C
          AA = DNORM*ETTP(IP,JMP)
C
          FEXTLA(IP,1)= FEXTLA(IP,1)+VAL*RA(INXPAR+IP-1) /AA
          FEXTLA(IP,2)= FEXTLA(IP,2)+VAL*RA(INYPAR+IP-1) /AA
          FEXTLA(IP,3)= FEXTLA(IP,3)+VAL*RA(INZPAR+IP-1) /AA
C
        ENDIF
C
      ENDDO
C
C**********************************************************************
C
C--------
C FORMATS
C--------
C
 9001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LA VALEUR DE LA FORCE IONIQUE EST NULLE                 ',/,
     &'@ OU LA VALEUR DE LA CONSTANTE DE FARADET EST NEGATIVE       ',/,
     &'@                                          OU NULLEE         ',/,
     &'@                                                            ',/,
     &'@       FORCE IONIQUE                : ',G15.7                ,/,
     &'@       CSTE DE FARADET              : ',G15.7                ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les valeurs de FION et de CSTFAR                 ',/,
     &'@                              dans la subroutine USLAG1.    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LA VALEUR DE L''EPAISSEUR DE LA DOUBLE COUCHE           ',/,
     &'@    EST NEGATIVE OU NULLE :                                 ',/,
     &'@       NUMERO DE PARTICULE          : ',I10                  ,/,
     &'@       EPAISSEUR                    : ',G15.7                ,/,
     &'@       TEMPERATURE                  : ',G15.7                ,/,
     &'@       FORCE IONIQUE                : ',G15.7                ,/,
     &'@       CSTE DIELECTIQUE DU VIDE     : ',G15.7                ,/,
     &'@       CSTE DIELECTIQUE DE L''EAU   : ',G15.7                ,/,
     &'@       CSTE DE FARADET              : ',G15.7                ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les valeurs des CSTES dans la subroutine USLAG1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LA VALEUR DE DISTANCE DE COUPURE EST                    ',/,
     &'@    EST NEGATIVE OU NULLE               : ',G15.7            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les valeurs des CSTES dans la subroutine USLAG1. ',/,
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
