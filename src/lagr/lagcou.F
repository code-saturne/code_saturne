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
                        SUBROUTINE LAGCOU
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &   NPRFML , NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  , INDEP  , IBORD  ,
     &   IDEVEL , ITUSER , IA     ,
     &   VOLUME , RTP    , PROPCE ,
     &   ETTP   , ETTPA  , TEPA   , TAUP   , TEMPCT , TSFEXT , TSLAGR ,
     &   CPGD1  , CPGD2  , CPGHT  ,
     &   TSLAG  , VOLP   , VOLM   ,
     &   AUXL1  , AUXL2  , AUXL3  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC       SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC       -----------------------------------
CFONC
CFONC      CALCUL DES TERMES SOURCES DU COUPLAGE RETOUR
CFONC
CFONC      Remarque : les termes sources sont calcules pour
CFONC                 la cellule de depart de la particule
CFONC                 lors de l'iteration courante. Attention, meme
CFONC                 si la particule est sortante du domaine de
CFONC                 calcul (peu importe la maniere) on doit calculer
CFONC                 un terme source qui correspond a ce qu'echange le
CFONC                 fluide porteur et la particule au debut du pas de
CFONC                 temps. Si NORDRE = 2 et que la particule est en
CFONC                 interaction avec la frontiere, alors les termes
CFONC                 source sont calcules comme si NORDRE=1
CFONC                 (on oublie le pre-remplissage de TSFEXT dans
CONFC                 LAGES2).
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
CARGU ! INDEP        ! TE !  -> ! POUR CHAQUE PARTICULE :              !
CARGU !  (NBPMAX)    !    !     !    NUMERO DE LA CELLULE DE DEPART    !
CARGU ! IBORD        ! TE !  -> ! CONTIENT LE NUMERO DE LA             !
CARGU !   (NBPMAX)   !    !     !   FACE D'INTERACTION PART/FRONTIERE  !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! TAUP(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE DYNAMIQUE      !
CARGU ! TSFEXT(NBPMAX! TR !  -> ! FORCES EXTERNES                      !
CARGU ! TEMPCT       ! TR !  -> ! TEMPS CARACTERISTIQUE THERMIQUE      !
CARGU !  (NBPMAX,2)  !    !     !                                      !
CARGU ! TSLAGR(NBPMAX! TR ! <-  ! TERMES SOURCES DE COUPLAGE RETOUR    !
CARGU !     NTERSL)  !    !     !                                      !
CARGU ! CPGD1,CPGD2, ! TR !  -> ! TERMES DE DEVOLATILISATION 1 ET 2 ET !
CARGU !  CPGHT(NBPMAX!    !     !   DE COMBUSION HETEROGENE (CHARBON   !
CARGU !              !    !     !   AVEC COUPLAGE RETOUR THERMIQUE)    !
CARGU ! TSLAG(NBPMAX,! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU !     NTERSL)  !    !     !                                      !
CARGU ! VOLP(NCELET) ! TR !  -  ! FRACTION VOLUMIQUE DES PARTICULES    !
CARGU ! VOLM(NCELET) ! TR !  -  ! FRACTION MASSIQUE DES PARTICULES     !
CARGU ! AUXL1(NBPMAX)! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! AUXL2(NBPMAX)! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! AUXL3(NBPMAX)! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "radiat.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ITEPA(NBPMAX,NIVEP), INDEP(NBPMAX), IBORD(NBPMAX)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION VOLUME(NCELET) , PROPCE(NCELET,*) , RTP(NCELET,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION TAUP(NBPMAX) , TEMPCT(NBPMAX,2)
      DOUBLE PRECISION TSFEXT(NBPMAX)
      DOUBLE PRECISION CPGD1(NBPMAX) , CPGD2(NBPMAX) , CPGHT(NBPMAX)
      DOUBLE PRECISION TSLAG(NCELET,NTERSL)
      DOUBLE PRECISION VOLP(NCELET) , VOLM(NCELET)
      DOUBLE PRECISION TSLAGR(NCELET,NTERSL)
      DOUBLE PRECISION AUXL1(NBPMAX) , AUXL2(NBPMAX) , AUXL3(NBPMAX)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          NPT , IEL , IVAR , ICHA , IPHAS
      DOUBLE PRECISION TVMAX , TAUV , TAUM , AUX1
      DOUBLE PRECISION UUF , VVF , WWF , MF
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IPHAS = ILPHAS
C
      TVMAX = 0.8D0
C
C   Nombre de passage pour les termes sources en stationnaire
C
      IF (ISTTIO.EQ.1 .AND. IPLAS.GE.NSTITS)  NPTS = NPTS + 1
C
      NTXERR = 0
      VMAX   = 0.D0
      TMAMAX = 0.D0
C
      DO IEL=1,NCEL
        VOLP(IEL) = 0.D0
        VOLM(IEL) = 0.D0
      ENDDO
C
      DO IVAR = 1,NTERSL
        DO IEL = 1,NCEL
          TSLAG(IEL,IVAR) = 0.D0
        ENDDO
      ENDDO
C
C=======================================================================
C 2. CALCULS PRELIMINAIRES
C=======================================================================
C
C Finalisation des forces externes (Si la particule a interagit avec
C une frontiere du domaine de calcul, on degenere a l'ordre 1).
C
C
      DO NPT = 1,NBPART
        AUX1 = DTP/TAUP(NPT)
        IF (NORDRE.EQ.1 .OR. IBORD(NPT).GT.0) THEN
          TSFEXT(NPT)= (1.D0-EXP(-AUX1)) *ETTP(NPT,JMP) *TAUP(NPT)
        ELSE
          TSFEXT(NPT) = TSFEXT(NPT)
     &                + (1.D0- (1.D0-EXP(-AUX1)) /AUX1 ) * TAUP(NPT)
     &                * ETTP(NPT,JMP)
        ENDIF
      ENDDO
C
      DO NPT = 1,NBPART
        AUXL1(NPT) = TEPA(NPT,JRPOI)*
     &        ( ETTP(NPT,JMP)  * ETTP(NPT,JUP)
     &         -ETTPA(NPT,JMP) * ETTPA(NPT,JUP)
     &         -GX*TSFEXT(NPT)  ) / DTP
        AUXL2(NPT) = TEPA(NPT,JRPOI)*
     &        ( ETTP(NPT,JMP)  * ETTP(NPT,JVP)
     &         -ETTPA(NPT,JMP)* ETTPA(NPT,JVP)
     &         -GY*TSFEXT(NPT) ) / DTP
        AUXL3(NPT) = TEPA(NPT,JRPOI)*
     &        ( ETTP(NPT,JMP)  * ETTP(NPT,JWP)
     &         -ETTPA(NPT,JMP)* ETTPA(NPT,JWP)
     &         -GZ*TSFEXT(NPT) ) / DTP
      ENDDO
C
C=======================================================================
C 3. TERMES SOURCES DE QUANTITE DE MOUVEMENT
C=======================================================================
C
      IF (LTSDYN.EQ.1) THEN
C
        DO NPT = 1,NBPART
C
          IEL = INDEP(NPT)
C
C Volume et masse des particules dans la maille
C
          VOLP(IEL) = VOLP(IEL)
     &              + TEPA(NPT,JRPOI)*PI*(ETTPA(NPT,JDP)**3)/6.D0
          VOLM(IEL) = VOLM(IEL)
     &              + TEPA(NPT,JRPOI)*ETTPA(NPT,JMP)
C
C TS de QM
C
          TSLAG(IEL,ITSVX) = TSLAG(IEL,ITSVX) - AUXL1(NPT)
          TSLAG(IEL,ITSVY) = TSLAG(IEL,ITSVY) - AUXL2(NPT)
          TSLAG(IEL,ITSVZ) = TSLAG(IEL,ITSVZ) - AUXL3(NPT)
          TSLAG(IEL,ITSLI) = TSLAG(IEL,ITSLI)
     &                     - 2.D0*TEPA(NPT,JRPOI)*ETTP(NPT,JMP)
     &                     / TAUP(NPT)
C
        ENDDO
C
C=======================================================================
C 4. TERMES SOURCES SUR LA TURBULENCE
C=======================================================================
C
        IF (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &       .OR. ITURB(IPHAS).EQ.60 ) THEN
C En v2f (ITURB=50) les TS lagrangiens influent uniquement sur k et eps
C (difficile d'ecrire quoi que ce soit sur v2, qui perd son sens de
C  "composante de Rij")
C
          DO NPT = 1,NBPART
C
            IEL = INDEP(NPT)
C
            UUF = 0.5D0 * ( ETTPA(NPT,JUF) + ETTP(NPT,JUF) )
            VVF = 0.5D0 * ( ETTPA(NPT,JVF) + ETTP(NPT,JVF) )
            WWF = 0.5D0 * ( ETTPA(NPT,JWF) + ETTP(NPT,JWF) )
C
            TSLAG(IEL,ITSKE) = TSLAG(IEL,ITSKE)
     &                       - UUF * AUXL1(NPT)
     &                       - VVF * AUXL2(NPT)
     &                       - WWF * AUXL3(NPT)
C
          ENDDO
C
          DO IEL = 1,NCEL
C
            TSLAG(IEL,ITSKE) = TSLAG(IEL,ITSKE)
     &                       - ETTP(NPT,JUF) * TSLAG(IEL,ITSVX)
     &                       - ETTP(NPT,JVF) * TSLAG(IEL,ITSVY)
     &                       - ETTP(NPT,JWF) * TSLAG(IEL,ITSVZ)
C
          ENDDO
C
        ELSE IF (ITYTUR(IPHAS).EQ.3) THEN
C
          DO NPT = 1,NBPART
C
            IEL = INDEP(NPT)
C
            UUF = 0.5D0 * ( ETTPA(NPT,JUF) + ETTP(NPT,JUF) )
            VVF = 0.5D0 * ( ETTPA(NPT,JVF) + ETTP(NPT,JVF) )
            WWF = 0.5D0 * ( ETTPA(NPT,JWF) + ETTP(NPT,JWF) )
C
            TSLAG(IEL,ITSR11) = TSLAG(IEL,ITSR11)
     &                        - 2.D0 * UUF * AUXL1(NPT)
C
            TSLAG(IEL,ITSR12) = TSLAG(IEL,ITSR12)
     &                        - UUF * AUXL2(NPT)
     &                        - VVF * AUXL1(NPT)
C
            TSLAG(IEL,ITSR13) = TSLAG(IEL,ITSR13)
     &                        - UUF * AUXL3(NPT)
     &                        - WWF * AUXL1(NPT)
C
            TSLAG(IEL,ITSR22) = TSLAG(IEL,ITSR22)
     &                        - 2.D0 * VVF * AUXL2(NPT)
C
            TSLAG(IEL,ITSR23) = TSLAG(IEL,ITSR23)
     &                        - VVF * AUXL3(NPT)
     &                        - WWF * AUXL2(NPT)
C
            TSLAG(IEL,ITSR33) = TSLAG(IEL,ITSR33)
     &                        - 2.D0 * WWF * AUXL3(NPT)
C
          ENDDO
C
          DO IEL = 1,NCEL
C
            TSLAG(IEL,ITSR11) = TSLAG(IEL,ITSR11)
     &                 - 2.D0 * RTP(IEL,IU(IPHAS)) * TSLAG(IEL,ITSVX)
C
            TSLAG(IEL,ITSR12) = TSLAG(IEL,ITSR12)
     &                        - RTP(IEL,IU(IPHAS)) * TSLAG(IEL,ITSVY)
     &                        - RTP(IEL,IV(IPHAS)) * TSLAG(IEL,ITSVX)
C
            TSLAG(IEL,ITSR13) = TSLAG(IEL,ITSR13)
     &                        - RTP(IEL,IU(IPHAS)) * TSLAG(IEL,ITSVZ)
     &                        - RTP(IEL,IW(IPHAS)) * TSLAG(IEL,ITSVX)
C
            TSLAG(IEL,ITSR22) = TSLAG(IEL,ITSR22)
     &                 - 2.D0 * RTP(IEL,IV(IPHAS)) * TSLAG(IEL,ITSVY)
C
            TSLAG(IEL,ITSR23) = TSLAG(IEL,ITSR23)
     &                        - RTP(IEL,IV(IPHAS)) * TSLAG(IEL,ITSVZ)
     &                        - RTP(IEL,IW(IPHAS)) * TSLAG(IEL,ITSVY)
C
            TSLAG(IEL,ITSR33) = TSLAG(IEL,ITSR33)
     &                 - 2.D0 * RTP(IEL,IW(IPHAS)) * TSLAG(IEL,ITSVZ)
C
          ENDDO
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C 5. TERME SOURCE MASSIQUES
C=======================================================================
C
      IF ( LTSMAS.EQ.1 .AND. (IMPVAR.EQ.1 .OR. IDPVAR.EQ.1) ) THEN
C
        DO NPT = 1,NBPART
C
C Dans saturne TSmasse > 0 ===> Apport de masse sur le fluide
C
          IEL = INDEP(NPT)
C
          TSLAG(IEL,ITSMAS) = TSLAG(IEL,ITSMAS) - TEPA(NPT,JRPOI)
     &     * ( ETTP(NPT,JMP) - ETTPA(NPT,JMP) ) /DTP
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C 6. TERMES SOURCES THERMIQUE
C=======================================================================
C
      IF (LTSTHE.EQ.1) THEN
C
        IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) THEN
C
          DO NPT = 1,NBPART
C
            IEL = INDEP(NPT)
C
            TSLAG(IEL,ITSTE) = TSLAG(IEL,ITSTE)
     &     -( ETTP(NPT,JMP)  *ETTP(NPT,JTP) *ETTP(NPT,JCP)
     &        -ETTPA(NPT,JMP) *ETTPA(NPT,JTP)
     &         *ETTPA(NPT,JCP) ) / DTP * TEPA(NPT,JRPOI)
C
            TSLAG(IEL,ITSTI) = TSLAG(IEL,ITSTI)
     &                       + TEMPCT(NPT,2) * TEPA(NPT,JRPOI)
C
          ENDDO
C
          IF (IIRAYO.EQ.1) THEN
C
            DO NPT = 1,NBPART
C
              IEL = INDEP(NPT)
C
              AUX1 = PI *ETTP(NPT,JDP) *ETTP(NPT,JDP) *TEPA(NPT,JREPS)
     &                *(PROPCE(IEL,IPPROC(ILUMN))
     &                -4.D0 *STEPHN *ETTP(NPT,JTP)**4 )
C
              TSLAG(IEL,ITSTE) =TSLAG(IEL,ITSTE)+AUX1*TEPA(NPT,JRPOI)
C
            ENDDO
C
          ENDIF
C
        ELSE IF (IPHYLA.EQ.2) THEN
C
          DO NPT = 1,NBPART
C
            IEL = INDEP(NPT)
            ICHA = ITEPA(NPT,JINCH)
C
            TSLAG(IEL,ITSTE) = TSLAG(IEL,ITSTE)
     &                -( ETTP(NPT,JMP)  *ETTP(NPT,JHP)
     &                    *ETTP(NPT,JCP)
     &                  -ETTPA(NPT,JMP)*ETTPA(NPT,JHP)
     &                    *ETTPA(NPT,JCP) )
     &                 /DTP*TEPA(NPT,JRPOI)
C
            TSLAG(IEL,ITSTI) = TSLAG(IEL,ITSTI)
     &                       + TEMPCT(NPT,2) * TEPA(NPT,JRPOI)
C
            TSLAG(IEL,ITSMV1(ICHA)) = TSLAG(IEL,ITSMV1(ICHA))
     &                              + TEPA(NPT,JRPOI) * CPGD1(NPT)
C
            TSLAG(IEL,ITSMV2(ICHA)) = TSLAG(IEL,ITSMV2(ICHA))
     &                              + TEPA(NPT,JRPOI) * CPGD2(NPT)
C
            TSLAG(IEL,ITSCO)  = TSLAG(IEL,ITSCO)
     &                        + TEPA(NPT,JRPOI) * CPGHT(NPT)
C
            TSLAG(IEL,ITSFP4) = 0.D0
C
          ENDDO
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C 7. Verif que le taux volumique maximal TVMAX admissible de particules
C    ne soit pas depasse dans quelques cellules.
C=======================================================================
C
      DO IEL = 1,NCEL
C
        MF   = VOLUME(IEL) * PROPCE(IEL,IPPROC(IROM(IPHAS)))
        TAUV = VOLP(IEL) / VOLUME(IEL)
        TAUM = VOLM(IEL) / MF
C
        IF (TAUV.GT.TVMAX) THEN
C
          NTXERR = NTXERR + 1
C
          DO IVAR =1,NTERSL
            TSLAGR(IEL,IVAR) = 0.D0
          ENDDO
C
        ENDIF
C
        VMAX   = MAX(TAUV,VMAX)
        TMAMAX = MAX(TMAMAX,TAUM)
C
      ENDDO
C
C=======================================================================
C 8. MOYENNE TEMPORELLE DES TERMES SOURCES
C=======================================================================
C
      IF (ISTTIO.EQ.1 .AND. NPTS.GT.0) THEN
C
        DO IVAR = 1,NTERSL
          DO IEL = 1,NCEL
            TSLAGR(IEL,IVAR) =
     &   ( TSLAG(IEL,IVAR) + (NPTS-1.D0)*TSLAGR(IEL,IVAR) ) / DBLE(NPTS)
          ENDDO
        ENDDO
C
      ELSE
C
        DO IVAR = 1,NTERSL
          DO IEL = 1,NCEL
            TSLAGR(IEL,IVAR) = TSLAG(IEL,IVAR)
          ENDDO
        ENDDO
C
      ENDIF
C
C=======================================================================
C
C----
C FIN
C----
C
      END
c@z
