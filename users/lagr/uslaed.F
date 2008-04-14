c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
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
                        SUBROUTINE USLAED
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITEPA  , IBORD  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT , TSVAR  ,
     &   AUXL1  , AUXL2  , AUXL3  , W1     , W2     , W3     ,
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
CFONC      SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)
CFONC
CFONC      INTEGRATION DES EDS POUR LES VARIABLES UTILISATEUR
CFONC        PAR DEFAUT LES VARIABLES SONT CONSTANTES
CFONC
CFONC
CFONC                                          d T       T - PIP
CFONC      LES EDS DOIVENT ETRE DE LA FORME : ----- = - ---------
CFONC                                          d t         Tca
CFONC
CFONC
CFONC      T : IIIIeme VARIABLE UTILISATEUR QUI EST REPEREE POUR LA
CFONC          PARTICULE IP PAR :
CFONC             T = ETTP(IP,JVLS(IIII))
CFONC             T = ETTPA(IP,JVLS(IIII))
CFONC
CFONC      Tca : TEMPS CARACTERISTIQUE DE L'EDS
CFONC               A DONNER DANS LE TABLEAU    AUXL1
CFONC
CFONC      PIP : COEFFICIENT DE L'EDS ("PSEUDO SECOND MEMBRE")
CFONC               A DONNER DANS LE TABLEAU    AUXL2
CFONC
CFONC               SI LE SCHEMA CHOISI EST D'ORDRE 1 (NORDRE=1)
CFONC               ALORS AU 1ER ET SEUL PASSAGE PIP S'EXPRIME
CFONC               EN FONCTION DE DES GRANDEURS DU PAS DE TEMPS
CFONC               PRECEDENT ETTPA
CFONC
CFONC               SI LE SCHEMA CHOISI EST D'ORDRE 2 (NORDRE=2)
CFONC               ALORS AU 1ER PASSAGE (NOR=1) PIP S'EXPRIME
CFONC               EN FONCTION DE DES GRANDEURS DU PAS DE TEMPS
CFONC               PRECEDENT ETTPA, ET AU 2E PASSAGE (NOR=2)
CFONC               PIP S'EXPRIME EN FONCTION DE DES GRANDEURS
CFONC               DU PAS DE TEMPS COURANT ETTP
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
CARGU ! IBORD        ! TE !  -> ! CONTIENT LE NUMERO DE LA             !
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
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR ! <-> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! TAUP(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE DYNAMIQUE      !
CARGU ! TLAG(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE FLUIDE         !
CARGU ! TEMPCT       ! TR !  -> ! TEMPS CARACTERISTIQUE THERMIQUE ET   !
CARGU !  (NBPMAX,2)  !    !     !   TS IMPLICITE DE COUPLAGE RETOUR    !
CARGU ! TSVAR        ! TR ! <-> ! PREDICTION 1ER SOUS-PAS POUR LA      !
CARGU ! (NBPMAX,NVP1)!    !     !   VARIABLE IVAR, UTILISE POUR LA     !
CARGU !              !    !     !   CORRECTION AU 2EME SOUS-PAS        !
CARGU ! AUXL1(NBPMAX)! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! AUXL2(NBPMAX)! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! AUXL3(NBPMAX)! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! W1..W3(NCELET! TR !  -  ! TABLEAUX DE TRAVAIL                  !
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
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "radiat.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
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
      INTEGER          ITEPA(NBPMAX,NIVEP)  , IBORD(NBPMAX)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION TAUP(NBPMAX) , TLAG(NBPMAX,3) , TEMPCT(NBPMAX,2)
      DOUBLE PRECISION TSVAR(NBPMAX,NVP1)
      DOUBLE PRECISION AUXL1(NBPMAX), AUXL2(NBPMAX), AUXL3(NBPMAX)
      DOUBLE PRECISION W1(NCELET) ,  W2(NCELET) ,  W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          NPT , IEL , IIII , IPL , IPHAS
C
C***********************************************************************
C
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C=======================================================================
C
C On ne rentre en effet dans ce sous programme que si on a defini des
C   variables supplementaires dans uslag1 ; il faut donc necessairement
C   dire comment on les resout (non?).
C
C
      IF(1.EQ.1) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS LE MODULE LAGRANGIEN             ',/,
     &'@    *********                                               ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR uslaed DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Des variables supplementaires ont ete declarees dans      ',/,
     &'@    uslag1 (NVLS=)                                          ',/,
     &'@  Le sous-programme uslaed doit etre complete pour preciser ',/,
     &'@    l''equation differentielle stochastique a resoudre.     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C FIXME : TODO ecrire un exemple utilisateur
C
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPHAS = ILPHAS
C
C=======================================================================
C 2. TEMPS CARACTERISTIQUE DE L'EDS COURANTE
C=======================================================================
C
C Boucle sur les variables supplementaires
C
      DO IIII = 1,NVLS
C
C      Numero dans ETTP de la variable traitee
C
        IPL = JVLS(IIII)
C
        DO NPT = 1,NBPART
C
          IF ( ITEPA(NPT,JISOR).GT.0 ) THEN
C
            IEL = ITEPA(NPT,JISOR)
C
C     Temps caracteristique Tca de l'equation differentielle
C     Cet exemple doit etre adapte au cas traite
C
            AUXL1(NPT) = 1.D0
C
C     Prediction au premier sous pas
C     Cet exemple doit etre adapte au cas traite
C
            IF (NOR.EQ.1) THEN
              AUXL2(NPT) = ETTPA(NPT,IPL)
            ELSE
C
C     Correction au deuxieme sous pas
C     Cet exemple doit etre adapte au cas traite
C
              AUXL2(NPT) = ETTP(NPT,IPL)
            ENDIF
C
          ENDIF
        ENDDO
C
C=======================================================================
C 3. INTEGRATION DE LA VARIABLE IPL
C=======================================================================
C
        CALL LAGITG
C       ===========
     &   ( NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &     IPL    ,
     &     ITEPA(1,JISOR)  , IBORD  ,
     &     ETTP   , ETTPA  , AUXL1  , AUXL2  , TSVAR  )
C
      ENDDO
C
C***********************************************************************
C
C----
C FIN
C----
C
      END
