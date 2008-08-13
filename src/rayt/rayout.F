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
                        SUBROUTINE RAYOUT
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RAYEXP , RAYIMP ,
     &   TPAROI , QINCID , XLAM   , EPA    , EPS    ,
     &   FLUNET , FLCONV , HFCONV ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE RAYONNEMENT :
CFONC    --------------------------------------
CFONC
CFONC   1) ECRITURE FICHIER SUITE,
CFONC   2) Ecriture des fichiers Ensight pour les sorties sur les
CFONC      frontieres du domaine de calcul
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
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDNOD       ! E  !  -> ! LONGUEUR DU TABLEAU ICOCEL (OPTIONNEL!
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU !              !    !     ! LE MODULE LAGRANGIEN                 !
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
CARGU ! RAYEXP(NCELET! TR ! <-  ! TERME SOURCE RADIATIF EXPLICITE      !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! RAYIMP(NCELET! TR ! <-  ! TERME SOURCE RADIATIF IMPLICITE      !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! TEMPK(NCELET)! TR ! <-  ! TEMPERATURE EN KELVIN                !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! TPAROI(NFABOR! TR !  -  ! TEMPERATURE DE PAROI EN KELVIN       !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! QINCID(NFABOR! TR ! <-  ! DENSITE DE FLUX RADIATIF AUX BORDS   !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! XLAM(NFABOR  ! TR ! <-  ! COEFFICIENT DE CONDUCTIVITE THERMIQUE!
CARGU !   ,NPHAST)   !    !     ! DES FACETTES DE PAROI (W/m/K)        !
CARGU ! EPA (NFABOR  ! TR ! <-  ! EPAISSEUR DES FACETTES DE PAROI (m)  !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! EPS (NFABOR  ! TR ! <-  ! EMISSIVITE DES FACETTES DE BORD      !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! FLUNET(NFABOR! TR ! <-  ! DENSITE DE FLUX NET RADIATIF AUX     !
CARGU !   ,NPHAST)   !    !     ! FACES DE BORD                        !
CARGU ! FLCONV(NFABOR! TR ! <-  ! DENSITE DE FLUX CONVECTIF AUX FACES  !
CARGU !   ,NPHAST)   !    !     ! DE BORD                              !
CARGU ! HFCONV(NFABOR! TR ! <-  ! COEFFICIENT D'ECHANGE FLUIDE AUX     !
CARGU !   ,NPHAST)   !    !     ! FACES DE BORD                        !
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
C-----------------------------------------------------------------------
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
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "radiat.h"
      INCLUDE "parall.h"
C
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          NVAR   , NSCAL  , NPHAS  , IPH
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION RAYEXP(NCELET,NPHAST), RAYIMP(NCELET,NPHAST)
      DOUBLE PRECISION TPAROI(NFABOR,NPHAST), QINCID(NFABOR,NPHAST)
      DOUBLE PRECISION XLAM(NFABOR,NPHAST), EPA(NFABOR,NPHAST)
      DOUBLE PRECISION EPS(NFABOR,NPHAST), FLUNET(NFABOR,NPHAST)
      DOUBLE PRECISION FLCONV(NFABOR,NPHAST), HFCONV(NFABOR,NPHAST)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER        RUBRIQ*64
      CHARACTER        CPHASE(NPHSMX)*2
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IFINIA , IFINRA
      INTEGER          ITRAV1 , IP
      INTEGER          IERROR , NBERRO , IRTYP , ITYSUP , NBVAL
      INTEGER          IVERS  , ILECEC , IPHAS
      INTEGER          IMPAVR
C
C***********************************************************************
C=======================================================================
C 0 - GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. ECRITURE DU FICHIER SUITE DU MODULE DE RAYONNEMENT
C=======================================================================
C
C
C ---> On code en chaine le numero des phases
C
C     Verifications pour les format et les numero de phase
C       et de scalaire en chaine.
      IF(NPHAS .GT.99) THEN
        WRITE(NFECRA,9012)NPHAS
        CALL CSEXIT (1)
      ENDIF
C
C     Codage en chaine de caracteres du numero de la phase
      DO IPHAS = 1, NPHAS
        WRITE(CPHASE(IPHAS),'(I2.2)')IPHAS
      ENDDO
C
C ---> Ouverture (et on saute si erreur)
C     ILECEC = 2 : ecriture
C
      WRITE(NFECRA,6010)
C
      ILECEC = 2
      CALL OPNSUI(FICAVR, LEN(FICAVR), ILECEC, IMPAVR, IERROR)
C     ===========
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,9020)
        GOTO 9998
      ENDIF
C
      WRITE(NFECRA,6011)
C
C Entete et Dimensions ou on saute si erreur
C     On inclut une rubrique destinee a distinguer ce fichier
C       d'un autre fichier suite
C     Pour le moment, IVERS n'est pas utilise
C
      NBERRO = 0
C
      IVERS  = 111
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      RUBRIQ = 'version_fichier_suite_rayonnement'
      CALL ECRSUI(IMPAVR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IVERS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
C
      RUBRIQ = 'nombre_phases_rayt'
      CALL ECRSUI(IMPAVR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NPHAST,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      IF(NBERRO.NE.0) THEN
        WRITE(NFECRA,9120)
        GOTO 9998
      ENDIF
C
      WRITE(NFECRA,6012)
C
C Temps (par securite)
C
      NBERRO = 0
C
      RUBRIQ = 'nbre_pas_de_temps'
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      CALL ECRSUI(IMPAVR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NTCABS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'instant_precedent'
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 2
      CALL ECRSUI(IMPAVR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,TTCABS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      IF(NBERRO.NE.0) THEN
        WRITE(NFECRA,8121)
      ENDIF
C
C Donnees
C
      NBERRO = 0
C
      DO IPHAS = 1, NPHAST
C
C     Aux faces de bord
C
        ITYSUP = 3
        NBVAL  = 1
        IRTYP  = 2
C
        RUBRIQ = 'tparoi_fb_phase'//CPHASE(IPHAS)
        CALL ECRSUI(IMPAVR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              TPAROI(1,IPHAS),IERROR)
        NBERRO=NBERRO+IERROR
C
        RUBRIQ = 'qincid_fb_phase'//CPHASE(IPHAS)
        CALL ECRSUI(IMPAVR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              QINCID(1,IPHAS),IERROR)
        NBERRO=NBERRO+IERROR
C
        RUBRIQ = 'hfconv_fb_phase'//CPHASE(IPHAS)
        CALL ECRSUI(IMPAVR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              HFCONV(1,IPHAS),IERROR)
        NBERRO=NBERRO+IERROR
C
        RUBRIQ = 'flconv_fb_phase'//CPHASE(IPHAS)
        CALL ECRSUI(IMPAVR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              FLCONV(1,IPHAS),IERROR)
        NBERRO=NBERRO+IERROR
C
C
C     Aux cellules
C
        ITYSUP = 1
        NBVAL  = 1
        IRTYP  = 2
C
        RUBRIQ = 'rayimp_ce_phase'//CPHASE(IPHAS)
        CALL ECRSUI(IMPAVR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              RAYIMP(1,IPHAS),IERROR)
        NBERRO=NBERRO+IERROR
C
        RUBRIQ = 'rayexp_ce_phase'//CPHASE(IPHAS)
        CALL ECRSUI(IMPAVR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              RAYEXP(1,IPHAS),IERROR)
        NBERRO=NBERRO+IERROR
C
      ENDDO
C
C  ---> Si pb : on saute
C
      IF(NBERRO.NE.0) THEN
        WRITE(NFECRA,9100)
        GOTO 9998
      ENDIF
C
      WRITE(NFECRA,6013)
C
C ---> Fermeture du fichier suite
      CALL CLSSUI(IMPAVR,IERROR)
C
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,8011) FICAVR
      ENDIF
C
      WRITE(NFECRA,6014)
C
C ---> En cas d'erreur, on continue quand meme
 9998 CONTINUE
C
C=======================================================================
C 2. ECRITURE DU FICHIER BORD DU MODULE DE RAYONNEMENT
C=======================================================================
C
C    Uniquement en monoprocesseur pour l'instant
C
      IF (NRANGP.EQ.1) THEN
C
        IFINIA = IDEBIA
C
        ITRAV1 = IDEBRA
        IFINRA = ITRAV1 + NFABOR
        CALL RASIZE ('RAYOUT',IFINRA)
C       ===========
C
        DO IP = 1,NPHAST
C
          IPH = IP
C
          CALL RAYBRD
C         ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IPH    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IITYPF) , IA(IITRIF)  , IA(IIZFRD) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RA(ISRFBN) ,
     &   TPAROI , QINCID , XLAM   , EPA    , EPS    ,
     &   FLUNET , FLCONV , HFCONV , RA(ITRAV1) ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDDO
C
      ELSE
C
        WRITE(NFECRA,7000)
C
      ENDIF
C
C
      RETURN
C
C
C--------
C FORMATS
C--------
C
 6010 FORMAT (/, 3X,'** INFORMATIONS SUR LE MODULE DE RAYONNEMENT ',/,
     &           3X,'   ------------------------------------------',/,
     &           3X,' Ecriture d''un fichier suite                ',/)
C
 6011 FORMAT (   3X,'   Debut de l''ecriture                      ',/)
 6012 FORMAT (   3X,'   Fin de l''ecriture des dimensions         ',/)
 6013 FORMAT (   3X,'   Fin de l''ecriture des donnees            ',/)
 6014 FORMAT (   3X,' Fin de l''ecriture du fichier suite         ',/)
C
 9020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    ERREUR A L''OUVERTURE DU FICHIER SUITE RAYONNEMENT      ',/,
     &'@                                                            ',/,
     &'@  Le calcul continue mais                                   ',/,
     &'@            ne fournira pas de fichier suite rayonnement.   ',/,
     &'@                                                            ',/,
     &'@  Verifier que le repertoire de travail est accessible en   ',/,
     &'@    ecriture et que le fichier suite peut y etre cree.      ',/,
     &'@  Voir le sous-programme rayout.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9012 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :ARRET A L''ECRITURE DU FICHIER SUITE         ',/,
     &'@    *********                                   RAYONNEMENT ',/,
     &'@                                                            ',/,
     &'@      Le nombre de phases            NPHAS  doit etre       ',/,
     &'@        inferieur ou egal a  99, il vaut ici ',I10           ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme rayout.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/)
C
 9120 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      ERREUR LORS DE L''ECRITURE DES DIMENSIONS             ',/,
     &'@                                                            ',/,
     &'@  Le calcul continue mais                                   ',/,
     &'@            ne fournira pas de fichier suite rayonnement.   ',/,
     &'@                                                            ',/,
     &'@  Voir le sous-programme rayout.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8121 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      ERREUR LORS DE L''ECRITURE DU PAS DE TEMPS ET DU TEMPS',/,
     &'@                                                            ',/,
     &'@    Le calcul continue...                                   ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme rayout.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      ERREUR LORS DE L''ECRITURE DES DONNEES                ',/,
     &'@                                                            ',/,
     &'@  Le calcul continue mais                                   ',/,
     &'@            ne fournira pas de fichier suite rayonnement.   ',/,
     &'@                                                            ',/,
     &'@  Voir le sous-programme rayout.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,
     &'@    *********                              AVAL RAYONNMEMENT',/,
     &'@                                                            ',/,
     &'@    Probleme sur le fichier de nom (',A13,')                ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7000 FORMAT(
     &'                                                             ',/,
     &' Module de rayonnement :                                     ',/,
     &'     Ecriture du fichier Ensight de bord non disponible      ',/,
     &'     en parallele.                                           ',/,
     &'                                                             ',/)
C
C----
C FIN
C----
C
      END
c@z
