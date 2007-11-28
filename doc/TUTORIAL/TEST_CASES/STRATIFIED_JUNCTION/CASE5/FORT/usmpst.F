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
                        SUBROUTINE USMPST
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , IPART  ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA ,
     &   NCELPS , NFACPS , NFBRPS ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IMODIF ,
     &   ITYPPS , IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   LSTCEL , LSTFAC , LSTFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , STATIS ,
     &   TRACEL , TRAFAC , TRAFBR , RDEVEL , RTUSER , RA     )
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR LA MODIFICATION DES LISTES DE CELLULES
CFONC  OU FACES INTERNES ET DE BORD DEFINISSANT UN MAILLAGE DE POST
CFONC  TRAITEMENT EXISTANT ; CETTE ROUTINE EST APPELEE AUX PAS DE
CFONC  TEMPS AUQUEL CE MAILLAGE EST ACTIF, ET UNIQUEMENT POUR LES
CFONC  MAILLAGES POST UTILISATEUR PRINCIPAUX (NON ALIAS), SI TOUS LES
CFONC  "WRITERS" ASSOCIES A CE MAILLAGE OU SES ALIAS PERMETTENT
CFONC  CETTE MODIFICATION
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! IPART        ! E  !  -> ! NUMERO DU MAILLAGE POST              !
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
CARGU ! NCELPS       ! E  ! <-> ! NOMBRE DE CELLULES DU MAILLAGE POST  !
CARGU ! NFACPS       ! E  ! <-> ! NOMBRE DE FACES INTERIEUR POST       !
CARGU ! NFBRPS       ! E  ! <-> ! NOMBRE DE FACES DE BORD POST         !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IMODIF       ! E  ! <-> ! 0 SI MAILLAGE NON MODIFIE PAR CETTE  !
CARGU !              !    !     ! FONCTION, 1 SI MODIFIE               !
CARGU ! ITYPPS(3)    ! TE !  -> ! INDICATEUR DE PRESENCE (0 OU 1) DE   !
CARGU !              !    !     ! CELLULES (1), FACES (2), OU FACES DE !
CARGU !              !    !     ! DE BORD (3) DANS LE MAILLAGE POST    !
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
CARGU ! LSTCEL(NCELPS! TE ! <-> ! LISTE DES CELLULES DU MAILLAGE POST  !
CARGU ! LSTFAC(NFACPS! TE ! <-> ! LISTE DES FACES INTERIEURES POST     !
CARGU ! LSTFBR(NFBRPS! TE ! <-> ! LISTE DES FACES DE BORD POST         !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TE !  -  ! MACRO TABLEAU ENTIER                 !
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
CARGU ! (NCELET)     !    !     !                                      !
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
CARGU ! STATIS       ! TR !  -> ! STATISTIQUES (LAGRANGIEN)            !
CARGU !NCELET,NVLSTA)!    !     !                                      !
CARGU ! TRACEL(*)    ! TR ! <-> ! TAB REEL VALEURS CELLULES POST       !
CARGU ! TRAFAC(*)    ! TR ! <-> ! TAB REEL VALEURS FACES INT. POST     !
CARGU ! TRAFBR(*)    ! TR ! <-> ! TAB REEL VALEURS FACES BORD POST     !
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
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "numvar.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          IPART
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , NVLSTA
      INTEGER          NCELPS , NFACPS , NFBRPS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE, IMODIF
C
      INTEGER          ITYPPS(3)
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          LSTCEL(NCELPS), LSTFAC(NFACPS), LSTFBR(NFBRPS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTPA(NCELET,*), RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION TRACEL(NCELPS*3)
      DOUBLE PRECISION TRAFAC(NFACPS*3), TRAFBR(NFBRPS*3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          II
C
C
C***********************************************************************
C
C     Remarque : le tableau ITYPPS permet de savoir si le maillage post
C                contient a l'origine des cellules, des faces internes,
C                ou des faces de bord (sur l'ensemble des processeurs).
C
C                Ceci permet d'avoir un traitement "generique" qui
C                peut fonctionner pour tous les numeros de maillage,
C                mais si le maillage post est vide a un instant de
C                post traitement donne, on ne saura plus s'il contenait
C                des cellules ou faces. Dans ce cas, il est preferable
C                d'utiliser explicitement le numero du maillage post
C                pour bien determiner s'il doit contenir des cellules
C                ou des faces.
C
C
C=======================================================================
C     1. TRAITEMENT DES MAILLAGES POST A REDEFINIR
C         A RENSEIGNER PAR L'UTILISATEUR aux endroits indiques
C=======================================================================
C
C     Exemple : maillage 2 : cellules ou T<21 degres
C              
C
C
      IF (IPART.EQ.2) THEN
C
        IMODIF = 1
C
        NCELPS = 0
        NFACPS = 0
        NFBRPS = 0
C
C       SI LE MAILLAGE POST CONTIENT DES CELLULES
C       -----------------------------------------

        IF (ITYPPS(1) .EQ. 1) THEN
C
          DO II = 1, NCEL
C
C
            IF (RTP(II,ISCA(1)) .LE. 21.D0) THEN
              NCELPS = NCELPS + 1
              LSTCEL(NCELPS) = II
            ENDIF
C
          ENDDO
C

        ENDIF
C
      ENDIF
C     Fin du test sur le numero de maillage post.
C
C
      RETURN
C
      END
c@z
