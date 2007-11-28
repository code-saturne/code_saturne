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
                        SUBROUTINE  USDPST
C                       ******************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   LSTCEL , LSTFAC , LSTFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR LOCALISER DES CELLULES, DES FACES
CFONC  INTERNES ET/OU DES FACES DE BORD DEFINISSANT UN MAILLAGE DE
CFONC  POST-TRAITEMENT.
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
CARGU ! (NFML,NPRFML)!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFAC+1)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (LNDFAC)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !  (NFABOR+1)  !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR  ) !    !     !  (OPTIONNEL)                         !
CARGU ! LSTCEL       ! TE !  -  ! TABLEAU DE TRAVAIL (LISTE DES        !
CARGU ! (NCELET)     !    !     !  CELLULES D'UN MAILLAGE DE SORTIE)   !
CARGU ! LSTFAC       ! TE !  -  ! TABLEAU DE TRAVAIL (LISTE DES FACES  !
CARGU ! (NFAC)       !    !     !  INTERNES D'UN MAILLAGE DE SORTIE)   !
CARGU ! LSTFBR       ! TE !  -  ! TABLEAU DE TRAVAIL (LISTE DES FACES  !
CARGU ! (NFABOR)     !    !     !  DE BORD D'UN MAILLAGE DE SORTIE)    !
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
CARGU ! (NCELET      !    !     !                                      !
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
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "optcal.h"
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
      INTEGER          NNOD   , LNDFAC , LNDFBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          LSTCEL(NCELET), LSTFAC(NFAC), LSTFBR(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          INDMOD, ICAS, NBCAS, IPART, NBPART, IPREF
      INTEGER          NTCHRL
C
      INTEGER          NLCEL, NLFAC , NLFBR
      INTEGER          IEL, IFAC  , II
      INTEGER          IDEBIA, IDEBRA
      INTEGER          ICOUL , ICOUL1, ICOUL2, IEL1  , IEL2
      CHARACTER*32     NOMCAS, NOMFMT, NOMMAI
      CHARACTER*96     NOMREP, OPTFMT
C
      DOUBLE PRECISION XFAC  , YFAC  , ZFAC
C
C
C***********************************************************************
C
C
C
      NBCAS  = 0
      NBPART = 0
C
C Entiers "pointeurs" sur la premiere case libre de IA et RA
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C     CREATION DES GESTIONNAIRES D'ECRITURE POUR LE POST TRAITEMENT
C         (UN PAR CAS ET PAR FORMAT, A RENSEIGNER PAR L'UTILISATEUR)
C=======================================================================
C
C     NOMBRE DE GESTIONNAIRES (case au sens EnSight, etude au sens MED,
C                              ou racine d'une arborescence CGNS)
C
      NBCAS = 1
C
      DO ICAS = 1, NBCAS
C
C       INITIALISATIONS DIVERSES
C
        DO II = 1, LEN(NOMCAS)
          NOMCAS (II:II) = ' '
        ENDDO
        DO II = 1, LEN(NOMREP)
          NOMREP (II:II) = ' '
        ENDDO
        DO II = 1, LEN(NOMFMT)
          NOMFMT (II:II) = ' '
        ENDDO
        DO II = 1, LEN(OPTFMT)
          OPTFMT (II:II) = ' '
        ENDDO
C
C       DEFINITION UTILISATEUR :
C
C       NOMCAS et NOMREP indiquent respectivement le prefixe du nom
C       des fichiers et le repertoire correspondant.
C       Si NOMREP est de la forme xxxx.ensight ou xxxx.med, le lanceur le
C       rapatriera automatiquement sous le nom XXXX.ENSIGHT.$DATE ou
C       XXXX.MED.$DATE dans le repertoire RESU. Si NOMREP est d'une autre
C       forme, il faudra gerer son rapatriement a la main.
C
C       NOMFMT permet de choisir le format de sortie
C       ("EnSight Gold" uniquement aujourd'hui, "MED_fichier" et "CGNS"
C       a terme, ou "text" a des fins de debug).
C
C       OPTFMT permet de fournir des options specifiques au format de
C       sortie (separees par des virgules) ;
C         Pour EnSight : "text" ou "binary" (defaut),
C         Pour EnSight, MED, ou CGNS :
C                        "discard_polygons" pour supprimer les polygones,
C                        "discard_polyhedra" pour supprimer les polyedres.
C         Pour EnSight  ou MED :
C                        "divide_polygons" pour découper les polygones,
C                        "divide_polyhedra" pour découper les polyedres.
C
C       INDMOD indique si les maillages ecrits seront :
C         0 : fixes,
C         1 : deformables a topologie constante,
C         2 : modifiables (pourront etre completement redefinis en
C             cours de calcul via le sous-programme USMPST).
C        10 : comme INDMOD = 0, avec champ de déplacement
C        11 : comme INDMOD = 1, avec champ de déplacement
C        12 : comme INDMOD = 2, avec champ de déplacement

C       NTCHRL donne la frequence de sortie par defaut associee,
C       (la sortie a un pas de temps donne pouvant etre forcee ou
C       empechee via le sous-programme utilisateur USNPST).
C
        IF (ICAS .EQ. 1) THEN
C
          NOMCAS = 'chr'
          NOMREP = 'TINF21.ensight'
          NOMFMT = 'EnSight Gold'
          OPTFMT = 'binary, discard_polygons'
          INDMOD = 2
          NTCHRL = 5
C
        ENDIF
C
C       DEFINITION EFFECTIVE

        CALL PSTCWR (ICAS  , NOMCAS, NOMREP, NOMFMT, OPTFMT,
C       ===========
     &               INDMOD, NTCHRL)
C
      ENDDO
C
C=======================================================================
C     NOMBRE DE MAILLAGES EXTRAITS POUR POST TRAITEMENT
C         A RENSEIGNER PAR L'UTILISATEUR
C=======================================================================
C
C   NBPART est le nombre de "parts" qui seront generees
C   (au sens EnSight ; les équivalents MED et CGNS sont le maillage
C    et la base respectivement)
C
C   Une "part" peut etre tout volume ou surface que l'on definira par
C   l'identification des cellules ou faces du maillage
C
C
      NBPART = 2
C
C=======================================================================
C     DEBUT DE LA BOUCLE SUR LES PARTS DEFINIES PAR L'UTILISATEUR
C=======================================================================
C
      DO IPART = 1, NBPART
C
C
C=======================================================================
C       INITIALISATIONS DIVERSES
C         PAS D'INTERVENTION UTILISATEUR REQUISE
C=======================================================================
C
        NLCEL = 0
        NLFAC = 0
        NLFBR = 0
        DO IEL = 1, NCELET
          LSTCEL(IEL) = 0
        ENDDO
        DO IFAC = 1, NFAC
          LSTFAC(IFAC) = 0
        ENDDO
        DO IFAC = 1, NFABOR
          LSTFBR(IFAC) = 0
        ENDDO
C
        DO II = 1, LEN(NOMMAI)
          NOMMAI(II:II) = ' '
        ENDDO
C
C=======================================================================
C       REPERAGE DES CELLULES OU FACES INCLUSES DANS LE MAILLAGE
C         A RENSEIGNER PAR L'UTILISATEUR
C=======================================================================
C
C       Ce sous programme est appele avant la definition des
C        conditions aux limites
C
C
C       POUR LA 1ere COUPE (PART 1) : coupe exemple
C
C         Exemple : on selectionne
C                   les faces internes sur le plan median
C
        IF (IPART .EQ. 1) THEN
C
          NOMMAI = 'Coupe'
C
C         Pour les faces internes
C
          DO IFAC = 1, NFAC
C
C           Determination si la face appartient a la coupe
C
            IF (ABS(CDGFAC(2,IFAC)).LT.1.D-4) THEN
              NLFAC = NLFAC+1
              LSTFAC(NLFAC)= IFAC
            ENDIF
C
          ENDDO
C
C
C
C
C       POUR LA 2eme COUPE (PART 2) : cellules a T<21 degres
C
C         Exemple : on initialise a NCELET, le choix
C                   des cellules sera fait dans usmpst.F
C
        ELSE IF (IPART .EQ. 2) THEN
C
          NOMMAI = 'celTinf21'
C
C         Si NLCEL = NCELET, pas la peine de definir LSTCEL
C
          NLCEL = NCELET
C
        ENDIF
C
C=======================================================================
C       CREATION DES STRUCTURES CONSERVANT LES DONNEES DES PARTS
C         PAS D'INTERVENTION UTILISATEUR REQUISE
C=======================================================================
C
        CALL PSTCMA (IPART, NOMMAI,
C       ===========
     &               NLCEL, NLFAC, NLFBR, LSTCEL, LSTFAC, LSTFBR)
C
C
C=======================================================================
C       IDENTIFICATION DU MAILLAGE EXTRAIT ET GESTION DE SORTIE
C         A RENSEIGNER PAR L'UTILISATEUR
C=======================================================================
C
        IF (IPART .EQ. 1) THEN
C
C         Le maillage 1 est associe a la sortie standard (ICAS=-1)
          ICAS = -1
          CALL PSTASS(IPART, ICAS)
C         ===========
C
        ELSE IF (IPART .EQ. 2) THEN
C
C         Le maillage 2 est associe au cas cree ici (ICAS=1)
          ICAS = 1
          CALL PSTASS(IPART, ICAS)
C         ===========
C
        ENDIF
C
C=======================================================================
C     FIN   DE LA BOUCLE SUR LES PARTS DEFINIES PAR L'UTILISATEUR
C=======================================================================
C
      ENDDO
C
      RETURN
C
C=======================================================================
C     FORMATS
C=======================================================================
C
      END
c@z

