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
                        SUBROUTINE NUMVEC
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NFAC   , NFABOR , NNOD   , LNDFAC , LNDFBR ,
     &   IRVECI , IRVECB ,
     &   IFACEL , IFABOR , IFMFBR , IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   INUMFI , INUMFB ,
     &   IWORKF , ISMBS  , ISMBV  ,
     &   IPNFAW , NODFAW , IPNFBW , NODFBW ,
     &   RWORKF , RSMBS  , RSMBV )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC CALCUL D UNE TABLE DE RENUMEROTATION DES FACES INTERNES ET DE BORD
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                 ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC  /NFABOR! E  !  -> ! NOMBRE TOTAL DE FACES INTERNES/DE BRD!
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! IRVECI       ! E  ! <-  ! INDICATEUR VECTORISATION FACE INTERN !
CARGU ! IRVECB       ! E  ! <-  ! INDICATEUR VECTORISATION FACE BORD   !
CARGU ! IFACEL       ! TE ! <-> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE ! <-> ! No DE L'ELT VOISIN D'UNE FACE DE BORD!
CARGU ! NFABOR  )    !    !     !                                      !
CARGU ! IFMFBR       ! TE ! <-> ! NUMERO DE LA FAMILLE DES FACES DE BRD!
CARGU ! (NFABOR  )   !    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFAC+1)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (LNDFAC)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFABOR+1) !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR)   !    !     !  (OPTIONNEL)                         !
CARGU ! INUMFI(NFAC) ! TE !  -  ! TABLE DE RENUM DES FACES INTERNES    !
CARGU ! INUMFB(NFABOR! TE !  -  ! TABLE DE RENUM DES FACES DE BORD     !
CARGU ! IWORKF(*     ! TE !  -  ! TAB DE TRAV DE DIM MAX(NFAC,NFABOR)  !
CARGU ! ISMBS (NCELET! TE !  -  ! TAB DE TRAV POUR ASSEMBLAGE SCALAIRE !
CARGU ! ISMBV (NCELET! TE !  -  ! TAB DE TRAV POUR ASSEMBLAGE VECTORIEL!
CARGU ! IPNFAW       ! TE !  -  ! TAB DE TRAV POUR IPNFAC              !
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! NODFAW       ! TE !  -  ! TAB DE TRAV POUR NODFAC              !
CARGU !   (LNDFAC)   !    !     !                                      !
CARGU ! IPNFBW       ! TE !  -  ! TAB DE TRAV POUR IPNFBR              !
CARGU !   (NFABOR+1) !    !     !                                      !
CARGU ! NODFBW       ! TE !  -  ! TAB DE TRAV POUR NODFBR              !
CARGU !   (LNDFBR)   !    !     !                                      !
CARGU ! RWORKF(*     ! TR !  -  ! TAB DE TRAV DE DIM MAX(NFAC,NFABOR)  !
CARGU ! RSMBS (NCELET! TR !  -  ! TAB DE TRAV POUR ASSEMBLAGE SCALAIRE !
CARGU ! RSMBV (NCELET! TR !  -  ! TAB DE TRAV POUR ASSEMBLAGE VECTORIEL!
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
      INCLUDE "vector.h"
      INCLUDE "entsor.h"
C
C ARGUMENTS
C
      INTEGER          NCELET, NCEL, NFAC, NFABOR, NNOD , LNDFAC, LNDFBR
      INTEGER          IRVECI, IRVECB
      INTEGER          IFACEL(2,NFAC),IFABOR(NFABOR), IFMFBR(NFABOR)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          INUMFI(NFAC), INUMFB(NFABOR)
      INTEGER          IWORKF(*), ISMBS(NCELET), ISMBV(NCELET)
      INTEGER          IPNFAW(NFAC+1), NODFAW(LNDFAC)
      INTEGER          IPNFBW(NFABOR+1), NODFBW(LNDFBR)
      DOUBLE PRECISION RWORKF(*), RSMBS(NCELET), RSMBV(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IRELII, NREGII, IRELIB, NREGIB
      INTEGER          ILOOP, IMODAV, IREGIP, IREGIC, JREGIC
      INTEGER          ILAST, INEXT
      INTEGER          IIECH
      INTEGER          II, JJ, IOK, IFAC, IVOIS, IFAC1, NFANP1
      INTEGER          IBLOC, IEL, IREG, ILIG, NFAMAX, ITMP, ISTOP
      INTEGER          JNOD, INOD, INOD1, NBNOD
C
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS COMMUNES
C=======================================================================
C
C --- Numerotation
      DO IFAC = 1, NFAC
        INUMFI(IFAC) = IFAC
      ENDDO
      DO IFAC = 1, NFABOR
        INUMFB(IFAC) = IFAC
      ENDDO
C
C
C=======================================================================
C 2. RANGEMENT DES FACES INTERNES
C     (pour le raisonnement, on place le reliquat a la fin)
C=======================================================================
C
C --- Si l'utilisateur a indique IVECTI = 0, il ne souhaite pas
C       vectoriser, sinon, IVECTI a ete initialise a -1
C
      IF(IVECTI.EQ.0) GOTO 400
C
C --- Indicateur de vectorisation possible = 1
      IVECTI = 0
C
C --- Determination du reliquat et du nbre de registres complets
      IRELII = MOD(NFAC,LREGIS)
      NREGII = NFAC/LREGIS
C
C --- Compteur de boucles
      ILOOP = 0
C
C --- Debut de la boucle externe
  100 CONTINUE
      ILOOP = ILOOP + 1
C
C --- IMODAV = 1 si on a echange un element et un element precedent
C                 dans la table INUMFI
      IMODAV = 0
C
C --- Registre precedent
C
      IREGIC = 0
C
C --- On parcourt les elements de INUMFI
C
      DO JJ = 1, NFAC
C
C --- Registre courant et position dans le registre
        IREGIP = IREGIC
        IREGIC = (JJ-1)/LREGIS+1
        JREGIC = MOD(JJ-1,LREGIS)+1
C
C --- On teste entre ILAST, debut du registre,
C       et la position courante
C
C     En prenant le cas le plus penalisant entre
C       reliquat au debut et reliquat a la fin : reliquat au debut :
C
        IF(IREGIC.EQ.1) THEN
          ILAST = 1
        ELSEIF (JREGIC.LE.IRELII) THEN
          ILAST = (IREGIC-2)*LREGIS+IRELII+1
        ELSE
          ILAST = (IREGIC-1)*LREGIS+1
        ENDIF
C
C     Avec reliquat a la fin :
C
C        ILAST = (IREGIC-1)*LREGIS+1
C
C --- On echange a partir de INEXT, debut du registre suivant
C
C     En prenant le cas le plus penalisant entre reliquat au debut et
C       reliquat a la fin : reliquat au debut :
C
        IF ((IREGIC.EQ.NREGII.AND.JREGIC.GT.IRELII).OR.
     &      (IREGIC.EQ.NREGII+1)                     ) THEN
          INEXT = 1
        ELSEIF (JREGIC.GT.IRELII) THEN
          INEXT = IREGIC*LREGIS+IRELII+1
        ELSE
          INEXT = IREGIC*LREGIS+1
        ENDIF
C
C     Sinon, reliquat a la fin :
C
C        IF ((IREGIC.EQ.NREGII.AND.IRELII.EQ.0) .OR.
C     &       IREGIC.EQ.NREGII+1                   ) THEN
C          INEXT = 1
C        ELSE
C          INEXT = IREGIC*LREGIS+1
C        ENDIF
C
        IF(IREGIC.NE.IREGIP) IIECH = INEXT-1
C
C --- Compteur pour ne pas echanger avec tous les elements de INUMFI
C       plus de n fois
        IBLOC = 0
C
C
C --- Test avec tous les elements precedents depuis ILAST
C      IIECH indique avec quel element de INUMFI on echange
C      IMODAV indique qu'on modifie un element deja vu
C      IBLOC indique qu'on a vu tous les elements et qu'il faut
C              melanger (il n'y a pas de solution)
C
 200    CONTINUE
C
        IFAC = INUMFI(JJ)
        DO II = ILAST, JJ-1
C
          IF ( (IFACEL(1,INUMFI(II)).EQ.IFACEL(1,IFAC)).OR.
     &         (IFACEL(2,INUMFI(II)).EQ.IFACEL(2,IFAC))   ) THEN
C
            IIECH           = IIECH+1
C
            IF(IIECH.GT.NFAC) THEN
              IIECH = 1
              IBLOC = IBLOC + 1
            ENDIF
            IF (IIECH.LT.JJ) IMODAV = 1
            IF (IBLOC.GE.2) THEN
              IBLOC = -1
              GOTO 450
            ENDIF
C
            ITMP             = INUMFI(IIECH )
            INUMFI(IIECH )   = INUMFI(JJ    )
            INUMFI(JJ    )   = ITMP
C
            GOTO 200
          ENDIF
        ENDDO
      ENDDO
C
C --- Si on n'a pas touche aux elements precedant le courant,
C       ca a marche
      IF(IMODAV.EQ.0) THEN
        IVECTI = 1
        GOTO 400
      ENDIF
C
C --- On melange s'il n'y a pas de solution ou si on a boucle 10 fois
 450  CONTINUE
      IF (ILOOP.LE.100.AND.(MOD(ILOOP,10).EQ.0.OR.IBLOC.EQ.-1)) THEN
        DO II = 1, (NFAC-4)/2, 2
          JJ = NFAC-II+1
          ITMP             = INUMFI(II    )
          INUMFI(II    )   = INUMFI(JJ    )
          INUMFI(JJ    )   = ITMP
        ENDDO
      ENDIF
C
C
C --- Et on recommence
      IF(ILOOP.LE.100) GOTO 100
C
 400  CONTINUE
C
C=======================================================================
C 3. RANGEMENT DES FACES DE BORD
C=======================================================================
C
C --- Si l'utilisateur a indique IVECTB = 0, il ne souhaite pas
C       vectoriser, sinon, IVECTB a ete initialise a -1
C
      IF(IVECTB.EQ.0) GOTO 900
C
C --- Indicateur de vectorisation possible = 1
      IVECTB = 0
C
C --- Determination du reliquat et du nbre de registres complets
      IRELIB = MOD(NFABOR,LREGIS)
      NREGIB = NFABOR/LREGIS
C
C --- Nombre max de faces de bord
C      Si > NREGIB : il n'y a pas de solution
      DO IEL = 1, NCEL
        ISMBS(IEL) = 0
      ENDDO
      DO IFAC = 1, NFABOR
        II = IFABOR(IFAC)
        ISMBS(II) = ISMBS(II) + 1
      ENDDO
      NFAMAX = 0
      NFANP1 = 0
      DO IEL = 1, NCEL
        NFAMAX = MAX(NFAMAX,ISMBS(IEL))
        IF(ISMBS(IEL).EQ.NREGIB+1) NFANP1 = NFANP1 + 1
      ENDDO
      IF ( NFAMAX.GT.NREGIB+1.OR.
     &    (NFAMAX.EQ.NREGIB+1.AND.NFANP1.GT.IRELIB)) THEN
        GOTO 900
      ENDIF
C
C --- On classe par nombre de faces de bord du voisin decroissant
C          et numero de voisin decroissant
C
      DO IFAC = 1, NFABOR
        IEL = IFABOR(IFAC)
        IFABOR(IFAC) = IEL + NCEL*ISMBS(IEL)
      ENDDO
      CALL ORDITA(NFABOR,IFABOR,IWORKF)
C     ===========
      DO IFAC = 1, NFABOR
        IEL = MOD(IFABOR(IFAC)-1,NCEL)+1
        IFABOR(IFAC) = IFABOR(IFAC) - NCEL*ISMBS(IEL)
      ENDDO
C
C --- On distribue les faces dans les registres
      DO IFAC = 1, NFABOR
        IF(IFAC.LE.IRELIB*(NREGIB+1)) THEN
          IREG=MOD(IFAC-1,NREGIB+1)+1
          ILIG=(IFAC-1)/(NREGIB+1)+1
          II = (IREG-1)*LREGIS+ILIG
        ELSE
          IFAC1 = IFAC-IRELIB*(NREGIB+1)
          IREG=MOD(IFAC1-1,NREGIB)+1
          ILIG=(IFAC1-1)/NREGIB+1+IRELIB
          II = (IREG-1)*LREGIS+ILIG
        ENDIF
        INUMFB(II)=IWORKF(IFAC)
      ENDDO
      IVECTB=1
C
  900 CONTINUE
C
C=======================================================================
C 4. VERIFICATIONS
C=======================================================================
C
C -----> Verif que toutes les faces se retrouvent une et une seule fois
C       dans INUMFB et INUMFI
C
      IF(IVECTI.EQ.1) THEN
C
        CALL ORDITA(NFAC  ,INUMFI,IWORKF)
C       ===========
C
        IOK = 0
        DO II = 1, NFAC
          IF(INUMFI(IWORKF(II)).NE.NFAC-II+1) IOK = IOK - 1
        ENDDO
        IF (IOK.NE.0) THEN
          WRITE(NFECRA,1100)IOK,1101
          IVECTI = 0
        ENDIF
C
      ENDIF
C
      IF(IVECTB.EQ.1) THEN
C
        CALL ORDITA(NFABOR,INUMFB,IWORKF)
C       ===========
C
        IOK = 0
        DO II = 1, NFABOR
          IF(INUMFB(IWORKF(II)).NE.NFABOR-II+1) IOK = IOK - 1
        ENDDO
        IF (IOK.NE.0) THEN
          WRITE(NFECRA,1200)IOK,1201
          IVECTB = 0
        ENDIF
C
      ENDIF
C
C
C -----> Test classique en balayant les faces precedentes
C
C --- Faces internes
C
      IF(IVECTI.EQ.1) THEN
C
        IOK = 0
        DO JJ = 1, NFAC
C
C --- Registre courant et position dans le registre
          IREGIC = (JJ-1)/LREGIS+1
          JREGIC = MOD(JJ-1,LREGIS)+1
C
C --- On teste entre ILAST, debut du registre,
C       et la position courante
C
C     En prenant le cas le plus penalisant entre
C       reliquat au debut et reliquat a la fin : reliquat au debut :
C
          IF(IREGIC.EQ.1) THEN
            ILAST = 1
          ELSEIF (JREGIC.LE.IRELII) THEN
            ILAST = (IREGIC-2)*LREGIS+IRELII+1
          ELSE
            ILAST = (IREGIC-1)*LREGIS+1
          ENDIF
C
C     Avec reliquat a la fin :
C
C          ILAST = (IREGIC-1)*LREGIS+1
C
C
C --- Test avec tous les elements precedents depuis ILAST
C
          DO II = ILAST, JJ-1
            IFAC = INUMFI(JJ)
            IF((IFACEL(1,INUMFI(II)).EQ.IFACEL(1,IFAC)).OR.
     &         (IFACEL(2,INUMFI(II)).EQ.IFACEL(2,IFAC))    )THEN
              IOK = IOK - 1
            ENDIF
          ENDDO
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,1100)IOK,1102
          IVECTI = 0
        ENDIF
C
      ENDIF
C
C
C --- Faces de bord
C
      IF(IVECTB.EQ.1) THEN
C
        IOK = 0
        DO JJ = 1, NFABOR
C
C --- Registre courant et position dans le registre
          IREGIC = (JJ-1)/LREGIS+1
          JREGIC = MOD(JJ-1,LREGIS)+1
C
C --- On teste entre ILAST, debut du registre,
C       et la position courante
C
C     En prenant le cas le plus penalisant entre
C       reliquat au debut et reliquat a la fin : reliquat au debut :
C
          IF(IREGIC.EQ.1) THEN
            ILAST = 1
          ELSEIF (JREGIC.LE.IRELIB) THEN
            ILAST = (IREGIC-2)*LREGIS+IRELIB+1
          ELSE
            ILAST = (IREGIC-1)*LREGIS+1
          ENDIF
C
C     Avec reliquat a la fin :
C
C          ILAST = (IREGIC-1)*LREGIS+1
C
C
C --- Test avec tous les elements precedents depuis ILAST
C
          DO II = ILAST, JJ-1
            IFAC = INUMFB(JJ)
            IF (IFABOR(INUMFB(II)).EQ.IFABOR(IFAC))THEN
              IOK = IOK - 1
            ENDIF
          ENDDO
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,1200)IOK,1202
          IVECTB = 0
        ENDIF
C
      ENDIF
C
C=======================================================================
C 5. ECHANGES
C=======================================================================
C
C -----> Faces internes
C
      IF(IVECTI.EQ.1) THEN
C
C --- IFACEL
        DO IVOIS = 1, 2
          DO IFAC = 1, NFAC
            IFAC1 = INUMFI(IFAC)
            IWORKF(IFAC)       = IFACEL(IVOIS,IFAC1)
          ENDDO
          DO IFAC = 1, NFAC
            IFACEL(IVOIS,IFAC) = IWORKF(IFAC)
          ENDDO
        ENDDO
C
C --- IPNFAC, NODFAC
        IF(LNDFAC.GT.0) THEN
          JNOD = 1
          IPNFAW(1) = JNOD
          DO IFAC = 1, NFAC
            IFAC1 = INUMFI(IFAC)
            INOD  = IPNFAC(IFAC1)
            INOD1 = IPNFAC(IFAC1+1)
            NBNOD = INOD1-INOD
            DO II = 1, NBNOD
              NODFAW(JNOD+II-1) = NODFAC(INOD+II-1)
            ENDDO
            JNOD = JNOD + NBNOD
            IPNFAW(IFAC+1) = JNOD
          ENDDO
          DO II = 1, LNDFAC
            NODFAC(II) = NODFAW(II)
          ENDDO
          DO II = 1, NFAC
            IPNFAC(II) = IPNFAW(II)
          ENDDO
C Test temporaire a virer
          IF(IPNFAC(NFAC+1).NE.IPNFAW(NFAC+1))THEN
            WRITE(NFECRA,*)' stop ds numvec IPNFAC '
            CALL CSEXIT (1)
          ENDIF
        ENDIF
C
      ENDIF
C
C
C -----> Faces de bord
C
      IF(IVECTB.EQ.1) THEN
C
C --- IFABOR
        DO IFAC = 1, NFABOR
          IFAC1 = INUMFB(IFAC)
          IWORKF(IFAC)       = IFABOR(IFAC1)
        ENDDO
        DO IFAC = 1, NFABOR
          IFABOR(IFAC) = IWORKF(IFAC)
        ENDDO
C
C --- IFMFBR
        DO IFAC = 1, NFABOR
          IFAC1 = INUMFB(IFAC)
          IWORKF(IFAC)         = IFMFBR(IFAC1)
        ENDDO
        DO IFAC = 1, NFABOR
          IFMFBR(IFAC)         = IWORKF(IFAC)
        ENDDO
C
C --- IPNFBR, NODFBR
        IF(LNDFBR.GT.0) THEN
          JNOD = 1
          IPNFBW(1) = JNOD
          DO IFAC = 1, NFABOR
            IFAC1 = INUMFB(IFAC)
            INOD  = IPNFBR(IFAC1)
            INOD1 = IPNFBR(IFAC1+1)
            NBNOD = INOD1-INOD
            DO II = 1, NBNOD
              NODFBW(JNOD+II-1) = NODFBR(INOD+II-1)
            ENDDO
            JNOD = JNOD + NBNOD
            IPNFBW(IFAC+1) = JNOD
          ENDDO
          DO II = 1, LNDFBR
            NODFBR(II) = NODFBW(II)
          ENDDO
          DO II = 1, NFABOR
            IPNFBR(II) = IPNFBW(II)
          ENDDO
        ENDIF
C
      ENDIF
C
      IF ((IVECTB.EQ.1).OR.(IVECTI.EQ.1)) THEN
C
C Stockage des tableaux de renumerotation pour suites
C
        CALL SAVNUM(IVECTI, IVECTB, INUMFI, INUMFB)
C       ===========
C
C Mise a jour des numeros d'elements parents pour les maillages post
C (doit etre appele apres SAVNUM)
C
        CALL PSTRNM(IVECTI, IVECTB, INUMFI, INUMFB)
C       ===========
C
      ENDIF
C
C=======================================================================
C 6. VERIF EN LIVE
C=======================================================================
C
C
C -----> Test d'assemblage sur des entiers
C
      ISTOP = 0
C
C --- Faces internes
C
      IF(IVECTI.EQ.1) THEN
C
        DO IFAC = 1, NFAC
          IWORKF(IFAC) = 1
        ENDDO
        DO IEL = 1, NCELET
          ISMBV(IEL) = 0
          ISMBS(IEL) = 0
        ENDDO
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          ISMBV(II) = ISMBV(II) + IWORKF(IFAC)
          ISMBV(JJ) = ISMBV(JJ) + IWORKF(IFAC)
        ENDDO
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          ISMBS(II) = ISMBS(II) + IWORKF(IFAC)
          ISMBS(JJ) = ISMBS(JJ) + IWORKF(IFAC)
        ENDDO
C
        IOK = 0
        DO IEL = 1, NCEL
          IF(ISMBS(IEL).NE.ISMBV(IEL)) THEN
            IOK = IOK - 1
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,3101)IOK
          ISTOP = 1
        ENDIF
C
      ENDIF
C
C --- Faces de bord
C
      IF(IVECTB.EQ.1) THEN
C
        DO IFAC = 1, NFABOR
          IWORKF(IFAC) = 1
        ENDDO
        DO IEL = 1, NCEL
          ISMBV(IEL) = 0
          ISMBS(IEL) = 0
        ENDDO
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          ISMBV(II) = ISMBV(II) + IWORKF(IFAC)
        ENDDO
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          ISMBS(II) = ISMBS(II) + IWORKF(IFAC)
        ENDDO
C
        IOK = 0
        DO IEL = 1, NCEL
          IF(ISMBS(IEL).NE.ISMBV(IEL)) THEN
            IOK = IOK - 1
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,3201)IOK
          ISTOP = 1
        ENDIF
C
      ENDIF
C
C
C -----> Test d'assemblage sur des reels
C
C --- Faces internes
C
      IF(IVECTI.EQ.1) THEN
C
        DO IFAC = 1, NFAC
          RWORKF(IFAC) = 1.D0
        ENDDO
        DO IEL = 1, NCELET
          RSMBV(IEL) = 0.D0
          RSMBS(IEL) = 0.D0
        ENDDO
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          RSMBV(II) = RSMBV(II) + RWORKF(IFAC)
          RSMBV(JJ) = RSMBV(JJ) + RWORKF(IFAC)
        ENDDO
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          RSMBS(II) = RSMBS(II) + RWORKF(IFAC)
          RSMBS(JJ) = RSMBS(JJ) + RWORKF(IFAC)
        ENDDO
C
        IOK = 0
        DO IEL = 1, NCEL
          IF(ABS(RSMBS(IEL)-RSMBV(IEL)).GT.1.D-20) THEN
            IOK = IOK - 1
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,3102)IOK
          ISTOP  = 1
        ENDIF
C
      ENDIF
C
C --- Faces de bord
C
      IF(IVECTB.EQ.1) THEN
C
        DO IFAC = 1, NFABOR
          RWORKF(IFAC) = 1.D0
        ENDDO
        DO IEL = 1, NCEL
          RSMBV(IEL) = 0.D0
          RSMBS(IEL) = 0.D0
        ENDDO
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          RSMBV(II) = RSMBV(II) + RWORKF(IFAC)
        ENDDO
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          RSMBS(II) = RSMBS(II) + RWORKF(IFAC)
        ENDDO
C
        IOK = 0
        DO IEL = 1, NCEL
          IF(ABS(RSMBS(IEL)-RSMBV(IEL)).GT.1.D-20) THEN
            IOK = IOK - 1
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,3202)IOK
          ISTOP  = 1
        ENDIF
C
      ENDIF
C
      IF(ISTOP.NE.0) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
C -----> Impression finale
C
      WRITE(NFECRA,2000)IVECTI,IVECTB
C
      IRVECI = IVECTI
      IRVECB = IVECTB
C
C=======================================================================
C 6. FORMATS
C=======================================================================
C
 1100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RENUMEROTATION DES FACES                    ',/,
     &'@    *********                                               ',/,
     &'@     PROBLEME DE RENUMEROTATION DES FACES INTERNES          ',/,
     &'@     ',I10   ,' OCCURRENCES DE L''EVENEMENT ',I5             ,/,
     &'@                                                            ',/,
     &'@  Le calcul peut etre execute.                              ',/,
     &'@  La vectorisation restera non forcee                       ',/,
     &'@                                                            ',/,
     &'@  Signaler ce message a l''assistance.                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RENUMEROTATION DES FACES                    ',/,
     &'@    *********                                               ',/,
     &'@     PROBLEME DE RENUMEROTATION DES FACES DE BORD           ',/,
     &'@     ',I10   ,' OCCURRENCES DE L''EVENEMENT ',I5             ,/,
     &'@                                                            ',/,
     &'@  Le calcul peut etre execute.                              ',/,
     &'@  La vectorisation restera non forcee                       ',/,
     &'@                                                            ',/,
     &'@  Signaler ce message a l''assistance.                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2000 FORMAT(
     &'                                                             ',/,
     &' ----------------------------------------------------------- ',/,
     &'                                                             ',/,
     &'                                                             ',/,
     & 3X,'** VECTORISATION (1 : POSSIBLE, 0 : IMPOSSIBLE)          ',/,
     & 3X,'   -------------                                       ',/,/,
     & 3X,'          Faces internes        Faces de bord            ',/,
     & 3X,'   IVECTI     ',I10   ,' IVECTB    ',I10                ,/,/,
     &' ----------------------------------------------------------- ',/)
C
 3101 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS numvec                           ',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME D''ASSEMBLAGE D''ENTIERS AUX FACES INTERNES    ',/,
     &'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Verfier le maillage.                                      ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3201 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS numvec                           ',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME D''ASSEMBLAGE D''ENTIERS AUX FACES DE BORD     ',/,
     &'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Verfier le maillage.                                      ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3102 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS numvec                           ',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME D''ASSEMBLAGE DE REELS   AUX FACES INTERNES    ',/,
     &'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Verfier le maillage.                                      ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3202 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS numvec                           ',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME D''ASSEMBLAGE DE REELS   AUX FACES DE BORD     ',/,
     &'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Verfier le maillage.                                      ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS numvec                           ',/,
     &'@    *********                                               ',/,
     &'@     CONFIGURATION IMPREVUE A LA RENUMEROTATION DES FACES   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Verfier le maillage.                                      ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
      END
c@z
