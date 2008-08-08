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
                        SUBROUTINE LETGEO
C                       *****************
C     -------------------------------------------------------------
     & ( NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR ,
     &   NTETRA , NPYRAM , NPRISM , NHEXAE , INODAL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   ICOTET , ICOPYR , ICOPRI , ICOHEX ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION :
C  --------
c@foncb
CFONC
CFONC  LECTURE DES TABLEAUX ENTITES GEOMETRIQUES
c@fonce
c@argub
CARGU
C-----------------------------------------------------------------------
C                                  ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
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
CARGU ! NTETRA       ! E  !  -> ! NOMBRE DE TETRAEDRES DU MAILLAGE     !
CARGU ! NPYRAM       ! E  !  -> ! NOMBRE DE PYRAMIDES  DU MAILLAGE     !
CARGU ! NPRISM       ! E  !  -> ! NOMBRE DE PRISMES    DU MAILLAGE     !
CARGU ! NHEXAE       ! E  !  -> ! NOMBRE DE HEXAEDRES  DU MAILLAGE     !
CARGU ! INODAL       ! E  ! <-  ! INDIQUE SI L'ON DOIT LIRE LA         !
CARGU !              !    !     !  CONNECTIVITE NODALE POUR LE         !
CARGU !              !    !     !  POST TRAITEMENT                     !
CARGU ! IFACEL       ! TE ! <-  ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE ! <-  ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE ! <-  ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE ! <-  ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE ! <-  ! PROPRIETES D'UNE FAMILLE             !
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! ICOTET       ! TE ! <-  ! CONNECTIVITE TETRAEDRES-NOEUDS       !
CARGU ! ICOPYR       ! TE ! <-  ! CONNECTIVITE PYRAMIDES-NOEUDS        !
CARGU ! ICOPRI       ! TE ! <-  ! CONNECTIVITE PRISMES-NOEUDS          !
CARGU ! ICOHEX       ! TE ! <-  ! CONNECTIVITE HEXAEDRES-NOEUDS        !
CARGU ! IPNFAC       ! TE !  -  ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -  ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -  ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -  ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! XYZCEN       ! TR ! <-  ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET !    !     !                                      !
CARGU ! SURFAC       ! TR ! <-  ! VECTEUR SURFACE DES FACES INTERNES   !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! SURFBO       ! TR ! <-  ! VECTEUR SURFACE DES FACES DE BORD    !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! CDGFAC       ! TR ! <-  ! CENTRE DE GRAVITE DES FACES INTERNES !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! CDGFBO       ! TR ! <-  ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! XYZNOD       ! TR !  -  ! COORDONNES DES NOEUDS                !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
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
      INCLUDE "entsor.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR
      INTEGER          NTETRA , NPYRAM , NPRISM , NHEXAE , INODAL
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          ICOTET(4,NTETRA), ICOPYR(5,NPYRAM)
      INTEGER          ICOPRI(6,NPRISM), ICOHEX(8,NHEXAE)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC)
      DOUBLE PRECISION SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC)
      DOUBLE PRECISION CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD)
C
C VARIABLES LOCALES
C
      INTEGER          KEL, IEL, IDIM, KFACE, IFACE, KFAFBR, KPRFFB
      INTEGER          IOK, IOK1, IOK2, IND
      INTEGER          JDIM,JCEL ,JFAC  ,JFABOR
      INTEGER          JPRFFB, JFAFBR
      INTEGER          JPOINT, JTETRA, JPYRAM, JPRISM, JHEXAE
      INTEGER          IFAC, NN, N1, N2, IP, II1, II2
      INTEGER          NBFAC1, NBFAC2
C
C***********************************************************************
C
      IOK = 0
C
      INODAL = 0
C
C=======================================================================
C 1.  OUVERTURE
C=======================================================================
C
      OPEN (FILE=FICGEO,UNIT=IMPGEO,FORM='formatted')
C
C
C=======================================================================
C 2.  LECTURE DES DIMENSIONS ET VERIFICATIONS
C=======================================================================
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      JDIM = 3
      READ (IMPGEO,1100)  JCEL, JFAC, JFABOR, JPOINT
C
      IF (JDIM.NE.NDIM.OR.JCEL .NE.NCEL ) THEN
        WRITE(NFECRA,2010) '  JDIM',JDIM,'  JCEL',JCEL ,
     &                     '  NDIM',NDIM,'  NCEL',NCEL
        IOK = IOK + 1
      ENDIF
      IF (JFAC  .NE.NFAC  .OR.JFABOR.NE.NFABOR) THEN
        WRITE(NFECRA,2010) '  JFAC',JFAC  ,'JFABOR',JFABOR,
     &                     '  NFAC',NFAC  ,'NFABOR',NFABOR
        IOK = IOK + 1
      ENDIF
      IF (JPOINT.NE.NNOD  ) THEN
        WRITE(NFECRA,2011) 'JPOINT',JPOINT,
     &                     'NNOD  ',NNOD
        IOK = IOK + 1
      ENDIF
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      READ (IMPGEO,1100) JTETRA,JPYRAM,JPRISM,JHEXAE
C
      IF (JTETRA.NE.NTETRA.OR.JPYRAM.NE.NPYRAM) THEN
        WRITE(NFECRA,2010) 'JTETRA',JTETRA,'JPYRAM',JPYRAM,
     &                     'NTETRA',NTETRA,'NPYRAM',NPYRAM
        IOK = IOK + 1
      ENDIF
      IF (JPRISM.NE.NPRISM.OR.JHEXAE.NE.NHEXAE) THEN
        WRITE(NFECRA,2010) 'JPRISM',JPRISM,'JHEXAE',JHEXAE,
     &                     'NPRISM',NPRISM,'NHEXAE',NHEXAE
        IOK = IOK + 1
      ENDIF
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      READ (IMPGEO,1100)  JPRFFB, JFAFBR
C
      IF (JPRFFB.NE.NPRFML.OR.JFAFBR.NE.NFML  ) THEN
        WRITE(NFECRA,2010) 'JPRFFB',JPRFFB,'JFAFBR',JFAFBR,
     &                     'NPRFML',NPRFML,'NFML  ',NFML
        IOK = IOK + 1
      ENDIF
C
      IF (IOK.NE.0) THEN
        WRITE(NFECRA,9999)
        CALL CSEXIT (1)
      ENDIF
C
C
C=======================================================================
C 3.  LECTURE DES TABLEAUX PRINCIPAUX
C=======================================================================
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      DO KFACE = 1, NFAC
        READ (IMPGEO,1100) IFACE,(IFACEL(IEL,KFACE),IEL=1,2)
      ENDDO
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      DO KFACE = 1, NFABOR
        READ (IMPGEO,1100) IFACE,IFABOR(KFACE)
      ENDDO
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      DO KEL = 1, NCEL
        READ (IMPGEO,1120) IEL,XYZCEN(1,KEL),XYZCEN(2,KEL),XYZCEN(3,KEL)
      ENDDO
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      DO KFACE = 1, NFAC
        READ (IMPGEO,1120) IFACE,(SURFAC(IDIM,KFACE),IDIM=1,NDIM)
      ENDDO
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      DO KFACE = 1, NFABOR
        READ (IMPGEO,1120) IFACE,(SURFBO(IDIM,KFACE),IDIM=1,NDIM)
      ENDDO
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      DO KFACE = 1, NFAC
        READ (IMPGEO,1120) IFACE,(CDGFAC(IDIM,KFACE),IDIM=1,NDIM)
      ENDDO
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      DO KFACE = 1, NFABOR
        READ (IMPGEO,1120) IFACE,(CDGFBO(IDIM,KFACE),IDIM=1,NDIM)
      ENDDO
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      DO KFACE = 1, NFABOR
        READ (IMPGEO,1100) IFACE,IFMFBR(KFACE)
      ENDDO
C
      READ (IMPGEO,   *)
      READ (IMPGEO,   *)
      DO KFAFBR = 1, NFML
        READ (IMPGEO,1100)
     &       JFAFBR,(IPRFML(JFAFBR,KPRFFB),KPRFFB=1,NPRFML)
      ENDDO
C
C     IFMCEL NON LU ICI, ON INITIALISE A LA PLACE
C
      DO KEL = 1, NCEL
        IFMCEL(KEL) = 0
      ENDDO
C
C=======================================================================
C 4.  LECTURE DES POINTS SUPPORT ET DE LA CONNECTIVITE NODALE
C=======================================================================
C
C     Selon si l'on dispose des connectivites faces->sommets ou non,
C     on pourra reconstruire une connectivite nodale pour le post
C     traitement. Si celui-ci est demande et que cette connectivite
C     n'est pas disponible, on utilise la connectivite nodale
C     fournie dans le fichier solcom.
C
      IF (LNDFAC.EQ.0 .AND. LNDFBR.EQ.0 .AND. ICHRVL.GT.0) THEN
        INODAL = 1
      ELSE
        INODAL = 0
      ENDIF
C
C     Coordonnees des sommets
C     -----------------------
C
C     Selon si l'on dispose des connectivites faces->sommets ou non,
C     le tableau XYZNOD passe en argument appartient a la connectivite
C     complete ou au maillage nodal destine au post traitement.
C
      READ(IMPGEO,*)
      READ(IMPGEO,*)
C
      IF ((LNDFAC.GT.0 .OR. LNDFBR.GT.0) .OR. INODAL.EQ.1) THEN
C
        DO N1 = 1,NNOD
          READ(IMPGEO,1120) NN,(XYZNOD(IP,NN),IP=1,3)
        ENDDO
C
      ELSE
C
        DO N1 = 1,NNOD
          READ(IMPGEO,*)
        ENDDO
C
      ENDIF
C
C     Connectivites cellules - sommets
C     --------------------------------
C
C     On n'utilise ces connectivites que lorsque l'on ne dispose pas
C     des connectivites faces -> sommets, et que l'on ne pourra donc
C     pas reconstruire une connectivite nodale pour le post traitement.
C
      READ (IMPGEO,*)
      READ (IMPGEO,*)
C
      IF (INODAL .EQ. 1) THEN
C
        DO N1 = 1, NTETRA
          READ(IMPGEO,1100) NN,(ICOTET(IP,N1),IP=1,4)
        ENDDO
        DO N1 = 1, NPYRAM
          READ(IMPGEO,1100) NN,(ICOPYR(IP,N1),IP=1,5)
        ENDDO
        DO N1 = 1, NPRISM
          READ(IMPGEO,1100) NN,(ICOPRI(IP,N1),IP=1,6)
        ENDDO
        DO N1 = 1, NHEXAE
          READ(IMPGEO,1100) NN,(ICOHEX(IP,N1),IP=1,8)
        ENDDO
C
      ELSE
C
        DO N1 = 1, NTETRA
          READ(IMPGEO,*)
        ENDDO
        DO N1 = 1, NPYRAM
          READ(IMPGEO,*)
        ENDDO
        DO N1 = 1, NPRISM
          READ(IMPGEO,*)
        ENDDO
        DO N1 = 1, NHEXAE
          READ(IMPGEO,*)
        ENDDO
C
      ENDIF
C
C=======================================================================
C 5.  LECTURE DES CONNECTIVITES FACE-SOMMET
C=======================================================================
C
C     On peut avoir besoin de la connectivite faces -> sommets pour
C       les modules rayonnement et Lagrangien. Les tableaux
C       correspondants se trouvent en fin de fichier, après
C       les coordonnees des sommets.
C
C     Si l'on n'a pas detecte la presence de ces tableaux
C       supplementaires, on sort immediatement.
C
      IF (LNDFAC.EQ.0 .AND. LNDFBR.EQ.0) RETURN
C
C     Remarque : ces rubriques ont deja ete lues dans LEDGEO, pour
C                determiner les dimensions. On ne prend donc pas de
C                precaution quant a leur existence.
C
C Connectivite faces de bord - sommets
      II2 = 0
      READ(IMPGEO,*,ERR=9000)
      READ(IMPGEO,*,ERR=9000)
      DO N1 = 1,NFABOR
        READ(IMPGEO,'(20I10)',ERR=9000)
     &       NN, IPNFBR(NN), (NODFBR(II1),II1=II2+1,II2+IPNFBR(NN))
        II2 = II2 + IPNFBR(NN)
      ENDDO
C
C Connectivite faces internes - sommets
      II2 = 0
      READ(IMPGEO,*,ERR=9000)
      READ(IMPGEO,*,ERR=9000)
      DO N1 = 1,NFAC
        READ(IMPGEO,'(20I10)',ERR=9000)
     &       NN, IPNFAC(NN), (NODFAC(II1),II1=II2+1,II2+IPNFAC(NN))
        II2 = II2 + IPNFAC(NN)
      ENDDO
C
C=======================================================================
C 6. CREATION DES TABLEAUX IPNFAC et IPNFBR
C=======================================================================
C
C --> On calcule les pointeurs des tableaux NODFAC et NODFBR
C
      NBFAC2 = IPNFAC(1)
      IPNFAC(1) = 1
      DO IFAC = 2,NFAC
        NBFAC1 = NBFAC2
        NBFAC2 = IPNFAC(IFAC)
        IPNFAC(IFAC) = IPNFAC(IFAC-1) + NBFAC1
      ENDDO
      IPNFAC(NFAC+1) = IPNFAC(NFAC) + NBFAC2
C
      NBFAC2 = IPNFBR(1)
      IPNFBR(1) = 1
      DO IFAC = 2, NFABOR
        NBFAC1 = NBFAC2
        NBFAC2 = IPNFBR(IFAC)
        IPNFBR(IFAC) = IPNFBR(IFAC-1) + NBFAC1
      ENDDO
      IPNFBR(NFABOR+1) = IPNFBR(NFABOR) + NBFAC2
C
C=======================================================================
C 7.  CONTROLES SIMPLES
C=======================================================================
C
C --> Verifications de IPNFAC, IPNFBR
C
      IOK = 0
      IF (IPNFAC(NFAC+1).NE.LNDFAC+1) THEN
        WRITE(NFECRA,2001)NFAC  ,LNDFAC+1,IPNFAC(NFAC+1)
        IOK = IOK + 1
      ENDIF
      IF (IPNFBR(NFABOR+1).NE.LNDFBR+1) THEN
        WRITE(NFECRA,2002)NFABOR,LNDFBR+1,IPNFBR(NFABOR+1)
        IOK = IOK + 1
      ENDIF
      IF(IOK.NE.0) THEN
        CALL CSEXIT (1)
      ENDIF
C
C
C --> connectivite faces internes - sommets
      IOK1=0
      DO IFAC = 1,NFAC
        IND = 0
        DO N1 = IPNFAC(IFAC),IPNFAC(IFAC+1)-1
          DO N2 = IPNFAC(IFAC),IPNFAC(IFAC+1)-1
            IF (N1.NE.N2 .AND.
     &         ((XYZNOD(1,NODFAC(N1)).EQ.XYZNOD(1,NODFAC(N2))) .AND.
     &          (XYZNOD(2,NODFAC(N1)).EQ.XYZNOD(2,NODFAC(N2))) .AND.
     &          (XYZNOD(3,NODFAC(N1)).EQ.XYZNOD(3,NODFAC(N2))) )) THEN
              IND=1
            ENDIF
          ENDDO
        ENDDO
        IF (IND.NE.0) THEN
          DO N1 = IPNFAC(IFAC),IPNFAC(IFAC+1)-1
            WRITE(NFECRA,3011)NODFAC(N1),
     &                        N1,XYZNOD(1,NODFAC(N1)),
     &                        N1,XYZNOD(2,NODFAC(N1)),
     &                        N1,XYZNOD(3,NODFAC(N1))
          ENDDO
          IOK1=IOK1+1
        ENDIF
      ENDDO
C
C --> connectivite faces de bord - sommets
      IOK2=0
      DO IFAC = 1,NFABOR
        IND = 0
        DO N1 = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
          DO N2 = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
            IF (N1.NE.N2 .AND.
     &         ((XYZNOD(1,NODFBR(N1)).EQ.XYZNOD(1,NODFBR(N2))) .AND.
     &          (XYZNOD(2,NODFBR(N1)).EQ.XYZNOD(2,NODFBR(N2))) .AND.
     &          (XYZNOD(3,NODFBR(N1)).EQ.XYZNOD(3,NODFBR(N2))) )) THEN
              IND=1
            ENDIF
          ENDDO
        ENDDO
        IF (IND.NE.0) THEN
          DO N1 = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
            WRITE(NFECRA,3012)NODFBR(N1),
     &                        N1,XYZNOD(1,NODFBR(N1)),
     &                        N1,XYZNOD(2,NODFBR(N1)),
     &                        N1,XYZNOD(3,NODFBR(N1))
          ENDDO
          IOK2=IOK2+1
        ENDIF
      ENDDO
C
C --> stop si erreur connectivites
C
      IF(IOK1.NE.0) THEN
        WRITE(NFECRA,3001)
      ENDIF
      IF(IOK2.NE.0) THEN
        WRITE(NFECRA,3002)
      ENDIF
      IF(IOK1.NE.0.OR.IOK2.NE.0) THEN
        CALL CSEXIT (1)
      ENDIF
C
C=======================================================================
C 8.  FERMETURE
C=======================================================================
C
      CLOSE(IMPGEO)
C
      RETURN
C
C=======================================================================
C 9. SORTIES SUR ERREUR DE LECTURE
C=======================================================================
C
C ---> Stop si le maillage lu est incomplet ou erreur a la lecture.
C
 9000 WRITE(NFECRA,9100)
      CALL CSEXIT (1)
C
C---> FORMATS
C
 1100 FORMAT(20I10)
 1120 FORMAT(I10, 5E23.15)
C
 2001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA CREATION DES CONNECTIVITES       ',/,
     &'@    *********   (LETGEO).                                   ',/,
     &'@                                                            ',/,
     &'@    Le nombre de faces internes est NFAC   = ',I10           ,/,
     &'@    IPNFAC(NFAC  +1)devrait valoir ',I10                     ,/,
     &'@                           il vaut ',I10                     ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@  Verifier LETGEO.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA CREATION DES CONNECTIVITES       ',/,
     &'@    *********   (LETGEO).                                   ',/,
     &'@                                                            ',/,
     &'@    Le nombre de faces de bord  est NFABOR = ',I10           ,/,
     &'@    IPNFBR(NFABOR+1)devrait valoir ',I10                     ,/,
     &'@                           il vaut ',I10                     ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@  Verifier LETGEO.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER GEOMETRIE     ',/,
     &'@    *********                                               ',/,
     &'@      INCOHERENCES RENCONTREES                              ',/,
     &'@                                                            ',/,
     &'@      Les dimensions lues dans letgeo ne sont pas en accord ',/,
     &'@        avec les dimensions lues dans ledgeo.               ',/,
     &'@      Lecture de  ',A6,' = ',I10   ,', ',A6,' = ',I10        ,/,
     &'@        alors que ',A6,' = ',I10   ,', ',A6,' = ',I10        ,/,
     &'@                                                            ',/,
     &'@      Le calcul ne peut etre execute.                       ',/,
     &'@                                                            ',/,
     &'@      Contacter l''assistance.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER GEOMETRIE     ',/,
     &'@    *********                                               ',/,
     &'@      INCOHERENCES RENCONTREES                              ',/,
     &'@                                                            ',/,
     &'@      Les dimensions lues dans letgeo ne sont pas en accord ',/,
     &'@        avec les dimensions lues dans ledgeo.               ',/,
     &'@      Lecture de  ',A6,' = ',I10                             ,/,
     &'@        alors que ',A6,' = ',I10                             ,/,
     &'@                                                            ',/,
     &'@      Le calcul ne peut etre execute.                       ',/,
     &'@                                                            ',/,
     &'@      Contacter l''assistance.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA CREATION DES CONNECTIVITES       ',/,
     &'@    *********   (LETGEO)                                    ',/,
     &'@                                                            ',/,
     &'@    Le test a echoue sur la connectivite                    ',/,
     &'@         faces internes -> points supports                  ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@  Verifier LETGEO.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA CREATION DES CONNECTIVITES       ',/,
     &'@    *********   (LETGEO)                                    ',/,
     &'@                                                            ',/,
     &'@    Le test a echoue sur la connectivite                    ',/,
     &'@         faces de bord  -> points supports                  ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@  Verifier LETGEO.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3011 FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LETGEO - CONNECTIVITE : NODFAC = ',I10       ,/,
     &'@                   X(',I10,') = ',E23.15                     ,/,
     &'@                   Y(',I10,') = ',E23.15                     ,/,
     &'@                   Z(',I10,') = ',E23.15                       )
 3012 FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LETGEO - CONNECTIVITE : NODFBR = ',I10       ,/,
     &'@                   X(',I10,') = ',E23.15                     ,/,
     &'@                   Y(',I10,') = ',E23.15                     ,/,
     &'@                   Z(',I10,') = ',E23.15                       )
C
 9100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU MAILLAGE              ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    LE MAILLAGE LU (.slc ou .tlc) EST INCOMPLET             ',/,
     &'@                                                            ',/,
     &'@  Le fichier de maillage ne contient pas les informations   ',/,
     &'@    relatives a la connectivite faces-noeuds, ou contient   ',/,
     &'@    des informations incompletes.                           ',/,
     &'@  Ces informations sont necessaires pour effectuer un calcul',/,
     &'@    en utilisant le module Lagrangien ou Rayonnement        ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier de maillage.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER GEOMETRIE     ',/,
     &'@    *********                                               ',/,
     &'@      INCOHERENCES RENCONTREES                              ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,
     &'@                                                            ',/,
     &'@  Se reporter aux impressions precedentes pour plus de      ',/,
     &'@    renseignements.                                         ',/,
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
