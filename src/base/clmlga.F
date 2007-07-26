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
                        SUBROUTINE CLMLGA
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ISYM   , IWARMG , NFECRA , IFINIA , IFINRA , NGR    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DA     , XA     , DAG    , XAG    ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC         MULTIGRILLE ALGEBRIQUE :
CFONC CONSTRUCTION D'UNE HIERARCHIE DE MAILLAGES A PARTIR
CFONC DU MAILLAGE FIN, POUR UNE RESOLUTION MULTIGRILLE
CFONC METHODE A C M (ADDITIVE CORRECTIVE MULTIGRID)
CFONC AVEC REGROUPEMENT DES CELLULES AU PLUS DEUX PAR DEUX
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
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYMETRIQUE    !
CARGU ! IWARMG       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! NFECRA       ! E  !  -> ! NUMERO DU FICHIER D'IMPRESSIONS      !
CARGU ! IFINIA/IFINRA! E  !  -> ! PREMIERE CAS LIBRE EN SORTIE         !
CARGU ! NGR          ! E  !  -> ! NOMBRE DE NIVEAUX DE GRILLES         !
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
CARGU ! IFACLG       ! TE !  -  ! MACRO VOISINS DES MACROS FACES INT   !
CARGU ! (2,NFAC)     !    !     !                                      !
CARGU ! IRESPR(NCELET! TE !  -  ! MACRO EL AUQUEL APPARTIENT UN MICROEL!
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
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DAG(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE (MGM)!
CARGU ! XAG(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE (MGM)!
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
      INCLUDE "mltgrd.h"
C
C***********************************************************************
C
C ARGUMENTS
C

      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ISYM   , IWARMG , NFECRA , IFINIA , IFINRA , NGR
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IFACLG(2,NFAC), IRESPR(NCELET)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DA (NCELET), XA (NFAC ,2)
      DOUBLE PRECISION DAG(NCELET), XAG(NFAC ,2)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          NCELG , NFACG , NCELF , NFACF
      INTEGER          IGR , IEL , IFAC
      INTEGER          JIVOIS , JIP , JIW , JINDIF , JIRESF , JICELF
      INTEGER          JICELC , JRW
C
C***********************************************************************
C
C=======================================================================
C 0.  INITIALISATION
C=======================================================================
C
      IF(IWARMG.GE.1) THEN
        WRITE(NFECRA,1000)
      ENDIF
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C=======================================================================
C 1.  CONSTRUCTION DES MAILLAGES PAR BOUCLE DE RECURSION
C=======================================================================
C
      NCELG = NCEL
      IGR   = 0
C
 100  CONTINUE
C
      IF (IGR.GE.NGRMMX) THEN
        WRITE(NFECRA,2000) NGRMMX
      ENDIF
C
      IF ((NCELG.GT.NCEGRM).AND.(IGR.LT.NGRMAX)) THEN
C
        IGR = IGR + 1
        IF (IWARMG.GE.2) THEN
          WRITE(NFECRA,1020) IGR
        ENDIF
C
        IF (IGR.EQ.1) THEN
C
C       RESERVATION DE TABLEAUX DE TRAVAIL POUR LA CREATION DES
C       MAILLAGES GROSSIERS
C
          CALL MEMGRG
C         ===========
     &    ( IDEBIA , IDEBRA ,
     &      NCEL   , NFAC   , ISYM   ,
     &      JIVOIS , JIP    , JIW    , JINDIF , JIRESF ,
     &      JICELF , JICELC , JRW    ,
     &      IFINIA , IFINRA )
C
C         CREATION DU PREMIER MAILLAGE GROSSIER A PARTIR
C         DU MAILLAGE FIN DE BASE
C
          CALL CALGRG
C         ===========
     &    ( IFINIA , IFINRA ,
     &      IDEBIA , IDEBRA ,
     &      NFECRA ,
     &      ISYM   , IGR    , NTGREN , NTGRRE ,
     &      JNCEL  , JNFAC  , JIFACE , JIRESP , JDA    , JXA    ,
     &      NCEL   , NFAC   ,
     &      IFACEL ,
     &      DA     , XA     ,
     &      NCELG  , NFACG  ,
     &      IPGREN(1,IGR) , IPGRRE(1,IGR) ,
     &      IRESPR , IFACLG ,
     &      DAG    , XAG    ,
     &      IA(JIVOIS) , IA(JIP) , IA(JIW) ,
     &      IA(JINDIF) , IA(JIRESF) , IA(JICELF) , IA(JICELC) ,
     &      RA(JRW)  )
C
        ELSE
C
C       RESERVATION DE TABLEAUX DE TRAVAIL POUR LA CREATION DES
C       MAILLAGES GROSSIERS
C
          CALL MEMGRG
C         ===========
     &    ( IDEBIA , IDEBRA ,
     &      IA(IPGREN(JNCEL,IGR-1))  , IA(IPGREN(JNFAC,IGR-1))  ,
     &      ISYM   ,
     &      JIVOIS , JIP    , JIW    , JINDIF , JIRESF ,
     &      JICELF , JICELC , JRW    ,
     &      IFINIA , IFINRA )
C
C         CREATION DU MAILLAGE GROSSIER DE NIVEAU IGR A PARTIR
C         DU MAILLAGE GROSSIER DE NIVEAU IGR-1
C
          CALL CALGRG
C         ===========
     &    ( IFINIA , IFINRA ,
     &      IDEBIA , IDEBRA ,
     &      NFECRA ,
     &      ISYM   , IGR    , NTGREN , NTGRRE ,
     &      JNCEL  , JNFAC  , JIFACE , JIRESP , JDA    , JXA    ,
     &      IA(IPGREN(JNCEL,IGR-1)) , IA(IPGREN(JNFAC,IGR-1)) ,
     &      IA(IPGREN(JIFACE,IGR-1)) ,
     &      RA(IPGRRE(JDA,IGR-1)) , RA(IPGRRE(JXA,IGR-1)) ,
     &      NCELG  , NFACG  ,
     &      IPGREN(1,IGR) , IPGRRE(1,IGR) ,
     &      IRESPR , IFACLG ,
     &      DAG    , XAG    ,
     &      IA(JIVOIS) , IA(JIP) , IA(JIW) ,
     &      IA(JINDIF) , IA(JIRESF) , IA(JICELF) , IA(JICELC) ,
     &      RA(JRW) )
C
        ENDIF
C
C       STOCKAGE DES DESCRIPTEURS DU MAILLAGE GROSSIER VENANT D'ETRE
C        CONSTRUIT DANS IA ET RA
C
C       VERIFICATION DE L'ESPACE DISPONIBLE
C
        IFINIA = IPGREN(NTGREN,IGR)
        IFINRA = IPGRRE(NTGRRE,IGR)
C
        CALL IASIZE('CLMLGA',IFINIA)
C       ===========
C
        CALL RASIZE('CLMLGA',IFINRA)
C       ===========
C
        IA(IPGREN(JNCEL,IGR)) = NCELG
        IA(IPGREN(JNFAC,IGR)) = NFACG
        IF(IGR.EQ.1) THEN
          NCELF = NCEL
          NFACF = NFAC
        ELSE
          NCELF = IA(IPGREN(JNCEL,IGR-1))
          NFACF = IA(IPGREN(JNFAC,IGR-1))
        ENDIF
        DO IEL = 1, NCELF
          IA(IPGREN(JIRESP,IGR)+IEL-1) = IRESPR(IEL)
        ENDDO
        DO IEL = 1, NCELG
          RA(IPGRRE(JDA,IGR)   +IEL-1) = DAG(IEL)
        ENDDO
        DO IFAC = 1, NFACG
          IA(IPGREN(JIFACE,IGR)+2*(IFAC-1)  ) = IFACLG(1,IFAC)
          IA(IPGREN(JIFACE,IGR)+2*(IFAC-1)+1) = IFACLG(2,IFAC)
        ENDDO
        DO IFAC = 1, NFACG
          RA(IPGRRE(JXA,IGR)+IFAC-1) = XAG(IFAC,1)
        ENDDO
        IF (ISYM.EQ.2) THEN
          DO IFAC = 1, NFACG
            RA(IPGRRE(JXA,IGR)+NFACG+IFAC-1) = XAG(NFACF+IFAC,1)
          ENDDO
        ENDIF
C
        IDEBIA = IFINIA
        IDEBRA = IFINRA
C
        IF (IWARMG.GE.2) THEN
          WRITE(NFECRA,1030) NCELG
          WRITE(NFECRA,1040) NFACG
        ENDIF
C
C ON TESTE : SI ON A TROP PEU AGGLOMERE DE CELLULES, ON S'ARRETE LA
C            SINON, ON CONTINUE
        IF (NCELG.GT.INT(0.8*NCELF)) THEN
          GOTO 200
        ELSE
          GOTO 100
        ENDIF
C
      ENDIF
C
 200  CONTINUE
C
      NGR = IGR
C
      IF (IWARMG.GE.1) THEN
        WRITE(NFECRA,1010) NGR , NCELG
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(/,' CONSTRUCTION DES MAILLAGES GROSSIERS      ',/)
 1010 FORMAT('NOMBRE DE MAILLAGES GROSSIERS                  : ',
     & I8,/, 'NOMBRE DE CELLULES DANS MAILLAGE LE + GROSSIER : ',
     & I8)
 1020 FORMAT('CONSTRUCTION DU MAILLAGE GROSSIER NUMERO: ',I8)
 1030 FORMAT('      NOMBRE DE MAILLES GROSSIERES : ',I8)
 1040 FORMAT('      NOMBRE DE FACES GROSSIERES   : ',I8)
 2000 FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION : clmlga  MULTIGRILLE ALGEBRIQUE              ',/,
     &'@    *********                                               ',/,
     &'@  Nombre de maillages grossiers maximal ',I10   ,' atteint. ',/,
     &'@                                                            '  )
C
C----
C FIN
C----
C
      END
c@z
