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
                        SUBROUTINE RAYBRD
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IPH    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , IZFRAD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   SURFBN ,
     &   TPAROI , QINCID , XLAM   , EPA    , EPS    ,
     &   FLUNET , FLCONV , HFCONV , W1     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
CFONC    -----------------------------------------
CFONC
CFONC  Sortie au format CASE de grandeurs
CFONC  sur les faces de frontieres du domaine de calcul.
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
CARGU ! LNDNOD       ! E  !  -> ! LONGUEUR DU TABLEAU ICOCEL (OPTIONNEL!
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! IPH          ! E  !  -> ! NUMERO DE LA PHASE COURANTE ASSOCIEE !
CARGU !              !    !     ! AU RAYONNEMENT                       !
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
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! ITRIFB(NFABOR! TE ! <-  ! TAB D'INDIRECTION POUR TRI DES FACES !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! IZFRAD(NFABOR! TE ! <-> ! NUMERO DE ZONE DES FACES DE BORD     !
CARGU ! ISOTHM(NFABOR! TE ! <-> ! TYPE DE CONDITION DE PAROI           !
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
CARGU ! SURFBN(NFABOR! TR !  -> ! SURFACE DES FACES DE BORD            !
CARGU ! TPAROI(NFABOR! TR !  -> ! TEMPERATURE DE PAROI EN KELVIN       !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! QINCID(NFABOR! TR !  -> ! DENSITE DE FLUX RADIATIF AUX BORDS   !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! XLAM(NFABOR  ! TR !  -> ! COEFFICIENT DE CONDUCTIVITE THERMIQUE!
CARGU !   ,NPHAST)   !    !     ! DES FACETTES DE PAROI (W/m/K)        !
CARGU ! EPA (NFABOR  ! TR !  -> ! EPAISSEUR DES FACETTES DE PAROI (m)  !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! EPS (NFABOR  ! TR !  -> ! EMISSIVITE DES FACETTES DE BORD      !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! FLUNET(NFABOR! TR !  -> ! DENSITE DE FLUX NET RADIATIF AUX     !
CARGU !   ,NPHAST)   !    !     ! FACES DE BORD                        !
CARGU ! FLCONV(NFABOR! TR !  -> ! DENSITE DE FLUX CONVECTIF AUX FACES  !
CARGU !   ,NPHAST)   !    !     ! DE BORD                              !
CARGU ! HFCONV(NFABOR! TR !  -> ! COEFFICIENT D'ECHANGE FLUIDE AUX     !
CARGU !   ,NPHAST)   !    !     ! FACES DE BORD                        !
CARGU ! W1(NFABOR)   ! TR ! <-> ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "radiat.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , IPH
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ITYPFB(NFABOR,NPHAS) , ITRIFB(NFABOR,NPHAS)
      INTEGER          IZFRAD(NFABOR,NPHAS)
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
      DOUBLE PRECISION SURFBN(NFABOR)
      DOUBLE PRECISION TPAROI(NFABOR,NPHAST), QINCID(NFABOR,NPHAST)
      DOUBLE PRECISION XLAM(NFABOR,NPHAST), EPA(NFABOR,NPHAST)
      DOUBLE PRECISION EPS(NFABOR,NPHAST), FLUNET(NFABOR,NPHAST)
      DOUBLE PRECISION FLCONV(NFABOR,NPHAST), HFCONV(NFABOR,NPHAST)
      DOUBLE PRECISION W1(NFABOR)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
C
      INTEGER          IFAC , NBP , IVAR , IPHAS
      INTEGER          II , JJ , IN , KPT
      INTEGER          IPART , NTRIA3 , NQUAD4 , NSIDED
C
      INTEGER          IZ , IOK , IZONE , NBZONE
      INTEGER          IZONES(NOZRDM) , IMPRAY
C
      INTEGER          II1 , II2 , LPOS , N1 , N2
      CHARACTER        FICH*80 , NAME*80 , ENTETE*7
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
C 1.  A t-on quelquechose a voir sur le maillage frontiere ?
C=======================================================================
C
      IOK = 0
      IPHAS = IRAPHA(IPH)
      DO IZ = 1,NBRAYF
        IF (IRAYVF(IZ,IPHAS).EQ.1) IOK = IOK + 1
      ENDDO
      IF (IOK.EQ.0) RETURN
C
      ENTETE = ' '
      WRITE(ENTETE,'(A5,I1,A1)')'bord_',IPHAS,'.'
C
      IMPRAY = IMPAVR
C
C=======================================================================
C 2.  ECRITURE DU MAILLAGE DE PEAU AU FORMAT ENSIGHT GOLD
C=======================================================================
C
C--> On stocke le nombre de zones et leur numero associe
C
      NBZONE = 0
      DO IZ = 1,NOZRDM
        IZONES(IZ) = -1
      ENDDO
C
      DO IFAC = 1,NFABOR
C
        IZONE = IZFRAD(IFAC,IPH)
C
        IF (IZONE.GT.0) THEN
C
          IOK = 0
          DO IZ = 1,NBZONE+1
            IF (IZONE.NE.IZONES(IZ) .AND. IOK.EQ.0) THEN
              IOK = 1
            ELSE IF (IZONES(IZ).EQ.IZONE) THEN
              IOK = 2
            ENDIF
          ENDDO
C
          IF (IOK.EQ.1) THEN
            NBZONE = NBZONE +1
            IF (NBZONE.GT.NOZRDM) THEN
              WRITE(NFECRA,9010)
              RETURN
            ENDIF
            IZONES(NBZONE) = IZONE
          ENDIF
C
        ENDIF
C
      ENDDO
C
      IF (NBZONE.EQ.0) THEN
        WRITE(NFECRA,9020)
        RETURN
      ENDIF
C
C=======================================================================
C 3.  ECRITURE DU MAILLAGE DE PEAU AU FORMAT ENSIGHT GOLD
C=======================================================================
C
C--> Ouverture
C
      OPEN (IMPRAY,FILE= ENTETE // 'geom',
     &      STATUS='UNKNOWN',FORM='FORMATTED',
     &      ACCESS='SEQUENTIAL')
C
C--> Entete
C
      WRITE(IMPRAY,'(A)') 'geometrie de la frontiere'
      WRITE(IMPRAY,'(A)') 'au format ensight gold'
      WRITE(IMPRAY,'(A)') 'node id given'
      WRITE(IMPRAY,'(A)') 'element id given'
C
      IPART = 0
C
      DO IZ = 1,NBZONE
C
        IZONE = IZONES(IZ)
C
        NQUAD4 = 0
        NTRIA3 = 0
        NSIDED = 0
C
        DO IFAC = 1,NFABOR
          IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
            NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
            IF (NBP.EQ.4) THEN
              NQUAD4 = NQUAD4 + 1
            ELSE IF (NBP.EQ.3) THEN
              NTRIA3 = NTRIA3 + 1
            ELSE IF (NBP.GE.5) THEN
              NSIDED = NSIDED + 1
            ENDIF
          ENDIF
        ENDDO
C
C
C-> Part des triangles
C
        IF (NTRIA3.GT.0) THEN
C
          IPART = IPART+1
          WRITE(IMPRAY,'(A)')    'part'
          WRITE(IMPRAY,'(I10)')  IPART
          WRITE(IMPRAY,'(A,I4)') 'faces triangles zone ',IZ
          WRITE(IMPRAY,'(A)')    'coordinates'
          WRITE(IMPRAY,'(I10)')  NTRIA3*3
C
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.EQ.3) THEN
                DO IN = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
                  WRITE(IMPRAY,'(I10)')  NODFBR(IN)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
C
          DO IN   = 1,3
            DO IFAC = 1,NFABOR
              IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
                NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
                IF (NBP.EQ.3) THEN
                  DO II = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
                    WRITE(IMPRAY,'(E12.5)') XYZNOD(IN,NODFBR(II))
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          WRITE(IMPRAY,'(A)')  'tria3'
          WRITE(IMPRAY,'(I10)') NTRIA3
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.EQ.3) WRITE(IMPRAY,'(I10)') IFAC
            ENDIF
          ENDDO
C
          KPT = 0
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.EQ.3) THEN
                WRITE(IMPRAY,'(3I10)') KPT+1,KPT+2,KPT+3
                KPT = KPT+3
              ENDIF
            ENDIF
          ENDDO
C
C Fin ecriture part des triangles
        ENDIF
C
C
C-> Part des quadrangles
C
        IF (NQUAD4.GT.0) THEN
C
          IPART = IPART+1
          WRITE(IMPRAY,'(A)')    'part'
          WRITE(IMPRAY,'(I10)')  IPART
          WRITE(IMPRAY,'(A,I4)') 'faces quadrangles zone ',IZ
          WRITE(IMPRAY,'(A)')    'coordinates'
          WRITE(IMPRAY,'(I10)')  NQUAD4*4
C
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.EQ.4) THEN
                DO IN = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
                  WRITE(IMPRAY,'(I10)') NODFBR(IN)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
C
          DO IN   = 1,3
            DO IFAC = 1,NFABOR
              IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
                NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
                IF (NBP.EQ.4) THEN
                  DO JJ = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
                    WRITE(IMPRAY,'(E12.5)') XYZNOD(IN,NODFBR(JJ))
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          WRITE(IMPRAY,'(A)')  'quad4'
          WRITE(IMPRAY,'(I10)') NQUAD4
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.EQ.4) WRITE(IMPRAY,'(I10)') IFAC
            ENDIF
          ENDDO
C
          KPT = 0
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.EQ.4) THEN
                WRITE(IMPRAY,'(4I10)') KPT+1,KPT+2,KPT+3,KPT+4
                KPT = KPT+4
              ENDIF
            ENDIF
          ENDDO
C
C Fin ecriture part des quadrangles
        ENDIF
C
C--> Geometrie : nsided
C
        IF (NSIDED.GT.0) THEN
C
          IPART = IPART+1
          WRITE(IMPRAY,'(A)')    'part'
          WRITE(IMPRAY,'(I10)')  IPART
          WRITE(IMPRAY,'(A,I4)') 'faces polygones zone ',IZ
          WRITE(IMPRAY,'(A)')    'coordinates'
C
          KPT = 0
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.GE.5) KPT = KPT + NBP
            ENDIF
          ENDDO
          WRITE(IMPRAY,'(I10)')  KPT
C
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.GE.5) THEN
                DO IN = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
                  WRITE(IMPRAY,'(I10)') NODFBR(IN)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
C
          DO IN   = 1,3
            DO IFAC = 1,NFABOR
              IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
                NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
                IF (NBP.GE.5) THEN
                  DO JJ = IPNFBR(IFAC),IPNFBR(IFAC+1)-1
                    WRITE(IMPRAY,'(E12.5)') XYZNOD(IN,NODFBR(JJ))
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
          ENDDO
C
          WRITE(IMPRAY,'(A)')  'nsided'
          WRITE(IMPRAY,'(I10)') NSIDED
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.GE.5) WRITE(IMPRAY,'(I10)') IFAC
            ENDIF
          ENDDO
C
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.GE.5) WRITE(IMPRAY,'(I10)') NBP
            ENDIF
          ENDDO
C
          KPT = 0
          DO IFAC = 1,NFABOR
            IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
              NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
              IF (NBP.GE.5) THEN
                WRITE(IMPRAY,'(100I10)') (KPT+II,II = 1,NBP)
                KPT = KPT+NBP
              ENDIF
            ENDIF
          ENDDO
C
C Fin ecriture Nsided
        ENDIF
C
C Fin boucle sur les zones
      ENDDO
C
      CLOSE(IMPRAY)
C
C=======================================================================
C 4.  ECRITURE FICHIERS VARIABLES AU FORMAT ENSIGHT GOLD
C=======================================================================
C
      IPHAS = IRAPHA(IPH)
C
      DO IVAR = 1,NBRAYF
C
        IF (IRAYVF(IVAR,IPHAS).EQ.1) THEN
C
          IF      (IVAR.EQ.ITPARP) THEN
            DO IFAC = 1,NFABOR
              W1(IFAC) = TPAROI(IFAC,IPH)
            ENDDO
          ELSE IF (IVAR.EQ.IQINCP) THEN
            DO IFAC = 1,NFABOR
              W1(IFAC) = QINCID(IFAC,IPH)
            ENDDO
          ELSE IF (IVAR.EQ.IXLAMP) THEN
            DO IFAC = 1,NFABOR
              W1(IFAC) = XLAM(IFAC,IPH)
            ENDDO
          ELSE IF (IVAR.EQ.IEPAP) THEN
            DO IFAC = 1,NFABOR
              W1(IFAC) = EPA(IFAC,IPH)
            ENDDO
          ELSE IF (IVAR.EQ.IEPSP) THEN
            DO IFAC = 1,NFABOR
              W1(IFAC) = EPS(IFAC,IPH)
            ENDDO
          ELSE IF (IVAR.EQ.IFNETP) THEN
            DO IFAC = 1,NFABOR
              W1(IFAC) = FLUNET(IFAC,IPH)
            ENDDO
          ELSE IF (IVAR.EQ.IFCONP) THEN
            DO IFAC = 1,NFABOR
              W1(IFAC) = FLCONV(IFAC,IPH)
            ENDDO
          ELSE IF (IVAR.EQ.IHCONP) THEN
            DO IFAC = 1,NFABOR
              W1(IFAC) = HFCONV(IFAC,IPH)
            ENDDO
          ENDIF
C
          FICH = ' '
          FICH = ENTETE
          CALL VERLON (FICH,II1,II2,LPOS)
          CALL VERLON (NBRVAF(IVAR,IPHAS),N1,N2,LPOS)
          FICH(II2+1:II2+LPOS) = NBRVAF(IVAR,IPHAS)(N1:N2)
          CALL VERLON (FICH,II1,II2,LPOS)
C
C--> Ouverture
C
          OPEN (IMPRAY,FILE=FICH(II1:II2),
     &          STATUS='UNKNOWN',FORM='FORMATTED',
     &          ACCESS='SEQUENTIAL')
C
C--> Entete
C
          WRITE(IMPRAY,'(A)')  NBRVAF(IVAR,IPHAS)(N1:N2)
C
          IPART = 0
C
          DO IZ = 1,NBZONE
C
            IZONE = IZONES(IZ)
C
            NQUAD4 = 0
            NTRIA3 = 0
            NSIDED = 0
C
            DO IFAC = 1,NFABOR
              IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
                NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
                IF (NBP.EQ.4) THEN
                  NQUAD4 = NQUAD4 + 1
                ELSE IF (NBP.EQ.3) THEN
                  NTRIA3 = NTRIA3 + 1
                ELSE IF (NBP.GE.5) THEN
                  NSIDED = NSIDED + 1
                ENDIF
              ENDIF
            ENDDO
C
C
C--> Variable scalaire
C
C
            IF (NTRIA3.GT.0) THEN
              IPART = IPART +1
              WRITE(IMPRAY,'(A)')    'part'
              WRITE(IMPRAY,'(I10)')  IPART
              WRITE(IMPRAY,'(A)')  'tria3'
              DO IFAC = 1,NFABOR
                IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
                  NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
                  IF(NBP.EQ.3) WRITE(IMPRAY,'(E12.5)') W1(IFAC)
                ENDIF
              ENDDO
            ENDIF
C
            IF (NQUAD4.GT.0) THEN
              IPART = IPART +1
              WRITE(IMPRAY,'(A)')    'part'
              WRITE(IMPRAY,'(I10)')  IPART
              WRITE(IMPRAY,'(A)')  'quad4'
              DO IFAC = 1,NFABOR
                IF (IZFRAD(IFAC,IPH).EQ.IZONE) THEN
                  NBP = IPNFBR(IFAC+1) - IPNFBR(IFAC)
                  IF(NBP.EQ.4) WRITE(IMPRAY,'(E12.5)') W1(IFAC)
                ENDIF
              ENDDO
            ENDIF
C
            IF (NSIDED.GT.0) IPART = IPART +1
C
C Fin boucle sur les zones
          ENDDO
C
          CLOSE(IMPRAY)
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 5.  ECRITURE FICHIER .CASE AU FORMAT ENSIGHT GOLD
C=======================================================================
C
      OPEN (IMPRAY,FILE = ENTETE // 'CASE',
     &      STATUS='UNKNOWN',FORM='FORMATTED',
     &      ACCESS='SEQUENTIAL')
C
      WRITE(IMPRAY,'(A)') 'FORMAT'
      WRITE(IMPRAY,'(A)') 'type:     ensight gold'
      WRITE(IMPRAY,'(A)') 'GEOMETRY'
      WRITE(IMPRAY,'(A)') 'model:    '//ENTETE//'geom'
      WRITE(IMPRAY,'(A)') 'VARIABLE'
C
      IPHAS = IRAPHA(IPH)
C
      DO IVAR = 1,NBRAYF
C
        IF (IRAYVF(IVAR,IPHAS).EQ.1) THEN
C
          FICH = ' '
          FICH = 'scalar per element:'
          CALL VERLON (FICH,II1,II2,LPOS)
          CALL VERLON (NBRVAF(IVAR,IPHAS),N1,N2,LPOS)
          FICH(II2+11:II2+11+LPOS) = NBRVAF(IVAR,IPHAS)(N1:N2)
C
          NAME = ' '
          NAME = ENTETE
          CALL VERLON (FICH,II1,II2,LPOS)
          CALL VERLON (NAME,N1,N2,LPOS)
          FICH(II2+2:II2+1+LPOS) = NAME(N1:N2)
C
          CALL VERLON (FICH,II1,II2,LPOS)
          CALL VERLON (NBRVAF(IVAR,IPHAS),N1,N2,LPOS)
          FICH(II2+1:II2+LPOS) = NBRVAF(IVAR,IPHAS)(N1:N2)
C
          CALL VERLON (FICH,II1,II2,LPOS)
C
          WRITE(IMPRAY,'(A)') FICH(II1:II2)
C
         ENDIF
C
       ENDDO
C
       CLOSE(IMPRAY)
C
C--------
C FORMATS
C--------
C
 9010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DE L''ECRITURE DES SORTIES ENSIGHT    ',/,
     &'@    *********   SUR LE MAILLAGE FRONTIERE DANS LE MODULE    ',/,
     &'@                DE TRANSFERT THERMIQUE RADIATIF             ',/,
     &'@                                                            ',/,
     &'@    SOUS DIMENSIONNEMENT DE NOZRDM                          ',/,
     &'@                                                            ',/,
     &'@    Il faut augmenter NOZRDM dans l''include radiat.h       ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit, mais les sorties ensight         ',/,
     &'@    ne sont pas effectuees                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DE L''ECRITURE DES SORTIES ENSIGHT    ',/,
     &'@    *********   SUR LE MAILLAGE FRONTIERE DANS LE MODULE    ',/,
     &'@                DE TRANSFERT THERMIQUE RADIATIF             ',/,
     &'@                                                            ',/,
     &'@    PAS DE ZONES DE FRONTIERES DETECTEES                    ',/,
     &'@                                                            ',/,
     &'@    Il faut verifier le contenu du tableau IZFRAD           ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit, mais les sorties ensight         ',/,
     &'@    ne sont pas effectuees                                  ',/,
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
