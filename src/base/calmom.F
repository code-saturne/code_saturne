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
                        SUBROUTINE CALMOM
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , NCEL   , NCELET ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDEVEL , ITUSER , IA     ,
     &   RTP    , DT     , PROPCE ,
     &   RDEVEL , RTUSER , RA )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  CALCUL DE LA MOYENNE TEMPORELLE DE CORRELATIONS
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
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE DE CELLULES                   !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! RTP          ! TR !  -> ! TABLEAUX DES VARIABLES AU PDT COURANT!
CARGU ! (NCELET,*   )!    !     !                                      !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
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
      INCLUDE "optcal.h"
      INCLUDE "numvar.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0 , NCEL   , NCELET
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
      DOUBLE PRECISION RTP(NCELET,*) , DT(NCELET), PROPCE(NCELET,*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          II, IEL, IMOM, ICMOM, IDTCM, ITRAVM, IFINRA
      INTEGER          JMOM, ICUMUT

C***********************************************************************
C
C     PHASE DE CUMUL DES TERMES
C     =========================
C
C Pour un tableau, on reserve la memoire ici.
C
      ITRAVM = IDBRA0
      IFINRA  = ITRAVM + NCELET
      CALL RASIZE('CALMOM',IFINRA)
C
C Boucle et test sur les moments a calculer
      DO IMOM = 1, NBMOMT
        IF(NTCABS.GE.NTDMOM(IMOM)) THEN
C
C   Position dans PROPCE du tableau de cumul des moments
          ICMOM = IPPROC(ICMOME(IMOM))
C
C   Tableau de travail pour recueillir la correlation par cellule
          DO IEL = 1, NCEL
            RA(ITRAVM+IEL-1) = 1.D0
          ENDDO
C
C   Calcul de la correlation (on suppose qu'il y a au moins une variable)
          DO II = 1, IDGMOM(IMOM)
C
C     Si la variable est dans RTP
            IF(IDFMOM(II,IMOM).GT.0) THEN
              DO IEL = 1, NCEL
                RA(ITRAVM+IEL-1) =
     &          RA(ITRAVM+IEL-1) * RTP(IEL,IDFMOM(II,IMOM))
              ENDDO
C     Si la variable est dans PROPCE
            ELSEIF(IDFMOM(II,IMOM).LT.0) THEN
              DO IEL = 1, NCEL
                RA(ITRAVM+IEL-1) =
     &          RA(ITRAVM+IEL-1) * PROPCE(IEL,IPPROC(-IDFMOM(II,IMOM)))
              ENDDO
            ENDIF
C
          ENDDO
C
C   Incrementation du tableau de cumul des correlations
          DO IEL = 1, NCEL
            PROPCE(IEL,ICMOM) =
     &      PROPCE(IEL,ICMOM) + DT(IEL)*RA(ITRAVM+IEL-1)
          ENDDO
C
C   Incrementation du tableau de cumul du temps (ncel ou constante)
C
C     Si plusieurs moyennes partagent le meme cumul temporel,
C     il ne faut incrementer que pour une seule (la premiere)
C
C     On decide si on doit incrementer
          ICUMUT = 1
          IF (IMOM.GT.1) THEN
            DO JMOM = 1, IMOM-1
              IF(IDTMOM(JMOM).EQ.IDTMOM(IMOM)) THEN
                ICUMUT = 0
              ENDIF
            ENDDO
          ENDIF
C
C     On incremente
          IF(ICUMUT.EQ.1) THEN
            IF(IDTMOM(IMOM).GT.0) THEN
              IDTCM = IPPROC(ICDTMO(IDTMOM(IMOM)))
              DO IEL = 1, NCEL
                PROPCE(IEL,IDTCM) = PROPCE(IEL,IDTCM) + DT(IEL)
              ENDDO
            ELSEIF(IDTMOM(IMOM).LT.0) THEN
              IDTCM = -IDTMOM(IMOM)
              DTCMOM(IDTCM) = DTCMOM(IDTCM) + DT(1)
            ENDIF
          ENDIF
C
        ENDIF
      ENDDO
C
      RETURN
C
      END
c@z
