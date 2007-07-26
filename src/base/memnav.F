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
                        SUBROUTINE MEMNAV
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVISCF , IVISCB , IVISFI , IVISBI ,
     &   IDAM   , IXAM   , IDAG   , IXAG   , IIFACL , IIRESP ,
     &   IDRTP  , IGRDP  , ISMBR  , IROVSD ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    , IW10   , IDFRCX , IFRCHY , IDFRHY ,
     &   ICOEFU , IESFLM , IESFLB , ITRAVA , IXIMPA ,
     &   IUVWK  ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   GESTION MEMOIRE NAVIER-STOKES
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
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
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IIFACL       ! E  ! <-  ! "POINTEUR" SUR IFACLG                !
CARGU ! IIRESP       ! E  ! <-  ! "POINTEUR" SUR IRESPR                !
CARGU ! IVISCF, B    ! E  ! <-  ! "POINTEUR" SUR VISCF, VISCB          !
CARGU ! IVISFI, BI   ! E  ! <-  ! "POINTEUR" SUR VISCFI, VISCBI        !
CARGU ! IDAM, IXAM   ! E  ! <-  ! "POINTEUR" SUR DAM, XAM              !
CARGU ! IDAG, IXAG   ! E  ! <-  ! "POINTEUR" SUR DAG, XAG              !
CARGU ! IDRTP        ! E  ! <-  ! "POINTEUR" SUR DRTP                  !
CARGU ! IGRDP        ! E  ! <-  ! "POINTEUR" SUR GRDP                  !
CARGU ! ISMBR        ! E  ! <-  ! "POINTEUR" SUR SMBR                  !
CARGU ! IROVSD       ! E  ! <-  ! "POINTEUR" SUR ROVSDT                !
CARGU ! IW1,2,...,10 ! E  ! <-  ! "POINTEUR" SUR W1 A W9               !
CARGU ! IDFRCX       ! E  ! <-  ! "POINTEUR" SUR DFRCXT                !
CARGU ! IFRCHY       ! E  ! <-  ! "POINTEUR" SUR FRCHY                 !
CARGU ! IDFRHY       ! E  ! <-  ! "POINTEUR" SUR DFRCHY                !
CARGU !IESFLM, IESFLB! E  ! <-  ! "POINTEUR" SUR ESFLUM ET ESFLUB      !
CARGU ! ITRAVA       ! E  ! <-  ! "POINTEUR" SUR TRAVA                 !
CARGU ! IXIMPA       ! E  ! <-  ! "POINTEUR" SUR IXMPA                 !
CARGU ! IUVWK        ! E  ! <-  ! "POINTEUR" SUR UVWK                  !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CASE LIBRE   !
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CASE LIBRE   !
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU .______________.____._____.______________________________________.
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
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IIFACL , IIRESP
      INTEGER          IVISCF , IVISCB , IVISFI , IVISBI
      INTEGER          IDAM   , IXAM   , IDAG   , IXAG
      INTEGER          IDRTP  , IGRDP  , ISMBR  , IROVSD
      INTEGER          IW1    , IW2    , IW3    , IW4    , IW5    , IW6
      INTEGER          IW7    , IW8    , IW9    , IW10
      INTEGER          IDFRCX , IFRCHY , IDFRHY
      INTEGER          ICOEFU , IESFLM , IESFLB
      INTEGER          ITRAVA , IXIMPA , IUVWK
      INTEGER          IFINIA , IFINRA
C
      INTEGER          IDEBIA , IDEBRA , IPHAS, IMLGRD, IRIJ
C
      INTEGER          IESCAT
      INTEGER          IIRNPN
C
C***********************************************************************
C
C---> INITIALISATION
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
C
      IRIJ = 0
      DO IPHAS = 1, NPHAS
        IF(ITYTUR(IPHAS).EQ.3.AND.IRIJNU(IPHAS).EQ.1) THEN
          IRIJ = 1
        ENDIF
      ENDDO
C
      IMLGRD = 0
      DO IPHAS = 1, NPHAS
        IF(IMGR(IPR(IPHAS)).EQ.1.OR.
     &     IMGR(IU (IPHAS)).EQ.1.OR.
     &     IMGR(IV (IPHAS)).EQ.1.OR.
     &     IMGR(IW (IPHAS)).EQ.1    ) THEN
          IMLGRD = 1
        ENDIF
      ENDDO
C
      IESCAT = 0
      DO IPHAS = 1, NPHAS
        IF(IESCAL(IESTOT,IPHAS).GT.0) THEN
          IESCAT = 1
        ENDIF
      ENDDO
C
      IIFACL =       IDEBIA
      IIRESP =       IIFACL + NFAC*2*IMLGRD
      IFINIA =       IIRESP + NCELET*IMLGRD
C
C Attention :
C   ci-dessous, on fait pointer VISCFI et VISCBI
C   sur VISCF et VISCB respectivement
C   dans les cas ou ils sont identiques.
C
C     Pour la norme de resolp
      IIRNPN = 0
      IF(IRNPNW.EQ.1) IIRNPN = 1
C
      IVISFI =       IDEBRA
      IVISCF =       IVISFI + NFAC  *IRIJ
      IVISBI =       IVISCF + NFAC
      IVISCB =       IVISBI + NFABOR*IRIJ
      IDAM   =       IVISCB + NFABOR
      IXAM   =       IDAM   + NCELET
      IDAG   =       IXAM   + NFAC*2
      IXAG   =       IDAG   + NCELET*IMLGRD
      IDRTP  =       IXAG   + NFAC*2*IMLGRD
      IGRDP  =       IDRTP  + NCELET
      ISMBR  =       IGRDP  + NCELET*3
      IROVSD =       ISMBR  + NCELET
      IW1    =       IROVSD + NCELET
      IW2    =       IW1    + NCELET
      IW3    =       IW2    + NCELET
      IW4    =       IW3    + NCELET
      IW5    =       IW4    + NCELET
      IW6    =       IW5    + NCELET
      IW7    =       IW6    + NCELET
      IW8    =       IW7    + NCELET
      IW9    =       IW8    + NCELET
      IW10   =       IW9    + NCELET
      IDFRCX =       IW10   + NCELET*IIRNPN
      IFRCHY =       IDFRCX + NCELET*3*NPHAS*IPHYDR
      IDFRHY =       IFRCHY + NCELET*NDIM*ICALHY
      ICOEFU =       IDFRHY + NCELET*NDIM*ICALHY
      IESFLM =       ICOEFU + NFABOR*NDIM
      IESFLB =       IESFLM + NFAC*IESCAT
      IXIMPA =       IESFLB + NFABOR*IESCAT
      IF(NTERUP.GT.1) THEN
        IUVWK  =       IXIMPA + NCELET*NDIM*NPHAS
        ITRAVA =       IUVWK  + NCELET*NDIM*NPHAS
        IFINRA =       ITRAVA + NCELET*NDIM*NPHAS
      ELSE
        IUVWK  =       IXIMPA
        ITRAVA =       IUVWK
        IFINRA =       ITRAVA
      ENDIF
C
C Dans une phase d'optimisation memoire, on pourra faire pointer
C  esflum et esflub sur fluint et flubrd dans le cas (rare)
C  avec estimateur total et reconstruction par moindres carres
C
C---> VERIFICATION
C
      CALL IASIZE('MEMNAV',IFINIA)
C     ===========
C
      CALL RASIZE('MEMNAV',IFINRA)
C     ===========
C
      RETURN
      END
c@z
