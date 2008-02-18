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
                        SUBROUTINE MEMDIS
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IIFACL , IIRESP ,
     &   IVISCF , IVISCB , IDAM   , IXAM   , ISMBR  , IROVSD ,
     &   IRTDP  , ICOEFA , ICOEFB ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   GESTION MEMOIRE DISTANCE A LA PAROI
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
CARGU ! IVISCF, B    ! E  ! <-  ! "POINTEUR" SUR VISCF, VISCB          !
CARGU ! IDAM, IXAM   ! E  ! <-  ! "POINTEUR" SUR DAM, XAM              !
CARGU ! IDRTP        ! E  ! <-  ! "POINTEUR" SUR DRTP                  !
CARGU ! ISMBR        ! E  ! <-  ! "POINTEUR" SUR SMBR                  !
CARGU ! IROVSD       ! E  ! <-  ! "POINTEUR" SUR ROVSDT                !
CARGU ! IW1,2,...,9  ! E  ! <-  ! "POINTEUR" SUR W1 A W9               !
CARGU ! ITIMEP       ! E  ! <-  ! "POINTEUR" SUR TIMEPA                !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
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
      INTEGER          IVISCF , IVISCB , IDAM   , IXAM
      INTEGER          ISMBR  , IROVSD
      INTEGER          IRTDP  , ICOEFA , ICOEFB
      INTEGER          IW1    , IW2    , IW3    , IW4    , IW5    , IW6
      INTEGER          IW7    , IW8    , IW9
      INTEGER          IFINIA , IFINRA
C
      INTEGER          IDEBIA, IDEBRA
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
      IIFACL =       IDEBIA
      IIRESP =       IIFACL + NFAC*2*IMGRPY
      IFINIA =       IIRESP + NCELET*IMGRPY
C
      IVISCF =       IDEBRA
      IVISCB =       IVISCF + NFAC
      IDAM   =       IVISCB + NFABOR
      IXAM   =       IDAM   + NCELET
      ISMBR  =       IXAM   + NFAC*2
      IROVSD =       ISMBR  + NCELET
      IRTDP  =       IROVSD + NCELET
      ICOEFA =       IRTDP  + NCELET
      ICOEFB =       ICOEFA + NFABOR
      IW1    =       ICOEFB + NFABOR
      IW2    =       IW1    + NCELET
      IW3    =       IW2    + NCELET
      IW4    =       IW3    + NCELET
      IW5    =       IW4    + NCELET
      IW6    =       IW5    + NCELET
      IW7    =       IW6    + NCELET
      IW8    =       IW7    + NCELET
      IW9    =       IW8    + NCELET
      IFINRA =       IW9    + NCELET
C
C---> VERIFICATION
C
      CALL IASIZE('MEMDIS',IFINIA)
C     ===========
C
      CALL RASIZE('MEMDIS',IFINRA)
C     ===========
C
      RETURN
      END
c@z
