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
                        SUBROUTINE MEMDTV
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IVISCF , IVISCB , IDAM   , ICOFBD , IW1   , IW2    , IW3     ,
     &   ICOFBR , IGRARX , IGRARY , IGRARZ , IWCF  ,
     &   IPTLRO , IPPMCF ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   GESTION MEMOIRE DTTVAR
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
CARGU ! IVISCF, B    ! E  ! <-  ! "POINTEUR" SUR VISCF, VISCB          !
CARGU ! IDAM         ! E  ! <-  ! "POINTEUR" SUR DAM                   !
CARGU ! ICOFBD       ! E  ! <-  ! "POINTEUR" SUR COFBDT                !
CARGU ! IW1          ! E  ! <-  ! "POINTEUR" SUR W1                    !
CARGU ! IW2          ! E  ! <-  ! "POINTEUR" SUR W2                    !
CARGU ! IW3          ! E  ! <-  ! "POINTEUR" SUR W3                    !
CARGU ! IGRARX,Y,Z   ! E  ! <-  ! "POINTEURS" SUR GRAROX,Y,Z (IPTLRO=1)!
CARGU ! IWCF         ! E  ! <-  ! "POINTEUR" SUR WCF                   !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              ! E  !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              ! E  !     !  DANS IA EN SORTIE                   !
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
C***********************************************************************
C
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL  , NFAC  , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR, NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          IVISCF , IVISCB , IDAM  , ICOFBD, ICOFBR
      INTEGER          IW1    , IW2    , IW3
      INTEGER          IGRARX, IGRARY, IGRARZ, IWCF
      INTEGER          IFINIA , IFINRA
      INTEGER          IPTLRO , IPPMCF
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IIPMCF, IIROCF
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
      IFINIA =       IDEBIA
C
C  Ceci, pour le compressible, n'est pas tres satisfaisant,
C  mais on se calque sur IPTLRO (certains tableaux n'existent
C  pas toujours)
C
      IIPMCF = 0
      IF(IPPMCF.GE.0) THEN
        IIPMCF = 1
      ENDIF
C
      IIROCF = 0
      IF(IPTLRO.EQ.1.OR.IPPMCF.GE.0) THEN
        IIROCF = 1
      ENDIF
C
      IVISCF =       IDEBRA
      IVISCB =       IVISCF + NFAC
      IDAM   =       IVISCB + NFABOR
      ICOFBD =       IDAM   + NCELET
      IW1    =       ICOFBD + NFABOR
      IW2    =       IW1    + NCELET
      IW3    =       IW2    + NCELET
      ICOFBR =       IW3    + NCELET
      IGRARX =       ICOFBR + NFABOR*IIROCF
      IGRARY =       IGRARX + NCELET*IIROCF
      IGRARZ =       IGRARY + NCELET*IIROCF
      IWCF   =       IGRARZ + NFABOR*IIROCF
      IFINRA =       IGRARZ + NCELET*IIPMCF
C
C---> VERIFICATION
C
      CALL IASIZE('MEMDTV',IFINIA)
C     ===========
C
      CALL RASIZE('MEMDTV',IFINRA)
C     ===========
C
      RETURN
      END
c@z
