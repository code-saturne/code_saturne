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
                        SUBROUTINE CFPROP
C                       *****************
C     -------------------------------------------------------------
     & ( IPROPP , IPPPST )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC      INIT DES POSITIONS DES VARIABLES D'ETAT POUR
CFONC               POUR LE COMPRESSIBLE SANS CHOC
CFONC          (DANS VECTEURS PROPCE, PROPFA, PROPFB)
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IPROPP       ! E  ! <-> ! NUMERO DE LA DERNIERE PROPRIETE      !
CARGU !              !    !     !  (LES PROPRIETES SONT DANS PROPCE,   !
CARGU !              !    !     !   PROPFA ou PRPFB)                   !
CARGU ! IPPPST       ! E  ! <-> ! POINTEUR INDIQUANT LE RANG DE LA     !
CARGU !              !    !     !  DERNIERE GRANDEUR DEFINIE AUX       !
CARGU !              !    !     !  CELLULES (RTP,PROPCE...) POUR LE    !
CARGU !              !    !     !  POST TRAITEMENT                     !
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
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER       IPROPP, IPPPST
C
C VARIABLES LOCALES
C
      INTEGER       IPROP, IPP, IPHAS
C
C***********************************************************************
C=======================================================================
C 1. POSITIONNEMENT DES PROPRIETES : PROPCE, PROPFA, PROPFB
C    Physique particuliere : Compressible sans choc
C=======================================================================
C
      IF ( IPPMOD(ICOMPF).GE.0 ) THEN
C
C ---> Definition des pointeurs relatifs aux variables d'etat
C
C
        IPROP = IPROPP
C
C  Proprietes des phases : CV s'il est variable
        DO IPHAS = 1, NPHAS
          IF(ICV(IPHAS).NE.0) THEN
            IPROP         = IPROP + 1
            ICV   (IPHAS) = IPROP
          ENDIF
        ENDDO
C
C  Proprietes des phases : Viscosite en volume
        DO IPHAS = 1, NPHAS
          IF(IVISCV(IPHAS).NE.0) THEN
            IPROP         = IPROP + 1
            IVISCV(IPHAS) = IPROP
          ENDIF
        ENDDO
C
C   Flux de masse specifique pour la vitesse (si on en veut un)
        DO IPHAS = 1, NPHAS
          IF(IFLMAU(IPHAS).GT.0) THEN
            IPROP         = IPROP + 1
            IFLUMA(IU  (IPHAS)) = IPROP
            IFLUMA(IV  (IPHAS)) = IPROP
            IFLUMA(IW  (IPHAS)) = IPROP
          ENDIF
        ENDDO
C
C    Flux de Rusanov au bord pour Qdm et E
        DO IPHAS = 1, NPHAS
          IPROP         = IPROP + 1
          IFBRHU(IPHAS) = IPROP
          IPROP         = IPROP + 1
          IFBRHV(IPHAS) = IPROP
          IPROP         = IPROP + 1
          IFBRHW(IPHAS) = IPROP
          IPROP         = IPROP + 1
          IFBENE(IPHAS) = IPROP
        ENDDO
C
C
C ----  Nb de variables algebriques (ou d'etat)
C         propre a la physique particuliere NSALPP
C         total NSALTO
C
        NSALPP = IPROP - IPROPP
        NSALTO = IPROP
C
C ----  On renvoie IPROPP au cas ou d'autres proprietes devraient
C         etre numerotees ensuite
C
        IPROPP = IPROP
C
C
C ---> Positionnement dans le tableau PROPCE
C      et reperage du rang pour le post-traitement
C
        IPROP = NPROCE
C
        DO IPHAS = 1, NPHAS
C
          IF(ICV(IPHAS).GT.0) THEN
            IPROP                 = IPROP + 1
            IPPROC(ICV   (IPHAS)) = IPROP
            IPPPST                = IPPPST + 1
            IPPPRO(IPROP)         = IPPPST
          ENDIF
C
          IF(IVISCV(IPHAS).GT.0) THEN
            IPROP                 = IPROP + 1
            IPPROC(IVISCV(IPHAS)) = IPROP
            IPPPST                = IPPPST + 1
            IPPPRO(IPROP)         = IPPPST
          ENDIF
C
        ENDDO
C
        NPROCE = IPROP
C
C
C ---> Positionnement dans le tableau PROPFB
C      Au centre des faces de bord
C
        IPROP = NPROFB
C
        DO IPHAS = 1, NPHAS
          IPROP                 = IPROP + 1
          IPPROB(IFBRHU(IPHAS)) = IPROP
          IPROP                 = IPROP + 1
          IPPROB(IFBRHV(IPHAS)) = IPROP
          IPROP                 = IPROP + 1
          IPPROB(IFBRHW(IPHAS)) = IPROP
          IPROP                 = IPROP + 1
          IPPROB(IFBENE(IPHAS)) = IPROP
        ENDDO
C
        NPROFB = IPROP
C
C
C ---> Positionnement dans le tableau PROPFA
C      Au centre des faces internes (flux de masse)
C
        IPROP = NPROFA
C
        DO IPHAS = 1, NPHAS
          IF(IFLMAU(IPHAS).GT.0) THEN
            IPROP                     = IPROP + 1
            IPPROF(IFLUMA(IU(IPHAS))) = IPROP
          ENDIF
        ENDDO
C
        NPROFA = IPROP
C
C=======================================================================
C 2. ENTREES SORTIES (entsor.h)
C=======================================================================
C
C     Comme pour les autres variables,
C       si l'on n'affecte pas les tableaux suivants,
C       les valeurs par defaut seront utilisees
C
C     NOMVAR( ) = nom de la variable
C     ICHRVR( ) = sortie chono (oui 1/non 0)
C     ILISVR( ) = suivi listing (oui 1/non 0)
C     IHISVR( ) = sortie historique (nombre de sondes et numeros)
C     si IHISVR(.,1)  = -1 sortie sur toutes les sondes
C
C     NB : Seuls les 8 premiers caracteres du nom seront repris dans le
C          listing le plus detaille
C
C-->  chaleur specifique a volume constant
        IF(ICV   (IPHAS).GT.0) THEN
          IPP = IPPPRO(IPPROC(ICV   (IPHAS)))
          NOMVAR(IPP)   = 'ch. sp volume constant'
          ICHRVR(IPP)   = 0
          ILISVR(IPP)   = 0
          IHISVR(IPP,1) = 0
        ENDIF
C
C-->  viscosite laminaire
        IF(IVISCV(IPHAS).GT.0) THEN
          IPP = IPPPRO(IPPROC(IVISCV(IPHAS)))
          NOMVAR(IPP)   = 'visc. volume'
          ICHRVR(IPP)   = 0
          ILISVR(IPP)   = 0
          IHISVR(IPP,1) = 0
        ENDIF
C
      ENDIF
C
C
      RETURN
      END
c@z
