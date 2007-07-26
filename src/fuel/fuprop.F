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
                        SUBROUTINE FUPROP
C                       *****************
C     -------------------------------------------------------------
     & ( IPROPP , IPPPST )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  INIT DES POSITIONS DES VARIABLES D'ETAT SELON
CFONC          COMBUSTION FUEL
CFONC    (DANS VECTEURS PROPCE, PROPFA, PROPFB)
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
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
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
      INTEGER       IPROP, IGE
C
C***********************************************************************
C
C ---> Definition des pointeurs relatifs aux variables d'etat
C
      IPROP = IPROPP
C
C    Phase continue (melange gazeux)
      IPROP   = IPROP + 1
      ITEMP1  = IPROP
      IPROP   = IPROP + 1
      IROM1   = IPROP
      DO IGE = 1, NGAZE
        IPROP     = IPROP +1
        IYM1(IGE) = IPROP
      ENDDO
      IPROP = IPROP + 1
      IMMEL = IPROP
C
C   Phase dispersee (classes de particules)
C
      IPROP  = IPROP
      ITEMP3 = IPROP
      IPROP        = IPROP + 1
      IROM3        = IPROP
      IPROP        = IPROP + 1
      IDIAM3       = IPROP
      IPROP        = IPROP + 1
      IH1HLF       = IPROP
      IPROP        = IPROP + 1
      IGMEVA       = IPROP
      IPROP        = IPROP + 1
      IGMHTF       = IPROP
C
C   Phase dispersee (Luminance)
      IPROP  = IPROP + 1
      ILUMI  = IPROP
C
C ---- Nb de variables algebriques (ou d'etat)
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
C ---> Positionnement dans le tableau PROPCE
C      et reperage du rang pour le post-traitement
C
      IPROP           = NPROCE
C
C    Phase continue (melange gazeux)
      IPROP           = IPROP + 1
      IPPROC(ITEMP1)  = IPROP
      IPPPST          = IPPPST + 1
      IPPPRO(IPROP)   = IPPPST
C
      IPROP           = IPROP + 1
      IPPROC(IROM1)   = IPROP
      IPPPST          = IPPPST + 1
      IPPPRO(IPROP)   = IPPPST
C
      DO IGE = 1, NGAZE
        IPROP              = IPROP + 1
        IPPROC(IYM1(IGE))  = IPROP
        IPPPST             = IPPPST + 1
        IPPPRO(IPROP)      = IPPPST
      ENDDO
C
      IPROP                 = IPROP + 1
      IPPROC(IMMEL)         = IPROP
      IPPPST                = IPPPST + 1
      IPPPRO(IPROP)         = IPPPST
C
C
C   Phase dispersee (classes de particules)
C
      IPROP                 = IPROP +1
      IPPROC(ITEMP3)        = IPROP
      IPPPST                = IPPPST + 1
      IPPPRO(IPROP)         = IPPPST
C
      IPROP                 = IPROP + 1
      IPPROC(IROM3)         = IPROP
      IPPPST                = IPPPST + 1
      IPPPRO(IPROP)         = IPPPST
C
      IPROP                 = IPROP + 1
      IPPROC(IDIAM3)        = IPROP
      IPPPST                = IPPPST + 1
      IPPPRO(IPROP)         = IPPPST
C
      IPROP                 = IPROP + 1
      IPPROC(IH1HLF)        = IPROP
      IPPPST                = IPPPST + 1
      IPPPRO(IPROP)         = IPPPST
C
      IPROP                 = IPROP + 1
      IPPROC(IGMEVA)        = IPROP
      IPPPST                = IPPPST + 1
      IPPPRO(IPROP)         = IPPPST
C
      IPROP                 = IPROP + 1
      IPPROC(IGMHTF)        = IPROP
      IPPPST                = IPPPST + 1
      IPPPRO(IPROP)         = IPPPST
C
C
C    Rayonnement : luminance
C
      IPROP         = IPROP + 1
      IPPROC(ILUMI) = IPROP
      IPPPST        = IPPPST + 1
      IPPPRO(IPROP) = IPPPST
C
      NPROCE = IPROP
C
C
C ---> Positionnement dans le tableau PROPFB
C      Au centre des faces de bord
C
      IPROP  = NPROFB
      NPROFB = IPROP
C
C ---> Positionnement dans le tableau PROPFA
C      Au centre des faces internes (flux de masse)
C
      IPROP  = NPROFA
      NPROFA = IPROP
C
      RETURN
      END
c@z
