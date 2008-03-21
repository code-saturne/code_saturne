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
                        SUBROUTINE CPPROP
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
CFONC          COMBUSTION CHARBON PULVERISE
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
      INCLUDE "ppincl.h"
      INCLUDE "ihmpre.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER       IPROPP, IPPPST
C
C VARIABLES LOCALES
C
      INTEGER       IPROP, IGE , ICLA, IPROP2
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
      IROM1  = IPROP
      DO IGE = 1, (NGAZE-2*NCHARB)
C ---- Cf. definition de NGAZE dans cplecd.F
        IPROP     = IPROP + 1
        IYM1(IGE) = IPROP
      ENDDO
      IPROP = IPROP + 1
      IMMEL = IPROP
C
      IPROP2 = IPROP
C
C   Phase dispersee (classes de particules)
      DO ICLA = 1, NCLACP
        IPROP        = IPROP2 + ICLA
        ITEMP2(ICLA) = IPROP
        IPROP        = IPROP2 + 1*NCLACP + ICLA
        IX2(ICLA)    = IPROP
        IPROP        = IPROP2 + 2*NCLACP + ICLA
        IROM2(ICLA)  = IPROP
        IPROP        = IPROP2 + 3*NCLACP + ICLA
        IDIAM2(ICLA) = IPROP
        IPROP        = IPROP2 + 4*NCLACP + ICLA
        IGMDCH(ICLA) = IPROP
        IPROP        = IPROP2 + 5*NCLACP + ICLA
        IGMDV1(ICLA) = IPROP
        IPROP        = IPROP2 + 6*NCLACP + ICLA
        IGMDV2(ICLA) = IPROP
        IPROP        = IPROP2 + 7*NCLACP + ICLA
        IGMHET(ICLA) = IPROP
        IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
          IPROP        = IPROP2 + 8*NCLACP + ICLA
          IGMSEC(ICLA) = IPROP
        ENDIF
      ENDDO
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
      IPROP         = NPROCE
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
      DO IGE = 1, (NGAZE-2*NCHARB)
C ---- Cf. definition de NGAZE dans cplecd.F
        IPROP                 = IPROP + 1
        IPPROC(IYM1(IGE))     = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
      ENDDO
C
      IPROP                 = IPROP + 1
      IPPROC(IMMEL)         = IPROP
      IPPPST                = IPPPST + 1
      IPPPRO(IPROP)         = IPPPST
C
      IPROP2 = IPROP
C
C   Phase dispersee (classes de particules)
      DO ICLA = 1, NCLACP
C
        IPROP                 = IPROP2 + ICLA
        IPPROC(ITEMP2(ICLA))  = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
C
        IPROP                 = IPROP2 + 1*NCLACP + ICLA
        IPPROC(IX2(ICLA))     = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
C
        IPROP                 = IPROP2 + 2*NCLACP + ICLA
        IPPROC(IROM2(ICLA))   = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
C
        IPROP                 = IPROP2 + 3*NCLACP + ICLA
        IPPROC(IDIAM2(ICLA))  = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
C
        IPROP                 = IPROP2 + 4*NCLACP + ICLA
        IPPROC(IGMDCH(ICLA))  = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
C
        IPROP                 = IPROP2 + 5*NCLACP + ICLA
        IPPROC(IGMDV1(ICLA))  = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
C
        IPROP                 = IPROP2 + 6*NCLACP + ICLA
        IPPROC(IGMDV2(ICLA))  = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
C
        IPROP                 = IPROP2 + 7*NCLACP + ICLA
        IPPROC(IGMHET(ICLA))  = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
C
        IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
          IPROP                 = IPROP2 + 8*NCLACP + ICLA
          IPPROC(IGMSEC(ICLA))  = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
        ENDIF
C
      ENDDO
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
      IPROP = NPROFB
      NPROFB = IPROP
C
C ---> Positionnement dans le tableau PROPFA
C      Au centre des faces internes (flux de masse)
C
      IPROP = NPROFA
      NPROFA = IPROP
C
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C     Construction de l'indirection entre la numerotation du noyau et XML
      IF (IIHMPR.EQ.1) THEN
        CALL UICPPR (NCLACP, NSALPP, NSALTO, IPPMOD, ICP3PL, IPPPRO,
     &               IPPROC, ITEMP1, IROM1, IYM1, IMMEL,
     &               ITEMP2, IX2, IROM2, IDIAM2, IGMDCH, IGMDV1,
     &               IGMDV2, IGMHET, IGMSEC, ILUMI)

      ENDIF

#endif
C
      RETURN
      END
c@z
