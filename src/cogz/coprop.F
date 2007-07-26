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
                        SUBROUTINE COPROP
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
CFONC               POUR LA COMBUSTION
CFONC         FLAMME DE DIFFUSION ET DE PREMELANGE
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
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
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
      INTEGER       IPROP, ICG, IDIRAC
C
C***********************************************************************
C=======================================================================
C 1. POSITIONNEMENT DES PROPRIETES : PROPCE, PROPFA, PROPFB
C    Physique particuliere : Flamme de diffusion chimie 3 points
C=======================================================================
C
      IF ( IPPMOD(ICOD3P).GE.0 ) THEN
C
C ---> Definition des pointeurs relatifs aux variables d'etat
C
C
        IPROP = IPROPP
C
        IPROP  = IPROP + 1
        ITEMP  = IPROP
        DO ICG = 1, NGAZG
          IPROP    = IPROP + 1
          IYM(ICG) = IPROP
        ENDDO
        IF ( IPPMOD(ICOD3P).EQ.1 .AND. IRAYPP.GT.0 ) THEN
          IPROP = IPROP + 1
          ICKABS= IPROP
          IPROP = IPROP + 1
          IT4M  = IPROP
          IPROP = IPROP + 1
          IT3M  = IPROP
        ENDIF
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
        IPROP                 = NPROCE
C
        IPROP                 = IPROP + 1
        IPPROC(ITEMP)         = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
C
        DO ICG = 1, NGAZG
          IPROP                 = IPROP + 1
          IPPROC(IYM(ICG))      = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
        ENDDO
C
        IF ( IPPMOD(ICOD3P).EQ.1 .AND. IRAYPP.GT.0 ) THEN
C
          IPROP                 = IPROP + 1
          IPPROC(ICKABS)        = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(IT4M)          = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(IT3M)          = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
        ENDIF
C
        NPROCE = IPROP
C
C
C ---> Positionnement dans le tableau PROPFB
C      Au centre des faces de bord
C
        IPROP = NPROFB
        IF ( IPPMOD(ICOD3P).EQ.1 .AND. IRAYPP.GT.0 ) THEN
          DO ICG = 1, NGAZG
            IPROP            = IPROP + 1
            IPPROB(IYM(ICG)) = IPROP
          ENDDO
        ENDIF
        NPROFB = IPROP
C
C
C ---> Positionnement dans le tableau PROPFA
C      Au centre des faces internes (flux de masse)
C
       IPROP = NPROFA
C Exemple INUTILE DANS NOTRE CAS
C       IPROP         = IPROP + 1
C       IPPROF(ITEMP) = IPROP
       NPROFA = IPROP
C
       ENDIF
C
C
C
C=======================================================================
C 2. POSITIONNEMENT DES PROPRIETES : PROPCE, PROPFA, PROPFB
C    Physique particuliere : Flamme de diffusion chimie equilibre
C=======================================================================
C
C      IF ( LPP(IDEQA).OR.LPP(IDEQR) ) THEN
C      ENDIF
C
C
C
C=======================================================================
C 3. POSITIONNEMENT DES PROPRIETES : PROPCE, PROPFA, PROPFB
C    Physique particuliere : Flamme de premelange - Modele EBU
C=======================================================================
C
      IF ( IPPMOD(ICOEBU).GE.0 ) THEN
C
C ---> Definition des pointeurs relatifs aux variables d'etat
C
        IPROP = IPROPP
C
        IPROP  = IPROP + 1
        ITEMP  = IPROP
        DO ICG = 1, NGAZG
          IPROP    = IPROP + 1
          IYM(ICG) = IPROP
        ENDDO
        IF ( ( IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 )
     &       .AND. ( IRAYPP.GE.1 )  ) THEN
          IPROP = IPROP + 1
          ICKABS= IPROP
          IPROP = IPROP + 1
          IT4M  = IPROP
          IPROP = IPROP + 1
          IT3M  = IPROP
        ENDIF
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
        IPROP         = NPROCE
C
        IPROP                 = IPROP + 1
        IPPROC(ITEMP)         = IPROP
        IPPPST                = IPPPST + 1
        IPPPRO(IPROP)         = IPPPST
C
        DO ICG = 1, NGAZG
          IPROP                 = IPROP + 1
          IPPROC(IYM(ICG))      = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
        ENDDO
C
        IF ( ( IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 )
     &       .AND. ( IRAYPP.GE.1 )  ) THEN
C
          IPROP                 = IPROP + 1
          IPPROC(ICKABS)        = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(IT4M)          = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(IT3M)          = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
        ENDIF
C
        NPROCE = IPROP
C
C
C ---> Positionnement dans le tableau PROPFB
C      Au centre des faces de bord
C
        IPROP = NPROFB
        IF ( ( IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 )
     &       .AND. ( IRAYPP.GE.1 )  ) THEN
          DO ICG = 1, NGAZG
            IPROP            = IPROP + 1
            IPPROB(IYM(ICG)) = IPROP
          ENDDO
        ENDIF
        NPROFB = IPROP
C
C
C ---> Positionnement dans le tableau PROPFA
C      Au centre des faces internes (flux de masse)
C
       IPROP = NPROFA
C Exemple INUTILE DANS NOTRE CAS
C       IPROP         = IPROP + 1
C       IPPROF(ITEMP) = IPROP
       NPROFA = IPROP
C
       ENDIF
C
C
C
C=======================================================================
C 4. POSITIONNEMENT DES PROPRIETES : PROPCE, PROPFA, PROPFB
C    Physique particuliere : Flamme de premelange - Modele BML
C=======================================================================
C
C      IF ( IPPMOD(ICOBML).GE.0 ) THEN
C      ENDIF
C
C
C
C=======================================================================
C 5. POSITIONNEMENT DES PROPRIETES : PROPCE, PROPFA, PROPFB
C    Physique particuliere : Flamme de premelange - Modele LWC
C=======================================================================
C
      IF ( IPPMOD(ICOLWC).GE.0 ) THEN
C
C ---> Definition des pointeurs relatifs aux variables d'etat
C
        IPROP = IPROPP
C
        IPROP  = IPROP + 1
        ITEMP  = IPROP
C
        IPROP  = IPROP + 1
        IMAM   = IPROP
C
        IPROP  = IPROP + 1
        ITSC   = IPROP
C
        DO ICG = 1, NGAZG
          IPROP    = IPROP + 1
          IYM(ICG) = IPROP
        ENDDO
C
        DO IDIRAC = 1, NDIRAC
C masse volumique Locale
          IPROP         = IPROP + 1
          IRHOL(IDIRAC) = IPROP
C Temperature L.
          IPROP         = IPROP + 1
          ITEML(IDIRAC) = IPROP
C Fraction de Melange L.
          IPROP         = IPROP + 1
          IFMEL(IDIRAC) = IPROP
C Fraction Massique L.
          IPROP         = IPROP + 1
          IFMAL(IDIRAC) = IPROP
C Amplitude L.
          IPROP         = IPROP + 1
          IAMPL(IDIRAC) = IPROP
C Terme Source Chimique L.
          IPROP         = IPROP + 1
          ITSCL(IDIRAC) = IPROP
C MAsse Molaire L.
          IPROP         = IPROP + 1
          IMAML(IDIRAC) = IPROP
        ENDDO
C
        IF ( ( IPPMOD(ICOLWC).EQ.1 .OR. IPPMOD(ICOLWC).EQ.3
     &                             .OR. IPPMOD(ICOLWC).EQ.5 )
     &       .AND. ( IRAYPP.GE.1 )  ) THEN
          IPROP = IPROP + 1
          ICKABS= IPROP
          IPROP = IPROP + 1
          IT4M  = IPROP
          IPROP = IPROP + 1
          IT3M  = IPROP
        ENDIF
C
C ----  Nb de variables algebriques (ou d'etat)
C         propre a la physique particuliere NSALPP
C         total NSALTO
C
        NSALPP = IPROP - IPROPP
        NSALTO = IPROP
C
C
C ---> Positionnement dans le tableau PROPCE
C
        IPROP         = NPROCE
C
        IPROP         = IPROP + 1
        IPPROC(ITEMP) = IPROP
        IPPPST        = IPPPST + 1
        IPPPRO(IPROP) = IPPPST
C
        IPROP         = IPROP + 1
        IPPROC(IMAM)  = IPROP
        IPPPST        = IPPPST + 1
        IPPPRO(IPROP) = IPPPST
C
        IPROP         = IPROP + 1
        IPPROC(ITSC)  = IPROP
        IPPPST        = IPPPST + 1
        IPPPRO(IPROP) = IPPPST
C
        DO ICG = 1, NGAZG
          IPROP            = IPROP + 1
          IPPROC(IYM(ICG)) = IPROP
          IPPPST           = IPPPST + 1
          IPPPRO(IPROP)    = IPPPST
        ENDDO
C
        DO IDIRAC = 1, NDIRAC
          IPROP                 = IPROP + 1
          IPPROC(IRHOL(IDIRAC)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(ITEML(IDIRAC)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(IFMEL(IDIRAC)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(IFMAL(IDIRAC)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(IAMPL(IDIRAC)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(ITSCL(IDIRAC)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(IMAML(IDIRAC)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
        ENDDO
C
        IF ( ( IPPMOD(ICOLWC).EQ.1 .OR. IPPMOD(ICOLWC).EQ.3
     &                             .OR. IPPMOD(ICOLWC).EQ.5 )
     &       .AND. ( IRAYPP.GE.1 )  ) THEN
C
          IPROP                 = IPROP + 1
          IPPROC(ICKABS)        = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(IT4M)          = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IPROP                 = IPROP + 1
          IPPROC(IT3M)          = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
        ENDIF
C
        NPROCE = IPROP
C
C
C ---> Positionnement dans le tableau PROPFB
C      Au centre des faces de bord
C
C
C
C ---> Positionnement dans le tableau PROPFB
C      Au centre des faces de bord
C
        IPROP = NPROFB
        IF ( ( IPPMOD(ICOLWC).EQ.1 .OR. IPPMOD(ICOLWC).EQ.3
     &                             .OR. IPPMOD(ICOLWC).EQ.5)
     &       .AND. ( IRAYPP.GE.1 )  ) THEN
          DO ICG = 1, NGAZG
            IPROP            = IPROP + 1
            IPPROB(IYM(ICG)) = IPROP
          ENDDO
        ENDIF
        NPROFB = IPROP
C
C ---> Positionnement dans le tableau PROPFA
C      Au centre des faces internes (flux de masse)
C
        IPROP = NPROFA
C Exemple INUTILE DANS NOTRE CAS
C       IPROP         = IPROP + 1
C       IPPROF(ITEMP) = IPROP
        NPROFA = IPROP
C
      ENDIF
C
C
      RETURN
C
      END
c@z
