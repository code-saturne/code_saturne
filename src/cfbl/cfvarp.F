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
                        SUBROUTINE CFVARP
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC               INIT DES POSITIONS DES VARIABLES
CFONC             POUR LE COMPRESSIBLE SANS CHOC SELON
CFONC  REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
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
C VARIABLES LOCALES
C
      INTEGER          II, IPHAS, IPROP, IOK, ICCFTH, IMODIF
      INTEGER          IIT(1)
      DOUBLE PRECISION DBLPRE(1)
C
C***********************************************************************
C=======================================================================
C 1. DEFINITION DES POINTEURS
C=======================================================================
C
C
      IF ( IPPMOD(ICOMPF).GE.0 ) THEN
C
        IPROP =0
C
        DO IPHAS = 1, NPHAS
C
C ---- Masse volumique
          IPROP = IPROP + 1
          IRHO(IPHAS) = ISCAPP(IPROP)
C     Alias pour les C.L.
          IRUN(IPHAS) = IRHO(IPHAS)
C
C ---- Energie totale
          IPROP = IPROP + 1
          IENERG(IPHAS) = ISCAPP(IPROP)
C     Alias pour les C.L.
          IRUNH(IPHAS) = IENERG(IPHAS)
C
C ---- Temperature (post)
          IPROP = IPROP + 1
          ITEMPK(IPHAS) = ISCAPP(IPROP)
C
C ---- Phase associee
          IPHSCA(IRHO  (IPHAS)) = IPHAS
          IPHSCA(IENERG(IPHAS)) = IPHAS
          IPHSCA(ITEMPK(IPHAS)) = IPHAS
C
C ---- Viscosite dynamique de reference relative au scalaire IRHO
          IVISLS(IRHO  (IPHAS)) = 0
          VISLS0(IRHO  (IPHAS)) = EPZERO
C
C ---- Viscosite dynamique de reference relative au scalaire ITEMPK
          IVISLS(ITEMPK(IPHAS)) = 0
          VISLS0(ITEMPK(IPHAS)) = EPZERO
C
C ---- Initialisation par defaut de la viscosite en volume (cste)
          IVISCV(IPHAS) = 0
          VISCV0(IPHAS) = 0.D0
C
C
C=======================================================================
C 2. OPTIONS DE CALCUL
C=======================================================================
C
C --> Cv constant ou variable (par defaut : constant)
          ICV(IPHAS) = 0
          CV0(IPHAS) = 0.D0
C
          ICCFTH = -1
          IMODIF = 0
          II     = 1
          IIT(1)    = 1
          DBLPRE(1) = 0.D0
          CALL USCFTH
C         ===========
     & ( II , II ,
     &   II , II , II , II , II , II , II ,
     &   II , II , II , II ,
     &   II , II , II ,
     &   ICCFTH , IMODIF  , IPHAS   ,
     &   II , II , II , II ,
     &   IIT , IIT , IIT , IIT , IIT ,
     &   IIT , IIT , IIT , IIT ,
     &   IIT , IIT , IIT , IIT ,
     &   DBLPRE , DBLPRE , DBLPRE , DBLPRE , DBLPRE , DBLPRE , DBLPRE ,
     &   DBLPRE , DBLPRE , DBLPRE , DBLPRE , DBLPRE , DBLPRE ,
     &   DBLPRE , DBLPRE ,
     &   DBLPRE , DBLPRE , DBLPRE , DBLPRE ,
     &   DBLPRE , DBLPRE , DBLPRE , DBLPRE )
C
C --> Utilisation d'un flux de masse specifique pour la vitesse
C
C     ATTENTION   PAS ENCORE IMPLEMENTE
C     *********   LAISSER IFLMAU(IPHAS) = 0
C
          IFLMAU(IPHAS) = 0
C
        ENDDO
C
C
C=======================================================================
C 3. ON REDONNE LA MAIN A L'UTILISATEUR
C=======================================================================
C
        CALL USCFX2
C       ===========
C
C
C=======================================================================
C 4. TRAITEMENT ET VERIFICATION DES DONNEES FOURNIES PAR L'UTILISATEUR
C=======================================================================
C
C ---- Viscosite dynamique de reference relative au scalaire IENERG
        DO IPHAS = 1, NPHAS
          IF(IVISLS(ITEMPK(IPHAS)).GT.0 .OR. ICV(IPHAS).GT.0) THEN
            IVISLS(IENERG(IPHAS)) = 1
          ELSE
            IVISLS(IENERG(IPHAS)) = 0
          ENDIF
C
          VISLS0(IENERG(IPHAS)) = EPZERO
C
        ENDDO
C
        IOK = 0
C
        DO IPHAS = 1, NPHAS
C
          IF(VISLS0(ITEMPK(IPHAS)).LE.0.D0) THEN
            WRITE(NFECRA,1000) IPHAS, VISLS0(ITEMPK(IPHAS))
            IOK = 1
          ENDIF
C
          IF(VISCV0(IPHAS).LT.0.D0) THEN
            WRITE(NFECRA,2000) IPHAS, VISCV0(IPHAS)
            IOK = 1
          ENDIF

        ENDDO
C
        IF(IOK.GT.0) CALL CSEXIT (1)
C
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********   MODULE COMPRESSIBLE                         ',/,
     &'@                                                            ',/,
     &'@    LA CONDUCTIVITE THERMIQUE (PHASE ',I6,') DOIT ETRE      ',/,
     &'@    UN REEL POSITIF STRICTEMENT                             ',/,
     &'@    ELLE A POUR VALEUR ',E12.4                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier uscfx2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********   MODULE COMPRESSIBLE                         ',/,
     &'@                                                            ',/,
     &'@    LA CONDUCTIVITE THERMIQUE (PHASE ',I6,') DOIT ETRE      ',/,
     &'@                                                            ',/,
     &'@    LA VISCOSITE EN VOLUME (PHASE ',I6,') DOIT ETRE         ',/,
     &'@    UN REEL POSITIF                                         ',/,
     &'@    ELLE A POUR VALEUR ',E12.4                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier uscfx2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
      END
C
c@z
