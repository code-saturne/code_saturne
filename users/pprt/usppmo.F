c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
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
                        SUBROUTINE USPPMO
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR
CFONC  UTILISATION OU NON D'UNE PHYSIQUE PARTICULIERE
CFONC
CFONC
CFONC        UNE SEULE PHYSIQUE PARTICULIERE A LA FOIS.
CFONC
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C
C
C***********************************************************************
C
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) THEN
        RETURN
      ENDIF
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 1.  DECLENCHEMENT DE L UTILISATION D'UNE PHYSIQUE PARTICULIERE
C=======================================================================
C
C
C ---- COD3P Flamme de diffusion en chimie complete rapide (3 points)
C        si = -1   modele non utilise
C        si =  0   modele utilise dans les conditions adiabatiques
C        si =  1   modele utilise dans les conditions permeatiques
C
      IPPMOD(ICOD3P) = -1
C
C
C----- CODEQ Flamme de diffusion en chimie rapide vers l'equilibre
C      ATTENTION : la version CODEQ n'EST PAS OPERATIONNELLE
C      ==========
C        si = -1   modele non utilise
      IPPMOD(ICODEQ) = -1
C
C
C----- COEBU Flamme premelangee en Eddy Break Up
C        si = -1   modele non utilise
C        si =  0   modele utilise dans les conditions adiabatiques
C        si =  1   modele utilise dans les conditions permeatiques (H)
C        si =  2   conditions adiabatiques    avec transport de f
C        si =  3   conditions permeatique (H) avec transport de f
C
      IPPMOD(ICOEBU) = -1
C
C
C----- COBML premelange avec le modele Bray - Moss - Libby
C      ATTENTION : la version COBML n'EST PAS OPERATIONNELLE
C      ==========
C        si = -1   modele non utilise
      IPPMOD(ICOBML) = -1
C
C
C----- COLWC non parfaitement premelange Libby Williams
C        si = -1   modele non utilise
C        si =  0   modele a 2 pics dans les conditions adiabatiques
C        si =  1   modele a 2 pics dans les conditions permeatiques
C                    (suppose rayonnement)
C        si =  2   modele a 3 pics dans les conditions adiabatiques
C        si =  3   modele a 3 pics dans les conditions permeatiques
C                    (suppose rayonnement)
C        si =  4   modele a 4 pics dans les conditions adiabatiques
C        si =  5   modele a 4 pics dans les conditions permeatiques
C                    (suppose rayonnement)
C
      IPPMOD(ICOLWC) = -1
C
C
C----- Charbon pulverise avec trois combustibles gazeux et
C        granulometrie.
C        CP3PL Combustible moyen local
C        IPPMOD(ICP...) = 0 : Transport d'H2
C        IPPMOD(ICP...) = 1 : Transport d'H2 + sechage
      IPPMOD(ICP3PL) = -1
C
C
C----- Charbon pulverise couple lagrangien avec trois combustibles
C      gazeux et granulometrie.
C        IPPMOD(ICPL3C) =-1 : Modele non utilise
C        IPPMOD(ICPL3C) = 0 : Transport d'H2
C        IPPMOD(ICPL3C) = 1 : Transport d'H2 + sechage (non operationnel)
C
      IPPMOD(ICPL3C) = -1
C
C
C----- Combustion Fuel
C        IPPMOD(ICFUEL) =-1 : modele non utilise
C        IPPMOD(ICFUEL) = 0 : modele active
      IPPMOD(ICFUEL) = -1
C
C----- Equation sur YCO2 : pour l'instant on le met ici mais il
C                          faudrait le deplacer, mais ou?
C      Valable pour le charbon et pour le Fuel mais non operationnel
      IEQCO2 = 0
C
C----- COMPF compressible sans choc
C      ==========
C        si = -1   modele non utilise
C        si = 0    modele active
      IPPMOD(ICOMPF) = -1
C
C----- VERSIONS ELECTRIQUES
C        Equation de l'energie obligatoire --> |IPPMOD(IEL...)| >= 1
C        + Possibilite de constituants
C
C       ELJOU : Effet Joule
C        IPPMOD(IELJOU) = 1 : Potentiel reel
C        IPPMOD(IELJOU) = 2 : Potentiel complexe
C        IPPMOD(IELJOU) = 3 : Potentiel reel     + CDL Transfo
C        IPPMOD(IELJOU) = 4 : Potentiel complexe + CDL Transfo
C
      IPPMOD(IELJOU) = -1
C
C       ELARC : Arc electrique
C        IPPMOD(IELARC) = 1 : Potentiel electrique
C        IPPMOD(IELARC) = 2 : Potentiel electrique +
C                             Potentiel vecteur (=>3D)
C
      IPPMOD(IELARC) = -1
C
C       ELION : Mobilite ionique
C        IPPMOD(IELION) = 1 : Potentiel electrique
C
C       ATTENTION : la version ELION n'EST PAS OPERATIONNELLE
C       ==========
C
      IPPMOD(IELION) = -1
C
C
C=======================================================================
C 2.  CHOIX DU FICHIER THERMOCHIMIE DANS LE CAS DE LA COMBUSTION GAZ
C=======================================================================
C
C
C-----Si INDJON=1 on utilise une tabulation ENTH-TEMP calculee par JANAF
C     sinon, l'utilisateur doit lui meme fournir sa propre tabulation
C
      INDJON = 1
C
C----
C FORMATS
C----
      RETURN
      END
c@z
