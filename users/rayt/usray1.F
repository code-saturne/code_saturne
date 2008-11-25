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
                        SUBROUTINE USRAY1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC
CFONC    ============================================================
CFONC    ROUTINE UTILISATEUR : INITIALISATION DES COMMONS RAYONNEMENT
CFONC    ============================================================
CFONC
CFONC
CFONC
CFONC   CALCUL DES FLUX ET DU TERME SOURCE RADIATIFS
CFONC
CFONC   METHODE DES VOLUMES FINIS.
CFONC
CFONC
CFONC  1/ DONNEES DES LUMINANCES ENTRANTES AUX LIMITES DU DOMAINE
CFONC        (C.L : REFLEXION ET EMISSION ISOTROPE)
CFONC
CFONC                               ->  ->           ->
CFONC  2/ CALCUL DE LA LUMINANCE L( X , S ) AU POINT X
CFONC
CFONC                                    D L
CFONC     PAR RESOLUTION DE L'EQUATION : --- = -TK.L +TS
CFONC                                    D S
CFONC                        ->                o
CFONC     OU ENCORE : DIV (L.S ) = -TK.L + TK.L
CFONC
CFONC                                  ->   /    ->  ->  ->
CFONC  3/ CALCUL DES DENSITES DE FLUX  Q = /  L( X , S ).S DOMEGA
CFONC                                     /4.PI
CFONC
CFONC                                       /    ->  ->
CFONC         ET DE L'ABSORPTION       SA= /  L( X , S ).  DOMEGA
CFONC                                     /4.PI
CFONC
CFONC     PAR INTEGRATION DES LUMINANCES SUR LES ANGLES SOLIDES.
CFONC
CFONC     N . B : CA SERT A CALCULER LE TAUX D'ECHAUFFEMENT
CFONC     -----
CFONC                                       /    ->  ->  ->  ->
CFONC  4/ CALCUL DU FLUX INCIDENT QINCID = /  L( X , S ).S . N DOMEGA
CFONC                                     /->->
CFONC        ->                          / S.N >0
CFONC        N NORMALE FLUIDE VERS PAROI
CFONC
CFONC
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
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "radiat.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
      INTEGER          IPHAS
      CHARACTER*2      NUM
C
C
C***********************************************************************
C
C
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C
C
C=======================================================================
C 1. UTILISATION DU MODULE DE TRANSFERTS RADIATIFS
C
C
C     DANS LE CAS DES PHYSIQUES PARTICULIERES (COMBUSTION/CHARBON/
C                                                            ELEC/FIOUL)
C
C         LES RENSEIGNEMENTS NE DOIVENT PAS ETRE FOURNIS
C                            ==============
C         (ils sont completes dans le fichiers de donnees)
C
C     DANS LES AUTRES CAS
C         LES RENSEIGNEMENTS SUIVANTS SONT OBLIGATOIRES
C
C
C
C=======================================================================
C
C--> Par defaut, il n'y a qu'une seule phase traitee.
C
      IPHAS = 1
C
C
      IF( IPPMOD(IPHPAR).LE.1 ) THEN
C
C-->  IRAYON = 0 : PAS DE TRANSFERTS RADIATIFS (PAR DEFAUT)
C            = 1 : TRANSFERTS RADIATIFS, METHODE DES ORDONNEES DISCRETES
C            = 2 : TRANSFERTS RADIATIFS, APPROXIMATION P-1
C
        IRAYON(IPHAS) = 1
C
      ENDIF
C
C=======================================================================
C 2. PARAMETRES DU MODULE DE TRANSFERTS RADIATIFS
C=======================================================================
C
C-->  INDICATEUR SUITE DE CALCUL (LECTURE DU FICHIER SUITE DE RAYONNEMENT)
C     (0      : PAS DE LECTURE D'UN FICHIER SUITE DE RAYONNEMENT
C      1      : RELECTURE D'UN FICHIER SUITE DE RAYONNEMENT
C      ISUITE : RELECTURE D'UN FICHIER SUITE DE RAYONNEMENT SI LE CALCUL FLUIDE EST
C               AUSSI UNE SUITE )
C
      ISUIRD = ISUITE
C
C-->  FREQUENCE DE PASSAGE DANS LE MODULE DE RAYONNEMENT
C
      NFREQR = 1
C
C-->  NOMBRE DE DIRECTIONS : 32 OU 128 (UTILE UNIQUEMENT SI IRAYON=1)
C
      NDIREC = 32
C
C-->  INITIALISATION DU MODE DE CALCUL DU TERME SOURCE RADIATIF EXPLICITE
C     IDIVER = 0 => CALCUL SEMI-ANALYTIQUE
C     IDIVER = 1 => CALCUL CONSERVATIF
C     IDIVER = 2 => CALCUL SEMI-ANALYTIQUE CORRIGE POUR ETRE CONSERVATIF
C     (EN RAYONNEMENT TRANSPARENT, LE CHOIX EST SANS INFLUENCE)
C
      IDIVER = 2
C
C--> NIVEAU D'AFFICHAGE (0,1,2) DES RENSEIGNEMENTS TEMPERATURE DE PAROI
C
      IIMPAR = 1
C
C--> NIVEAU D'AFFICHAGE (0,1,2) DES RENSEIGNEMENTS SOLVEUR
C
      IIMLUM = 0
C
C=======================================================================
C 3. RENSEIGNEMENTS POUR LE POST-PROCESSING
C=======================================================================
C
C ATTENTION : les grandeurs physiques liees au milieu semi-transparent,
C             controlees par le mot-cle IRAYVP, respectent la frequence
C             de sortie chronologique imposée via NTCHR (dans USINI1.F),
C             alors que les grandeurs liees a la frontiere ne sont
C             disponibles que pour le dernier pas de temps.
C             En cas de besoin de visualisation de l'evolution en temps
C             de grandeurs surfaciques contacter l'equipe de
C             developpement.
C
      DO IPHAS = 1, NPHAS
        WRITE(NUM,'(I1)') IPHAS
C
C=======================================================================
C 3.1 VARIABLE DU MILIEU SEMI-TRANSPARENT
C=======================================================================
C
C    ****************************************************************
C    * IL FAUT METTRE LA VALEUR DE IRAYVP A 1 POUR LA VISUALISATION *
C    ****************************************************************
C
C--> TERME SOURCE RADIATIF (ANALYTIQUE/CONSERVATIF/SEMI-ANALYTIQUE)
C
        NBRVAP(ITSRAY,IPHAS) = 'Srad_'//NUM
        IRAYVP(ITSRAY,IPHAS) = -1
C
C--> VECTEUR DENSITE DE FLUX RADIATIF
C
        NBRVAP(IQRAYP,IPHAS) = 'Qrad_'//NUM
        IRAYVP(IQRAYP,IPHAS) = -1
C
C--> PART DE L'ABSORPTION DANS LE TERME SOURCE RADIATIF
C
        NBRVAP( IABSP,IPHAS) = 'Absorp_'//NUM
        IRAYVP( IABSP,IPHAS) = -1
C
C--> PART DE L'EMISSION DANS LE TERME SOURCE RADIATIF
C
        NBRVAP( IEMIP,IPHAS) = 'Emiss_'//NUM
        IRAYVP( IEMIP,IPHAS) = -1
C
C--> COEFFICIENT D'ABSORPTION DU MILIEU SEMI-TRANSPARENT
C
        NBRVAP(  ICAKP,IPHAS) = 'CoefAb_'//NUM
        IRAYVP(  ICAKP,IPHAS) = -1
C
C=======================================================================
C 3.2 VARIABLE SUR LES FRONTIERES DE TYPE PAROI DU DOMAINE DE CALCUL
C=======================================================================
C
C    ****************************************************************
C    * IL FAUT METTRE LA VALEUR DE IRAYVF A 1 POUR LA VISUALISATION *
C    ****************************************************************
C
C--> TEMPERATURE DES FACES FRONTIERES DE PAROI
C
        NBRVAF(ITPARP,IPHAS) = 'Temp_paroi_'//NUM
        IRAYVF(ITPARP,IPHAS) = -1
C
C--> FLUX INCIDENT RADIATIF RECU PAR LES FACES FRONTIERES DE PAROI
C
        NBRVAF(IQINCP,IPHAS) = 'Flux_incident_'//NUM
        IRAYVF(IQINCP,IPHAS) = -1
C
C--> CONDUCTIVITE THERMIQUES DES FACES FRONTIERES DE PAROIS
C
        NBRVAF(IXLAMP,IPHAS) = 'Conductivite_th_'//NUM
        IRAYVF(IXLAMP,IPHAS) = -1
C
C--> EPAISSEUR DES FACES FRONTIERES DE PAROIS
C
        NBRVAF( IEPAP,IPHAS) = 'Epaisseur_'//NUM
        IRAYVF( IEPAP,IPHAS) = -1
C
C--> EMISSIVITE DES FACES FRONTIERES DE PAROIS
C
        NBRVAF( IEPSP,IPHAS) = 'Emissivite_'//NUM
        IRAYVF( IEPSP,IPHAS) = -1
C
C--> FLUX NET RADIATIF AUX FACES FRONTIERES DE PAROIS
C
        NBRVAF(IFNETP,IPHAS) = 'Flux_net_'//NUM
        IRAYVF(IFNETP,IPHAS) = -1
C
C--> FLUX CONVECTIF AUX FACES FRONTIERES DE PAROIS
C
        NBRVAF(IFCONP,IPHAS) = 'Flux_convectif_'//NUM
        IRAYVF(IFCONP,IPHAS) = -1
C
C--> COEFFICIENT D'ECHANGE CONVECTIF AUX FACES FRONTIERES DE PAROIS
C
        NBRVAF(IHCONP,IPHAS) = 'Coef_ech_convectif_'//NUM
        IRAYVF(IHCONP,IPHAS) = -1
C
C=======================================================================
C
      ENDDO
C
      RETURN
C
      END
c@z
