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
                        SUBROUTINE ELVERI
C                       *****************
C     -------------------------------------------------------------
     & ( IOK    )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  VERIFICATION DES PARAMETRES DE CALCUL POUR LE MODULE ELECTRIQUE
CFONC      APRES INTERVENTION UTILISATEUR (COMMONS)
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
      INCLUDE "elincl.h"
      INCLUDE "radiat.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IOK
C
C VARIABLES LOCALES
C
      INTEGER          IPHAS , IMONO , NTF
C
C***********************************************************************
C
C --> Verifications tardives ...
      IF( IPPMOD(IELARC).NE.-1.AND.IPPMOD(IELARC).NE.2) THEN
        WRITE(NFECRA,1000)IPPMOD(IELARC)
        IOK = IOK + 1
      ENDIF
      IF( IPPMOD(IELJOU).NE.-1.AND.IPPMOD(IELJOU).NE.1.AND.
     &                             IPPMOD(IELJOU).NE.2.AND.
     &                             IPPMOD(IELJOU).NE.3.AND.
     &                             IPPMOD(IELJOU).NE.4     ) THEN
        WRITE(NFECRA,1001)IPPMOD(IELJOU)
        IOK = IOK + 1
      ENDIF
      IF( IPPMOD(IELION).NE.-1) THEN
        WRITE(NFECRA,1002)IPPMOD(IELION)
        IOK = IOK + 1
      ENDIF
C
C=======================================================================
C 2. OPTIONS DU CALCUL : TABLEAUX DE ppincl.h
C                                    elincl.h
C=======================================================================
C
C --> Coefficient de relaxation de la masse volumique
      IF( SRROM.LT.0.D0 .OR. SRROM.GE.1.D0) THEN
        WRITE(NFECRA,2000)'SSROM ', SRROM
        IOK = IOK + 1
      ENDIF
C
C --> Recalage des variables electriques
      IF( IELCOR.NE.0 .AND. IELCOR.NE.1) THEN
        WRITE(NFECRA,2001)'IELCOR', IELCOR
        IOK = IOK + 1
      ENDIF
C
C --> Si on recale
      IF( IELCOR.EQ.1) THEN
C
C   -- en arc
        IF(IPPMOD(IELARC).GE.1) THEN
C        Intensite de courant imposee
          IF( COUIMP.LT.0.D0 ) THEN
            WRITE(NFECRA,2002)'COUIMP', COUIMP
            IOK = IOK + 1
          ENDIF
C        Difference de potentiel initiale imposee
          IF( DPOT  .LT.0.D0 ) THEN
            WRITE(NFECRA,2002)'DPOT  ', DPOT
            IOK = IOK + 1
          ENDIF
        ENDIF
C
C   -- en Joule
        IF( IPPMOD(IELJOU).GE.1) THEN
C         Puissance dissipee imposee
          IF( PUISIM.LT.0.D0 ) THEN
            WRITE(NFECRA,2002)'PUISIM', PUISIM
            IOK = IOK + 1
          ENDIF
C         Coefficient multiplicatif correcteur initial
          IF( COEJOU.LT.0.D0 ) THEN
            WRITE(NFECRA,2002)'COEJOU', COEJOU
            IOK = IOK + 1
          ENDIF
c         IF( ABS(COEJOU-1.D0).GT.EPZERO ) THEN
c           WRITE(NFECRA,2003)'COEJOU', COEJOU
c           IOK = IOK + 1
c         ENDIF
C        Difference de potentiel initiale imposee
          IF( DPOT  .LT.0.D0 ) THEN
            WRITE(NFECRA,2002)'DPOT  ', DPOT
            IOK = IOK + 1
          ENDIF
        ENDIF
C
      ENDIF
C
C --> Rho et viscosite variables
      IPHAS = 1
      IF(IROVAR(IPHAS).NE.1) THEN
        WRITE(NFECRA,2010)IROVAR(IPHAS)
        IOK = IOK + 1
      ENDIF
      IF(IVIVAR(IPHAS).NE.1) THEN
        WRITE(NFECRA,2011)IVIVAR(IPHAS)
        IOK = IOK + 1
      ENDIF
C
C --> Cp et visls (ihm) variables
      IPHAS = 1
      IF(ICP(IPHAS).LE.0) THEN
        WRITE(NFECRA,2012)'       ICP(1)',    ICP(IPHAS),
     &                    '       ICP(1)','la chaleur massique     ',
     &                    '       ICP(1)'
        IOK = IOK + 1
      ENDIF
      IF(IVISLS(IHM).LE.0) THEN
        WRITE(NFECRA,2012)'  IVISLS(IHM)',   IVISLS(IHM),
     &                    '  IVISLS(IHM)','la propriete lambda/Cp  ',
     &                    '  IVISLS(IHM)'
        IOK = IOK + 1
      ENDIF
C
C --> Conductivites electriques variables
C     Pour le potentiel reel
      IF(IVISLS(IPOTR).LE.0) THEN
        WRITE(NFECRA,2012)'IVISLS(IPOTR)', IVISLS(IPOTR),
     &                    'IVISLS(IPOTR)','la cond. elec. reelle   ',
     &                    'IVISLS(IPOTR)'
        IOK = IOK + 1
      ENDIF
C     Pour le potentiel imaginaire (Joule)
      IF(IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4 ) THEN
        IF(IVISLS(IPOTI).LE.0) THEN
        WRITE(NFECRA,2012)'IVISLS(IPOTI)', IVISLS(IPOTI),
     &                    'IVISLS(IPOTI)','la cond. elec. imag.    ',
     &                    'IVISLS(IPOTI)'
          IOK = IOK + 1
        ENDIF
      ENDIF
C
C --> verif option sur les transfos
C
      IF(IPPMOD(IELJOU).EQ.3 .OR. IPPMOD(IELJOU).EQ.4 ) THEN
C
        IMONO = 1
        DO NTF=1,NBTRF
C
          IF ( IBRPR(NTF) .NE. 2 .OR.
     &         IBRSEC(NTF).NE. 2       ) THEN
            IMONO = 0
          ENDIF
        ENDDO
C
        IF ( IMONO .EQ. 1 .AND. IPPMOD(IELJOU).NE.3 ) THEN
          WRITE(NFECRA,2020)
          IOK = IOK + 1
        ENDIF
C
        IF ( IMONO .EQ. 0 .AND. IPPMOD(IELJOU).NE.4 ) THEN
          WRITE(NFECRA,2021)
          IOK = IOK + 1
        ENDIF
      ENDIF
C
C ---> verif rayonnement/Arc electrique
C
C
      IF ( IPPMOD(IELARC).GE.1 .AND. IXKABE .EQ. 2
     &                         .AND. IIRAYO .GT. 0 ) THEN
        WRITE(NFECRA,2022)
        IOK = IOK + 1
      ENDIF
C
      IF ( IPPMOD(IELARC).GE.1 .AND. IXKABE .NE. 1
     &                         .AND. IIRAYO .GT. 0 ) THEN
        WRITE(NFECRA,2023)
        IOK = IOK + 1
      ENDIF


C
C=======================================================================
C 5. FORMATS VERIFICATION
C=======================================================================
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    IPPMOD(IELARC) NE PEUT PRENDRE QUE LES VALEURS -1 ET 2  ',/,
     &'@                                                            ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  IPPMOD(IELARC) = -1 : module arc electrique desactive     ',/,
     &'@                 =  2 : module arc electrique 3D            ',/,
     &'@                  avec equation sur le potentiel vecteur A  ',/,
     &'@                                                            ',/,
     &'@                  Aucune autre valeur n''est permise. En    ',/,
     &'@                    particulier, la version avec calcul de B',/,
     &'@                    par le theoreme d''Ampere n''est pas    ',/,
     &'@                    disponible.                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usppmo.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@  IPPMOD(IELJOU) NE PEUT PRENDRE QUE LES VALEURS -1, 1 ET 2 ',/,
     &'@                                                            ',/,
     &'@  IL VAUT ICI ',I10                                          ,/,
     &'@                                                            ',/,
     &'@  IPPMOD(IELJOU) = -1 : module Joule desactive              ',/,
     &'@                 =  1 : module Joule avec potentiel reel    ',/,
     &'@                 =  2 : module Joule avec potentiel complexe',/,
     &'@                                                            ',/,
     &'@                  Aucune autre valeur n''est permise.       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usppmo.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@  IPPMOD(IELION) NE PEUT PRENDRE QUE LA VALEUR -1           ',/,
     &'@                                                            ',/,
     &'@  IL VAUT ICI ',I10                                          ,/,
     &'@                                                            ',/,
     &'@  IPPMOD(IELION) = -1 : module mobilite ionique desactive   ',/,
     &'@                                                            ',/,
     &'@                  Aucune autre valeur n''est permise.       ',/,
     &'@                  Le module n''est pas disponible.          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usppmo.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    ',A6,' DOIT ETRE UN REEL SUPERIEUR OU EGAL A ZERO ET    ',/,
     &'@                             INFERIEUR STRICTEMENT A 1.     ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier useli1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER COMPRIS EGAL A 0 OU 1        ',/,
     &'@                                                            ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier useli1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    ',A6,' DOIT ETRE UN REEL STRICTEMENT POSITIF            ',/,
     &'@                                                            ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier useli1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
c2003 FORMAT(
c    &'@                                                            ',/,
c    &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
c    &'@                                                            ',/,
c    &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
c    &'@    *********                                               ',/,
c    &'@                MODULE ELECTRIQUE                           ',/,
c    &'@                                                            ',/,
c    &'@    ',A6,' DOIT ETRE UN REEL EGAL A 1                       ',/,
c    &'@                                                            ',/,
c    &'@    IL VAUT ICI ',E14.5                                      ,/,
c    &'@                                                            ',/,
c    &'@  Le calcul ne peut etre execute.                           ',/,
c    &'@                                                            ',/,
c    &'@  Verifier useli1.                                          ',/,
c    &'@                                                            ',/,
c    &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
c    &'@                                                            ',/)
C
 2010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    IROVAR(1) DOIT ETRE UN ENTIER EGAL A 1                  ',/,
     &'@                                                            ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  IROVAR = 1 indique que la masse volumique       est       ',/,
     &'@    variable. Elle est donnee par une loi dans uselph ou    ',/,
     &'@    par fichier de donnees.                                 ',/,
     &'@  IROVAR ne doit pas etre modifie par l''utilisateur        ',/,
     &'@    lorsque le module electrique est active.                ',/,
     &'@                                                            ',/,
     &'@  Verifier useli1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    IVIVAR(1) DOIT ETRE UN ENTIER EGAL A 1                  ',/,
     &'@                                                            ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  IVIVAR = 1 indique que la viscosite moleculaire est       ',/,
     &'@    variable. Elle est donnee par une loi dans uselph ou    ',/,
     &'@    par fichier de donnees.                                 ',/,
     &'@  IVIVAR ne doit pas etre modifie par l''utilisateur        ',/,
     &'@    lorsque le module electrique est active.                ',/,
     &'@                                                            ',/,
     &'@  Verifier useli1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2012 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    ',A13,      ' DOIT ETRE UN ENTIER STRICTEMENT POSITIF   ',/,
     &'@                                                            ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  ',A13,      ' > 0 indique que ',A26                        ,/,
     &'@    est variable. Elle est donnee par une loi dans uselph ou',/,
     &'@    par fichier de donnees.                                 ',/,
     &'@  ',A13,      ' ne doit pas etre modifie par                ',/,
     &'@    l''utilisateur lorsque le module electrique est active. ',/,
     &'@                                                            ',/,
     &'@  Verifier useli1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    Tout vos branchements sont mono et vous avez choisi     ',/,
     &'@    un Type d''ecoulement avec transport du potentiel       ',/,
     &'@    imaginaire                                              ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usppmo et votre fichier dp_elec                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2021 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    Certains de vos branchements ne sont pas mono et        ',/,
     &'@    vous avez choisis un Type d''ecoulement sans            ',/,
     &'@    transport du potentiel imaginaire                       ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usppmo et votre fichier dp_elec                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2022 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    VOUS NE POUVEZ PAS ACTIVER LE MODULE RAYONNEMENT        ',/,
     &'@    ET INTRODUIRE UN TERME SOURCE RADIATIF DANS LE          ',/,
     &'@    FICHIER DP_ELEC                                         ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@  Verifier votre fichier dp_elec ou desactiver le module    ',/,
     &'@  rayonnement                                               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2023 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                MODULE ELECTRIQUE                           ',/,
     &'@                                                            ',/,
     &'@    POUR ACTIVER LE MODULE RAYONNEMENT AVEC LE MODULE       ',/,
     &'@    ARC ELECTRIQUE IL FAUT DONNER LE COEFFICIENT            ',/,
     &'@    D''ABSORPTION DANS LE FICHIER DP_ELEC                   ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@  Verifier votre fichier dp_elec ou desactiver le module    ',/,
     &'@  rayonnement                                               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C=======================================================================
C 6. SORTIE
C=======================================================================
C
      RETURN
      END
c@z
