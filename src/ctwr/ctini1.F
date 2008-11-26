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
                        SUBROUTINE CTINI1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC    INIT DES OPTIONS DES VARIABLES POUR LE MODULE AEROREFRIGERANTS
CFONC       EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USINI1
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
      INCLUDE "ctincl.h"
C
C***********************************************************************
C
C VARIABLES LOCALES
C
      INTEGER IPHAS, II, JJ, ISC, IPP
C
C***********************************************************************
C
C=======================================================================
C 0. VERIFICATION ISCALT, ISCSTH
C=======================================================================
C     L'utilisateur ne doit pas y avoir touche.
C
      DO IPHAS = 1, NPHAS
        IF(ISCALT(IPHAS).NE.-1) THEN
          WRITE(NFECRA,1000)IPHAS,ISCALT(IPHAS)
          CALL CSEXIT (1)
C         ===========
        ENDIF
      ENDDO
      DO II = 1, NSCAPP
        IF(ISCSTH(ISCAPP(II)).NE.-10) THEN
          WRITE(NFECRA,1001)II,ISCAPP(II),ISCAPP(II),ISCSTH(ISCAPP(II))
          CALL CSEXIT (1)
C         ===========
        ENDIF
      ENDDO
C
C
C=======================================================================
C 1. VARIABLES TRANSPORTEES
C=======================================================================
C
      ISCSTH(ITEMP4) = 1
      ISCSTH(IHUMID) = 0
C
      IPHAS = 1
      ISCALT(IPHAS) = ITEMP4
C
      IROVAR(IPHAS) = 1
      IVIVAR(IPHAS) = 0
C
C --> Donnees physiques ou numeriques propres aux scalaires
C
      DO ISC = 1, NSCAPP
C
        JJ = ISCAPP(ISC)
C
        IF (ISCAVR(JJ).LE.0) THEN
          VISLS0(JJ) = VISCL0(IPHSCA(JJ))
        ENDIF
C
        BLENCV(ISCA(JJ)) = 1.D0
C
      ENDDO
C
      DO IPHAS = 1, NPHAS
C
        IPP = IPPRTP(ISCA(ITEMP4))
        NOMVAR(IPP)  = 'Temperature'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
        IPP = IPPRTP(ISCA(IHUMID))
        NOMVAR(IPP)  = 'Humidite'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
      ENDDO
C
C=======================================================================
C 2. ON DONNE LA MAIN A L'UTLISATEUR
C=======================================================================
C
      CALL USCTI1
C     ===========
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
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (AEROREFRIGERANTS) DEMANDEE       ',/,
     &'@                                                            ',/,
     &'@  La valeur de ISCALT est renseignee automatiquement.       ',/,
     &'@                                                            ',/,
     &'@  L''utilisateur ne doit pas la renseigner dans usini1, or  ',/,
     &'@    elle a ete affectee comme suit :                        ',/,
     &'@    ISCALT(',I10   ,') = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (AEROREFRIGERANTS) DEMANDEE       ',/,
     &'@                                                            ',/,
     &'@  Les valeurs de ISCSTH sont renseignees automatiquement.   ',/,
     &'@                                                            ',/,
     &'@  L''utilisateur ne doit pas les renseigner dans usini1, or ',/,
     &'@    pour le scalaire ',I10   ,' correspondant au scalaire   ',/,
     &'@    physique particuliere ',I10   ,' on a                   ',/,
     &'@    ISCSTH(',I10   ,') = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
      END
c@z
