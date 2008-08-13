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
                        SUBROUTINE RAYLEC
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDEVEL , ITUSER , IA     ,
     &   RAYEXP , RAYIMP , TPAROI , QINCID , FLCONV , HFCONV ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE RAYONNEMENT :
CFONC    --------------------------------------
CFONC
CFONC          Lecture du fichier suite au 1er passage
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
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
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! RAYEXP(NCELET! TR ! <-  ! TERME SOURCE RADIATIF EXPLICITE      !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! RAYIMP(NCELET! TR ! <-  ! TERME SOURCE RADIATIF IMPLICITE      !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! TPAROI(NFABOR! TR !  -  ! TEMPERATURE DE PAROI EN KELVIN       !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! QINCID(NFABOR! TR ! <-  ! DENSITE DE FLUX RADIATIF AUX BORDS   !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! FLCONV(NFABOR! TR ! <-  ! DENSITE DE FLUX CONVECTIF AUX FACES  !
CARGU !   ,NPHAST)   !    !     ! DE BORD                              !
CARGU ! HFCONV(NFABOR! TR ! <-  ! COEFFICIENT D'ECHANGE FLUIDE AUX     !
CARGU !   ,NPHAST)   !    !     ! FACES DE BORD                        !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
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
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "radiat.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION RAYEXP(NCELET,NPHAST), RAYIMP(NCELET,NPHAST)
      DOUBLE PRECISION TPAROI(NFABOR,NPHAST), QINCID(NFABOR,NPHAST)
      DOUBLE PRECISION FLCONV(NFABOR,NPHAST), HFCONV(NFABOR,NPHAST)
C
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER        RUBRIQ*64
      CHARACTER        CPHASE(NPHSMX)*2
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IPHAS, IOK
C
      INTEGER          JPHAST
      INTEGER          NCELOK , NFAIOK , NFABOK , NSOMOK
      INTEGER          IERROR , IRTYP  , ITYSUP , NBVAL
      INTEGER          ILECEC , NBERRO , IVERS
      INTEGER          IMPAMR
C
C***********************************************************************
C=======================================================================
C 0 - GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. LECTURE DU FICHIER SUITE
C=======================================================================
C
      IF (ISUIRD.EQ.1) THEN
C
C ---> On code en chaine le numero des phases
C
C     Verifications pour les format et les numero de phase
C       et de scalaire en chaine.
          IF(NPHSMX.GT.99) THEN
            WRITE(NFECRA,9001)NPHSMX
            CALL CSEXIT (1)
          ENDIF
C
C     Codage en chaine de caracteres du numero de la phase
C       Aller jusqu'a NPHAST suffirait
          DO IPHAS = 1, NPHSMX
            WRITE(CPHASE(IPHAS),'(I2.2)')IPHAS
          ENDDO
C
C  ---> Ouverture
C
          WRITE(NFECRA,6000)
C
C     (ILECEC=1:lecture)
          ILECEC = 1
          CALL OPNSUI(FICAMR,LEN(FICAMR),ILECEC,IMPAMR,IERROR)
C         ===========
          IF (IERROR.NE.0) THEN
            WRITE(NFECRA,9011) FICAMR
            CALL CSEXIT (1)
          ENDIF
C
          WRITE(NFECRA,6010)
C
C
C  ---> Type de fichier suite
C        Pourrait porter le numero de version si besoin.
C        On ne se sert pas de IVERS pour le moment
C
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'version_fichier_suite_rayonnement'
          CALL LECSUI(IMPAMR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IVERS,IERROR)
C
          IF (IERROR.NE.0) THEN
            WRITE(NFECRA,9200)FICAMR
            CALL CSEXIT (1)
          ENDIF
C
C
C  ---> Tests
C
          IOK = 0

C     Dimensions des supports
C
          CALL TSTSUI(IMPAMR,NCELOK,NFAIOK,NFABOK,NSOMOK)
C         ===========
          IF (NCELOK.EQ.0) THEN
            WRITE(NFECRA,9210)
            IOK = IOK + 1
          ENDIF
          IF (NFABOK.EQ.0) THEN
            WRITE(NFECRA,9211)
            IOK = IOK + 1
          ENDIF
C
C     Nombre de phases
C
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'nombre_phases_rayt'
          JPHAST = -1
          CALL LECSUI(IMPAMR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                JPHAST,IERROR)
C
          IF(JPHAST.NE.NPHAST.OR.IERROR.NE.0) THEN
            WRITE(NFECRA,1110) JPHAST,NPHAST
            IOK = IOK + 1
          ENDIF
C
C     Arret
          IF (IOK.NE.0) THEN
            CALL CSEXIT (1)
          ENDIF
C
C  ---> Pour test ulterieur si pb : arret
C
          NBERRO = 0
C
C  ---> Lecture des donnees
C
          DO IPHAS = 1, NPHAST
C
C     Aux faces de bord
C
            ITYSUP = 3
            NBVAL  = 1
            IRTYP  = 2
C
            RUBRIQ = 'tparoi_fb_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           TPAROI(1,IPHAS),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'qincid_fb_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           QINCID(1,IPHAS),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'hfconv_fb_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           HFCONV(1,IPHAS),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'flconv_fb_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           FLCONV(1,IPHAS),IERROR)
            NBERRO=NBERRO+IERROR
C
C
C     Aux cellules
C
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
C
            RUBRIQ = 'rayimp_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           RAYIMP(1,IPHAS),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'rayexp_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMR,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           RAYEXP(1,IPHAS),IERROR)
            NBERRO=NBERRO+IERROR
C
         ENDDO
C
C  ---> Si pb : arret
C
          IF(NBERRO.NE.0) THEN
            WRITE(NFECRA,9100)
            CALL CSEXIT (1)
          ENDIF
C
          WRITE(NFECRA,6011)

C  ---> Fermeture du fichier suite
C
          CALL CLSSUI(IMPAMR,IERROR)
C
          IF (IERROR.NE.0) THEN
            WRITE(NFECRA,8011) FICAMR
          ENDIF
C
          WRITE(NFECRA,6099)
C
C Fin détection suite rayonnement
      ENDIF
C
C--------
C FORMATS
C--------
C
 1110 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LECTURE DU FICHIER SUITE RAYONNEMENT        ',/,
     &'@    *********                                               ',/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Le de phases qui rayonnent a ete modifie                ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite rayonnement.                  ',/,
     &'@                                                            ',/,
     &'@               NPHAS                                        ',/,
     &'@  AMONT : ',I10                                              ,/,
     &'@  ACTUEL: ',I10                                              ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 6000 FORMAT (   3X,'** INFORMATIONS SUR LE MODULE DE RAYONNEMENT   ',/,
     &           3X,'   ------------------------------------------  ',/,
     &           3X,' Lecture d''un fichier suite                   '  )
 6010 FORMAT (   3X,'   Debut de la lecture                         '  )
 6011 FORMAT (   3X,'   Fin   de la lecture                         '  )
 6099 FORMAT (   3X,' Fin de la lecture du fichier suite            ',/,
     &'                                                             ',/,
     &'-------------------------------------------------------------',/)
C
 8011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,
     &'@    *********   RAYONNEMENT                                 ',/,
     &'@                                                            ',/,
     &'@    Probleme sur le fichier de nom (',A13,')                ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********   RAYONNEMENT                                 ',/,
     &'@      ERREUR A L''OUVERTURE DU FICHIER SUITE                ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier l''existence et le nom (',A13,') du            ',/,
     &'@        fichier suite dans le repertoire de travail.        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/)
 9200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                    RAYONNEMENT',/,
     &'@      TYPE DE FICHIER INCORRECT                             ',/,
     &'@                                                            ',/,
     &'@    Le fichier ',A13      ,' ne semble pas etre un fichier  ',/,
     &'@      suite rayonnement.                                    ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise correspond bien   ',/,
     &'@        a un fichier suite rayonnement.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********   RAYONNEMENT                                 ',/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    Le nombre de cellules a ete modifie                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9211 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********   RAYONNEMENT                                 ',/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    Le nombre de faces de bord a ete modifie                ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :ARRET A LA LECTURE DU FICHIER SUITE          ',/,
     &'@    *********   RAYONNEMENT                                 ',/,
     &'@      Le nombre de phases    maximal NPHSMX doit etre       ',/,
     &'@        inferieur ou egal a  99, il vaut ici ',I10           ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme raycli.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: ARRET A LA LECTURE DU FICHIER SUITE          ',/,
     &'@    *********   RAYONNEMENT                                 ',/,
     &'@      ERREUR LORS DE LA LECTURE DES DONNEES                 ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
