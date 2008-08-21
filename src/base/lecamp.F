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
                        SUBROUTINE LECAMP
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCELET , NCEL   ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDEVEL , ITUSER , IA     ,
     &   RTP    ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
c@foncb
CFONC
CFONC FONCTION :
CFONC ----------
CFONC  LECTURE DU FICHIER SUITE PRINCIPAL
CFONC
CFONC  ON S'ARRETE SI ERREUR DE LECTURE OU
CFONC              OU NCEL DIFFERENT
CFONC              OU NTPABS.GE.NTMABS
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                       ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! JPHAS        ! E  !  -> ! NOMBRE DE PHASES DU CALCUL PRECEDENT !
CARGU ! LJTU         ! E  !  -> ! LONGUEUR DE JTURB                    !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! JTURB        ! TE !  -> ! MODELES DE TURB CALCUL PRECEDENT     !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! RTP          ! TR ! <-  ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT        )!
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
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "albase.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NCELET , NCEL
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION RTP(NCELET,*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER        RUBRIQ*64,RUBRIK*64,CAR4*4
      CHARACTER        CINDFP*2
      CHARACTER        CPHASE(NPHSMX)*2
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL
      INTEGER          IPHAS , IVAR  , ISCAL , II    ,  IVERS
      INTEGER          JPHAS , JVAR  , JSCAL , JSCAUS, JSCAPP
      INTEGER          IKIPH , IEIPH , IOMGIP
      INTEGER          IR11IP, IR22IP, IR33IP
      INTEGER          IR12IP, IR13IP, IR23IP
      INTEGER          IOK   , NCELOK, NFAIOK, NFABOK, NSOMOK
      INTEGER          IERROR, IRTYP,  ITYSUP, NBVAL
      INTEGER          NBERRO, ILECEC
      INTEGER          ITURPH, JTURPH, ITYTPH, JTYTPH
      INTEGER          NFMTPH, NFMTSC, NFMTRU
      INTEGER          JTURB(NPHSMX), JTYTUR(NPHSMX), JALE
      INTEGER          IMPAMO
      DOUBLE PRECISION D2S3, D2S3XK
C
C***********************************************************************
C
C=======================================================================
C 0. INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C  ---> Banniere
      WRITE(NFECRA,1000)
C
C     Longueur pour format de print
      NFMTRU = 36
C
C  --->  On code en chaine le numero des phases et scalaires
C
C     Nombre de phases, de scalaires max pour les formats choisis
      NFMTPH = 99
      NFMTSC = 9999
C
C     Indefini a 2 caracteres
      CINDFP='YY'
C
C     Codage en chaine de caracteres du numero de la phase
C       Aller jusqu'a NPHAS suffirait
C       On suppose que les phases ont le meme numero en suite de calcul
      DO IPHAS = 1, MIN(NPHAS,NFMTPH)
        WRITE(CPHASE(IPHAS),'(I2.2)') IPHAS
      ENDDO
      DO IPHAS = MIN(NPHAS,NFMTPH)+1,NPHAS
       CPHASE(IPHAS) = CINDFP
      ENDDO
C
C     Avertissement
      IF(NPHSMX.GT.NFMTPH) THEN
        WRITE(NFECRA,8000)NFMTPH,NPHSMX
      ENDIF
      IF(NSCAMX.GT.NFMTSC) THEN
        WRITE(NFECRA,8001)NFMTSC,NSCAMX
      ENDIF
C
C=======================================================================
C 1. OUVERTURE DU FICHIER OU STOP
C=======================================================================
C
C  ---> Ouverture du fichier (ILECEC = 1 : lecture)
      ILECEC = 1
      CALL OPNSUI(FICAMO,LEN(FICAMO),ILECEC,IMPAMO,IERROR)
C     ===========
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,9100) FICAMO
        CALL CSEXIT (1)
      ENDIF
C
C ---> Debut de la lecture
      WRITE(NFECRA,1100)
C
C=======================================================================
C 2. ENTETES DU FICHIER SUITE OU STOP
C=======================================================================
C
C  --->  Rubrique "fichier suite ppal"
C        Pourrait porter le numero de version si besoin.
C        On ne se sert pas de IVERS pour le moment.
C
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      RUBRIQ = 'version_fichier_suite_principal'
      CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IVERS,
     &            IERROR)
C
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,9200)FICAMO
        CALL CSEXIT (1)
      ENDIF

C  --->  Tests sur les supports
C
      CALL TSTSUI(IMPAMO,NCELOK,NFAIOK,NFABOK,NSOMOK)
C     ===========
C
      IF (NCELOK.EQ.0) THEN
        WRITE(NFECRA,9201)
        CALL CSEXIT (1)
      ENDIF
C
C     Inutile de tester les supports "faces internes" et "faces de bord"
C     ils ne sont pas utilises ici
C      IF (NFAIOK.EQ.0) WRITE(NFECRA,8200)
C      IF (NFABOK.EQ.0) WRITE(NFECRA,8205)
C
C
C  --->  Tests sur les nombres de variables
C
      NBERRO = 0
C
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
C
      RUBRIQ = 'nombre_variables'
      CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,JVAR,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'nombre_scalaires'
      CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,JSCAL,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'nombre_scalaires_us'
      CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,JSCAUS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'nombre_scalaires_pp'
      CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,JSCAPP,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'nombre_phases'
      CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,JPHAS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
C
C  ---> On s'arrete si erreur
C
      IF (NBERRO.NE.0) THEN
        WRITE(NFECRA,9210)
        CALL CSEXIT (1)
      ENDIF
C
C  ---> On previent si des parametres sont differents
C
      IF ( JVAR  .NE.NVAR   .OR. JSCAL .NE.NSCAL  .OR.
     &     JSCAUS.NE.NSCAUS .OR. JSCAPP.NE.NSCAPP .OR.
     &     JPHAS .NE.NPHAS                                 ) THEN
        WRITE(NFECRA,8210)
     &         JVAR, JSCAL, JSCAUS, JSCAPP, JPHAS,
     &         NVAR, NSCAL, NSCAUS, NSCAPP, NPHAS
      ENDIF
C
C ---> Fin de la lecture des dimensions
      WRITE(NFECRA,1299)
C
C=======================================================================
C 3. CORRESPONDANCE DES SCALAIRES OU STOP
C=======================================================================
C
C  ---> Correspondance avec les scalaires du calcul precedent.
C       Le tableau ISCOLD manipule des numero de scalaire "globaux"
C         i.e. de 1 a NSCAL (et non pas de 1 a NSCAUS ou 1 a NSCAPP)
C
C     On complete si besoin le tableau ISCOLD. Si ses valeurs sont
C       egales a -999, c'est que l'utilisateur ne les a pas modifiees,
C       et qu'il ne souhaite donc rien de particulier quant a la
C       correspondance des nouveaux et des anciens scalaires.
C       On etablit alors la correspondance
C       avec l'ancien scalaire (qui est peut etre une variance) de meme
C       numero s'il existe ou a rien du tout sinon.
C
C     Pour bien faire la difference entre scalaires utilisateurs et
C       physiques particulieres, on part du principe suivant :
C         par defaut
C         pour le nouveau scalaire utilisateur j
C           l'ancien scalaire correspondant est
C             l'ancien scalaire utilisateur j s'il existe
C                                             et rien sinon
C         pour le nouveau scalaire physique particuliere j
C           l'ancien scalaire correspondant est pour le moment
C             l'ancien scalaire physique particuliere j s'il existe
C                                                       et rien sinon
C             a terme, on pourra preciser les correspondances plus
C             finement, en particulier si on souhaite faire des suites
C             entre deux modeles differents.
C
C     Attention, on suppose ici que les scalaires sont numerotes dans
C       l'ordre scalaires utilisateurs (NSCAUS) puis scalaires
C       physiques particulieres (NSCAPP)
C
      IF(NSCAUS.GT.0) THEN
        DO II = 1, NSCAUS
          ISCAL = II
          IF(ISCOLD(ISCAL).EQ.-999) THEN
            IF(II.LE.JSCAUS) THEN
              ISCOLD(ISCAL) = II
            ELSE
              ISCOLD(ISCAL) = 0
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
      IF(NSCAPP.GT.0) THEN
        DO II = 1, NSCAPP
          ISCAL = ISCAPP(II)
          IF(ISCOLD(ISCAL).EQ.-999) THEN
            IF(II.LE.JSCAPP) THEN
              ISCOLD(ISCAL) = II+JSCAUS
            ELSE
              ISCOLD(ISCAL) = 0
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C     On verifie les valeurs du tableau ISCOLD. Elles doivent toutes
C       etre inferieures a JSCAL et positives ou nulles.
      IF(NSCAL.GT.0) THEN
        IOK = 0
        DO ISCAL = 1, NSCAL
          IF(ISCOLD(ISCAL).GT.JSCAL) THEN
            WRITE(NFECRA,9320)ISCAL, ISCOLD(ISCAL), JSCAL
            IOK = IOK+1
          ENDIF
          IF(ISCOLD(ISCAL).LT.0) THEN
            WRITE(NFECRA,9321)ISCAL, ISCOLD(ISCAL), JSCAL
            IOK = IOK+1
          ENDIF
        ENDDO
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,9330) IOK
          CALL CSEXIT (1)
        ENDIF
      ENDIF
C
C=======================================================================
C 4. OPTIONS OU STOP
C=======================================================================
C
C  --->  Lecture des options
C
      NBERRO = 0
C
C     Nombre de pas de temps, instant precedent
C
      RUBRIQ = 'nbre_pas_de_temps'
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NTPABS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'instant_precedent'
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 2
      CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,TTPABS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
C
C     Modeles de turbulence
C
      DO IPHAS = 1, MIN(NPHAS,JPHAS)
        RUBRIQ = 'modele_turbulence_phase'//CPHASE(IPHAS)
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
        CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              JTURB(IPHAS),IERROR)
        NBERRO=NBERRO+IERROR
        JTYTUR(IPHAS)=JTURB(IPHAS)/10
      ENDDO
C
C --->  Stop si erreur
      IF (NBERRO.NE.0) THEN
        WRITE(NFECRA,9400)
        CALL CSEXIT (1)
      ENDIF
C
C     Methode ALE
C
      NBERRO = 0
C
      RUBRIQ = 'methode_ALE'
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,JALE,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
C --->  Message si erreur (pas de stop pour compatibilite avec les fichiers anterieurs)
C       -> on n'affiche le message que si IALE=1 (sinon RAS)
      IF (NBERRO.NE.0) THEN
        IF (IALE.EQ.1) WRITE(NFECRA,9401)
        JALE = 0
      ENDIF
C
C --->  Stop si pas de temps incoherent
      IF(NTPABS.GE.NTMABS.AND.INPDT0.EQ.0) THEN
        WRITE(NFECRA,9410)NTPABS,NTMABS
        CALL CSEXIT (1)
      ENDIF
C
C --->  Informations
      WRITE(NFECRA,2410) NTPABS
      WRITE(NFECRA,2411) TTPABS
C
C --->  Donnees modifiees
      DO IPHAS = 1, MIN(NPHAS,JPHAS)
         IF (ITURB(IPHAS) .NE. JTURB(IPHAS))
     &        WRITE(NFECRA,8410) IPHAS,ITURB(IPHAS),
     &                           IPHAS,JTURB(IPHAS),
     &                           IPHAS
      ENDDO
C
C --->  Si le calcul precedent etait en ALE, on DOIT relire les
C         coordonnees des noeuds dans le fichier auxiliaire
      IF (IALE.EQ.1 .AND. JALE.EQ.1) THEN
        IF (ILEAUX.NE.1) THEN
          WRITE(NFECRA,9402)JALE,IALE,ILEAUX
          CALL CSEXIT(1)
        ENDIF
      ENDIF
C
C --->  Fin de la lecture des options
      WRITE(NFECRA,1499)
C
C========================================================
C 5. LECTURE DES VARIABLES
C========================================================
C
      NBERRO = 0
C
C --->  Pression
C     (a priori une seule (non fonction du nombre de phases))
C
      RUBRIQ = 'pression_ce_phase'//CPHASE(1)
      ITYSUP = 1
      NBVAL  = 1
      IRTYP  = 2
      CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &            RTP(1,IPR(1)),IERROR)
      NBERRO=NBERRO+IERROR
C
C
C --->  Vitesse et turbulence
C     On boucle sur les phases communes afin de recuperer
C       les infos existantes dans le fichier suite.
C       Pour les nouvelles phases : deja fait par defaut dans INIVA0
C
      DO IPHAS = 1, MIN(JPHAS,NPHAS)
C
         ITURPH = ITURB(IPHAS)
         JTURPH = JTURB(IPHAS)
         ITYTPH = ITYTUR(IPHAS)
         JTYTPH = JTYTUR(IPHAS)
C
C     Vitesse
C
         ITYSUP = 1
         NBVAL  = 1
         IRTYP  = 2

         RUBRIQ = 'vitesse_u_ce_phase'//CPHASE(IPHAS)
         CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &               RTP(1,IU(IPHAS)),IERROR)
         NBERRO=NBERRO+IERROR
C
         RUBRIQ = 'vitesse_v_ce_phase'//CPHASE(IPHAS)
         CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &               RTP(1,IV(IPHAS)),IERROR)
         NBERRO=NBERRO+IERROR
C
         RUBRIQ = 'vitesse_w_ce_phase'//CPHASE(IPHAS)
         CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &               RTP(1,IW(IPHAS)),IERROR)
         NBERRO=NBERRO+IERROR
C
C
C
C     Turbulence (ke -> ke, ke -> rij, rij -> ke, rij -> rij ...)
C       Quand on ne sait pas deduire les grandeurs turbulence des infos
C         disponibles dans le fichier suite, on ne fait rien :
C         les variables gardent alors les initialisations par defaut
C         faites dans INIVA0
C

C   -- Le nouveau calcul est en k-epsilon

C
         IF (ITYTPH.EQ.2) THEN
C
C     * k-e ou v2f -> k-e
C
          IF(JTYTPH.EQ.2 .OR. JTURPH.EQ.50) THEN
C
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2

            RUBRIQ = 'k_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IK(IPHAS)),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'eps_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IEP(IPHAS)),IERROR)
            NBERRO=NBERRO+IERROR
C
C     * rij -> k-e
C
          ELSEIF(JTYTPH.EQ.3) THEN
C
             IKIPH  = IK (IPHAS)
             IEIPH  = IEP(IPHAS)

             ITYSUP = 1
             NBVAL  = 1
             IRTYP  = 2
C
             RUBRIQ = 'R11_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IKIPH),IERROR)
             NBERRO=NBERRO+IERROR
C
C            La variable epsilon sert de tableau de travail
             RUBRIQ = 'R22_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IEIPH),IERROR)
             NBERRO=NBERRO+IERROR
C
             DO IEL = 1, NCEL
               RTP(IEL,IKIPH) = RTP(IEL,IKIPH) + RTP(IEL,IEIPH)
             ENDDO
C
             RUBRIQ = 'R33_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IEIPH),IERROR)
             NBERRO=NBERRO+IERROR
C
             DO IEL = 1, NCEL
               RTP(IEL,IKIPH) = 0.5D0*(RTP(IEL,IKIPH)+RTP(IEL,IEIPH))
             ENDDO
C
             RUBRIQ = 'eps_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IEIPH),IERROR)
             NBERRO=NBERRO+IERROR
C
C     * k-omega -> k-e
C
          ELSE IF(JTURPH.EQ.60) THEN
C
            IKIPH = IK (IPHAS)
            IEIPH = IEP(IPHAS)
C
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2

            RUBRIQ = 'k_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IKIPH),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'omega_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IEIPH),IERROR)
            NBERRO=NBERRO+IERROR
C           On transforme ensuite omega en epsilon
            DO IEL = 1, NCEL
              RTP(IEL,IEIPH) = CMU*RTP(IEL,IEIPH)*RTP(IEL,IKIPH)
            ENDDO
C
          ENDIF
C         Rq : laminaire -> k-e  (On ne fait rien, deja fait dans iniva0)
C
C
C   -- Le nouveau calcul est en Rij-epsilon
C
        ELSEIF(ITYTPH.EQ.3) THEN
C
C     * k-e ou v2f -> rij
C
          IF (JTYTPH.EQ.2 .OR. JTURPH.EQ.50) THEN
C
             IR11IP=IR11(IPHAS)
             IR22IP=IR22(IPHAS)
             IR33IP=IR33(IPHAS)
             IR12IP=IR12(IPHAS)
             IR13IP=IR13(IPHAS)
             IR23IP=IR23(IPHAS)
C
             ITYSUP = 1
             NBVAL  = 1
             IRTYP  = 2

             RUBRIQ = 'k_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IR11IP),IERROR)
             NBERRO=NBERRO+IERROR
C
             D2S3 = 2.D0/3.D0
             DO IEL = 1, NCEL
                D2S3XK = RTP(IEL,IR11IP)*D2S3
                RTP(IEL,IR11IP) = D2S3XK
                RTP(IEL,IR22IP) = D2S3XK
                RTP(IEL,IR33IP) = D2S3XK
                RTP(IEL,IR12IP) = 0.D0
                RTP(IEL,IR13IP) = 0.D0
                RTP(IEL,IR23IP) = 0.D0
             ENDDO

             RUBRIQ = 'eps_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IEP(IPHAS)),IERROR)
             NBERRO=NBERRO+IERROR
C
C     * rij -> rij
C
          ELSEIF(JTYTPH.EQ.3) THEN
C
           DO II   = 1, 7
            IF (II  .EQ.1) THEN
               RUBRIQ = 'R11_ce_phase'//CPHASE(IPHAS)
               IVAR = IR11(IPHAS)
            ELSEIF (II  .EQ.2) THEN
               RUBRIQ = 'R22_ce_phase'//CPHASE(IPHAS)
               IVAR = IR22(IPHAS)
            ELSEIF (II  .EQ.3) THEN
               RUBRIQ = 'R33_ce_phase'//CPHASE(IPHAS)
               IVAR = IR33(IPHAS)
            ELSEIF (II  .EQ.4) THEN
               RUBRIQ = 'R12_ce_phase'//CPHASE(IPHAS)
               IVAR = IR12(IPHAS)
            ELSEIF (II  .EQ.5) THEN
               RUBRIQ = 'R13_ce_phase'//CPHASE(IPHAS)
               IVAR = IR13(IPHAS)
            ELSEIF (II  .EQ.6) THEN
               RUBRIQ = 'R23_ce_phase'//CPHASE(IPHAS)
               IVAR = IR23(IPHAS)
            ELSEIF (II  .EQ.7) THEN
               RUBRIQ = 'eps_ce_phase'//CPHASE(IPHAS)
               IVAR = IEP(IPHAS)
            ENDIF
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IVAR),IERROR)
            NBERRO=NBERRO+IERROR
           ENDDO
C
C     * k-omega -> rij
C
          ELSE IF (JTURPH.EQ.60) THEN
C
            IR11IP=IR11(IPHAS)
            IR22IP=IR22(IPHAS)
            IR33IP=IR33(IPHAS)
            IR12IP=IR12(IPHAS)
            IR13IP=IR13(IPHAS)
            IR23IP=IR23(IPHAS)
            IEIPH =IEP (IPHAS)
C
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2

            RUBRIQ = 'k_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           RTP(1,IR11IP),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'omega_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           RTP(1,IEIPH),IERROR)
            NBERRO=NBERRO+IERROR
C     On transforme ensuite omega en epsilon
            DO IEL = 1, NCEL
              RTP(IEL,IEIPH) = CMU*RTP(IEL,IEIPH)*RTP(IEL,IR11IP)
            ENDDO
C
            D2S3 = 2.D0/3.D0
            DO IEL = 1, NCEL
              D2S3XK = RTP(IEL,IR11IP)*D2S3
              RTP(IEL,IR11IP) = D2S3XK
              RTP(IEL,IR22IP) = D2S3XK
              RTP(IEL,IR33IP) = D2S3XK
              RTP(IEL,IR12IP) = 0.D0
              RTP(IEL,IR13IP) = 0.D0
              RTP(IEL,IR23IP) = 0.D0
            ENDDO
C
          ENDIF
C       Rq : laminaire -> rij   (On ne fait rien, deja fait dans iniva0)
C
C   -- Le nouveau calcul est en v2f
C
        ELSEIF(ITURPH.EQ.50) THEN
C
C     * k-e -> v2f
C
          IF(JTYTPH.EQ.2) THEN
C
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2

            RUBRIQ = 'k_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IK(IPHAS)),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'eps_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IEP(IPHAS)),IERROR)
            NBERRO=NBERRO+IERROR
C     On laisse pour phi et fb les initialisations de iniva0
C
C     * rij -> v2f
C
          ELSEIF(JTYTPH.EQ.3) THEN
C
             IKIPH  = IK (IPHAS)
             IEIPH  = IEP(IPHAS)

             ITYSUP = 1
             NBVAL  = 1
             IRTYP  = 2
C
             RUBRIQ = 'R11_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IKIPH),IERROR)
             NBERRO=NBERRO+IERROR
C
C            La variable epsilon sert de tableau de travail
             RUBRIQ = 'R22_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IEIPH),IERROR)
             NBERRO=NBERRO+IERROR
C
             DO IEL = 1, NCEL
               RTP(IEL,IKIPH) = RTP(IEL,IKIPH) + RTP(IEL,IEIPH)
             ENDDO
C
             RUBRIQ = 'R33_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IEIPH),IERROR)
             NBERRO=NBERRO+IERROR
C
             DO IEL = 1, NCEL
               RTP(IEL,IKIPH) = 0.5D0*(RTP(IEL,IKIPH)+RTP(IEL,IEIPH))
             ENDDO
C
             RUBRIQ = 'eps_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IEIPH),IERROR)
             NBERRO=NBERRO+IERROR
C     On laisse pour phi et fb l'initialisations de iniva0
C     (le v2 qui intervient dans phi n'est pas vraiment une composante de Rij
C      et qui plus est, on ne saurait pas quelle composante prendre ...)
C
C     * v2f -> v2f
C
          ELSEIF(JTURPH.EQ.50) THEN
C
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2

            RUBRIQ = 'k_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IK(IPHAS)),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'eps_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IEP(IPHAS)),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'phi_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IPHI(IPHAS)),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'fb_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IFB(IPHAS)),IERROR)
            NBERRO=NBERRO+IERROR
C     * k-omega -> v2f
C
          ELSE IF(JTURPH.EQ.60) THEN
C
            IKIPH = IK (IPHAS)
            IEIPH = IEP(IPHAS)
C
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2

            RUBRIQ = 'k_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IKIPH),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'omega_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IEIPH),IERROR)
            NBERRO=NBERRO+IERROR
C           On transforme ensuite omega en epsilon
            DO IEL = 1, NCEL
              RTP(IEL,IEIPH) = CMU*RTP(IEL,IEIPH)*RTP(IEL,IKIPH)
            ENDDO
C     On laisse pour phi et fb l'initialisations de iniva0
C
C
          ENDIF
C         Rq : laminaire -> v2f  (On ne fait rien, deja fait dans iniva0)
C
C
C   -- Le nouveau calcul est en k-omega
C
         ELSE IF (ITURPH.EQ.60) THEN
C
C     * k-e ou v2f -> k-omega
C
          IF(JTYTPH.EQ.2 .OR. JTURPH.EQ.50) THEN
C
            IKIPH  = IK  (IPHAS)
            IOMGIP = IOMG(IPHAS)
C
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2

            RUBRIQ = 'k_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IKIPH),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'eps_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IOMGIP),IERROR)
            NBERRO=NBERRO+IERROR
C           On transforme ensuite epsilon en omega
            DO IEL = 1, NCEL
              RTP(IEL,IOMGIP) = RTP(IEL,IOMGIP)/CMU/RTP(IEL,IKIPH)
            ENDDO
C
C     * rij -> k-omega
C
          ELSEIF(JTYTPH.EQ.3) THEN
C
            IKIPH  = IK  (IPHAS)
            IOMGIP = IOMG(IPHAS)
C
             ITYSUP = 1
             NBVAL  = 1
             IRTYP  = 2
C
             RUBRIQ = 'R11_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IKIPH),IERROR)
             NBERRO=NBERRO+IERROR
C
C            La variable omega sert de tableau de travail
             RUBRIQ = 'R22_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IOMGIP),IERROR)
             NBERRO=NBERRO+IERROR
C
             DO IEL = 1, NCEL
               RTP(IEL,IKIPH) = RTP(IEL,IKIPH) + RTP(IEL,IOMGIP)
             ENDDO
C
             RUBRIQ = 'R33_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IOMGIP),IERROR)
             NBERRO=NBERRO+IERROR
C
             DO IEL = 1, NCEL
               RTP(IEL,IKIPH) = 0.5D0*(RTP(IEL,IKIPH)+RTP(IEL,IOMGIP))
             ENDDO
C
             RUBRIQ = 'eps_ce_phase'//CPHASE(IPHAS)
             CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                   RTP(1,IOMGIP),IERROR)
             NBERRO=NBERRO+IERROR
C           On transforme ensuite epsilon en omega
            DO IEL = 1, NCEL
              RTP(IEL,IOMGIP) = RTP(IEL,IOMGIP)/CMU/RTP(IEL,IKIPH)
            ENDDO
C
C     * k-omega -> k-omega
C
          ELSE IF(JTURPH.EQ.60) THEN
C
            IKIPH = IK  (IPHAS)
            IOMGIP= IOMG(IPHAS)
C
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2

            RUBRIQ = 'k_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IKIPH),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'omega_ce_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RTP(1,IOMGIP),IERROR)
            NBERRO=NBERRO+IERROR
C
          ENDIF
C         Rq : laminaire -> k-omega  (On ne fait rien, deja fait dans iniva0)
C
C
C   -- Le nouveau calcul est en laminaire, longueur de melange ou LES
C           --> rien a lire
C
        ENDIF
C
      ENDDO
C    Fin de la boucle sur les phases communes
C
      IF (NBERRO.NE.0) THEN
         WRITE(NFECRA,9510)
         CALL CSEXIT (1)
      ENDIF
C
C --->  Fin de la lecture des variables pression, vitesse, turbulence
      WRITE(NFECRA,1598)
C
C
C  ---> Scalaires
C
      NBERRO = 0
C
      IF(NSCAL.GT.0) THEN
        DO ISCAL = 1, NSCAL
          IVAR = ISCA(ISCAL)
C         Si le scalaire existait precedemment on le lit
C         sinon on ne fait rien (init par defaut dans INIVA0)
          IF (ISCOLD(ISCAL).GT.0) THEN
            IF(ISCOLD(ISCAL).LE.NFMTSC) THEN
              WRITE(CAR4,'(I4.4)')ISCOLD(ISCAL)
              RUBRIQ = 'scalaire_ce_'//CAR4
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    RTP(1,IVAR),IERROR)
              NBERRO=NBERRO+IERROR
            ELSE
              IERROR= -1
              NBERRO=NBERRO+IERROR
            ENDIF
            IF(IERROR.NE.0) THEN
              RUBRIK = RUBRIQ(1:1)
              DO II = 2, MIN(LEN(RUBRIQ),NFMTRU)
                RUBRIK = RUBRIK(1:II-1)//RUBRIQ(II:II)
              ENDDO
              DO II = MIN(LEN(RUBRIQ),NFMTRU)+1,NFMTRU
                RUBRIK = RUBRIK(1:II-1)//' '
              ENDDO
              WRITE(NFECRA,9511)RUBRIK
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
      IF (NBERRO.NE.0) THEN
        WRITE(NFECRA,9512)
        CALL CSEXIT (1)
      ENDIF
C
      WRITE(NFECRA,1599)
C
C
C  ---> Vitesse de maillage en ALE
C
      IF (IALE.EQ.1 .AND. JALE.EQ.1) THEN
C
        NBERRO = 0
C
        ITYSUP = 1
        NBVAL  = 1
        IRTYP  = 2

        RUBRIQ = 'vit_maillage_u_ce'
        CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              RTP(1,IUMA),IERROR)
        NBERRO=NBERRO+IERROR
C
        RUBRIQ = 'vit_maillage_v_ce'
        CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              RTP(1,IVMA),IERROR)
        NBERRO=NBERRO+IERROR
C
        RUBRIQ = 'vit_maillage_w_ce'
        CALL LECSUI(IMPAMO,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              RTP(1,IWMA),IERROR)
        NBERRO=NBERRO+IERROR
C
        IF (NBERRO.NE.0) THEN
          WRITE(NFECRA,9513)
          DO IEL = 1, NCEL
            RTP(IEL,IUMA) = 0.D0
            RTP(IEL,IVMA) = 0.D0
            RTP(IEL,IWMA) = 0.D0
          ENDDO
        ENDIF
C
        WRITE(NFECRA,1599)
C
      ENDIF
C
C========================================================
C 6. LECTURE D'INFORMATIONS COMPLEMENTAIRES LEGERES
C========================================================
C
C     Pour le moment, il n'y en a pas.
C
C========================================================
C 7. FERMETURE DU FICHIER SUITE PRINCIPAL
C========================================================
C
      CALL CLSSUI(IMPAMO,IERROR)
C     ===========
C
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,8711) FICAMO
      ENDIF
C
      WRITE(NFECRA,1799)

C=======================================================================
C 8. SORTIE
C=======================================================================
C
      RETURN
C
C=======================================================================
C 9. FORMATS
C=======================================================================
C
C --- ETAPES
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(/,
     &     3X,'   LECTURE DU FICHIER SUITE PRINCIPAL                ',/)
 1100 FORMAT(' Debut de la lecture                                  '  )
 1299 FORMAT('   Fin de la lecture des dimensions                   '  )
 1499 FORMAT('   Fin de la lecture des options                      '  )
 1598 FORMAT('   Fin de la lecture des variables                    ',/,
     &       '        pression, vitesse, turbulence                 '  )
 1599 FORMAT('   Fin de la lecture des scalaires                    '  )
 1600 FORMAT('   Fin de la lecture de la vitesse de maillage        ',/,
     &       '                                  (methode ALE)       '  )
 1799 FORMAT(' Fin de la lecture                                    '  )
C
#else
C
 1000 FORMAT(/,
     &     3X,'   READING THE MAIN RESTART FILE                     ',/)
 1100 FORMAT(' Start reading                                        '  )
 1299 FORMAT('   Reading dimensions complete                        '  )
 1499 FORMAT('   Reading options complete                           '  )
 1598 FORMAT('   Reading the pressure, velocity, turbulent          ',/,
     &       '        variables complete                            '  )
 1599 FORMAT('   Reading scalars complete                           '  )
 1600 FORMAT('   Reading the mesh velocity (ALE method) complete    '  )
 1799 FORMAT(' Reading complete                                     '  )
C
#endif
C
C --- INFORMATIONS
C
#if defined(_CS_LANG_FR)
C
 2410 FORMAT
     & ('   Lecture du pas de temps precedent (suite) ',
     &                                                  'NTPABS = ',I10)
 2411 FORMAT
     & ('   Lecture du pas de temps precedent (suite) ',
     &                                                'TTPABS = ',E12.4)
C
#else
C
 2410 FORMAT
     & ('   Reading the previous time step number ',
     &                      '(restarting computation)  NTPABS =   ',I10)
 2411 FORMAT
     & ('   Reading the previous time step number ',
     &                      '(restarting computation)  TTPABS = ',E12.4)
C
#endif
C
C
C --- MISES EN GARDE
C
#if defined(_CS_LANG_FR)
C
 8000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@      Le nombre de phases    maximal NPHSMX supporte pas le ',/,
     &'@        format d''ecriture du fichier suite est             ',/,
     &'@        NFMTPH = ',I10                                       ,/,
     &'@      On a ici un nombre de phases    maximal superieur     ',/,
     &'@        NPHSMX = ',I10                                       ,/,
     &'@      Si le nombre de phases effectif est superieur, elles  ',/,
     &'@        ne seront pas relues.                               ',/,
     &'@                                                            ',/,
     &'@    Le calcul sera execute.                                 ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme lecamp.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@      Le nombre de scalaires maximal NSCAMX supporte pas le ',/,
     &'@        format d''ecriture du fichier suite est             ',/,
     &'@        NFMTSC = ',I10                                       ,/,
     &'@      On a ici un nombre de scalaires maximal superieur     ',/,
     &'@        NSCAMX = ',I10                                       ,/,
     &'@      On ne pourra pas relire les scalaires dont le numero  ',/,
     &'@        est superieur                                       ',/,
     &'@                                                            ',/,
     &'@    Le calcul sera execute.                                 ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme lecamp.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LECTURE DU FICHIER SUITE PRINCIPAL          ',/,
     &'@    *********                                               ',/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Le nombre de variables, de scalaires ou de phases a ete ',/,
     &'@      modifie.                                              ',/,
     &'@    Le calcul peut etre execute.                            ',/,
     &'@                                                            ',/,
     &'@    Il est cependant conseille de verifier                  ',/,
     &'@      les dimensions suivantes dans usini1 :                ',/,
     &'@                                                            ',/,
     &'@                NVAR     NSCAL    NSCAUS    NSCAPP     NPHAS',/,
     &'@  AMONT : ',5I10                                             ,/,
     &'@  ACTUEL: ',5I10                                             ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8410 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LECTURE DU FICHIER SUITE PRINCIPAL          ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      REPRISE  DE CALCUL           AVEC ITURB(',I3,') = ',I4 ,/,
     &'@      A PARTIR D''UN CALCUL REALISE AVEC ITURB(',I3,') = ',I4,/,
     &'@                                                            ',/,
     &'@    Le modele de turbulence en temps a ete modifie.         ',/,
     &'@    Le calcul peut etre execute.                            ',/,
     &'@                                                            ',/,
     &'@    Il est conseille cependant de                           ',/,
     &'@      verifier la valeur de ITURB(',I2,') dans usini1       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8711 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@    Probleme sur le fichier de nom (',A13,')                ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 8000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING :   WHEN READING THE MAIN RESTART FILE          ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      The maximum number of phases NPHSMX supported         ',/,
     &'@        by the restart file writing format is               ',/,
     &'@        NFMTPH = ',I10                                       ,/,
     &'@      Here is the maximum number of phases larger           ',/,
     &'@        NPHSMX = ',I10                                       ,/,
     &'@      If the actual number of phases is more than NFMTPH,   ',/,
     &'@        some of them might not be read.                     ',/,
     &'@                                                            ',/,
     &'@    The calculation will be run.                            ',/,
     &'@                                                            ',/,
     &'@    See the sub-routine lecamp.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING :  WHEN READING THE MAIN RESTART FILE        ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      The maximum number of scalars NSCAMX supported        ',/,
     &'@        by the restart file writing format is               ',/,
     &'@        NFMTSC = ',I10                                       ,/,
     &'@      Here is the maximum number of scalars larger          ',/,
     &'@        NSCAMX = ',I10                                       ,/,
     &'@      Scalars whose number is larger wont be read.          ',/,
     &'@                                                            ',/,
     &'@    The calculation will be run.                            ',/,
     &'@                                                            ',/,
     &'@    See the sub-routine lecamp                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : WHEN READING THE MAIN RESTARTING FILE       ',/,
     &'@    *********                                               ',/,
     &'@      THE RESTART AND CHECKPOINT DATA ARE DIFFERENT         ',/,
     &'@                                                            ',/,
     &'@    The number of variables, scalars or phases has changed. ',/,
     &'@    The calculation will be run.                            ',/,
     &'@                                                            ',/,
     &'@    However, it is strongly advised to check                ',/,
     &'@      the following dimensions in usini1 :                  ',/,
     &'@                                                            ',/,
     &'@                NVAR     NSCAL    NSCAUS    NSCAPP     NPHAS',/,
     &'@ PREVIOUS:',5I10                                             ,/,
     &'@ CURRENT :',5I10                                             ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8410 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : WHEN READING THE MAIN RESTART FILE            ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    THE CURRENT CALCULATION USES ITURB(',I3,') = ',I4        ,/,
     &'@    BUT RESTARTS FROM ANOTHER ONE USING ITURB(',I3,') = ',I4 ,/,
     &'@                                                            ',/,
     &'@    The turbulence model has changed.                       ',/,
     &'@    The computation can be executed.                        ',/,
     &'@                                                            ',/,
     &'@    However, it is strongly advised to check                ',/,
     &'@      the value of the variable ITURB(',I2,') in usini1     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8711 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : ERROR AT CLOSING THE MAIN RESTART FILE        ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    Problem on file called (',A13,')                        ',/,
     &'@                                                            ',/,
     &'@    The calculation will be run.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#endif
C
C --- ERREURS
C
#if defined(_CS_LANG_FR)
C
 9100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                               ',/,
     &'@      ERREUR A L''OUVERTURE DU FICHIER SUITE PRINCIPAL      ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier l''existence et le nom (',A13,') du            ',/,
     &'@        fichier suite dans le repertoire de travail.        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@      TYPE DE FICHIER INCORRECT                             ',/,
     &'@                                                            ',/,
     &'@    Le fichier ',A13      ,' ne semble pas etre un fichier  ',/,
     &'@      suite principal.                                      ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise correspond bien   ',/,
     &'@        a un fichier suite principal.                       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9201 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    Le nombre de cellules a ete modifie                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise correspond bien   ',/,
     &'@        au cas traite.                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@      ERREUR A LA LECTURE DES DIMENSIONS                    ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise n''a pas ete      ',/,
     &'@        endommage.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9320 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@    On souhaite faire correspondre le scalaire  ',I10        ,/,
     &'@            du present calcul avec le scalaire  ',I10        ,/,
     &'@            du calcul precedent, or, le nombre maximal de   ',/,
     &'@            scalaires dans le fichier suite est ',I10        ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier les valeurs de ISCOLD dans usini1.             ',/,
     &'@    Verifier que le fichier suite utilise correspond bien   ',/,
     &'@        au cas traite.                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9321 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@    On souhaite faire correspondre le scalaire  ',I10        ,/,
     &'@            du present calcul avec le scalaire  ',I10        ,/,
     &'@            du calcul precedent, or, le numero des anciens  ',/,
     &'@            scalaires doit etre strictement positif et      ',/,
     &'@            inferieur ou egal a                 ',I10        ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier les valeurs de ISCOLD dans usini1.             ',/,
     &'@    Verifier que le fichier suite utilise correspond bien   ',/,
     &'@        au cas traite.                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9330 FORMAT (
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@      CORRESPONDANCE AVEC LES ANCIENS SCALAIRES IMPOSSIBLE  ',/,
     &'@      ',I10,   ' ERREURS REPORTEES CI-DESSUS                ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@      Les lignes precedentes identifient les scalaires      ',/,
     &'@        definis dans le present calcul pour lesquels la     ',/,
     &'@        la correspondance fournie dans usini1 par ISCOLD    ',/,
     &'@        est incorrecte.                                     ',/,
     &'@                                                            ',/,
     &'@      Le calcul ne peut etre execute.                       ',/,
     &'@                                                            ',/,
     &'@      Verifier que le fichier suite UTILISE correspond bien ',/,
     &'@        au cas traite.                                      ',/,
     &'@      Verifier ISCOLD dans usini1.                          ',/,
     &'@        ISCOLD(ISCAL) = 0 indique que le scalaire ISCAL     ',/,
     &'@          du present calcul ne correspond a aucun scalaire  ',/,
     &'@          du calcul precedent                               ',/,
     &'@        ISCOLD(ISCAL) > 0 indique le numero du scalaire     ',/,
     &'@          du calcul precedent auquel correspond ISCAL       ',/,
     &'@      Si les correspondances avec les anciens scalaires     ',/,
     &'@        ne sont pas necessaires, eliminer ISCOLD de usini1. ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@      ERREUR A LA LECTURE DES INFORMATIONS TEMPORELLES OU   ',/,
     &'@        DES MODELES DE TURBULENCE                           ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise n''a pas ete      ',/,
     &'@        endommage.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9401 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A LA LECTURE DU FICHIER SUITE        ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@      ERREUR A LA LECTURE DE L''INDICATEUR DE METHODE ALE   ',/,
     &'@                                                            ',/,
     &'@    Il se peut que le fichier suite relu corresponde a une  ',/,
     &'@      version anterieure de Code_Saturne, sans methode ALE. ',/,
     &'@    Le calcul sera execute en reinitialisant toutes les     ',/,
     &'@      donnees ALE.                                          ',/,
     &'@    Verifier neanmoins que le fichier suite utilise n''a    ',/,
     &'@        pas ete endommage.                                  ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9402 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@      INDICATEUR IALE DU CALCUL PRECEDENT = ',I10            ,/,
     &'@      INDICATEUR IALE DU CALCUL ACTUEL    = ',I10            ,/,
     &'@                                                            ',/,
     &'@    Les coordonnees des noeuds du maillage doivent etre     ',/,
     &'@      relues. Elles sont stockees dans le fichier suite     ',/,
     &'@      auxiliaire.                                           ',/,
     &'@    L''indicateur ILEAUX doit donc etre positionne a 1.     ',/,
     &'@    Il vaut ici ILEAUX = ',I10                               ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@    Verifier ILEAUX dans usini1 ou dans l''Interface        ',/,
     &'@                                                  Graphique.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9410 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@      NUMERO DU PAS DE TEMPS PRECEDENT NTPABS = ',I10        ,/,
     &'@      NUMERO DU PAS DE TEMPS VISE      NTMABS = ',I10        ,/,
     &'@                                                            ',/,
     &'@    Le nombre de pas de temps (absolu) vise, NTMABS,        ',/,
     &'@      doit etre superieur au nombre de pas de temps         ',/,
     &'@      (absolu) deja effectues, NTPABS.                      ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier (augmenter) NTMABS dans usini1.                ',/,
     &'@    Verifier que le fichier suite utilise correspond bien   ',/,
     &'@        au cas traite.                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9510 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@      ERREUR LORS DE LA LECTURE DES VARIABLES PRINCIPALES   ',/,
     &'@        PRESSION, VITESSE, TURBULENCE                       ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise n''a pas ete      ',/,
     &'@        endommage.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C NFMTRU = 36 pour A36
 9511 FORMAT(
     &'@ Erreur a la lecture de ',A36                                  )
 9512 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@      ERREUR LORS DE LA LECTURE DES VARIABLES PRINCIPALES   ',/,
     &'@        SCALAIRES                                           ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise n''a pas ete      ',/,
     &'@        endommage.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9513 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                      PRINCIPAL',/,
     &'@                                                            ',/,
     &'@      ERREUR LORS DE LA LECTURE DE LA VITESSE DE MAILLAGE   ',/,
     &'@        (METHODE ALE)                                       ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise n''a pas ete      ',/,
     &'@        endommage.                                          ',/,
     &'@                                                            ',/,
     &'@    Le calcul peut neanmoins etre execute (vitesse de       ',/,
     &'@                                     maillage reinitialisee)',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 9100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE RESTART FILE READING              ',/,
     &'@    *********                                               ',/,
     &'@      ERROR AT OPENING THE MAIN RESTART FILE                ',/,
     &'@                                                            ',/,
     &'@    The computation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@    Please check the existence and the name (',A13,')       ',/,
     &'@        of the restart file in the working directory.       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@      WRONG FILE TYPE                                       ',/,
     &'@                                                            ',/,
     &'@    The file ',A13      ,' does not look like a proper      ',/,
     &'@      main restart file.                                    ',/,
     &'@                                                            ',/,
     &'@    The calculation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@    Please make sure the file used as a restart file        ',/,
     &'@        actually is a correct main restart file.            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9201 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@      INCONSISTANT RESTART AND CHECKPOINT DATA              ',/,
     &'@                                                            ',/,
     &'@    The number of cells has changed                         ',/,
     &'@                                                            ',/,
     &'@    The calculation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@    Please make sure the file used as restart file does     ',/,
     &'@        correspond to your case                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE RESTART FILE READING              ',/,
     &'@    *********                                               ',/,
     &'@      ERROR AT READING DIMENSIONS                           ',/,
     &'@                                                            ',/,
     &'@    The calculation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@    Please check the integrity of the file used as          ',/,
     &'@        restart file                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9320 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    The scalar ',I10   ,' of the current calculation        ',/,
     &'@        is attempted to match the scalar      ',I10          ,/,
     &'@        of the previous one. But the maximum number         ',/,
     &'@        of scalar read in the restart file is ',I10          ,/,
     &'@                                                            ',/,
     &'@    The calculation cannot be excuted.                      ',/,
     &'@                                                            ',/,
     &'@    Please check the value of ISCOLD in usini1.             ',/,
     &'@    Please make sure the file used as restart file does     ',/,
     &'@        correspond to your case                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9321 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    The scalar ',I10   ,' of the current calculation        ',/,
     &'@        is attempted to match the scalar ',I10               ,/,
     &'@        of the previous one. But the number                 ',/,
     &'@        of the previous scalar must be stricly larger       ',/,
     &'@        than zero and smaller than or equal to     ',I10     ,/,
     &'@                                                            ',/,
     &'@    The calculation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@    Please check the value of ISCOLD in usini1.             ',/,
     &'@    Please make sure the file used as restart file does     ',/,
     &'@        correspond to your case                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9330 FORMAT (
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      IMPOSSIBLE MATCHING WITH PREVIOUS SCALAR              ',/,
     &'@      ',I10,   ' ERRORS REPORTED IN THE ABOVE MESSAGES      ',/,
     &'@                                                            ',/,
     &'@    The calculation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@      The previous lines identify the scalars defined      ',/,
     &'@        in the current computation for which the matching   ',/,
     &'@        defined in usini1 by ISCOLD is not correct          ',/,
     &'@                                                            ',/,
     &'@      The computation cannot be executed.                   ',/,
     &'@                                                            ',/,
     &'@      Please make sure the file used as restart file does   ',/,
     &'@          correspond to your case                           ',/,
     &'@      Please check the value of ISCOLD in usini1.           ',/,
     &'@        ISCOLD(ISCAL) = 0 means that the scalar ISCAL       ',/,
     &'@          of the current calculation does not correspond to ',/,
     &'@          any of the scalars of the previous calculation    ',/,
     &'@        ISCOLD(ISCAL) > 0 defines the scalar number in the  ',/,
     &'@          previous calculation ISCAL corresponds to         ',/,
     &'@      If matching with the previous scalar is unnecessary   ',/,
     &'@        please delete ISCOLD from usini1.                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      ERROR AT READING THE TEMPORAL INFORMATION OR          ',/,
     &'@        TURBULENCE MODELS                                   ',/,
     &'@                                                            ',/,
     &'@    The computation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@    Please check the integrity of the file used as          ',/,
     &'@        restart file                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9401 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : ERROR AT THE MAIN RESTART FILE READING        ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      ERROR AT READING THE INDICATOR OF ALE METHOD          ',/,
     &'@                                                            ',/,
     &'@    The read restart file might come from a previous        ',/,
     &'@      version of Code Saturne, without ALE.                 ',/,
     &'@    The calculation will be executed but                    ',/,
     &'@      ALE data will be reset.                               ',/,
     &'@    Please check the integrity of the file used as          ',/,
     &'@        restart file, however.                              ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9402 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@      IALE INDICATOR OF THE PREVIOUS CALCULATION = ',I10     ,/,
     &'@      IALE INDICATOR OF THE CURRECT CALCULATION  = ',I10     ,/,
     &'@                                                            ',/,
     &'@    The coordinates of the mesh nodes need to be read again.',/,
     &'@      They are stored in the auxiliary restart file.        ',/,
     &'@    Therefore the ILEAUX indicator needs to be equal to 1.  ',/,
     &'@    Its current value is ILEAUX = ',I10                     ,/,
     &'@                                                            ',/,
     &'@    The calculation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@    Please check the value of ILEAUX in usini1 or in the    ',/,
     &'@                                   Graphical User Interface ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9410 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@      NUMBER OF THE PREVIOUS TIME STEP  NTPABS = ',I10       ,/,
     &'@      NUMBER OF TIME STEPS WANTED       NTMABS = ',I10       ,/,
     &'@                                                            ',/,
     &'@    The number of time steps (absolute) wanted, NTMABS,     ',/,
     &'@      has to be larger than the number of time steps        ',/,
     &'@      (absolute) already done, NTPABS.                      ',/,
     &'@                                                            ',/,
     &'@    The calculation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@    Please check (increase) NTMABS in usini1.               ',/,
     &'@    Please make sure the file used as restart file does     ',/,
     &'@          correspond to your case                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9510 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      ERROR AT THE READING OF MAIN VARIABLES                ',/,
     &'@        PRESSURE, VELOCITY, TURBULENCE                      ',/,
     &'@                                                            ',/,
     &'@    The calculation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@    Please check the integrity of the file used as          ',/,
     &'@        restart file                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9511 FORMAT(
     &'@ Error at the reading of ',A36                                 )
 9512 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      ERROR AT THE READING OF MAIN SCALAR VARIABLES         ',/,
     &'@                                                            ',/,
     &'@    The calculation cannot be executed.                     ',/,
     &'@                                                            ',/,
     &'@    Please check the integrity of the file used as          ',/,
     &'@        restart file                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9513 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      ERROR AT THE READING OF MESH VELOCITY                 ',/,
     &'@        (ALE METHOD)                                       ',/,
     &'@                                                            ',/,
     &'@    Please check the integrity of the file used as          ',/,
     &'@        restart file                                        ',/,
     &'@                                                            ',/,
     &'@    The calculation can be executed however (mesh           ',/,
     &'@                                     velocity reset)        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#endif
C
      END
c@z
