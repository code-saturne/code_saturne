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
                        SUBROUTINE LECAMX
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NNOD   ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDEVEL , ITUSER , IA     ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , FRCXT  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
c@foncb
CFONC
CFONC FONCTION :
CFONC ----------
CFONC  LECTURE DU FICHIER SUITE AUXILIAIRE
CFONC
CFONC  ON PASSE ICI SI ILEAUX = 1
CFONC  ON S'ARRETE SI ON NE PEUT PAS OUVRIR LE FICHIER
CFONC           OU SI CE N'EST PAS UN FICHIER AUXILIAIRE
CFONC           OU SI NCEL N'EST PAS CORRECT
CFONC           OU SI ON NE PEUT PAS LIRE JPHAS, JTURB, JDTVAR
CFONC           OU SI ON NE PEUT PAS LIRE UNE MOYENNE QU'ON VEUT POURSUIVRE
CFONC
CFONC
CFONC
CFONC
CFONC
CFONC
CFONC
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
CARGU ! NDIM         ! E  !  -> ! DIMENSION DU CALCUL                  !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE NOEUDS                     !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! JTURB        ! TE !  -> ! MODELES DE TURB CALCUL PRECEDENT     !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! DT(NCELET)   ! TR ! <-  ! PAS DE TEMPS                         !
CARGU ! RTP          ! TR ! <-  ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT        )!
CARGU ! PROPCE       ! TR ! <-  ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR ! <-  ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR ! <-  ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR ! <-  ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! FRCXT(NCELET,! TR ! <-  ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   3,NPHAS)   !    !     !  HYDROSTATIQUE                       !
CARGU ! RACELL(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! RAFACL(NFAC  ! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! RAFABL(NFABOR! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "albase.h"
      INCLUDE "alstru.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR, NNOD
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
C
      CHARACTER        RUBRIQ*64,CAR4*4,CAR2*2
      CHARACTER        CNUM4*4, CAR54*54
      CHARACTER        CINDFP*2,CINDFS*4,CINDFF*4,CINDFM*4
      CHARACTER        CINDFC*2,CINDFL*4
      CHARACTER        CPHASE(NPHSMX)*2
      CHARACTER        NOMFLU(NVARMX)*18,NOMCLI(NVARMX)*18
      CHARACTER        CSTRUC(NSTRMX)*2, CINDST*2
      INTEGER          IEL   , IFAC, II, ISTR
      INTEGER          IPHAS , IVAR  , ISCAL , JPHAS , ISCO
      INTEGER          IDECAL, ICLAPC, ICHA  , ICLA
      INTEGER          IMOM  , IMOLD
      INTEGER          JDTVAR
      INTEGER          JORTVM, IPCVMX, IPCVMY, IPCVMZ
      INTEGER          ICLVAR, ICLVAF
      INTEGER          IDTCM
      INTEGER          IPTSNA, IPTSTA, IPTSCA
      INTEGER          NUMERO, IPCEFJ, IPCLA1, IPCLA2, IPCLA3
      INTEGER          IOK   , INIFOK, IOKDT
      INTEGER          NCELOK, NFAIOK, NFABOK, NSOMOK
      INTEGER          IERROR, IRTYP,  ITYSUP, NBVAL
      INTEGER          NBERRO, INIERR, IVERS
      INTEGER          ILU   , ILECEC, IDEBLU, IANNUL, IERRCH
      INTEGER          IMPAMX
      INTEGER          NFMTPH, NFMTSC, NFMTFL, NFMTMO, NFMTCH, NFMTCL
      INTEGER          NFMTST
      INTEGER          NUMFLU(NVARMX)
      INTEGER          JTURB (NPHSMX), JTYTUR(NPHSMX), JALE
      DOUBLE PRECISION D2S3  , TSRII , CDTCM
      DOUBLE PRECISION TMPSTR(27)

C
C***********************************************************************
C
C=======================================================================
C 0. INITIALISATION
C=======================================================================
C
C  ---> Banniere
      WRITE(NFECRA,1000)
C
C  --->  On code en chaine le numero des phases et scalaires
C
C     Nombre max pour les formats choisis
      NFMTPH = 99
      NFMTSC = 9999
      NFMTFL = 9999
      NFMTMO = 9999
      NFMTCH = 99
      NFMTCL = 9999
C
C     Indefini a 2 et 4 caracteres
      CINDFP='YY'
      CINDFS='YYYY'
      CINDFF='YYYY'
      CINDFM='YYYY'
      CINDFC='YY'
      CINDFL='YYYY'
C
C     Codage en chaine de caracteres du numero de la phase
C       Aller jusqu'a NPHAS suffirait
      DO IPHAS = 1, MIN(NPHAS,NFMTPH)
        WRITE(CPHASE(IPHAS),'(I2.2)') IPHAS
      ENDDO
      DO IPHAS = MIN(NPHAS,NFMTPH)+1,NPHAS
       CPHASE(IPHAS) = CINDFP
      ENDDO
C
C
C     Avertissement
      IF(NPHSMX.GT.NFMTPH) THEN
        WRITE(NFECRA,8000)NFMTPH,NPHSMX
      ENDIF
      IF(NSCAMX.GT.NFMTSC) THEN
        WRITE(NFECRA,8001)NFMTSC,NSCAMX
      ENDIF
      IF(NVARMX.GT.NFMTFL) THEN
        WRITE(NFECRA,8002)NFMTFL,NVARMX
      ENDIF
      IF(NBMOMX.GT.NFMTMO) THEN
        WRITE(NFECRA,8003)NFMTMO,NBMOMX
      ENDIF
C
C=======================================================================
C 1. OUVERTURE DU FICHIER SUITE AUXILIAIRE
C=======================================================================
C
C  ---> Ouverture du fichier (ILECEC = 1 : lecture)
      ILECEC = 1
      CALL OPNSUI(FICAMX,LEN(FICAMX),ILECEC,IMPAMX,IERROR)
C     ===========
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,9000) FICAMX
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
C  Rubrique "fichier suite aux"
C        Pourrait porter le numero de version si besoin.
C        On ne se sert pas de IVERS pour le moment.
C
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      RUBRIQ = 'version_fichier_suite_auxiliaire'
      CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IVERS,
     &            IERROR)
C
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,9100)FICAMX
        CALL CSEXIT (1)
      ENDIF
C
C     Supports
C
      CALL TSTSUI(IMPAMX,NCELOK,NFAIOK,NFABOK,NSOMOK)
C     ===========
C
      IF (NCELOK.EQ.0) THEN
        WRITE(NFECRA,9101)
        CALL CSEXIT (1)
      ENDIF
C
      IF (NFAIOK.EQ.0) WRITE(NFECRA,8200)'internes','internes'
C
      IF (NFABOK.EQ.0) WRITE(NFECRA,8200)'de bord ','de bord '
C
C
C     On n'a besoin que
C       du nombre de phases et du modele de turbulence
C
      NBERRO = 0
C
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
C
      RUBRIQ = 'nombre_phases'
      CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,JPHAS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
C  ---> On s'arrete si erreur
C     Si on ne peut pas relire un entier, c'est que le fichier n'est pas bon
      IF (NBERRO.NE.0) THEN
        CAR54 ='ERREUR A lA LECTURE DES DIMENSIONS (NPHAS)            '
        WRITE(NFECRA,9200)CAR54
        CALL CSEXIT (1)
      ENDIF
C
C  ---> On previent si des parametres sont differents
      IF ( JPHAS .NE.NPHAS ) THEN
        WRITE(NFECRA,8210) JPHAS, NPHAS
      ENDIF
C
C     Modele de turbulence
C
      NBERRO = 0
C
      DO IPHAS = 1, MIN(NPHAS,JPHAS)
        RUBRIQ = 'modele_turbulence_phase'//CPHASE(IPHAS)
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              JTURB(IPHAS),IERROR)
        NBERRO=NBERRO+IERROR
        JTYTUR(IPHAS)=JTURB(IPHAS)/10
      ENDDO
C
C --->  Stop si erreur
C     Si on ne peut pas relire un entier, c'est que le fichier n'est pas bon
      IF (NBERRO.NE.0) THEN
        CAR54 ='ERREUR A lA LECTURE DES MODELES DE TURBULENCE         '
        WRITE(NFECRA,9200)CAR54
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
      CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,JALE,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
C --->  Message si erreur (pas de stop pour compatibilite avec les fichiers anterieurs)
C       -> on n'affiche le message que si IALE=1 (sinon RAS)
      IF (NBERRO.NE.0) THEN
        IF (IALE.EQ.1) WRITE(NFECRA,9210)
        JALE = 0
      ENDIF
C
C ---> Pas d'iteration d'initialisation si suite de calcul ALE
      IF (ITALIN.EQ.-999) THEN
        IF (IALE.EQ.1 .AND. JALE.EQ.1) THEN
          ITALIN = 0
        ELSE IF (IALE.EQ.1) THEN
          ITALIN = 1
        ELSE
          ITALIN = 0
        ENDIF
      ENDIF
C
C --->  Donnees modifiees
      DO IPHAS = 1, MIN(NPHAS,JPHAS)
        IF (ITURB(IPHAS) .NE. JTURB(IPHAS))
     &       WRITE(NFECRA,8220) IPHAS,ITURB(IPHAS),
     &                          IPHAS,JTURB(IPHAS),
     &                          IPHAS
      ENDDO
C
      CAR54 =' Fin de la lecture des options                        '
      WRITE(NFECRA,1110)CAR54
C
C=======================================================================
C 3. PROPRIETES PHYSIQUES
C=======================================================================
C
      NBERRO = 0
C
C     On lit les infos des phases communes
C     Pour les phases supplementaires, les initialisations par defaut
C     ont ete faites dans INIVA0
C
      DO IPHAS = 1, MIN(NPHAS,JPHAS)
C
C ---> Point de reference de pression
C     On lit les coordonnees si l'utilisateur n'a rien specifie, i.e.
C       si IXYZP0=-1, et on met IXYZP0 a 1 pour eviter de le changer ensuite.
        IF (IXYZP0(IPHAS).EQ.-1) THEN
          RUBRIQ = 'ref_presstot'//CPHASE(IPHAS)
          ITYSUP = 0
          NBVAL  = 3
          IRTYP  = 2
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                XYZP0(1,IPHAS),IERROR)
          NBERRO = NBERRO+IERROR
          IF (IERROR.EQ.0) THEN
            WRITE(NFECRA,7000) IPHAS,(XYZP0(II,IPHAS),II=1,3)
            IXYZP0(IPHAS) = 1
          ENDIF
        ENDIF
C
C ---> Masse volumique
C     On la lit, qu'elle soit extrapolee ou pas,
C       pour permettre les sous-relaxations
C     Pour les suites a rho constant, cependant, on ne la lit pas,
C       afin de ne pas ecraser la valeur RO0 fixee par l'utilisateur
C       et eventuellement modifiee.
C
        INIERR = 0
C
        IF(IROVAR(IPHAS).EQ.1) THEN
C       Masse volumique - cellules
          RUBRIQ = 'rho_ce_phase'//CPHASE(IPHAS)
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPPROC(IROM(IPHAS))),IERROR)
          NBERRO = NBERRO+IERROR
          INIERR = INIERR+IERROR
C
C       Masse volumique - faces de bord
          IF (NFABOK.EQ.1) THEN
            RUBRIQ = 'rho_fb_phase'//CPHASE(IPHAS)
            ITYSUP = 3
            NBVAL  = 1
            IRTYP  = 2
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPFB(1,IPPROB(IROM(IPHAS))),IERROR)
            NBERRO = NBERRO+IERROR
            INIERR = INIERR+IERROR
          ENDIF
C
C     Si on a reussi a initialiser la masse volumique aux cellules ET
C       aux faces de bord, on l'indique (pour schtmp)
          IF (NFABOK.EQ.1.AND.INIERR.EQ.0) THEN
            INITRO(IPHAS) = 1
          ENDIF
C
        ELSE
C     Si la masse volumique est constante, elle est initialisee
C       correctement
          INITRO(IPHAS) = 1
        ENDIF
C
C ---> Viscosite moleculaire et "turbulente" ou de "sous-maille"
C     Si elle est extrapolee en temps, on la lit
C     Si on reussit, on l'indique
C
C     On cherche a lire uniquement si on doit extrapoler en temps
        IF(IVIEXT(IPHAS).GT.0) THEN
C
          INIERR = 0
C
C         Viscosite moleculaire - cellules
C         Uniquement si elle est variable
          IF(IVIVAR(IPHAS).EQ.1) THEN
            RUBRIQ = 'viscl_ce_phase'//CPHASE(IPHAS)
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPPROC(IVISCL(IPHAS))),IERROR)
            NBERRO = NBERRO+IERROR
            INIERR = INIERR+IERROR
          ENDIF
C
C         Viscosite turbulente ou de sous-maille - cellules
          RUBRIQ = 'visct_ce_phase'//CPHASE(IPHAS)
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPPROC(IVISCT(IPHAS))),IERROR)
          NBERRO = NBERRO+IERROR
          INIERR = INIERR+IERROR
C
C     Si on a initialise les viscosites, on l'indique (pour schtmp)
          IF (INIERR.EQ.0) THEN
            INITVI(IPHAS) = 1
          ENDIF
C
        ENDIF
C
C
C ---> Chaleur massique
C     On cherche a la lire si elle est variable
C       et qu'on l'extrapole ou qu'on est en Joule
C     Si on reussit, on l'indique
C       (Ca sert quand elle est extrapolee en temps
C        et quand l'utilisateur peut s'en servir pour passer
C        de H a T, comme en effet Joule par exemple).
C
        IF((ICPEXT(IPHAS ).GT.0.AND.ICP(IPHAS).GT.0).OR.
     &     (IPPMOD(IELJOU).GE.1.AND.ICP(IPHAS).GT.0)) THEN
C
          INIERR = 0
C
C         Chaleur massique - cellules
          RUBRIQ = 'cp_ce_phase'//CPHASE(IPHAS)
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPPROC(ICP(IPHAS))),IERROR)
          NBERRO = NBERRO+IERROR
          INIERR = INIERR+IERROR
C
C     Si on a initialise Cp, on l'indique (pour schtmp)
          IF (INIERR.EQ.0) THEN
            INITCP(IPHAS) = 1
          ENDIF
C
        ENDIF
C
      ENDDO
C
C
C     Si on a des scalaires, on lit a diffusivite
C       si le scalaire a  un correspondant et
C       si on doit extrapoler la diffusivite
C       (et qu'elle est variable, et que le scalaire n'est pas une variance)
C
      IF(NSCAL.GT.0) THEN
C
        DO ISCAL = 1, NSCAL
          ISCO = ISCOLD(ISCAL)
          IF ( ISCO         .GT.0.AND.
     &         IVSEXT(ISCAL).GT.0.AND.
     &         IVISLS(ISCAL).GT.0.AND.
     &        (ISCAVR(ISCAL).LE.0.OR.ISCAVR(ISCAL).GT.NSCAL) ) THEN
C
            INIERR = 0
C
C         Diffusivite - cellules
            IF(ISCO.LE.NFMTSC) THEN
              WRITE(CAR4,'(I4.4)')ISCO
              RUBRIQ = 'visls_ce_scalaire'//CAR4
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPPROC(IVISLS(ISCAL))),IERROR)
              NBERRO = NBERRO+IERROR
              INIERR = INIERR+IERROR
            ELSE
              NBERRO = NBERRO-1
              INIERR = INIERR-1
            ENDIF
C
C     Si on a initialise visls, on l'indique (pour schtmp)
            IF (INIERR.EQ.0) THEN
              INITVS(ISCAL) = 1
            ENDIF
          ENDIF
        ENDDO
C
C     Pour les variances, il suffit de dire qu'on a initialise ou non
C       selon ce qu'on a fait au scalaire correspondant
        DO ISCAL = 1, NSCAL
          IF(ISCAVR(ISCAL).GT.0.AND.ISCAVR(ISCAL).LE.NSCAL) THEN
            INITVS(ISCAL) = INITVS(ISCAVR(ISCAL))
          ENDIF
        ENDDO
C
      ENDIF
C
C     Si erreur, on previent mais pas stop :
C       auparavant on n'avait pas stocke les prop phy
C         si on n'etait pas a l'ordre 2
C       c'est discutable pour rho
C
      IF (NBERRO.NE.0) THEN
        CAR54 = 'LECTURE DES PROPRIETES PHYSIQUES                      '
        WRITE(NFECRA,8300)CAR54
      ENDIF
C
      CAR54 = ' Fin de la lecture des proprietes physiques           '
      WRITE(NFECRA,1110)CAR54
C
C
C=======================================================================
C 4. PAS DE TEMPS
C=======================================================================
C
C
C  ---> Indicateur de pas de temps variable
      RUBRIQ = 'indic_dt_variable'
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,JDTVAR,
     &            IERROR)
C
C  ---> On s'arrete si erreur
C     Si on ne peut pas relire un entier, c'est que le fichier n'est pas bon
      IF (IERROR.NE.0) THEN
        CAR54 ='ERREUR A LA LECTURE DU MODE DE MARCHE EN TEMPS        '
        WRITE(NFECRA,9200)CAR54
        CALL CSEXIT(1)
      ENDIF
C
C  ---> Pas de temps
C     Rq : si marche en temps differente, on conserve la valeur par defaut
C          DTREF imposee dans INIVA0
C
      NBERRO = 0
C
      IF (IDTVAR.NE.JDTVAR) THEN
        WRITE(NFECRA,8400)IDTVAR,JDTVAR,DTREF
C
      ELSEIF (IDTVAR.EQ.1) THEN
        RUBRIQ = 'dt_variable_temps'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,DT(1),
     &              IERROR)
        NBERRO=NBERRO+IERROR
        DO IEL = 1, NCEL
          DT(IEL) = DT(1)
        ENDDO
C
      ELSEIF (IDTVAR.EQ.2) THEN
        RUBRIQ = 'dt_variable_espace_ce'
        ITYSUP = 1
        NBVAL  = 1
        IRTYP  = 2
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,DT,
     &              IERROR)
        NBERRO=NBERRO+IERROR
      ENDIF
C
C     Si erreur, on previent mais pas stop :
C       le pas de temps n'est pas une donnee majeure
C       c'est discutable
C
      IF (NBERRO.NE.0) THEN
        CAR54 = 'LECTURE DU PAS DE TEMPS                               '
        WRITE(NFECRA,8300)CAR54
      ENDIF
C
      CAR54 = ' Fin de la lecture du pas de temps                    '
      WRITE(NFECRA,1110)CAR54
C
C=======================================================================
C 5. FLUX DE MASSE
C=======================================================================
C     Pour retrouver la correspondance entre les variables
C     et les flux de masse, on utilise le nom de chaque variable
C     (NOMFLU(I)= nom de la ieme variable)
C     Ensuite, pour chaque variable, si on a deja le flux, on ne fait
C       rien, sinon on lit quel est le numero local du
C       flux qui lui est associe (en pratique 1 ou 2) et le flux lui meme.
C
C     Les flux de masse ne sont a lire que si les supports des faces
C     de bord ou des faces internes coincident
C
C     On lit d'abord le flux de masse (a l'instant n) et
C       ensuite a l'instant precedent si on est en schema en temps
C       particulier (ISTMPF NE 1)
C
      IF (NFAIOK.EQ.1 .OR. NFABOK.EQ.1) THEN
C
        NBERRO=0
C
C        Initialisation des tableaux de travail
        DO IVAR = 1, NVARMX
          NUMFLU(IVAR) = 0
        ENDDO
C
C     Nom du flux associe a la variable dans le calcul precedent
        NOMFLU(IPR(1))='fm_p_phase'//CPHASE(1)
        DO IPHAS = 1, MIN(NPHAS,JPHAS)
          NOMFLU(IU(IPHAS))='fm_u_phase'//CPHASE(IPHAS)
          NOMFLU(IV(IPHAS))='fm_v_phase'//CPHASE(IPHAS)
          NOMFLU(IW(IPHAS))='fm_w_phase'//CPHASE(IPHAS)
          IF (ITYTUR(IPHAS).EQ.2.AND.
     &         (JTYTUR(IPHAS).EQ.2.OR.JTURB(IPHAS).EQ.50) ) THEN
            NOMFLU(IK(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
          ELSEIF(ITYTUR(IPHAS).EQ.2.AND.JTYTUR(IPHAS).EQ.3) THEN
            NOMFLU(IK(IPHAS))='fm_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
          ELSEIF(ITYTUR(IPHAS).EQ.2.AND.JTURB(IPHAS).EQ.60) THEN
            NOMFLU(IK(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_omega_phase'//CPHASE(IPHAS)
          ELSEIF (ITYTUR(IPHAS).EQ.3.AND.
     &         (JTYTUR(IPHAS).EQ.2.OR.JTURB(IPHAS).EQ.50) ) THEN
            NOMFLU(IR11(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR22(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR33(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR12(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR13(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR23(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITYTUR(IPHAS).EQ.3.AND.JTYTUR(IPHAS).EQ.3) THEN
            NOMFLU(IR11(IPHAS))='fm_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IR22(IPHAS))='fm_R22_phase'//CPHASE(IPHAS)
            NOMFLU(IR33(IPHAS))='fm_R33_phase'//CPHASE(IPHAS)
            NOMFLU(IR12(IPHAS))='fm_R12_phase'//CPHASE(IPHAS)
            NOMFLU(IR13(IPHAS))='fm_R13_phase'//CPHASE(IPHAS)
            NOMFLU(IR23(IPHAS))='fm_R23_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITYTUR(IPHAS).EQ.3.AND.JTURB(IPHAS).EQ.60) THEN
            NOMFLU(IR11(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR22(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR33(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR12(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR13(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR23(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_omega_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50.AND.JTYTUR(IPHAS).EQ.2) THEN
            NOMFLU(IK(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
            NOMFLU(IPHI(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IFB(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50.AND.JTYTUR(IPHAS).EQ.3) THEN
            NOMFLU(IK(IPHAS))='fm_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
            NOMFLU(IPHI(IPHAS))='fm_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IFB(IPHAS))='fm_R11_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50.AND.JTURB(IPHAS).EQ.50) THEN
            NOMFLU(IK(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
            NOMFLU(IPHI(IPHAS))='fm_phi_phase'//CPHASE(IPHAS)
            NOMFLU(IFB(IPHAS))='fm_fb_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50.AND.JTURB(IPHAS).EQ.60) THEN
            NOMFLU(IK(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_omega_phase'//CPHASE(IPHAS)
            NOMFLU(IPHI(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IFB(IPHAS))='fm_k_phase'//CPHASE(IPHAS)

          ELSEIF (ITURB(IPHAS).EQ.60.AND.
     &         (JTYTUR(IPHAS).EQ.2.OR.JTURB(IPHAS).EQ.50) ) THEN
            NOMFLU(IK  (IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IOMG(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.60.AND.JTYTUR(IPHAS).EQ.3) THEN
            NOMFLU(IK  (IPHAS))='fm_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IOMG(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.60.AND.JTURB(IPHAS).EQ.60) THEN
            NOMFLU(IK  (IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IOMG(IPHAS))='fm_omega_phase'//CPHASE(IPHAS)
          ENDIF
        ENDDO
        IF(NSCAL.GT.0) THEN
          DO ISCAL = 1, NSCAL
            IF(ISCOLD(ISCAL).GT.0) THEN
              IF(ISCOLD(ISCAL).LE.NFMTSC) THEN
                WRITE(CAR4,'(I4.4)')ISCOLD(ISCAL)
              ELSE
                CAR4 = CINDFS
              ENDIF
              NOMFLU(ISCA(ISCAL))='fm_scalaire'//CAR4
            ENDIF
          ENDDO
        ENDIF
        IF (IALE.EQ.1) THEN
          NOMFLU(IUMA)='fm_vit_maill_u'
          NOMFLU(IVMA)='fm_vit_maill_v'
          NOMFLU(IWMA)='fm_vit_maill_w'
        ENDIF
C
C     --Pour les variables
        DO IVAR = 1, NVAR
C        S'il y a un flux
          IF(IFLUMA(IVAR).GT.0) THEN
C          Si le flux n'a pas deja ete rempli, on le fait
            IF (NUMFLU(IFLUMA(IVAR)).EQ.0) THEN
C            Lecture du numero local du flux correspondant
              RUBRIQ = NOMFLU(IVAR)
              ITYSUP = 0
              NBVAL  = 1
              IRTYP  = 1
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    NUMERO,IERROR)
              NBERRO = NBERRO+IERROR
C            S'il existe, on le lit
              IF (IERROR.EQ.0) THEN
                INIFOK = 0
                IF(NUMERO.GT.0) THEN
                  IF(NUMERO.LE.NFMTFL) THEN
                    WRITE(CAR4,'(I4.4)')NUMERO
                  ELSE
                    CAR4 = CINDFF
                  ENDIF
                  IF (NFAIOK.EQ.1) THEN
                    RUBRIQ = 'flux_masse_fi_'//CAR4
                    ITYSUP = 2
                    NBVAL  = 1
                    IRTYP  = 2
                    CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &                          IRTYP,PROPFA(1,IPPROF(IFLUMA(IVAR))),
     &                          IERROR)
                    NBERRO = NBERRO+IERROR
                    IF(IERROR.EQ.0) INIFOK = INIFOK+1
                  ENDIF
                  IF (NFABOK.EQ.1) THEN
                    RUBRIQ = 'flux_masse_fb_'//CAR4
                    ITYSUP = 3
                    NBVAL  = 1
                    IRTYP  = 2
                    CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &                          IRTYP,PROPFB(1,IPPROB(IFLUMA(IVAR))),
     &                          IERROR)
                    NBERRO = NBERRO+IERROR
                    IF(IERROR.EQ.0) INIFOK = INIFOK+1
                  ENDIF
                ENDIF
C             Si ca a marche, on indique que ce flux a ete rempli
                IF(INIFOK.EQ.2) THEN
                  NUMFLU(IFLUMA(IVAR))=1
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C
C
C     ET MAINTENANT : pareil, mais FLUMAS au pas de temps precedent
C
C        Initialisation des tableaux de travail
        DO IVAR = 1, NVARMX
          NUMFLU(IVAR) = 0
        ENDDO
C
C     Nom du flux associe a la variable dans le calcul precedent
        NOMFLU(IPR(1))='fm_a_p_phase'//CPHASE(1)
        DO IPHAS = 1, MIN(NPHAS,JPHAS)
          NOMFLU(IU(IPHAS))='fm_a_u_phase'//CPHASE(IPHAS)
          NOMFLU(IV(IPHAS))='fm_a_v_phase'//CPHASE(IPHAS)
          NOMFLU(IW(IPHAS))='fm_a_w_phase'//CPHASE(IPHAS)
          IF (ITYTUR(IPHAS).EQ.2.AND.
     &         (JTYTUR(IPHAS).EQ.2.OR.JTURB(IPHAS).EQ.50) ) THEN
            NOMFLU(IK(IPHAS))='fm_a_k_a_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
          ELSEIF(ITYTUR(IPHAS).EQ.2.AND.JTYTUR(IPHAS).EQ.3) THEN
            NOMFLU(IK(IPHAS))='fm_a_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
          ELSEIF(ITYTUR(IPHAS).EQ.2.AND.JTURB(IPHAS).EQ.60) THEN
            NOMFLU(IK(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_omega_phase'//CPHASE(IPHAS)
          ELSEIF (ITYTUR(IPHAS).EQ.3.AND.
     &         (JTYTUR(IPHAS).EQ.2.OR.JTURB(IPHAS).EQ.50) ) THEN
            NOMFLU(IR11(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR22(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR33(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR12(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR13(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR23(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITYTUR(IPHAS).EQ.3.AND.JTYTUR(IPHAS).EQ.3) THEN
            NOMFLU(IR11(IPHAS))='fm_a_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IR22(IPHAS))='fm_a_R22_phase'//CPHASE(IPHAS)
            NOMFLU(IR33(IPHAS))='fm_a_R33_phase'//CPHASE(IPHAS)
            NOMFLU(IR12(IPHAS))='fm_a_R12_phase'//CPHASE(IPHAS)
            NOMFLU(IR13(IPHAS))='fm_a_R13_phase'//CPHASE(IPHAS)
            NOMFLU(IR23(IPHAS))='fm_a_R23_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITYTUR(IPHAS).EQ.3.AND.JTURB(IPHAS).EQ.60) THEN
            NOMFLU(IR11(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR22(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR33(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR12(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR13(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IR23(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_omega_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50.AND.JTYTUR(IPHAS).EQ.2) THEN
            NOMFLU(IK(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
            NOMFLU(IPHI(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IFB(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50.AND.JTYTUR(IPHAS).EQ.3) THEN
            NOMFLU(IK(IPHAS))='fm_a_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
            NOMFLU(IPHI(IPHAS))='fm_a_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IFB(IPHAS))='fm_a_R11_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50.AND.JTURB(IPHAS).EQ.50) THEN
            NOMFLU(IK(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
            NOMFLU(IPHI(IPHAS))='fm_a_phi_phase'//CPHASE(IPHAS)
            NOMFLU(IFB(IPHAS))='fm_a_fb_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50.AND.JTURB(IPHAS).EQ.60) THEN
            NOMFLU(IK(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_omega_phase'//CPHASE(IPHAS)
            NOMFLU(IPHI(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IFB(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.60.AND.
     &         (JTYTUR(IPHAS).EQ.2.OR.JTURB(IPHAS).EQ.50) ) THEN
            NOMFLU(IK  (IPHAS))='fm_a_k_a_phase'//CPHASE(IPHAS)
            NOMFLU(IOMG(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
          ELSEIF(ITURB(IPHAS).EQ.60.AND.JTYTUR(IPHAS).EQ.3) THEN
            NOMFLU(IK  (IPHAS))='fm_a_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IOMG(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
          ELSEIF(ITURB(IPHAS).EQ.60.AND.JTURB(IPHAS).EQ.60) THEN
            NOMFLU(IK  (IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IOMG(IPHAS))='fm_a_omega_phase'//CPHASE(IPHAS)
          ENDIF
        ENDDO
        IF(NSCAL.GT.0) THEN
          DO ISCAL = 1, NSCAL
            IF(ISCOLD(ISCAL).GT.0) THEN
              IF(ISCOLD(ISCAL).LE.NFMTSC) THEN
                WRITE(CAR4,'(I4.4)')ISCOLD(ISCAL)
              ELSE
                CAR4 = CINDFS
              ENDIF
              NOMFLU(ISCA(ISCAL))='fm_a_scalaire'//CAR4
            ENDIF
          ENDDO
        ENDIF
        IF (IALE.EQ.1) THEN
          NOMFLU(IUMA)='fm_a_vit_maill_u'
          NOMFLU(IVMA)='fm_a_vit_maill_v'
          NOMFLU(IWMA)='fm_a_vit_maill_w'
        ENDIF
C
C     --Pour les variables
        DO IVAR = 1, NVAR
C        S'il y a un flux
          IF(IFLUAA(IVAR).GT.0) THEN
C          Si le flux n'a pas deja ete rempli, on le fait
            IF (NUMFLU(IFLUAA(IVAR)).EQ.0) THEN
C            Lecture du numero local du flux correspondant
              RUBRIQ = NOMFLU(IVAR)
              ITYSUP = 0
              NBVAL  = 1
              IRTYP  = 1
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    NUMERO,IERROR)
              NBERRO = NBERRO+IERROR
C            S'il existe, on le lit
              IF (IERROR.EQ.0) THEN
                INIFOK = 0
                IF(NUMERO.GT.0) THEN
                  IF(NUMERO.LE.NFMTFL) THEN
                    WRITE(CAR4,'(I4.4)')NUMERO
                  ELSE
                    CAR4 = CINDFF
                  ENDIF
                  IF (NFAIOK.EQ.1) THEN
                    RUBRIQ = 'flux_masse_a_fi_'//CAR4
                    ITYSUP = 2
                    NBVAL  = 1
                    IRTYP  = 2
                    CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &                          IRTYP,PROPFA(1,IPPROF(IFLUAA(IVAR))),
     &                          IERROR)
                    NBERRO = NBERRO+IERROR
                    IF(IERROR.EQ.0) INIFOK = INIFOK+1
                  ENDIF
                  IF (NFABOK.EQ.1) THEN
                    RUBRIQ = 'flux_masse_a_fb_'//CAR4
                    ITYSUP = 3
                    NBVAL  = 1
                    IRTYP  = 2
                    CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &                          IRTYP,PROPFB(1,IPPROB(IFLUAA(IVAR))),
     &                          IERROR)
                    NBERRO = NBERRO+IERROR
                    IF(IERROR.EQ.0) INIFOK = INIFOK+1
                  ENDIF
                ENDIF
C              Si ca a marche, on indique que ce flux a ete rempli
                IF(INIFOK.EQ.2) THEN
                  NUMFLU(IFLUAA(IVAR))=1
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C
C
C     Si erreur, on previent mais pas stop :
        IF (NBERRO.NE.0) THEN
          CAR54 =
     &         'LECTURE DES FLUX DE MASSE                             '
          WRITE(NFECRA,8300)CAR54
        ENDIF
C
      CAR54 = ' Fin de la lecture des flux de masse                  '
      WRITE(NFECRA,1110)CAR54
C
      ENDIF
C     fin de "s'il faut lire les flux de masse (ie. supports coincidents)"
C
C=======================================================================
C 6. CONDITIONS AUX LIMITES
C=======================================================================
C
C     A ne relire que si les supports sont identiques (faces de bord)
C
      ILU = 0
C
      IF (NFABOK.EQ.1) THEN
C
        ILU = 1
C
        NBERRO=0
C
C     Nom de variable associe a la variable dans le calcul precedent
        NOMCLI(IPR(1))='_p_phase'//CPHASE(1)
        DO IPHAS = 1, MIN(NPHAS,JPHAS)
          NOMCLI(IU(IPHAS))='_u_phase'//CPHASE(IPHAS)
          NOMCLI(IV(IPHAS))='_v_phase'//CPHASE(IPHAS)
          NOMCLI(IW(IPHAS))='_w_phase'//CPHASE(IPHAS)
C     Pour un calcul k-eps, on peut recuperer les CL d'un calcul k-eps
C     ou d'un calcul v2f
          IF (ITYTUR(IPHAS).EQ.2.AND.
     &         (JTYTUR(IPHAS).EQ.2.OR.JTURB(IPHAS).EQ.50)) THEN
            NOMCLI(IK(IPHAS))='_k_phase'//CPHASE(IPHAS)
            NOMCLI(IEP(IPHAS))='_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITYTUR(IPHAS).EQ.3.AND.JTYTUR(IPHAS).EQ.3) THEN
            NOMCLI(IR11(IPHAS))='_R11_phase'//CPHASE(IPHAS)
            NOMCLI(IR22(IPHAS))='_R22_phase'//CPHASE(IPHAS)
            NOMCLI(IR33(IPHAS))='_R33_phase'//CPHASE(IPHAS)
            NOMCLI(IR12(IPHAS))='_R12_phase'//CPHASE(IPHAS)
            NOMCLI(IR13(IPHAS))='_R13_phase'//CPHASE(IPHAS)
            NOMCLI(IR23(IPHAS))='_R23_phase'//CPHASE(IPHAS)
            NOMCLI(IEP(IPHAS))='_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50.AND.JTURB(IPHAS).EQ.50) THEN
            NOMCLI(IK(IPHAS))='_k_phase'//CPHASE(IPHAS)
            NOMCLI(IEP(IPHAS))='_eps_phase'//CPHASE(IPHAS)
            NOMCLI(IPHI(IPHAS))='_phi_phase'//CPHASE(IPHAS)
            NOMCLI(IFB(IPHAS))='_fb_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.60.AND.JTURB(IPHAS).EQ.60) THEN
            NOMCLI(IK(IPHAS))='_k_phase'//CPHASE(IPHAS)
            NOMCLI(IOMG(IPHAS))='_omega_phase'//CPHASE(IPHAS)
C     On peut aussi recuperer les CL de k et eps pour un calcul v2f suite
C     d'un calcul k-eps
          ELSEIF (ITURB(IPHAS).EQ.50.AND.JTYTUR(IPHAS).EQ.2) THEN
            NOMCLI(IK(IPHAS))='_k_phase'//CPHASE(IPHAS)
            NOMCLI(IEP(IPHAS))='_eps_phase'//CPHASE(IPHAS)
C     On peut aussi recuperer les CL de k pour un calcul k-eps, v2f ou k-omega suite
C     d'un calcul k-eps, v2f ou k-omega et qui n'est pas deja un des cas ci-dessus.
          ELSEIF ( (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &           .OR. ITURB(IPHAS).EQ.60)  .AND. (JTYTUR(IPHAS).EQ.2
     &           .OR. JTURB(IPHAS).EQ.50 .OR. JTURB(IPHAS).EQ.60) ) THEN
            NOMCLI(IK(IPHAS))='_k_phase'//CPHASE(IPHAS)
          ENDIF
        ENDDO
        IF(NSCAL.GT.0) THEN
          DO ISCAL = 1, NSCAL
            IF(ISCOLD(ISCAL).GT.0) THEN
              IF(ISCOLD(ISCAL).LE.NFMTSC) THEN
                WRITE(CAR4,'(I4.4)')ISCOLD(ISCAL)
              ELSE
                CAR4 = CINDFS
              ENDIF
              NOMCLI(ISCA(ISCAL))='_scalaire'//CAR4
            ENDIF
          ENDDO
        ENDIF
        IF (IALE.EQ.1) THEN
          NOMCLI(IUMA)='_vit_maillage_u'
          NOMCLI(IVMA)='_vit_maillage_v'
          NOMCLI(IWMA)='_vit_maillage_w'
        ENDIF
C
C     --Pour les variables
        DO IVAR = 1, NVAR
C
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 2
C
C              Coefficients numeros 1
          ICLVAR = ICLRTP(IVAR,ICOEF)
          RUBRIQ = 'cla1'//NOMCLI(IVAR)
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                COEFA(1,ICLVAR),IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'clb1'//NOMCLI(IVAR)
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                COEFB(1,ICLVAR),IERROR)
          NBERRO=NBERRO+IERROR
C
C              Coefficients numeros 2
          ICLVAF = ICLRTP(IVAR,ICOEFF)
          IF (ICLVAR.NE.ICLVAF) THEN
C
            RUBRIQ = 'cla2'//NOMCLI(IVAR)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &                  IRTYP,COEFA(1,ICLVAF),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'clb2'//NOMCLI(IVAR)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &                  IRTYP,COEFB(1,ICLVAF),IERROR)
            NBERRO=NBERRO+IERROR
          ENDIF
C
        ENDDO
C
C     Type symétrie (utilisé pour les gradients par moindres carrés
C       sur support étendu, avec extrapolation du gradient au bord).
C
        DO IPHAS = 1, MIN(NPHAS,JPHAS)
C
          RUBRIQ = 'isympa_fb_phase'//CPHASE(IPHAS)
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IISYMP+NFABOR*(IPHAS-1)),IERROR)
          NBERRO = NBERRO+IERROR
C
        ENDDO
C
C
      ENDIF
C     fin du test "si les supports des faces de bord sont identiques"
C
      IF (ILU.EQ.1) THEN
C
C     Si erreur, on previent mais pas stop :
C       (on n'a pas forcement les coefs 2, on n'a pas forcément isympa
C        si on prend les fichiers d'une version antérieure)
        IF (NBERRO.NE.0) THEN
          CAR54 =
     &         'LECTURE DES CONDITIONS AUX LIMITES                    '
          WRITE(NFECRA,8300)CAR54
        ENDIF
C
        CAR54 = ' Fin de la lecture des conditions aux limites         '
        WRITE(NFECRA,1110)CAR54
C
      ENDIF
C
C
C=======================================================================
C 7. TERMES SOURCES EXTRAPOLES
C=======================================================================
C
C
      NBERRO=0
      ILU = 0
C
C ---> Termes sources Navier-Stokes
C
C     Boucle sur les phases communes
      DO IPHAS = 1, MIN(JPHAS,NPHAS)
C     Si le terme est a l'ordre 2
        IF(ISNO2T(IPHAS).GT.0) THEN
C
          IPTSNA = IPPROC(ITSNSA(IPHAS))
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
C
          RUBRIQ = 'tsource_ns_ce_x_phase'//CPHASE(IPHAS)
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPTSNA  ),IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'tsource_ns_ce_y_phase'//CPHASE(IPHAS)
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPTSNA+1),IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'tsource_ns_ce_z_phase'//CPHASE(IPHAS)
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPTSNA+2),IERROR)
          NBERRO=NBERRO+IERROR
C
          ILU = ILU + 1
C
        ENDIF
      ENDDO
C
C ---> Termes sources turbulence
C
C     Boucle sur les phases communes
      DO IPHAS = 1, MIN(JPHAS,NPHAS)
C     Si le terme est a l'ordre 2
        IF(ISTO2T(IPHAS).GT.0) THEN
C
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
C
C     Keps suite de keps ou v2f
          IF(ITYTUR(IPHAS).EQ.2.AND.
     &         (JTYTUR(IPHAS).EQ.2.OR.JTURB(IPHAS).EQ.50)) THEN
            IPTSTA = IPPROC(ITSTUA(IPHAS))
C
            RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA  ),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+1),IERROR)
            NBERRO=NBERRO+IERROR
C
            ILU = ILU + 1
C
C     Keps suite de Rij, on on utilise 1/2 de tsr11+tsr22+ts33
          ELSEIF(ITYTUR(IPHAS).EQ.2.AND.JTYTUR(IPHAS).EQ.3) THEN
C
C       Ici on veut lire tout ou rien, comme on fait des operations
C          pour reconstruire le terme source sur k (discutable...)
C          IDEBLU = 1 indique qu'on a commence a lire.
C          Si erreur quand on a commence a lire, on remet tout a zero
            IDEBLU = 0
            IANNUL = 0
C
            IPTSTA = IPPROC(ITSTUA(IPHAS))
C
            RUBRIQ = 'tsource_tu_ce_R11_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA  ),IERROR)
            NBERRO=NBERRO+IERROR
            IF(IERROR.EQ.0) IDEBLU = 1
C
            IF(IDEBLU.EQ.1) THEN
              RUBRIQ = 'tsource_tu_ce_R22_phase'//CPHASE(IPHAS)
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+1),IERROR)
              NBERRO=NBERRO+IERROR
              IF(IERROR.NE.0) IANNUL = 1
            ENDIF
C
            IF(IDEBLU.EQ.1.AND.IANNUL.EQ.0) THEN
              DO IEL = 1, NCEL
                PROPCE(IEL,IPTSTA  ) =
     &               PROPCE(IEL,IPTSTA  )+PROPCE(IEL,IPTSTA+1)
              ENDDO
              RUBRIQ = 'tsource_tu_ce_R33_phase'//CPHASE(IPHAS)
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+1),IERROR)
              NBERRO=NBERRO+IERROR
              IF(IERROR.NE.0) IANNUL = 1
            ENDIF
C
            IF(IDEBLU.EQ.1.AND.IANNUL.EQ.0) THEN
              DO IEL = 1, NCEL
                PROPCE(IEL,IPTSTA  ) = 0.5D0*(
     &             PROPCE(IEL,IPTSTA  )+PROPCE(IEL,IPTSTA+1))
              ENDDO
              RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE(IPHAS)
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+1),IERROR)
              NBERRO=NBERRO+IERROR
              IF(IERROR.NE.0) IANNUL = 1
            ENDIF
C
            IF(IANNUL.EQ.1) THEN
              DO IEL = 1, NCEL
                PROPCE(IEL,IPTSTA  ) = 0.D0
                PROPCE(IEL,IPTSTA+1) = 0.D0
              ENDDO
            ENDIF
C
            ILU = ILU + 1
C
C     Rij suite de keps ou de v2f, on utilise 2/3 de tsk
          ELSEIF(ITYTUR(IPHAS).EQ.3.AND.
     &           (JTYTUR(IPHAS).EQ.2.OR.JTURB(IPHAS).EQ.50)) THEN
            IPTSTA = IPPROC(ITSTUA(IPHAS))
C
            RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA  ),IERROR)
            NBERRO=NBERRO+IERROR
C
            D2S3 = 2.D0/3.D0
            DO IEL = 1, NCEL
              TSRII = D2S3*PROPCE(IEL,IPTSTA  )
              PROPCE(IEL,IPTSTA  ) = TSRII
              PROPCE(IEL,IPTSTA+1) = TSRII
              PROPCE(IEL,IPTSTA+2) = TSRII
              PROPCE(IEL,IPTSTA+3) = 0.D0
              PROPCE(IEL,IPTSTA+4) = 0.D0
              PROPCE(IEL,IPTSTA+5) = 0.D0
            ENDDO
C
            RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+6),IERROR)
            NBERRO=NBERRO+IERROR
C
            ILU = ILU + 1
C
C     Rij suite de Rij
          ELSEIF(ITYTUR(IPHAS).EQ.3.AND.JTYTUR(IPHAS).EQ.3) THEN
            IPTSTA = IPPROC(ITSTUA(IPHAS))
C
            RUBRIQ = 'tsource_tu_ce_R11_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA  ),IERROR)
            NBERRO=NBERRO+IERROR
            RUBRIQ = 'tsource_tu_ce_R22_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+1),IERROR)
            NBERRO=NBERRO+IERROR
            RUBRIQ = 'tsource_tu_ce_R33_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+2),IERROR)
            NBERRO=NBERRO+IERROR
            RUBRIQ = 'tsource_tu_ce_R12_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+3),IERROR)
            NBERRO=NBERRO+IERROR
            RUBRIQ = 'tsource_tu_ce_R13_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+4),IERROR)
            NBERRO=NBERRO+IERROR
            RUBRIQ = 'tsource_tu_ce_R23_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+5),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+6),IERROR)
            NBERRO=NBERRO+IERROR
C
            ILU = ILU + 1
C
C     v2f suite de keps : On ne relit que les TS de k et eps
C     on laisse 0 pour phi et f_barre (valeur mise dans iniva0)
          ELSEIF(ITURB(IPHAS).EQ.50.AND.JTYTUR(IPHAS).EQ.2) THEN
            IPTSTA = IPPROC(ITSTUA(IPHAS))
C
            RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA  ),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+1),IERROR)
            NBERRO=NBERRO+IERROR
C
            ILU = ILU + 1
C
C     v2f suite de Rij : On ne relit que les TS de k et eps
C     on laisse 0 pour phi et f_barre (valeur mise dans iniva0)
C     on utilise 1/2 de tsr11+tsr22+ts33 pour le ts de k
          ELSEIF(ITURB(IPHAS).EQ.50.AND.JTYTUR(IPHAS).EQ.3) THEN
C
C       Ici on veut lire tout ou rien, comme on fait des operations
C          pour reconstruire le terme source sur k (discutable...)
C          IDEBLU = 1 indique qu'on a commence a lire.
C          Si erreur quand on a commence a lire, on remet tout a zero
            IDEBLU = 0
            IANNUL = 0
C
            IPTSTA = IPPROC(ITSTUA(IPHAS))
C
            RUBRIQ = 'tsource_tu_ce_R11_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA  ),IERROR)
            NBERRO=NBERRO+IERROR
            IF(IERROR.EQ.0) IDEBLU = 1
C
            IF(IDEBLU.EQ.1) THEN
              RUBRIQ = 'tsource_tu_ce_R22_phase'//CPHASE(IPHAS)
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+1),IERROR)
              NBERRO=NBERRO+IERROR
              IF(IERROR.NE.0) IANNUL = 1
            ENDIF
C
            IF(IDEBLU.EQ.1.AND.IANNUL.EQ.0) THEN
              DO IEL = 1, NCEL
                PROPCE(IEL,IPTSTA  ) =
     &               PROPCE(IEL,IPTSTA  )+PROPCE(IEL,IPTSTA+1)
              ENDDO
              RUBRIQ = 'tsource_tu_ce_R33_phase'//CPHASE(IPHAS)
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+1),IERROR)
              NBERRO=NBERRO+IERROR
              IF(IERROR.NE.0) IANNUL = 1
            ENDIF
C
            IF(IDEBLU.EQ.1.AND.IANNUL.EQ.0) THEN
              DO IEL = 1, NCEL
                PROPCE(IEL,IPTSTA  ) = 0.5D0*(
     &             PROPCE(IEL,IPTSTA  )+PROPCE(IEL,IPTSTA+1))
              ENDDO
              RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE(IPHAS)
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+1),IERROR)
              NBERRO=NBERRO+IERROR
              IF(IERROR.NE.0) IANNUL = 1
            ENDIF
C
            IF(IANNUL.EQ.1) THEN
              DO IEL = 1, NCEL
                PROPCE(IEL,IPTSTA  ) = 0.D0
                PROPCE(IEL,IPTSTA+1) = 0.D0
              ENDDO
            ENDIF
C
            ILU = ILU + 1
C
C     v2f suite de v2f
          ELSEIF(ITURB(IPHAS).EQ.50.AND.JTURB(IPHAS).EQ.50) THEN
            IPTSTA = IPPROC(ITSTUA(IPHAS))
C
            RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA  ),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+1),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'tsource_tu_ce_phi_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+2),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'tsource_tu_ce_fb_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+3),IERROR)
            NBERRO=NBERRO+IERROR
C
            ILU = ILU + 1
C
C     k-omega suite de k-omega
C     Pour le k-omega, on ne relit les termes sources que si on est deja en k-omega.
C     Il n'est en effet pas forcement tres judicieux de relire le TS de k et de mettre
C     0 pour le TS de eps.
C     Quant a essayer de transformer le TS de omega en TS de eps ... no comment !
          ELSEIF(ITURB(IPHAS).EQ.60.AND.JTURB(IPHAS).EQ.60) THEN
            IPTSTA = IPPROC(ITSTUA(IPHAS))
C
            RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA  ),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'tsource_tu_ce_omega_phase'//CPHASE(IPHAS)
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSTA+1),IERROR)
            NBERRO=NBERRO+IERROR
C
            ILU = ILU + 1
C
          ENDIF
C
        ENDIF
      ENDDO
C
C ---> Termes sources scalaires
C
C     Boucle sur les scalaires
      DO ISCAL = 1, NSCAL
        ISCO = ISCOLD(ISCAL)
C     Si il y a un correspondant
        IF(ISCO.GT.0) THEN
C     Si ce correspondant est ok pour le format
          IF(ISCO.LE.NFMTSC) THEN
C     Si le terme est a l'ordre 2
            IF(ISSO2T(ISCAL).GT.0) THEN
              IPTSCA = IPPROC(ITSSCA(ISCAL))
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              WRITE(CAR4,'(I4.4)')ISCO
              RUBRIQ = 'tsource_sc_ce_scalaire'//CAR4
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSCA  ),IERROR)
              NBERRO=NBERRO+IERROR
              ILU = ILU + 1
            ENDIF
          ELSE
            NBERRO=NBERRO-1
            ILU = ILU + 1
          ENDIF
        ENDIF
      ENDDO
C
C     Si erreur, on previent mais pas stop :
C       (on n'a pas forcement les ts a l'ordre 2 avant)
      IF (NBERRO.NE.0) THEN
          CAR54 =
     &         'LECTURE DES TERMES SOURCES                            '
          WRITE(NFECRA,8300)CAR54
      ENDIF
C
      IF(ILU.NE.0) THEN
        CAR54 =' Fin de la lecture des termes sources                 '
        WRITE(NFECRA,1110)CAR54
      ENDIF
C
C=======================================================================
C 8. MOYENNES
C=======================================================================
C
C     Indicateur ok (=0) ou non (>0)
      IOK = 0
C
      ILU = 0
C
      DO IMOM = 1, NBMOMT
C
        IMOLD = IMOOLD(IMOM)
C     Si on doit lire la moyenne
        IF(IMOLD.GT.0) THEN
          ILU = 1
C     Si ce correspondant est ok pour le format
          IF(IMOLD.LE.NFMTMO) THEN
            WRITE(CAR4,'(I4.4)')IMOLD
          ELSE
            CAR4 = CINDFM
          ENDIF
C
C       On la lit
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
          RUBRIQ = 'cumul_ce_moment'//CAR4
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPPROC(ICMOME(IMOM))),IERROR)
          NBERRO=NBERRO+IERROR
C
C       Si ca ne marche pas, on s'arrete
C         (on pourrait tenter de prendre des mesures correctives,
C          mais le cumul de la duree associee sert peut etre a d'autres
C            moyennes
C          en outre, si l'utilisateur a indique qu'il voulait relire
C            une moyenne, c'est qu'il veut faire un calcul propre
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,9300)IMOLD,IMOM,IMOM
            IOK = IOK + 1
          ENDIF
C
C
C
C       Si on a reussi a lire la moyenne, il faut obligatoirement disposer
C         du cumul de duree associe, sinon, on devra s'arreter.
C       IOKDT different de 0 indiquera qu'on n'a pas pu l'obtenir.
          IOKDT = 0
C
C       On cherche le numero (local au fichier suite)
C         du cumul de temps de la moyenne du calcul precedent
C         correspondant a la moyenne IMOM du calcul courant (ouf!)
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'numero_cumul_temps_moment'//CAR4
          NUMERO = 0
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                NUMERO,IERROR)
          NBERRO=NBERRO+IERROR
C
C       Si on n'a pas trouve, on a echoue
          IF(NUMERO.EQ.0.OR.IERROR.NE.0) THEN
            IOKDT = 1
C
C       Sinon, on a trouve un cumul en temps correspondant
          ELSE
C
C         Si NUMERO > 0, il s'agissait d'un cumul variable en espace
C             et le cumul courant est forcement variable en espace
            IF(NUMERO.GT.0) THEN
C
              IF(NUMERO.LE.NFMTMO) THEN
                WRITE(CNUM4,'(I4.4)')NUMERO
              ELSE
                CNUM4 = CINDFM
              ENDIF
C
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              RUBRIQ = 'cumul_temps_ce_'//CNUM4
              IDTCM  = IPPROC(ICDTMO(IDTMOM(IMOM)))
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IDTCM),IERROR)
              NBERRO=NBERRO+IERROR
C
C           Si on n'a pas lu, on a echoue
              IF(IERROR.NE.0) THEN
                IOKDT = 1
              ENDIF
C
C         Si NUMERO < 0, il s'agissait d'un cumul uniforme
C             et le cumul courant est variable en espace ou non
            ELSEIF(NUMERO.LT.0) THEN
C
              NUMERO = -NUMERO
              IF(NUMERO.LE.NFMTMO) THEN
                WRITE(CNUM4,'(I4.4)')NUMERO
              ELSE
                CNUM4 = CINDFM
              ENDIF
C
              ITYSUP = 0
              NBVAL  = 1
              IRTYP  = 2
              RUBRIQ = 'cumul_temps_'//CNUM4
              CDTCM  = 0.D0
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    CDTCM,IERROR)
              NBERRO=NBERRO+IERROR
C
C           Si on n'a pas lu, on a echoue
              IF(IERROR.NE.0) THEN
                IOKDT = 1
C
C           Sinon, selon la nature du cumul actuel, on l'affecte
              ELSE
                IF(IDTMOM(IMOM).GT.0) THEN
                  IDTCM  = IPPROC(ICDTMO(IDTMOM(IMOM)))
                  DO IEL = 1, NCEL
                    PROPCE(IEL,IDTCM) = CDTCM
                  ENDDO
                ELSEIF(IDTMOM(IMOM).LT.0) THEN
                  IDTCM  = -IDTMOM(IMOM)
                  DTCMOM(IDTCM) = CDTCM
                ENDIF
              ENDIF
C
            ENDIF
C
          ENDIF
C
          IF(IOKDT.NE.0) THEN
            WRITE(NFECRA,9310)IMOLD,IMOM,IMOM
            IOK = IOK + 1
          ENDIF
C
        ENDIF
      ENDDO
C
C     Si pb on s'arrete car si on cherche a faire des moyennes
C       c'est qu'on veut vraiement les relire
      IF(IOK.NE.0) THEN
        CALL CSEXIT(1)
      ENDIF
C
      IF(ILU.GT.0) THEN
        CAR54 = ' Fin de la lecture des moyennes temporelles           '
        WRITE(NFECRA,1110)CAR54
      ENDIF
C
C=======================================================================
C 9. DISTANCE A LA PAROI
C=======================================================================
C
      ILU = 0
      NBERRO = 0
C
C     MODE DE CALCUL DIRECT (NON COMPATIBLE PARALLELISME ET PERIODICITE)
C
      IF(ABS(ICDPAR).EQ.2) THEN
C
C     On la lit si on en a besoin uniquement.
C
C     Si l'utilisateur a force le recalcul, on ne la lit pas
C       il faudra la mettre a jour (sauf si zero pas de temps).
C
C     Sinon, on cherche a la lire.
C
C     On ne relit les numeros des faces de bord que si on a toujours
C       le meme nombre de faces de bord (sinon, la numerotation a change)
        IF(INEEDY.EQ.1) THEN
          IF(ICDPAR.EQ.2.OR.INPDT0.EQ.1) THEN
            IF(NFABOK.EQ.1) THEN
C
              DO IPHAS = 1, MIN(NPHAS,JPHAS)
                IF(IIFAPA(IPHAS).GT.0) THEN
                  ITYSUP = 1
                  NBVAL  = 1
                  IRTYP  = 1
                  RUBRIQ = 'num_fac_par_ce_phase'//CPHASE(IPHAS)
                  CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &                    IRTYP,IA(IIFAPA(IPHAS)),IERROR)
                  NBERRO=NBERRO+IERROR
                  ILU   = ILU + 1
                ENDIF
              ENDDO
C
            ENDIF
          ENDIF
        ENDIF
C
C
C     MODE DE CALCUL PAR EQUATION DE DIFFUSION
C
      ELSEIF(ABS(ICDPAR).EQ.1) THEN
C
C     On la lit si on en a besoin uniquement.
C
C     Si l'utilisateur a force le recalcul, on ne la lit pas
C       il faudra la mettre a jour (sauf si zero pas de temps).
C
C     Sinon, on cherche a la lire.
C       On pourrait la relire aussi quand le nombre de faces a
C       change, mais il vaut mieux la recalculer au cas ou des faces de
C       paroi auraient disparu
C       Si on arrive a la lire, on note qu'elle est a jour (sauf si ALE).
C
        IF(INEEDY.EQ.1) THEN
          IF(ICDPAR.EQ.1.OR.INPDT0.EQ.1) THEN
            IF(NFABOK.EQ.1) THEN
              IPHAS  = 1
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              RUBRIQ = 'dist_fac_par_ce_phase'//CPHASE(IPHAS)
              CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    RA(IDIPAR),IERROR)
              NBERRO=NBERRO+IERROR
              IF(IERROR.EQ.0 .AND. IALE.EQ.0 ) THEN
                IMAJDY = 1
              ENDIF
              ILU   = ILU + 1
            ENDIF
          ENDIF
        ENDIF
C
      ENDIF
C
      IF (NBERRO.NE.0) THEN
        CAR54 =
     &         'LECTURE DE LA DISTANCE A LA PAROI                     '
        WRITE(NFECRA,8300)CAR54
      ENDIF
C
      IF(ILU.NE.0) THEN
        CAR54=' Fin de la lecture de la distance a la paroi          '
        WRITE(NFECRA,1110)CAR54
      ENDIF
C
C
C=======================================================================
C 10.  FORCE EXTERIEURE
C=======================================================================
C
      IF(IPHYDR.EQ.1) THEN
        NBERRO=0
        DO IPHAS = 1, MIN(NPHAS,JPHAS)
C
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
C
          RUBRIQ = 'force_ext_ce_x_phase'//CPHASE(IPHAS)
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                FRCXT(1,1,IPHAS),IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'force_ext_ce_y_phase'//CPHASE(IPHAS)
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                FRCXT(1,2,IPHAS),IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'force_ext_ce_z_phase'//CPHASE(IPHAS)
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                FRCXT(1,3,IPHAS),IERROR)
          NBERRO=NBERRO+IERROR
        ENDDO
C
       IF (NBERRO.NE.0) THEN
          CAR54 =
     &         'LECTURE DES FORCES EXTERIEURES                        '
          WRITE(NFECRA,8300)CAR54
        ENDIF
C
        CAR54 =' Fin de la lecture des forces exterieures             '
        WRITE(NFECRA,1110)CAR54
C
      ENDIF
C
C
C=======================================================================
C 11.  DEPLACEMENT AUX NOEUDS EN ALE
C=======================================================================
C
      IF (IALE.EQ.1 .AND. JALE.EQ.1) THEN
        NBERRO = 0
C
        ITYSUP = 4
        NBVAL  = 1
        IRTYP  = 2
C
        RUBRIQ = 'deplact_x_no'
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              RA(IDEPAL),IERROR)
        NBERRO=NBERRO+IERROR
        RUBRIQ = 'deplact_y_no'
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              RA(IDEPAL+NNOD),IERROR)
        NBERRO=NBERRO+IERROR
        RUBRIQ = 'deplact_z_no'
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              RA(IDEPAL+2*NNOD),IERROR)
        NBERRO=NBERRO+IERROR
C
C Si JALE=1, on doit avoir le deplacement dans le fichier suite, sinon
C   les resultats relus n'ont pas de sens -> on s'arrete si pb
        IF (NBERRO.NE.0) THEN
          WRITE(NFECRA,9320)
          CALL CSEXIT(1)
        ENDIF
C
        NBERRO = 0
        RUBRIQ = 'type_visc_mail'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &       JORTVM,IERROR)

        RUBRIQ = 'visc_maillage_x'
        ITYSUP = 1
        NBVAL  = 1
        IRTYP  = 2
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &       PROPCE(1,IPPROC(IVISMA(1))),IERROR)
C
        IF (IORTVM.EQ.1) THEN
          IF (JORTVM.EQ.1) THEN
            RUBRIQ = 'visc_maillage_y'
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           PROPCE(1,IPPROC(IVISMA(2))),IERROR)
            RUBRIQ = 'visc_maillage_z'
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           PROPCE(1,IPPROC(IVISMA(3))),IERROR)
          ELSE
            DO IEL = 1, NCEL
              IPCVMX = IPPROC(IVISMA(1))
              IPCVMY = IPPROC(IVISMA(2))
              IPCVMZ = IPPROC(IVISMA(3))
              PROPCE(IEL,IPCVMY) = PROPCE(IEL,IPCVMX)
              PROPCE(IEL,IPCVMZ) = PROPCE(IEL,IPCVMX)
            ENDDO
          ENDIF
        ENDIF
C
        CAR54 =' Fin de la lecture des donnees ALE                    '
        WRITE(NFECRA,1110)CAR54
C
        NBERRO=0
        RUBRIQ = 'nombre_structures'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &       NBSTRU,IERROR)
        NBERRO=NBERRO+IERROR
C
        IF (NBSTRU.GT.0) THEN
C
          NFMTST = 99
          CINDST='YY'
          DO ISTR = 1, MIN(NBSTRU,NSTRMX)
            WRITE(CSTRUC(ISTR),'(I2.2)') ISTR
          ENDDO
          DO ISTR = MIN(NBSTRU,NFMTST)+1,NBSTRU
            CSTRUC(ISTR) = CINDST
          ENDDO
          IF(NSTRMX.GT.NFMTST) THEN
            WRITE(NFECRA,8004)NFMTST,NSTRMX
          ENDIF
C
          DO ISTR = 1, NBSTRU
C
            RUBRIQ = 'donnees_structure_'//CSTRUC(ISTR)
            ITYSUP = 0
            NBVAL  = 27
            IRTYP  = 2
C
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &           TMPSTR,IERROR)
            NBERRO=NBERRO+IERROR
C
            DO II = 1, 3
              XSTR  (II,ISTR) = TMPSTR(   II)
              XPSTR (II,ISTR) = TMPSTR(3 +II)
              XPPSTR(II,ISTR) = TMPSTR(6 +II)
              XSTA  (II,ISTR) = TMPSTR(9 +II)
              XPSTA (II,ISTR) = TMPSTR(12+II)
              XPPSTA(II,ISTR) = TMPSTR(15+II)
              XSTP  (II,ISTR) = TMPSTR(18+II)
              FORSTR(II,ISTR) = TMPSTR(21+II)
              FORSTA(II,ISTR) = TMPSTR(24+II)
            ENDDO
C
          ENDDO
C
          CAR54 =' Fin de la lecture des donnees des structures ALE   '
          WRITE(NFECRA,1110)CAR54
C
        ENDIF
C
        IF (NBERRO.NE.0) THEN
          WRITE(NFECRA,9321)
          CALL CSEXIT(1)
        ENDIF
C
      ENDIF
C
C=======================================================================
C 12. LECTURE DES INFORMATIONS COMPLEMENTAIRES COMBUSTION GAZ, CP ET
C                                                                  FUEL
C=======================================================================
C
      NBERRO = 0
      ILU = 0
C
C     Modele COD3P :
C     ============
C
      IF ( IPPMOD(ICOD3P).GE.0 ) THEN
C
        RUBRIQ = 'hinfue_cod3p'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,HINFUE,
     &              IERROR)
        NBERRO=NBERRO+IERROR
        ILU = ILU + 1
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9400)
        ENDIF
C
        RUBRIQ = 'hinoxy_cod3p'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        ILU    = ILU + 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,HINOXY,
     &              IERROR)
        NBERRO=NBERRO+IERROR
        ILU = ILU + 1
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9400)
        ENDIF
C
        RUBRIQ = 'tinfue_cod3p'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        ILU    = ILU + 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,TINFUE,
     &              IERROR)
        NBERRO=NBERRO+IERROR
        ILU = ILU + 1
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9400)
        ENDIF
C
        RUBRIQ = 'tinoxy_cod3p'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        ILU    = ILU + 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,TINOXY,
     &              IERROR)
        NBERRO=NBERRO+IERROR
        ILU = ILU + 1
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9400)
        ENDIF
C
C       Il faut le meme nbr de faces de bord, sinon on ne lit pas
        IF(NFABOK.EQ.1) THEN
C
          ILU = ILU + 1
C
          IERRCH = 0
C
C       Numero des zones
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'num_zone_fb_cod3p'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IIZFPP), IERROR)
          NBERRO=NBERRO+IERROR
C
C       Type entree Fuel
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientfu_zone_bord_cod3p'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTFU, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C       Type entree Oxydant
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientox_zone_bord_cod3p'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTOX, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C     Par securite, si on ne parvient pas a lire la
C       IENTCPFU ou IENTOX, on remet a zero le numero des zones IZFPPP
C       car il a peut etre ete lu.
C       Ceci permettra d'eviter de se servir des valeurs par defaut
C
          IF(IERRCH.NE.0) THEN
            DO IFAC = 1, NFABOR
              IA(IIZFPP+IFAC-1) = 0
            ENDDO
          ENDIF
C
        ENDIF
C
      ENDIF
C
C     Modele EBU :
C     ==========
C
      IF ( IPPMOD(ICOEBU).GE.0 ) THEN
C
        RUBRIQ = 'temperature_gaz_frais_ebu'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        ILU    = 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,TGF,
     &              IERROR)
        NBERRO=NBERRO+IERROR
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9500)
        ENDIF
C
        RUBRIQ = 'frmel_ebu'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        ILU    = 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,FRMEL,
     &              IERROR)
        NBERRO=NBERRO+IERROR
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9500)
        ENDIF
C
C       Il faut le meme nbr de faces de bord, sinon on ne lit pas
        IF(NFABOK.EQ.1) THEN
C
          ILU = ILU + 1
C
          IERRCH = 0
C
C       Numero des zones
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'num_zone_fb_ebu'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IIZFPP), IERROR)
          NBERRO=NBERRO+IERROR
C
C       Type entree Gaz brulee
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientgb_zone_bord_ebu'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTGB, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C       Type entree gaz frais
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientgf_zone_bord_ebu'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTGF, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C       FMENT
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'fment_zone_bord_ebu'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                FMENT , IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C       TKENT
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'tkent_zone_bord_ebu'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TKENT, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C     Par securite, si on ne parvient pas a lire la
C       IENTCPFU ou IENTOX, on remet a zero le numero des zones IZFPPP
C       car il a peut etre ete lu.
C       Ceci permettra d'eviter de se servir des valeurs par defaut
C
          IF(IERRCH.NE.0) THEN
            DO IFAC = 1, NFABOR
              IA(IIZFPP+IFAC-1) = 0
            ENDDO
          ENDIF
C
        ENDIF
C
      ENDIF
C
C     Modele LWC :
C     ==========
C
      IF ( IPPMOD(ICOLWC).GE.0 ) THEN
C
        RUBRIQ = 'fmin_lwc'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        ILU    = 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,FMIN,
     &              IERROR)
        NBERRO=NBERRO+IERROR
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9600)
        ENDIF
C
        RUBRIQ = 'fmax_lwc'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        ILU    = 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,FMAX,
     &              IERROR)
        NBERRO=NBERRO+IERROR
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9600)
        ENDIF
C
        RUBRIQ = 'hmin_lwc'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        ILU    = 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,HMIN,
     &              IERROR)
        NBERRO=NBERRO+IERROR
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9600)
        ENDIF
C
        RUBRIQ = 'hmax_lwc'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        ILU    = 1
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,HMAX,
     &              IERROR)
        NBERRO=NBERRO+IERROR
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,9600)
        ENDIF
C
C       Il faut le meme nbr de faces de bord, sinon on ne lit pas
        IF(NFABOK.EQ.1) THEN
C
          ILU = ILU + 1
C
          IERRCH = 0
C
C       Numero des zones
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'num_zone_fb_lwc'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IIZFPP), IERROR)
          NBERRO=NBERRO+IERROR
C
C       Type entree Gaz brulee
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientgb_zone_bord_lwc'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTGB, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C       Type entree gaz frais
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientgf_zone_bord_lwc'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTGF, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C       FMENT
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'fment_zone_bord_lwc'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                FMENT , IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C       TKENT
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'tkent_zone_bord_lwc'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TKENT, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C     Par securite, si on ne parvient pas a lire la
C       IENTCPFU ou IENTOX, on remet a zero le numero des zones IZFPPP
C       car il a peut etre ete lu.
C       Ceci permettra d'eviter de se servir des valeurs par defaut
C
          IF(IERRCH.NE.0) THEN
            DO IFAC = 1, NFABOR
              IA(IIZFPP+IFAC-1) = 0
            ENDDO
          ENDIF
C
        ENDIF
C
      ENDIF
C
C     Charbon PuLVerise : masse vol des charbons
      IF ( IPPMOD(ICP3PL).GE.0.OR. IPPMOD(ICPL3C).GE.0 ) THEN
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 2
        IERRCH = 0
        DO ICHA = 1, NCHARB
          IF(ICHA.LE.NFMTCH) THEN
            WRITE(CAR2,'(I2.2)')ICHA
          ELSE
            CAR2 = CINDFC
          ENDIF
          RUBRIQ = 'masse_volumique_charbon'//CAR2
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                RHOCK(ICHA), IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
          ILU = ILU + 1
        ENDDO
        IF (IERRCH.NE.0) THEN
          WRITE(NFECRA,8611)
          DO ICHA = 1, NCHARB
            WRITE(NFECRA,8612)ICHA,RHOCK(ICHA)
          ENDDO
          WRITE(NFECRA,8613)
        ENDIF
      ENDIF
C
C
C     Charbon PuLVerise : type de zones de bord, ientat, ientcp, timpat
C       et x20 pour le calcul de rho au bord en entree
      IF ( IPPMOD(ICP3PL).GE.0.OR. IPPMOD(ICPL3C).GE.0 ) THEN
C
C       Il faut le meme nbr de faces de bord, sinon on ne lit pas
        IF(NFABOK.EQ.1) THEN
C
          ILU = ILU + 1
C
          IERRCH = 0
C
C       Numero des zones
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'num_zone_fb_charbon_pulverise'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IIZFPP), IERROR)
          NBERRO=NBERRO+IERROR
C
C       Type entree air ou cp (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientat_zone_bord_charbon_pulverise'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTAT, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C         ientcp et x20 ne servent pas pour le CP couple Lagrangien (cplphy)
          IF ( IPPMOD(ICP3PL).GE.0 ) THEN
C
            ITYSUP = 0
            NBVAL  = NOZPPM
            IRTYP  = 1
            RUBRIQ = 'ientcp_zone_bord_charbon_pulverise'
            CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  IENTCP, IERROR)
            IERRCH=IERRCH+IERROR
            NBERRO=NBERRO+IERROR
C
C
            ITYSUP = 0
            NBVAL  = NOZPPM
            IRTYP  = 2
C
            IDECAL = 0
            DO ICHA = 1, NCHARB
              DO ICLAPC = 1, NCLPCH(ICHA)
                ICLA = ICLAPC + IDECAL
                IF(ICHA.LE.NFMTCH.AND.ICLAPC.LE.NFMTCL) THEN
                  WRITE(CAR2,'(I2.2)')ICHA
                  WRITE(CAR4,'(I4.4)')ICLAPC
                ELSE
                  CAR2 = CINDFC
                  CAR4 = CINDFL
                ENDIF
                RUBRIQ = 'x20_zone_bord_charbon'//CAR2//'_classe'//CAR4
                CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &               IRTYP,X20(1,ICLA), IERROR)
                IERRCH=IERRCH+IERROR
                NBERRO=NBERRO+IERROR
C
              ENDDO
            ENDDO
C
          ENDIF
C
C       Temperature
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'timpat_zone_bord_charbon_pulverise'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TIMPAT, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C     Par securite, si on ne parvient pas a lire la temperature TIMPAT,
C       IENTCP ou IENTAT, on remet a zero le numero des zones IZFPPP
C       car il a peut etre ete lu.
C       Ceci permettra d'eviter de se servir des valeurs par defaut (=0)
C       de TIMPAT dans cpphyv et cplphy.
          IF(IERRCH.NE.0) THEN
            DO IFAC = 1, NFABOR
              IA(IIZFPP+IFAC-1) = 0
            ENDDO
          ENDIF
C
        ENDIF
C
      ENDIF
C
C
C     FUEL : type de zones de bord, ientat, ientfl, timpat
C       qimpat et qimpfl pour le calcul de rho au bord en entree
      IF ( IPPMOD(ICFUEL).GE.0 ) THEN
C
C       Il faut le meme nbr de faces de bord, sinon on ne lit pas
        IF(NFABOK.EQ.1) THEN
C
          ILU = ILU + 1
C
          IERRCH = 0
C
C       Numero des zones
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'num_zone_fb_fuel'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IIZFPP), IERROR)
          NBERRO=NBERRO+IERROR
C
C       Type entree air ou fuel (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientat_zone_bord_fuel'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTAT, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientfl_zone_bord_fuel'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTFL, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C       TIMPAT
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'timpat_zone_bord_fuel'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TIMPAT, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C       QIMPAT
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'qimpat_zone_bord_fuel'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                QIMPAT, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C       QIMPFL
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'qimpfl_zone_bord_fuel'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                QIMPFL, IERROR)
          IERRCH=IERRCH+IERROR
          NBERRO=NBERRO+IERROR
C
C     Par securite, si on ne parvient pas a lire la temperature TIMPAT,
C       IENTCP ou IENTAT, on remet a zero le numero des zones IZFPPP
C       car il a peut etre ete lu.
C       Ceci permettra d'eviter de se servir des valeurs par defaut (=0)
C       de TIMPAT dans cpphyv et cplphy.
          IF(IERRCH.NE.0) THEN
            DO IFAC = 1, NFABOR
              IA(IIZFPP+IFAC-1) = 0
            ENDDO
          ENDIF
C
        ENDIF
C
      ENDIF
C
      IF (NBERRO.NE.0) THEN
        CAR54 =
     &       'LECTURE DES INFORMATIONS COMBUSTION                   '
        WRITE(NFECRA,8300)CAR54
      ENDIF
C
      IF(ILU.NE.0) THEN
        CAR54=' Fin de la lecture des informations combustion        '
        WRITE(NFECRA,1110)CAR54
      ENDIF
C
C=======================================================================
C 13. LECTURE DES INFORMATIONS COMPLEMENTAIRES ELECTRIQUES
C=======================================================================
C
      NBERRO=0
      ILU  = 0
C
C     Recalage des CL pot des versions electriques
C
      IF ( IPPMOD(IELJOU).GE.1       ) THEN
        IF(IELCOR.EQ.1) THEN
          ILU = ILU + 1
          RUBRIQ = 'coeff_recalage_joule'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                COEJOU,IERROR)
          NBERRO=NBERRO+IERROR
        ENDIF
      ENDIF
      IF ( IPPMOD(IELARC).GE.1  .OR. IPPMOD(IELJOU).GE.1 ) THEN
        IF(IELCOR.EQ.1) THEN
          ILU = 1
          RUBRIQ = 'ddpot_recalage_arc_elec'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                DPOT  ,IERROR)
          NBERRO=NBERRO+IERROR
        ENDIF
      ENDIF
C
C ---> Termes sources des versions electriques
C
      IF ( IPPMOD(IELJOU).GE.1 .OR.
     &     IPPMOD(IELARC).GE.1 .OR.
     &     IPPMOD(IELION).GE.1       ) THEN
C
        IPCEFJ = IPPROC(IEFJOU)
        ITYSUP = 1
        NBVAL  = 1
        IRTYP  = 2
C
        RUBRIQ = 'tsource_sc_ce_joule'
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              PROPCE(1,IPCEFJ  ),IERROR)
        NBERRO=NBERRO+IERROR
        ILU = ILU + 1
C
      ENDIF
C
      IF( IPPMOD(IELARC).GE.1 ) THEN
C
        IPCLA1 = IPPROC(ILAPLA(1))
        IPCLA2 = IPPROC(ILAPLA(2))
        IPCLA3 = IPPROC(ILAPLA(3))
        ITYSUP = 1
        NBVAL  = 1
        IRTYP  = 2
C
        RUBRIQ = 'tsource_ns_ce_x_laplace'
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              PROPCE(1,IPCLA1  ),IERROR)
        NBERRO=NBERRO+IERROR
        ILU = ILU + 1
C
        RUBRIQ = 'tsource_ns_ce_y_laplace'
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              PROPCE(1,IPCLA2  ),IERROR)
        NBERRO=NBERRO+IERROR
        ILU = ILU + 1
C
        RUBRIQ = 'tsource_ns_ce_z_laplace'
        CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              PROPCE(1,IPCLA3  ),IERROR)
        NBERRO=NBERRO+IERROR
        ILU = ILU + 1
C
      ENDIF
C
      IF (NBERRO.NE.0) THEN
        CAR54 =
     &       'LECTURE DES INFORMATIONS ELECTRIQUES                  '
        WRITE(NFECRA,8300)CAR54
      ENDIF
C
      IF(ILU.NE.0) THEN
        CAR54=' Fin de la lecture des informations electriques       '
        WRITE(NFECRA,1110)CAR54
      ENDIF
C
C=======================================================================
C 14.  FERMETURE DU FICHIER SUITE AUXILAIRE
C=======================================================================
C
C
C     Fermeture du fichier suite auxilaire
      CALL CLSSUI(IMPAMX,IERROR)
C
      IF (IERROR.NE.0) THEN
         WRITE(NFECRA,8900) FICAMX
      ENDIF
C
      WRITE(NFECRA,1200)
C
C=======================================================================
C 15. SORTIE
C=======================================================================
C
      RETURN
C
C=======================================================================
C 16. FORMATS
C=======================================================================
C
C --- ETAPES
C
 1000 FORMAT(/,
     &     3X,'   LECTURE DU FICHIER SUITE AUXILIAIRE               ',/)
 1100 FORMAT(' Debut de la lecture                                    ')
 1110 FORMAT('  ',A54                                                  )
 1200 FORMAT(' Fin de la lecture                                      ')
C
C --- INFORMATIONS
C
 7000 FORMAT(/,
     &'   PHASE ',I4,' :                                            ',/,
     &'   Mise a jour du point de reference pour la pression totale ',/,
     &'     par relecture du fichier suite                          ',/,
     &'    XYZP0 = ',       E14.5,        E14.5,        E14.5        ,/)
C
C --- MISES EN GARDE
C
 8000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@      Le nombre de phases    maximal NPHSMX supporte par le ',/,
     &'@        format d''ecriture du fichier suite est             ',/,
     &'@        NFMTPH = ',I10                                       ,/,
     &'@      On a ici un nombre de phases    maximal superieur     ',/,
     &'@        NPHSMX = ',I10                                       ,/,
     &'@      Si le nombre de phases effectif est superieur, elles  ',/,
     &'@        ne seront pas relues.                               ',/,
     &'@                                                            ',/,
     &'@    Le calcul sera execute.                                 ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme lecamx.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@      Le nombre de scalaires maximal NSCAMX supporte par le ',/,
     &'@        format d''ecriture du fichier suite est             ',/,
     &'@        NFMTSC = ',I10                                       ,/,
     &'@      On a ici un nombre de scalaires maximal superieur     ',/,
     &'@        NSCAMX = ',I10                                       ,/,
     &'@      On ne pourra pas relire les scalaires dont le numero  ',/,
     &'@        est superieur                                       ',/,
     &'@                                                            ',/,
     &'@    Le calcul sera execute.                                 ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme lecamx.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@      Le nombre de flux de masse max NVARMX supporte par le ',/,
     &'@        format d''ecriture du fichier suite est             ',/,
     &'@        NFMTFL = ',I10                                       ,/,
     &'@      On a ici un nombre de flux      maximal superieur     ',/,
     &'@        NVARMX = ',I10                                       ,/,
     &'@      On ne pourra pas relire les flux      dont le numero  ',/,
     &'@        est superieur                                       ',/,
     &'@                                                            ',/,
     &'@    Le calcul sera execute.                                 ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme lecamx.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8003 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@      Le nombre de moments mx        NBMOMX supporte par le ',/,
     &'@        format d''ecriture du fichier suite est             ',/,
     &'@        NFMTMO = ',I10                                       ,/,
     &'@      On a ici un nombre de moments   maximal superieur     ',/,
     &'@        NBMOMX = ',I10                                       ,/,
     &'@      On ne pourra pas relire les moments   dont le numero  ',/,
     &'@        est superieur                                       ',/,
     &'@                                                            ',/,
     &'@    Le calcul sera execute.                                 ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme lecamx.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8004 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@      Le nombre de structures maximal NSTRMX supporte par le',/,
     &'@        format d''ecriture du fichier suite est             ',/,
     &'@        NFMTST = ',I10                                       ,/,
     &'@      On a ici un nombre de structures maximal superieur    ',/,
     &'@        NSTRMX = ',I10                                       ,/,
     &'@      Si le nombre de structures effectif est superieur,    ',/,
     &'@        elles ne seront pas relues.                         ',/,
     &'@                                                            ',/,
     &'@    Le calcul sera execute.                                 ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme lecamx.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,
     &'@    *********                                               ',/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Le nombre de faces ',A8  ,' a ete modifie.              ',/,
     &'@                                                            ',/,
     &'@    Le calcul peut etre execute mais les donnees            ',/,
     &'@      sur les faces ',A8  ,' ne seront pas relues           ',/,
     &'@      dans le fichier suite.                                ',/,
     &'@    Elles seront initialisees par des valeurs par defaut.   ',/,
     &'@                                                            ',/,
     &'@    Cette situation peut se produire lorsque le fichier     ',/,
     &'@      suite est issu d''un calcul realise avec des options  ',/,
     &'@      de recollement differentes ou lorsque l''on modifie   ',/,
     &'@      la prise en compte de periodicite.                    ',/,
     &'@    Cette situation peut egalement se produire lorsque l''on',/,
     &'@      realise une suite sur une machine de calcul differente',/,
     &'@      et que le jeu de la precision machine modifie le      ',/,
     &'@      nombre de faces issues des recollements.              ',/,
     &'@                                                            ',/,
     &'@    Cette situation peut enfin se produire lorsque le       ',/,
     &'@      fichier suite auxiliaire ne correspond pas au cas     ',/,
     &'@      traite.                                               ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite auxiliaire utilise        ',/,
     &'@      correspond bien au cas traite                         ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,
     &'@    *********                                               ',/,
     &'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,
     &'@                                                            ',/,
     &'@    Le nombre de phases a ete modifie.                      ',/,
     &'@    Le calcul peut etre execute.                            ',/,
     &'@                                                            ',/,
     &'@    Il est cependant conseille de verifier                  ',/,
     &'@      les dimensions suivantes dans usini1 :                ',/,
     &'@                                                            ',/,
     &'@               NPHAS                                        ',/,
     &'@  AMONT : ',I10                                              ,/,
     &'@  ACTUEL: ',I10                                              ,/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite auxiliaire utilise        ',/,
     &'@      correspond bien au cas traite                         ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8220 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,
     &'@    *********                                               ',/,
     &'@      REPRISE  DE CALCUL           AVEC ITURB(',I3,') = ',I4 ,/,
     &'@      A PARTIR D''UN CALCUL REALISE AVEC ITURB(',I3,') = ',I4,/,
     &'@                                                            ',/,
     &'@    Le modele de turbulence a ete modifie.                  ',/,
     &'@                                                            ',/,
     &'@    Il est conseille cependant de                           ',/,
     &'@      verifier la valeur de ITURB(',I2,') dans usini1       ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite auxiliaire utilise        ',/,
     &'@      correspond bien au cas traite                         ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8300 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,
     &'@    *********                                               ',/,
     &'@      ', A54                                                 ,/,
     &'@                                                            ',/,
     &'@    Certaines grandeurs n''ont pas pu etre lues dans le     ',/,
     &'@      fichier suite auxiliaire.                             ',/,
     &'@    Elles seront initialisees par des valeurs par defaut.   ',/,
     &'@                                                            ',/,
     &'@    Cette situation peut se produire lorsque le fichier     ',/,
     &'@      suite est issu d''un calcul realise avec des options  ',/,
     &'@      differentes ou lorsqu''il a ete endommage.            ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,
     &'@    *********                                               ',/,
     &'@      REPRISE  DE CALCUL           AVEC IDTVAR = ',I10       ,/,
     &'@      A PARTIR D''UN CALCUL REALISE AVEC IDTVAR = ',I10      ,/,
     &'@                                                            ',/,
     &'@    Le mode de marche en temps a ete modifie.               ',/,
     &'@    La valeur (uniforme) du pas de temps est                ',/,
     &'@      DTREF = ',E12.4   ,' fournie dans usini1.             ',/,
     &'@                                                            ',/,
     &'@    Il est conseille cependant de                           ',/,
     &'@      verifier la valeur de IDTVAR dans usini1.             ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite auxiliaire utilise        ',/,
     &'@      correspond bien au cas traite                         ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8611 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      Modele de combustion charbon pulverise.               ',/,
     &'@      On ne trouve pas la masse volumique des charbons dans ',/,
     &'@        le fichier suite. C''est naturel si le calcul       ',/,
     &'@        precedent n''etait pas un calcul charbon pulverise. ',/,
     &'@        La valeur par defaut est utilisee comme valeur      ',/,
     &'@        initiale :                                          ',/,
     &'@         Charbon        rho                                 '  )
 8612 FORMAT(
     &'@        ',I10   ,'  ',E14.5                                    )
 8613 FORMAT(
     &'@                                                            ',/,
     &'@    Le calcul peut etre execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8900 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@    Probleme sur le fichier de nom (',A13,')                ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C --- ERREURS
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@      ERREUR A L''OUVERTURE DU FICHIER SUITE                ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier l''existence et le nom (',A13,') du            ',/,
     &'@        fichier suite dans le repertoire de travail.        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@      TYPE DE FICHIER INCORRECT                             ',/,
     &'@                                                            ',/,
     &'@    Le fichier ',A13      ,' ne semble pas etre un fichier  ',/,
     &'@      suite auxiliaire.                                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise correspond bien   ',/,
     &'@        a un fichier suite auxiliaire.                      ',/,
     &'@    Si necessaire, il est possible de desactiver la lecture ',/,
     &'@      du fichier suite auxiliaire par ILEAUX = 0 (usini1).  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9101 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    Le nombre de cellules a ete modifie                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise correspond bien   ',/,
     &'@        au cas traite.                                      ',/,
     &'@    Si necessaire, il est possible de desactiver la lecture ',/,
     &'@      du fichier suite auxiliaire par ILEAUX = 0 (usini1).  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION: ARRET A LA LECTURE DU FICHIER SUITE          ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@      ', A54                                                 ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise n''a pas ete      ',/,
     &'@      endommage.                                            ',/,
     &'@    Si necessaire, il est possible de desactiver la lecture ',/,
     &'@      du fichier suite auxiliaire par ILEAUX = 0 (usini1).  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A LA LECTURE DU FICHIER SUITE        ',/,
     &'@    *********                                     AUXILIAIRE',/,
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
 9300 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@    La relecture de l''ancienne moyenne ',I10  ,' qui doit  ',/,
     &'@      permettre d''initialiser la nouvelle moyenne ',I10     ,/,
     &'@      a echoue                                              ',/
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas  execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite auxiliaire ou                 ',/,
     &'@      specifier dans usini1 que la moyenne doit etre        ',/,
     &'@      reinitialisee, en indiquant : IMOOLD(',I10   ,') = -1 ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9310 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@    La relecture du cumul de la duree associee a            ',/,
     &'@                    l''ancienne moyenne ',I10  ,' qui doit  ',/,
     &'@      permettre d''initialiser la nouvelle moyenne ',I10     ,/,
     &'@      a echoue                                              ',/
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas  execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite auxiliaire ou                 ',/,
     &'@      specifier dans usini1 que la moyenne doit etre        ',/,
     &'@      reinitialisee, en indiquant : IMOOLD(',I10   ,') = -1 ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9320 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@      ERREUR LORS DE LA LECTURE DU DEPLACEMENT AUX NOEUDS   ',/,
     &'@        DU MAILLAGE (METHODE ALE)                           ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise n''a pas ete      ',/,
     &'@        endommage.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9321 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@      ERREUR LORS DE LA LECTURE DES DONNEES DES STRUCTURES  ',/,
     &'@        MOBILES (METHODE ALE)                               ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise n''a pas ete      ',/,
     &'@        endommage.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@    La relecture des variable CO3DP                          ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas  execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite auxiliaire                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9500 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@    La relecture des variable EBU                           ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas  execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite auxiliaire                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9600 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@    La relecture des variable LWC                           ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas  execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite auxiliaire                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)

C
      END
c@z
