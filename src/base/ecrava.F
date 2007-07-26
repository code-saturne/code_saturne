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
                        SUBROUTINE ECRAVA
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NNOD   ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , FRCXT  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
c@foncb
CFONC
CFONC FONCTION :
CFONC ----------
CFONC  ECRITURE D'UN FICHIER SUITE
CFONC
CFONC  PAS DE STOP (RETURN SI ON ECHOUE)
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
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE NOEUDS                     !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT        )!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! FRCXT(NCELET,! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   3,NPHAS)   !    !     !  HYDROSTATIQUE                       !
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
      INCLUDE "numvar.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "optcal.h"
      INCLUDE "albase.h"
      INCLUDE "alstru.h"
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
C
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
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
      INTEGER          NBMOM2
      PARAMETER       (NBMOM2=NBMOMX*2)
C
      CHARACTER        RUBRIQ*64,CAR2*2,CAR4*4,CAR54*54
      CHARACTER        CINDFP*2,CINDFS*4,CINDFF*4,CINDFM*4
      CHARACTER        CINDFC*2,CINDFL*4
      CHARACTER        CPHASE(NPHSMX)*2 , CSCAL(NSCAMX)*4
      CHARACTER        CFLU  (NVARMX)*4 , CMOM (NBMOMX)*4
      CHARACTER        NOMFLU(NVARMX)*18, NOMRTP(NVARMX)*20
      CHARACTER        NOMCLI(NVARMX)*18
      CHARACTER        CSTRUC(NSTRMX)*2, CINDST*2
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IVAR  , IPHAS , ISCAL , IMOM
      INTEGER          IDECAL, ICLAPC, ICHA  , ICLA
      INTEGER          II    , IVERS , IDTM  , IDTCM
      INTEGER          ICLVAR, ICLVAF, IPTSNA, IPTSTA, IPTSCA
      INTEGER          IERROR, NBERRO, IRTYP , ITYSUP, NBVAL
      INTEGER          NBCTM , IPCEFJ, IPCLA1, IPCLA2, IPCLA3
      INTEGER          NFMTPH, NFMTSC, NFMTFL, NFMTMO, NFMTCH, NFMTCL
      INTEGER          NFMTST
      INTEGER          NBFLU , ILECEC, IECR
      INTEGER          IFAIT (NVARMX)
      INTEGER          ICDTVU(NBMOM2)
      INTEGER          IW1, IFINRA, IFAC, IEL, ISTR
      DOUBLE PRECISION TMPSTR(27)
C
C***********************************************************************
C     A noter :
C        Lorsque qu'il est necessaire d'utiliser un ordre implicite
C        de rangement des variables, on a choisi :
C          P,
C          (U, V, W, turbulence, scalaires)_phase1,
C          (...)_phase2,
C          (...)_...
C          scalaires
C
C          avec turbulence = k, epsilon
C                       ou   R11, R22, R33, R12, R13, R23, epsilon
C                       ou   k, epsilon, phi, f_barre
C                       ou   k, omega
C
C        Ceci est par exemple utilise pour relier les flux de masse aux
C        variables
C
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. VERIFICATIONS DE BASE ET CODAGE DES CHAINES DE CARACTERES
C=======================================================================
C
C  --->  On code en chaine le numero des phases et scalaires
C        ----------------------------------------------------
C
C     Nombre de phases, de scalaires, de flux, de moments et de charbons
C       max pour les formats choisis
      NFMTPH = 99
      NFMTSC = 9999
      NFMTFL = 9999
      NFMTMO = 9999
      NFMTCH = 99
      NFMTCL = 9999
C
C     Indefini (on met qqch de different de lecamo (pour generer une
C       erreur a la lecture)
      CINDFP = 'XX'
      CINDFS = 'XXXX'
      CINDFF = 'XXXX'
      CINDFM = 'XXXX'
      CINDFC = 'XX'
      CINDFL = 'XXXX'
C
C
C     Codage en chaine de caracteres du numero de la phase
      DO IPHAS = 1, MIN(NPHAS ,NFMTPH)
        WRITE(CPHASE(IPHAS),'(I2.2)') IPHAS
      ENDDO
      DO IPHAS = MIN(NPHAS ,NFMTPH)+1,NPHAS
       CPHASE(IPHAS) = CINDFP
      ENDDO
C
C     Codage en chaine de caracteres du numero du scalaire
      DO ISCAL = 1, MIN(NSCAL ,NFMTSC)
        WRITE(CSCAL(ISCAL),'(I4.4)')ISCAL
      ENDDO
      DO ISCAL = MIN(NSCAL ,NFMTSC)+1,NSCAL
        CSCAL(ISCAL) = CINDFS
      ENDDO
C
C     Codage en chaine de caracteres du numero du flux de masse
      DO IVAR = 1, MIN(NVAR  ,NFMTFL)
        WRITE(CFLU(IVAR),'(I4.4)')IVAR
      ENDDO
      DO IVAR = MIN(NVAR  ,NFMTFL)+1,NVAR
        CFLU(IVAR) = CINDFF
      ENDDO
C
C     Codage en chaine de caracteres du numero du moment
      DO IMOM = 1, MIN(NBMOMT,NFMTMO)
        WRITE(CMOM(IMOM),'(I4.4)')IMOM
      ENDDO
      DO IMOM = MIN(NBMOMT,NFMTMO)+1,NBMOMT
        CMOM(IMOM) = CINDFM
      ENDDO
C
C     Verifications pour les formats et les numero de phase
C       et de scalaire en chaine.
C     Avertissement (certaines infos passent a la trappe)
      IF(NPHSMX.GT.NFMTPH) THEN
        WRITE(NFECRA,7000)NFMTPH,NPHSMX
      ENDIF
      IF(NSCAMX.GT.NFMTSC) THEN
        WRITE(NFECRA,7001)NFMTSC,NSCAMX
      ENDIF
      IF(NVARMX.GT.NFMTFL) THEN
        WRITE(NFECRA,7002)NFMTFL,NVARMX
      ENDIF
      IF(NBMOMX.GT.NFMTMO) THEN
        WRITE(NFECRA,7003)NFMTMO,NBMOMX
      ENDIF
      IF(NCHARM.GT.NFMTCH) THEN
        WRITE(NFECRA,7004)NFMTCH,NCHARM
      ENDIF
      IF(NCPCMX.GT.NFMTCL) THEN
        WRITE(NFECRA,7005)NFMTCL,NCPCMX
      ENDIF
C
C
C=======================================================================
C 2. OUVERTURE FICHIER SUITE DE BASE
C=======================================================================
C IFORMA       =0 --> non formate       =1 --> formate
C ILECEC = 2 : ecriture
C
      WRITE(NFECRA,1000)
C
      IF(IFOAVA.EQ.0) THEN
        WRITE(NFECRA,1010)
      ELSEIF(IFOAVA.EQ.1) THEN
        WRITE(NFECRA,1020)
      ENDIF
C
      ILECEC = 2
C
      CALL OPNSUI(FICAVA, LEN(FICAVA), ILECEC, IFOAVA, IMPAVA, IERROR)
C     ===========
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,8000) FICAVA
        RETURN
      ENDIF
C
C
C=======================================================================
C 3. ECRITURE FICHIER SUITE DE BASE
C=======================================================================
C
      WRITE(NFECRA,1100)
C
C
C 3.0 VERSION  : Rubrique "fichier suite ppal"
C=============   Pourrait porter le numero de version si besoin.
C                On ne se sert pas de IVERS (=1.2.0) pour le moment.
C
      IVERS  = 120
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      RUBRIQ = 'version_fichier_suite_principal'
      CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IVERS,
     &            IERROR)
      IF (IERROR.NE.0) THEN
        CAR54='ERREUR A L''ECRITURE DE L''ENTETE                     '
        WRITE(NFECRA,8100)CAR54
      ENDIF
C
C
C 3.1 DIMENSIONS : les dimensions geometriques sont ecrites
C===============   automatiquement lors de l'ouverture du fichier
C                  on ecrit ici les nombres de variables
C
      NBERRO = 0
C
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
C
      RUBRIQ = 'nombre_variables'
      CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NVAR,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'nombre_scalaires'
      CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NSCAL,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'nombre_scalaires_us'
      CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NSCAUS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'nombre_scalaires_pp'
      CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NSCAPP,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'nombre_phases'
      CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NPHAS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      IF (NBERRO.NE.0) THEN
        CAR54='ERREUR A L''ECRITURE DES DIMENSIONS                   '
        WRITE(NFECRA,8100)CAR54
      ENDIF
C
      CAR54 =' Fin de l''ecriture des dimensions                    '
      WRITE(NFECRA,1110)CAR54
C
C 3.2 OPTIONS (Celles servant a donner le nombre de tableaux a lire)
C====================================================================
C Remarque : ces variables ne sont pas toutes utilies pour la relecture
C            en revanche, elles peuvent completer les infos au sein
C            du fichier suite
C
      NBERRO = 0
C
C  ---> Nombre de pas de temps, instant precedent
      RUBRIQ = 'nbre_pas_de_temps'
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NTCABS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
      RUBRIQ = 'instant_precedent'
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 2
      CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,TTCABS,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
C  ---> Modeles de turbulence
      DO IPHAS = 1, NPHAS
        RUBRIQ = 'modele_turbulence_phase'//CPHASE(IPHAS)
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
        CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              ITURB(IPHAS),IERROR)
        NBERRO=NBERRO+IERROR
      ENDDO
C
C
C  ---> Methode ALE
      RUBRIQ = 'methode_ALE'
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IALE,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
C
      IF (NBERRO.NE.0) THEN
        CAR54='ERREUR A L''ECRITURE DES OPTIONS                      '
        WRITE(NFECRA,8100)CAR54
      ENDIF
C
      CAR54 =' Fin de l''ecriture des options                       '
      WRITE(NFECRA,1110)CAR54
C
C 3.3 VARIABLES "PRINCIPALES"
C============================
C
      NBERRO = 0
C
      NOMRTP(IPR(1))='pression_ce_phase'//CPHASE(1)
      DO IPHAS = 1, NPHAS
        NOMRTP(IU(IPHAS))='vitesse_u_ce_phase'//CPHASE(IPHAS)
        NOMRTP(IV(IPHAS))='vitesse_v_ce_phase'//CPHASE(IPHAS)
        NOMRTP(IW(IPHAS))='vitesse_w_ce_phase'//CPHASE(IPHAS)
        IF (ITYTUR(IPHAS).EQ.2) THEN
          NOMRTP(IK(IPHAS))='k_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IEP(IPHAS))='eps_ce_phase'//CPHASE(IPHAS)
        ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
          NOMRTP(IR11(IPHAS))='R11_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IR22(IPHAS))='R22_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IR33(IPHAS))='R33_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IR12(IPHAS))='R12_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IR13(IPHAS))='R13_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IR23(IPHAS))='R23_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IEP(IPHAS))='eps_ce_phase'//CPHASE(IPHAS)
        ELSEIF (ITURB(IPHAS).EQ.50) THEN
          NOMRTP(IK(IPHAS))='k_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IEP(IPHAS))='eps_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IPHI(IPHAS))='phi_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IFB(IPHAS))='fb_ce_phase'//CPHASE(IPHAS)
        ELSEIF (ITURB(IPHAS).EQ.60) THEN
          NOMRTP(IK(IPHAS))='k_ce_phase'//CPHASE(IPHAS)
          NOMRTP(IOMG(IPHAS))='omega_ce_phase'//CPHASE(IPHAS)
        ENDIF
      ENDDO
      IF(NSCAL.GT.0) THEN
        DO ISCAL = 1, NSCAL
          NOMRTP(ISCA(ISCAL))='scalaire_ce_'//CSCAL(ISCAL)
        ENDDO
      ENDIF
      IF (IALE.EQ.1) THEN
        NOMRTP(IUMA)='vit_maillage_u_ce'
        NOMRTP(IVMA)='vit_maillage_v_ce'
        NOMRTP(IWMA)='vit_maillage_w_ce'
      ENDIF
C
C     Dans le cas ou il y a plusieurs phases,
C       on ne veut ecrire la pression qu'une seule fois
C       mais ca tombe bien, car on ne la verra qu'une seule fois
C       dans la liste des variables
C
      DO IVAR = 1, NVAR
C
        ITYSUP = 1
        NBVAL  = 1
        IRTYP  = 2
        RUBRIQ = NOMRTP(IVAR)
        CALL ECRSUI(IMPAVA,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              RTP(1,IVAR),IERROR)
        NBERRO=NBERRO+IERROR
C
      ENDDO
C
      IF (NBERRO.NE.0) THEN
        CAR54='ERREUR A L''ECRITURE DES VARIABLES PRINCIPALES        '
        WRITE(NFECRA,8100)CAR54
      ENDIF
C
      CAR54 =' Fin de l''ecriture des variables principales         '
      WRITE(NFECRA,1110)CAR54
C
C
C 3.4 INFORMATIONS COMPLEMENTAIRES LEGERES (IE UN ENTIER, UN REEL...)
C==========================================
C
      NBERRO = 0
C
      IECR = 0
C

C
      IF (NBERRO.NE.0) THEN
        CAR54='ERREUR A L''ECRITURE DES INFORMATIONS COMPLEMENTAIRES '
        WRITE(NFECRA,8100)CAR54
      ENDIF
C
      IF(IECR.NE.0) THEN
        CAR54 =' Fin de l''ecriture des informations complementaires  '
        WRITE(NFECRA,1110)CAR54
      ENDIF
C
C
C=======================================================================
C 4. FERMETURE FICHIER SUITE DE BASE
C=======================================================================
C
C     Fermeture du fichier suite principal
      CALL CLSSUI(IMPAVA,IERROR)
C
      IF (IERROR.NE.0) THEN
        WRITE(NFECRA,8010) FICAVA
      ENDIF
C
      WRITE(NFECRA,1200)
C
C=======================================================================
C 5. ECRITURE FICHIER SUITE AUXILIAIRE
C=======================================================================
C
C     Si  l'ecriture du fichier suite auxiliaire est demandee
      IF (IECAUX.EQ.1) THEN
C
C 5.0. OUVERTURE FICHIER SUITE AUXILIAIRE
C========================================
C IFORMA       =0 --> non formate       =1 --> formate
C
        WRITE(NFECRA,2000)
C
        IF(IFOAVX.EQ.0) THEN
          WRITE(NFECRA,1010)
        ELSEIF(IFOAVX.EQ.1) THEN
          WRITE(NFECRA,1020)
        ENDIF
C
        ILECEC = 2
        CALL OPNSUI(FICAVX, LEN(FICAVX), ILECEC, IFOAVX, IMPAVX, IERROR)
C       ===========
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,8001) FICAVX
          RETURN
        ENDIF
C
        WRITE(NFECRA,1100)
C
C 5.1 VERSION  : Rubrique "fichier suite auxiliaire"
C=============   Pourrait porter le numero de version si besoin.
C                On ne se sert pas de IVERS (=1.2.0) pour le moment.
C
        IVERS  = 120
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
        RUBRIQ = 'version_fichier_suite_auxiliaire'
        CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &              IVERS,IERROR)
        IF (IERROR.NE.0) THEN
          CAR54='ERREUR A L''ECRITURE DE L''ENTETE                     '
          WRITE(NFECRA,8101)CAR54
        ENDIF
C
C 5.2 DIMENSIONS : les dimensions geometriques sont ecrites
C===============   automatiquement lors de l'ouverture du fichier
C
C
        NBERRO=0
C
C  ---> Nombre de phases
C       On les reecrit ici car on en aura besoin a la relecture
        RUBRIQ = 'nombre_phases'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
        CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NPHAS,
     &              IERROR)
        NBERRO=NBERRO+IERROR
C
C  ---> Nombre de pas de temps, instant precedent
C       On les reecrit ici car on en aura besoin a la relecture
        RUBRIQ = 'nbre_pas_de_temps'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
        CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NTCABS,
     &              IERROR)
        NBERRO=NBERRO+IERROR
C
C  ---> Indicateur de pas de temps variable
        RUBRIQ = 'indic_dt_variable'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
        CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IDTVAR,
     &        IERROR)
        NBERRO=NBERRO+IERROR
C
C  ---> Modeles de turbulence
C       On les reecrit ici car on en aura besoin a la relecture
        DO IPHAS = 1, NPHAS
          RUBRIQ = 'modele_turbulence_phase'//CPHASE(IPHAS)
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 1
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                ITURB(IPHAS),IERROR)
          NBERRO=NBERRO+IERROR
        ENDDO
C
      RUBRIQ = 'methode_ALE'
      ITYSUP = 0
      NBVAL  = 1
      IRTYP  = 1
      CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,IALE,
     &            IERROR)
      NBERRO=NBERRO+IERROR
C
        IF (NBERRO.NE.0) THEN
          CAR54='ERREUR A L''ECRITURE DES DIMENSIONS ET DES OPTIONS    '
          WRITE(NFECRA,8101)CAR54
        ENDIF
C
        CAR54 =' Fin de l''ecriture des dimensions et des options     '
        WRITE(NFECRA,1110)CAR54
C
C 5.3 ECRITURE DES VARIABLES
C===========================
C
C --->  Proprietes physiques
C
        NBERRO=0
C
        DO IPHAS = 1, NPHAS
C
C     Point de reference pour la pression totale
C     On n'ecrit que si XYZP0 a ete specifie par l'utilisateur ou
C       calcule a partir de faces de sorties ou de Dirichlet
          IF (IXYZP0(IPHAS).EQ.1) THEN
            RUBRIQ = 'ref_presstot'//CPHASE(IPHAS)
            ITYSUP = 0
            NBVAL  = 3
            IRTYP  = 2
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  XYZP0(1,IPHAS),IERROR)
            NBERRO=NBERRO+IERROR
          ENDIF
C
C     Masse volumique si elle est variable uniquement
          IF(IROVAR(IPHAS).EQ.1) THEN
C          Masse volumique - cellules
            RUBRIQ = 'rho_ce_phase'//CPHASE(IPHAS)
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPPROC(IROM(IPHAS))),IERROR)
            NBERRO=NBERRO+IERROR
C
C          Masse volumique - faces de bord
            RUBRIQ = 'rho_fb_phase'//CPHASE(IPHAS)
            ITYSUP = 3
            NBVAL  = 1
            IRTYP  = 2
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPFB(1,IPPROB(IROM(IPHAS))),IERROR)
            NBERRO=NBERRO+IERROR
          ENDIF
C
C     On n'ecrit les proprietes physiques que si on les extrapole.
C       On pourrait les ecrire a tous les coups en prevision d'une
C       suite avec extrapolation, mais
C          - c'est rare
C          - si on demarre un calcul a l'ordre deux a partir d'un calcul
C            a l'ordre 1, on peut estimer que les premiers pas de temps
C            sont a jeter de toute facon.
C       Une exception : on ecrit egalement Cp en effet joule pour
C         pouvoir calculer la temperature H/Cp en debut de calcul
C
          IF(IVIEXT(IPHAS).GT.0) THEN
C         Viscosite moleculaire - cellules (si variable)
            IF(IVIVAR(IPHAS).EQ.1) THEN
              RUBRIQ = 'viscl_ce_phase'//CPHASE(IPHAS)
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPPROC(IVISCL(IPHAS))),IERROR)
              NBERRO = NBERRO+IERROR
            ENDIF
C
C         Viscosite turbulente ou de sous-maille - cellules
            RUBRIQ = 'visct_ce_phase'//CPHASE(IPHAS)
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPPROC(IVISCT(IPHAS))),IERROR)
            NBERRO = NBERRO+IERROR
          ENDIF
C
          IF((ICPEXT(IPHAS ).GT.0.AND.ICP(IPHAS).GT.0).OR.
     &       (IPPMOD(IELJOU).GE.1.AND.ICP(IPHAS).GT.0))  THEN
C         Chaleur massique - cellules
            RUBRIQ = 'cp_ce_phase'//CPHASE(IPHAS)
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPPROC(ICP(IPHAS))),IERROR)
            NBERRO = NBERRO+IERROR
          ENDIF
C
        ENDDO
C
C     Si on a des scalaires, on ecrit leur diffusivite
C        (on ne l'ecrit pas pour les variances)
        IF(NSCAL.GT.0) THEN
          DO ISCAL = 1, NSCAL
            IF(IVSEXT(ISCAL).GT.0.AND.IVISLS(ISCAL).GT.0.AND.
     &         (ISCAVR(ISCAL).LE.0.OR.ISCAVR(ISCAL).GT.NSCAL) ) THEN
C         Diffusivite - cellules
              RUBRIQ = 'visls_ce_scalaire'//CSCAL(ISCAL)
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPPROC(IVISLS(ISCAL))),IERROR)
              NBERRO = NBERRO+IERROR
            ENDIF
          ENDDO
        ENDIF
C
        IF (NBERRO.NE.0) THEN
          CAR54='ERREUR A L''ECRITURE DES PROPRIETES PHYSIQUES         '
          WRITE(NFECRA,8101)CAR54
        ENDIF
C
        CAR54 =' Fin de l''ecriture des proprietes physiques          '
        WRITE(NFECRA,1110)CAR54
C
C ---> Pas de temps
C
        NBERRO = 0
C
        IF(IDTVAR.EQ.2) THEN
          RUBRIQ = 'dt_variable_espace_ce'
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,DT,
     &                IERROR)
          NBERRO=NBERRO+IERROR
        ELSEIF(IDTVAR.EQ.1) THEN
          RUBRIQ = 'dt_variable_temps'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,DT,
     &                IERROR)
          NBERRO=NBERRO+IERROR
        ENDIF
C
        IF (NBERRO.NE.0) THEN
          CAR54='ERREUR A L''ECRITURE DU PAS DE TEMPS                  '
          WRITE(NFECRA,8101)CAR54
        ENDIF
C
        CAR54 =' Fin de l''ecriture du pas de temps                   '
        WRITE(NFECRA,1110)CAR54
C
C ---> Flux de masse
C
C     Pour garder la memoire de la correspondance entre les variables
C     et les flux de masse, on memorise le nom de chaque variable
C     (NOMFLU(I)= nom de la ieme variable)
C     Ensuite, pour chaque variable, on ecrit son nom et le numero
C     local du flux de masse correspondant (en pratique 1 ou 2)
C
        NBERRO=0
C
C       Initialisation des tableaux de travail
        NBFLU = 0
        DO II = 1, NVARMX
          IFAIT(II) = 0
        ENDDO
C
        IPHAS = 1
        NOMFLU(IPR(IPHAS))='fm_p_phase'//CPHASE(IPHAS)
        DO IPHAS=1,NPHAS
          NOMFLU(IU(IPHAS))='fm_u_phase'//CPHASE(IPHAS)
          NOMFLU(IV(IPHAS))='fm_v_phase'//CPHASE(IPHAS)
          NOMFLU(IW(IPHAS))='fm_w_phase'//CPHASE(IPHAS)
          IF (ITYTUR(IPHAS).EQ.2) THEN
            NOMFLU(IK(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
            NOMFLU(IR11(IPHAS))='fm_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IR22(IPHAS))='fm_R22_phase'//CPHASE(IPHAS)
            NOMFLU(IR33(IPHAS))='fm_R33_phase'//CPHASE(IPHAS)
            NOMFLU(IR12(IPHAS))='fm_R12_phase'//CPHASE(IPHAS)
            NOMFLU(IR13(IPHAS))='fm_R13_phase'//CPHASE(IPHAS)
            NOMFLU(IR23(IPHAS))='fm_R23_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50) THEN
            NOMFLU(IK(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_eps_phase'//CPHASE(IPHAS)
            NOMFLU(IPHI(IPHAS))='fm_phi_phase'//CPHASE(IPHAS)
C     On n'utilise pas le flux de masse pour fb en fait mais on le laisse ici, car ca
C     ne change rien (le flux n'est ecrit qu'une seule fois)
            NOMFLU(IFB(IPHAS))='fm_fb_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.60) THEN
            NOMFLU(IK(IPHAS))='fm_k_phase'//CPHASE(IPHAS)
            NOMFLU(IOMG(IPHAS))='fm_omega_phase'//CPHASE(IPHAS)
          ENDIF
        ENDDO
        IF(NSCAL.GT.0) THEN
          DO ISCAL = 1, NSCAL
            NOMFLU(ISCA(ISCAL))='fm_scalaire'//CSCAL(ISCAL)
          ENDDO
        ENDIF
        IF (IALE.EQ.1) THEN
          NOMFLU(IUMA)='fm_vit_maill_u'
          NOMFLU(IVMA)='fm_vit_maill_v'
          NOMFLU(IWMA)='fm_vit_maill_w'
        ENDIF
C
C
        DO IVAR = 1, NVAR
C
C        Si la variable n'est pas associee a un flux de masse
C          on ne fait rien
          IF (IFLUMA(IVAR).GT.0) THEN
C          Si le flux de masse n'a pas encore ete ecrit
            IF (IFAIT(IFLUMA(IVAR)).EQ.0) THEN
C            C'est un nouveau flux a ecrire
              NBFLU=NBFLU+1
C            On note que ce flux est alors traite (c'est le NBFLUieme)
              IFAIT(IFLUMA(IVAR))=NBFLU
C
C            Ecriture du flux de masse sur les faces internes
              RUBRIQ = 'flux_masse_fi_'//CFLU(NBFLU)
              ITYSUP = 2
              NBVAL  = 1
              IRTYP  = 2
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPFA(1,IPPROF(IFLUMA(IVAR))),IERROR)
              NBERRO=NBERRO+IERROR
C
C            Ecriture du flux de masse sur les faces de bord
              RUBRIQ = 'flux_masse_fb_'//CFLU(NBFLU)
              ITYSUP = 3
              NBVAL  = 1
              IRTYP  = 2
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPFB(1,IPPROB(IFLUMA(IVAR))),IERROR)
              NBERRO=NBERRO+IERROR
            ENDIF
C          Que le flux de masse ait deja ete ecrit ou pas,
C            on ecrit  le lien variable/numero de flux
            RUBRIQ = NOMFLU(IVAR)
            ITYSUP = 0
            NBVAL  = 1
            IRTYP  = 1
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  IFAIT(IFLUMA(IVAR)),IERROR)
            NBERRO=NBERRO+IERROR
          ENDIF
        ENDDO
C
C
C     On fait de meme pour les flux de masse au temps precedent.
C       Initialisation des tableaux de travail
C
        NBFLU = 0
        DO II = 1, NVARMX
          IFAIT(II) = 0
        ENDDO
C
        IPHAS = 1
        NOMFLU(IPR(IPHAS))='fm_a_p_phase'//CPHASE(IPHAS)
        DO IPHAS=1,NPHAS
          NOMFLU(IU(IPHAS))='fm_a_u_phase'//CPHASE(IPHAS)
          NOMFLU(IV(IPHAS))='fm_a_v_phase'//CPHASE(IPHAS)
          NOMFLU(IW(IPHAS))='fm_a_w_phase'//CPHASE(IPHAS)
          IF (ITYTUR(IPHAS).EQ.2) THEN
            NOMFLU(IK(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
            NOMFLU(IR11(IPHAS))='fm_a_R11_phase'//CPHASE(IPHAS)
            NOMFLU(IR22(IPHAS))='fm_a_R22_phase'//CPHASE(IPHAS)
            NOMFLU(IR33(IPHAS))='fm_a_R33_phase'//CPHASE(IPHAS)
            NOMFLU(IR12(IPHAS))='fm_a_R12_phase'//CPHASE(IPHAS)
            NOMFLU(IR13(IPHAS))='fm_a_R13_phase'//CPHASE(IPHAS)
            NOMFLU(IR23(IPHAS))='fm_a_R23_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50) THEN
            NOMFLU(IK(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IEP(IPHAS))='fm_a_eps_phase'//CPHASE(IPHAS)
            NOMFLU(IPHI(IPHAS))='fm_a_phi_phase'//CPHASE(IPHAS)
            NOMFLU(IFB(IPHAS))='fm_a_fb_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.60) THEN
            NOMFLU(IK(IPHAS))='fm_a_k_phase'//CPHASE(IPHAS)
            NOMFLU(IOMG(IPHAS))='fm_a_omega_phase'//CPHASE(IPHAS)
          ENDIF
        ENDDO
        IF(NSCAL.GT.0) THEN
          DO ISCAL = 1, NSCAL
            NOMFLU(ISCA(ISCAL))='fm_a_scalaire'//CSCAL(ISCAL)
          ENDDO
        ENDIF
        IF (IALE.EQ.1) THEN
          NOMFLU(IUMA)='fm_a_vit_maill_u'
          NOMFLU(IVMA)='fm_a_vit_maill_v'
          NOMFLU(IWMA)='fm_a_vit_maill_w'
        ENDIF
C
C
        DO IVAR = 1, NVAR
C
C        Si la variable n'est pas associee a un flux de masse
C          on ne fait rien
          IF (IFLUAA(IVAR).GT.0) THEN
C          Si le flux de masse n'a pas encore ete ecrit
            IF (IFAIT(IFLUAA(IVAR)).EQ.0) THEN
C            C'est un nouveau flux a ecrire
              NBFLU=NBFLU+1
C            On note que ce flux est alors traite (c'est le NBFLUieme)
              IFAIT(IFLUAA(IVAR))=NBFLU
C
C            Ecriture du flux de masse sur les faces internes
              RUBRIQ = 'flux_masse_a_fi_'//CFLU(NBFLU)
              ITYSUP = 2
              NBVAL  = 1
              IRTYP  = 2
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPFA(1,IPPROF(IFLUAA(IVAR))),IERROR)
              NBERRO=NBERRO+IERROR
C
C            Ecriture du flux de masse sur les faces de bord
              RUBRIQ = 'flux_masse_a_fb_'//CFLU(NBFLU)
              ITYSUP = 3
              NBVAL  = 1
              IRTYP  = 2
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPFB(1,IPPROB(IFLUAA(IVAR))),IERROR)
              NBERRO=NBERRO+IERROR
            ENDIF
C          Que le flux de masse ait deja ete ecrit ou pas,
C            on ecrit  le lien variable/numero de flux
            RUBRIQ = NOMFLU(IVAR)
            ITYSUP = 0
            NBVAL  = 1
            IRTYP  = 1
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  IFAIT(IFLUAA(IVAR)),IERROR)
            NBERRO=NBERRO+IERROR
          ENDIF
        ENDDO
C
        IF (NBERRO.NE.0) THEN
          CAR54='ERREUR A L''ECRITURE DES FLUX DE MASSE                '
          WRITE(NFECRA,8101)CAR54
        ENDIF
C
        CAR54 =' Fin de l''ecriture des flux de masse                 '
        WRITE(NFECRA,1110)CAR54
C
C ---> Conditions aux limites
C
        NBERRO=0
C
        NOMCLI(IPR(1))='_p_phase'//CPHASE(1)
        DO IPHAS = 1, NPHAS
          NOMCLI(IU(IPHAS))='_u_phase'//CPHASE(IPHAS)
          NOMCLI(IV(IPHAS))='_v_phase'//CPHASE(IPHAS)
          NOMCLI(IW(IPHAS))='_w_phase'//CPHASE(IPHAS)
          IF (ITYTUR(IPHAS).EQ.2) THEN
            NOMCLI(IK(IPHAS))='_k_phase'//CPHASE(IPHAS)
            NOMCLI(IEP(IPHAS))='_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
            NOMCLI(IR11(IPHAS))='_R11_phase'//CPHASE(IPHAS)
            NOMCLI(IR22(IPHAS))='_R22_phase'//CPHASE(IPHAS)
            NOMCLI(IR33(IPHAS))='_R33_phase'//CPHASE(IPHAS)
            NOMCLI(IR12(IPHAS))='_R12_phase'//CPHASE(IPHAS)
            NOMCLI(IR13(IPHAS))='_R13_phase'//CPHASE(IPHAS)
            NOMCLI(IR23(IPHAS))='_R23_phase'//CPHASE(IPHAS)
            NOMCLI(IEP(IPHAS))='_eps_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.50) THEN
            NOMCLI(IK(IPHAS))='_k_phase'//CPHASE(IPHAS)
            NOMCLI(IEP(IPHAS))='_eps_phase'//CPHASE(IPHAS)
            NOMCLI(IPHI(IPHAS))='_phi_phase'//CPHASE(IPHAS)
            NOMCLI(IFB(IPHAS))='_fb_phase'//CPHASE(IPHAS)
          ELSEIF (ITURB(IPHAS).EQ.60) THEN
            NOMCLI(IK(IPHAS))='_k_phase'//CPHASE(IPHAS)
            NOMCLI(IOMG(IPHAS))='_omega_phase'//CPHASE(IPHAS)
          ENDIF
        ENDDO
        IF(NSCAL.GT.0) THEN
          DO ISCAL = 1, NSCAL
            NOMCLI(ISCA(ISCAL))='_scalaire'//CSCAL(ISCAL)
          ENDDO
        ENDIF
        IF (IALE.EQ.1) THEN
          NOMCLI(IUMA)='_vit_maillage_u'
          NOMCLI(IVMA)='_vit_maillage_v'
          NOMCLI(IWMA)='_vit_maillage_w'
        ENDIF
C
C     Dans le cas ou il y a plusieurs phases,
C       on ne veut ecrire la pression qu'une seule fois
C       mais ca tombe bien, car on ne la verra qu'une seule fois
C       dans la liste des variables
C
        DO IVAR = 1, NVAR
C
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 2
C
C          Coefficients numeros 1
          ICLVAR = ICLRTP(IVAR,ICOEF)
          RUBRIQ = 'cla1'//NOMCLI(IVAR)
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                COEFA(1,ICLVAR),IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'clb1'//NOMCLI(IVAR)
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                COEFB(1,ICLVAR),IERROR)
          NBERRO=NBERRO+IERROR
C
C          Coefficients numeros 2
          ICLVAF = ICLRTP(IVAR,ICOEFF)
          IF (ICLVAR.NE.ICLVAF) THEN
C
            RUBRIQ = 'cla2'//NOMCLI(IVAR)
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  COEFA(1,ICLVAF),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'clb2'//NOMCLI(IVAR)
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  COEFB(1,ICLVAF),IERROR)
            NBERRO=NBERRO+IERROR
          ENDIF
C
        ENDDO
C
C
C     Type symétrie (utilisé pour les gradients par moindres carrés
C       sur support étendu, avec extrapolation du gradient au bord).
C
        DO IPHAS = 1, NPHAS
C
          RUBRIQ = 'isympa_fb_phase'//CPHASE(IPHAS)
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IISYMP+NFABOR*(IPHAS-1)),IERROR)
          NBERRO=NBERRO+IERROR
C
        ENDDO
C
C
        IF (NBERRO.NE.0) THEN
          CAR54='ERREUR A L''ECRITURE DES CONDITIONS AUX LIMITES       '
          WRITE(NFECRA,8101)CAR54
        ENDIF
C
        CAR54 =' Fin de l''ecriture des conditions aux limites        '
        WRITE(NFECRA,1110)CAR54
C
C ---> Termes sources
C      Lorsqu'ils sont extrapoles (pour les versions elec, voir plus bas)
C
        NBERRO=0
C
        IECR = 0
C
C ---> Termes sources Navier-Stokes
        DO IPHAS = 1, NPHAS
C
C     Si les termes sont a l'ordre 2
          IF(ISNO2T(IPHAS).GT.0) THEN
C
            IECR = 1
C
            IPTSNA = IPPROC(ITSNSA(IPHAS))
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
C
            RUBRIQ = 'tsource_ns_ce_x_phase'//CPHASE(IPHAS)
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSNA  ),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'tsource_ns_ce_y_phase'//CPHASE(IPHAS)
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSNA+1),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'tsource_ns_ce_z_phase'//CPHASE(IPHAS)
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSNA+2),IERROR)
            NBERRO=NBERRO+IERROR
C
          ENDIF
C
        ENDDO
C
C ---> Termes sources turbulence
        DO IPHAS = 1, NPHAS
C
C        Si les termes sont a l'ordre 2
          IF(ISTO2T(IPHAS).GT.0) THEN
C
            IECR = 1
C
            IPTSTA = IPPROC(ITSTUA(IPHAS))
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
C
C          En k-eps
            IF(ITYTUR(IPHAS).EQ.2) THEN
C
              RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA  ),IERROR)
              NBERRO=NBERRO+IERROR
C
              RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+1),IERROR)
              NBERRO=NBERRO+IERROR
C
C          En Rij
            ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
              RUBRIQ = 'tsource_tu_ce_R11_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA  ),IERROR)
              NBERRO=NBERRO+IERROR
              RUBRIQ = 'tsource_tu_ce_R22_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+1),IERROR)
              NBERRO=NBERRO+IERROR
              RUBRIQ = 'tsource_tu_ce_R33_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+2),IERROR)
              NBERRO=NBERRO+IERROR
              RUBRIQ = 'tsource_tu_ce_R12_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+3),IERROR)
              NBERRO=NBERRO+IERROR
              RUBRIQ = 'tsource_tu_ce_R13_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+4),IERROR)
              NBERRO=NBERRO+IERROR
              RUBRIQ = 'tsource_tu_ce_R23_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+5),IERROR)
              NBERRO=NBERRO+IERROR
C
              RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+6),IERROR)
              NBERRO=NBERRO+IERROR
C
C          En v2f
            ELSEIF(ITURB(IPHAS).EQ.50) THEN
C
              RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA  ),IERROR)
              NBERRO=NBERRO+IERROR
C
              RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+1),IERROR)
              NBERRO=NBERRO+IERROR
C
              RUBRIQ = 'tsource_tu_ce_phi_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+2),IERROR)
              NBERRO=NBERRO+IERROR
C
              RUBRIQ = 'tsource_tu_ce_fb_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+3),IERROR)
              NBERRO=NBERRO+IERROR
C
C          En k-omega
            ELSEIF(ITURB(IPHAS).EQ.60) THEN
C
              RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA  ),IERROR)
              NBERRO=NBERRO+IERROR
C
              RUBRIQ = 'tsource_tu_ce_omega_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IPTSTA+1),IERROR)
              NBERRO=NBERRO+IERROR
C
            ENDIF
          ENDIF
        ENDDO
C
C
C ---> Termes sources scalaires
C
C     Boucle sur les scalaires
        DO ISCAL = 1, NSCAL
C     Si le terme est a l'ordre 2
          IF(ISSO2T(ISCAL).GT.0) THEN
C
            IECR = 1
C
            IPTSCA = IPPROC(ITSSCA(ISCAL))
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
            RUBRIQ = 'tsource_sc_ce_scalaire'//CSCAL(ISCAL)
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  PROPCE(1,IPTSCA  ),IERROR)
            NBERRO=NBERRO+IERROR
C
          ENDIF
        ENDDO
C
C
        IF (NBERRO.NE.0) THEN
          CAR54='ERREUR A L''ECRITURE DES TERMES SOURCES               '
          WRITE(NFECRA,8101)CAR54
        ENDIF
C
        IF (IECR.NE.0) THEN
          CAR54=' Fin de l''ecriture des termes sources                '
          WRITE(NFECRA,1110)CAR54
        ENDIF
C
C
C ---> Moyennes (cumuls)
C
        NBERRO = 0
C
C  ---> Nombre de moyennes
        RUBRIQ = 'nombre_moyennes_temps'
        ITYSUP = 0
        NBVAL  = 1
        IRTYP  = 1
        CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,NBMOMT,
     &              IERROR)
        NBERRO=NBERRO+IERROR
C
C     Cumuls des moyennes
        DO IMOM = 1, NBMOMT
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
          RUBRIQ = 'cumul_ce_moment'//CMOM(IMOM)
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPPROC(ICMOME(IMOM))),IERROR)
          NBERRO=NBERRO+IERROR
        ENDDO
C
C     Cumuls des durees
C
C     On determine un numero unique de duree, local au fichier suite.
C       les durees variable ou non sont considerees sans distinction
C     Le tableau ICDTVU(2*NBMOMX) renvoie aux cumuls temporels variables
C       en espace pour les indices 1->NBMOMX et aux cumuls temporels
C       uniformes pour les indices NBMOMX+1->2*NBMOMX
C
C     Initialisation des elements de travail
        NBCTM = 0
        DO II = 1, NBMOM2
          ICDTVU(II) = 0
        ENDDO
C
C
        DO IMOM = 1, NBMOMT
          IDTM = IDTMOM(IMOM)
C        Si cumul variable en espace
          IF(IDTM.GT.0) THEN
C          Si cumul pas encore vu
            IF(ICDTVU(IDTM).EQ.0) THEN
C            C'est un nouveau, le NBCTM ieme
              NBCTM = NBCTM+1
C            On le marque
              ICDTVU(IDTM) = NBCTM
C            On ecrit ses valeurs
              IF(NBCTM.LE.NFMTMO) THEN
                WRITE(CAR4,'(I4.4)')NBCTM
              ELSE
                CAR4 = CINDFM
              ENDIF
              IDTCM  = IPPROC(ICDTMO(IDTM))
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              RUBRIQ = 'cumul_temps_ce_'//CAR4
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    PROPCE(1,IDTCM),IERROR)
              NBERRO=NBERRO+IERROR
            ENDIF
C          Cumul vu ou pas, on ecrit son numero, >0 pour dire qu'il est
C            variable en espace
            RUBRIQ = 'numero_cumul_temps_moment'//CMOM(IMOM)
            ITYSUP = 0
            NBVAL  = 1
            IRTYP  = 1
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  ICDTVU(IDTM),IERROR)
            NBERRO=NBERRO+IERROR
C
C        Sinon, si cumul uniforme en espace
          ELSEIF(IDTM.LT.0) THEN
C          Si cumul pas encore vu
            IF(ICDTVU(NBMOMX-IDTM).EQ.0) THEN
C            C'est un nouveau, le NBCTM ieme
              NBCTM = NBCTM+1
C            On le marque
              ICDTVU(NBMOMX-IDTM) = -NBCTM
C            On ecrit ses valeurs
              IF(NBCTM.LE.NFMTMO) THEN
                WRITE(CAR4,'(I4.4)')NBCTM
              ELSE
                CAR4 = CINDFM
              ENDIF
              IDTCM  = -IDTM
              ITYSUP = 0
              NBVAL  = 1
              IRTYP  = 2
              RUBRIQ = 'cumul_temps_'//CAR4
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                    DTCMOM(IDTCM),IERROR)
              NBERRO=NBERRO+IERROR
            ENDIF
C          Cumul vu ou pas, on ecrit son numero, <0 pour dire qu'il est
C            uniforme en espace
            RUBRIQ = 'numero_cumul_temps_moment'//CMOM(IMOM)
            ITYSUP = 0
            NBVAL  = 1
            IRTYP  = 1
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  ICDTVU(NBMOMX-IDTM),IERROR)
            NBERRO=NBERRO+IERROR
          ENDIF
        ENDDO
C
        IF (NBERRO.NE.0) THEN
          CAR54='ERREUR A L''ECRITURE DES MOYENNES TEMPORELLES         '
          WRITE(NFECRA,8101)CAR54
        ENDIF
C
        IF(NBMOMT.GT.0) THEN
          CAR54=' Fin de l''ecriture des moyennes temporelles          '
          WRITE(NFECRA,1110)CAR54
        ENDIF
C
C ---> Distance a la paroi
C      On pourra ecrire ici la distance a la paroi
C
        NBERRO = 0
C
        IECR = 0
C
        IF(INEEDY.EQ.1) THEN
C     Ancien mode de calcul. On ecrit aussi la distance a la paroi,
C       au cas ou on fait une suite en ICDPAR=1.
          IF(ABS(ICDPAR).EQ.2) THEN
            DO IPHAS = 1, NPHAS
              IF(IIFAPA(IPHAS).GT.0) THEN
                IECR   = 1
                ITYSUP = 1
                NBVAL  = 1
                IRTYP  = 1
                RUBRIQ = 'num_fac_par_ce_phase'//CPHASE(IPHAS)
                CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &                    IRTYP,IA(IIFAPA(IPHAS)),IERROR)
                NBERRO=NBERRO+IERROR
              ENDIF
            ENDDO
C     Pour la distance reelle, on a besoin d'un tableau provisoire
C     on ne prend que la phase 1
            IPHAS  = 1
            IF(IIFAPA(IPHAS).GT.0) THEN
              IW1    = IDEBRA
              IFINRA = IW1 + NCELET
              CALL RASIZE('ECRAVA',IFINRA)
              DO IEL = 1, NCEL
                IFAC = IA(IIFAPA(IPHAS)+IEL-1)
                RA(IW1+IEL-1) =
     &               SQRT((CDGFBO(1,IFAC)-XYZCEN(1,IEL))**2
     &               +       (CDGFBO(2,IFAC)-XYZCEN(2,IEL))**2
     &               +       (CDGFBO(3,IFAC)-XYZCEN(3,IEL))**2)
              ENDDO
              IECR   = 1
              ITYSUP = 1
              NBVAL  = 1
              IRTYP  = 2
              RUBRIQ = 'dist_fac_par_ce_phase'//CPHASE(IPHAS)
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &             RA(IW1),IERROR)
              NBERRO=NBERRO+IERROR
            ENDIF
C
C     Nouveau mode de calcul
          ELSEIF(ABS(ICDPAR).EQ.1) THEN
            IPHAS  = 1
            IECR   = 1
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
            RUBRIQ = 'dist_fac_par_ce_phase'//CPHASE(IPHAS)
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RA(IDIPAR),IERROR)
            NBERRO=NBERRO+IERROR
          ENDIF
        ENDIF
C
        IF (NBERRO.NE.0) THEN
          CAR54='ERREUR A L''ECRITURE DE LA DISTANCE A LA PAROI        '
          WRITE(NFECRA,8101)CAR54
        ENDIF
C
        IF(IECR.NE.0) THEN
          CAR54=' Fin de l''ecriture de la distance a la paroi         '
          WRITE(NFECRA,1110)CAR54
        ENDIF
C
C ---> Force exterieure
C
        IF(IPHYDR.EQ.1) THEN
          NBERRO=0
C
          DO IPHAS = 1, NPHAS
C
            ITYSUP = 1
            NBVAL  = 1
            IRTYP  = 2
C
            RUBRIQ = 'force_ext_ce_x_phase'//CPHASE(IPHAS)
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  FRCXT(1,1,IPHAS),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'force_ext_ce_y_phase'//CPHASE(IPHAS)
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  FRCXT(1,2,IPHAS),IERROR)
            NBERRO=NBERRO+IERROR
C
            RUBRIQ = 'force_ext_ce_z_phase'//CPHASE(IPHAS)
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  FRCXT(1,3,IPHAS),IERROR)
            NBERRO=NBERRO+IERROR
          ENDDO
C
          IF (NBERRO.NE.0) THEN
            CAR54=
     &          'ERREUR A L''ECRITURE DES FORCES EXTERIEURES           '
            WRITE(NFECRA,8101)CAR54
          ENDIF
C
          CAR54=' Fin de l''ecriture des forces exterieures            '
          WRITE(NFECRA,1110)CAR54
C
        ENDIF
C
C ---> Methode ALE
C
        IF(IALE.EQ.1) THEN
          NBERRO=0
C
          ITYSUP = 4
          NBVAL  = 1
          IRTYP  = 2
C
          RUBRIQ = 'deplact_x_no'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                RA(IDEPAL),IERROR)
          NBERRO=NBERRO+IERROR
          RUBRIQ = 'deplact_y_no'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                RA(IDEPAL+NNOD),IERROR)
          NBERRO=NBERRO+IERROR
          RUBRIQ = 'deplact_z_no'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                RA(IDEPAL+2*NNOD),IERROR)
          NBERRO=NBERRO+IERROR
C
          IF (NBERRO.NE.0) THEN
            CAR54=
     &          'ERREUR A L''ECRITURE DU DEPLACEMENT AUX NOEUDS (ALE)  '
            WRITE(NFECRA,8101)CAR54
          ENDIF
C
C     Viscosite de maillage (elle est souvent definie geometriquement sur le
C       maillage initial ... il est donc plus facile de la relire ensuite)

          NBERRO = 0
          RUBRIQ = 'type_visc_mail'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 1
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &         IORTVM,IERROR)
          NBERRO = NBERRO+IERROR
C
          NBERRO = 0
          RUBRIQ = 'visc_maillage_x'
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &         PROPCE(1,IPPROC(IVISMA(1))),IERROR)
          NBERRO = NBERRO+IERROR
C
          IF (IORTVM.EQ.1) THEN
            RUBRIQ = 'visc_maillage_y'
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &         PROPCE(1,IPPROC(IVISMA(2))),IERROR)
            RUBRIQ = 'visc_maillage_z'
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &         PROPCE(1,IPPROC(IVISMA(3))),IERROR)
          ENDIF
C
          IF (NBERRO.NE.0) THEN
            CAR54=
     &          'ERREUR A L''ECRITURE DE LA VISCOSITE DE MAILLAGE (ALE)'
            WRITE(NFECRA,8101)CAR54
          ENDIF
C
          CAR54=' Fin de l''ecriture des donnees ALE    '
          WRITE(NFECRA,1110)CAR54
C
          NBERRO=0
          RUBRIQ = 'nombre_structures'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 1
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &         NBSTRU,IERROR)
          NBERRO=NBERRO+IERROR
C
          IF (NBSTRU.GT.0) THEN
C
            NFMTST = 99
            CINDST = 'XX'
C     Codage en chaine de caracteres du numero de la structure
            DO ISTR = 1, MIN(NBSTRU ,NFMTST)
              WRITE(CSTRUC(ISTR),'(I2.2)') ISTR
            ENDDO
            DO ISTR = MIN(NBSTRU ,NFMTST)+1,NBSTRU
              CSTRUC(ISTR) = CINDST
            ENDDO

            DO ISTR = 1, NBSTRU
              RUBRIQ = 'donnees_structure_'//CSTRUC(ISTR)
              ITYSUP = 0
              NBVAL  = 27
              IRTYP  = 2
C
              DO II = 1, 3
                TMPSTR(   II) = XSTR  (II,ISTR)
                TMPSTR(3 +II) = XPSTR (II,ISTR)
                TMPSTR(6 +II) = XPPSTR(II,ISTR)
                TMPSTR(9 +II) = XSTA  (II,ISTR)
                TMPSTR(12+II) = XPSTA (II,ISTR)
                TMPSTR(15+II) = XPPSTA(II,ISTR)
                TMPSTR(18+II) = XSTP  (II,ISTR)
                TMPSTR(21+II) = FORSTR(II,ISTR)
                TMPSTR(24+II) = FORSTA(II,ISTR)
              ENDDO
C
              CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &             TMPSTR,IERROR)
              NBERRO=NBERRO+IERROR
            ENDDO
C
            IF (NBERRO.NE.0) THEN
              CAR54=
     &          'ERREUR A L''ECRITURE DES DONNEES DES STRUCTURES (ALE) '
              WRITE(NFECRA,8101)CAR54
            ENDIF
C
            CAR54=' Fin de l''ecriture des donnees des strucutres (ALE)'
            WRITE(NFECRA,1110)CAR54
C
          ENDIF
        ENDIF
C
C
C ---> Grandeurs complementaires pour la combustion gaz
C
C     Modele COD3P :
C     ============
C
        IF ( IPPMOD(ICOD3P).GE.0 ) THEN
C
          NBERRO=0
C
          RUBRIQ = 'hinfue_cod3p'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                HINFUE,IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'hinoxy_cod3p'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                HINOXY,IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'tinfue_cod3p'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TINFUE,IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'tinoxy_cod3p'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TINOXY,IERROR)
          NBERRO=NBERRO+IERROR
C
C       Numero des zones
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'num_zone_fb_cod3p'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IIZFPP), IERROR)
          NBERRO=NBERRO+IERROR
C
C       Entree Fuel (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientfu_zone_bord_cod3p'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTFU, IERROR)
          NBERRO=NBERRO+IERROR
C
C       Entree oxydant (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientox_zone_bord_cod3p'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTOX, IERROR)
          NBERRO=NBERRO+IERROR
C
          IF (NBERRO.NE.0) THEN
            CAR54=
     &          'ERREUR A L''ECRITURE DES INFORMATIONS COMBUSTION COD3P'
            WRITE(NFECRA,8101)CAR54
          ENDIF
C
          CAR54=' Fin de l''ecriture des informations combustion COD3P'
          WRITE(NFECRA,1110)CAR54
C
        ENDIF
C
C      Modele EBU :
C      ==========
C
        IF ( IPPMOD(ICOEBU).GE.0 ) THEN
C
          NBERRO=0
C
          RUBRIQ = 'temperature_gaz_frais_ebu'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TGF,IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'frmel_ebu'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                FRMEL,IERROR)
          NBERRO=NBERRO+IERROR
C
C       Numero des zones
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'num_zone_fb_ebu'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IIZFPP), IERROR)
          NBERRO=NBERRO+IERROR
C
C       Entree Gaz brule(si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientgb_zone_bord_ebu'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTGB, IERROR)
          NBERRO=NBERRO+IERROR
C
C       Entree gaz frais (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientgf_zone_bord_ebu'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTGF, IERROR)
          NBERRO=NBERRO+IERROR
C
C       FMENT (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'fment_zone_bord_ebu'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                FMENT, IERROR)
          NBERRO=NBERRO+IERROR
C
C       TKENT (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'tkent_zone_bord_ebu'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TKENT, IERROR)
          NBERRO=NBERRO+IERROR
C
          IF (NBERRO.NE.0) THEN
            CAR54=
     &          'ERREUR A L''ECRITURE DES INFORMATIONS COMBUSTION EBU'
            WRITE(NFECRA,8101)CAR54
          ENDIF
C
          CAR54=' Fin de l''ecriture des informations combustion EBU '
          WRITE(NFECRA,1110)CAR54
C
        ENDIF
C
C      Modele LWC :
C      ==========
C
        IF ( IPPMOD(ICOLWC).GE.0 ) THEN
C
          NBERRO=0
C
          RUBRIQ = 'fmin_lwc'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                FMIN,IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'fmax_lwc'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                FMAX,IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'hmin_lwc'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                HMIN,IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'hmax_lwc'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                HMAX,IERROR)
          NBERRO=NBERRO+IERROR
C
C       Numero des zones
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'num_zone_fb_lwc'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IIZFPP), IERROR)
          NBERRO=NBERRO+IERROR
C
C       Entree Gaz brule(si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientgb_zone_bord_lwc'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTGB, IERROR)
          NBERRO=NBERRO+IERROR
C
C       Entree gaz frais (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientgf_zone_bord_lwc'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTGF, IERROR)
          NBERRO=NBERRO+IERROR
C
C       FMENT (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'fment_zone_bord_lwc'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                FMENT, IERROR)
          NBERRO=NBERRO+IERROR
C
C       TKENT (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'tkent_zone_bord_lwc'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TKENT, IERROR)
          NBERRO=NBERRO+IERROR
C
          IF (NBERRO.NE.0) THEN
            CAR54=
     &          'ERREUR A L''ECRITURE DES INFORMATIONS COMBUSTION LWC'
            WRITE(NFECRA,8101)CAR54
          ENDIF
C
          CAR54=' Fin de l''ecriture des informations combustion LWC '
          WRITE(NFECRA,1110)CAR54
C
        ENDIF
C
C ---> Grandeurs complementaires pour la combustion CP
C
        IF ( IPPMOD(ICP3PL).GE.0 .OR. IPPMOD(ICPL3C).GE.0 ) THEN
          NBERRO=0
C
C
C     Charbon PuLVerise : masse vol des charbons
C
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 2
          DO ICHA = 1, NCHARB
            IF(ICHA.LE.NFMTCH) THEN
              WRITE(CAR2,'(I2.2)')ICHA
            ELSE
              CAR2 = CINDFC
            ENDIF
            RUBRIQ = 'masse_volumique_charbon'//CAR2
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  RHOCK(ICHA), IERROR)
            NBERRO=NBERRO+IERROR
          ENDDO
C
C
C     Charbon PuLVerise : type de zones de bord, ientat, ientcp, timpat
C       x20, pour le calcul de rho au bord en entree
C
C       Numero des zones
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'num_zone_fb_charbon_pulverise'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IIZFPP), IERROR)
          NBERRO=NBERRO+IERROR
C
C
C       Type entree air ou cp (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientat_zone_bord_charbon_pulverise'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTAT, IERROR)
          NBERRO=NBERRO+IERROR
C
C       ientat et x20 ne servent pas pour le CP couple Lagrangien (cplphy)
          IF ( IPPMOD(ICP3PL).GE.0 ) THEN
C
            ITYSUP = 0
            NBVAL  = NOZPPM
            IRTYP  = 1
            RUBRIQ = 'ientcp_zone_bord_charbon_pulverise'
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  IENTCP, IERROR)
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
                CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &               IRTYP,X20(1,ICLA), IERROR)
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
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TIMPAT, IERROR)
          NBERRO=NBERRO+IERROR
C
C
          IF (NBERRO.NE.0) THEN
            CAR54=
     &          'ERREUR A L''ECRITURE DES INFORMATIONS COMBUSTION CP   '
            WRITE(NFECRA,8101)CAR54
          ENDIF
C
          CAR54=' Fin de l''ecriture des informations combustion CP    '
          WRITE(NFECRA,1110)CAR54
C
        ENDIF
C
C ---> Grandeurs complementaires pour la FUEL
C
        IF ( IPPMOD(ICFUEL).GE.0 ) THEN
          NBERRO=0
C
C
C     Charbon PuLVerise : type de zones de bord, ientat, ientfl, timpat
C       qimpat et qimpfl  pour le calcul de rho au bord en entree
C
C       Numero des zones
          ITYSUP = 3
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'num_zone_fb_fuel'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IA(IIZFPP), IERROR)
          NBERRO=NBERRO+IERROR
C
C
C       Type entree air ou fuel (si ce n'est pas NOZPPM, erreur)
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientat_zone_bord_fuel'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTAT, IERROR)
          NBERRO=NBERRO+IERROR
C
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 1
          RUBRIQ = 'ientfl_zone_bord_fuel'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IENTFL, IERROR)
          NBERRO=NBERRO+IERROR
C
C       Timpat
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'timpat_zone_bord_fuel'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                TIMPAT, IERROR)
          NBERRO=NBERRO+IERROR
C
C       Qimpat
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'qimpat_zone_bord_fuel'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                QIMPAT, IERROR)
          NBERRO=NBERRO+IERROR
C
C       Qimpfl
          ITYSUP = 0
          NBVAL  = NOZPPM
          IRTYP  = 2
          RUBRIQ = 'qimpfl_zone_bord_fuel'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                QIMPFL, IERROR)
          NBERRO=NBERRO+IERROR
C
          IF (NBERRO.NE.0) THEN
            CAR54=
     &          'ERREUR A L''ECRITURE DES INFORMATIONS COMBUSTION FUEL '
            WRITE(NFECRA,8101)CAR54
          ENDIF
C
          CAR54=' Fin de l''ecriture des informations combustion FUEL  '
          WRITE(NFECRA,1110)CAR54
C
        ENDIF
C
C ---> Grandeurs complementaires pour les versions electriques
C
        NBERRO=0
        IECR = 0
C
C     Recalage des CL pot des versions electriques
C
        IF ( IPPMOD(IELJOU).GE.1       ) THEN
          IF(IELCOR.EQ.1) THEN
C
            IECR   = 1
            RUBRIQ = 'coeff_recalage_joule'
            ITYSUP = 0
            NBVAL  = 1
            IRTYP  = 2
C
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  COEJOU,IERROR)
            NBERRO=NBERRO+IERROR
C
          ENDIF
        ENDIF
C
        IF ( IPPMOD(IELARC).GE.1 .OR. IPPMOD(IELJOU).GE.1 ) THEN
          IF(IELCOR.EQ.1) THEN
C
            IECR   = 1
            RUBRIQ = 'ddpot_recalage_arc_elec'
            ITYSUP = 0
            NBVAL  = 1
            IRTYP  = 2
C
            CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                  DPOT  ,IERROR)
            NBERRO=NBERRO+IERROR
C
          ENDIF
        ENDIF
C
C     Termes sources des versions electriques
C
        IF ( IPPMOD(IELJOU).GE.1 .OR.
     &       IPPMOD(IELARC).GE.1 .OR.
     &       IPPMOD(IELION).GE.1       ) THEN
C
          IECR   = 1
          IPCEFJ = IPPROC(IEFJOU)
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
C
          RUBRIQ = 'tsource_sc_ce_joule'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPCEFJ  ),IERROR)
          NBERRO=NBERRO+IERROR
C
        ENDIF
C
        IF( IPPMOD(IELARC).GE.1 ) THEN
C
          IECR   = 1
          IPCLA1 = IPPROC(ILAPLA(1))
          IPCLA2 = IPPROC(ILAPLA(2))
          IPCLA3 = IPPROC(ILAPLA(3))
          ITYSUP = 1
          NBVAL  = 1
          IRTYP  = 2
C
          RUBRIQ = 'tsource_ns_ce_x_laplace'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPCLA1  ),IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'tsource_ns_ce_y_laplace'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPCLA2  ),IERROR)
          NBERRO=NBERRO+IERROR
C
          RUBRIQ = 'tsource_ns_ce_z_laplace'
          CALL ECRSUI(IMPAVX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                PROPCE(1,IPCLA3  ),IERROR)
          NBERRO=NBERRO+IERROR
C
        ENDIF
C
        IF (NBERRO.NE.0) THEN
          CAR54='ERREUR A L''ECRITURE DES INFORMATIONS ELECTRIQUES     '
          WRITE(NFECRA,8101)CAR54
        ENDIF
C
        IF (IECR.NE.0) THEN
          CAR54=' Fin de l''ecriture des informations electriques      '
          WRITE(NFECRA,1110)CAR54
        ENDIF
C
C
C       Fermeture du fichiers suite auxiliaire
        CALL CLSSUI(IMPAVX,IERROR)
C
        IF (IERROR.NE.0) THEN
          WRITE(NFECRA,8011) FICAVX
        ENDIF
C
        WRITE(NFECRA,1200)
C
      ENDIF
C     Fin de l'ecriture du fichier suite auxiliaire

C
      RETURN
C
C=======================================================================
C 6. FORMATS
C=======================================================================
C
 1000 FORMAT(3X,'** ECRITURE DU FICHIER SUITE PRINCIPAL',/,
     &       3X,'   ----------------------------------- ',/)
 1010 FORMAT(3X,'   Ouverture d un fichier binaire',/)
 1020 FORMAT(3X,'   Ouverture d un fichier ASCII  ',/)
 1100 FORMAT(' Debut de l''ecriture                                   ')
 1110 FORMAT('  ',A54                                                  )
 1200 FORMAT(' Fin de l''ecriture                                     ')
 2000 FORMAT(/,3X,'** ECRITURE DU FICHIER SUITE AUXILIAIRE',/,
     &         3X,'   ------------------------------------ ',/)
C
 7000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,
     &'@    *********                                               ',/,
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
     &'@    Voir le sous-programme ecrava.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,
     &'@    *********                                               ',/,
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
     &'@    Voir le sous-programme ecrava.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,
     &'@    *********                                               ',/,
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
     &'@    Voir le sous-programme ecrava.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7003 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      Le nombre de moments       max NBMOMX supporte par le ',/,
     &'@        format d''ecriture du fichier suite est             ',/,
     &'@        NFMTMO = ',I10                                       ,/,
     &'@      On a ici un nombre de moments   maximal superieur     ',/,
     &'@        NBMOMX = ',I10                                       ,/,
     &'@      On ne pourra pas relire les moments   dont le numero  ',/,
     &'@        est superieur                                       ',/,
     &'@                                                            ',/,
     &'@    Le calcul sera execute.                                 ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme ecrava.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7004 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      Le nombre de charbons      max NCHARM supporte par le ',/,
     &'@        format d''ecriture du fichier suite est             ',/,
     &'@        NFMTCH = ',I10                                       ,/,
     &'@      On a ici un nombre de charbons  maximal superieur     ',/,
     &'@        NCHARM = ',I10                                       ,/,
     &'@      On ne pourra pas relire certaines informations        ',/,
     &'@        realtives aux charbons dont le numero               ',/,
     &'@        est superieur                                       ',/,
     &'@                                                            ',/,
     &'@    Le calcul sera execute.                                 ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme ecrava.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7005 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@      Le nombre de classes par charbon max NCPCMX supporte  ',/,
     &'@        par le format d''ecriture du fichier suite est      ',/,
     &'@        NFMTCL = ',I10                                       ,/,
     &'@      On a ici un nombre de classes par charbon superieur   ',/,
     &'@        NCPCMX = ',I10                                       ,/,
     &'@      On ne pourra pas relire certaines informations        ',/,
     &'@        realtives aux classes  dont le numero               ',/,
     &'@        est superieur                                       ',/,
     &'@                                                            ',/,
     &'@    Le calcul sera execute.                                 ',/,
     &'@                                                            ',/,
     &'@    Voir le sous-programme ecrava.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 8000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A L''OUVERTURE DU FICHIER SUITE      ',/,
     &'@    *********                                 AVAL PRINCIPAL',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier ',A13,'peut etre                ',/,
     &'@            cree dans le repertoire de travail.             ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A L''OUVERTURE DU FICHIER SUITE      ',/,
     &'@    *********                                AVAL AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier ',A13,'peut etre                ',/,
     &'@            cree dans le repertoire de travail.             ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 8010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,
     &'@    *********                                 AVAL PRINCIPAL',/,
     &'@                                                            ',/,
     &'@    Probleme sur le fichier de nom (',A13,')                ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,
     &'@    *********                                AVAL AUXILIAIRE',/,
     &'@                                                            ',/,
     &'@    Probleme sur le fichier de nom (',A13,')                ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 8100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''ECRITURE DU FICHIER SUITE              ',/,
     &'@    *********                                 AVAL PRINCIPAL',/,
     &'@      ', A54                                                 ,/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8101 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''ECRITURE DU FICHIER SUITE              ',/,
     &'@    *********                                AVAL AUXILIAIRE',/,
     &'@      ', A54                                                 ,/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
      END
c@z
