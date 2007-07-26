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
                        SUBROUTINE ECRLIS
C                       *****************
C     ---------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NVAR   , NPHAS  , NDIM   , NCELET , NCEL   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IRTP   ,
     &   IDEVEL , ITUSER , IA     ,
     &   RTP    , RTPA   , DT     , VOLUME , XYZCEN ,
     &   RDEVEL , RTUSER , RA     )
C     ---------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE D'ECRITURE DES INFOS LISTING
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
CARGU ! NVAR         ! E  !  -> ! NOMBRE DE VARIABLES                  !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IRTP         ! E  !  -> ! INDICE DE RTP DANS RA                !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! RTP          ! TR !  -> ! TABLEAUX DES VARIABLES AU PDT COURANT!
CARGU ! (NCELET,NVAR)!    !     !                                      !
CARGU ! RTPA         ! TR !  -> ! TABLEAUX DES VARIABLES AU PDT PREC   !
CARGU ! (NCELET,NVAR)!    !     !                                      !
CARGU ! DT   (NCELET)! TR !  -> ! VALEUR DU PAS DE TEMPS               !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET !    !     !                                      !
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
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
      INCLUDE "cstphy.h"
      INCLUDE "albase.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
      INTEGER          IDBIA0, IDBRA0, NVAR, NPHAS, NDIM, NCELET, NCEL
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IRTP
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
      DOUBLE PRECISION RTPA(NCELET,NVAR), RTP(NCELET,NVAR)
      DOUBLE PRECISION DT(NCELET), VOLUME(NCELET)
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ, IC, ICEL, IPP, IRA, IVRTP, IOK
      INTEGER          IPHAS, KPHAS, IPRNEW, IPUVW
      INTEGER          ICMIN, ICMAX
      INTEGER          NBRVAL
      INTEGER          IDIVDT, IXMSDT, IDEBIA, IDEBRA, IFINRA, IEL
      DOUBLE PRECISION PETIT,XYZMIN(3),XYZMAX(3)
      CHARACTER*200    CHAIN, CHAINC
C
C
C=======================================================================
C 0. INITIALISATIONS LOCALES
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      PETIT  =-GRAND
C
C========================================
C 1. CALCUL DES MIN ET MAX DES VARIABLES
C========================================
C
      DO IPP = 2, NVPPMX
        IF(ILISVR(IPP).EQ.1) THEN
          VARMIN(IPP) = GRAND
          VARMAX(IPP) = PETIT
          IRA = ABS(IPP2RA(IPP))
C
C     Pour les moments, il faut eventuellement diviser par le temps cumule
          IDIVDT = IPPMOM(IPP)
          IF(IDIVDT.EQ.0) THEN
            IXMSDT = IRA
          ELSE
            IXMSDT = IDEBRA
            IFINRA = IXMSDT + NCEL
            CALL RASIZE ('ECRLIS', IFINRA)
C           ===========
          ENDIF
          IF(IDIVDT.GT.0) THEN
            DO IEL = 1, NCEL
              RA(IXMSDT+IEL-1) = RA(IRA+IEL-1)/
     &             MAX(RA(IDIVDT+IEL-1),EPZERO)
            ENDDO
          ELSEIF(IDIVDT.LT.0) THEN
            DO IEL = 1, NCEL
              RA(IXMSDT+IEL-1) = RA(IRA+IEL-1)/
     &             MAX(DTCMOM(-IDIVDT),EPZERO)
            ENDDO
C         ELSE
C           RA(IXMSDT+IEL-1) = RA(IRA+IEL-1)
C           inutile car on a pose IXMSDT = IRA
          ENDIF

          DO ICEL = 1, NCEL
            IF(RA(IXMSDT+ICEL-1).LT.VARMIN(IPP))
     &         VARMIN(IPP) = RA(IXMSDT+ICEL-1)
            IF(RA(IXMSDT+ICEL-1).GT.VARMAX(IPP))
     &         VARMAX(IPP) = RA(IXMSDT+ICEL-1)
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMIN (VARMIN(IPP))
C           ===========
            CALL PARMAX (VARMAX(IPP))
C           ===========
          ENDIF
        ENDIF
      ENDDO
C
C==========================================================
C 2. DERIVE POUR LES VARIABLES TRANSPORTEES (sauf pression)
C==========================================================
C
      DO IPP = 2, NVPPMX
        IOK = 1
        DO IPHAS = 1, NPHAS
          IF(IPP.EQ.IPPRTP(IPR(IPHAS))) THEN
            IOK = 0
          ENDIF
        ENDDO
        IF(ILISVR(IPP).EQ.1.AND.ITRSVR(IPP).GE.1) THEN
          IF(IOK.EQ.1) THEN
            IRA = ABS(IPP2RA(IPP))
            IVRTP = (IRA-IRTP)/NCELET+1
            DERVAR(IPP) = 0.D0
            DO ICEL = 1, NCEL
              DERVAR(IPP) = DERVAR(IPP)
     &                + (RTP(ICEL,IVRTP)-RTPA(ICEL,IVRTP))**2
     &                *  VOLUME(ICEL)/DT(ICEL)
            ENDDO
            IF(IRANGP.GE.0) CALL PARSOM (DERVAR(IPP))
C           ===========
            DERVAR(IPP) = DERVAR(IPP) / VOLTOT
          ENDIF
        ENDIF
      ENDDO
C
      DO IPHAS = 1, NPHAS
        IPRNEW = 1
        IF(IPHAS.GT.1) THEN
          DO KPHAS = 1, IPHAS-1
            IF(IPR(IPHAS).EQ.IPR(KPHAS)) THEN
              IPRNEW = 0
            ENDIF
          ENDDO
        ENDIF
        IF(IPRNEW.EQ.1) THEN
          IPP = IPPRTP(IPR(IPHAS))
          IF(DERVAR(IPP).LT.EPZERO) THEN
            DERVAR(IPP) = -1.D0
          ENDIF
          DERVAR(IPP) = RNSMBR(IPP) / DERVAR(IPP)
        ENDIF
      ENDDO

C
C==========================================================
C 3. MIN MAX DT + LOCALISATION
C==========================================================
C
      IF (IDTVAR.GE.0) THEN
C
C     -> min max dt
        PTPLOC(7,1) = GRAND
        PTPLOC(8,1) = PETIT
        ICMIN = 1
        ICMAX = 1
        DO ICEL = 1, NCEL
          IF(DT(ICEL).LT.PTPLOC(7,1)) THEN
            PTPLOC(7,1) = DT(ICEL)
            ICMIN = ICEL
          ENDIF
          IF(DT(ICEL).GT.PTPLOC(8,1)) THEN
            PTPLOC(8,1) = DT(ICEL)
            ICMAX = ICEL
          ENDIF
        ENDDO
C
        XYZMIN(1) = XYZCEN(1,ICMIN)
        XYZMIN(2) = XYZCEN(2,ICMIN)
        XYZMIN(3) = XYZCEN(3,ICMIN)
        XYZMAX(1) = XYZCEN(1,ICMAX)
        XYZMAX(2) = XYZCEN(2,ICMAX)
        XYZMAX(3) = XYZCEN(3,ICMAX)
C
        IF (IRANGP.GE.0) THEN
          NBRVAL = 3
          CALL PARMNL (NBRVAL, PTPLOC(7,1), XYZMIN)
C         ===========
          CALL PARMXL (NBRVAL, PTPLOC(8,1), XYZMAX)
C         ===========
        ENDIF
C
        PTPLOC(7,2) = XYZMIN(1)
        PTPLOC(7,3) = XYZMIN(2)
        PTPLOC(7,4) = XYZMIN(3)
        PTPLOC(8,2) = XYZMAX(1)
        PTPLOC(8,3) = XYZMAX(2)
        PTPLOC(8,4) = XYZMAX(3)
C
      ENDIF
C
C
C=======================================================================
C 4. ECRITURE DES CRITERES DE CONVERGENCE
C=======================================================================
C
C
      WRITE(NFECRA,1000)
      WRITE(NFECRA,1010)
      WRITE(NFECRA,1011)
      WRITE(NFECRA,1010)
C
      DO IPP = 2, NVPPMX
        IF(ILISVR(IPP).EQ.1.AND.ITRSVR(IPP).GE.1) THEN
C
          CHAINC = 'c'
          CHAIN = ' '
          IC = 4
          CHAIN = NOMVAR(IPP)
          CHAINC(IC:IC+12) = CHAIN(1:12)
          IC=IC+12
          CHAIN = ' '
          WRITE(CHAIN,3000) RNSMBR(IPP)
          CHAINC(IC:IC+12) = CHAIN(1:12)
          IC=IC+13
          CHAIN = ' '
          WRITE(CHAIN,4000) NBIVAR(IPP)
          CHAINC(IC:IC+7) = CHAIN(1:7)
          IC=IC+9
          CHAIN = ' '
          WRITE(CHAIN,3000) RESVAR(IPP)
          CHAINC(IC:IC+12) = CHAIN(1:12)
          IC=IC+14
          CHAIN = ' '
          WRITE(CHAIN,3000) DERVAR(IPP)
          CHAINC(IC:IC+12) = CHAIN(1:12)
          IC=IC+12
C
          WRITE(NFECRA,'(A)') CHAINC(1:IC)
C
        ENDIF
      ENDDO
C
      WRITE(NFECRA,1010)
      WRITE(NFECRA,*) ' '
      WRITE(NFECRA,*) ' '
C
C
C=======================================================================
C 5.  SUIVI DES VARIABLES , MIN - MAX - CLIPPING
C=======================================================================
C
C
      WRITE(NFECRA,1100)
      WRITE(NFECRA,1110)
      WRITE(NFECRA,1111)
      WRITE(NFECRA,1110)
C
      DO IPP = 2, NVPPMX
        IF(ILISVR(IPP).EQ.1) THEN
C
          CHAINC = 'v'
          CHAIN = ' '
          IC = 4
          CHAIN = NOMVAR(IPP)
          CHAINC(IC:IC+12) = CHAIN(1:12)
          IC=IC+12
          CHAIN = ' '
          WRITE(CHAIN,3000) VARMIN(IPP)
          CHAINC(IC:IC+12) = CHAIN(1:12)
          IC=IC+14
          CHAIN = ' '
          WRITE(CHAIN,3000) VARMAX(IPP)
          CHAINC(IC:IC+12) = CHAIN(1:12)
          IC=IC+16
          IPUVW = 0
          DO IPHAS = 1, NPHAS
            IF(IPP.EQ.IPPRTP(IPR(IPHAS)) .OR.
     &         IPP.EQ.IPPRTP(IU (IPHAS)) .OR.
     &         IPP.EQ.IPPRTP(IV (IPHAS)) .OR.
     &         IPP.EQ.IPPRTP(IW (IPHAS)) ) THEN
              IPUVW = 1
            ENDIF
C   En v2f on ne clippe jamais f_barrre, on ne l'affiche donc pas
            IF (ITURB(IPHAS).EQ.50) THEN
              IF (IPP.EQ.IPPRTP(IFB(IPHAS))) IPUVW = 1
            ENDIF
          ENDDO
C   En ALE on ne clippe pas la vitesse de maillage
          IF (IALE.EQ.1) THEN
            IF (IPP.EQ.IPPRTP(IUMA) .OR.
     &          IPP.EQ.IPPRTP(IVMA) .OR.
     &          IPP.EQ.IPPRTP(IWMA)) IPUVW = 1
          ENDIF
          IF(IPUVW.EQ.1) THEN
            CHAIN = '     --         --'
            CHAINC(IC:IC+18) = CHAIN(1:18)
            IC = IC+18
          ELSE
            CHAIN = ' '
            WRITE(CHAIN,4000) ICLPMN(IPP)
            CHAINC(IC:IC+7) = CHAIN(1:7)
            IC=IC+11
            CHAIN = ' '
            WRITE(CHAIN,4000) ICLPMX(IPP)
            CHAINC(IC:IC+7) = CHAIN(1:7)
            IC=IC+7
          ENDIF
C
          WRITE(NFECRA,'(A)') CHAINC(1:IC)
CMO          ICLPMN(IPP) = 0
CMO          ICLPMX(IPP) = 0
        ENDIF
      ENDDO
C
      WRITE(NFECRA,1110)
      WRITE(NFECRA,*) ' '
      WRITE(NFECRA,*) ' '
C
C
C
C
      WRITE(NFECRA,1150)
      WRITE(NFECRA,1160)
      WRITE(NFECRA,1161)
      WRITE(NFECRA,1160)
C
      DO IPP = 2, NVPPMX
        IPUVW = 0
        DO IPHAS = 1, NPHAS
          IF(IPP.EQ.IPPRTP(IPR(IPHAS)) .OR.
     &       IPP.EQ.IPPRTP(IU (IPHAS)) .OR.
     &       IPP.EQ.IPPRTP(IV (IPHAS)) .OR.
     &       IPP.EQ.IPPRTP(IW (IPHAS)) ) THEN
            IPUVW = 1
          ENDIF
C   En v2f on ne clippe jamais f_barrre, on ne l'affiche donc pas
          IF (ITURB(IPHAS).EQ.50) THEN
            IF (IPP.EQ.IPPRTP(IFB(IPHAS))) IPUVW = 1
          ENDIF
C   En ALE on ne clippe pas la vitesse de maillage
          IF (IALE.EQ.1) THEN
            IF (IPP.EQ.IPPRTP(IUMA) .OR.
     &          IPP.EQ.IPPRTP(IVMA) .OR.
     &          IPP.EQ.IPPRTP(IWMA)) IPUVW = 1
          ENDIF
C   Compressible
          IF(IPPMOD(ICOMPF).GE.0) THEN
            IF(IPP.EQ.IPPRTP(ISCA(ITEMPK(IPHAS)))) THEN
              IPUVW = 1
            ENDIF
          ENDIF
        ENDDO
        IF(ILISVR(IPP).EQ.1.AND.ITRSVR(IPP).GT.0.AND.IPUVW.EQ.0) THEN
          CHAINC = 'a'
          CHAIN = ' '
          IC = 4
          CHAIN = NOMVAR(IPP)
          CHAINC(IC:IC+12) = CHAIN(1:12)
          IC=IC+12
          CHAIN = ' '
          WRITE(CHAIN,3000) VARMNA(IPP)
          CHAINC(IC:IC+12) = CHAIN(1:12)
          IC=IC+14
          CHAIN = ' '
          WRITE(CHAIN,3000) VARMXA(IPP)
          CHAINC(IC:IC+12) = CHAIN(1:12)
          IC=IC+16
C
          CHAIN = ' '
          WRITE(CHAIN,4000) ICLPMN(IPP)
          CHAINC(IC:IC+7) = CHAIN(1:7)
          IC=IC+11
          CHAIN = ' '
          WRITE(CHAIN,4000) ICLPMX(IPP)
          CHAINC(IC:IC+7) = CHAIN(1:7)
          IC=IC+7
C
          WRITE(NFECRA,'(A)') CHAINC(1:IC)
          ICLPMN(IPP) = 0
          ICLPMX(IPP) = 0
        ENDIF
      ENDDO
C
      WRITE(NFECRA,1160)
      WRITE(NFECRA,*) ' '
      WRITE(NFECRA,*) ' '
C
C
C=======================================================================
C 6.  PAS DE TEMPS LOCAL
C=======================================================================
C
      IF (IDTVAR.GE.0) THEN
C
        WRITE(NFECRA,1200)
        WRITE(NFECRA,1210)
        WRITE(NFECRA,1211)
        WRITE(NFECRA,1210)
        DO II = 1, 6
          CHAINC = ' '
          IF(II.EQ.1) CHAINC = 'Courant min'
          IF(II.EQ.2) CHAINC = 'Courant max'
          IF(II.EQ.3) CHAINC = 'Fourier min'
          IF(II.EQ.4) CHAINC = 'Fourier max'
          IF(II.EQ.5) CHAINC = 'Cou_Fou min'
          IF(II.EQ.6) CHAINC = 'Cou_Fou max'
C    Compressible
          IF(IPPMOD(ICOMPF).GE.0) THEN
            IF(II.EQ.5) CHAINC = 'CFL/Mas min'
            IF(II.EQ.6) CHAINC = 'CFL/Mas max'
          ENDIF
C
          IC = 13
          DO JJ = 1, 4
            CHAIN = ' '
            WRITE(CHAIN,3000) PTPLOC(II,JJ)
            CHAINC(IC:IC+13) = CHAIN(1:13)
            IC = IC+13
          ENDDO
          WRITE(NFECRA,'(A)') CHAINC(1:IC)
        ENDDO
        WRITE(NFECRA,1210)
C
        IF(IDTVAR.GT.0) THEN
          DO II = 7, 8
            CHAINC = ' '
            IF(II.EQ.7) CHAINC = 'Dt min'
            IF(II.EQ.8) CHAINC = 'Dt max'
            IC = 13
            DO JJ = 1, 4
              CHAIN = ' '
              WRITE(CHAIN,3000) PTPLOC(II,JJ)
              CHAINC(IC:IC+13) = CHAIN(1:13)
              IC = IC+13
            ENDDO
            WRITE(NFECRA,'(A)') CHAINC(1:IC)
          ENDDO
          WRITE(NFECRA,1210)
        ENDIF
C
        IF (IPTLRO.EQ.1) THEN
          CHAINC = ' '
          CHAINC = 'Dt/Dtrho max'
          IC = 13
          DO JJ = 1, 4
            CHAIN = ' '
            WRITE(CHAIN,3000) RPDTRO(JJ)
            CHAINC(IC:IC+13) = CHAIN(1:13)
            IC = IC+13
          ENDDO
          WRITE(NFECRA,'(A)') CHAINC(1:IC)
          CHAINC = ' '
          CHAINC = 'Clips a Dtrho'
          IC = 57
          CHAIN = ' '
          WRITE(CHAIN,4000) NCLPTR
          CHAINC(IC:IC+7) = CHAIN(1:7)
          IC = IC+7
          WRITE(NFECRA,'(A)') CHAINC(1:IC)
          WRITE(NFECRA,1210)
        ENDIF
C
        WRITE(NFECRA,*) ' '
        WRITE(NFECRA,*) ' '
C
      ENDIF
C
C=======================================================================
C 3.
C=======================================================================
C--------
C FORMATS
C--------
C
 1000 FORMAT (/,3X,'** INFORMATIONS SUR LA CONVERGENCE',/,
     &          3X,'   -------------------------------')
 1011 FORMAT ('   Variable    Norm 2nd mb.',
     &        '  Nbiter  Residu norme        derive')
 1010 FORMAT ('---------------------------',
     &        '------------------------------------')
C
C
 1100 FORMAT (/,3X,'** INFORMATIONS SUR LES VARIABLES',/,
     &         3X,'   ------------------------------')
 1111 FORMAT ('   Variable      Valeur min    Valeur max',
     &        '   Clip min   Clip max')
 1110 FORMAT ('-----------------------------------------',
     &        '----------------------')
C
 1150 FORMAT (/,3X,'** INFORMATIONS SUR LES CLIPPINGS',/,
     &         3X,'   ------------------------------')
 1161 FORMAT ('   Variable    Min ss clips  Max ss clips',
     &        '   Clip min   Clip max')
 1160 FORMAT ('-----------------------------------------',
     &        '----------------------')
C
C
 1200 FORMAT (/,3X,'** INFORMATIONS SUR LE PAS DE TEMPS',/,
     &         3X,'   --------------------------------')
 1210 FORMAT ('---------------------------',
     &        '------------------------------------')
 1211 FORMAT ('Critere      Valeur     en         xc',
     &        '           yc           zc')
C
 3000 FORMAT (E12.5)
 4000 FORMAT (I7)
C
      RETURN
      END
c@z
