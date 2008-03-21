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
                        SUBROUTINE COU1DO
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCP    , NFPT1D ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IENTHA , IFPT1D , ICLT1D ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   TPPT1D , TEPT1D , HEPT1D , FEPT1D ,
     &   XLMBT1 , RCPT1D , DTPT1D , DT     , RTPA   ,
     &   PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   CPCST  , CP     , HBORD  , TBORD  ,
     &   RDEVEL , RTUSER , RA     )
C
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC ECRITURE DE DONNEES RELATIVES A UN COUPLAGE AVEC SYRTHES
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCP          ! E  !  -> ! DIMENSION DE CP (NCELET OU 1)        !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NFPT1D       ! E  !  -> ! NOMBRE DE FACES AVEC MODULE THERM 1D !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IENTHA       ! E  !  -> ! 1 SI TPAROI EST UNE ENTHALPIE        !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFPT1D       ! TE !  -> ! NUMERO DE LA FACE EN TRAITEMENT      !
CARGU !              !    !     ! THERMIQUE EN PAROI                   !
CARGU ! ICLT1D       ! TE !  -> ! TYPE DE CONDITION LIMITE             !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! CPCST        ! R  !  -> ! CHALEUR SPECIFIQUE SI CONSTANTE      !
CARGU ! CP(NCP)      ! TR !  -> ! CHALEUR SPECIFIQUE SI VARIABLE       !
CARGU ! HBORD        ! TR !  -> ! COEFFICIENTS D'ECHANGE AUX BORDS     !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! TBORD        ! TR !  -> ! TEMPERATURES AUX BORDS               !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! TPPT1D       ! TR !  -> ! TEMPERATURE DE PAROI                 !
CARGU ! TEPT1D       ! TR !  -> ! TEMPERATURE EXTERIEURE               !
CARGU ! HEPT1D       ! TR !  -> ! COEFFICIENT D'ECHANGE EXTERIEUR      !
CARGU ! FEPT1D       ! TR !  -> ! FLUX EXTERIEUR                       !
CARGU ! XLMBT1       ! TR !  -> ! DIFFUSIVITE THERMIQUE                !
CARGU ! RCPT1D       ! TR !  -> ! ROCP                                 !
CARGU ! DTPT1D       ! TR !  -> ! PAS DE TEMPS                         !
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
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
C
C***********************************************************************
C
C     ARGUMENTS
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFPT1D
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , NCP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IFPT1D(NFPT1D), ICLT1D(NFPT1D)
      INTEGER          IENTHA
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION HBORD(NFABOR),TBORD(NFABOR)
      DOUBLE PRECISION CPCST, CP(NCP)
      DOUBLE PRECISION TPPT1D(NFPT1D)
      DOUBLE PRECISION TEPT1D(NFPT1D), HEPT1D(NFPT1D), FEPT1D(NFPT1D)
      DOUBLE PRECISION XLMBT1(NFPT1D), RCPT1D(NFPT1D), DTPT1D(NFPT1D)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C     VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, MODE
      INTEGER          IPHAS  , IAPPEL
      INTEGER          IFAC, IEL , II
      INTEGER          MAXELT, ILS, IDBIA1
      DOUBLE PRECISION ENTHAL, TEMPER

C***********************************************************************
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C     SI ENTHALPIE, ON TRANSFORME EN TEMPERATURE
C     Il est necessaire de transmettre a SYRTHES des Temperatures
C     Afin de conserver le flux Phi = (lambda/d     ) Delta T
C     ou Phi = (lambda/(d Cp)) Delta H
C     on multiplie HBORD = lambda/(d Cp) par Cp pris dans la
C     cellule adjacente.
C     Le resultat n'est pas garanti (conservation en particulier),
C     on ajoute donc un avertissement.
C
C     On ne change les TBORD et HBORD que sur les faces couplees. Du coup ces
C     tableaux contiennent des choses differentes suivant les faces.
C     C'est dangereux mais pas trop grave car on les jette juste apres
C     (COUPBO passe avant).
C
      IF(IENTHA.EQ.1) THEN
         WRITE(NFECRA,1000)
         MODE = 1
         DO II = 1, NFPT1D
            IFAC = IFPT1D(II)
            IEL  = IFABOR(IFAC)
            ENTHAL = TBORD(IFAC)
            CALL USTHHT (MODE   , ENTHAL , TEMPER  )
C           ===========
            TBORD(IFAC) = TEMPER
            IF(NCP.EQ.NCELET) THEN
               HBORD(IFAC) = HBORD(IFAC)*CP(IEL)
            ELSE
               HBORD(IFAC) = HBORD(IFAC)*CPCST
            ENDIF
         ENDDO
      ENDIF
C
C     Pour l'instant on bloque le couplage si la variable est l'energie
C     -> on pourra le calquer sur coupbo si necessaire.
      IF (IENTHA.EQ.2) THEN
        WRITE(NFECRA,2000)
        CALL CSEXIT(1)
      ENDIF
C
C     Mise a jour des conditions aux limites externes du module 1D
      IPHAS = 1
      IAPPEL = 3
C
      MAXELT = MAX(NCELET,NFAC,NFABOR)
      ILS    = IDEBIA
      IDBIA1 = ILS + MAXELT
      CALL IASIZE('COU1DO',IDBIA1)
C
      CALL  USPT1D
C     ============
     & ( IDBIA1 , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NFPT1D , IPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFPT1D , IA(IDEBIA), ICLT1D ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   TPPT1D , RA(IDEBRA), RA(IDEBRA),
     &   TEPT1D , HEPT1D , FEPT1D ,
     &   XLMBT1 , RCPT1D , DTPT1D ,
     &   DT     , RTPA   ,
     &   PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
C
      IAPPEL = 3
      CALL VERT1D
C     ===========
     & (IDEBIA     , IDEBRA     ,
     &  NFABOR     , NFPT1D     , IAPPEL    ,
     &  IFPT1D     , IA(IDEBIA) , ICLT1D    , IA     ,
     &  RA(IDEBRA) , RA(IDEBRA) ,
     &  XLMBT1     , RCPT1D     , DTPT1D    , RA      )
C
      DO II = 1, NFPT1D
C
         IFAC = IFPT1D(II)
C
         CALL TPAR1D
C        ===========
     & ( II-1      , ICLT1D(II), TBORD(IFAC), HBORD(IFAC),
     &   TEPT1D(II), HEPT1D(II), FEPT1D(II) , XLMBT1(II) ,
     &   RCPT1D(II), DTPT1D(II), TPPT1D(II) )
C
C
      ENDDO
C
C --> FORMATS
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION : COUPLAGE SYRTHES/MODULE 1D AVEC CALCUL EN   ',/,
     &'@                ENTHALPIE                                   ',/,
     &'@    *********                                               ',/,
     &'@      OPTION NON VALIDEE - CONTACTER L''EQUIPE DE DVPT      ',/,
     &'@  ')
 2000 FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION : COUPLAGE SYRTHES/MODULE 1D AVEC CALCUL EN   ',/,
     &'@                ENERGIE                                   ',/,
     &'@    *********                                               ',/,
     &'@      OPTION NON PERMISE - CONTACTER L''EQUIPE DE DVPT      ',/,
     &'@                                                            ',/,
     &'@      Le calcul s''arrete                                   ',/,
     &'@  ')
C
      RETURN
      END
c@z
