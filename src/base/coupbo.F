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
                        SUBROUTINE COUPBO
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , ISVTB  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NCP , NCV , IENTHA ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   CPCST  , CP     , CVCST  , CV     ,
     &   HBORD  , TBORD  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
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
CARGU ! NCV          ! E  !  -> ! DIMENSION DE CV (NCELET OU 1)        !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IENTHA       ! E  !  -> ! 1 SI TPAROI EST UNE ENTHALPIE        !
CARGU !              ! E  !  -> ! 2 SI TPAROI EST UNE ENERGIE          !
CARGU !              !    !     !    (Compressible)                    !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! CPCST        ! R  !  -> ! CHALEUR SPECIFIQUE SI CONSTANTE      !
CARGU ! CVCST        ! R  !  -> ! CHALEUR SPECIFIQUE SI CONSTANTE      !
CARGU ! CP(NCP)      ! TR !  -> ! CHALEUR SPECIFIQUE SI VARIABLE       !
CARGU ! CV(NCP)      ! TR !  -> ! CHALEUR SPECIFIQUE SI VARIABLE       !
CARGU ! HBORD        ! TR !  -> ! COEFFICIENTS D'ECHANGE AUX BORDS     !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! TBORD        ! TR !  -> ! TEMPERATURES AUX BORDS               !
CARGU ! (NFABOR)     !    !     !                                      !
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
      INCLUDE "cstphy.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0, IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          ISVTB
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          NCP    , NCV    , IENTHA
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION CPCST  , CVCST
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*),PROPFA(NFAC,*),PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION CP(NCP), CV(NCV)
      DOUBLE PRECISION HBORD(NFABOR),TBORD(NFABOR)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          NBCCOU, INBCOU, INBCOO, NBFCOU, IFAC, ILOC, IEL
      INTEGER          ITFLUI, IHPARO
      INTEGER          IDEBIA, IDEBRA, IFINIA, IFINRA, IPFCOU, MODE
      INTEGER          ICCFTH, IMODIF, IPHAS
      INTEGER          IEPSEL, IEPSFA, IGAMAG, IXMASM
      DOUBLE PRECISION ENTHAL, TEMPER, ENERG, CVT
C
C***********************************************************************
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C     COUPLAGE SYRTHES : CALCUL DE LA TEMPERATURE FLUIDE ET DU
C                        COEFFICIENT D'ECHANGE
C=======================================================================
C
C     RECUPERATION DU NOMBRE DE CAS DE COUPLAGE
C
      CALL NBCSYR (NBCCOU)
C     ===========
C
C---> BOUCLE SUR LES CAS DE COUPLAGE DE TYPE SYRTHES
C
      DO INBCOU = 1, NBCCOU
C
C        NOMBRE DE FACES DE BORD PAR CAS DE COUPLAGE
        INBCOO = INBCOU
        CALL NBFSYR (INBCOO, NBFCOU)
C       ===========
C
C        GESTION MEMOIRE POUR CONSTRUIRE LES TABLEAUX
C
        IPFCOU = IDEBIA
        IFINIA = IPFCOU + NBFCOU
C
        ITFLUI = IDEBRA
        IHPARO = ITFLUI + NBFCOU
        IFINRA = IHPARO + NBFCOU
C
C Compressible : couplage avec l'energie
C
        IF(IENTHA .EQ. 2) THEN
          IEPSEL = IFINRA
          IEPSFA = IEPSEL + NCELET
          IGAMAG = IEPSFA + NFABOR
          IXMASM = IGAMAG + NCELET
          IFINRA = IXMASM + NCELET
        ENDIF
C
C Fin Compressible
C
        CALL RASIZE('COUPBO',IFINIA)
C       ===========
        CALL RASIZE('COUPBO',IFINRA)
C       ===========
C
C        BOUCLE SUR LES FACES DE COUPLAGE ET CALCUL DES COEFFICIENTS
C
        INBCOO = INBCOU
        CALL LFASYR(INBCOO, NBFCOU, IA(IPFCOU))
C       ===========
C
        DO ILOC = 1, NBFCOU
C
          IFAC = IA(IPFCOU+ILOC-1)
C
C           TEMPERATURES FLUIDES SAUVEGARDEES
          RA(ITFLUI+ILOC-1) = TBORD(IFAC)
C
C           COEFFICIENTS D'ECHANGE SAUVEGARDES
          RA(IHPARO+ILOC-1) = HBORD(IFAC)
C
        ENDDO
C
C        SI ENTHALPIE, ON TRANSFORME EN TEMPERATURE
C          Il est necessaire de transmettre a SYRTHES des Temperatures
C          Afin de conserver le flux Phi = (lambda/d     ) Delta T
C                                 ou Phi = (lambda/(d Cp)) Delta H
C            on multiplie HBORD = lambda/(d Cp) par Cp pris dans la
C              cellule adjacente.
C          Le resultat n'est pas garanti (conservation en particulier),
C             on ajoute donc un avertissement.
C
C
        IF(IENTHA.EQ.1) THEN
C
          WRITE(NFECRA,1000)
          MODE = 1
          DO ILOC = 1, NBFCOU
            IFAC = IA(IPFCOU+ILOC-1)
            IEL  = IFABOR(IFAC)
            ENTHAL = RA(ITFLUI+ILOC-1)
            CALL USTHHT (MODE   , ENTHAL , TEMPER  )
C           ===========
            RA(ITFLUI+ILOC-1) = TEMPER
            IF(NCP.EQ.NCELET) THEN
              RA(IHPARO+ILOC-1) = RA(IHPARO+ILOC-1)*CP(IEL)
            ELSE
              RA(IHPARO+ILOC-1) = RA(IHPARO+ILOC-1)*CPCST
            ENDIF
          ENDDO
C
        ELSE IF(IENTHA.EQ.2) THEN
C
C        SI ENERGIE, ON TRANSFORME EN TEMPERATURE
C          Il est necessaire de transmettre a SYRTHES des Temperatures
C          Afin de conserver le flux Phi = (lambda/d     ) Delta T
C                                 ou Phi = (lambda/(d Cv)) (Ei - Ep)
C            on multiplie HBORD = lambda/(d Cv) par Cv pris dans la
C              cellule adjacente.
C            noter que Ei = Cv Ti + 1/2 Ui*Ui + Epsilon_sup_i
C               et que Ep = Cv Tp + 1/2 Ui*Ui + Epsilon_sup_i
C               (l'ecart est donc bien Cv Delta T)
C
C     Modif temperature et coeff d'echange
C
C       Calcul de e - CvT
C
          ICCFTH = 7
          IMODIF = 0
          IPHAS  = IPHSCA(ISVTB)
C
          CALL USCFTH
C         ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IMODIF , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RA(IEPSEL) , RA(IEPSFA) , RA(IGAMAG) , RA(IXMASM) ,
C        ----------   ---------
     &   RDEVEL , RTUSER , RA     )
C
          DO ILOC = 1, NBFCOU
            IFAC  = IA(IPFCOU+ILOC-1)
            IEL   = IFABOR(IFAC)
            ENERG = RA(ITFLUI+ILOC-1)
            CVT   = ENERG
     &             -(0.5D0*( RTP(IEL,IU(IPHAS))**2
     &                      +RTP(IEL,IV(IPHAS))**2
     &                      +RTP(IEL,IW(IPHAS))**2)
     &               + RA(IEPSEL+IEL-1)           )
             IF(NCV.EQ.NCELET) THEN
               RA(ITFLUI+ILOC-1) = CVT/CV(IEL)
               RA(IHPARO+ILOC-1) = RA(IHPARO+ILOC-1)*CV(IEL)
             ELSE
               RA(ITFLUI+ILOC-1) = CVT/CVCST
               RA(IHPARO+ILOC-1) = RA(IHPARO+ILOC-1)*CVCST
             ENDIF
           ENDDO
C
        ENDIF
C
C       ENVOI DE LA TEMPERATURE FLUIDE ET DU COEFFICIENT D'ECHANGE
C
        INBCOO = INBCOU
        CALL VARSYO (INBCOO, NBFCOU, RA(ITFLUI), RA(IHPARO))
C       ===========
C
      ENDDO
C
C======================================================================
C     FIN DES COUPLAGES DE BORD
C======================================================================
C
      RETURN
C
C FORMATS
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION : COUPLAGE SYRTHES AVEC CALCUL EN ENTHALPIE   ',/,
     &'@    *********                                               ',/,
     &'@      OPTION NON VALIDEE - CONTACTER L''EQUIPE DE DVPT      ',/,
     &'@                                                            ')
C
      END
c@z
