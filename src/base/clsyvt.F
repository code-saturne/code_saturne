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
                        SUBROUTINE CLSYVT
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ISYMPA ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , RCODCL ,
     &   COEFU  , RIJIPB , COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CONDITIONS LIMITES EN SYMETRIE POUR LES VECTEURS ET TENSEURS
CFONC
CFONC ON SUPPOSE QUE ICODCL(IU) = 4 =>
CFONC                      SYMETRIE POUR LA VITESSE ET RIJ
CFONC   (A PRIORI PEU RESTRICTIF EN MONOPHASIQUE)
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
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE PHASE                      !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ICODCL       ! TE ! <-  ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU ! ISYMPA       ! TE ! <-> ! ZERO POUR ANNULER LE FLUX DE MASSE   !
CARGU ! (NFABOR     )!    !     !(SYMETRIES ET PAROIS AVEC CL COUPLEES)!
CARGU !              !    !     ! UN SINON                             !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET !    !     !                                      !
CARGU ! SURFAC       ! TR !  -> ! VECTEUR SURFACE DES FACES INTERNES   !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! SURFBO       ! TR !  -> ! VECTEUR SURFACE DES FACES DE BORD    !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! CDGFAC       ! TR !  -> ! CENTRE DE GRAVITE DES FACES INTERNES !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS (OPTIONNEL)    !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! RCODCL       ! TR ! <-  ! VALEUR DES CONDITIONS AUX LIMITES    !
CARGU !  (NFABOR,NVAR!    !     !  AUX FACES DE BORD                   !
CARGU !              !    !     ! RCODCL(1) = VALEUR DU DIRICHLET      !
CARGU !              !    !     ! RCODCL(2) = VALEUR DU COEF. D'ECHANGE!
CARGU !              !    !     !  EXT. (INFINIE SI PAS D'ECHANGE)     !
CARGU !              !    !     ! RCODCL(3) = VALEUR DE LA DENSITE DE  !
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2 OU      !
CARGU !              !    !     !  HAUTEUR DE RUGOSITE (m) si ICODCL=6 !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
CARGU ! COEFU        ! TR !  -> ! TAB DE TRAV POUR VALEURS EN IPRIME   !
CARGU ! (NFABOR,3   )!    !     !  DES COMP DE LA VITESSE AU BORD      !
CARGU ! RIJIPB       ! TR !  -> ! TAB DE TRAV POUR VALEURS EN IPRIME   !
CARGU ! (NFABOR,6   )!    !     !  DES RIJ AU BORD                     !
CARGU ! COEFA, COEFB ! TR ! <-> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
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
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "albase.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IPHAS
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR), ISYMPA(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION COEFU(NFABOR,NDIM), RIJIPB(NFABOR,6)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, II, ISOU
      INTEGER          IUIPH , IVIPH , IWIPH
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          ICLU  , ICLV  , ICLW
      INTEGER          ICL11 , ICL22 , ICL33 , ICL12 , ICL13 , ICL23
      INTEGER          ICLUF , ICLVF , ICLWF
      INTEGER          ICLVAR
      DOUBLE PRECISION RNX, RNY, RNZ, RXNN
      DOUBLE PRECISION UPX, UPY, UPZ, USN
      DOUBLE PRECISION TX, TY, TZ, TXN, T2X, T2Y, T2Z
      DOUBLE PRECISION CLSYME
      DOUBLE PRECISION ELOGLO(3,3), ALPHA(6,6)
      DOUBLE PRECISION SRFBNF, RCODCN
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C --- Variables
      IUIPH  = IU  (IPHAS)
      IVIPH  = IV  (IPHAS)
      IWIPH  = IW  (IPHAS)
      IF(ITYTUR(IPHAS).EQ.3) THEN
        IR11IP = IR11(IPHAS)
        IR22IP = IR22(IPHAS)
        IR33IP = IR33(IPHAS)
        IR12IP = IR12(IPHAS)
        IR13IP = IR13(IPHAS)
        IR23IP = IR23(IPHAS)
      ENDIF
C
C --- Conditions aux limites
      ICLU   = ICLRTP(IUIPH ,ICOEF)
      ICLV   = ICLRTP(IVIPH ,ICOEF)
      ICLW   = ICLRTP(IWIPH ,ICOEF)
      IF(ITYTUR(IPHAS).EQ.3) THEN
        ICL11  = ICLRTP(IR11IP,ICOEF)
        ICL22  = ICLRTP(IR22IP,ICOEF)
        ICL33  = ICLRTP(IR33IP,ICOEF)
        ICL12  = ICLRTP(IR12IP,ICOEF)
        ICL13  = ICLRTP(IR13IP,ICOEF)
        ICL23  = ICLRTP(IR23IP,ICOEF)
      ENDIF
C
      ICLUF  = ICLRTP(IUIPH ,ICOEFF)
      ICLVF  = ICLRTP(IVIPH ,ICOEFF)
      ICLWF  = ICLRTP(IWIPH ,ICOEFF)
C
C
C --- Boucle sur les faces de bord : debut
      DO IFAC = 1, NFABOR
C
C --- Test sur la presence d'une condition de symetrie vitesse : debut
        IF( ICODCL(IFAC,IUIPH).EQ.4 ) THEN
C
C --- Pour annuler le flux de masse
          ISYMPA(IFAC) = 0
C
C --- Grandeurs geometriques
          SRFBNF = RA(ISRFBN-1+IFAC)
C
C=======================================================================
C 1. REPERE LOCAL
C      POUR LA VITESSE, SEULE EST NECESSAIRE LA NORMALE
C      POUR RIJ, IL FAUT LE REPERE COMPLET
C=======================================================================
C
C ---> NORMALE UNITAIRE
C
          RNX = SURFBO(1,IFAC)/SRFBNF
          RNY = SURFBO(2,IFAC)/SRFBNF
          RNZ = SURFBO(3,IFAC)/SRFBNF
C
C     En ALE, on a eventuellement une vitesse de deplacement de la face
C       donc seule la composante normale importe (on continue a determiner
C       TX a partir de la vitesse tangentielle absolue car l'orientation
C       de TX et T2X est sans importance pour les symetries)
          RCODCN = 0.D0
          IF (IALE.EQ.1) THEN
            RCODCN = RCODCL(IFAC,IUIPH,1)*RNX
     &             + RCODCL(IFAC,IVIPH,1)*RNY
     &             + RCODCL(IFAC,IWIPH,1)*RNZ
          ENDIF
C
          UPX = COEFU(IFAC,1)
          UPY = COEFU(IFAC,2)
          UPZ = COEFU(IFAC,3)
C
          IF (ITYTUR(IPHAS).EQ.3) THEN
C
C ---> VITESSE TANGENTIELLE RELATIVE
C
            USN = UPX*RNX+UPY*RNY+UPZ*RNZ
            TX  = UPX -USN*RNX
            TY  = UPY -USN*RNY
            TZ  = UPZ -USN*RNZ
            TXN = SQRT( TX**2 +TY**2 +TZ**2 )
C
C ---> TANGENTE UNITAIRE
C
            IF( TXN.GE.EPZERO) THEN
C
              TX  = TX/TXN
              TY  = TY/TXN
              TZ  = TZ/TXN
C
            ELSE
C
C      SI LA VITESSE EST NULLE, LE VECTEUR T EST NORMAL ET QCQUE
C
              IF(ABS(RNY).GE.EPZERO.OR.ABS(RNZ).GE.EPZERO)THEN
                RXNN = SQRT(RNY**2+RNZ**2)
                TX  =  0.D0
                TY  =  RNZ/RXNN
                TZ  = -RNY/RXNN
              ELSEIF(ABS(RNX).GE.EPZERO.OR.ABS(RNZ).GE.EPZERO)THEN
                RXNN = SQRT(RNX**2+RNZ**2)
                TX  =  RNZ/RXNN
                TY  =  0.D0
                TZ  = -RNX/RXNN
              ELSE
                WRITE(NFECRA,1000)IFAC,RNX,RNY,RNZ
                CALL CSEXIT (1)
              ENDIF
C
            ENDIF
C
C
C ---> T2 = RN X T (OU X EST LE PRODUIT VECTORIEL)
C
            T2X = RNY*TZ - RNZ*TY
            T2Y = RNZ*TX - RNX*TZ
            T2Z = RNX*TY - RNY*TX
C
C     --> MATRICE ORTHOGONALE DE CHANGEMENT DE BASE ELOGLOij
C         (DE LA BASE LOCALE VERS LA BASE GLOBALE)
C
C                            |TX  -RNX  T2X|
C                   ELOGLO = |TY  -RNY  T2Y|
C                            |TZ  -RNZ  T2Z|
C
C         SA TRANSPOSEE ELOGLOt EST SON INVERSE
C
C
            ELOGLO(1,1) =  TX
            ELOGLO(1,2) = -RNX
            ELOGLO(1,3) =  T2X
            ELOGLO(2,1) =  TY
            ELOGLO(2,2) = -RNY
            ELOGLO(2,3) =  T2Y
            ELOGLO(3,1) =  TZ
            ELOGLO(3,2) = -RNZ
            ELOGLO(3,3) =  T2Z
C
C     --> ON CALCULE ALPHA(6,6)
C
C       SOIT f LE CENTRE DE LA FACE DE BORD ET
C            I LE CENTRE DE LA CELLULE CORRESPONDANTE
C
C       EN NOTE RG (RESP RL) INDICE PAR f OU PAR I
C          LE TENSEUR DE REYNOLDS DANS LA BASE GLOBALE (RESP LOCALE)
C
C       LA MATRICE ALPHA APPLIQUEE AU VECTEUR GLOBAL EN I'
C         (RG11,I'|RG22,I'|RG33,I'|RG12,I'|RG13,I'|RG23,I')t
C         DOIT DONNER LES VALEURS A IMPOSER A LA FACE
C         (RG11,f |RG22,f |RG33,f |RG12,f |RG13,f |RG23,f )t
C         AUX CONDITIONS LIMITES DE DIRICHLET PRES (AJOUTEES ENSUITE)
C
C       ON LA DEFINIT EN CALCULANT RG,f EN FONCTION DE RG,I' COMME SUIT
C
C         RG,f = ELOGLO.RL,f.ELOGLOt (PRODUITS MATRICIELS)
C
C                          | RL,I'(1,1)     B*U*.Uk     C*RL,I'(1,3) |
C           AVEC    RL,f = | B*U*.Uk       RL,I'(2,2)       0        |
C                          | C*RL,I'(1,3)     0         RL,I'(3,3)   |
C
C                  AVEC    RL,I = ELOGLOt.RG,I'.ELOGLO
C                          B = 0
C                    ET    C = 0 EN PAROI (1 EN SYMETRIE)
C
C
C
C          ON CALCULE EN FAIT   ELOGLO.PROJECTEUR.ELOGLOt
C
C
            CLSYME=1.D0
            CALL CLCA66 ( CLSYME , ELOGLO , ALPHA )
C           ===========

          ENDIF
C
C=======================================================================
C 2. CONDITIONS SUR LES VITESSES (PARTIELLEMENT IMPLICITES)
C=======================================================================
C
          COEFA(IFAC,ICLU) = RCODCN*RNX - RNX*(RNY*UPY+RNZ*UPZ)
          COEFB(IFAC,ICLU) = 1.D0-RNX**2
          COEFA(IFAC,ICLV) = RCODCN*RNY - RNY*(RNZ*UPZ+RNX*UPX)
          COEFB(IFAC,ICLV) = 1.D0-RNY**2
          COEFA(IFAC,ICLW) = RCODCN*RNZ - RNZ*(RNX*UPX+RNY*UPY)
          COEFB(IFAC,ICLW) = 1.D0-RNZ**2
C
C=======================================================================
C 3. CONDITIONS SUR RIJ (PARTIELLEMENT IMPLICITES)
C=======================================================================
C
          IF (ITYTUR(IPHAS).EQ.3) THEN
C
            DO ISOU = 1, 6
C
              IF(ISOU.EQ.1) ICLVAR = ICL11
              IF(ISOU.EQ.2) ICLVAR = ICL22
              IF(ISOU.EQ.3) ICLVAR = ICL33
              IF(ISOU.EQ.4) ICLVAR = ICL12
              IF(ISOU.EQ.5) ICLVAR = ICL13
              IF(ISOU.EQ.6) ICLVAR = ICL23
C
              COEFA(IFAC,ICLVAR) = 0.0D0
              COEFB(IFAC,ICLVAR) = 0.0D0
C
            ENDDO
C
            DO ISOU = 1,6
C
              IF(ISOU.EQ.1) THEN
                ICLVAR = ICL11
              ELSE IF(ISOU.EQ.2) THEN
                ICLVAR = ICL22
              ELSE IF(ISOU.EQ.3) THEN
                ICLVAR = ICL33
              ELSE IF(ISOU.EQ.4) THEN
                ICLVAR = ICL12
              ELSE IF(ISOU.EQ.5) THEN
                ICLVAR = ICL13
              ELSE IF(ISOU.EQ.6) THEN
                ICLVAR = ICL23
              ENDIF
C
C     IMPLICITATION PARTIELLE EVENTUELLE DES CL
              IF (ICLSYR(IPHAS).EQ.1) THEN
                DO II = 1, 6
                  IF (II.NE.ISOU) THEN
                    COEFA(IFAC,ICLVAR) = COEFA(IFAC,ICLVAR) +
     &                   ALPHA(ISOU,II) * RIJIPB(IFAC,II)
                  ENDIF
                ENDDO
                COEFB(IFAC,ICLVAR) = ALPHA(ISOU,ISOU)
              ELSE
                DO II = 1, 6
                  COEFA(IFAC,ICLVAR) = COEFA(IFAC,ICLVAR) +
     &                 ALPHA(ISOU,II) * RIJIPB(IFAC,II)
                ENDDO
                COEFB(IFAC,ICLVAR) = 0.D0
              ENDIF
C
            ENDDO
C
          ENDIF
C
        ENDIF
C --- Test sur la presence d'une condition de symetrie vitesse : fin
C
      ENDDO
C ---  Boucle sur les faces de bord : fin
C
C=======================================================================
C 4. COEFAF et COEFBF BIDONS POUR LES VITESSES
C=======================================================================
C
        IF(ICLU.NE.ICLUF) THEN
          DO IFAC = 1, NFABOR
            IF( ICODCL(IFAC,IUIPH).EQ.4) THEN
              COEFA(IFAC,ICLUF) = COEFA(IFAC,ICLU)
              COEFB(IFAC,ICLUF) = COEFB(IFAC,ICLU)
              COEFA(IFAC,ICLVF) = COEFA(IFAC,ICLV)
              COEFB(IFAC,ICLVF) = COEFB(IFAC,ICLV)
              COEFA(IFAC,ICLWF) = COEFA(IFAC,ICLW)
              COEFB(IFAC,ICLWF) = COEFB(IFAC,ICLW)
            ENDIF
          ENDDO
        ENDIF
C
C=======================================================================
C 7.  FORMATS
C=======================================================================
C
C
 1000 FORMAT(/,' LA NORMALE A LA FACE DE BORD DE SYMETRIE ',I10,/,
     &         ' EST NULLE ; COORDONNEES : ',3E12.5)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
