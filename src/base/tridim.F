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
                        SUBROUTINE TRIDIM
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , ITRALE ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   TSLAGR , COEFA  , COEFB  , FRCXT  ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  RESOLUTION DES EQUATIONS N-S MONOPHASIQUES INCOMPRESSIBLE ET DES
CFONC  EQUATIONS SCALAIRES POUR UN PAS DE TEMPS
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
CARGU ! ITRALE       ! E  !  -> ! NUMERO D'ITERATION POUR L'ALE        !
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
CARGU !   (NFAC+1)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (LNDFAC)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFABOR+1) !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR)   !    !     !  (OPTIONNEL)                         !
CARGU ! ISOSTD       ! TE !  -> ! INDICATEUR DE SORTIE STANDARD        !
CARGU !    (NFABOR+1)!    !     !  +NUMERO DE LA FACE DE REFERENCE     !
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
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! TSLAGR       ! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !(NCELET,*)    !    !     !     LAGRANGIEN                       !
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
C
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
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "albase.h"
      INCLUDE "alstru.h"
      INCLUDE "radiat.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "coincl.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagdim.h"
      INCLUDE "lagran.h"
      INCLUDE "vortex.h"
      INCLUDE "ihmpre.h"
      INCLUDE "matiss.h"
C
C les includes pp* ne servent que pour recuperer le pointeur IIZFPP
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0 , ITRALE
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ISOSTD(NFABOR+1,NPHAS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION TSLAGR(NCELET,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IDBIA1, IDBRA1
      INTEGER          IFINIA, IFINRA, IFNIA1, IFNRA1
      INTEGER          IEL   , IFAC  , IVAR  , ISCAL , IAPPEL
      INTEGER          NCP   , NCV   , IOK
      INTEGER          IICODC, IRCODC
      INTEGER          ICOEFU, IRIJIP, IHBORD, ITBORD
      INTEGER          IDTR  , IVISCF, IVISCB, IVISFI, IVISBI, IIPTOT
      INTEGER          IDAM  , IXAM
      INTEGER          ICOFBD
      INTEGER          IDRTP , IGRDP , ISMBR , IROVSD
      INTEGER          ITINSK, ITINSE, IDIVU , IPRV2F
      INTEGER          IW1   , IW2   , IW3   , IW4   , IW5   , IW6
      INTEGER          IW7   , IW8   , IW9   , IW10  , IW11  , IW12
      INTEGER          IXMIJ
      INTEGER          IFRCHY, IDFRHY, IDFRCX
      INTEGER          IGRDVT, IPRODU, IGRARX, IGRARY, IGRARZ
      INTEGER          IESFLM, IESFLB
      INTEGER          ITRAVA, IXIMPA, IUVWK
      INTEGER          NBCCOU
      INTEGER          ICOFBR
C
      INTEGER          ISVHB , ISVTB
      INTEGER          IPHAS , KPHAS , II    , IPPCP , IENTHA, IPPCV
      INTEGER          IKIPH , IEIPH , IOMIPH
      INTEGER          IUIPH , IVIPH , IWIPH , IPRIPH, IPHIPH, IPHASS
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          IPCROM, IPCROA
      INTEGER          IPRNEW, IDIMTE, ITENSO
      INTEGER          IFINIB, IFINRB, IIIFAP
      INTEGER          IFLUA , IFLUB
      INTEGER          ICOAX , ICOBX , ICOAY , ICOBY , ICOAZ , ICOBZ
      INTEGER          IQFX  , IQFY  , IQFZ  , ICOEFQ
      INTEGER          IIRHO , IIRHOB, ICOEFX
      INTEGER          IRTDP , IDRTDP, ICOFAY, ICOFBY, IISMPH
      INTEGER          IISOTH, ITEXT , ITINT , ITEK
      INTEGER          ICORUA, ICORUB, IFLXMA, IFLXMB
      INTEGER          ITERNS, INSLST, ICVRGE, IUETBO, IVSVDR
      INTEGER          IWFLMS, IWFLMB
      INTEGER          IWCF  , IPH   , IFLMAS, IFLMAB
      INTEGER          ITALIM, ITRFIN, INEEFL
      INTEGER          IFLALF, IFLALB, IPRALE, ICOALE
      INTEGER          MAXELT, ILS, IILZFB, NBZFMX, NOZFMX, IQCALC

      DOUBLE PRECISION CPCST , TDITOT, TDIST2, TDIST1, CVCST
      DOUBLE PRECISION RO0IPH, P0IPH, PR0IPH, XXP0, XYP0, XZP0
      DOUBLE PRECISION RELAXK, RELAXE, RELAXW
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
      INTEGER          INFPAR
      SAVE             INFPAR
C
C***********************************************************************
C
C
C=======================================================================
C 1.  INITIALISATION
C=======================================================================
C
      IF(IWARNI(IU(1)).GE.1) THEN
        WRITE(NFECRA,1000)
      ENDIF
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      MAXELT = MAX(NCELET, NFAC, NFABOR)
C
      IPASS = IPASS + 1
C
C --- Indicateur de stockage d'un scalaire et de son coef
C      d'echange associe.
C     Pour le moment, on stocke uniquement dans le cas couplage SYRTHES.
C     ISVTB donne le numero du scalaire (on suppose qu'il n'y en a
C       qu'un).
C     Dans le cas ou on a un couplage avec le module thermique 1D en paroi,
C     on utilise le meme scalaire que celui qui sert a Syrthes (s'il y a
C     couplage Syrthes), sinon on stocke le scalaire thermique de la phase 1.
C
      CALL NBCSYR (NBCCOU)
C     ===========
      ISVHB = 0
      ISVTB = 0
      IF (NBCCOU .GE. 1) THEN
        DO ISCAL = 1, NSCAL
          IF(ICPSYR(ISCAL).EQ.1) THEN
            ISVHB = ISCAL
            ISVTB = ISCAL
          ENDIF
        ENDDO
      ENDIF
C
      IF ((NFPT1T.GT.0).AND.(NBCCOU.LE.0)) THEN
        IPHAS = 1
        ISVHB = ISCALT(IPHAS)
        ISVTB = ISCALT(IPHAS)
      ENDIF
C
C     Si la distance a la paroi doit etre mise a jour, on l'initialise a GRAND
C     des maintenant (pour le premier passage dans phyvar en k-omega)
      IF(IPASS.EQ.1.AND.INEEDY.EQ.1.AND.ABS(ICDPAR).EQ.1.AND.
     &                                  IMAJDY.EQ.0) THEN
        DO IEL = 1, NCEL
          RA(IDIPAR+IEL-1) = GRAND
        ENDDO
      ENDIF
C
C=======================================================================
C 2.  AU DEBUT DU CALCUL ON REINITIALISE LA PRESSION
C=======================================================================
C
C On le fait sur 2 pas de temps, car souvent, le champ de flux de masse
C   initial n'est pas a divergence nulle (CL incluses) et l'obtention
C   d'un flux a divergence nulle coherent avec la contrainte stationnaire
C   peut prendre quelques pas de temps.
C Remarquer que la pression est rattrapee dans l'etape de stokes.
C On ne le fait pas dans le cas de la prise en compte de la pression
C   hydrostatique, ni dans le cas du compressible
C
      IF( NTCABS.LE.2 .AND. ISUITE.EQ.0 .AND. IPHYDR.EQ.0
     &                .AND. IPPMOD(ICOMPF).LT.0           ) THEN
C
        IF(IWARNI(IPR(1)).GE.2) THEN
          WRITE(NFECRA,2000) NTCABS
        ENDIF
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
            IIPTOT = IPPROC(IPRTOT(IPHAS))
            RO0IPH = RO0  (IPHAS)
            P0IPH  = P0   (IPHAS)
            PR0IPH = PRED0(IPHAS)
            XXP0   = XYZP0(1,IPHAS)
            XYP0   = XYZP0(2,IPHAS)
            XZP0   = XYZP0(3,IPHAS)
            DO IEL = 1, NCEL
              RTP(IEL,IPR(IPHAS)) = PR0IPH
              PROPCE(IEL,IIPTOT) = P0IPH
     &             + RO0IPH*( GX*(XYZCEN(1,IEL)-XXP0)
     &                      + GY*(XYZCEN(2,IEL)-XYP0)
     &                      + GZ*(XYZCEN(3,IEL)-XZP0) )
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
 2000 FORMAT(
     &  ' REINITIALISATION DE LA PRESSION A L''ITERATION ',I10)
C
C=======================================================================
C 3.  COMMUNICATIONS
C=======================================================================
C
C ---> On echange ici les variables RTP en debut de pas de temps.
C         Ce n'est pas fait dans caltri pour faciliter la lecture
C         (manipulation des tableaux)
C      Ceci a l'avantage d'echanger egalement ce qui provient de
C         inivar et lecamo (et par rapport a une solution d'echange
C         en fin de pas de temps, on evite une communication)
C      L'inconvenient est que l'on ne dispose pas des halos sitot les
C         RTP calcules (fin de navsto, turbke, turrij, boucle dans scalai)
C      On pourrait penser a echanger aussi les PROPCE et le DT
C         Pour le moment ca ne s'impose pas : on a besoin d'avoir
C         echange RHO dans navsto pour un affichage et les viscosites
C         dans vissec. On fera les transferts localement quand necessaire.
C      Par principe, on suppose que
C         c'est a la charge de celui qui utilise des valeurs voisines de
C           s'assurer qu'elles existent (ss pgm utilisateurs en
C           particulier)
C         seules les RTPA sont echangees de maniere systematique
C           (eviter donc d'utiliser RTP)
C         le calcul du gradient fournit aussi les valeurs sur le halo
C           (utile pour les reconstructions)
C         les seules boucles sur NCELET sont des boucles d'initialisation
C           (on n'a pas a faire de calcul sur les cellules halo, elles
C            elles sont remplies par les routines de communication et
C            on y preleve seulement de l'information)
C
C
C
C ---> Parallelisme
C
      IF(IRANGP.GE.0) THEN
C
        DO IVAR = 1, NVAR
          CALL PARCVE (RTP(1,IVAR))
C         ===========
        ENDDO
C
      ENDIF
C
C ---> Periodicite
C
      IF(IPERIO.EQ.1) THEN
C
C  -- Vitesse
C
        DO IPHAS = 1, NPHAS
C
          IUIPH = IU(IPHAS)
          IVIPH = IV(IPHAS)
          IWIPH = IW(IPHAS)
          IDIMTE = 1
          ITENSO = 0
          CALL PERCVE
C         ===========
     &    ( IDIMTE , ITENSO ,
     &      RTP(1,IUIPH), RTP(1,IUIPH), RTP(1,IUIPH),
     &      RTP(1,IVIPH), RTP(1,IVIPH), RTP(1,IVIPH),
     &      RTP(1,IWIPH), RTP(1,IWIPH), RTP(1,IWIPH))
C
        ENDDO
C
C  -- Tenseur de Reynolds
C
        DO IPHAS = 1, NPHAS
C
          IF(ITYTUR(IPHAS).EQ.3) THEN
            IDIMTE = 2
            ITENSO = 0
            IR11IP = IR11(IPHAS)
            IR22IP = IR22(IPHAS)
            IR33IP = IR33(IPHAS)
            IR12IP = IR12(IPHAS)
            IR13IP = IR13(IPHAS)
            IR23IP = IR23(IPHAS)
            CALL PERCVE
C           ===========
     &    ( IDIMTE , ITENSO ,
     &      RTP(1,IR11IP), RTP(1,IR12IP), RTP(1,IR13IP),
     &      RTP(1,IR12IP), RTP(1,IR22IP), RTP(1,IR23IP),
     &      RTP(1,IR13IP), RTP(1,IR23IP), RTP(1,IR33IP) )
          ENDIF
C
C  -- Remarque pour le v2f
C     v2 (donc phi) est lie a une orientation locale, on peut donc le traiter
C     comme un scalaire dans la periodicite de rotation
C
C
        ENDDO
C
C  -- Variables scalaires
C
        DO IVAR = 1, NVAR
          DO IPHAS = 1, NPHAS
            IF(IVAR.NE.IU(IPHAS).AND.IVAR.NE.IV(IPHAS).AND.
     &                               IVAR.NE.IW(IPHAS).AND.
     &         (ITYTUR(IPHAS).NE.3.OR.
     &          (IVAR.NE.IR11(IPHAS).AND.IVAR.NE.IR22(IPHAS).AND.
     &           IVAR.NE.IR33(IPHAS).AND.IVAR.NE.IR12(IPHAS).AND.
     &           IVAR.NE.IR13(IPHAS).AND.IVAR.NE.IR23(IPHAS)))) THEN
C
              IDIMTE = 0
              ITENSO = 0
              CALL PERCVE
C             ===========
     &      ( IDIMTE , ITENSO ,
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR) )
C
            ENDIF
          ENDDO
        ENDDO
C
      ENDIF
C
C=======================================================================
C 4.  POUR IPHYDR ON DOIT COMMUNIQUER FRCXT AU PREMIER PASSAGE
C     (FRCXT SERT DANS TYPECL)
C     SI ICALHY=1, ON COMMUNIQUE AUSSI RHO POUR REMPLIR
C     PROPCE(1,IPPROC(IROMA(IPHAS)))
C=======================================================================
C
      IF(IPASS.EQ.1) THEN
C
C --- Communication de FRCXT
        IF (IPHYDR.EQ.1) THEN
C
          DO IPHAS = 1, NPHAS
            IF(IRANGP.GE.0) THEN
              CALL PARCVE (FRCXT(1,1,IPHAS))
C             ===========
              CALL PARCVE (FRCXT(1,2,IPHAS))
C             ===========
              CALL PARCVE (FRCXT(1,3,IPHAS))
C             ===========
            ENDIF
            IF(IPERIO.EQ.1) THEN
              IDIMTE = 1
              ITENSO = 0
              CALL PERCVE
C             ===========
     &    ( IDIMTE , ITENSO ,
     &      FRCXT(1,1,IPHAS),FRCXT(1,1,IPHAS),FRCXT(1,1,IPHAS),
     &      FRCXT(1,2,IPHAS),FRCXT(1,2,IPHAS),FRCXT(1,2,IPHAS),
     &      FRCXT(1,3,IPHAS),FRCXT(1,3,IPHAS),FRCXT(1,3,IPHAS) )
            ENDIF
C
          ENDDO
C
        ENDIF
C
C --- Communication de RHO
        IF (ICALHY.EQ.1) THEN
C
          DO IPHAS = 1, NPHAS
            IPCROM = IPPROC(IROM  (IPHAS))
            IF(IRANGP.GE.0) THEN
              CALL PARCVE (PROPCE(1,IPCROM))
C             ===========
            ENDIF
            IF(IPERIO.EQ.1) THEN
              IDIMTE = 0
              ITENSO = 0
              CALL PERCVE
C             ===========
     &    ( IDIMTE , ITENSO ,
     &      PROPCE(1,IPCROM),PROPCE(1,IPCROM),PROPCE(1,IPCROM),
     &      PROPCE(1,IPCROM),PROPCE(1,IPCROM),PROPCE(1,IPCROM),
     &      PROPCE(1,IPCROM),PROPCE(1,IPCROM),PROPCE(1,IPCROM) )
            ENDIF
          ENDDO
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C 5.  LES VALEURS COURANTES ECRASENT LES VALEURS ANTERIEURES
C=======================================================================
C
C --- Noter que exceptionnellement, on fait un calcul avec NCELET,
C       pour eviter une nouvelle communication sur RTPA et les autres
C       tableaux du pas de temps precedent
C
      DO IVAR = 1, NVAR
        DO IEL = 1, NCELET
          RTPA (IEL,IVAR) = RTP (IEL,IVAR)
        ENDDO
      ENDDO
C
      IF (ICALHY.EQ.1) THEN
        DO IPHAS = 1, NPHAS
          IPCROM = IPPROC(IROM  (IPHAS))
          IPCROA = IPPROC(IROMA (IPHAS))
          DO IEL = 1, NCELET
            PROPCE(IEL,IPCROA) = PROPCE(IEL,IPCROM)
          ENDDO
        ENDDO
      ENDIF
C
C=======================================================================
C 6. DANS LE CAS  "zero pas de temps" EN "SUITE" DE CALCUL
C      ON SORT ICI
C=======================================================================
C  on sort avant SCHTMP car sinon a l'ordre 2 en temps la valeur du
C  flux de masse au pas de temps precedent est ecrasee par la valeur
C  au pas de temps actuel et la valeur au pas de temps actuel est
C  remplacee par une extrapolation qui n'a pas lieu d'etre puisque
C  NTCABS n'est pas incremente. Dans le cas INPDT0=1 sans suite, il
C  n'y a pas de probleme puisque tous les flux de masse sont a 0 !
C  Si ITRALE=0, on est a l'iteration d'initialisation de l'ALE,
C  on ne touche pas au flux de masse non plus
C
C
      IF(INPDT0.EQ.1.AND.ISUITE.EQ.1) GOTO 200
C
      IF (ITRALE.GT.0) THEN
        IAPPEL = 1
        CALL SCHTMP
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER ,
     &   RA     )
      ENDIF
C
C
C=======================================================================
C 7.  CALCUL DES PROPRIETES PHYSIQUES VARIABLES
C      SOIT VARIABLES AU COURS DU TEMPS
C      SOIT VARIABLES LORS D'UNE REPRISE DE CALCUL
C        (VISCOSITES ET MASSE VOLUMIQUE)
C=======================================================================
C
      IF(IWARNI(IU(1)).GE.1) THEN
        WRITE(NFECRA,1010)
      ENDIF
C
      CALL MEMPHY
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NVAR   , NCELET , NCEL   , NFAC   , NFABOR , NPHAS  ,
     &   IW1    , IW2    , IW3    , IW4    ,
     &   IW5    , IW6    , IW7    , IW8    ,
     &   IW9    , IW10   , IXMIJ  ,
     &   IFINIA , IFINRA )
C
      CALL PHYVAR
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RA(IW1), RA(IW2), RA(IW3), RA(IW4),
     &   RA(IW5), RA(IW6), RA(IW7), RA(IW8), RA(IW9), RA(IW10),
     &   RA(IXMIJ) ,
     &   RDEVEL , RTUSER , RA     )
C
      IF (ITRALE.GT.0) THEN
        IAPPEL = 2
        CALL SCHTMP
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER ,
     &   RA     )
      ENDIF
C
C
C REMPLISSAGE DES COEFS DE PDC
C    ON Y PASSE MEME S'IL N'Y A PAS DE PDC SUR LE PROC COURANT AU CAS OU
C    UN UTILISATEUR DECIDERAIT D'AVOIR UN COEFF DE PDC DEPENDANT DE
C    LA VITESSE MOYENNE OU MAX.
C
C
      DO IPHAS = 1, NPHAS
C
        IF (NCPDCT(IPHAS).GT.0) THEN
C
          IAPPEL = 3
          IF (IMATIS.EQ.1) THEN
C
            CALL MTKPDC
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCEPDC(IPHAS) , IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IICEPD(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)) ,
     &   RDEVEL , RTUSER , RA     )
C
          ELSE
C
            ILS    = IDEBIA
            IDBIA1 = ILS + MAXELT
            CALL IASIZE('TRIDIM',IDBIA1)
C
            CALL USKPDC
C           ===========
     & ( IDBIA1 , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCEPDC(IPHAS) , IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IICEPD(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)) ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
C
C REMPLISSAGE DES COEFS DE TERME SOURCE DE MASSE
C
      DO IPHAS = 1, NPHAS
C
C    ON Y PASSE MEME S'IL N'Y A PAS DE TSM SUR LE PROC COURANT AU CAS OU
C    UN UTILISATEUR DECIDERAIT D'AVOIR UN TSM DEPENDANT DE
C    VALEURS GLOBALES OU MAX.
        IF(NCTSMT(IPHAS).GT.0) THEN
C
          ILS    = IDEBIA
          IDBIA1 = ILS + MAXELT
          CALL IASIZE('TRIDIM',IDBIA1)
C
C     Mise a zero du tableau de type de TS masse et source
          DO II = 1, NCETSM(IPHAS)*NVAR
            IA(IITPSM(IPHAS)+II-1) = 0
            RA(ISMACE(IPHAS)+II-1) = 0.D0
          ENDDO
C
          IAPPEL = 3
          CALL  USTSMA
C         ============
     & ( IDBIA1 , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDC(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCETSM(IPHAS)   , IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) ,
     &   IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)), RA(ISMACE(IPHAS)),
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 8.  CALCUL DU NOMBRE DE COURANT ET DE FOURIER
C     CALCUL DU PAS DE TEMPS SI VARIABLE
C=======================================================================
C
      IF(IWARNI(IU(1)).GE.1) THEN
        WRITE(NFECRA,1020)
      ENDIF
C
      CALL MEMDTV
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IVISCF , IVISCB , IDAM   , ICOFBD , IW1    , IW2    , IW3    ,
     &   ICOFBR , IGRARX , IGRARY , IGRARZ , IWCF   ,
     &   IPTLRO , IPPMOD(ICOMPF) ,
     &   IFINIA , IFINRA )
C
      IPHAS = 1
C
      CALL DTTVAR
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS)   , NCETSM(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IWARNI(IU(1))   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)), IA(IICESM(IPHAS)), IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)) , RA(ISMACE(IPHAS)),
     &   RA(IVISCF)      , RA(IVISCB)      , RA(IDAM)      ,
     &   RA(ICOFBD)      , RA(IW1)         , RA(IW2)       , RA(IW3)   ,
     &   RA(ICOFBR)      , RA(IGRARX)      , RA(IGRARY)    , RA(IGRARZ),
     &   RA(IWCF),
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 9.  CHARGEMENT ET TRADUCTION DES CONDITIONS AUX LIMITES
C=======================================================================
C
      IF(IWARNI(IU(1)).GE.1) THEN
        WRITE(NFECRA,1030)
      ENDIF
C
C  -- Methode des vortex en LES :
C     Definition ou modification eventuelle des parametres
C     Mise a jour des vortex
C
      IF (IVRTEX.EQ.1) THEN
C
        ILS    = IFINIA
        IFNIA1 = ILS + MAXELT
        CALL IASIZE('TRIDIM',IFNIA1)
C
        IPHAS  = 1
        IAPPEL = 2
        CALL USVORT
C       ===========
     & ( IFNIA1 , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML, MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IIREPV)      ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   DT     , RTPA   ,
     &   PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
C
C     Verification des donnees entrees par l'utilisateur
C       (au premier passage seulement)
        IF (IPASS.EQ.1) THEN
          CALL VORVER ( NFABOR , IA(IIREPV)  , IAPPEL )
C         ===========
        ENDIF
C
        IF(IRANGP.LE.0) THEN
          CALL VORTEX
C         ===========
     & ( IA(IIVRCE) , RA(IVISV)  , RA(IXYZV) ,
     &   RA(IYZCEL) , RA(IUVORT) , RA(IVVORT),
     &   RA(IWVORT) ,
     &   RA(IYZVOR) , RA(IYZVOA) , RA(ISIGNV),
     &   RA(IXSIGM) , RA(IXGAMM) , RA(IXTMP) ,
     &   RA(IXTMPL) )
        ENDIF
C
C -- Fin de zone Methode des vortex
C
      ENDIF
C
C --- Methode ALE : debut de boucle d'implicitation du deplacement des
C       structures. ITRFIN=0 indique qu'on a besoin de refaire une iteration
C       pour Syrthes, T1D ou rayonnement.
      ITALIM = 1
      ITRFIN = 1
      IDBIA1 = IDEBIA
      IDBRA1 = IDEBRA
      IFLALF = 1
      IFLALB = 1
      IPRALE = 1
      INEEFL = 0
      IF (IALE.EQ.1 .AND. NALIMX.GT.1 .AND. ITRALE.GT.NALINF) THEN
C     On reserve certains tableaux pour permettre le retour a l'etat
C       initial en fin d'iteration ALE
C       - flux de masse
C       - conditions aux limites de gradient de P et U (car on a un appel
C         a GDRCEL pour les non orthogonalites pour calculer les CL reelles)
C         -> n'est peut-etre pas reellement necessaire
C       - la pression initiale (car RTPA est aussi ecrase dans le cas
C         ou NTERUP>1) -> on pourrait optimiser en ne reservant que si
C         necessaire ...
C       Pas la peine de tester les depassements car on passe dans
C       memcli juste apres.
        IFLALF = IDEBRA
        IFLALB = IFLALF + NFAC
        ICOALE = IFLALB + NFABOR
        IPRALE = ICOALE + 8*NFABOR
        IDBRA1 = IPRALE + NCELET
        INEEFL = 1
C
        IF (NBCCOU.GT.0 .OR. NFPT1T.GT.0 .OR. IIRAYO.GT.0) ITRFIN = 0
C
      ENDIF
C
 300  CONTINUE
C
C
      CALL MEMCLI
C     ===========
     & ( IDBIA1 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ISVHB  , ISVTB  ,
     &   IICODC , IRCODC ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    ,
     &   ICOEFU , IRIJIP , IUETBO , IVSVDR , IHBORD , ITBORD ,
     &   IFINIA , IFINRA )
C
      CALL PRECLI
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICODC)      , IA(IIZFPP)      ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RA(IRCODC) , RA(ICOEFU) ,
     &   RA(IW1   ) , RA(IW2   ) , RA(IW3   ) ,
     &   RA(IW4   ) , RA(IW5   ) , RA(IW6   ) ,
     &   RDEVEL , RTUSER , RA     )
C
C
      IF (IMATIS.EQ.0) THEN
C
C ON NE FAIT PAS DE MATISSE
C
C     - Interface Code_Saturne
C       ======================
C
        IF(IIHMPR.EQ.1) THEN
C
C N.B. Zones de face de bord : on utilise provisoirement les zones des
C    physiques particulieres, meme sans physique particuliere
C    -> sera modifie lors de la restructuration des zones de bord
C
          CALL UICLIM
C         ===========
     & ( NFABOR,
     &   NOZPPM, NCHARM, NCHARB, NCLPCH,
     &   IINDEF, IENTRE, IPAROI, IPARUG, ISYMET, ISOLIB,
     &   IQIMP,  ICALKE, IENTAT, IENTCP,
     &   IA(IITYPF), IA(IIZFPP), IA(IICODC),
     &   SURFBO,
     &   QIMP,   QIMPAT, QIMPCP, DH,     XINTUR,
     &   TIMPAT, TIMPCP, DISTCH, RA(IRCODC) )
C
          IF (IPPMOD(IPHPAR).EQ.0) THEN
C
C ON NE FAIT PAS DE LA PHYSIQUE PARTICULIERE NI DE MATISSE
C
            NBZFMX = NBZPPM
            NOZFMX = NOZPPM
            IILZFB = IFINIA
            IFNIA1 = IILZFB + NBZFMX
            IQCALC = IFINRA
            IFNRA1 = IQCALC + NOZFMX
            CALL IASIZE('TRIDIM',IFNIA1)
            CALL RASIZE('TRIDIM',IFNRA1)
C
            CALL STDTCL
C           ===========
     & ( IFNIA1 , IFNRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NBZFMX , NOZFMX ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IQIMP  , ICALKE , QIMP   , DH , XINTUR,
     &   IA(IICODC)      , IA(IITRIF)   , IA(IITYPF)   , IA(IIZFPP)   ,
     &   IA(IILZFB)      , IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(IRCODC)      ,
     &   RA(IW1), RA(IW2), RA(IW3), RA(IW4), RA(IW5), RA(IW6),
     &   RA(ICOEFU)      , RA(IQCALC)      ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C
        ENDIF
C
C     - Sous-programme utilisateur
C       ==========================
C
        IF (NSCAPP.EQ.0) THEN
C
          ILS    = IFINIA
          IFNIA1 = ILS + MAXELT
          CALL IASIZE('TRIDIM',IFNIA1)
C
          CALL USCLIM
C         ===========
     & ( IFNIA1 , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICODC)      , IA(IITRIF)   , IA(IITYPF)   ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(IRCODC)      ,
     &   RA(IW1), RA(IW2), RA(IW3), RA(IW4), RA(IW5), RA(IW6),
     &   RA(ICOEFU)      ,
     &   RDEVEL , RTUSER , RA     )
C
        ELSE
C
C ON FAIT DE LA PHYSIQUE PARTICULIERE (MAIS PAS DE MATISSE)
C
          CALL PPCLIM
C         ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICODC)      , IA(IITRIF)   , IA(IITYPF)   , IA(IIZFPP) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(IRCODC)      ,
     &   RA(IW1), RA(IW2), RA(IW3), RA(IW4), RA(IW5), RA(IW6),
     &   RA(ICOEFU)      ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
C     - Interface Code_Saturne
C       ======================
C
        IF(IIHMPR.EQ.1) THEN
C
          CALL UICLVE
C         ===========
     & ( NFABOR, IINDEF, IENTRE, IPAROI, IPARUG, ISYMET, ISOLIB,
     &   IA(IITYPF), IA(IIZFPP) )
C
        ENDIF
C
      ENDIF
C
C -- Methode des vortex en L.E.S. :
C    (Transfert des vortex dans les tableaux RCODCL)
C
      IF(IVRTEX.EQ.1) THEN
C
        CALL VOR2CL
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICODC)      , IA(IITRIF)      , IA(IITYPF)   ,
     &   IDEVEL , IA(IIREPV)      , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(IRCODC)      ,
     &   RA(IW1), RA(IW2), RA(IW3), RA(IW4), RA(IW5), RA(IW6),
     &   RA(ICOEFU)       ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C
C -- Methode ALE (CL de vitesse de maillage et deplacement aux noeuds)
C
      IF (IALE.EQ.1) THEN
C
        DO II = 1, NNOD
          IA(IIMPAL+II-1) = 0
        ENDDO
C
        ILS    = IFINIA
        IFNIA1 = ILS + MAXELT
        CALL IASIZE('TRIDIM',IFNIA1)
C
        CALL USALCL
C       ===========
     & ( IFNIA1 , IFINRA , ITRALE ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICODC)      , IA(IITYPF)      , IA(IIALTY)      ,
     &   IA(IIMPAL)      ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(IRCODC)      ,
     &   RA(IXYZN0)      , RA(IDEPAL)      ,
     &   RDEVEL , RTUSER , RA     )
C
C     Au cas ou l'utilisateur aurait touche DEPALE sans mettre IMPALE=1, on
C       remet le deplacement initial
        DO II  = 1, NNOD
          IF (IA(IIMPAL+II-1).EQ.0) THEN
            RA(IDEPAL+II-1       ) =
     &           XYZNOD(1,II)-RA(IXYZN0+(II-1)*NDIM  )
            RA(IDEPAL+II-1+  NNOD) =
     &           XYZNOD(2,II)-RA(IXYZN0+(II-1)*NDIM+1)
            RA(IDEPAL+II-1+2*NNOD) =
     &           XYZNOD(3,II)-RA(IXYZN0+(II-1)*NDIM+2)
          ENDIF
        ENDDO
C
C     En cas de couplage de structures, on calcule un deplacement predit
        IF (NBSTRU.GT.0) THEN
C
          CALL STRPRE
C         ===========
     & ( IFINIA , IFINRA , ITRALE , ITALIM , INEEFL ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IIMPAL)      ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RA(IFLALF), RA(IFLALB), RA(IPRALE), RA(ICOALE), RA(IDEPAL),
     &   RDEVEL , RTUSER , RA  )
C
        ENDIF
C
      ENDIF
C
C     UNE FOIS CERTAINS CODES DE CONDITIONS LIMITES INITIALISES PAR
C     L'UTILISATEUR, ON PEUT COMPLETER CES CODES PAR LES COUPLAGES
C     AUX BORDS (TYPE SYRTHES), SAUF SI ON DOIT Y REPASSER ENSUITE
C
      IF (ITRFIN.EQ.1) THEN
C
        CALL COUPBI
C       ===========
     & ( IFINIA , IFINRA ,
     &   NFABOR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IA(IICODC) ,
     &   IDEVEL , ITUSER , IA     ,
     &   RA(IRCODC) ,
     &   RDEVEL , RTUSER , RA     )
C
        IF (NFPT1T.GT.0) THEN
          CALL COU1DI
C         ===========
     & ( IFINIA , IFINRA ,
     &   NFABOR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ISVTB  , IA(IICODC) ,
     &   IDEVEL , ITUSER , IA     ,
     &   RA(IRCODC) ,
     &   RDEVEL , RTUSER , RA     )
        ENDIF
C
      ENDIF
C
C
      IF(IIRAYO.GT.0 .AND. ITRFIN.EQ.1) THEN
C
        CALL MEMRA3
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IISOTH , ITEK   , ITEXT  , ITINT  ,
     &   IFINIB , IFINRB )
C
        CALL RAYCLI
C       ===========
     & ( IFINIB , IFINRB ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISVHB  , ISVTB  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICODC) , IA(IITRIF)   , IA(IITYPF)   ,
C
     &   IA(IIZFRD) , IA(IISOTH)      ,
C
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &                                                RA(IRCODC)      ,
     &   COEFA  , COEFB  , RA(IHBORD)      , RA(ITBORD)      ,
     &   RA(IW1), RA(IW2), RA(IW3), RA(IW4), RA(IW5), RA(IW6),
C
     &   RA(ITSRE)  , RA(ITSRI)  ,
     &   RA(ITPARO) , RA(IQINCI) , RA(IXLAM) , RA(IEPA)   , RA(IEPS)  ,
     &   RA(IFNET)  , RA(IFCONV) , RA(IHCONV) , RA(ITEXT)  , RA(ITINT) ,
     &   RA(ITEK)   ,
C
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C     ON CALCULE LES COEFFICIENTS ASSOCIES AUX CONDITIONS LIMITES
C
      CALL CONDLI
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISVHB  , ISVTB  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICODC)      , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &                                                RA(IRCODC)      ,
     &   COEFA  , COEFB  , RA(IUETBO) , RA(IVSVDR)  ,
     &   RA(IHBORD)  , RA(ITBORD)      ,
     &   FRCXT  ,
     &   RA(IW1), RA(IW2), RA(IW3), RA(IW4), RA(IW5), RA(IW6),
     &                              RA(ICOEFU)      , RA(IRIJIP)      ,
     &   RDEVEL , RTUSER , RA     )
C
C
C     UNE FOIS LES COEFFICIENTS CALCULES, ON PEUT EN DEDUIRE PLUS
C     FACILEMENT (I.E. SANS RECALCULS INUTILES) LES TERMES A
C     ENVOYER POUR LES COUPLAGES AUX BORDS (TYPE SYRTHES)
C
C
C On indique si la variable couplee est l'enthalpie
      IENTHA = 0
      IF(ISVTB.GT.0) THEN
        IF(ISCSTH(ISVTB).EQ.2) THEN
          IENTHA = 1
        ENDIF
      ENDIF
C
C Compressible : on indique si la variable couple est l'energie
C
      IF ( IPPMOD(ICOMPF).GE.0 ) THEN
        IF(ISVTB.GT.0) THEN
          IF(ISCSTH(ISVTB).EQ.3) THEN
            IENTHA = 2
          ENDIF
        ENDIF
      ENDIF
C
C On recupere le Cp de la phase couplee
C  (ou de la phase 1, si pas de couplage)
      IPHAS = 1
      IF(ISVTB.GT.0) THEN
        IPHAS = IPHSCA(ISVTB)
      ENDIF
      IF(ICP(IPHAS).GT.0) THEN
        IPPCP = IPPROC(ICP(IPHAS))
        NCP   = NCELET
        CPCST = 0.D0
      ELSE
        IPPCP = IFINRA
        NCP   = 1
        CPCST = CP0(IPHAS)
      ENDIF
C
C En compressible et si on couple ave l'energie
C on recupere de Cv de la phase couplee
C
      IF ( IPPMOD(ICOMPF).GE.0 .AND. IENTHA .EQ. 2 ) THEN
C
        IPHAS = IPHSCA(ISVTB)
        IF(ICV(IPHAS).GT.0) THEN
          IPPCV = IPPROC(ICV(IPHAS))
          NCV   = NCELET
          CVCST = 0.D0
        ELSE
          IPPCV = IFINRA
          NCV   = 1
          CVCST = CV0(IPHAS)
        ENDIF
      ELSE
        IPPCV = IFINRA
        NCV   = 1
        CVCST = 0.D0
      ENDIF
C
C On envoie le tout vers SYRTHES, en distinguant CP
C  constant ou variable
      IF (ITRFIN.EQ.1) THEN
C
        CALL COUPBO
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , ISVTB  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NCP  , NCV , IENTHA ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   CPCST  , PROPCE(1,IPPCP) , CVCST  , PROPCE(1,IPPCV),
     &   RA(IHBORD)      , RA(ITBORD)      ,
     &   RDEVEL , RTUSER , RA     )
C
C
        IF (NFPT1T.GT.0) THEN
          CALL COU1DO
C         ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCP    , NFPT1D ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IENTHA , IA(IIFPT1), IA(IICLT1),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(ITPPT1), RA(ITEPT1), RA(IHEPT1),
     &   RA(IFEPT1), RA(IXLMT1), RA(IRCPT1), RA(IDTPT1),
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   CPCST  , PROPCE(1,IPPCP) , RA(IHBORD) , RA(ITBORD)  ,
     &   RDEVEL , RTUSER , RA     )
        ENDIF
      ENDIF
C
C     ON N'A PLUS BESOIN DE ISVHB OU ISVHT (POUR HBORD ET TBORD)
C     A PARTIR D'ICI
C
C
C
C     CALCUL DE LA DISTANCE A LA PAROI
C       (Nouvel algorithme. L'ancien est dans condli)
C     En ALE on ne fait ce calcul qu'a la premiere des
C       sous-iterations d'implicitation ITALIM, car le maillage
C       n'est pas modifie dans les sous-iterations suivantes
C
      IF (ITALIM.EQ.1) THEN
C
        IF(INEEDY.EQ.1.AND.IWARNY.GE.1) THEN
          CALL DMTMPS(TDIST1)
        ENDIF
C
C
C     On ne fait le calcul que s'il y a des parois (RA(IDIPAR) est reserve
C       et initialise a GRAND avant. S'il n'y a pas de paroi, il restera = GRAND)
C
C     Pour le moment, on suppose que l'on peut se contenter de faire
C       cela au premier passage, sauf avec les maillages mobiles. Attention donc
C       aux conditions aux limites variables (faces qui deviennent des parois ou
C       parois qui deviennent autre chose)
C
C     Nombre de faces de paroi
        IF(IPASS.EQ.1) THEN
          IF(INEEDY.EQ.1) THEN
            INFPAR = 0
            DO IFAC = 1, NFABOR
              IF(IA(IITYPF-1+IFAC).EQ.IPAROI .OR.
     &           IA(IITYPF-1+IFAC).EQ.IPARUG) THEN
                INFPAR = INFPAR+1
              ENDIF
            ENDDO
            IF(IRANGP.GE.0) THEN
              CALL PARCPT(INFPAR)
            ENDIF
          ENDIF
        ENDIF
C
C
C     On calcule la distance a la paroi
C          si elle doit etre mise a jour
C       et si on en a besoin,
C       et si on a choisi ce mode de calcul,
      IF( IMAJDY.EQ.0.AND.INEEDY.EQ.1.AND.ABS(ICDPAR).EQ.1) THEN
C
C     S'il n'y a pas de paroi, on garde l'initialisation a GRAND
        IF(INFPAR.EQ.0) THEN
          IMAJDY = 1
C
C     S'il y a des parois, il faut calculer
        ELSE
C
C
C     On doit conserver la memoire de memcli a cause de RA(IUETBO)
C       dans DISTYP (uniquement en LES avec van Driest mais tant pis)
C
            CALL MEMDIS
C           ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVISCF , IVISCB , IDAM   , IXAM   , ISMBR  , IROVSD ,
     &   IRTDP  , ICOFAY , ICOFBY ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    ,
     &   IFINIB , IFINRB )
C
            CALL DISTPR
C           ===========
     & ( IFINIB , IFINRB ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IITYPF)      ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(IDIPAR)      ,
     &   RA(IVISCF)      , RA(IVISCB)      ,
     &   RA(IDAM)        , RA(IXAM)        ,
     &   RA(ISMBR )      , RA(IROVSD)      ,
     &   RA(IRTDP)       , RA(ICOFAY)      , RA(ICOFBY)      ,
     &   RA(IW1), RA(IW2), RA(IW3), RA(IW4), RA(IW5), RA(IW6), RA(IW7),
     &   RA(IW8), RA(IW9),
     &   RDEVEL , RTUSER , RA     )
C
C     La distance n'a plus a etre mise a jour sauf en ALE
            IF (IALE.EQ.0) IMAJDY = 1
C
          ENDIF
        ENDIF
C
      ENDIF
C
C
C     CALCUL DE L'AMORTISSEMENT DE VAN DRIEST
C     OU CALCUL DE Y+ POUR LE LAGRANGIEN
C
C
      DO IPHAS = 1, NPHAS
C
C       Pour passer en argument
        IPHASS = IPHAS
C
C     On calcule y+ si on en a besoin
C
        IF( (ITYTUR(IPHAS).EQ.4.AND.IDRIES(IPHAS).EQ.1)
     &       .OR. (IILAGR.GE.1 .AND. IROULE.EQ.2) ) THEN
C
C       On calcule si on a demande ce mode de calcul
C               et s'il y a des parois (si pas de paroi, pas de y+)
          IF(ABS(ICDPAR).EQ.1.AND.INFPAR.GT.0) THEN
C
            IISMPH = IISYMP+NFABOR*(IPHAS-1)
C
C     On doit conserver la memoire de memcli a cause de RA(IUETBO)
C       dans DISTYP
C
            CALL MEMDYP
C           ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDAM   , IXAM   , ISMBR  , IROVSD ,
     &   IRTDP  , IDRTDP ,
     &   IQFX   , IQFY   , IQFZ   , ICOEFQ , IIRHO  , IIRHOB ,
     &   IFLUA  , IFLUB  ,
     &   ICOAX  , ICOBX  , ICOAY  , ICOBY  , ICOAZ  , ICOBZ  ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    ,
     &   IFINIB , IFINRB )
C
            CALL DISTYP
C           ===========
     & ( IFINIB , IFINRB ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IPHASS ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IITYPF) , IA(IISMPH),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(IDIPAR), PROPCE    , RA(IUETBO), RA(IYPPAR),
     &   RA(IDAM  ), RA(IXAM  ), RA(ISMBR ), RA(IROVSD),
     &   RA(IRTDP ), RA(IDRTDP),
     &   RA(IQFX  ), RA(IQFY  ), RA(IQFZ  ), RA(ICOEFQ),
     &               RA(IIRHO ), RA(IIRHOB),
     &   RA(IFLUA ), RA(IFLUB ),
     &   RA(ICOAX ), RA(ICOBX ), RA(ICOAY ), RA(ICOBY ),
     &               RA(ICOAZ ), RA(ICOBZ ),
     &   RA(IW1), RA(IW2), RA(IW3), RA(IW4), RA(IW5), RA(IW6), RA(IW7),
     &   RA(IW8), RA(IW9),
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C
        ENDIF
C
        IF (ITYTUR(IPHAS).EQ.4 .AND. IDRIES(IPHAS).EQ.1) THEN
C
C     Pas d'amortissement si pas de paroi
          IF(INFPAR.GT.0) THEN
            IF(IIFAPA(IPHAS).GT.0) THEN
              IIIFAP = IIFAPA(IPHAS)
            ELSE
              IIIFAP = IFINIB
            ENDIF
            CALL VANDRI
C           ===========
     &  (  NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NPHAS ,
     &     NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHASS ,
     &     IA(IITYPF) , IFABOR, IA(IIIFAP),
     &     IDEVEL , ITUSER , IA     ,
     &     XYZCEN , CDGFBO , RA(IUETBO) , RA(IVSVDR) , RA(IYPPAR) ,
     &     PROPCE , RDEVEL , RTUSER , RA     )
          ENDIF
C
        ENDIF
C
      ENDDO
C
      IF(INEEDY.EQ.1.AND.IWARNY.GE.1) THEN
        CALL DMTMPS(TDIST2)
        TDITOT = TDIST2-TDIST1
        IF (IRANGP.GE.0) CALL PARSOM (TDITOT)
C                        ===========
        WRITE(NFECRA,4010)TDITOT
      ENDIF
C
C
C=======================================================================
C 10. DANS LE CAS  "zero pas de temps" EN "NON SUITE" DE CALCUL
C      ON SORT ICI
C=======================================================================
C
      IF (INPDT0.EQ.1.AND.ISUITE.EQ.0) GOTO 200
C
      IF (IILAGR.EQ.3) GOTO 200
C
C=======================================================================
C 11. RESOLUTION DE LA VITESSE DE MAILLAGE EN ALE
C=======================================================================
C
      IF (IALE.EQ.1) THEN
C
        IF (ITRALE.EQ.0 .OR. ITRALE.GT.NALINF) THEN
C
          CALL MEMALE
C         ===========
     & ( IDBIA1 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVISCF , IVISCB , IDAM   , IXAM   ,
     &   IDRTP  , ISMBR  , IROVSD ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    ,
     &   IFINIA , IFINRA )
C
          CALL ALELAP
C         ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RA(IVISCF) , RA(IVISCB) ,
     &   RA(IDAM  ) , RA(IXAM  ) ,
     &   RA(IDRTP ) , RA(ISMBR ) , RA(IROVSD) ,
     &   RA(IW1   ) , RA(IW2   ) , RA(IW3   ) , RA(IW4   ),
     &   RA(IW5   ) , RA(IW6   ) , RA(IW7   ) , RA(IW8   ) , RA(IW9   ),
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
        IF (ITRALE.EQ.0) GOTO 200
C
      ENDIF
C
C=======================================================================
C 11. CALCUL A CHAMP DE VITESSE NON FIGE :
C      ON RESOUT VITESSE ET TURBULENCE
C    ON SUPPOSE QUE TOUTES LES PHASES SONT FIGEES OU AUCUNE
C=======================================================================
C
      IF (ICCVFG.EQ.0) THEN
C    ========================
C
C=======================================================================
C 12. RESOLUTION MASSE (MODULE COMPRESSIBLE UNIQUEMENT)
C=======================================================================
        IF ( IPPMOD(ICOMPF).GE.0 ) THEN
C
          IF(IWARNI(IU(1)).GE.1) THEN
            WRITE(NFECRA,1080)
          ENDIF
C
          CALL MEMCFM
C         ===========
     & ( IDBIA1 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTR   , IVISCF , IVISCB , IDAM   , IXAM   ,
     &   IDRTP  , ISMBR  , IROVSD ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    ,
     &   IW7    , IW8    , IW9    , IW10   , IW11   , IW12   ,
     &   IWFLMS , IWFLMB ,
     &   ICOEFU ,
     &   IFINIA , IFINRA )
C
          DO IPHAS = 1, NPHAS
C
            ISCAL = IRHO(IPHAS)
C
            CALL CFMSVL
C           ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS)   , NCETSM(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS))        , IA(IICESM(IPHAS))       ,
     &   IA(IITPSM(IPHAS))        ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RA(ICKUPD(IPHAS))        , RA(ISMACE(IPHAS))        ,
     &   RA(IVISCF) , RA(IVISCB)  ,
     &   RA(IDAM  ) , RA(IXAM  )  ,
     &   RA(IDRTP ) , RA(ISMBR )  , RA(IROVSD) ,
     &   RA(IW1   ) , RA(IW2   )  , RA(IW3   ) ,
     &   RA(IW4   ) , RA(IW5   )  , RA(IW6   ) ,
     &   RA(IW7   ) , RA(IW8   )  , RA(IW9   ) ,
     &   RA(IW10  ) , RA(IW11  )  , RA(IW12  ) ,
     &   RA(IWFLMS) , RA(IWFLMB)  ,
     &   RA(ICOEFU) ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDDO
C
        ENDIF
C
C=======================================================================
C 13. RESOLUTION QUANTITE DE MOUVEMENT ET MASSE (INCOMPRESSIBLE)
C=======================================================================
C
        IF(IWARNI(IU(1)).GE.1) THEN
          WRITE(NFECRA,1040)
        ENDIF
C
        CALL MEMNAV
C       ===========
     & ( IDBIA1 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVISCF , IVISCB , IVISFI , IVISBI ,
     &   IDAM   , IXAM   ,
     &   IDRTP  , IGRDP  , ISMBR  , IROVSD ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    , IW10   , IDFRCX , IFRCHY , IDFRHY ,
     &   ICOEFU , IESFLM , IESFLB , ITRAVA , IXIMPA ,
     &   IUVWK  ,
     &   IFINIA , IFINRA )
C
C     SI LE COMPRESSIBLE SANS CHOC EST ACTIF, ON RESOUT AVEC CFQDMV
        IF ( IPPMOD(ICOMPF).GE.0 ) THEN
C
          DO IPHAS = 1, NPHAS
C
            IUIPH  = IU(IPHAS)
            IFLMAS = IPPROF(IFLUMA(IUIPH))
            IFLMAB = IPPROB(IFLUMA(IUIPH))
            IPH    = IPHAS
C
            CALL CFQDMV
C           ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS)   , NCETSM(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPH    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS))        , IA(IICESM(IPHAS))       ,
     &   IA(IITPSM(IPHAS))        ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   COEFA  , COEFB  ,
     &   RA(ICKUPD(IPHAS))        , RA(ISMACE(IPHAS))        ,
     &   FRCXT  , RA(IDFRCX)      , RA(ITPUCO)      , RA(IGRDP)       ,
     &   RA(IVISCF) , RA(IVISCB)  , RA(IVISFI) , RA(IVISBI)  ,
     &   RA(IDAM  ) , RA(IXAM  )  ,
     &   RA(IDRTP ) , RA(ISMBR )  , RA(IROVSD) ,
     &   RA(IW1   ) , RA(IW2   )  , RA(IW3   ) ,
     &   RA(IW4   ) , RA(IW5   )  , RA(IW6   ) ,
     &   RA(IW7   ) , RA(IW8   )  , RA(IW9   ) ,
     &   RA(ICOEFU) ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDDO
C
        ELSE
C
C --- Boucle sur navsto pour couplage vitesse/pression
C     on s'arrete a NTERUP ou quand TOUTES les phases on converge
C
          ICVRGE = 0
          INSLST = 0
          DO ITERNS = 1, NTERUP
C
            CALL NAVSTO
C           ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , ITERNS , ICVRGE ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   TSLAGR , COEFA  , COEFB  , FRCXT  ,
     &   RA(ITRAVA) , RA(IXIMPA) , RA(IUVWK ) ,
     &   RA(IVISCF) , RA(IVISCB) , RA(IVISFI) , RA(IVISBI) ,
     &   RA(IDAM  ) , RA(IXAM  ) ,
     &   RA(IDRTP ) , RA(IGRDP ) , RA(ISMBR ) , RA(IROVSD) ,
     &   RA(IW1   ) , RA(IW2   ) , RA(IW3   ) ,
     &   RA(IW4   ) , RA(IW5   ) , RA(IW6   ) ,
     &   RA(IW7   ) , RA(IW8   ) , RA(IW9   ) , RA(IW10  ) ,
     &   RA(IDFRCX) , RA(IFRCHY) , RA(IDFRHY) ,
     &   RA(ICOEFU) , RA(IESFLM) , RA(IESFLB),
     &   RDEVEL , RTUSER , RA     )
C
C     Mise a jour de la pression si on utilise un couplage vitesse/pression
C       par point fixe
C     En parallele, l'echange est fait au debut de navsto.
            IF(NTERUP.GT.1) THEN
              DO IPHAS = 1, NPHAS
                IPRIPH  = IPR(IPHAS)
                DO IEL = 1, NCEL
                  RTPA(IEL,IPRIPH) = RTP(IEL,IPRIPH)
                ENDDO
              ENDDO
            ENDIF
C
C     Si c'est la derniere iteration : INSLST = 1
            IF((ICVRGE.EQ.1).OR.(ITERNS.EQ.NTERUP)) THEN
              INSLST = 1
            ENDIF
C
C     Si ISTMPF(IPHAS).EQ.0 (explicite) on ne traite pas le flux de
C       masse a la derniere iteration
C     Sinon on traite le flux de masse a toutes les iterations
C
C     On teste le flux de masse de la phase 1 (toutes les phases sont
C     necessairement traitees de la meme facon, cf. VERINI)
            IF( (ISTMPF(1).EQ.0.AND.INSLST.EQ.0).OR.
     &           (ISTMPF(1).NE.0                ) ) THEN
              IAPPEL = 3
              CALL SCHTMP
C             ===========
     & ( IDBIA1 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER ,
     &   RA     )
            ENDIF
C
            IF(INSLST.EQ.1) GOTO 100
C
          ENDDO
C
        ENDIF
C
 100    CONTINUE
C
C=======================================================================
C 14.  DEPLACEMENT DES STRUCTURES EN ALE ET TEST DE BOUCLAGE IMPLICITE
C=======================================================================
C
      IF (NBSTRU.GT.0) THEN
C
        CALL STRDEP
C       ===========
     & ( IDBIA1 , IDBRA1 , ITRALE , ITALIM , ITRFIN ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NVAR   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RA(IFLALF), RA(IFLALB), RA(ICOALE), RA(IPRALE),
     &   RDEVEL , RTUSER ,
     &   RA     )
C
C     On boucle eventuellement sur de deplacement des structures
        IF (ITRFIN.NE.-1) THEN
          ITALIM = ITALIM + 1
          GOTO 300
        ENDIF
C
      ENDIF
C
C --- On libere les tableaux IFLALF, IFLALB ICOALE et IPRALE
C
C     On ne passe dans SCHTMP que si ISTMPF(IPHAS).EQ.0 (explicite)
C     On teste le flux de masse de la phase 1 (toutes les phases sont
C     necessairement traitees de la meme facon, cf. VERINI)
C     pour conserver
      IF( ISTMPF(1).EQ.0 ) THEN
        IAPPEL = 4
        CALL SCHTMP
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER ,
     &   RA     )
      ENDIF
C
C=======================================================================
C 15. RESOLUTION TURBULENCE
C=======================================================================
C
        IOK = 0
        IF(IWARNI(IU(1)).GE.1) THEN
          DO IPHAS = 1, NPHAS
            IF( ITYTUR(IPHAS).EQ.2 .OR. ITYTUR(IPHAS).EQ.3
     &           .OR. ITURB(IPHAS).EQ.50 .OR. ITURB(IPHAS).EQ.60 ) THEN
              IOK = 1
            ENDIF
          ENDDO
          IF(IOK.EQ.1) THEN
            WRITE(NFECRA,1050)
          ENDIF
        ENDIF
C
        DO IPHAS = 1, NPHAS
C
C     Si on est en v2f, on reserve un tableau de taille NCELET pour
C     eviter de recalculer la production dans RESV2F (trois appels
C     a GRDCEL)
          IDBIA1 = IDEBIA
          IDBRA1 = IDEBRA
          IPRV2F = IDEBRA
          IF (ITURB(IPHAS).EQ.50) THEN
            IDBRA1 = IPRV2F + NCELET
c     Pas la peine de tester les depassements de tableaux puisqu'on
c     passe juste apres dans MEMKEP
          ENDIF
C
          IF( (ITYTUR(IPHAS).EQ.2) .OR. (ITURB(IPHAS).EQ.50) ) THEN
C
            IKIPH  = IK  (IPHAS)
C
            CALL MEMKEP
C           ===========
     & ( IDBIA1 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTR   , IVISCF , IVISCB , IDAM   , IXAM   ,
     &   IDRTP  , ISMBR  , IROVSD , ITINSK , ITINSE , IDIVU ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    ,
     &   IFINIA , IFINRA )
C
            IF(CDTVAR(IKIPH).NE.1.D0) THEN
              DO IEL = 1, NCEL
                RA(IDTR-1+IEL) = DT(IEL)*CDTVAR(IKIPH)
              ENDDO
            ELSE
              DO IEL = 1, NCEL
                RA(IDTR-1+IEL) = DT(IEL)
              ENDDO
            ENDIF
C
            CALL TURBKE
C           ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS) , NCETSM(IPHAS) ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) , IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(IDTR) , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   TSLAGR   ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)) , RA(ISMACE(IPHAS)) ,
     &   RA(IVISCF) , RA(IVISCB) , RA(IPRV2F),
     &   RA(IDAM  ) , RA(IXAM  ) ,
     &   RA(IDRTP ) , RA(ISMBR ) , RA(IROVSD) , RA(ITINSK) , RA(ITINSE),
     &   RA(IDIVU ) , RA(IW1   ) , RA(IW2   ) , RA(IW3   ) , RA(IW4   ),
     &   RA(IW5   ) , RA(IW6   ) , RA(IW7   ) , RA(IW8   ) , RA(IW9   ),
     &   RDEVEL , RTUSER , RA     )
C
            IF( ITURB(IPHAS).EQ.50 )  THEN
C
              IPHIPH  = IPHI(IPHAS)
C
              CALL MEMV2F
C             ===========
     & ( IDBIA1 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTR   , IVISCF , IVISCB , IDAM   , IXAM   ,
     &   IDRTP  , ISMBR  , IROVSD ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    , IW10   ,
     &   IFINIA , IFINRA )
C
              IF(CDTVAR(IPHIPH).NE.1.D0) THEN
                DO IEL = 1, NCEL
                  RA(IDTR-1+IEL) = DT(IEL)*CDTVAR(IPHIPH)
                ENDDO
              ELSE
                DO IEL = 1, NCEL
                  RA(IDTR-1+IEL) = DT(IEL)
                ENDDO
              ENDIF
C
              CALL RESV2F
C             ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS) , NCETSM(IPHAS) ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) , IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(IDTR) , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)) , RA(ISMACE(IPHAS)) ,
     &   RA(IVISCF) , RA(IVISCB) , RA(IPRV2F),
     &   RA(IDAM  ) , RA(IXAM  ) ,
     &   RA(IDRTP ) , RA(ISMBR ) , RA(IROVSD) ,
     &   RA(IW1   ) , RA(IW2   ) , RA(IW3   ) , RA(IW4   ),
     &   RA(IW5   ) , RA(IW6   ) , RA(IW7   ) , RA(IW8   ) , RA(IW9   ),
     &   RA(IW10  ) , RDEVEL , RTUSER , RA     )
C
            ENDIF
C
C  RELAXATION DE K ET EPSILON SI IKECOU=0 EN INSTATIONNAIRE
            IF (IKECOU(IPHAS).EQ.0 .AND. IDTVAR.GE.0) THEN
              IKIPH  = IK (IPHAS)
              IEIPH  = IEP(IPHAS)
              RELAXK = RELAXV(IKIPH)
              RELAXE = RELAXV(IEIPH)
              DO IEL = 1,NCEL
                RTP(IEL,IKIPH) = RELAXK*RTP(IEL,IKIPH)
     &                             + (1.D0-RELAXK)*RTPA(IEL,IKIPH)
                RTP(IEL,IEIPH) = RELAXE*RTP(IEL,IEIPH)
     &                             + (1.D0-RELAXE)*RTPA(IEL,IEIPH)
              ENDDO
            ENDIF
C
          ELSE IF(ITYTUR(IPHAS).EQ.3) THEN
C
            IR11IP = IR11(IPHAS)
C
            CALL MEMRIJ
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ITURB(IPHAS) ,
     &   IDTR   , IVISCF , IVISCB , ICOEFX ,
     &   IDAM   , IXAM   , IDRTP  ,
     &   ISMBR  , IROVSD , IGRDVT , IPRODU , IGRARX , IGRARY , IGRARZ ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    ,
     &   IFINIA , IFINRA )
C
            IF(CDTVAR(IR11IP).NE.1.D0) THEN
              DO IEL = 1, NCEL
                RA(IDTR-1+IEL) = DT(IEL)*CDTVAR(IR11IP)
              ENDDO
            ELSE
              DO IEL = 1, NCEL
                RA(IDTR-1+IEL) = DT(IEL)
              ENDDO
            ENDIF
C
            CALL TURRIJ
C           ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS) , NCETSM(IPHAS) ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) , IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(IDTR) , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   TSLAGR   ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)) , RA(ISMACE(IPHAS)) ,
     &   RA(IVISCF) , RA(IVISCB) , RA(ICOEFX),
     &   RA(IDAM  ) , RA(IXAM  ) ,
     &   RA(IDRTP ) , RA(ISMBR ) , RA(IROVSD) , RA(IGRDVT) ,
     &   RA(IPRODU) , RA(IGRARX) , RA(IGRARY) , RA(IGRARZ) ,
     &   RA(IW1   ) , RA(IW2   ) , RA(IW3   ) , RA(IW4   ) ,
     &   RA(IW5   ) , RA(IW6   ) , RA(IW7   ) , RA(IW8   ) , RA(IW9   ),
     &   RDEVEL , RTUSER , RA     )
C
          ELSE IF( ITURB(IPHAS).EQ.60 ) THEN
C
            IKIPH  = IK  (IPHAS)
C
            CALL MEMKOM
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTR   , IVISCF , IVISCB , IDAM   , IXAM   ,
     &   IDRTP  , ISMBR  , IROVSD , ITINSK , ITINSE , IDIVU  ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    ,
     &   IFINIA , IFINRA )
C
            IF(CDTVAR(IKIPH).NE.1.D0) THEN
              DO IEL = 1, NCEL
                RA(IDTR-1+IEL) = DT(IEL)*CDTVAR(IKIPH)
              ENDDO
            ELSE
              DO IEL = 1, NCEL
                RA(IDTR-1+IEL) = DT(IEL)
              ENDDO
            ENDIF
C
            CALL TURBKW
C           ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS) , NCETSM(IPHAS) ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) , IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(IDTR) , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   TSLAGR   ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)) , RA(ISMACE(IPHAS)) ,
     &   RA(IS2KW(IPHAS)), RA(IDVUKW(IPHAS)),
     &   RA(IVISCF) , RA(IVISCB) ,
     &   RA(IDAM  ) , RA(IXAM  ) ,
     &   RA(IDRTP ) , RA(ISMBR ) , RA(IROVSD) , RA(ITINSK) , RA(ITINSE),
     &   RA(IDIVU ) , RA(IW1   ) , RA(IW2   ) , RA(IW3   ) , RA(IW4   ),
     &   RA(IW5   ) , RA(IW6   ) , RA(IW7   ) , RA(IW8   ) , RA(IW9   ),
     &   RDEVEL , RTUSER , RA     )
C
C  RELAXATION DE K ET OMEGA SI IKECOU=0
            IF (IKECOU(IPHAS).EQ.0 .AND. IDTVAR.GE.0) THEN
              IKIPH   = IK  (IPHAS)
              IOMIPH  = IOMG(IPHAS)
              RELAXK = RELAXV(IKIPH )
              RELAXW = RELAXV(IOMIPH)
              DO IEL = 1,NCEL
                RTP(IEL,IKIPH)  = RELAXK*RTP(IEL,IKIPH)
     &                             + (1.D0-RELAXK)*RTPA(IEL,IKIPH)
                RTP(IEL,IOMIPH) = RELAXW*RTP(IEL,IOMIPH)
     &                             + (1.D0-RELAXW)*RTPA(IEL,IOMIPH)
              ENDDO
            ENDIF
C
          ENDIF
C
        ENDDO
C
C
C  FIN DU TEST SUR ICCVFG : CI DESSOUS
C ========================
C
      ENDIF
C
C
C     Ici on peut liberer les eventuels tableaux SKW et DIVUKW
C
C=======================================================================
C 16.  RESOLUTION DES SCALAIRES
C=======================================================================
C
C
      IF (NSCAL.GE.1 .AND. IIRAYO.GT.0) THEN
C
        IF(IWARNI(IU(1)).GE.1 .AND. MOD(NTCABS,NFREQR).EQ.0) THEN
          WRITE(NFECRA,1070)
        ENDIF
C
        CALL MEMRA2
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTR   , IVISCF , IVISCB , IDAM   , IXAM   ,
     &   IDRTP  , ISMBR  , IROVSD ,
C
     &   ICORUA , ICORUB , IFLXMA , IFLXMB , ITEK   ,
C
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    , IW10   ,
     &   IFINIA , IFINRA )
C
C
        CALL RAYDOM
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , IA(IITYPF) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IIZFRD) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
C
     &   RA(ICORUA) , RA(ICORUB) ,
     &   RA(IFLXMA) , RA(IFLXMB) ,
C
     &   RA(IDTR)   , RA(IVISCF) , RA(IVISCB) ,
     &   RA(IDAM  ) , RA(IXAM  ) ,
     &                RA(IDRTP ) , RA(ISMBR ) , RA(IROVSD) ,
     &   RA(IW1   ) , RA(IW2   ) , RA(IW3   ) , RA(IW4   ) , RA(IW5   ),
     &   RA(IW6   ) , RA(IW7   ) , RA(IW8   ) , RA(IW9   ) ,
     &   RA(IW10  ) ,
     &   RDEVEL , RTUSER ,
C
     &   RA(ITSRE)  , RA(ITSRI)  , RA(IQX)    , RA(IQY)    , RA(IQZ)   ,
     &   RA(IABS)   , RA(IEMI)   , RA(ICAK)   , RA(ITEK)   ,
     &   RA(ITPARO) , RA(IQINCI) , RA(IXLAM)  , RA(IEPA)   , RA(IEPS)  ,
     &   RA(IFNET)  , RA(IFCONV) , RA(IHCONV) ,
C
     &   RA     )
C
      ENDIF
C
C
      IF (NSCAL.GE.1) THEN
C
        IF(IWARNI(IU(1)).GE.1) THEN
          WRITE(NFECRA,1060)
        ENDIF
C
        CALL MEMSCA
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTR   , IVISCF , IVISCB , IDAM   , IXAM   ,
     &   IDRTP  , ISMBR  , IROVSD ,
     &   IW1    , IW2    , IW3    , IW4    , IW5    , IW6    , IW7    ,
     &   IW8    , IW9    ,
     &   IFINIA , IFINRA )
C
        CALL SCALAI
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   TSLAGR , COEFA  , COEFB  ,
     &   RA(IDTR)   , RA(IVISCF) , RA(IVISCB) ,
     &   RA(IDAM  ) , RA(IXAM  ) ,
     &                RA(IDRTP ) , RA(ISMBR ) , RA(IROVSD) ,
     &   RA(IW1   ) , RA(IW2   ) , RA(IW3   ) , RA(IW4   ) , RA(IW5   ),
     &   RA(IW6   ) , RA(IW7   ) , RA(IW8   ) , RA(IW9   ) ,
     &   RDEVEL     , RTUSER ,
     &   RA     )
C
      ENDIF
C
C=======================================================================
C 17.  TRAITEMENT DU FLUX DE MASSE, DE LA VISCOSITE,
C      DE LA MASSE VOLUMIQUE ET DE LA CHALEUR SPECIFIQUE POUR
C      UN THETA SCHEMA
C=======================================================================
C
C
      IAPPEL = 5
      CALL SCHTMP
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C
C=======================================================================
C 18.  SORTIE DANS LE CAS DE "zero pas de temps" ET INIT ALE
C=======================================================================
C
 200  CONTINUE
C
C
C=======================================================================
C
C--------
C FORMATS
C--------
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  INITIALISATIONS                                            ',/,
     &'  ===============                                            ',/)
 1010 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  CALCUL DES GRANDEURS PHYSIQUES                             ',/,
     &'  ==============================                             ',/)
 1020 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  CALCUL DU CFL, DU FOURIER ET DU DT VARIABLE                ',/,
     &'  ===========================================                ',/)
 1030 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  MISE EN PLACE DES CONDITIONS AUX LIMITES                   ',/,
     &'  ========================================                   ',/)
 1040 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  RESOLUTION DES EQUATIONS DE NAVIER-STOKES                  ',/,
     &'  =========================================                  ',/)
 1050 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  RESOLUTION DES EQUATIONS DES VARIABLES TURBULENTES         ',/,
     &'  ==================================================         ',/)
 1060 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  RESOLUTION DES EQUATIONS SUR L''ENERGIE ET LES SCALAIRES   ',/,
     &'  ========================================================   ',/)
 1070 FORMAT(/,
     & '------------------------------------------------------------',/,
     &                                                              /,/,
     & ' RESOLUTION DES TRANSFERTS THERMIQUES RADIATIFS             ',/,
     &'  ==============================================             ',/)
 1080 FORMAT(/,
     & '------------------------------------------------------------',/,
     &                                                              /,/,
     & ' RESOLUTION DE L''EQUATION DE MASSE                         ',/,
     &'  ==================================                         ',/)
C
 4010 FORMAT(/,
     &' ** TEMPS CPU TOTAL POUR LA DISTANCE A LA PAROI : ',E14.5     ,/,
     &'    -------------------------------------------              ',/)
C
#else
C
 1000 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  INITIALISATIONS                                            ',/,
     &'  ===============                                            ',/)
 1010 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  COMPUTATION OF PHYSICAL QUANTITIES                         ',/,
     &'  ==================================                         ',/)
 1020 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  COMPUTATION OF CFL, FOURIER AND VARIABLE DT                ',/,
     &'  ===========================================                ',/)
 1030 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  SETTING UP THE BOUNDARY CONDITIONS                         ',/,
     &'  ==================================                         ',/)
 1040 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  SOLVING NAVIER-STOKES EQUATIONS                            ',/,
     &'  ===============================                            ',/)
 1050 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  SOLVING TURBULENT VARIABLES EQUATIONS                      ',/,
     &'  =====================================                      ',/)
 1060 FORMAT(/,
     &' ------------------------------------------------------------',/,
     &                                                              /,/,
     &'  SOLVING ENERGY AND SCALARS EQUATIONS                       ',/,
     &'  ====================================                       ',/)
 1070 FORMAT(/,
     & '------------------------------------------------------------',/,
     &                                                              /,/,
     & ' SOLVING THERMAL RADIATIVE TRANSFER                         ',/,
     &'  ==================================                         ',/)
 1080 FORMAT(/,
     & '------------------------------------------------------------',/,
     &                                                              /,/,
     & ' SOLVING MASS EQUATION                                      ',/,
     &'  =====================                                      ',/)
C
 4010 FORMAT(/,
     &' ** TOTAL CPU TIME FOR THE WALL DISTANCE: ',E14.5             ,/,
     &'    -------------------------------------                    ',/)
C
#endif
C
C----
C FIN
C----
C
      END
c@z
