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
                        SUBROUTINE INIINI
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  INITIALISATION PAR DEFAUT DES COMMONS
CFONC    AVANT DE PASSER LA MAIN A L'UTILISATEUR
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
      INCLUDE "cstnum.h"
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "vector.h"
      INCLUDE "albase.h"
      INCLUDE "alstru.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "mltgrd.h"
      INCLUDE "ihmpre.h"
      INCLUDE "matiss.h"
C
C***********************************************************************
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ, IVAR, ISCAL, IPHAS, IPROP, IEST, IMOM
      INTEGER          ISTR
C
C***********************************************************************
C
C=======================================================================
C 1. STOCKAGE DES ARGUMENTS ET IMPRESSIONS INITIALES
C=======================================================================
C
      WRITE(NFECRA, 900)
C
#if defined(_CS_LANG_FR)
C
 900  FORMAT(/,
     &'***************************************************************',
     &/,/,
     &'                   PREPARATION DU CALCUL                     ',/,
     &'                   =====================                     ',/,
     &                                                                /,
     &                                                                /,
     &' =========================================================== ',/,
     &                                                                /,
     &                                                                /)
C
#else
C
 900  FORMAT(/,
     &'***************************************************************',
     &/,/,
     &'                   CALCULATION PREPARATION'                   ,/,
     &'                   ======================='                   ,/,
     &                                                                /,
     &                                                                /,
     &' ===========================================================' ,/,
     &                                                                /,
     &                                                                /)
C
#endif
C
C=======================================================================
C 2. ENTREES SORTIES entsor.h
C=======================================================================
C
C ---> Impressions standard
C         -10000 : non initialise : voir modini
      DO II = 1, NVARMX
        IWARNI(II) = -10000
      ENDDO
C
C ---> NFECRA vaut 6 par defaut ou 9 en parallele (CSINIT)
C
C ---> Geometrie
C
      IMPGEO = 10
      FICGEO = 'geomet'
C
C ---> Fichier amont
C
      FICAMO = 'suiamo'
      FICAMX = 'suiamx'
C
C    Module thermique 1D en paroi
      FICMT1 = 't1damo'
C
C    Methode des vortex : on utilise la meme unite
C                  Pour le fichier de donnees, on specifie l'unite ici, mais le
C                  nom est laisse dans usvort. On utilise 20 et pas 11
C                  car en cas d'entree multiple, les deux fichiers doivent etre
C                  ouverts en meme temps dans VORINI
      IMPMVO = 11
      FICMVO = 'voramo'
      IMPDVO = 20
C
C    Rayonnement
      FICAMR = 'rayamo'
C
C    Lagrangien
C      Fichier suite calcul
      FICAML = 'lagamo'
C      Fichier suite statistique
      FICMLS = 'lasamo'
C
C ---> Fichier stop
C
      IMPSTP = 12
      FICSTP = 'ficstp'
C
C ---> Fichier aval
C
C     NTSUIT : Periode de sauvegarde du fichier suite
C            (il est de toutes facons ecrit en fin de calcul)
C              -1   : a la fin seulement
C              0    : par defaut (4 fois par calcul)
C              > 0  : periode
C
C
      FICAVA = 'suiava'
      FICAVX = 'suiavx'
C
      NTSUIT = 0
C
C    Module thermique 1D en paroi
      FICVT1 = 't1dava'
C
C    Methode des vortex : on utilise la meme unite
C                  et le format ascii obligatoirement
C                  (pas de detection automatique, fichiers de taille faible)
      IMPVVO = 20
      FICVVO = 'vorava'
C
C    Rayonnement : on utilise la meme unite (pour le moment)
C                  et le meme format (voir modini)
      FICAVR = 'rayava'
C
C    Lagrangien
C      Fichier suite calcul
      FICAVL = 'lagava'
C      Fichier suite statistique
      FICVLS = 'lasava'
C
C    Fichier listing Lagrangien
C
      IMPLAL = 80
      FICLAL = 'listla'
      NTLAL  = 1
C
C    Fichier historique Lagrangien
C
      IMPLI1 = 81
      IMPLI2 = 82
C
C
C    Post traitement
C
C       ICHRVL : Post traitement du domaine fluide
C         (1 oui, 0 non)
C       ICHRBO : Post traitement du bord du domaine
C         (1 oui, 0 non)
C       ICHRSY : Post traitement des zones couplees avec Syrthes
C         (1 oui, 0 non)
C
C       ICHRMD : indique si les maillages ecrits seront :
C         0 : fixes,
C         1 : deformables a topologie constante,
C         2 : modifiables (pourront etre completement redefinis en
C             cours de calcul via le sous-programme USMPST).
C        10 : comme INDMOD = 0, avec champ de déplacement
C        11 : comme INDMOD = 1, avec champ de déplacement
C        12 : comme INDMOD = 2, avec champ de déplacement
C
C       NTCHR  : Periode   de sortie Post
C         -1 : une seule sortie a la fin
C         >0 : periode
C       ICHRVR  : Variables a sortir
C         (1 oui, 0 non, sinon : non initialise)
C
      ICHRVL = 1
      ICHRBO = 0
      ICHRSY = 0
C
      ICHRMD = 0
      NTCHR  = -1
C
      DO II = 1, NVPPMX
        ICHRVR(II )   = -999
      ENDDO
C
C       FMTCHR : format ('EnSight Gold', 'MED_fichier', ou 'CGNS')
C       OPTCHR : options associees au format de sortie
C
      FMTCHR = 'EnSight Gold'
      OPTCHR = 'binary'
C
C ---> Fichier thermochinie
C        FPP : utilisateur
C        JNF : Janaf
C        Les deux fichiers peuvent partager la meme unite
C          puisqu'ils sont lus l'un a pres l'autre.
C      En prime, INDJON (janaf=1 ou non=0)
C
      IMPFPP = 25
      FICFPP = 'dp_tch'
C
      IMPJNF = IMPFPP
      FICJNF = 'JANAF'
C
      INDJON = 1
C
C ---> Fichiers module atmospherique
      IMPMET = 26
      FICMET = 'meteo'
C
C ---> Fichier resuMatisse
      IMPMAT = 28
C
C
C ---> Fichiers historiques
C
C     IMPHIS : fichier stock + unite d'ecriture des variables
C     EMPHIS : EMPlacement
C     EXTHIS : EXTension
C     IMPUSH : Unite fichiers specifiques ushist
C     FICUSH : Nom   fichiers specifiques ushist
C     IMPSTH : fichier stock + unite d'ecriture des variables
C              des structures mobiles
C
      IMPHIS(1) = 30
      IMPHIS(2) = 31
      IMPSTH(1) = 32
C     Les fichiers d'ecriture ne sont pas ouverts en meme temps
      IMPSTH(2) = IMPHIS(2)
C
      EMPHIS = './'
      EXTHIS = 'hst'
C
      DO II = 1, NUSHMX
        IMPUSH(II) = 32+II
        IF (IRANGP .LE. 0) THEN
           WRITE(FICUSH(II),'(1A3,I3.3)')'ush',II
        ELSE
           WRITE(FICUSH(II),'(1A3,I3.3,1A3,I4.4)')'ush',II,
     &                                            '.n_',IRANGP+1
        ENDIF
      ENDDO
C
C     NCAPT  : nombre de sondes total (limite a NCAPTM)
C     NTHIST : Periode de sortie
C         ( > 0 ou -1 (jamais))
C     NTHSAV : Periode de sauvegarde
C         ( > 0 ou -1 (a la fin) ou 0 (NTMABS/4))
C     IHISVR : nb de sonde et numero par variable
C         (-999 non initialise)
C     IHISTR : indicateur d'ecriture des historiques des structures
C              mobiles internes (=0 ou 1)
C     NCAPT  : nombre de sondes total (limite a NCAPTM)
C     NODCAP : element correspondant aux sondes
C     NDRCAP : rang du processus contenant NODCAP (parallelisme)
C     XYZCAP : position demandee des sondes
C
      NCAPT = 0
C
      NTHIST = 1
      NTHSAV = 0
C
      DO II = 1, NVPPMX
        DO JJ = 1, NCAPTM+1
          IHISVR(II ,JJ) = -999
        ENDDO
      ENDDO
C
      IHISTR = 0
C
      DO II = 1, NCAPTM
        NODCAP(II) = 1
        NDRCAP(II) = -1
      ENDDO
C
      DO II = 1, NCAPTM
        XYZCAP(1,II) = 0.D0
        XYZCAP(2,II) = 0.D0
        XYZCAP(3,II) = 0.D0
      ENDDO
C
C
C
C ---> Fichiers Lagrangiens
C
C     IMPLA1 : Unite fichier specifique Lagrangien
C     IMPLA2 : Unite fichier specifique Lagrangien
C     IMPLA3 : Unite fichier SCRATCH pour stockage temporaire
C     IMPLA4 : Unite fichier SCRATCH pour stockage temporaire
C     IMPLA5 : Unite d'ecriture des variables associees
C
      IMPLA1 = 50
      IMPLA2 = 51
      IMPLA3 = 52
      IMPLA4 = 53
C
      IMPLA5(1)  = 54
      IMPLA5(2)  = 55
      IMPLA5(3)  = 56
      IMPLA5(4)  = 57
      IMPLA5(5)  = 58
      IMPLA5(6)  = 59
      IMPLA5(7)  = 60
      IMPLA5(8)  = 61
      IMPLA5(9)  = 62
      IMPLA5(10) = 63
      IMPLA5(11) = 64
      IMPLA5(12) = 65
      IMPLA5(13) = 66
      IMPLA5(14) = 67
      IMPLA5(15) = 68
C
C ---> Fichiers utilisateurs
C
      DO II = 1, NUSRMX
        IMPUSR(II) = 69+II
        IF (IRANGP .LE. 0) THEN
           WRITE(FICUSR(II),'(1A4,I2.2)')'usrf',II
        ELSE
           WRITE(FICUSR(II),'(1A4,I2.2,1A3,I4.4)')'usrf',II,
     &                                            '.n_',IRANGP+1
        ENDIF
      ENDDO
C
C ---> Sorties listing
C
C   COMMUNES
C     IPP*   : Pointeurs de reperage des variables pour les sorties
C              1 pointe sur une case poubelle et constitue donc une
C              bonne initialisation
C     NOMVAR : Nom des variables
C     ILISVR : on suit la variable (1) ou non (0) ou non initialise
C     ITRSVR : numero de variable si IPP correspond a une variable resolue (p,u,k...)
C              0 si IPP correspond a une variable annexe (cp, mut...)ou a rien
C     NTLIST : periode d'ecriture
C       ( -1 : dernier pas de temps : > 0 : periode)
C
      DO II = 1, NVARMX
        IPPRTP(II) = 1
      ENDDO
      DO II = 1, NPROMX
        IPPPRO(II) = 1
      ENDDO
      IPPDT        = 1
      IPPTX        = 1
      IPPTY        = 1
      IPPTZ        = 1
      DO II = 1, NVPPMX
        IPP2RA(II) = 0
      ENDDO
C
      DO II = 1, NVPPMX
        NOMVAR(II)    = ' '
        ILISVR(II)    = -999
        ITRSVR(II)    = 0
      ENDDO
C
      NTLIST = 1
C
C   PARAMETRES DE SUIVI DE CALCUL, MIN-MAX, CLIPMIN, CLIPMAX
C
      DO II = 1, NVPPMX
        ICLPMN(II) = 0
        ICLPMX(II) = 0
      ENDDO
C
      DO II = 1, NVPPMX
        VARMIN(II) = 0.D0
        VARMAX(II) = 0.D0
        VARMNA(II) = 0.D0
        VARMXA(II) = 0.D0
      ENDDO
C
C   PARAMETRES DE CONVERGENCE, NORME DU SECOND MEMBRE, NOMBRE ITERATIONS
C                                RESIDU NORME, DERIVE
      DO II = 1, NVPPMX
        NBIVAR(II) = 0
      ENDDO
      DO II = 1, NVPPMX
        RNSMBR(II) = 0.D0
        RESVAR(II) = 0.D0
        DERVAR(II) = 0.D0
      ENDDO
C
C   PARAMETRES DU PAS DE TEMPS LOCAL
C
      DO II = 1, 8
        DO JJ = 1, 4
          PTPLOC(II,JJ) = 0.D0
        ENDDO
      ENDDO
C
C
C ---> Post traitement automatique (bord)
C
      IF (IFOENV .EQ. 0) THEN
        IPSTDV = 0
      ELSE
        IPSTDV = IPSTYP*IPSTCL*IPSTFT*IPSTFO
      ENDIF
C
C
C ---> CPU
C      TMARUS : marge (Arret du calcul avant limite CPU)
C        Si TMARUS negatif, le code calcule une marge seul
C        Sinon, il utilise TMARUS (donnee en secondes)
C
      TMARUS = -1.D0
C
C
C Ici entsor.h est completement initialise
C
C=======================================================================
C 3. DIMENSIONS DE dimens.h (GEOMETRIE, sauf NCELBR)
C=======================================================================
C
C---> GEOMETRIE
C
C---> LECTURE SELON UTILISATION PREPROCESSEUR OU ANCIEN FICHIER
C
         NDIM   = 0
         NCEL   = 0
         NCELET = 0
         NFAC   = 0
         NFABOR = 0
         NFML   = 0
         NPRFML = 0
         NNOD   = 0
         LNDFAC = 0
         LNDFBR = 0
         NCELBR = 0
C
         NCELGB = 0
         NFACGB = 0
         NFBRGB = 0
         NSOMGB = 0
C
C Par defaut, on suppose qu'il n'y a pas de periodicite
         IPERIO = 0
         IPEROT = 0
C
      IF (IFOENV .EQ. 0) THEN
C
C        NECESSITE FICGEO
C
         CALL LEDGEO
C        ===========
     &        ( NDIM   , NCELET , NCEL   , NFAC   , NFABOR ,
     &          NPRFML , NFML   , NNOD   , LNDFAC , LNDFBR )
C
C
C --- Pas de decoupage en parallele prevu sans preprocesseur.
C     NCELET=NCEL est fait dans ledgeo
C
      ELSE IF (IFOENV .EQ. 1) THEN
C
C --- PASSAGE DANS LEDEVI :
C        On ne lit ici que les dimensions.
C        Les tableaux seront lus apres allocation.
C        On initialise les dimensions a 0 avant la lecture, certaines
C          rubriques etant optionnelles.
C
         CALL LEDEVI
C        ===========
     &   ( NDIM   , NFML   , NPRFML , IPERIO , IPEROT )
C
      ENDIF
C
C
C=======================================================================
C 4. DIMENSIONS de dimens.h (PHYSIQUE)
C=======================================================================
C
C --- Nombre de phases, de scalaires, de scalaires a diffusivite
C            variable, de variables
C
      NPHAS  = 0
      NSCAL  = 0
      NSCAUS = 0
      NSCAPP = 0
      NVAR   = 0
C
C --- Nombre de proprietes physiques (utile ?)
C
      NPROCE = 0
      NPROFA = 0
      NPROFB = 0
      NFLUMA = 0
C
C --- Nombre de tableaux NFABOR de type COEFA (ou COEFB)
C
      NCOFAB = 0
C
C
C=======================================================================
C 5. POSITION DES VARIABLES DE numvar.h
C=======================================================================
C
C --- Variables de calcul resolues (RTP, RTPA)
C
      DO IPHAS = 1, NPHSMX
        IPR   (IPHAS) = 0
        IU    (IPHAS) = 0
        IV    (IPHAS) = 0
        IW    (IPHAS) = 0
        IK    (IPHAS) = 0
        IEP   (IPHAS) = 0
        IR11  (IPHAS) = 0
        IR22  (IPHAS) = 0
        IR33  (IPHAS) = 0
        IR12  (IPHAS) = 0
        IR13  (IPHAS) = 0
        IR23  (IPHAS) = 0
        IPHI  (IPHAS) = 0
        IFB   (IPHAS) = 0
      ENDDO
C
      DO ISCAL = 1, NSCAMX
        ISCA  (ISCAL) = 0
        ISCAPP(ISCAL) = 0
      ENDDO
C
C --- Initialisation par defaut des commons pour la physique particuliere
C
      CALL PPINII
C     ===========
C
C --- Proprietes physiques au sens large (PROPCE, PROFA, PROPFB)
C
      DO IPROP  = 1, NPROMX
        IPPROC(IPROP) = 0
        IPPROF(IPROP) = 0
        IPPROB(IPROP) = 0
      ENDDO
C
      DO IPHAS = 1, NPHSMX
        IROM  (IPHAS) = 0
        IVISCL(IPHAS) = 0
        IVISCT(IPHAS) = 0
        ICOUR (IPHAS) = 0
        IFOUR (IPHAS) = 0
        ICP   (IPHAS) = 0
        IPRTOT(IPHAS) = 0
      ENDDO
C
      DO ISCAL = 1, NSCAMX
        IVISLS(ISCAL) = -1
      ENDDO
C
      DO IVAR  = 1, NVARMX
        IFLUMA(IVAR) = 0
        IFLUAA(IVAR) = 0
      ENDDO
C
C --- Conditions aux limites (COEFA, COEFB)
C
      ICOEF  = 1
      ICOEFF = 2
C
      DO IVAR = 1, NVARMX
        ICLRTP(IVAR,ICOEF ) = 0
        ICLRTP(IVAR,ICOEFF) = 0
      ENDDO
C
C --- Ici tout numvar est initialise.
C
C=======================================================================
C 6. POSITION DES VARIABLES DE pointe.h
C=======================================================================
C
C --- Geometrie
C
      IICELB = 0
      ISRFAN = 0
      ISRFBN = 0
      IDIST  = 0
      IDISTB = 0
      IPOND  = 0
      IDIJPF = 0
      IDIIPB = 0
      IDOFIJ = 0
C
C --- Auxiliaires independants du nb phases
C
      ICOCG  = 0
      ICOCGB = 0
      ITPUCO = 0
      IDIPAR = 0
      IYPPAR = 0
      NFPT1D = 0
      NMXT1D = 0
      INPPT1 = 0
      IIFPT1 = 0
      IICLT1 = 0
      IEPPT1 = 0
      IRGPT1 = 0
      ITPPT1 = 0
      ITEPT1 = 0
      IHEPT1 = 0
      IFEPT1 = 0
      IXLMT1 = 0
      IRCPT1 = 0
      IDTPT1 = 0
C
C --- Autres auxiliaires dependant du nb de phases
C
      IITYPF = 0
      IITRIF = 0
      IISYMP = 0
C
      DO IPHAS = 1, NPHSMX
        IIFAPA(IPHAS) = 0
        NCEPDC(IPHAS) = 0
        IICEPD(IPHAS) = 0
        ICKUPD(IPHAS) = 0
        NCETSM(IPHAS) = 0
        IICESM(IPHAS) = 0
        IITPSM(IPHAS) = 0
        ISMACE(IPHAS) = 0
        IS2KW (IPHAS) = 0
        IDVUKW(IPHAS) = 0
      ENDDO
C
C --- Auxiliaires pour la periodicite
C
      IDUDXY =0
      IDRDXY =0
      IWDUDX =0
      IWDRDX =0
C
C --- Ici tout pointe.h est initialise
C
C=======================================================================
C 7. OPTIONS DU CALCUL : TABLEAUX DE optcal.h
C=======================================================================
C
C --- Definition des equations
C       (convection-diffusion instationnaire, avec dirichlet
C        sauf pour la pression, diffusion instationnaire)
C        IDIFFT en particulier multiplie la diffusion turbulente
C           (quand elle est activee par le modele)
C
      DO II = 1, NVARMX
        ICONV (II) = 1
        ISTAT (II) = 1
        IDIFF (II) = 1
        IDIFFT(II) = 1
      ENDDO
C
C --- Schema en temps
C
C     NOTER BIEN que les valeurs de THETA pour
C       rho, visc, cp, flux de masse, termes sources
C       sont renseignees a partir des indicateurs I..EXT et ISTMPF
C       Elles ne sont pas completees directement par l'utilisateur
C       Les tests dans l'algo portent sur ces indicateurs
C
C
C   -- Schema en temps (regroupera les options suivantes)
C     = 1 Standard ordre 1
C     = 2 Standard ordre 2
C     si on veut, on peut en rajouter.
C     Cette variable conditionne toutes les autres de maniere automatique,
C     pour definir des schemas coherents dans modini.
      DO IPHAS = 1, NPHSMX
        ISCHTP(IPHAS) = -999
      ENDDO
C
C   -- Variables
C     Pour un schema centre (en n+1/2) il faut prendre theta = 0.5
C     On devrait pouvoir faire de l'explicite pur avec theta = 0 mais
C       ce point reste a voir
C     Ce theta sert aux termes de convection diffusion (partie implicitee
C       d'ordinaire)
C     Il est applique sous la forme (1-theta) ancien + theta nouveau
C       (ce n'est pas une extrapolation, contrairement aux termes sources)
      DO II = 1, NVARMX
        THETAV(II) =-999.D0
      ENDDO
C
      DO IPHAS = 1, NPHSMX
C
C   -- Flux de masse (-999 = non initialise)
C     = 1 Standard d'ordre 1 (THETFL = -999 inutile)
C     = 0 Explicite (THETFL = 0)
C     = 2 Ordre deux (THETFL = 0.5)
        ISTMPF(IPHAS) = -999
        THETFL(IPHAS) =-999.D0
C
C   -- Termes sources Navier Stokes
C     Pour les termes sources explicites en std, I..EXT definit
C       l'extrapolation -theta ancien + (1+theta) nouveau
C     = 0 explicite
C     = 1 extrapolation avec theta = 1/2
C     = 2 extrapolation avec theta = 1
C       0 implique pas de reservation de tableaux
C       1 et 2 sont deux options equivalentes, la difference etant faite
C       uniquement au moment de fixer theta
C     Pour les termes sources implicites en std, I..EXT definit
C       la mise a l'ordre 2 ou non avec le thetav de la variable associee
C     = 0 implicite (std)
C     > 0 utilisation du thetav
C     Noter cpdt que le TS d'acc. masse n'est pas regi par I..EXT
C       (il suit bilsc2)
        ISNO2T(IPHAS) = -999
        THETSN(IPHAS) =-999.D0
C   -- Termes sources Grandeurs turbulentes
C     Pour les termes sources explicites en std, I..EXT definit
C       l'extrapolation -theta ancien + (1+theta) nouveau
C     = 0 explicite
C     = 1 extrapolation avec theta = 1/2
C     = 2 extrapolation avec theta = 1
C       0 implique pas de reservation de tableaux
C       1 et 2 sont deux options equivalentes, la difference etant faite
C       uniquement au moment de fixer theta
C     Pour les termes sources implicites en std, I..EXT definit
C       la mise a l'ordre 2 ou non avec le thetav de la variable associee
C     = 0 implicite (std)
C     > 0 utilisation du thetav
C     Noter cpdt que le TS d'acc. masse n'est pas regi par I..EXT
C       (il suit bilsc2)
        ISTO2T(IPHAS) = -999
        THETST(IPHAS) =-999.D0
C
C    -- Proprietes physiques
C     I..EXT definit l'extrapolation -theta ancien + (1+theta) nouveau
C     = 0 explicite
C     = 1 extrapolation avec theta = 1/2
C     = 2 extrapolation avec theta = 1
C       0 implique pas de reservation de tableaux
C       1 et 2 sont deux options equivalentes, la difference etant faite
C       uniquement au moment de fixer theta
C     INIT.. =1 indique que la variable a ete proprement initialisee (dans un
C       fichier suite portant les valeurs adaptees)
C
C     Masse volumique
        IROEXT(IPHAS) = -999
        THETRO(IPHAS) = -999.D0
        INITRO(IPHAS) = 0
C     Viscosite totale
        IVIEXT(IPHAS) = -999
        THETVI(IPHAS) = -999.D0
        INITVI(IPHAS) = 0
C     Chaleur specifique
        ICPEXT(IPHAS) = -999
        THETCP(IPHAS) = -999.D0
        INITCP(IPHAS) = 0
C
C   -- Convergence point fixe vitesse pression
        EPSUP (IPHAS) = 1.D-5
C
C   -- Tab de travail pour normes de navsto
        XNRMU0(IPHAS) = 0.D0
        XNRMU (IPHAS) = 0.D0
C
      ENDDO
C
C   -- Nb d'iter point fixe vitesse pression
      NTERUP = 1
C
      DO ISCAL = 1, NSCAMX
C
C   -- Termes sources des scalaires
C     Pour les termes sources explicites en std, I..EXT definit
C       l'extrapolation -theta ancien + (1+theta) nouveau
C     = 0 explicite
C     = 1 extrapolation avec theta = 1/2
C     = 2 extrapolation avec theta = 1
C       0 implique pas de reservation de tableaux
C       1 et 2 sont deux options equivalentes, la difference etant faite
C       uniquement au moment de fixer theta
C     Pour les termes sources implicites en std, I..EXT definit
C       la mise a l'ordre 2 ou non avec le thetav de la variable associee
C     = 0 implicite (std)
C     > 0 utilisation du thetav
C     Noter cpdt que le TS d'acc. masse n'est pas regi par I..EXT
C       (il suit bilsc2)
        ISSO2T(ISCAL) = -999
        THETSS(ISCAL) =-999.D0
C
C    -- Proprietes physiques
C     I..EXT definit l'extrapolation -theta ancien + (1+theta) nouveau
C     = 0 explicite
C     = 1 extrapolation avec theta = 1/2
C     = 2 extrapolation avec theta = 1
C       0 implique pas de reservation de tableaux
C       1 et 2 sont deux options equivalentes, la difference etant faite
C       uniquement au moment de fixer theta
C     INIT.. =1 indique que la variable a ete proprement initialisee (dans un
C       fichier suite portant les valeurs adaptees)
C
        IVSEXT(ISCAL) = -999
        THETVS(ISCAL) = -999.D0
        INITVS(ISCAL) = 0
C
      ENDDO

C
C
C --- Schema convectif
C       (a decider, centre-upwind si ordre 2,  test de pente a decider)
C
      DO II = 1, NVARMX
        BLENCV(II) = -999.D0
      ENDDO
      DO II = 1, NVARMX
        ISCHCV(II) = 1
        ISSTPC(II) = -999
      ENDDO
C
C --- Reconstruction des gradients
C       On donne les valeurs par defaut
C       Pour la methode de limitation, on decidera plus tard
C         selon les choix de l'utilisateur
C       On n'active pas l'extrapolation des gradients par defaut
C         meme pour la pression (par securite : moins stable sur certains cas)
C
      IMRGRA = 0
      ANOMAX = -GRAND*10.D0
C
      DO II = 1, NVARMX
        NSWRGR(II) = 100
        NSWRSM(II) = -999
        IMLIGR(II) = -999
        IRCFLU(II) = 1
      ENDDO
C
      DO II = 1, NVARMX
        EPSRGR(II) = 1.D-5
        CLIMGR(II) = 1.5D0
        EXTRAG(II) = 0.D0
      ENDDO
C
C --- Solveurs iteratifs
C       La methode de resolution sera choisie selon les equations
C       La valeur de epsilon relatif est tres faible
C         (1.D-5 pourrait suffire)
C       On met IDIRCL a 1 pour toutes les variables. Pour toutes les
C       variables sauf pression et fb (en v2f) on a ISTAT=1, le
C       decalage de diagonale ne sera pas active. Pour la pression
C       on decalera la diagonale si necessaire. Pour fb, on sait qu'il
C       y a un autre terme diagonal (meme si ISTAT=0), donc IDIRCL
C       sera mis a 0 dans varpos.
C
      DO II = 1, NVARMX
        NITMAX(II) = 10000
        IRESOL(II) = -1
        IDIRCL(II) = 1
        NDIRCL(II) = 0
      ENDDO
      DO II = 1, NVARMX
        EPSILO(II) = -999.D0
        EPSRSM(II) = -999.D0
      ENDDO
C
C --- Multigrille
C       On donne ici les valeurs par defaut des options de base.
C
      DO II = 1, NVARMX
        IMGR(II) = 0
      ENDDO
C
      DO II = 1, NVARMX
        NCYMAX(II) = 100
        NITMGF(II) = 10
        NAGMX0(II)= 8
        IAGMX0(II)= 1
        NCPMGR(II)= 0
      ENDDO
C
C --- Nombre max de cellules et de niveaux
C
      NCEGRM = 30
      NGRMAX = 25
C
      RLXP1 = 1.D0
C
C --- Suite de calcul
C       Calcul non suite par defaut
C       Ecriture du fichier suite auxiliaire par defaut
C       Lecture du fichier suite auxiliaire par defaut (si suite de calcul)
C       La correspondance nouveaux scalaires -> anciens scalaires
C         sera etablie plus tard (usini1 et lecamo)
C       L'indicateur de suite du module thermique 1D de paroi est initialise
C         par defaut a -1, pour obliger l'utilisateur a le remplir dans
C         uspt1d.
C       Idem pour l'indicateur de suite de la methode des vortex.
C
      ISUITE = 0
      IECAUX = 1
      ILEAUX = 1
      DO II = 1, NSCAMX
        ISCOLD(II) = -999
      ENDDO
      ISUIT1 = -1
      ISUIVO = -1
C
C --- Reperage du temps
C
      NTPABS = 0
      NTCABS = NTPABS
      NTMABS = 10
C
      INPDT0 = 0
C
      TTPABS = 0.D0
      TTCABS = TTPABS
C
C --- Marche en temps
C       Par defaut pas de temps uniforme et constant,
C         sans coef multiplicatif
C       Dans les cas a pas de temps variable, par defaut
C         COUMAX = FOUMAX = 0.5 et variation max 10%
C       Les autres grandeurs sont a renseigner par l'utilisateur
C       Pour DTMIN et DTMAX, si on impose DTMIN > DTMAX,
C         les bornes sont prises egales a +/-1D12
C
      IDTVAR = 0
C
      COUMAX = 1.D0
      FOUMAX = 10.D0
      DTMIN  = -GRAND*10.D0
      DTMAX  = -GRAND*10.D0
      VARRDT = 0.1D0
C
      DTREF  = -GRAND*10.D0
C
      DO II = 1, NVARMX
        CDTVAR(II) = 1.D0
        RELAXV(II) =-999.D0
      ENDDO
      RELXST = 0.9D0
C
C
C     Par defaut, pas de limitation du pas de temps liee aux
C     effets de densite
C
      IPTLRO = 0
C
C --- Turbulence
C     Le modele de turbulence devra etre choisi par l'utilisateur
C     En fait on n'a pas besoin d'initialiser ITYTUR (cf. varpos)
C     On suppose qu'il n'y a pas de temperature
C
      DO IPHAS = 1, NPHSMX
        ITURB (IPHAS) =-999
        ITYTUR(IPHAS) =-999
        ISCALT(IPHAS) =-1
C Parfois, IGRHOK=1 donne des vecteurs non physiques en paroi
C        IGRHOK(IPHAS) = 1
        IGRHOK(IPHAS) = 0
        IGRAKE(IPHAS) = 1
        IDEUCH(IPHAS) =-999
        ILOGPO(IPHAS) = 1
        ICLKEP(IPHAS) = 0
        IKECOU(IPHAS) =-999
        IRIJNU(IPHAS) = 0
        IRIJRB(IPHAS) = 0
        IRIJEC(IPHAS) = 0
        IGRARI(IPHAS) = 1
        IDIFRE(IPHAS) = 1
        ICLSYR(IPHAS) = 0
        ICLPTR(IPHAS) = 0
        IDRIES(IPHAS) =-1
      ENDDO
C
C
C --- Viscosite secondaire
C
      DO IPHAS = 1, NPHSMX
        IVISSE(IPHAS) = 1
      ENDDO
C
C --- Stokes
C     On suppose que l'on traite l'etape de Stokes
C
      IPRCO  = 1
      DO IPHAS = 1, NPHSMX
        IREVMC(IPHAS) = 0
        ARAK  (IPHAS) = 1.D0
      ENDDO
C     indicateur developpeur temporaire sur le test de conv resolp
      IRNPNW = 1
C
C --- Couplage U-P
C       Non active par defaut
C
      IPUCOU = 0
C
C --- Prise en compte de l'equilibre gradient de pression
C       termes sources de gravite et de pertes de charge
C       Non active par defaut
C       ICALHY=1 permet de calculer la pression
C       hydrostatique pour les Dirichlet de pression en sortie
C       Sera modifie dans modini
C
      IPHYDR = 0
      ICALHY = -1
C
C --- Champ de vitesse fige (non fige par defaut)
C
      ICCVFG = 0
C
C --- Interpolation face des viscosites
C     = 0 ARITHMETIQUE
C     = 1 HARMONIQUE
C
      IMVISF = 0
C
C --- Type des CL, tables de tri
C       Sera calcule apres usclim.
C
      DO II = 1, NTYPMX
        DO IPHAS = 1, NPHSMX
          IDEBTY(II,IPHAS) = 0
          IFINTY(II,IPHAS) = 0
        ENDDO
      ENDDO
C
C --- Traitement de la temperature pour couplage SYRTHES
C
C     TRAITEMENT PRECIS DE LA TEMPERATURE AU BORD, VOIR CONDLI
C     (UTILISE POUR COUPLAGE SYRTHES et le RAYONNEMENT)
C     Disons 0 pour eviter que des temperatures issues du traitement
C       des flux aux parois ne viennent polluer le champ de T.
C
      ITBRRB = 0
      DO ISCAL = 1, NSCAMX
        ICPSYR(ISCAL) = -999
      ENDDO
C
C
C --- Estimateurs d'erreur pour Navier-Stokes
C       En attendant un retour d'experience et pour l'avoir,
C       on active les estimateurs par defaut.
C
C     Le numero d'estimateur IEST prend les valeurs suivantes
C        IESPRE : prediction
C                 L'estimateur est base sur la grandeur
C                 I = rho_n (u*-u_n)/dt + rho_n u_n grad u*
C                   - rho_n div (mu+mu_t)_n grad u* + grad P_n
C                   - reste du smb(u_n, P_n, autres variables_n)
C                 Idealement nul quand les methodes de reconstruction
C                   sont parfaites et le systeme est resolu exactement
C        IESDER : derive
C                 L'estimateur est base sur la grandeur
C                 I = div (flux de masse corrige apres etape pression)
C                 Idealement nul quand l'equation de Poisson est resolue
C                   exactement
C        IESCOR : correction
C                 L'estimateur est base sur la grandeur
C                 I = div (rho_n u_(n+1))
C                 Idealement nul quand IESDER est nul et que le passage
C                   des flux de masse aux faces vers les vitesses au centre
C                   se fait dans un espace de fonctions a divergence nulle.
C        IESTOT : total
C                 L'estimateur est base sur la grandeur
C                 I = rho_n (u_(n+1)-u_n)/dt + rho_n u_(n+1) grad u_(n+1)
C                   - rho_n div (mu+mu_t)_n grad u_(n+1) + gradP_(n_+1)
C                   - reste du smb(u_(n+1), P_(n+1), autres variables_n)
C                 Le flux du terme convectif est calcule a partir de u_(n+1)
C                   pris au centre des cellules (et non pas a partir du flux
C                   de masse aux faces actualise)
C
C     On evalue l'estimateur IEST selon les valeurs de IESCAL
C
C        IESCAL(IEST,IPHAS) = 0 : l'estimateur IEST n'est pas calcule
C        IESCAL(IEST,IPHAS) = 1 : l'estimateur IEST   est     calcule,
C                         sans contribution du volume  (on prend abs(I))
C        IESCAL(IEST,IPHAS) = 2 : l'estimateur IEST   est     calcule,
C                         avec contribution du volume ("norme L2")
C                         soit abs(I)*SQRT(Volume_cellule),
C                         sauf pour IESCOR : on calcule abs(I)*Volume_cellule
C                         pour mesurer l'ecart en kg/s
C
C
      DO IPHAS = 1, NPHSMX
        DO IEST = 1, NESTMX
          IESCAL(IEST,IPHAS) = 0
        ENDDO
      ENDDO
C
C --- Somme de NCEPDC (pour les calculs paralleles)
C
      DO IPHAS = 1, NPHSMX
        NCPDCT(IPHAS) = 0
      ENDDO
C
C --- Somme de NCETSM (pour les calculs paralleles)
C
      DO IPHAS = 1, NPHSMX
        NCTSMT(IPHAS) = 0
      ENDDO
C
C --- Somme de NFPT1D (pour les calculs paralleles)
C
      NFPT1T = 0
C
C --- Calcul des moyennes temporelles
C
C
      DO IMOM = 1, NBMOMX
C       Pas de temps de depart (-1 : jamais)
        NTDMOM(IMOM) = -1
C       Ancien moment a relire ou -1 pour (re)initialisation
        IMOOLD(IMOM) = -2
      ENDDO
C
C     Variables composant les moments
      DO II = 1, NDGMOX
        DO IMOM = 1, NBMOMX
          IDFMOM(II,IMOM) = 0
        ENDDO
      ENDDO
C
C     Non utilisateur
C
C       Nombre de moments
      NBMOMT = 0
C       Nombre de tableaux ncel pourle temps cumule
      NBDTCM = 0
C
      DO IMOM = 1, NBMOMX
C       Pointeur sur les moments
        ICMOME(IMOM) = 0
C       Pointeur sur le pas de temps cumule
        IDTMOM(IMOM) = 0
C       Degre
        IDGMOM(IMOM) = 0
      ENDDO
C
C       Repere pour  diviser la variable par le temps cumule au post.
      DO II = 1, NVPPMX
        IPPMOM(II) = 0
      ENDDO
C
C
C --- Calcul de la distance a la paroi
C     Seules variables utilisateur : ICDPAR, IWARNY
C
      INEEDY = 0
      IMAJDY = 0
      ICDPAR = -999
      NITMAY = 10000
      NSWRSY = 1
      NSWRGY = 100
      IMLIGY = -999
      IRCFLY = 1
      ISCHCY = 1
      ISSTPY = 0
      IMGRPY = 0
      IWARNY = -999
      NTCMXY = 1000
C
C
      BLENCY = 0.0D0
      EPSILY = 1.0D-8
      EPSRGY = 1.0D-5
      CLIMGY = 1.5D0
      EXTRAY = 0.0D0
      COUMXY = 5000.D0
      EPSCVY = 1.0D-8
      YPLMXY = 200.D0
C
C --- Methode des vortex
      IVRTEX = 0
C
C --- Calcul des efforts aux parois
      INEEDF = 0
C
C --- Ici tout optcal.h est initialise
C
C=======================================================================
C 8. TABLEAUX DE cstphy.h
C=======================================================================
C
C --- Gravite
C
      GX = 0.D0
      GY = 0.D0
      GZ = 0.D0
C
C --- Constantes physiques de chaque phase
C       RO0,VISCL0 et CP0 devront etre initialises par l'utilisateur
C       P0 est donne par phase, mais seul P0(1) est utilise, idem
C        pour PRED0 et XYZREF
C       T0 ne sert a rien, sauf a etre dispo pour l'utilisateur
C
C       IROVAR : en attendant mieux, IROVAR indique si rho est constant
C         (et c'est donc RO0) : utilise uniquement pour les suites de
C         calcul, il indique qu'il ne faut pas relire la valeur de rho
C         qui a ete stockee dans le fichier suite. Ceci permet de ne pas
C         ecraser la valeur de rho0 (eventuellement modifiee par
C         l'utilisateur) par la valeur de l'ancien calcul.
C         Le pb ne se pose pas pour la viscosite car elle n'est relue
C         que si elle est extrapolee et elle est alors forcement variable
C         (du moins, on l'espere...) : on adopte la meme methode pour la
C         symetrie.
C
      DO IPHAS = 1, NPHSMX
        IROVAR(IPHAS) = -1
        IVIVAR(IPHAS) = -1
        RO0   (IPHAS) = -GRAND*10.D0
        VISCL0(IPHAS) = -GRAND*10.D0
        P0    (IPHAS) = 1.013D5
        PRED0 (IPHAS) = 0.D0
        XYZP0(1,IPHAS)= -RINFIN
        XYZP0(2,IPHAS)= -RINFIN
        XYZP0(3,IPHAS)= -RINFIN
        IXYZP0(IPHAS) = -1
        T0    (IPHAS) = 0.D0
        CP0   (IPHAS) = -GRAND*10.D0
      ENDDO
C
C --- Turbulence
C     YPLULI est mis a -GRAND*10. Si l'utilisateur ne l'a pas specifie dans usini1, on
C     modifie sa valeur dans modini (10.88 avec les lois de paroi invariantes,
C     1/kappa sinon)
      DO IPHAS = 1, NPHSMX
        YPLULI(IPHAS) = -GRAND*10.D0
      ENDDO
      XKAPPA  = 0.42D0
      CSTLOG  = 5.2D0
C
      APOW    = 8.3D0
      BPOW    = 1.D0/7.D0
      CPOW    = APOW**(2.D0/(1.D0-BPOW))
      DPOW    = 1.D0/(1.D0+BPOW)
C
      CMU     = 0.09D0
      CMU025  = CMU**0.25D0
C
C   pour le k-epsilon
      CE1     = 1.44D0
      CE2     = 1.92D0
      CE4     = 1.20D0
      SIGMAK  = 1.00D0
      SIGMAE  = 1.30D0
C
C   pour le Rij-epsilon standard (et SSG pour CRIJ3)
      CRIJ1  = 1.80D0
      CRIJ2  = 0.60D0
      CRIJ3  = 0.55D0
      CRIJEP = 0.18D0
      CSRIJ  = 0.22D0
      CRIJP1 = 0.50D0
      CRIJP2 = 0.30D0
C
C   pour le Rij-epsilon SSG
      CSSGS1  = 1.70D0
      CSSGS2  =-1.05D0
      CSSGR1  = 0.90D0
      CSSGR2  = 0.80D0
      CSSGR3  = 0.65D0
      CSSGR4  = 0.625D0
      CSSGR5  = 0.20D0
      CSSGE2  = 1.83D0
C
C   pour la LES
      DO IPHAS = 1, NPHSMX
        XLESFL(IPHAS) = 2.D0
        ALES(IPHAS)   = 1.D0
        BLES(IPHAS)   = 1.D0/3.D0
        CSMAGO(IPHAS) = 0.065D0
        XLESFD(IPHAS) = 1.5D0
        SMAGMX(IPHAS) = 10.D0*CSMAGO(IPHAS)
        CDRIES(IPHAS) = 26.D0
      ENDDO
C
C   pour le v2f phi-model
      CV2FA1 = 0.05D0
      CV2FE2 = 1.85D0
      CV2FMU = 0.22D0
      CV2FC1 = 1.4D0
      CV2FC2 = 0.3D0
      CV2FCT = 6.D0
      CV2FCL = 0.25D0
      CV2FET = 110.D0
C
C   pour le modele k-omega sst
      CKWSK1 = 1.D0/0.85D0
      CKWSK2 = 2.D0
      CKWSW1 = 2.D0
      CKWSW2 = 1.D0/0.856D0
      CKWBT1 = 0.075D0
      CKWBT2 = 0.0828D0
      CKWGM1 = CKWBT1/CMU - XKAPPA**2/(CKWSW1*SQRT(CMU))
      CKWGM2 = CKWBT2/CMU - XKAPPA**2/(CKWSW2*SQRT(CMU))
      CKWA1  = 0.31D0
      CKWC1  = 10.D0
C
C   echelle de longueur negative, recalculee par la suite
C    ou entree par l'utilisateur
      DO IPHAS = 1, NPHSMX
        ALMAX(IPHAS)   = -GRAND*10.D0
      ENDDO
C
C   vitesse de reference pour l'initialisation de la turbulence
C    doit etre entree par l'utilisateur, sauf s'il initialise lui-meme
C    la turbulence.
      DO IPHAS = 1, NPHSMX
        UREF(IPHAS)    = -GRAND*10.D0
      ENDDO
C
C   longueur caracteristique pour le modele de longueur de melange
C    doit etre entree par l'utilisateur
      DO IPHAS = 1, NPHSMX
        XLOMLG(IPHAS)    = -GRAND*10.D0
      ENDDO
C
C --- Scalaires
C       L'utilisateur devra remplir VISLS0
C       On remplira plus tard, selon les modifs utilisateur,
C         ISCSTH, IPHSCA
C       On modifiera eventuellement plus tard, selon modifs utilisateur,
C         IVISLS, ce dernier a ete initialise plus haut
C       On donne la valeur par defaut pour les autres
C       En particulier, on suppose qu'on n'a pas de variance (ISCAVR=0)
C         qu'on clippe les variances a zero seulement,
C         qu'on ne clippe pas les scalaires (sauf a +/-GRAND)
C
      DO ISCAL = 1, NSCAMX
        ISCSTH(ISCAL) =-10
        ICLVFL(ISCAL) = -1
        ISCAVR(ISCAL) = 0
        IPHSCA(ISCAL) = 0
        SCAMIN(ISCAL) =-GRAND
        SCAMAX(ISCAL) =+GRAND
        VISLS0(ISCAL) =-GRAND*10.D0
        SIGMAS(ISCAL) = 1.0D0
        RVARFL(ISCAL) = 0.8D0
      ENDDO
C
C --- Ici tout cstphy a ete initialise
C
C=======================================================================
C 9. INDICATEURS DE VECTORISATION
C=======================================================================
C
C  On les prend ici egaux a -1 ; si l'utilisateur ne les positionne pas
C    a zero, ils seront ensuite renseignes dans numvec, apres des tests
C    de vectorisation.
C
C  Sur machine vectorielle, on effectue les tests de vectorisation
C  (exemple ici pour le VPP 5000, obsolete, mais on pourrait porter
C  le code sur NEC SX8 ou Cray X1 de maniere analogue).
C  Sur machine non vectorielle, inutile de chercher a vectoriser.
C
#if defined(__uxpvp__)
      IVECTI = -1
      IVECTB = -1
#else
      IVECTI = 0
      IVECTB = 0
#endif
C
C=======================================================================
C 10. INITIALISATION DES PARAMETRES DE PERIODICITE de period.h
C=======================================================================
C
      IGUPER = 0
      IGRPER = 0
C
C=======================================================================
C 11. INITIALISATION DES PARAMETRES DE IHM de ihmpre.h
C=======================================================================
C
C     Par defaut, pas de fichier IHM consulte (on regarde ensuite si on
C       en trouve un a partir de la ligne de commande)
C
      IIHMPR = 0
C
C=======================================================================
C 12. INITIALISATION DES PARAMETRES ALE de albase.h et alstru.h
C=======================================================================
C
C --- Methode ALE
      IALE = 0
C
C --- Iterations d'initialisation fluide seul
      NALINF = -999
C
C --- Type de viscosite de maillage (isotrope par defaut)
      IORTVM = 0
C
C --- Nombre de structures
C     (sera relu ou recalcule)
      NBSTRU = -999
C
C --- Parametres du couplage implicite
      NALIMX = -999
      EPALIM = 1.D-5
C
C --- Iteration d'initialisation de l'ALE
      ITALIN = -999
C
C --- Tableaux des structures
      DO ISTR = 1, NSTRMX
        DTSTR(ISTR) = DTREF
        DO II = 1, 3
          XSTR(II,ISTR)   = 0.D0
          XPSTR(II,ISTR)  = 0.D0
          XPPSTR(II,ISTR) = 0.D0
          XSTA(II,ISTR)   = 0.D0
          XPSTA(II,ISTR)  = 0.D0
          XPPSTA(II,ISTR) = 0.D0
          XSTP(II,ISTR)   = 0.D0
          FORSTR(II,ISTR) = 0.D0
          FORSTA(II,ISTR) = 0.D0
          FORSTP(II,ISTR) = 0.D0
          XSTREQ(II,ISTR) = 0.D0
          DO JJ = 1, 3
            XMSTRU(II,JJ,ISTR) = 0.D0
            XCSTRU(II,JJ,ISTR) = 0.D0
            XKSTRU(II,JJ,ISTR) = 0.D0
          ENDDO
        ENDDO
      ENDDO
C
C --- Schema de couplage des structures
      AEXXST = -GRAND
      BEXXST = -GRAND
      CFOPRE = -GRAND
C
C --- Methode de Newmark HHT
      ALPNMK = 0.D0
      BETNMK = -GRAND
      GAMNMK = -GRAND
C
C=======================================================================
C 14. SORTIE
C=======================================================================
C
      RETURN
      END
c@z
