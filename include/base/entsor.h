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
C                              entsor.h
C***********************************************************************
C
C --- SORTIE STD
C     NFECRA : UNITE SORTIE STD
C     IWARNI : NIVEAU D'IMPRESSION
C
      INTEGER           NFECRA,IWARNI(NVARMX)
      COMMON / ICONTR / NFECRA,IWARNI
C
C --- FICHIER GEOMETRIE
C
      CHARACTER*6       FICGEO
      COMMON / AGEOTL / FICGEO
C
      INTEGER                   IMPGEO
      COMMON / IGEOTL /         IMPGEO
C
C --- FICHIER SUITE AMONT
C     -------------------
C     FICAMO        --> fichier suite de base
C     FICAMX        --> fichier suite auxiliaire
C     FICMT1        --> fichier suite module thermique 1D
C     FICMVO,IMPMVO --> fichier suite methode des vortex (ASCII
C                       obligatoirement, structure specifique)
C     IMPDVO        --> fichier de donnees de la methode des vortex
C                       (nom FICDAT laisse a l'utilisateur dans usvort)
C
      CHARACTER*13      FICAMO,FICAMX,FICMT1,FICMVO
      COMMON / AFAMON / FICAMO,FICAMX,FICMT1,FICMVO
C
      INTEGER                   IMPMVO, IMPDVO
      COMMON / IFAMON /         IMPMVO, IMPDVO
C
C --- FICHIERS SUITE AVAL
C     -------------------
C     FICAVA        --> fichier suite de base
C     FICAVX        --> fichier suite auxiliaire
C     FICVT1        --> fichier suite module thermique 1D
C     FICVVO,IMPVVO --> fichier suite methode des vortex (ASCII
C                       obligatoirement, structure specifique)
C
C
      CHARACTER*13      FICAVA,FICAVX,FICVT1,FICVVO
      COMMON / AFAVAL / FICAVA,FICAVX,FICVT1,FICVVO
C
      INTEGER           NTSUIT, IMPVVO
      COMMON / IFAVAL / NTSUIT, IMPVVO
C
C --- FICHIER SUITE AMONT RAYONNEMENT
C
      CHARACTER*13      FICAMR
      COMMON / AFAMRD / FICAMR
C
C --- FICHIER SUITE AVAL RAYONNEMENT
C
      CHARACTER*13      FICAVR
      COMMON / AFAVRD / FICAVR
C
      INTEGER                   NTSUIR
      COMMON / IFAVRD /         NTSUIR
C
C
C --- FICHIER NSTOP
C
      CHARACTER*6       FICSTP
      COMMON / AFARRE / FICSTP
C
      INTEGER                   IMPSTP
      COMMON / IFARRE /         IMPSTP
C
C --- LECTURE PREPROCESSEUR
C
      INTEGER           IFOENV
      COMMON / IFENVP / IFOENV
C
C --- SORTIES POST TRAITEMENT (via FVM)
C
C     ICHRVL : Post traitement du domaine fluide
C     ICHRBO : Post traitement du bord du domaine
C     ICHRSY : Post traitement des zones de bord couplees avec Syrthes
C     ICHRMD : Indique si les maillages ecrits seront :
C               0 : fixes,
C               1 : deformables a topologie constante,
C               2 : modifiables (pourront etre completement redefinis
C                   en cours de calcul via le sous-programme USMPST),
C              10 : comme INDMOD = 0, avec champ de déplacement,
C              11 : comme INDMOD = 1, avec champ de déplacement,
C              12 : comme INDMOD = 2, avec champ de déplacement.
C     NTCHR  : frequence de sortie par defaut ( > 0 ou -1 (a la fin) )
C     ICHRVR : on sort la variable (1) ou non (0) ou non initialise
C     FMTCHR : format de sortie ('EnSight Gold', 'MED_fichier', 'CGNS')
C     OPTCHR : options associees au format de sortie
C
      INTEGER           ICHRVL, ICHRBO, ICHRSY,
     &                  ICHRMD, NTCHR , ICHRVR(NVPPMX)
      COMMON / IEPOST / ICHRVL, ICHRBO, ICHRSY,
     &                  ICHRMD, NTCHR, ICHRVR
C
      CHARACTER*32      FMTCHR
      CHARACTER*96      OPTCHR
      COMMON / AEPOST / FMTCHR, OPTCHR
C
C --- FICHIER SPECIFIQUE PHYSIQUE PARTICULIERE
C
C     IMP    --> Unite logique du fichier
C     FPP    --> Fichier utilisateur lorsqu'on utilise Janaf
C     JNF    --> Janaf
C     JON    --> Utilisation de Janaf ou non
C
      CHARACTER*6       FICFPP
      CHARACTER*5       FICJNF
      COMMON / AFCPPP / FICFPP, FICJNF
C
      INTEGER           IMPFPP, IMPJNF, INDJON
      COMMON / IFCPPP / IMPFPP, IMPJNF, INDJON
C
C --- FICHIERS HISTORIQUES
C
C     IMPHIS : fichier stock + unite d'ecriture des variables
C     EMPHIS : EMPlacement
C     EXTHIS : EXTension des fichiers
C     IMPUSH : Unite fichiers specifiques ushist
C     FICUSH : Nom   fichiers specifiques ushist
C     IMPSTH : fichier stock + unite d'ecriture des variables
C              des structures mobiles
C
      CHARACTER*80      EMPHIS, EXTHIS
      COMMON / AVHIST / EMPHIS, EXTHIS
C
      CHARACTER*13      FICUSH(NUSHMX)
      COMMON / AFHIST / FICUSH
C
      INTEGER                   IMPHIS(2), IMPUSH(NUSHMX), IMPSTH(2)
      COMMON / IFHIST /         IMPHIS, IMPUSH, IMPSTH
C
C     NCAPT  : nombre de sondes total (limite a NCAPTM)
C     NTHIST : Frequence de sortie
C         ( > 0 ou -1 (jamais) ou non initialise -999)
C     NTHSAV : Frequence de sauvegarde
C         ( > 0 ou -1 (a la fin) ou non initialise -999)
C     IHISVR : nb de sonde et numero par variable
C         (-999 non initialise)
C     IHISTR : indicateur d'ecriture des historiques des structures
C              mobiles internes (=0 ou 1)
C     NCAPT  : nombre de sondes total (limite a NCAPTM)
C     NODCAP : element correspondant aux sondes
C     NDRCAP : rang processus contenant NODCAP (parallelisme)
C     XYZCAP : position demandee des sondes
C
      INTEGER           NCAPT, NTHIST, NTHSAV,
     &                  IHISVR(NVPPMX,NCAPTM+1), IHISTR,
     &                  NODCAP(NCAPTM), NDRCAP(NCAPTM)
C
      COMMON / IVHIST / NCAPT , NTHIST, NTHSAV,
     &                  IHISVR, IHISTR,
     &                  NODCAP, NDRCAP
C
      DOUBLE PRECISION  XYZCAP(3,NCAPTM)
      COMMON / RVHIST / XYZCAP
C
C
C --- FICHIERS LAGRANGIENS
C
C   - FICHIER SUITE ET SUITE STATISTISQUE  AMONT LAGRANGIEN
C
      CHARACTER*13      FICAML, FICMLS
      COMMON / AFAMLA / FICAML, FICMLS
C
C   - FICHIER SUITE ET SUITE STATISTISQUE AVAL LAGRANGIEN
C
      CHARACTER*13      FICAVL, FICVLS
      COMMON / AFAVLA / FICAVL, FICVLS
C
C
C   - FICHIER LISTING LAGRANGIEN
C
C     FICLAL : Nom du fichier
C     IMPLAL : Unite du fichier
C     NTLAL  : Periode de sortie du listing
      CHARACTER*6       FICLAL
      COMMON / AFALAL / FICLAL
C
      INTEGER                   IMPLAL, NTLAL
      COMMON / IFALAL /         IMPLAL, NTLAL
C
C   - FICHIER HISTORIQUE LAGRANGIEN
C
      INTEGER                   IMPLI1, IMPLI2
      COMMON / IFALAH /         IMPLI1, IMPLI2
C
C   - AUTRES FICHIERS LAGRANGIEN
C
      INTEGER           IMPLA1 , IMPLA2 , IMPLA3 , IMPLA4 , IMPLA5(15)
      COMMON / IFALAG / IMPLA1 , IMPLA2 , IMPLA3 , IMPLA4 , IMPLA5
C
C
C --- FICHIERS UTILISATEURS
C
      CHARACTER*13      FICUSR(NUSRMX)
      COMMON / AFUSER / FICUSR
C
      INTEGER           IMPUSR(NUSRMX)
      COMMON / IFUSER / IMPUSR
C
C
C --- SORTIES LISTING
C
C   COMMUNES
C     IPP*   : Pointeurs de reperage des variables pour les sorties
C     NOMVAR : Nom des variables
C     ILISVR : on suit la variable (1) ou non (0) ou non initialise
C     ITRSVR : numero de variable si IPP correspond a une variable resolue (p,u,k...)
C              0 si IPP correspond a une variable annexe (cp, mut...)ou a rien
C     NTLIST : periode  ecriture
C       ( -1 : dernier pas de temps : > 0 : periode)
C
      INTEGER           IPPRTP(NVARMX),
     &                  IPPPRO(NPROMX),
     &                  IPPDT         ,
     &                  IPPTX         , IPPTY         , IPPTZ  ,
     &                  IPP2RA(NVPPMX)
      COMMON / IPNTPP / IPPRTP        ,
     &                  IPPPRO        ,
     &                  IPPDT         ,
     &                  IPPTX         , IPPTY         , IPPTZ  ,
     &                  IPP2RA
C
      CHARACTER*80      NOMVAR(NVPPMX)
      COMMON / ANAMPP / NOMVAR
C
      INTEGER           ILISVR(NVPPMX),ITRSVR(NVPPMX)
      COMMON / IPOSTP / ILISVR        ,ITRSVR
C

      INTEGER           NTLIST
      COMMON / ILISTI / NTLIST
C
C   PARAMETRES DE SUIVI DE CALCUL, MIN-MAX, CLIPMIN, CLIPMAX
C
      INTEGER           ICLPMN(NVPPMX) , ICLPMX(NVPPMX)
      COMMON / ISUIVI / ICLPMN         , ICLPMX
C
      DOUBLE PRECISION  VARMIN(NVPPMX) , VARMAX(NVPPMX),
     &                  VARMNA(NVPPMX) , VARMXA(NVPPMX)
      COMMON / RSUIVI / VARMIN ,         VARMAX,
     &                  VARMNA ,         VARMXA
C
C
C   PARAMETRES DE CONVERGENCE, NORME DU SECOND MEMBRE, NOMBRE ITERATIONS
C                                RESIDU NORME, DERIVE
C
      INTEGER           NBIVAR(NVPPMX)
      COMMON / ICONVG / NBIVAR
C
      DOUBLE PRECISION  RNSMBR(NVPPMX) ,
     &                  RESVAR(NVPPMX) , DERVAR(NVPPMX)
      COMMON / RCONVG / RNSMBR ,
     &                  RESVAR         , DERVAR
C
C
C   PARAMETRES DU PAS DE TEMPS LOCAL =
C     NB COURANT, FOURIER ET COMBINE MIN ET MAX + POINTS ASSOCIES
C     PTPLOC(.,1)   = NOMBRE , PTPLOC(.,2 3 et 4) = POINT ASSOCIE
C     PTPLOC(1 2,.) = COURANT MIN/MAX
C     PTPLOC(3 4,.) = FOURIER MIN/MAX
C     PTPLOC(5 6,.) = COU/FOU MIN/MAX
C     PTPLOC(7 8,.) = DT      MIN/MAX
C     NCLPTR        = NB DE CLIPPINGS PAR LE PDT MAX LIE AUX EFFETS DE DENSITE
C     RPDTRO        = RAPPORT MAX ENTRE DT ET DTmax LIE AUX EFFETS DE DENSITE
C                      +LOCALISATION
C
      INTEGER           NCLPTR
      COMMON / IPTLOC /  NCLPTR

      DOUBLE PRECISION  PTPLOC(8,4), RPDTRO(4)
      COMMON / RPTLOC /  PTPLOC, RPDTRO
C
C
C   PARAMETRES DES SORTIES AU BORD =
C
C     IPSTDV = PROPRIETES POST TRAITEES
C        IPSTDV EST LE PRODUIT DES VALEURS ENTIERES SUIVANTES (NOMBRES PREMIERS) :
C          IPSTYP => YPLUS
C          IPSTCL => VARIABLES NON RECONSTRUITES (suffisant pour Dirichlet)
C          IPSTFT => FLUX THERMIQUE RECONSTRUIT
C
C     SI IPSTDV = 1 = PAS DE SORTIE
C
C
      INTEGER            IPSTDV
      COMMON / IIPSTD /  IPSTDV
C
C
C --- TEMPS CPU
C
      DOUBLE PRECISION  TMARUS
      COMMON / TEMCPU / TMARUS
C
C FIN
c@z
