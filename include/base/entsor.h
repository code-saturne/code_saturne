!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

!                              entsor.h
!===============================================================================

! --- SORTIE STD
!     NFECRA : UNITE SORTIE STD
!     IWARNI : NIVEAU D'IMPRESSION

integer           nfecra,iwarni(nvarmx)
common / icontr / nfecra,iwarni

! --- FICHIER GEOMETRIE

character*6       ficgeo
common / ageotl / ficgeo

integer                   impgeo
common / igeotl /         impgeo

! --- FICHIER SUITE AMONT
!     -------------------
!     FICAMO        --> fichier suite de base
!     FICAMX        --> fichier suite auxiliaire
!     FICMT1        --> fichier suite module thermique 1D
!     FICMVO,IMPMVO --> fichier suite methode des vortex (ASCII
!                       obligatoirement, structure specifique)
!     IMPDVO        --> fichier de donnees de la methode des vortex
!                       (nom FICDAT laisse a l'utilisateur dans usvort)

character*13      ficamo,ficamx,ficmt1,ficmvo
common / afamon / ficamo,ficamx,ficmt1,ficmvo

integer                   impmvo, impdvo
common / ifamon /         impmvo, impdvo

! --- FICHIERS SUITE AVAL
!     -------------------
!     FICAVA        --> fichier suite de base
!     FICAVX        --> fichier suite auxiliaire
!     FICVT1        --> fichier suite module thermique 1D
!     FICVVO,IMPVVO --> fichier suite methode des vortex (ASCII
!                       obligatoirement, structure specifique)


character*13      ficava,ficavx,ficvt1,ficvvo
common / afaval / ficava,ficavx,ficvt1,ficvvo

integer           ntsuit, impvvo
common / ifaval / ntsuit, impvvo

! --- FICHIER SUITE AMONT RAYONNEMENT

character*13      ficamr
common / afamrd / ficamr

! --- FICHIER SUITE AVAL RAYONNEMENT

character*13      ficavr
common / afavrd / ficavr

integer                   ntsuir
common / ifavrd /         ntsuir


! --- FICHIER NSTOP

character*6       ficstp
common / afarre / ficstp

integer                   impstp
common / ifarre /         impstp

! --- LECTURE PREPROCESSEUR

integer           ifoenv
common / ifenvp / ifoenv

! --- SORTIES POST TRAITEMENT (via FVM)

!     ICHRVL : Post traitement du domaine fluide
!     ICHRBO : Post traitement du bord du domaine
!     ICHRSY : Post traitement des zones de bord couplees avec SYRTHES
!     ICHRZE : Post traitement des zones d'echange aerorefrigerants
!     ICHRMD : Indique si les maillages ecrits seront :
!               0 : fixes,
!               1 : deformables a topologie constante,
!               2 : modifiables (pourront etre completement redefinis
!                   en cours de calcul via le sous-programme USMPST),
!              10 : comme INDMOD = 0, avec champ de deplacement,
!              11 : comme INDMOD = 1, avec champ de deplacement,
!              12 : comme INDMOD = 2, avec champ de deplacement.
!     NTCHR  : frequence de sortie par defaut ( > 0 ou -1 (a la fin) )
!     ICHRVR : on sort la variable (1) ou non (0) ou non initialise
!     FMTCHR : format de sortie ('EnSight Gold', 'MED_fichier', 'CGNS')
!     OPTCHR : options associees au format de sortie

integer           ichrvl, ichrbo, ichrsy, ichrze,                 &
                  ichrmd, ntchr , ichrvr(nvppmx)
common / iepost / ichrvl, ichrbo, ichrsy, ichrze,                 &
                  ichrmd, ntchr, ichrvr

character*32      fmtchr
character*96      optchr
common / aepost / fmtchr, optchr

! --- FICHIER THERMOPHYSIQUE SPECIFIQUE PHYSIQUE PARTICULIERE

!     IMP    --> Unite logique du fichier
!     FPP    --> Fichier utilisateur lorsqu'on utilise Janaf
!     JNF    --> Janaf
!     JON    --> Utilisation de Janaf ou non

character*6       ficfpp
character*5       ficjnf
common / afcppp / ficfpp, ficjnf

integer           impfpp, impjnf, indjon
common / ifcppp / impfpp, impjnf, indjon

! --- INPUT FILES FOR THE ATMOSPHERIC SPECIFIC PHYSICS
!        IMPMET --> logical unit of the meteo profile file
!        FICMET --> name of the meteo profile file

  integer impmet
  common / ifcmet / impmet

  character*10 ficmet
  common / afcmet / ficmet

! --- FICHIERS HISTORIQUES

!     IMPHIS : fichier stock + unite d'ecriture des variables
!     EMPHIS : EMPlacement
!     EXTHIS : EXTension des fichiers
!     IMPUSH : Unite fichiers specifiques ushist
!     FICUSH : Nom   fichiers specifiques ushist
!     IMPSTH : fichier stock + unite d'ecriture des variables
!              des structures mobiles

character*80      emphis, exthis
common / avhist / emphis, exthis

character*13      ficush(nushmx)
common / afhist / ficush

integer                   imphis(2), impush(nushmx), impsth(2)
common / ifhist /         imphis, impush, impsth

!     NCAPT  : nombre de sondes total (limite a NCAPTM)
!     NTHIST : Frequence de sortie
!         ( > 0 ou -1 (jamais) ou non initialise -999)
!     NTHSAV : Frequence de sauvegarde
!         ( > 0 ou -1 (a la fin) ou non initialise -999)
!     IHISVR : nb de sonde et numero par variable
!         (-999 non initialise)
!     IHISTR : indicateur d'ecriture des historiques des structures
!              mobiles internes (=0 ou 1)
!     NCAPT  : nombre de sondes total (limite a NCAPTM)
!     NODCAP : element correspondant aux sondes
!     NDRCAP : rang processus contenant NODCAP (parallelisme)
!     XYZCAP : position demandee des sondes

integer           ncapt, nthist, nthsav,                          &
                  ihisvr(nvppmx,ncaptm+1), ihistr,                &
                  nodcap(ncaptm), ndrcap(ncaptm)

common / ivhist / ncapt , nthist, nthsav,                         &
                  ihisvr, ihistr,                                 &
                  nodcap, ndrcap

double precision  xyzcap(3,ncaptm)
common / rvhist / xyzcap


! --- FICHIERS LAGRANGIENS

!   - FICHIER SUITE ET SUITE STATISTISQUE  AMONT LAGRANGIEN

character*13      ficaml, ficmls
common / afamla / ficaml, ficmls

!   - FICHIER SUITE ET SUITE STATISTISQUE AVAL LAGRANGIEN

character*13      ficavl, ficvls
common / afavla / ficavl, ficvls


!   - FICHIER LISTING LAGRANGIEN

!     FICLAL : Nom du fichier
!     IMPLAL : Unite du fichier
!     NTLAL  : Periode de sortie du listing
character*6       ficlal
common / afalal / ficlal

integer                   implal, ntlal
common / ifalal /         implal, ntlal

!   - FICHIER HISTORIQUE LAGRANGIEN

integer                   impli1, impli2
common / ifalah /         impli1, impli2

!   - AUTRES FICHIERS LAGRANGIEN

integer           impla1 , impla2 , impla3 , impla4 , impla5(15)
common / ifalag / impla1 , impla2 , impla3 , impla4 , impla5


! --- FICHIERS UTILISATEURS

character*13      ficusr(nusrmx)
common / afuser / ficusr

integer           impusr(nusrmx)
common / ifuser / impusr


! --- SORTIES LISTING

!   COMMUNES
!     IPP*   : Pointeurs de reperage des variables pour les sorties
!     NOMVAR : Nom des variables
!     ILISVR : on suit la variable (1) ou non (0) ou non initialise
!     ITRSVR : numero de variable si IPP correspond a une variable resolue (p,u,k...)
!              0 si IPP correspond a une variable annexe (cp, mut...)ou a rien
!     NTLIST : periode  ecriture
!       ( -1 : dernier pas de temps : > 0 : periode)

integer           ipprtp(nvarmx),                                 &
                  ipppro(npromx),                                 &
                  ippdt         ,                                 &
                  ipptx         , ippty         , ipptz  ,        &
                  ipp2ra(nvppmx)
common / ipntpp / ipprtp        ,                                 &
                  ipppro        ,                                 &
                  ippdt         ,                                 &
                  ipptx         , ippty         , ipptz  ,        &
                  ipp2ra

character*80      nomvar(nvppmx)
common / anampp / nomvar

integer           ilisvr(nvppmx),itrsvr(nvppmx)
common / ipostp / ilisvr        ,itrsvr


integer           ntlist
common / ilisti / ntlist

!   PARAMETRES DE SUIVI DE CALCUL, MIN-MAX, CLIPMIN, CLIPMAX

integer           iclpmn(nvppmx) , iclpmx(nvppmx)
common / isuivi / iclpmn         , iclpmx

double precision  varmin(nvppmx) , varmax(nvppmx),                &
                  varmna(nvppmx) , varmxa(nvppmx)
common / rsuivi / varmin ,         varmax,                        &
                  varmna ,         varmxa


!   PARAMETRES DE CONVERGENCE, NORME DU SECOND MEMBRE, NOMBRE ITERATIONS
!                                RESIDU NORME, DERIVE

integer           nbivar(nvppmx)
common / iconvg / nbivar

double precision  rnsmbr(nvppmx) ,                                &
                  resvar(nvppmx) , dervar(nvppmx)
common / rconvg / rnsmbr ,                                        &
                  resvar         , dervar


!   PARAMETRES DU PAS DE TEMPS LOCAL =
!     NB COURANT, FOURIER ET COMBINE MIN ET MAX + POINTS ASSOCIES
!     PTPLOC(.,1)   = NOMBRE , PTPLOC(.,2 3 et 4) = POINT ASSOCIE
!     PTPLOC(1 2,.) = COURANT MIN/MAX
!     PTPLOC(3 4,.) = FOURIER MIN/MAX
!     PTPLOC(5 6,.) = COU/FOU MIN/MAX
!     PTPLOC(7 8,.) = DT      MIN/MAX
!     NCLPTR        = NB DE CLIPPINGS PAR LE PDT MAX LIE AUX EFFETS DE DENSITE
!     RPDTRO        = RAPPORT MAX ENTRE DT ET DTmax LIE AUX EFFETS DE DENSITE
!                      +LOCALISATION

integer           nclptr
common / iptloc /  nclptr

double precision  ptploc(8,4), rpdtro(4)
common / rptloc /  ptploc, rpdtro


!   PARAMETRES DES SORTIES AU BORD =

!     IPSTDV = PROPRIETES POST TRAITEES
!        IPSTDV EST LE PRODUIT DES VALEURS ENTIERES SUIVANTES (NOMBRES PREMIERS) :
!          IPSTYP => YPLUS
!          IPSTCL => VARIABLES NON RECONSTRUITES (suffisant pour Dirichlet)
!          IPSTFT => FLUX THERMIQUE RECONSTRUIT

!     SI IPSTDV = 1 = PAS DE SORTIE


integer            ipstdv
common / iipstd /  ipstdv


! --- TEMPS CPU

double precision  tmarus
common / temcpu / tmarus

! FIN
