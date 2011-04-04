!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine iniini
!================


!===============================================================================
!  FONCTION  :
!  ---------

! INITIALISATION PAR DEFAUT DES COMMONS
!   AVANT DE PASSER LA MAIN A L'UTILISATEUR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use pointe
use albase
use alstru
use alaste
use parall
use period
use mltgrd
use ihmpre
use cplsat
use mesh

!===============================================================================

implicit none

! Local variables

integer          ii, jj, ivar, iscal, iphas, iprop, iest, imom
integer          istr

!===============================================================================

!===============================================================================
! 1. STOCKAGE DES ARGUMENTS ET IMPRESSIONS INITIALES
!===============================================================================

write(nfecra, 900)

#if defined(_CS_LANG_FR)

 900  format(/,                                                   &
'===============================================================',&
/,/,                                                        &
'                   PREPARATION DU CALCUL                     ',/,&
'                   =====================                     ',/,&
                                                                /,&
                                                                /,&
' =========================================================== ',/,&
                                                                /,&
                                                                /)

#else

 900  format(/,                                                   &
'===============================================================',&
/,/,                                                        &
'                   CALCULATION PREPARATION'                   ,/,&
'                   ======================='                   ,/,&
                                                                /,&
                                                                /,&
' ===========================================================' ,/,&
                                                                /,&
                                                                /)

#endif

!===============================================================================
! 2. ENTREES SORTIES entsor.f90
!===============================================================================

! ---> Impressions standard
!         -10000 : non initialise : voir modini
do ii = 1, nvarmx
  iwarni(ii) = -10000
enddo

! ---> NFECRA vaut 6 par defaut ou 9 en parallele (CSINIT)

! ---> Geometrie

impgeo = 10
ficgeo = 'geomet'

!    Methode des vortex : on utilise la meme unite
!                  Pour le fichier de donnees, on specifie l'unite ici, mais le
!                  nom est laisse dans usvort. On utilise 20 et pas 11
!                  car en cas d'entree multiple, les deux fichiers doivent etre
!                  ouverts en meme temps dans VORINI
impmvo = 11
impdvo = 20

! ---> Fichier stop

impstp = 12
ficstp = 'ficstp'

! ---> Fichier aval

!     NTSUIT : Periode de sauvegarde du fichier suite
!            (il est de toutes facons ecrit en fin de calcul)
!              -1   : a la fin seulement
!              0    : par defaut (4 fois par calcul)
!              > 0  : periode

ntsuit = 0

!    Methode des vortex : on utilise la meme unite
!                  et le format ascii obligatoirement
!                  (pas de detection automatique, fichiers de taille faible)
impvvo = 20

!    Fichier listing Lagrangien

implal = 80
ficlal = 'listla'
ntlal  = 1

!    Post traitement

!       ICHRVL : Post traitement du domaine fluide
!         (1 oui, 0 non)
!       ICHRBO : Post traitement du bord du domaine
!         (1 oui, 0 non)
!       ICHRSY : Post traitement des zones couplees avec SYRTHES
!         (1 oui, 0 non)
!       ICHRZE: Post traitement des zones d'echange aerorefrigerants
!         (1 oui, 0 non)

!       ICHRMD : indique si les maillages ecrits seront :
!         0 : fixes,
!         1 : deformables a topologie constante,
!         2 : modifiables (pourront etre completement redefinis en
!             cours de calcul via le sous-programme USMPST).
!        10 : comme INDMOD = 0, avec champ de deplacement
!        11 : comme INDMOD = 1, avec champ de deplacement
!        12 : comme INDMOD = 2, avec champ de deplacement

!       NTCHR  : Periode   de sortie Post
!         -1 : une seule sortie a la fin
!         >0 : periode
!       FRCHR  : frequence de sortie (en secondes)
!       ICHRVR  : Variables a sortir
!         (1 oui, 0 non, sinon : non initialise)

ichrvl = 1
ichrbo = 0
ichrsy = 0
ichrze = 0

ichrmd = 0
ntchr  = -1
frchr  = -1.d0

do ii = 1, nvppmx
  ichrvr(ii )   = -999
enddo

!       FMTCHR : format ('EnSight Gold', 'MED_fichier', ou 'CGNS')
!       OPTCHR : options associees au format de sortie

fmtchr = 'EnSight Gold'
optchr = 'binary'

! ---> Fichier thermochinie
!        FPP : utilisateur
!        JNF : Janaf
!        Les deux fichiers peuvent partager la meme unite
!          puisqu'ils sont lus l'un a pres l'autre.
!      En prime, INDJON (janaf=1 ou non=0)

impfpp = 25
ficfpp = 'dp_tch'

impjnf = impfpp
ficjnf = 'JANAF'

indjon = 1

! ---> Fichiers module atmospherique
impmet = 26
ficmet = 'meteo'

! ---> Fichiers historiques

!     EMPHIS : EMPlacement
!     PREHIS : PREfixe
!     EXTHIS : EXTension
!     IMPUSH : Unite fichiers specifiques ushist
!     FICUSH : Nom   fichiers specifiques ushist
!     IMPSTH : fichier stock + unite d'ecriture des variables
!              des structures mobiles

impsth(1) = 30
impsth(2) = 31

emphis = 'monitoring/'
prehis = 'probes_'

do ii = 1, nushmx
  impush(ii) = 32+ii
  if (irangp .le. 0) then
    write(ficush(ii),'(1a3,i3.3)')'ush',ii
  else
    write(ficush(ii),'(1a3,i3.3,1a3,i4.4)')'ush',ii,'.n_',irangp+1
  endif
enddo

! tplfmt : time plot format (1: .dat, 2: .csv, 3: both)
! ncapt  : nombre de sondes total (limite a ncaptm)
! nthist : periode de sortie (> 0 ou -1 (jamais))
! frhist : frequence de sortie, en secondes (prioritaire sur nthist si > 0)
! nthsav : periode de sauvegarde (> 0 (fichiers ouverts et refermes) ou -1 )
! ihisvr : nb de sonde et numero par variable (-999 non initialise)
! ihistr : indicateur d'ecriture des historiques des structures
!          mobiles internes (=0 ou 1)
! ncapt  : nombre de sondes total (limite a ncaptm)
! nodcap : element correspondant aux sondes
! ndrcap : rang du processus contenant nodcap (parallelisme)
! xyzcap : position demandee des sondes
! tplflw : time plot flush wall-time interval (none if <= 0)

tplfmt = 1
ncapt = 0

nthist = 1
frhist = -1.d0
nthsav = -1

do ii = 1, nvppmx
  do jj = 1, ncaptm+1
    ihisvr(ii ,jj) = -999
  enddo
enddo

ihistr = 0

do ii = 1, ncaptm
  nodcap(ii) = 1
  ndrcap(ii) = -1
enddo

do ii = 1, ncaptm
  xyzcap(1,ii) = 0.d0
  xyzcap(2,ii) = 0.d0
  xyzcap(3,ii) = 0.d0
enddo

tplflw = -1.d0

! ---> Fichiers Lagrangiens

!     IMPLA1 : Unite fichier specifique Lagrangien
!     IMPLA2 : Unite fichier specifique Lagrangien
!     IMPLA3 : Unite fichier SCRATCH pour stockage temporaire
!     IMPLA4 : Unite fichier SCRATCH pour stockage temporaire
!     IMPLA5 : Unite d'ecriture des variables associees

impla1 = 50
impla2 = 51
impla3 = 52
impla4 = 53

impla5(1)  = 54
impla5(2)  = 55
impla5(3)  = 56
impla5(4)  = 57
impla5(5)  = 58
impla5(6)  = 59
impla5(7)  = 60
impla5(8)  = 61
impla5(9)  = 62
impla5(10) = 63
impla5(11) = 64
impla5(12) = 65
impla5(13) = 66
impla5(14) = 67
impla5(15) = 68

! ---> Fichiers utilisateurs

do ii = 1, nusrmx
  impusr(ii) = 69+ii
  if (irangp .le. 0) then
    write(ficusr(ii),'(1a4,i2.2)')'usrf',ii
  else
    write(ficusr(ii),'(1a4,i2.2,1a3,i4.4)')'usrf',ii,'.n_',irangp+1
  endif
enddo

! ---> Sorties listing

!   COMMUNES
!     IPP*   : Pointeurs de reperage des variables pour les sorties
!              1 pointe sur une case poubelle et constitue donc une
!              bonne initialisation
!     NOMVAR : Nom des variables
!     ILISVR : on suit la variable (1) ou non (0) ou non initialise
!     ITRSVR : numero de variable si IPP correspond a une variable resolue (p,u,k...)
!              0 si IPP correspond a une variable annexe (cp, mut...)ou a rien
!     NTLIST : periode d'ecriture
!       ( -1 : dernier pas de temps : > 0 : periode)

do ii = 1, nvarmx
  ipprtp(ii) = 1
enddo
do ii = 1, npromx
  ipppro(ii) = 1
enddo
ippdt        = 1
ipptx        = 1
ippty        = 1
ipptz        = 1
do ii = 1, nvppmx
  ipp2ra(ii) = 0
enddo

do ii = 1, nvppmx
  nomvar(ii)    = ' '
  ilisvr(ii)    = -999
  itrsvr(ii)    = 0
enddo

ntlist = 1

!   PARAMETRES DE SUIVI DE CALCUL, MIN-MAX, CLIPMIN, CLIPMAX

do ii = 1, nvppmx
  iclpmn(ii) = 0
  iclpmx(ii) = 0
enddo

do ii = 1, nvppmx
  varmin(ii) = 0.d0
  varmax(ii) = 0.d0
  varmna(ii) = 0.d0
  varmxa(ii) = 0.d0
enddo

!   PARAMETRES DE CONVERGENCE, NORME DU SECOND MEMBRE, NOMBRE ITERATIONS
!                                RESIDU NORME, DERIVE
do ii = 1, nvppmx
  nbivar(ii) = 0
enddo
do ii = 1, nvppmx
  rnsmbr(ii) = 0.d0
  resvar(ii) = 0.d0
  dervar(ii) = 0.d0
enddo

!   PARAMETRES DU PAS DE TEMPS LOCAL

do ii = 1, 8
  do jj = 1, 4
    ptploc(ii,jj) = 0.d0
  enddo
enddo


! ---> Post traitement automatique (bord)

if (ifoenv .eq. 0) then
  ipstdv = 0
else
  ipstdv = ipstyp*ipstcl*ipstft*ipstfo
endif


! ---> CPU
!      TMARUS : marge (Arret du calcul avant limite CPU)
!        Si TMARUS negatif, le code calcule une marge seul
!        Sinon, il utilise TMARUS (donnee en secondes)

tmarus = -1.d0


! Ici entsor.f90 est completement initialise

!===============================================================================
! 3. DIMENSIONS DE dimens.f90 (GEOMETRIE, sauf NCELBR)
!===============================================================================

!---> GEOMETRIE

!---> LECTURE SELON UTILISATION PREPROCESSEUR OU ANCIEN FICHIER

   ncel   = 0
   ncelet = 0
   nfac   = 0
   nfabor = 0
   nfml   = 0
   nprfml = 0
   nnod   = 0
   lndfac = 0
   lndfbr = 0
   ncelbr = 0

   ncelgb = 0
   nfacgb = 0
   nfbrgb = 0
   nsomgb = 0

! Par defaut, on suppose qu'il n'y a pas de periodicite
   iperio = 0
   iperot = 0

if (ifoenv .eq. 0) then

!        NECESSITE FICGEO

   call ledgeo                                                    &
   !==========
        ( ndim   , ncelet , ncel   , nfac   , nfabor ,            &
          nprfml , nfml   , nnod   , lndfac , lndfbr )


! --- Pas de decoupage en parallele prevu sans preprocesseur.
!     NCELET=NCEL est fait dans ledgeo

else if (ifoenv .eq. 1) then

! --- PASSAGE DANS LEDEVI :
!        On ne lit ici que les dimensions.
!        Les tableaux seront lus apres allocation.
!        On initialise les dimensions a 0 avant la lecture, certaines
!          rubriques etant optionnelles.

   call ledevi(ndim, nfml, nprfml, iperio, iperot)
   !==========

   call tstjpe(iperio, iperot)
   !==========

endif


!===============================================================================
! 4. DIMENSIONS de dimens.f90 (PHYSIQUE)
!===============================================================================

! --- Nombre de phases, de scalaires, de scalaires a diffusivite
!            variable, de variables

nphas  = 0
nscal  = 0
nscaus = 0
nscapp = 0
nvar   = 0

! --- Nombre de proprietes physiques (utile ?)

nproce = 0
nprofa = 0
nprofb = 0
nfluma = 0

! --- Nombre de tableaux NFABOR de type COEFA (ou COEFB)

ncofab = 0


!===============================================================================
! 5. POSITION DES VARIABLES DE numvar.f90
!===============================================================================

! --- Variables de calcul resolues (RTP, RTPA)

do iphas = 1, nphsmx
  ipr   (iphas) = 0
  iu    (iphas) = 0
  iv    (iphas) = 0
  iw    (iphas) = 0
  ik    (iphas) = 0
  iep   (iphas) = 0
  ir11  (iphas) = 0
  ir22  (iphas) = 0
  ir33  (iphas) = 0
  ir12  (iphas) = 0
  ir13  (iphas) = 0
  ir23  (iphas) = 0
  iphi  (iphas) = 0
  ifb   (iphas) = 0
enddo

do iscal = 1, nscamx
  isca  (iscal) = 0
  iscapp(iscal) = 0
enddo

! --- Initialisation par defaut des commons pour la physique particuliere

call ppinii
!==========

! --- Proprietes physiques au sens large (PROPCE, PROFA, PROPFB)

do iprop  = 1, npromx
  ipproc(iprop) = 0
  ipprof(iprop) = 0
  ipprob(iprop) = 0
enddo

do iphas = 1, nphsmx
  irom  (iphas) = 0
  iviscl(iphas) = 0
  ivisct(iphas) = 0
  icour (iphas) = 0
  ifour (iphas) = 0
  icp   (iphas) = 0
  iprtot(iphas) = 0
enddo

do iscal = 1, nscamx
  ivisls(iscal) = -1
enddo

do ivar  = 1, nvarmx
  ifluma(ivar) = 0
  ifluaa(ivar) = 0
enddo

! --- Conditions aux limites (COEFA, COEFB)

icoef  = 1
icoeff = 2

do ivar = 1, nvarmx
  iclrtp(ivar,icoef ) = 0
  iclrtp(ivar,icoeff) = 0
enddo

! --- Ici tout numvar est initialise.

!===============================================================================
! 6. POSITION DES VARIABLES DE pointe.f90
!===============================================================================

! --- Geometrie

idist  = 0
idistb = 0
ipond  = 0
idijpf = 0
idiipb = 0
idofij = 0

! --- Auxiliaires independants du nb phases

icocg  = 0
icocgb = 0
itpuco = 0
idipar = 0
iyppar = 0
nfpt1d = 0
nmxt1d = 0
inppt1 = 0
iifpt1 = 0
iiclt1 = 0
ieppt1 = 0
irgpt1 = 0
itppt1 = 0
itept1 = 0
ihept1 = 0
ifept1 = 0
ixlmt1 = 0
ircpt1 = 0
idtpt1 = 0

! --- Autres auxiliaires dependant du nb de phases

iitypf = 0
iitrif = 0
iisymp = 0

do iphas = 1, nphsmx
  iifapa(iphas) = 0
  ncepdc(iphas) = 0
  iicepd(iphas) = 0
  ickupd(iphas) = 0
  ncetsm(iphas) = 0
  iicesm(iphas) = 0
  iitpsm(iphas) = 0
  ismace(iphas) = 0
  is2kw (iphas) = 0
  idvukw(iphas) = 0
enddo

! --- Auxiliaires pour la periodicite

idudxy =0
idrdxy =0
iwdudx =0
iwdrdx =0

! --- Ici tout pointe.f90 est initialise

!===============================================================================
! 7. OPTIONS DU CALCUL : TABLEAUX DE optcal.f90
!===============================================================================

! --- Definition des equations
!       (convection-diffusion instationnaire, avec dirichlet
!        sauf pour la pression, diffusion instationnaire)
!        IDIFFT en particulier multiplie la diffusion turbulente
!           (quand elle est activee par le modele)

do ii = 1, nvarmx
  iconv (ii) = 1
  istat (ii) = 1
  idiff (ii) = 1
  idifft(ii) = 1
enddo

! --- Schema en temps

!     NOTER BIEN que les valeurs de THETA pour
!       rho, visc, cp, flux de masse, termes sources
!       sont renseignees a partir des indicateurs I..EXT et ISTMPF
!       Elles ne sont pas completees directement par l'utilisateur
!       Les tests dans l'algo portent sur ces indicateurs


!   -- Schema en temps (regroupera les options suivantes)
!     = 1 Standard ordre 1
!     = 2 Standard ordre 2
!     si on veut, on peut en rajouter.
!     Cette variable conditionne toutes les autres de maniere automatique,
!     pour definir des schemas coherents dans modini.
do iphas = 1, nphsmx
  ischtp(iphas) = -999
enddo

!   -- Variables
!     Pour un schema centre (en n+1/2) il faut prendre theta = 0.5
!     On devrait pouvoir faire de l'explicite pur avec theta = 0 mais
!       ce point reste a voir
!     Ce theta sert aux termes de convection diffusion (partie implicitee
!       d'ordinaire)
!     Il est applique sous la forme (1-theta) ancien + theta nouveau
!       (ce n'est pas une extrapolation, contrairement aux termes sources)
do ii = 1, nvarmx
  thetav(ii) =-999.d0
enddo

do iphas = 1, nphsmx

!   -- Flux de masse (-999 = non initialise)
!     = 1 Standard d'ordre 1 (THETFL = -999 inutile)
!     = 0 Explicite (THETFL = 0)
!     = 2 Ordre deux (THETFL = 0.5)
  istmpf(iphas) = -999
  thetfl(iphas) =-999.d0

!   -- Termes sources Navier Stokes
!     Pour les termes sources explicites en std, I..EXT definit
!       l'extrapolation -theta ancien + (1+theta) nouveau
!     = 0 explicite
!     = 1 extrapolation avec theta = 1/2
!     = 2 extrapolation avec theta = 1
!       0 implique pas de reservation de tableaux
!       1 et 2 sont deux options equivalentes, la difference etant faite
!       uniquement au moment de fixer theta
!     Pour les termes sources implicites en std, I..EXT definit
!       la mise a l'ordre 2 ou non avec le thetav de la variable associee
!     = 0 implicite (std)
!     > 0 utilisation du thetav
!     Noter cpdt que le TS d'acc. masse n'est pas regi par I..EXT
!       (il suit bilsc2)
  isno2t(iphas) = -999
  thetsn(iphas) =-999.d0
!   -- Termes sources Grandeurs turbulentes
!     Pour les termes sources explicites en std, I..EXT definit
!       l'extrapolation -theta ancien + (1+theta) nouveau
!     = 0 explicite
!     = 1 extrapolation avec theta = 1/2
!     = 2 extrapolation avec theta = 1
!       0 implique pas de reservation de tableaux
!       1 et 2 sont deux options equivalentes, la difference etant faite
!       uniquement au moment de fixer theta
!     Pour les termes sources implicites en std, I..EXT definit
!       la mise a l'ordre 2 ou non avec le thetav de la variable associee
!     = 0 implicite (std)
!     > 0 utilisation du thetav
!     Noter cpdt que le TS d'acc. masse n'est pas regi par I..EXT
!       (il suit bilsc2)
  isto2t(iphas) = -999
  thetst(iphas) =-999.d0

!    -- Proprietes physiques
!     I..EXT definit l'extrapolation -theta ancien + (1+theta) nouveau
!     = 0 explicite
!     = 1 extrapolation avec theta = 1/2
!     = 2 extrapolation avec theta = 1
!       0 implique pas de reservation de tableaux
!       1 et 2 sont deux options equivalentes, la difference etant faite
!       uniquement au moment de fixer theta
!     INIT.. =1 indique que la variable a ete proprement initialisee (dans un
!       fichier suite portant les valeurs adaptees)

!     Masse volumique
  iroext(iphas) = -999
  thetro(iphas) = -999.d0
  initro(iphas) = 0
!     Viscosite totale
  iviext(iphas) = -999
  thetvi(iphas) = -999.d0
  initvi(iphas) = 0
!     Chaleur specifique
  icpext(iphas) = -999
  thetcp(iphas) = -999.d0
  initcp(iphas) = 0

!   -- Convergence point fixe vitesse pression
  epsup (iphas) = 1.d-5

!   -- Tab de travail pour normes de navsto
  xnrmu0(iphas) = 0.d0
  xnrmu (iphas) = 0.d0

enddo

!   -- Nb d'iter point fixe vitesse pression
nterup = 1

do iscal = 1, nscamx

!   -- Termes sources des scalaires
!     Pour les termes sources explicites en std, I..EXT definit
!       l'extrapolation -theta ancien + (1+theta) nouveau
!     = 0 explicite
!     = 1 extrapolation avec theta = 1/2
!     = 2 extrapolation avec theta = 1
!       0 implique pas de reservation de tableaux
!       1 et 2 sont deux options equivalentes, la difference etant faite
!       uniquement au moment de fixer theta
!     Pour les termes sources implicites en std, I..EXT definit
!       la mise a l'ordre 2 ou non avec le thetav de la variable associee
!     = 0 implicite (std)
!     > 0 utilisation du thetav
!     Noter cpdt que le TS d'acc. masse n'est pas regi par I..EXT
!       (il suit bilsc2)
  isso2t(iscal) = -999
  thetss(iscal) =-999.d0

!    -- Proprietes physiques
!     I..EXT definit l'extrapolation -theta ancien + (1+theta) nouveau
!     = 0 explicite
!     = 1 extrapolation avec theta = 1/2
!     = 2 extrapolation avec theta = 1
!       0 implique pas de reservation de tableaux
!       1 et 2 sont deux options equivalentes, la difference etant faite
!       uniquement au moment de fixer theta
!     INIT.. =1 indique que la variable a ete proprement initialisee (dans un
!       fichier suite portant les valeurs adaptees)

  ivsext(iscal) = -999
  thetvs(iscal) = -999.d0
  initvs(iscal) = 0

enddo



! --- Schema convectif
!       (a decider, centre-upwind si ordre 2,  test de pente a decider)

do ii = 1, nvarmx
  blencv(ii) = -999.d0
enddo
do ii = 1, nvarmx
  ischcv(ii) = 1
  isstpc(ii) = -999
enddo

! --- Reconstruction des gradients
!       On donne les valeurs par defaut
!       Pour la methode de limitation, on decidera plus tard
!         selon les choix de l'utilisateur
!       On n'active pas l'extrapolation des gradients par defaut
!         meme pour la pression (par securite : moins stable sur certains cas)

imrgra = 0
anomax = -grand*10.d0

do ii = 1, nvarmx
  nswrgr(ii) = 100
  nswrsm(ii) = -999
  imligr(ii) = -999
  ircflu(ii) = 1
enddo

do ii = 1, nvarmx
  epsrgr(ii) = 1.d-5
  climgr(ii) = 1.5d0
  extrag(ii) = 0.d0
enddo

! --- Solveurs iteratifs
!       La methode de resolution sera choisie selon les equations
!       La valeur de epsilon relatif est tres faible
!         (1.D-5 pourrait suffire)
!       On met IDIRCL a 1 pour toutes les variables. Pour toutes les
!       variables sauf pression et fb (en v2f) on a ISTAT=1, le
!       decalage de diagonale ne sera pas active. Pour la pression
!       on decalera la diagonale si necessaire. Pour fb, on sait qu'il
!       y a un autre terme diagonal (meme si ISTAT=0), donc IDIRCL
!       sera mis a 0 dans varpos.

do ii = 1, nvarmx
  nitmax(ii) = 10000
  iresol(ii) = -1
  idircl(ii) = 1
  ndircl(ii) = 0
enddo
do ii = 1, nvarmx
  epsilo(ii) = -999.d0
  epsrsm(ii) = -999.d0
enddo

! --- Multigrille
!       On donne ici les valeurs par defaut des options de base.

do ii = 1, nvarmx
  imgr(ii) = 0
enddo

do ii = 1, nvarmx
  ncymax(ii) = 100
  nitmgf(ii) = 10
  nagmx0(ii)= 3
  iagmx0(ii)= 1
  ncpmgr(ii)= 0
enddo

! --- Nombre max de cellules et de niveaux

ncegrm = 30
ngrmax = 25

rlxp1 = 0.95d0

! --- Suite de calcul
!       Calcul non suite par defaut
!       Ecriture du fichier suite auxiliaire par defaut
!       Lecture du fichier suite auxiliaire par defaut (si suite de calcul)
!       La correspondance nouveaux scalaires -> anciens scalaires
!         sera etablie plus tard (usini1 et lecamo)
!       L'indicateur de suite du module thermique 1D de paroi est initialise
!         par defaut a -1, pour obliger l'utilisateur a le remplir dans
!         uspt1d.
!       Idem pour l'indicateur de suite de la methode des vortex.

isuite = 0
iecaux = 1
ileaux = 1
do ii = 1, nscamx
  iscold(ii) = -999
enddo
isuit1 = -1
isuivo = -1

! --- Reperage du temps

ntpabs = 0
ntcabs = ntpabs
ntmabs = 10

inpdt0 = 0

ttpabs = 0.d0
ttcabs = ttpabs

! --- Marche en temps
!       Par defaut pas de temps uniforme et constant,
!         sans coef multiplicatif
!       Dans les cas a pas de temps variable, par defaut
!         COUMAX = FOUMAX = 0.5 et variation max 10%
!       Les autres grandeurs sont a renseigner par l'utilisateur
!       Pour DTMIN et DTMAX, si on impose DTMIN > DTMAX,
!         les bornes sont prises egales a +/-1D12

idtvar = 0

coumax = 1.d0
foumax = 10.d0
dtmin  = -grand*10.d0
dtmax  = -grand*10.d0
varrdt = 0.1d0

dtref  = -grand*10.d0

do ii = 1, nvarmx
  cdtvar(ii) = 1.d0
  relaxv(ii) =-999.d0
enddo
relxst = 0.9d0


!     Par defaut, pas de limitation du pas de temps liee aux
!     effets de densite

iptlro = 0

! --- Turbulence
!     Le modele de turbulence devra etre choisi par l'utilisateur
!     En fait on n'a pas besoin d'initialiser ITYTUR (cf. varpos)
!     On suppose qu'il n'y a pas de temperature

do iphas = 1, nphsmx
  iturb (iphas) =-999
  itytur(iphas) =-999
  iscalt(iphas) =-1
! Parfois, IGRHOK=1 donne des vecteurs non physiques en paroi
!        IGRHOK(IPHAS) = 1
  igrhok(iphas) = 0
  igrake(iphas) = 1
  ideuch(iphas) =-999
  ilogpo(iphas) = 1
  iclkep(iphas) = 0
  ikecou(iphas) =-999
  irijnu(iphas) = 0
  irijrb(iphas) = 0
  irijec(iphas) = 0
  igrari(iphas) = 1
  idifre(iphas) = 1
  iclsyr(iphas) = 0
  iclptr(iphas) = 0
  idries(iphas) =-1
enddo


! --- Viscosite secondaire

do iphas = 1, nphsmx
  ivisse(iphas) = 1
enddo

! --- Stokes
!     On suppose que l'on traite l'etape de Stokes

iprco  = 1
do iphas = 1, nphsmx
  irevmc(iphas) = 0
  arak  (iphas) = 1.d0
enddo
!     indicateur developpeur temporaire sur le test de conv resolp
irnpnw = 1

! --- Couplage U-P
!       Non active par defaut

ipucou = 0

! --- Prise en compte de l'equilibre gradient de pression
!       termes sources de gravite et de pertes de charge
!       Non active par defaut
!       ICALHY=1 permet de calculer la pression
!       hydrostatique pour les Dirichlet de pression en sortie
!       Sera modifie dans modini

iphydr = 0
icalhy = -1

! --- Champ de vitesse fige (non fige par defaut)

iccvfg = 0

! --- Interpolation face des viscosites
!     = 0 ARITHMETIQUE
!     = 1 HARMONIQUE

imvisf = 0

! --- Type des CL, tables de tri
!       Sera calcule apres usclim.

do ii = 1, ntypmx
  do iphas = 1, nphsmx
    idebty(ii,iphas) = 0
    ifinty(ii,iphas) = 0
  enddo
enddo

! --- Traitement de la temperature pour couplage SYRTHES

!     TRAITEMENT PRECIS DE LA TEMPERATURE AU BORD, VOIR CONDLI
!     (UTILISE POUR COUPLAGE SYRTHES et le RAYONNEMENT)
!     Disons 0 pour eviter que des temperatures issues du traitement
!       des flux aux parois ne viennent polluer le champ de T.

itbrrb = 0
do iscal = 1, nscamx
  icpsyr(iscal) = -999
enddo


! --- Estimateurs d'erreur pour Navier-Stokes
!       En attendant un retour d'experience et pour l'avoir,
!       on active les estimateurs par defaut.

!     Le numero d'estimateur IEST prend les valeurs suivantes
!        IESPRE : prediction
!                 L'estimateur est base sur la grandeur
!                 I = rho_n (u*-u_n)/dt + rho_n u_n grad u*
!                   - rho_n div (mu+mu_t)_n grad u* + grad P_n
!                   - reste du smb(u_n, P_n, autres variables_n)
!                 Idealement nul quand les methodes de reconstruction
!                   sont parfaites et le systeme est resolu exactement
!        IESDER : derive
!                 L'estimateur est base sur la grandeur
!                 I = div (flux de masse corrige apres etape pression)
!                 Idealement nul quand l'equation de Poisson est resolue
!                   exactement
!        IESCOR : correction
!                 L'estimateur est base sur la grandeur
!                 I = div (rho_n u_(n+1))
!                 Idealement nul quand IESDER est nul et que le passage
!                   des flux de masse aux faces vers les vitesses au centre
!                   se fait dans un espace de fonctions a divergence nulle.
!        IESTOT : total
!                 L'estimateur est base sur la grandeur
!                 I = rho_n (u_(n+1)-u_n)/dt + rho_n u_(n+1) grad u_(n+1)
!                   - rho_n div (mu+mu_t)_n grad u_(n+1) + gradP_(n_+1)
!                   - reste du smb(u_(n+1), P_(n+1), autres variables_n)
!                 Le flux du terme convectif est calcule a partir de u_(n+1)
!                   pris au centre des cellules (et non pas a partir du flux
!                   de masse aux faces actualise)

!     On evalue l'estimateur IEST selon les valeurs de IESCAL

!        IESCAL(IEST,IPHAS) = 0 : l'estimateur IEST n'est pas calcule
!        IESCAL(IEST,IPHAS) = 1 : l'estimateur IEST   est     calcule,
!                         sans contribution du volume  (on prend abs(I))
!        IESCAL(IEST,IPHAS) = 2 : l'estimateur IEST   est     calcule,
!                         avec contribution du volume ("norme L2")
!                         soit abs(I)*SQRT(Volume_cellule),
!                         sauf pour IESCOR : on calcule abs(I)*Volume_cellule
!                         pour mesurer l'ecart en kg/s


do iphas = 1, nphsmx
  do iest = 1, nestmx
    iescal(iest,iphas) = 0
  enddo
enddo

! --- Somme de NCEPDC (pour les calculs paralleles)

do iphas = 1, nphsmx
  ncpdct(iphas) = 0
enddo

! --- Somme de NCETSM (pour les calculs paralleles)

do iphas = 1, nphsmx
  nctsmt(iphas) = 0
enddo

! --- Somme de NFPT1D (pour les calculs paralleles)

nfpt1t = 0

! --- Calcul des moyennes temporelles


do imom = 1, nbmomx
!       Pas de temps de depart (-1 : jamais)
  ntdmom(imom) = -1
!       Ancien moment a relire ou -1 pour (re)initialisation
  imoold(imom) = -2
enddo

!     Variables composant les moments
do ii = 1, ndgmox
  do imom = 1, nbmomx
    idfmom(ii,imom) = 0
  enddo
enddo

!     Non utilisateur

!       Nombre de moments
nbmomt = 0
!       Nombre de tableaux ncel pourle temps cumule
nbdtcm = 0

do imom = 1, nbmomx
!       Pointeur sur les moments
  icmome(imom) = 0
!       Pointeur sur le pas de temps cumule
  idtmom(imom) = 0
!       Degre
  idgmom(imom) = 0
enddo

!       Repere pour  diviser la variable par le temps cumule au post.
do ii = 1, nvppmx
  ippmom(ii) = 0
enddo


! --- Calcul de la distance a la paroi
!     Seules variables utilisateur : ICDPAR, IWARNY

ineedy = 0
imajdy = 0
icdpar = -999
nitmay = 10000
nswrsy = 1
nswrgy = 100
imligy = -999
ircfly = 1
ischcy = 1
isstpy = 0
imgrpy = 0
iwarny = -999
ntcmxy = 1000


blency = 0.0d0
epsily = 1.0d-8
epsrgy = 1.0d-5
climgy = 1.5d0
extray = 0.0d0
coumxy = 5000.d0
epscvy = 1.0d-8
yplmxy = 200.d0

! --- Methode des vortex
ivrtex = 0

! --- Calcul des efforts aux parois
ineedf = 0

! --- Ici tout optcal.f90 est initialise

!===============================================================================
! 8. TABLEAUX DE cstphy.f90
!===============================================================================

! --- Gravite

gx = 0.d0
gy = 0.d0
gz = 0.d0

! --- Vecteur rotation

icorio = 0

omegax = 0.d0
omegay = 0.d0
omegaz = 0.d0

! --- Constantes physiques de chaque phase
!       RO0,VISCL0 et CP0 devront etre initialises par l'utilisateur
!       P0 est donne par phase, mais seul P0(1) est utilise, idem
!        pour PRED0 et XYZREF
!       T0 ne sert a rien, sauf a etre dispo pour l'utilisateur

!       IROVAR : en attendant mieux, IROVAR indique si rho est constant
!         (et c'est donc RO0) : utilise uniquement pour les suites de
!         calcul, il indique qu'il ne faut pas relire la valeur de rho
!         qui a ete stockee dans le fichier suite. Ceci permet de ne pas
!         ecraser la valeur de rho0 (eventuellement modifiee par
!         l'utilisateur) par la valeur de l'ancien calcul.
!         Le pb ne se pose pas pour la viscosite car elle n'est relue
!         que si elle est extrapolee et elle est alors forcement variable
!         (du moins, on l'espere...) : on adopte la meme methode pour la
!         symetrie.

do iphas = 1, nphsmx
  irovar(iphas) = -1
  ivivar(iphas) = -1
  ro0   (iphas) = -grand*10.d0
  viscl0(iphas) = -grand*10.d0
  p0    (iphas) = 1.013d5
  pred0 (iphas) = 0.d0
  xyzp0(1,iphas)= -rinfin
  xyzp0(2,iphas)= -rinfin
  xyzp0(3,iphas)= -rinfin
  ixyzp0(iphas) = -1
  t0    (iphas) = 0.d0
  cp0   (iphas) = -grand*10.d0
enddo

! --- Turbulence
!     YPLULI est mis a -GRAND*10. Si l'utilisateur ne l'a pas specifie dans usini1, on
!     modifie sa valeur dans modini (10.88 avec les lois de paroi invariantes,
!     1/kappa sinon)
do iphas = 1, nphsmx
  ypluli(iphas) = -grand*10.d0
enddo
xkappa  = 0.42d0
cstlog  = 5.2d0

apow    = 8.3d0
bpow    = 1.d0/7.d0
cpow    = apow**(2.d0/(1.d0-bpow))
dpow    = 1.d0/(1.d0+bpow)

cmu     = 0.09d0
cmu025  = cmu**0.25d0

!   pour le k-epsilon
ce1     = 1.44d0
ce2     = 1.92d0
ce4     = 1.20d0
sigmak  = 1.00d0
sigmae  = 1.30d0

!   pour le Rij-epsilon standard (et SSG pour CRIJ3)
crij1  = 1.80d0
crij2  = 0.60d0
crij3  = 0.55d0
crijep = 0.18d0
csrij  = 0.22d0
crijp1 = 0.50d0
crijp2 = 0.30d0

!   pour le Rij-epsilon SSG
cssgs1  = 1.70d0
cssgs2  =-1.05d0
cssgr1  = 0.90d0
cssgr2  = 0.80d0
cssgr3  = 0.65d0
cssgr4  = 0.625d0
cssgr5  = 0.20d0
cssge2  = 1.83d0

!   pour la LES
do iphas = 1, nphsmx
  xlesfl(iphas) = 2.d0
  ales(iphas)   = 1.d0
  bles(iphas)   = 1.d0/3.d0
  csmago(iphas) = 0.065d0
  cwale(iphas)  = 0.25d0
  xlesfd(iphas) = 1.5d0
  smagmx(iphas) = 10.d0*csmago(iphas)
  cdries(iphas) = 26.d0
enddo

!   pour le v2f phi-model
cv2fa1 = 0.05d0
cv2fe2 = 1.85d0
cv2fmu = 0.22d0
cv2fc1 = 1.4d0
cv2fc2 = 0.3d0
cv2fct = 6.d0
cv2fcl = 0.25d0
cv2fet = 110.d0

!   pour le modele k-omega sst
ckwsk1 = 1.d0/0.85d0
ckwsk2 = 1.d0
ckwsw1 = 2.d0
ckwsw2 = 1.d0/0.856d0
ckwbt1 = 0.075d0
ckwbt2 = 0.0828d0
ckwgm1 = ckwbt1/cmu - xkappa**2/(ckwsw1*sqrt(cmu))
ckwgm2 = ckwbt2/cmu - xkappa**2/(ckwsw2*sqrt(cmu))
ckwa1  = 0.31d0
ckwc1  = 10.d0

!   echelle de longueur negative, recalculee par la suite
!    ou entree par l'utilisateur
do iphas = 1, nphsmx
  almax(iphas)   = -grand*10.d0
enddo

!   vitesse de reference pour l'initialisation de la turbulence
!    doit etre entree par l'utilisateur, sauf s'il initialise lui-meme
!    la turbulence.
do iphas = 1, nphsmx
  uref(iphas)    = -grand*10.d0
enddo

!   longueur caracteristique pour le modele de longueur de melange
!    doit etre entree par l'utilisateur
do iphas = 1, nphsmx
  xlomlg(iphas)    = -grand*10.d0
enddo

! --- Scalaires
!       L'utilisateur devra remplir VISLS0
!       On remplira plus tard, selon les modifs utilisateur,
!         ISCSTH, IPHSCA
!       On modifiera eventuellement plus tard, selon modifs utilisateur,
!         IVISLS, ce dernier a ete initialise plus haut
!       On donne la valeur par defaut pour les autres
!       En particulier, on suppose qu'on n'a pas de variance (ISCAVR=0)
!         qu'on clippe les variances a zero seulement,
!         qu'on ne clippe pas les scalaires (sauf a +/-GRAND)

do iscal = 1, nscamx
  iscsth(iscal) =-10
  iclvfl(iscal) = -1
  iscavr(iscal) = 0
  iphsca(iscal) = 0
  scamin(iscal) =-grand
  scamax(iscal) =+grand
  visls0(iscal) =-grand*10.d0
  sigmas(iscal) = 1.0d0
  rvarfl(iscal) = 0.8d0
enddo

! --- Ici tout cstphy a ete initialise

!===============================================================================
! 9. INDICATEURS DE VECTORISATION
!===============================================================================

!  On les prend ici egaux a -1 ; si l'utilisateur ne les positionne pas
!    a zero, ils seront ensuite renseignes dans numvec, apres des tests
!    de vectorisation.

!  Sur machine vectorielle, on effectue les tests de vectorisation
!  (exemple ici pour le VPP 5000, obsolete, mais on pourrait porter
!  le code sur NEC SX8 ou Cray X1 de maniere analogue).
!  Sur machine non vectorielle, inutile de chercher a vectoriser.

ivecti = -1
ivectb = -1

!===============================================================================
! 10. INITIALISATION DES PARAMETRES DE PERIODICITE de period.f90
!===============================================================================

iguper = 0
igrper = 0

!===============================================================================
! 11. INITIALISATION DES PARAMETRES DE IHM de ihmpre.f90
!===============================================================================

!     Par defaut, pas de fichier IHM consulte (on regarde ensuite si on
!       en trouve un a partir de la ligne de commande)

iihmpr = 0

!===============================================================================
! 12. INITIALISATION DES PARAMETRES ALE de albase.f90 et alstru.f90
!===============================================================================

! --- Methode ALE
iale = 0

! --- Iterations d'initialisation fluide seul
nalinf = -999

! --- Type de viscosite de maillage (isotrope par defaut)
iortvm = 0

! --- Nombre de structures internes
!     (sera relu ou recalcule)
nbstru = -999

! --- Nombre de structures externes
!     (sera relu ou recalcule)
nbaste = -999

! --- Numero d'iteration de couplage externe
ntcast = 0

! --- Parametres du couplage implicite
nalimx = -999
epalim = 1.d-5

! --- Iteration d'initialisation de l'ALE
italin = -999

! --- Tableaux des structures
do istr = 1, nstrmx
  dtstr(istr) = dtref
  do ii = 1, 3
    xstr(ii,istr)   = 0.d0
    xpstr(ii,istr)  = 0.d0
    xppstr(ii,istr) = 0.d0
    xsta(ii,istr)   = 0.d0
    xpsta(ii,istr)  = 0.d0
    xppsta(ii,istr) = 0.d0
    xstp(ii,istr)   = 0.d0
    forstr(ii,istr) = 0.d0
    forsta(ii,istr) = 0.d0
    forstp(ii,istr) = 0.d0
    xstreq(ii,istr) = 0.d0
    do jj = 1, 3
      xmstru(ii,jj,istr) = 0.d0
      xcstru(ii,jj,istr) = 0.d0
      xkstru(ii,jj,istr) = 0.d0
    enddo
  enddo
enddo

! --- Schema de couplage des structures
aexxst = -grand
bexxst = -grand
cfopre = -grand

! --- Methode de Newmark HHT
alpnmk = 0.d0
betnmk = -grand
gamnmk = -grand

!===============================================================================
! 14. INITIALISATION DES PARAMETRES DE COUPLAGE CS/CS
!===============================================================================

! --- Nombre de couplage
nbrcpl = 0

! --- Couplage uniquement par les faces
ifaccp = 0

!===============================================================================
! 15. SORTIE
!===============================================================================

return
end subroutine
