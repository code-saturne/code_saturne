!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

subroutine mttssc &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iscal  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   iconra ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   crvexp , crvimp ,                                              &
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , w10    , w11    ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------


! CALCUL DES TERMES SOURCES POUR LE SCALAIRE ISCAL

!    SOUS-PROGRAMME SPECIFIQUE A MATISSE (COPIE DE USTSSC)

!    LES SEULS SCALAIRES TRAITES SONT
!      ISCA(ITAAMT) TEMPERATURE DE L'AIR
!      ISCA(ITPCMT) TEMPERATURE DE PEAU DES COLIS
!      ISCA(ITPPMT) TEMPERATURE DE PEAU DES MURS


! ON RESOUT RHO*VOLUME*D(VAR)/DT = CRVIMP*VAR + CRVEXP

! ON FOURNIT ICI CRVIMP ET CRVEXP (ILS CONTIENNENT RHO*VOLUME)
!    CRVEXP en kg variable/s :
!     ex : pour les temperatures      kg degres/s
!          pour les enthalpies        Joules/s
!    CRVIMP en kg /s :

! VEILLER A UTILISER UN CRVIMP NEGATIF
! (ON IMPLICITERA CRVIMP
!  IE SUR LA DIAGONALE DE LA MATRICE, LE CODE AJOUTERA :
!   MAX(-CRVIMP,0) EN SCHEMA STANDARD EN TEMPS
!       -CRVIMP    SI LES TERMES SOURCES SONT A L'ORDRE 2

! CES TABLEAUX SONT INITIALISES A ZERO AVANT APPEL A CE SOUS
!   PROGRAMME ET AJOUTES ENSUITE AUX TABLEAUX PRIS EN COMPTE
!   POUR LA RESOLUTION

! EN CAS D'ORDRE 2 DEMANDE SUR LES TERMES SOURCES, ON DOIT
!   FOURNIR CRVEXP A L'INSTANT N     (IL SERA EXTRAPOLE) ET
!           CRVIMP A L'INSTANT N+1/2 (IL EST  DANS LA MATRICE,
!                                     ON LE SUPPOSE NEGATIF)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iscal            ! i  ! <-- ! scalar number                                  !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! iconra           ! te ! <-- ! tab de connectivite pour                       !
! (ncelet+1)       !    !     !  le rayonnement et les panaches                !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! crvexp(ncelet    ! tr ! --> ! tableau de travail pour part explicit          !
! crvimp(ncelet    ! tr ! --> ! tableau de travail pour terme instat           !
! viscf(nfac)      ! tr ! --- ! tableau de travail    faces internes           !
! viscb(nfabor     ! tr ! --- ! tableau de travail    faces de bord            !
! xam(nfac,2)      ! tr ! --- ! tableau de travail    faces de bord            !
! w1..11(ncelet    ! tr ! --- ! tableau de travail    cellules                 !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "cstnum.f90"
include "paramx.f90"
include "pointe.f90"
include "numvar.f90"
include "entsor.f90"
include "optcal.f90"
include "cstphy.f90"
include "parall.f90"
include "period.f90"
include "matiss.f90"


!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse
integer          iscal

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          iconra(ncelet)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision xam(nfac,2)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet), w11(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          ivar  , iiscvr, ipcrom, iphas
integer          icoul , ifml  , iel   , jel
integer          modntl
integer          nzones, izone
double precision xnu0  , xlamb0, prand0, phi0
double precision gmat  , untier
double precision fform , emis
double precision hplus , hmoin , hmax  , difx  , dify  , dhmur
double precision fragh , xxav  , tbulk , un0   , ptot
double precision raleiz, reynol
double precision xnsraz, xnsrey, xnsalv
double precision hraz  , hrey  , h0    , hray
double precision hmrey , hmraz , hm0   , hmray , hnb
double precision cpcont, cpparo
double precision xlimin, xlimax, yrgmin, yrgmax, zalmin, zalmax
double precision sqz
double precision dhe   , dhs   , ptcn
double precision fptsto, pvccsc, pccsc

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

! --- Gestion memoire
idebia = idbia0
idebra = idbra0


! --- Numero du scalaire a traiter : ISCAL (argument)

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
chaine = nomvar(ipprtp(ivar))

! --- Indicateur de variance
!         Si ISCAVR = 0 :
!           le scalaire ISCAL n'est pas une variance
!         Si ISCAVR > 0 et ISCAVR < NSCAL + 1 :
!           le scalaire ISCAL est une variance associee
!           au scalaire ISCAVR
iiscvr = iscavr(iscal)

! --- Numero de phase associee au scalaire ISCAL
iphas = iphsca(iscal)

! --- Masse volumique
ipcrom = ipproc(irom(iphas))


! --- Reperage du pas de temps pour les impressions listing
if(ntlist.gt.0) then
  modntl = mod(ntcabs,ntlist)
elseif(ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
  modntl = 0
else
  modntl = 1
endif

! --- Impression
if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif
if((irangp.le.0).and.(modntl.eq.0)) then
  write(nfecra,2001) chaine(1:8)
endif

! --- Numerique
untier = 1.d0/3.d0

!===============================================================================
! 1. DIMENSIONS GEOMETRIQUES
!===============================================================================


!===============================================================================
! 2. PROPRIETES PHYSIQUES
!===============================================================================

! --- Proprietes de l'air

!     Nombre de Prandtl moleculaire
prand0 = 0.7d0
!     Viscosite cinematique moleculaire (nu en m2.s-1)
xnu0   = xmumat/ro0(iphas)
!     Condutivite thermique moleculaire (lambda en W.m-1.K-1)
xlamb0 = (xnu0/prand0)*cp0(iphas)*ro0(iphas)


! --- Capacite thermique fictive pour le solide

!     On cherche un etat stationnaire comme limite d'un transitoire
!       Les equations portant sur la temperature T des solides
!         comportent un terme en CP dT/dt qui tend vers zero a
!         convergence. La valeur de CP est utilisee pour accelerer
!         la convergence : elle correspond a l'utilisation d'un
!         pas de temps adapte pour la temperature des solides (comme
!         on pourrait le faire avec CDTVAR)

cpcont = 1.d0
cpparo = 1.d0


!===============================================================================
! 2. EROSION DES PANACHES DE CONVECTION NATURELLE
!===============================================================================

!     Le but est de reporter dans la zone de ciel une partie de la
!       puissance transmise a l'air dans la zone de stockage.
!       Ceci permet de simuler l'effet des panaches, le debit
!       enthalpique de ceux-ci etant prescrit par l'utilisateur.

!     Calcul de DHS, DHE

! --- Debit enthalpique de convection naturelle calcule (initialisation)
dhs = 0.d0

! --- Debit enthalpique de convection naturelle impose pour les panaches

!     Nul par defaut (pas de modelisation des panaches)
dhe = 0.d0
!     Calcul, si modelisation des panaches
if(imdcnt.eq.1) then
  dhe = dhpcnt/frdtra
endif

!===============================================================================
! 2. RAYONNEMENT
!===============================================================================

! --- Flux d'un conteneur (puissance ramenee a la surface en W.m-2)

!       En cas de stockage en alveole, il y a un plenum inferieur libre
if(ialveo.eq.1) then
  phi0 = puicon/((hreso-hplen)*dmcont*pi)
else
  phi0 = puicon/(hreso*dmcont*pi)
endif

! --- Facteur de forme (alveole ou mur) / conteneur

!     En alveole, le conteneur voit l'alveole entiere et c'est tout
if(ialveo.eq.1) then

  fform = 1.d0

!     Sans alveole, le facteur de forme depend de la structure du reseau
else

!       Pour un reseau en ligne, considerer une rangee suffit
!         (valable jusqu'a PTRRES/DMCONT=4)
  fform =  sqrt( sqrt(dmcont/ptrres) + (dmcont/ptrres) )

!       Pour un reseau en quinconce (reseau a pas triangulaire), on
!         complete par la contribution d'une seconde rangee
  if(iconlg.eq.0) then
    fform = max( (fform + 0.22d0), 1.d0 )
  endif

endif


!  --- Emissivite equivalente (conteneur / mur)
emis = 1.d0/(1.d0/emicon + 1.d0/emimur)

!  --- Densite de surface d'echange (par unite de hauteur)
xxav = pi*dmcont/(ptrres*plgres)


!  --- Connectivite pour le rayonnement colis / plafond
!                et pour la modelisation des panaches
!        Une connectivite est necessaire pour relier la couche
!          superieure de la zone occupee par des colis et les
!          mailles situees directement sous le plafond

!     Cote du haut des mailles du haut des conteneurs
hplus = hreso
!     Cote du bas des mailles du haut des conteneurs
hmoin = hreso - epchel
!     Cote du bas des mailles touchant le plafond
hmax  = epchel*nchest - epchel

!     Calcul de la connectivite
!       Attention, c'est quadratique ...
!       mais on ne le fait qu'au premier passage

!     Calcul uniquement si pas deja fait et si utile
!       (utile si les colis ne touchent pas le plafond)
if( (icnrok.eq.0).and.(hreso.lt.epchel*nchest) ) then
!     On repere les mailles de plafond dans la zone de stockage
  do iel = 1, ncel
    ifml  = ifmcel(iel   )
    icoul = iprfml(ifml,1)
    if(icoul.eq.icmtst.and.xyzcen(3,iel).gt.hmax)then
!     On repere la maille du haut des colis situee dessous
      do jel = 1, ncel
        difx = abs(xyzcen(1,jel)-xyzcen(1,iel))
        dify = abs(xyzcen(2,jel)-xyzcen(2,iel))
        if(difx.lt.ptrres*1.d-2.and.dify.lt.plgres*1.d-2) then
          if ( xyzcen(3,jel).lt.hplus.and.                        &
               xyzcen(3,jel).gt.hmoin ) then
            iconra(iel) = jel
          endif
        endif
      enddo
    endif
  enddo
!     On indique que la connectivite a ete calculee
  icnrok = 1
endif

!===============================================================================

!  ON VA COMPLETER CI-DESSOUS LES TERMES SOURCES
!     - PUREMENT EXPLICITES (PUISSANCE)
!     - PARTIELLEMENT EXPLICITES ET IMPLICITES (ECHANGES)


!  ON UTILISE LE MODELE DE TERME SOURCE SUIVANT, POUR UN SCALAIRE F :

!                             S = A * F + B

!     APPARAISSANT DANS LES EQUATIONS DE TRANSPORT DE F SOUS LA FORME :

!                       RHO VOLUME D(F)/Dt = VOLUME*S


!   CE TERME A UNE PARTIE QU'ON VEUT IMPLICITER         : A
!           ET UNE PARTIE QU'ON VA TRAITER EN EXPLICITE : B


!   PAR EXEMPLE SI ON A :
!     A = - RHO / TAUF
!     B =   RHO * PRODF
!        AVEC
!     TAUF   = 10.D0  [secondes  ] (TEMPS DE DISSIPATION DE F)
!     PRODF  = 100.D0 [variable/s] (PRODUCTION DE F PAR UNITE DE TEMPS)

!   ON A ALORS
!     CRVIMP(IEL) = VOLUME(IEL)* A = - VOLUME(IEL) (RHO / TAUF )
!     CRVEXP(IEL) = VOLUME(IEL)* B =   VOLUME(IEL) (RHO * PRODF)

!===============================================================================

! === Initialisation
! ===========================================================

! --- Compteurs
!       pour moyenne et affichage, mais pas seulement

!     Coefficient d'echange moyen de convection forcee
hmrey = 0.d0
!     Coefficient d'echange moyen de convection naturelle
hmraz = 0.d0
!     Coefficient d'echange convectif moyen efficace
hm0   = 0.d0
!     Coefficient d'echange radiatif moyen pour les parois
hmray = 0.d0
!     Compteur pour moyenne
hnb   = 0.d0
!     Puissance totale sans panache
ptot  = 0.d0
!     Puissance totale avec panache
ptcn  = 0.d0

!     Gravite
gmat = sqrt(gx**2+gy**2+gz**2)


! === Boucle sur les cellules et selection selon couleur
! ===========================================================

do iel = 1, ncel

!     Vitesse horizontale
  un0 = sqrt(rtp(iel,iu(iphas))**2+rtp(iel,iv(iphas))**2)

!     Temperature de reference Kelvin
  tbulk = (tinit+rtp(iel,isca(itaamt)))*0.5d0 + tkelvi

!     Raleigh/DeltaT = FRAGH = g*H**3*beta/(nu a)
!       . g    = GMAT
!       . H    = XYZCEN(3,IEL) : on suppose que le sol est a z=0
!       . beta = 1/TBULK
!       . nu   = XNU0
!       . a    = nu/Pr = XNU0/PRAND0

  fragh = gmat*xyzcen(3,iel)**3*prand0/(tbulk*xnu0**2)

!     Couleur pour selection de la zone de stockage
  ifml  = ifmcel(iel   )
  icoul = iprfml(ifml,1)

  if(icoul.eq.icmtst) then


! === Calcul du coefficient d'echange par rayonnement
! ===========================================================

!     Coeff d'echange = flux / (S (Tcolis-Tparoi))
!       = Sigma Epsilon FactForme (Tcolis+Tparoi)(Tcolis**2+Tparoi**2)

!     La paroi consideree est le plafond ou les murs ou les alveoles

! --- Zone de ciel
!       On utilise la connectivite calculee precedemment
    if(icnrok.eq.1.and.xyzcen(3,iel).gt.hmax) then

!     En presence d'alveoles,
!       le plafond echange avec les alveoles face a face
!       (seule l'alveole directement en face est vue,
!        l'alveole voit un disque de plafond)
!       La temperature du plafond et des alveoles est representee
!       par le meme scalaire
!       Ce coefficient est utilise pour le calcul des termes sources
!       permettant de determiner la temperature de peau des parois
!       (plafond, mur, alveoles)
!       On le verra aussi intervenir pour un calcul de temperature
!       de colis, mais dans la zone de ciel ou la temperature de colis
!       n'a pas de signification (et comme il n'y a pas de diffusion
!       en alveoles, c'est ok)
      if(ialveo.eq.1) then
        hray = stephn*emimur                                      &
             * (0.25d0*pi*dmcont**2)/(ptrres*plgres)              &
             * ( rtp(iconra(iel),isca(itppmt))+tkelvi             &
                +rtp(       iel ,isca(itppmt))+tkelvi)            &
             * ((rtp(iconra(iel),isca(itppmt))+tkelvi)**2         &
               +(rtp(       iel ,isca(itppmt))+tkelvi)**2)

!     Sans alveoles,
!       le plafond echange avec les colis face a face
!       Ce coefficient sera utilise pour le calcul des termes sources
!       permettant de determiner la temperature de peau du plafond
      else
        hray = stephn*emis                                        &
             * (0.25d0*pi*dmcont**2)/(ptrres*plgres)              &
             * ( rtp(iconra(iel),isca(itpcmt))+tkelvi             &
                +rtp(       iel ,isca(itppmt))+tkelvi)            &
             * ((rtp(iconra(iel),isca(itpcmt))+tkelvi)**2         &
               +(rtp(       iel ,isca(itppmt))+tkelvi)**2)
      endif

! --- Zones de murs ou d'alveoles
!       les murs et les alveoles echangent avec les colis
!       le facteur de forme calcule precedemment s'applique
!       il prend en compte le cas avec ou sans alveoles
!       Ce coefficient sera utilise pour le calcul des termes sources
!       permettant de determiner
!         - la temperature de peau des parois
!         - la temperature de peau des colis, en alveole (FFORM=1)

    else
      hray = stephn*emis                                          &
           * fform                                                &
           * ( rtp(iel,isca(itpcmt))+tkelvi                       &
              +rtp(iel,isca(itppmt))+tkelvi)                      &
           * ((rtp(iel,isca(itpcmt))+tkelvi)**2                   &
             +(rtp(iel,isca(itppmt))+tkelvi)**2)
    endif


! === Puissance transmise a l'air et a la peau des colis par conduction
! ===========================================================

!     Pour l'air :
!       en stationnaire, toute la puissance des colis se retrouve tot ou
!       tard dans l'air (elle peut transiter par la peau des colis,
!       par les alveoles, par les murs, mais elle est extraite par
!       l'air).

!     Pour la peau des colis :
!       toute la puissance des colis est transmise par conduction a
!       la peau des colis.


    if (iscal.eq.itaamt.or.iscal.eq.itpcmt) then

! --- Puissance des colis par unite de surface au sol
!       On prend en compte la carte de remplissage ensuite

!     On doit fournit un CRVEXP en Watt/CP, soit de la dimension de
!       PUICON/CP.

!     Comme la puissance integree sur la hauteur du colis doit etre
!       la puissance totale du colis, PUICON, on doit imposer dans
!       chaque maille une fraction de PUICON.
!     Pour une repartition uniforme en hauteur, cette fraction est
!       la hauteur de la maille VOLUME/(PTRRES*PLGRES) divisee par
!       la hauteur totale du colis.
!     Les repartitions non uniformes en hauteur sont traitees par
!       les cartes de repartition en altitude, plus bas.
!     On impose donc CRVEXP = PUICON/CP * VOLUME/(PTRRES*PLGRES)
!       puis on divise par la hauteur de colis (plus bas, lors du
!       traitement des cartes de repartition en altitude).

!     Noter que comme on multiplie par VOLUME/(PTRRES*PLGRES) ici
!       et qu'on divise plus bas par EPCHEL, on pourrait simplifier,
!       puisque, a priori, EPCHEL = VOLUME/(PTRRES*PLGRES)


!     Air ambiant
      if (iscal.eq.itaamt) then
        crvexp(iel) = puicon/cp0(iphas)                           &
             * volume(iel)/(ptrres*plgres)
!     Peau des colis
      elseif(iscal.eq.itpcmt) then
        crvexp(iel) = puicon/cpcont                               &
             * volume(iel)/(ptrres*plgres)
      endif

! --- Prise en compte des cartes de remplissage (adimensionnelles)
!     Les entrepots ne sont pas uniformement remplis.

!     Pour représenter la carte d'encombrement des lignes de colis, on
!       repere les lignes occupees (coordonnee X, indicateur ILIGNE)
!       par leurs bornes XLIMIN et XLIMAX qui sont stockees dans VIZCAR
!       (avec l'indicateur ICPUIS faisant reference a la puissance)
      nzones = nzocar(iligne,icpuis)
      do izone = 1, nzones
        xlimin = ptrres*vizcar(1,izone,iligne,icpuis)
        xlimax = ptrres*vizcar(2,izone,iligne,icpuis)
        if ((xyzcen(1,iel).ge.xlimin).and.                        &
             (xyzcen(1,iel).lt.xlimax)) then
          crvexp(iel) = crvexp(iel) * vcarth(izone,iligne)
        endif
      enddo

!     Pour représenter la carte d'encombrement des rangees de colis, on
!       repere les rangees occupees (coordonnee Y, indicateur IRANGE)
!       par leurs bornes YRGMIN et YRGMAX qui sont stockees dans VIZCAR
!       (avec l'indicateur ICPUIS faisant reference a la puissance)
      nzones = nzocar(irange,icpuis)
      do izone = 1, nzones
        yrgmin = plgres*vizcar(1,izone,irange,icpuis)
        yrgmax = plgres*vizcar(2,izone,irange,icpuis)
        if ((xyzcen(2,iel).ge.yrgmin).and.                        &
             (xyzcen(2,iel).lt.yrgmax)) then
          crvexp(iel) = crvexp(iel) * vcarth(izone,irange)
        endif
      enddo

!     Pour représenter la carte d'encombrement en hauteur de colis, on
!       repere les couches occupees (coordonnee Z, indicateur IALTIT)
!       par leurs bornes ZALMIN et ZALMAX qui sont stockees dans VIZCAR
!       (avec l'indicateur ICPUIS faisant reference a la puissance)
!     SQZ permet d'integrer VCARTH sur les zones afin que
!       de normer VCARTH de telle sorte que la multiplication par
!       VCARTH/SQZ permette de retrouver, en integrant sur la hauteur,
!       toute la puissance du colis.
      sqz = 0.d0
      nzones = nzocar(ialtit,icpuis)
      do izone = 1, nzones
        zalmin = epchel*vizcar(1,izone,ialtit,icpuis)
        zalmax = epchel*vizcar(2,izone,ialtit,icpuis)
        sqz = sqz + vcarth(izone,ialtit) *                        &
             ( vizcar(2,izone,ialtit,icpuis) -                    &
             vizcar(1,izone,ialtit,icpuis) )
      enddo
      do izone = 1, nzones
        zalmin = epchel*vizcar(1,izone,ialtit,icpuis)
        zalmax = epchel*vizcar(2,izone,ialtit,icpuis)
        if ((xyzcen(3,iel).ge.zalmin).and.                        &
             (xyzcen(3,iel).lt.zalmax)) then
          crvexp(iel) = crvexp(iel)*                              &
               vcarth(izone,ialtit)/(sqz*epchel)
        endif
      enddo

!       On calcule la puissance totale injectee (Watt)

!     Air ambiant
      if (iscal.eq.itaamt) then
        ptot = ptot + crvexp(iel)*cp0(iphas)
      endif

    endif


! === Echange convectif conteneur / air et radiatif conteneur / mur
! ===========================================================

!     Pour la determination de la temperature de peau des colis

!     On calcule le coefficient d'echange convectif efficace selon
!       la configuration (convection naturelle, forcee, mixte).
!     On ajoute le terme de rayonnement a partir du coefficient
!       equivalent HRAY calcule precedemment.

    if (iscal.eq.itpcmt) then

! --- Convection naturelle verticale le long des colis
!       Van Vliet&Ross et Vliet&Liu

!     Rayleigh a partir du flux d'un conteneur
      raleiz = fragh*phi0*xyzcen(3,iel)/xlamb0
!     Nusselt associe
      xnsraz = 0.17d0 * raleiz**0.25d0
      xnsraz = max(xnsraz,(0.6d0 * raleiz**0.2d0))
!     Coefficient d'echange de convection naturelle       conteneur / air
!       Hypothese que le zero est a la base des colis
      hraz   = xnsraz*xlamb0/xyzcen(3,iel)
!     Coefficient d'echange moyen de convection naturelle
      hmraz  = hmraz + hraz


! --- Convection forcee transversale dans le reseau de colis
!       Zukauskas

!     Reynolds (vitesse horizontale gap, diametre
      reynol = (un0*ptrres/(ptrres-dmcont)) * dmcont / xnu0
!     Nusselt associe
      xnsrey = 0.35d0*(ptrres/plgres)**0.2d0                      &
           * 0.89d0*reynol**0.6d0
!     Coefficient d'echange de convection forcee          conteneur / air
      hrey   = xnsrey*xlamb0/dmcont
!     Coefficient d'echange moyen de convection forcee
      hmrey = hmrey + hrey

! --- Compteur pour les moyennes

      hnb = hnb + 1.d0


! --- Prise en compte du type de reseau (en ligne ou en quinconce)
!       et de la presence d'alveoles

!   - Alveoles
!       convection naturelle uniquement

      if(ialveo.eq.1) then

!     Nusselt associe
        xnsalv = 0.17d0 * raleiz**0.25d0
!     Coefficient d'echange de convection naturelle       conteneur / air
!       Hypothese que le zero est a la base des colis
        h0 = xnsalv*xlamb0/xyzcen(3,iel)


!   - Reseau en ligne (pas carre) sans alveoles
!       2/3 conv.force (face+cotes) + 1/3 conv.nat (sillage arriere)

      elseif(iconlg.eq.1)then
        h0 = (hraz+2.d0*hrey)/3.d0

!   - Reseau en quinconce (pas triangulaire) sans alveoles
!       conv.force (pas de protection du colis amont)

      else
        h0 = hrey
      endif

!   - Sans alveoles
!       le coefficient d'echange est toujours au moins egal a celui
!       de la convection naturelle
      if(ialveo.ne.1) then
        h0 = max(hraz,h0)
      endif

!    - Coefficient d'echange convectif moyen efficace
      hm0 = hm0 + h0


! --- Ajout du terme conductif et du terme radiatif

!       CRVEXP contient deja la puissance des colis transmise
!         par conduction
!       On veut ajouter dans CRVEXP des Watt/CP soit H0*Surface*T/CP.
!         On a ici XXAV*VOLUME = PI*DMCONT*EPCHEL, qui est la surface
!         d'échange dans la cellule.


!   - Alveoles
!       echange par convection naturelle et
!       les colis rayonnent avec le plafond et les alveoles

      if(ialveo.eq.1) then

        crvexp(iel) = crvexp(iel)                                 &
             +h0  *(xxav*volume(iel))                             &
                  *rtp(iel,isca(itaamt))/cpcont                   &
             +hray*(xxav*volume(iel))                             &
                  *rtp(iel,isca(itppmt))/cpcont
        crvimp(iel) = -(h0+hray)*(xxav*volume(iel))/cpcont


!   - Sans alveole
!       echange par convection mixte seul
!       les colis ne rayonnent qu'entre eux, et ceci est pris en compte
!         par diffusion

      else

        crvexp(iel) = crvexp(iel)                                 &
             +h0  *(xxav*volume(iel))                             &
                  *rtp(iel,isca(itaamt))/cpcont
        crvimp(iel) = - h0      *(xxav*volume(iel))/cpcont

      endif

    endif



! === Echange convectif mur / air et radiatif mur / conteneur
! ===========================================================

!     Pour la determination de la temperature de peau des parois
!       (murs et alveoles)

    if (iscal.eq.itppmt) then

! --- Echange par convection paroi / air

!   - Alveoles
!       convection naturelle uniquement

      if(ialveo.eq.1) then

!     Rayleigh a partir du flux d'un conteneur
        raleiz = fragh*phi0*xyzcen(3,iel)/xlamb0
!     Nusselt associe
        xnsalv = 0.17d0 * raleiz**0.25d0
!     Coefficient d'echange de convection naturelle       alveole / air
!       Hypothese que le zero est a la base des colis
        hraz = xnsalv*xlamb0/xyzcen(3,iel)
!     Coefficient d'echange efficace alveole / air
!       Convection naturelle seule
        h0 = hraz


!   - Sans alveole
!       convection naturelle (Mac Adams) ou convection forcee

      else

!    .Convection naturelle

!     Rayleigh a partir de l'ecart de temperature paroi - air
        raleiz =                                                  &
             fragh*(rtp(iel,isca(itppmt))-rtp(iel,isca(itaamt)))
!     Nusselt associe (Mac Adams)
        xnsraz = 0.12d0*(abs(raleiz))**untier
!     Coefficient d'echange de convection naturelle       paroi / air
        hraz   = xnsraz*xlamb0/xyzcen(3,iel)
!     Coefficient d'echange moyen de convection naturelle
        hmraz  = hmraz + hraz


!    .Convection forcee

!     Espace mur reseau
        dhmur  = ptrres-dmcont
!     Reynolds associé (vitesse horizontale hors reseau, pres des murs)
        reynol = dhmur*un0/xnu0
!     Nusselt associe (Colburn)
        xnsrey = 0.023d0*reynol**0.8d0*prand0**untier
!     Coefficient d'echange de convection forcee          paroi / air
        hrey   = xnsrey*xlamb0/dhmur
!     Coeficient d'echange moyen de convection forcee
        hmrey  = hmrey + hrey


!     .Coefficient d'echange efficace (le plus efficace est dominant)
        h0     = max(hraz,hrey)

      endif


!   - Coefficient d'echange convectif moyen efficace
      hm0   = hm0 + h0
!   - Coefficient d'echange moyen de convection naturelle
      hmraz = hmraz + hraz
!   - Compteur pour les moyennes
      hnb   = hnb + 1.d0


! --- Ajout du terme conductif et du terme radiatif

!       On veut ajouter dans CRVEXP des Watt/CP soit H0*Surface*T/CP.
!         On calcule H0*VOLUME*T/CP et on divise par une distance
!         plus bas. On pourrait faire les choses plus directement.

!       Les echanges radiatifs sont les echanges avec les colis.


!   - Coefficient d'echange radiatif moyen pour les parois
      hmray = hmray + hray


!   - Echanges convectifs
      crvexp(iel) = h0*volume(iel)*rtp(iel,isca(itaamt))/cpparo
      crvimp(iel) =-h0*volume(iel)/cpparo


!   - Echange radiatif dans le ciel
      if(icnrok.eq.1.and.xyzcen(3,iel).gt.hmax) then

!    .Alveoles : le plafond echange avec les alveoles
        if(ialveo.eq.1) then
          crvexp(iel) = crvexp(iel)                               &
               +hray*volume(iel)                                  &
               *rtp(iconra(iel),isca(itppmt))/cpparo

!    .Sans alveole : le plafond echange avec le dessus des colis
        else
          crvexp(iel) = crvexp(iel)                               &
               +hray*volume(iel)                                  &
               *rtp(iconra(iel),isca(itpcmt))/cpparo
        endif


!   - Echange radiatif avec les colis pour les autres parois
!       (murs et alveoles)
      else
        crvexp(iel) = crvexp(iel)                                 &
             +hray*volume(iel)*rtp(iel,isca(itpcmt))/cpparo

      endif


!   - La partie implicite est la meme pour tous
      crvimp(iel) = crvimp(iel) - hray*volume(iel)/cpparo


!   - On divise par une distance
!       Plus exactement, on multiplie par la surface d'echange et
!       on divise par le volume (on aurait donc pu eviter de multiplier
!       par le volume).


!    .Au dessus ou au dessous des colis, le plafond ou le sol echangent
!       sur une surface de PTRRES*PLGRES
!       dans des cellules de volume EPCHEL*PTRRES*PLGRES.
!       Le rapport Surface/Volume est donc 1/EPCHEL

      if(xyzcen(3,iel).lt.epchel.or.xyzcen(3,iel).gt.hmax) then
        crvexp(iel) = crvexp(iel)/epchel
        crvimp(iel) = crvimp(iel)/epchel


!    .En stockage en alveoles, les alveoles echangent
!       sur une surface de PI*DHALVE*EPCHEL
!       dans des cellules de volume EPCHEL*PTRRES*PLGRES.
!       Le rapport Surface/Volume est donc PI*DHALVE/(PTRRES*PLGRES)
!     A priori DHALVE est positif, mais on laisse le ABS (pour le
!       cas, ou il serait defini comme la difference entre le
!       diametre des colis et des alveoles)

      elseif(ialveo.eq.1) then
        crvexp(iel) =                                             &
             crvexp(iel)*pi*abs(dhalve)/(ptrres*plgres)
        crvimp(iel) =                                             &
             crvimp(iel)*pi*abs(dhalve)/(ptrres*plgres)


!    .Sans alveoles, les murs echangent
!       sur une surface de PLGRES*EPCHEL
!       dans des cellules de volume EPCHEL*PTRRES*PLGRES.
!       Le rapport Surface/Volume est donc 1/PTRRES
!       La ou il n'y a pas de mur, la temperature de mur ne signifie
!       rien, mais on peut quand meme l'y calculer (il n'y a pas de
!       phenomene de diffusion, tout est local)

      else
        crvexp(iel) = crvexp(iel)/ptrres
        crvimp(iel) = crvimp(iel)/ptrres

      endif

    endif


!     Fin du test sur la zone de stockage
  endif
enddo


! === Pour l'air, sans modelisation des panaches de convection naturelle
! ===========================================================

if(iscal.eq.itaamt) then

  if(imdcnt.eq.0)then

! --- Calcul du debit enthalpique de convection naturelle
!       sortant en haut de la zone des colis : rho*W*S*Cp*DeltaT
!       (il aurait ete plus correct d'utiliser FLUMAS)

    do iel = 1, ncel

      ifml  = ifmcel(iel   )
      icoul = iprfml(ifml,1)
      if(icoul.eq.icmtst) then

        if(icnrok.eq.1.and.xyzcen(3,iel).gt.hmax) then
          jel = iconra(iel)
          dhs = dhs +                                             &
               propce(jel,ipcrom)*rtp(jel,iw(iphas))              &
               *(volume(jel)/epchel)                              &
               *cp0(iphas)*(rtp(jel,isca(itaamt))-tinit)
        endif

      endif
    enddo


! --- Calcul de la puissance en kW,
!       avec correction si representation 3D en 2D

    puitot = ptot*frdtra*1.d-3


! --- Calcul du debit enthalpique de convection naturelle en kW,
!       avec correction si representation 3D en 2D

    debcon = dhs*frdtra*1.d-3

! --- Impression de PUITOT et DEBCON
    if (irangp.le.0.and.modntl.eq.0) then
      write(nfecra,2002) puitot, debcon
    endif


! === Pour l'air, modelisation des panaches de convection naturelle
! ===========================================================

  elseif(imdcnt.eq.1) then


! --- Le debit enthalpique des panaches ne peut pas etre plus grand
!       que la puissance totale disponible ; sinon, c'est qu'il
!       a ete mal predit ou que la puissance a ete mal imposee :
!       les cartes ou la puissance d'un conteneur peuvent etre en
!       cause
!     On fait ce test ici, dans la mesure ou PTOT est calcule juste
!       au dessus.

    if(dhpcnt.ge.ptot*frdtra)then
      write(nfecra,9001) imdcnt, dhpcnt, ptot*frdtra,             &
           dhpcnt, ptot*frdtra
      call csexit (1)
      !==========
    endif

! --- Modification du terme source pour la modelisation des panaches

!       On redistribue la puissance transmise a l'air uniquement :
!         une fraction de la puissance totale est retiree de la zone
!         de stockage et injectee directement dans le ciel au dessus
!         des conteneurs. La repartition se fait a priori de maniere
!         homogene (le panache de tous les conteneurs presents recoit
!         la meme fraction de la puissance totale), puis on applique
!         la carte de repartition thermique XY.

!   - Calcul des grandeurs independantes de la maille

!     Fraction de la puissance totale a conserver
!       dans la zone de stockage
    fptsto = 1.d0-dhe/ptot

!     Puissance Volumique a ajouter dans la zone de Ciel au dessus
!       d'un Conteneur de puissance PUICON, divise par CP
    pvccsc = ( (dhe/ptot)*puicon / (hercnt*ptrres*plgres) )       &
         / cp0(iphas)

    do iel = 1, ncel

      ifml  = ifmcel(iel   )
      icoul = iprfml(ifml,1)

      if(icoul.eq.icmtst) then

!   - Reduction dans la zone de stockage

!     La puissance est reduite pour tous les conteneurs de maniere
!       homogene, en fonction du rapport entre le debit enthalpique
!       des panaches impose par l'utilisateur et la puissance totale
!       presente dans l'entrepot. Les cartes de repartition thermique
!       ayant deja ete appliquees a CRVEXP, c'est bien le rapport de
!       reduction (et non pas la valeur absolue de la reduction) qui
!       est homogene.

        crvexp(iel) = crvexp(iel)*fptsto


        if ( xyzcen(3,iel).lt.(hercnt+hreso).and.                 &
             xyzcen(3,iel).gt.hreso) then

!   - Report dans la zone des panaches (1/2)

!     Calcul de la puissance a ajouter dans la zone de ciel, au dessus
!       d'un conteneur de puissance PUICON (deja divise par CP)

          pccsc = pvccsc * volume(iel)

!     Les entrepots ne sont pas uniformement remplis.
!     Application de la carte thermique horizontale XY pour
!       selectionner les positions auxquelles un conteneur (ou
!       une fraction de conteneur est presente, puisque l'on
!       peut avoir des demi-conteneurs, par exemple sur un plan
!       de symetrie)

!     Pour représenter la carte d'encombrement des lignes de colis, on
!       repere les lignes occupees (coordonnee X, indicateur ILIGNE)
!       par leurs bornes XLIMIN et XLIMAX qui sont stockees dans VIZCAR
!       (avec l'indicateur ICPUIS faisant reference a la puissance)
          nzones = nzocar(iligne,icpuis)
          do izone = 1, nzones
            xlimin = ptrres*vizcar(1,izone,iligne,icpuis)
            xlimax = ptrres*vizcar(2,izone,iligne,icpuis)
            if ((xyzcen(1,iel).ge.xlimin).and.                    &
                 (xyzcen(1,iel).lt.xlimax)) then
              pccsc = pccsc *                                     &
                   vcarth(izone,iligne)
            endif
          enddo

!     Pour représenter la carte d'encombrement des rangees de colis, on
!       repere les rangees occupees (coordonnee Y, indicateur IRANGE)
!       par leurs bornes YRGMIN et YRGMAX qui sont stockees dans VIZCAR
!       (avec l'indicateur ICPUIS faisant reference a la puissance)
          nzones = nzocar(irange,icpuis)
          do izone = 1, nzones
            yrgmin = plgres*vizcar(1,izone,irange,icpuis)
            yrgmax = plgres*vizcar(2,izone,irange,icpuis)
            if ((xyzcen(2,iel).ge.yrgmin).and.                    &
                 (xyzcen(2,iel).lt.yrgmax)) then
              pccsc = pccsc *                                     &
                   vcarth(izone,irange)
            endif
          enddo

!     Ajout de la puissance des panaches a la puissance deja existante

          crvexp(iel) = crvexp(iel) + pccsc

!     Fin si zone de panache
        endif

!     Fin si zone de stockage
      endif


! --- Calcul de la puissance totale incluant la modelisation
!       des panaches de convection natuelle
!       (normalement la puissance totale n'a pas ete modifiee)

      ptcn = ptcn + crvexp(iel)*cp0(iphas)

!     Fin boucle NCEL
    enddo

! --- Calcul et impression de la puissance totale incluant la modelisation
!       des panaches de convection naturelle
!       en kW, avec correction si representation 3D en 2D

    puitot = ptcn*frdtra*1.d-3
    if ((irangp.le.0).and.(modntl.eq.0)) then
      write(nfecra,2003) puitot
    endif

!     Fin si on modelise les panaches
  endif
!     Fin si scalaire = air ambiant
endif



! === Impressions pour la temperature de peau des colis et des parois
! ===========================================================

! --- Valeur fictive : pas de convection forcee en alveole
!     (valeur bizarre, mais uniquement pour affichage)
if (ialveo.eq.1) then
  hmrey = -1.00001d10
endif

! --- Temperature de peau des colis

if(iscal.eq.itpcmt) then

!     Si rien a moyenner, on prend une valeur negative
!       (mais les numerateurs seront nuls, normalement)
  if (hnb.lt.epzero) then
    hnb = -epzero
  endif

!     Stockage en common pour impression ulterieure
  cfecca = hm0/hnb

!     Impressions des coefficients d'echange moyens
  if ((irangp.le.0).and.(modntl.eq.0)) then
    write(nfecra,2004) hmrey/hnb, hmraz/hnb, hm0/hnb
  endif

endif

! --- Temperature de peau des parois

if(iscal.eq.itppmt) then

!     Si rien a moyenner, on prend une valeur negative
!       (mais les numerateurs seront nuls, normalement)
  if (hnb.lt.epzero) then
    hnb = -epzero
  endif

!     Stockage en common pour impression ulterieure
  cfecma = hm0/hnb

!     Impressions des coefficients d'echange moyens
  if ((irangp.le.0).and.(modntl.eq.0)) then
    write(nfecra,2005)                                            &
         hmrey/hnb, hmraz/hnb, hmray/hnb, hm0/hnb
  endif

endif


!--------
! FORMATS
!--------

 1000 format(' TERMES SOURCES MATISSE POUR LA VARIABLE ',A8,/)

 2001 FORMAT (/,3X,'** INFORMATIONS SUR MATISSE, VARIABLE ',A8,/, &
          3X,'   -------------------------------------------')
 2002 format (                                                          &
'mati --------------------------------------------------------',/,&
'mati Puissance (sans modelisation panaches)                  ',/,&
'mati --------------------------------------                  ',/,&
'mati Puissance totale                      ', E12.5  ,'    kW',/,&
'mati Debit enthalpique de conv. naturelle  ', E12.5  ,'    kW',/,&
'mati --------------------------------------------------------',/)
 2003 format (                                                          &
'mati --------------------------------------------------------',/,&
'mati Puissance (avec modelisation panaches)                  ',/,&
'mati --------------------------------------                  ',/,&
'mati Puissance totale avec erosion panaches', E12.5  ,'    kW',/,&
'mati --------------------------------------------------------',/)
 2004 format (                                                          &
'mati --------------------------------------------------------',/,&
'mati Coeff. d echange moyens (calcul T° de peau des colis)   ',/,&
'mati -----------------------------------------------------   ',/,&
'mati Convection forcee                   ', E12.5  ,' W/m2/°C',/,&
'mati Convection naturelle                ', E12.5  ,' W/m2/°C',/,&
'mati Efficace                            ', E12.5  ,' W/m2/°C',/,&
'mati --------------------------------------------------------',/)
 2005 format (                                                          &
'mati --------------------------------------------------------',/,&
'mati Coeff. d echange moyens (calcul T° de peau des parois)  ',/,&
'mati ------------------------------------------------------  ',/,&
'mati Convection forcee                   ', E12.5  ,' W/m2/°C',/,&
'mati Convection naturelle                ', E12.5  ,' W/m2/°C',/,&
'mati Rayonnement                         ', E12.5  ,' W/m2/°C',/,&
'mati Efficace                            ', E12.5  ,' W/m2/°C',/,&
'mati --------------------------------------------------------',/)

 9001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET MATISSE (MTTSSC)                      ',/,&
'@    =========                                               ',/,&
'@      DEBIT ENTHALPIQUE DES PANACHES INADAPTE.              ',/,&
'@                                                            ',/,&
'@    La modelisation des panaches est activee                ',/,&
'@      avec IMDCNT = ',I10   ,'                              ',/,&
'@    Le debit enthalpique des panaches de convection         ',/,&
'@      naturelle doit etre strictement inferieur a la        ',/,&
'@      puissance totale des colis presents dans l''entrepot  ',/,&
'@      puisqu''il en represente une partie.                  ',/,&
'@                                                            ',/,&
'@    Le debit enthalpique prescrit est  DHPCNT  =',E12.5      ,/,&
'@    La puissance totale reelle est PTOT*FRDTRA =',E12.5      ,/,&
'@                                                            ',/,&
'@    L''inegalite suivante n''est pas verifiee :             ',/,&
'@          DHPCNT    <    PTOT*FRDTRA                        ',/,&
'@       ',E12.5   ,'     ',E12.5   ,'                        ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Desactiver la modelisation des panaches ou              ',/,&
'@      modifier la valeur du debit enthalpique prescrit ou   ',/,&
'@      modifier la puissance totale (carte de remplissage ou ',/,&
'@      puissance des colis).                                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return

end subroutine
