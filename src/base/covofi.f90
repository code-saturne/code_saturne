!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine covofi &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  , itspdv ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb , tslagr , &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt )

!===============================================================================
! FONCTION :
! ----------

! Solving the advection/diffusion equation (with source terms) for a scalar
! quantity over a time step

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! itspdv           ! e  ! <-- ! calcul termes sources prod et dissip           !
!                  !    !     !  (0 : non , 1 : oui)                           !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
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
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use ppthch
use pointe, only: porosi, visten
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use lagpar
use lagran
use radiat
use field
use ihmpre, only: iihmpr
use mesh
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal  , itspdv

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)

! Local variables

character*80     chaine
integer          ivar
integer          ifac  , iel
integer          init  , inc   , iccocg, isqrt, iii, iiun, ibcl
integer          ivarsc
integer          iiscav, iicp
integer          iclvar, iclvaf
integer          ipcrom, ipcroa, ipcrho, ipcvst, ipcvsl, iflmas, iflmab
integer          ippvar, ipp   , iptsca, ipcvso
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          imucpp, idftnp, iswdyp
integer          f_id

double precision epsrgp, climgp, extrap, relaxp, blencp, epsilp
double precision epsrsp
double precision rhovst, xk    , xe    , sclnor
double precision thetv , thets , thetap, thetp1
double precision smbexp, drtp
double precision trrij , csteps

double precision rvoid(1)

character*80     fname

double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:) :: xcpp
double precision, allocatable, dimension(:) :: srcmas

double precision, dimension(:,:), pointer :: xut

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays
allocate(w1(ncelet))
allocate(dpvar(ncelet))

! Initialize variables to avoid compiler warnings

xe = 0.d0
xk = 0.d0

! --- Numero de variable de calcul et de post associe au scalaire traite
ivar   = isca(iscal)
ippvar = ipprtp(ivar)

! --- Numero du scalaire eventuel associe dans le cas fluctuation
!         et numero de variable de calcul
iiscav = iscavr(iscal)
if (iiscav.gt.0.and.iiscav.le.nscal) then
  ivarsc = isca(iiscav)
else
  ivarsc = 0
endif

! --- Numero des conditions aux limites
iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

! --- Numero des grandeurs physiques
ipcrom = ipproc(irom)
if (iroma .gt. 0) then
  ipcroa = ipproc(iroma)
else
  ipcroa = 0
endif
ipcvst = ipproc(ivisct)
iflmas = ipprof(ifluma(ivar ))
iflmab = ipprob(ifluma(ivar ))
if (ivisls(iscal).gt.0) then
  ipcvsl = ipproc(ivisls(iscal))
else
  ipcvsl = 0
endif

! --- Numero du terme source dans PROPCE si extrapolation
if (isso2t(iscal).gt.0) then
  iptsca = ipproc(itssca(iscal))
else
  iptsca = 0
endif

! S pour Source, V pour Variable
thets  = thetss(iscal)
thetv  = thetav(ivar )

chaine = nomvar(ippvar)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:16)
endif

! When solving the Temperature, we solve:
!  cp*Vol*dT/dt + ...
if (iscalt.gt.0) then
  if (ivar.eq.isca(iscalt) .or. iscavr(iscal).eq.iscalt) then
    if (abs(iscsth(iscalt)).eq.1) then
      imucpp = 1
    else
      imucpp = 0
    endif
  else
    imucpp = 0
  endif
else
  imucpp = 0
endif

allocate(xcpp(ncelet))

if (imucpp.eq.0) then
  do iel = 1, ncel
    xcpp(iel) = 1.d0
  enddo
elseif (imucpp.eq.1) then
  if (icp.gt.0) then
    do iel = 1, ncel
      xcpp(iel) = propce(iel,ipproc(icp))
    enddo
  else
    do iel = 1, ncel
      xcpp(iel) = cp0
    enddo
  endif
endif

! Handle parallelism and periodicity
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(xcpp)
endif

!===============================================================================
! 2. Source terms
!===============================================================================

! --> Initialization

do iel = 1, ncel
  rovsdt(iel) = 0.d0
  smbrs(iel) = 0.d0
enddo

if (iihmpr.eq.1) then

  if (iscal.ne.iscalt) then
    call uitssc &
    ( iscal  , rtp(1,ivar), smbrs  , rovsdt )
  else
    call uitsth &
    ( iscal  , rtp(1,ivar), smbrs  , rovsdt )
  endif
endif

call ustssc &
!==========
( nvar   , nscal  , ncepdp , ncesmp ,                            &
  iscal  ,                                                       &
  icepdc , icetsm , itypsm ,                                     &
  dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
  ckupdc , smacel , smbrs  , rovsdt )

! Si on extrapole les TS :
!   SMBRS recoit -theta PROPCE du pas de temps precedent
!     (on aurait pu le faire avant ustssc, mais avec le risque que
!      l'utilisateur l'ecrase)
!   SMBRS recoit la partie du terme source qui depend de la variable
!   A l'ordre 2, on suppose que le ROVSDT fourni par l'utilisateur est <0
!     on implicite le terme (donc ROVSDT*RTPA va dans SMBRS)
!   En std, on adapte le traitement au signe de ROVSDT, mais ROVSDT*RTPA va
!     quand meme dans SMBRS (pas d'autre choix)
if (isso2t(iscal).gt.0) then
  do iel = 1, ncel
    ! Stockage temporaire pour economiser un tableau
    smbexp = propce(iel,iptsca)
    ! Terme source utilisateur explicite
    propce(iel,iptsca) = smbrs(iel)
    ! Terme source du pas de temps precedent et
    ! On suppose -ROVSDT > 0 : on implicite
    !    le terme source utilisateur (le reste)
    smbrs(iel) = rovsdt(iel)*rtpa(iel,ivar) - thets*smbexp
    ! Diagonale
    rovsdt(iel) = - thetv*rovsdt(iel)
  enddo

! Si on n'extrapole pas les TS :
else
  do iel = 1, ncel
    ! Terme source utilisateur
    smbrs(iel) = smbrs(iel) + rovsdt(iel)*rtpa(iel,ivar)
    ! Diagonale
    rovsdt(iel) = max(-rovsdt(iel),zero)
  enddo
endif

! Add thermodynamic pressure variation for the low-Mach algorithm:
! NB: iscalt is the Enthalpy
if (idilat.eq.3 .and. iscalt.gt.0) then
  if (ivar.eq.isca(iscalt)) then
    ! unsteady thermodynamic source term added
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) + (pther - pthera)/dt(iel)*volume(iel)
    enddo
  endif
endif

! --> Couplage volumique avec Syrthes
!     Ordre 2 non pris en compte

if (iscal.eq.iscalt) then
  call cptssy &
  !==========
( nvar   , nscal  ,                                              &
  iscal  ,                                                       &
  dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
  smbrs  , rovsdt )
endif

! --> Physique particulieres
!     Ordre 2 non pris en compte

if (ippmod(iphpar).ge.1) then
  call pptssc &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt , tslagr )
endif

! --> Rayonnement
!     Ordre 2 non pris en compte

if (iirayo.ge.1) then

  if (iscal.eq.iscalt) then
    call raysca &
    !==========
  ( iscalt,ncelet,ncel,     &
    smbrs, rovsdt,volume,propce )

    ! Store the explicit radiative source term
    if (idilat.eq.4) then
      do iel = 1, ncel
        propce(iel,ipproc(iustdy(ihm))) = &
        propce(iel,ipproc(iustdy(ihm)))   &
        + propce(iel,ipproc(itsre(1)))*volume(iel)
      enddo
    endif
  endif

  !-> Charbon pulverise
  !   Ordre 2 non pris en compte
  ! old model
  if ( ippmod(icp3pl) .ge. 0 ) then
    if ( isca(iscal).ge.isca(ih2(1)) .and.       &
         isca(iscal).le.isca(ih2(nclacp)) ) then

      call cprays &
      !==========
    ( ivar  ,ncelet, ncel  ,       &
      volume,propce,smbrs,rovsdt)

    endif
  endif
  ! new model
  if (ippmod(iccoal) .ge. 0) then
    if (isca(iscal).ge.isca(ih2(1)) .and.       &
        isca(iscal).le.isca(ih2(nclacp))) then

      call cs_coal_radst &
      !=================
      ( ivar   , ncelet , ncel  ,                &
        volume , propce , smbrs , rovsdt )

    endif
  endif

  ! -> Fuel
  !    Ordre 2 non pris en compte
  !    Pour l'instant rayonnement non compatible avec Fuel

  if (ippmod(icfuel) .ge. 0) then
    if (isca(iscal).ge.isca(ih2(1)) .and.       &
        isca(iscal).le.isca(ih2(nclafu))) then

      call cs_fuel_radst &
     !==================
    ( ivar  ,ncelet, ncel  ,                      &
      volume,rtpa  , propce,smbrs,rovsdt)

    endif
  endif

endif

! --> Lagrangien (couplage retour thermique)
!     Ordre 2 non pris en compte

if (iilagr.eq.2 .and. ltsthe.eq.1)  then

  if ((iscsth(iscal).eq.2).or.(abs(iscsth(iscal)).eq.1)) then

    do iel = 1, ncel
      smbrs (iel) = smbrs(iel)  + tslagr(iel,itste)
      rovsdt(iel) = rovsdt(iel) + xcpp(iel)*max(tslagr(iel,itsti),zero)
    enddo

  endif

endif

! Mass source term
if (ncesmp.gt.0) then

  ! Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

  allocate(srcmas(ncesmp))

  ! When treating the Temperature, the equation is multiplied by Cp
  do iel = 1, ncesmp
    if (smacel(iel,ipr).gt.0.d0 .and.itypsm(iel,ivar).eq.1) then
      srcmas(iel) = smacel(iel,ipr)*xcpp(icetsm(iel))
    else
      srcmas(iel) = 0.d0
    endif
  enddo

  ! On incremente SMBRS par -Gamma RTPA et ROVSDT par Gamma (*theta)
  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   , isso2t(iscal) , thetv  ,   &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , srcmas   ,            &
   smbrs  , rovsdt , w1)

  deallocate(srcmas)

  ! Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if (isso2t(iscal).gt.0) then
    do iel = 1, ncel
      propce(iel,iptsca) = propce(iel,iptsca) + w1(iel)
    enddo
  ! Sinon on le met directement dans SMBRS
  else
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) + w1(iel)
    enddo
  endif

endif

! If the current scalar is the variance of an other scalar,
! production and dissipation terms are added.
if (itspdv.eq.1) then

  if (itytur.eq.2 .or. itytur.eq.3 .or. itytur.eq.5 .or. iturb.eq.60) then

    ! Allocate a temporary array for the gradient reconstruction
    allocate(grad(ncelet,3))

    ! Remarque : on a prevu la possibilite de scalaire associe non
    !  variable de calcul, mais des adaptations sont requises

    if (ivarsc.gt.0) then
      iii = ivarsc
    else
      write(nfecra,9000)ivarsc
      call csexit(1)
    endif

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(iii)
    imligp = imligr(iii)
    iwarnp = iwarni(iii)
    epsrgp = epsrgr(iii)
    climgp = climgr(iii)
    extrap = extrag(iii)

    call grdcel &
    !==========
  ( iii    , imrgra , inc    , iccocg , nswrgp , imligp ,          &
    iwarnp , nfecra ,                                              &
    epsrgp , climgp , extrap ,                                     &
    rtpa(1,iii) , coefa(1,iclrtp(iii,icoef)) ,                     &
                  coefb(1,iclrtp(iii,icoef)) ,                     &
    grad   )

    ! Traitement de la production
    ! On utilise MAX(PROPCE,ZERO) car en LES dynamique on fait un clipping
    ! tel que (mu + mu_t)>0, donc mu_t peut etre negatif et donc
    ! potentiellement (lambda/Cp + mu_t/sigma) aussi
    ! Ceci ne pose probleme que quand on resout une equation de variance
    ! de scalaire avec un modele LES ... ce qui serait curieux mais n'est
    ! pas interdit par le code.
    !   Si extrapolation : dans PROPCE
    if (isso2t(iscal).gt.0) then
      ! On prend la viscosite a l'instant n, meme si elle est extrapolee
      ipcvso = ipcvst
      if (iviext.gt.0) ipcvso = ipproc(ivista)


      ! iscal is the variance of the scalar iiscav
      ! with modelized turbulent fluxes GGDH or AFM or DFM
      if (ityturt(iiscav).ge.1) then

        ! Name of the scalar iiscav associated to the variance iscal
        call field_get_name(ivarfl(ivarsc), fname)

        ! Index of the corresponding turbulent flux
        call field_get_id(trim(fname)//'_turbulent_flux', f_id)

        call field_get_val_v(f_id, xut)

        do iel = 1, ncel
          propce(iel,iptsca) = propce(iel,iptsca) -2.d0*xcpp(iel)*volume(iel) &
                                                  *(xut(1,iel)*grad(iel,1)    &
                                                   +xut(2,iel)*grad(iel,2)    &
                                                   +xut(3,iel)*grad(iel,3) )
        enddo
      ! SGDH model
      else
        do iel = 1, ncel
          propce(iel,iptsca) = propce(iel,iptsca)                             &
               + 2.d0*xcpp(iel)*max(propce(iel,ipcvso),zero)                  &
               *volume(iel)/sigmas(iscal)                                     &
               *(grad(iel,1)**2 + grad(iel,2)**2 + grad(iel,3)**2)
        enddo
      endif
    ! Sinon : dans SMBRS
    else
      ipcvso = ipcvst

      ! iscal is the variance of the scalar iiscav
      ! with modelized turbulent fluxes GGDH or AFM or DFM
      if (ityturt(iiscav).ge.1) then


        ! Name of the scalar ivarsc associated to the variance iscal
        call field_get_name(ivarfl(ivarsc), fname)

        ! Index of the corresponding turbulent flux
        call field_get_id(trim(fname)//'_turbulent_flux', f_id)

        call field_get_val_v(f_id, xut)

        do iel = 1, ncel
          smbrs(iel) = smbrs(iel) -2.d0*xcpp(iel)*volume(iel)   &
                                  *(xut(1,iel)*grad(iel,1)      &
                                   +xut(2,iel)*grad(iel,2)      &
                                   +xut(3,iel)*grad(iel,3) )
        enddo

      ! SGDH model
      else
        do iel = 1, ncel
          smbrs(iel) = smbrs(iel)                                            &
                     + 2.d0*xcpp(iel)*max(propce(iel,ipcvso),zero)           &
                     * volume(iel)/sigmas(iscal)                             &
                     * (grad(iel,1)**2 + grad(iel,2)**2 + grad(iel,3)**2)
        enddo
      endif

      ! Production term for a variance  TODO compute ustdy when isso2t >0
      if (idilat.eq.4) then
        do iel = 1, ncel
          propce(iel,ipproc(iustdy(iscal))) =                     &
          propce(iel,ipproc(iustdy(iscal))) +                     &
               2.d0*xcpp(iel)*max(propce(iel,ipcvso),zero)        &
             *volume(iel)/sigmas(iscal)                           &
             *(grad(iel,1)**2 + grad(iel,2)**2 + grad(iel,3)**2)
        enddo
      endif
    endif

    ! Free memory
    deallocate(grad)

    ! Traitement de la dissipation
    if (isso2t(iscal).gt.0) then
      thetap = thetv
    else
      thetap = 1.d0
    endif
    do iel = 1, ncel
      if (itytur.eq.2 .or. itytur.eq.5) then
        xk = rtpa(iel,ik)
        xe = rtpa(iel,iep)
      elseif (itytur.eq.3) then
        xk = 0.5d0*(rtpa(iel,ir11)+rtpa(iel,ir22)+rtpa(iel,ir33))
        xe = rtpa(iel,iep)
      elseif(iturb.eq.60) then
        xk = rtpa(iel,ik)
        xe = cmu*xk*rtpa(iel,iomg)
      endif
      rhovst = xcpp(iel)*propce(iel,ipcrom)*xe/(xk * rvarfl(iscal))       &
             *volume(iel)

      ! La diagonale recoit eps/Rk, (*theta eventuellement)
      rovsdt(iel) = rovsdt(iel) + rhovst*thetap
      ! SMBRS recoit la dissipation
      smbrs(iel) = smbrs(iel) - rhovst*rtpa(iel,ivar)
    enddo

  endif

endif

if (isso2t(iscal).gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + thetp1 * propce(iel,iptsca)
  enddo
endif

! Low Mach compressible algos (conservative in time)
if (idilat.gt.1) then
  ipcrho = ipcroa

! Standard algo
else
  ipcrho = ipcrom
endif

idftnp = idften(ivar)

! "VITESSE" DE DIFFUSION FACETTE

! On prend le MAX(mu_t,0) car en LES dynamique mu_t peut etre negatif
! (clipping sur (mu + mu_t)). On aurait pu prendre
! MAX(K + K_t,0) mais cela autoriserait des K_t negatif, ce qui est
! considere ici comme non physique.
if (idiff(ivar).ge.1) then
  ! Scalar diffusivity
  if (idftnp.eq.1) then
    if (ipcvsl.eq.0) then
      do iel = 1, ncel
        w1(iel) = visls0(iscal)                                     &
           + idifft(ivar)*xcpp(iel)*max(propce(iel,ipcvst),zero)/sigmas(iscal)
      enddo
    else
      do iel = 1, ncel
        w1(iel) = propce(iel,ipcvsl)                                &
           + idifft(ivar)*xcpp(iel)*max(propce(iel,ipcvst),zero)/sigmas(iscal)
      enddo
    endif

    call viscfa &
    !==========
   ( imvisf ,                      &
     w1     ,                      &
     viscf  , viscb  )

  ! Symmetric tensor diffusivity (GGDH)
  elseif (idftnp.eq.6) then

    ! Allocate temporary arrays
    allocate(viscce(6,ncelet))
    allocate(weighf(2,nfac))
    allocate(weighb(nfabor))

    if (ipcvsl.eq.0) then
      do iel = 1, ncel

        viscce(1,iel) = visls0(iscal)                                      &
                      + idifft(ivar)*xcpp(iel)*visten(1,iel)*ctheta(iscal)
        viscce(2,iel) = visls0(iscal)                                      &
                      + idifft(ivar)*xcpp(iel)*visten(2,iel)*ctheta(iscal)
        viscce(3,iel) = visls0(iscal)                                      &
                      + idifft(ivar)*xcpp(iel)*visten(3,iel)*ctheta(iscal)
        viscce(4,iel) = idifft(ivar)*xcpp(iel)*visten(4,iel)*ctheta(iscal)
        viscce(5,iel) = idifft(ivar)*xcpp(iel)*visten(5,iel)*ctheta(iscal)
        viscce(6,iel) = idifft(ivar)*xcpp(iel)*visten(6,iel)*ctheta(iscal)

      enddo
    else
      do iel = 1, ncel

        viscce(1,iel) = propce(iel,ipcvsl)                                 &
                      + idifft(ivar)*xcpp(iel)*visten(1,iel)*ctheta(iscal)
        viscce(2,iel) = propce(iel,ipcvsl)                                 &
                      + idifft(ivar)*xcpp(iel)*visten(2,iel)*ctheta(iscal)
        viscce(3,iel) = propce(iel,ipcvsl)                                 &
                      + idifft(ivar)*xcpp(iel)*visten(3,iel)*ctheta(iscal)
        viscce(4,iel) = idifft(ivar)*xcpp(iel)*visten(4,iel)*ctheta(iscal)
        viscce(5,iel) = idifft(ivar)*xcpp(iel)*visten(5,iel)*ctheta(iscal)
        viscce(6,iel) = idifft(ivar)*xcpp(iel)*visten(6,iel)*ctheta(iscal)

      enddo
    endif

    iwarnp = iwarni(ivar)

    call vitens &
    !==========
   ( imvisf ,                      &
     viscce , iwarnp ,             &
     weighf , weighb ,             &
     viscf  , viscb  )

  endif

  ! AFM model or DFM models: add div(Cp*rho*T'u') to smbrs
  ! Compute T'u' for GGDH
  if (ityturt(iscal).ge.1) then

    call divrit &
    !==========
    ( nvar   , nscal  ,                                              &
      iscal  , itspdv ,                                              &
      dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
      coefa  , coefb  ,                                              &
      xcpp   ,                                                       &
      smbrs  )

  endif

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

! Without porosity
if (iporos.eq.0) then

  ! --> Non stationnary term and mass aggregation term
  do iel = 1, ncel
    rovsdt(iel) = rovsdt(iel)                                                 &
                + istat(ivar)*xcpp(iel)*propce(iel,ipcrho)*volume(iel)/dt(iel)
  enddo

! With porosity
else

  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)*porosi(iel)
  enddo

  ! --> Non stationnary term and mass aggregation term
  do iel = 1, ncel
    rovsdt(iel) = ( rovsdt(iel)                                        &
                  + istat(ivar)*xcpp(iel)*propce(iel,ipcrho)*volume(iel)/dt(iel) &
                  ) * porosi(iel)
  enddo

endif

!===============================================================================
! 3. Solving
!===============================================================================
iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
ndircp = ndircl(ivar)
nitmap = nitmax(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
iescap = 0
iswdyp = iswdyn(ivar)
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
ipp    = ippvar
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)

call codits &
!==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   propfa(1,iflmas), propfb(1,iflmab),                            &
   viscf  , viscb  , viscce , viscf  , viscb  , viscce ,          &
   weighf , weighb ,                                              &
   rovsdt , smbrs  , rtp(1,ivar)     , dpvar  ,                   &
   xcpp   , rvoid  )

!===============================================================================
! 4. Writting and clipping
!===============================================================================

if (ivarsc.gt.0) then
  iii = ivarsc
else
! Valeur bidon
  iii = 1
endif

call clpsca &
!==========
 ( ncelet , ncel   , nvar   , nscal  , iscal  ,                   &
   propce , rtp(1,iii)      , rtp    )

! Store the dissipation term for a variance
if (idilat.eq.4.and.itspdv.eq.1) then

  do iel = 1, ncel

    if (itytur.eq.2 .or. itytur.eq.5) then
      xk = rtpa(iel,ik)
      xe = rtpa(iel,iep)
    elseif (itytur.eq.3) then
      xk = 0.5d0*(rtpa(iel,ir11)+rtpa(iel,ir22)+rtpa(iel,ir33))
      xe = rtpa(iel,iep)
    elseif(iturb.eq.60) then
      xk = rtpa(iel,ik)
      xe = cmu*xk*rtpa(iel,iomg)
    endif

    rhovst = xcpp(iel)*propce(iel,ipcrom)*xe/(xk * rvarfl(iscal))     &
           *volume(iel)
    propce(iel,ipproc(iustdy(iscal))) =                               &
      propce(iel,ipproc(iustdy(iscal))) - rhovst*rtp(iel,ivar)

  enddo

endif

! Store the implicit part of the radiative source term
if (idilat.eq.4.and.iirayo.ge.1.and.iscal.eq.iscalt) then
  do iel = 1, ncel
    ivar = isca(iscalt)
    drtp = rtp(iel,ivar)-rtpa(iel,ivar)
    propce(iel,ipproc(iustdy(ihm))) = &
    propce(iel,ipproc(iustdy(ihm)))   &
    - propce(iel,ipproc(itsri(1)))*drtp*volume(iel)
  enddo
endif

! BILAN EXPLICITE (VOIR CODITS : ON ENLEVE L'INCREMENT)
! Ceci devrait etre valable avec le theta schema sur les Termes source

if (iwarni(ivar).ge.2) then
  if (nswrsm(ivar).gt.1) then
    ibcl = 1
  else
    ibcl = 0
  endif
  ! Low Mach compressible algos (conservative in time)
  if (idilat.gt.1) then
    ipcrho = ipcroa

    ! Standard algo
  else
    ipcrho = ipcrom
  endif
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                                 &
            - istat(ivar)*xcpp(iel)*(propce(iel,ipcrom)/dt(iel))*volume(iel)&
                *(rtp(iel,ivar)-rtpa(iel,ivar))*ibcl
  enddo
  isqrt = 1
  call prodsc(ncel,isqrt,smbrs,smbrs,sclnor)
  write(nfecra,1200)chaine(1:16) ,sclnor
endif

! Free memory
deallocate(w1)
if (allocated(viscce)) deallocate(viscce)
if (allocated(weighf)) deallocate(weighf, weighb)
deallocate(dpvar)
deallocate(xcpp)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A16                       ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A16,' : BILAN EXPLICITE = ',E14.5)
 9000 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS COVOFI                          ',/,&
'@    =========                                               ',/,&
'@    IVARSC DOIT ETRE UN ENTIER POSITIF STRICTEMENT          ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(/,                                                   &
'   ** SOLVING VARIABLE ',A16                                  ,/,&
'      ----------------'                                       ,/)
 1200 format(1X,A16,' : EXPLICIT BALANCE = ',E14.5)
 9000 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ERROR IN COVOFI'                                ,/,&
'@    ========'                                                ,/,&
'@    IVARSC MUST BE A STRICTLY POSITIVE INTEGER'              ,/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run.'                          ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! End
!----

return

end subroutine
