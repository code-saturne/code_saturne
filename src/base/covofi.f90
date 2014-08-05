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

!===============================================================================
! Function:
! ---------

!> \file covofi.f90
!>
!> \brief This subroutine performs the solving the convection/diffusion
!> equation (with eventually source terms and/or drift) for a scalar quantity
!> over a time step.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     iscal         scalar number
!> \param[in]     itspdv        indicator to compute production/dissipation
!>                              terms for a variance:
!>                               - 0: no
!>                               - 1: yes
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        type of mass source term for the variables
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     tslagr        coupling term for the Lagrangian module
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!> \param[in]     frcxt         external forces making hydrostatic pressure
!> \param[in]     dfrcxt        variation of the external forces
!> \param[in]                    making the hydrostatic pressure
!> \param[in]     tpucou        non scalar time step in case of
!>                               velocity pressure coupling
!> \param[in]     trav          right hand side for the normalizing
!>                               the residual
!> \param[in]     viscf         visc*surface/dist aux faces internes
!> \param[in]     viscb         visc*surface/dist aux faces de bord
!> \param[in]     smbrs         tableau de travail
!> \param[in]     rovsdt        tableau de travail
!_______________________________________________________________________________

subroutine covofi &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  , itspdv ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , tslagr ,                   &
   ckupdc , smacel ,                                              &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use lagpar
use lagran
use radiat
use field
use field_operator
use ihmpre, only: iihmpr
use mesh
use parall
use period
use cs_f_interfaces
use atchem

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal  , itspdv

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision tslagr(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)

! Local variables

logical          lprev
character(len=80) :: chaine, fname
integer          ivar
integer          ifac  , iel
integer          iprev , inc   , iccocg, isqrt, iii, iiun, ibcl
integer          ivarsc
integer          iiscav
integer          ipcvsl, iflmas, iflmab
integer          ippvar, ipp   , ipcvso
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          imucpp, idftnp, iswdyp
integer          iflid , f_id, st_prv_id,  keydri, iscdri
integer          icvflb

integer          ivoid(1)

double precision epsrgp, climgp, extrap, relaxp, blencp, epsilp
double precision epsrsp
double precision rhovst, xk    , xe    , sclnor
double precision thetv , thets , thetap, thetp1
double precision smbexp
double precision temp, idifftp

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:) :: xcpp
double precision, allocatable, dimension(:) :: srcmas

double precision, dimension(:,:), pointer :: xut, visten
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom, croma, pcrom
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: porosi
double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_omg
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: visct, cpro_cp, c_st_scal

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Index of the field
iflid = ivarfl(isca(iscal))

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Id of the mass flux
call field_get_key_int(iflid, kimasf, iflmas) ! interior mass flux
! Pointer to the internal mass flux
call field_get_val_s(iflmas, imasfl)

! Id of the mass flux
call field_get_key_int(iflid, kbmasf, iflmab) ! boundary mass flux
! Pointer to the Boundary mass flux
call field_get_val_s(iflmab, bmasfl)

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

! --- Numero des grandeurs physiques
call field_get_val_s(icrom, crom)
call field_have_previous(icrom, lprev)
if (lprev) then
  call field_get_val_prev_s(icrom, croma)
endif
call field_get_val_s(iprpfl(ivisct), visct)

if (ivisls(iscal).gt.0) then
  ipcvsl = ipproc(ivisls(iscal))
else
  ipcvsl = 0
endif

! --- Numero du terme source dans PROPCE si extrapolation
call field_get_key_int(iflid, kstprv, st_prv_id)
if (st_prv_id .ge.0) then
  call field_get_val_s(st_prv_id, c_st_scal)
else
  c_st_scal => null()
endif

! S pour Source, V pour Variable
thets  = thetss(iscal)
thetv  = thetav(ivar)

call field_get_name(ivarfl(ivar), chaine)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:16)
endif

! When solving the Temperature, we solve:
!  cp*Vol*dT/dt + ...

imucpp = 0
if (iscavr(iscal).gt.0) then
  if (abs(iscacp(iscavr(iscal))).eq.1) then
    imucpp = 1
  endif
else
  if (abs(iscacp(iscal)).eq.1) then
    imucpp = 1
  endif
endif

allocate(xcpp(ncelet))

if (imucpp.eq.0) then
  do iel = 1, ncel
    xcpp(iel) = 1.d0
  enddo
elseif (imucpp.eq.1) then
  if (icp.gt.0) then
    call field_get_val_s(iprpfl(icp), cpro_cp)
    do iel = 1, ncel
      xcpp(iel) = cpro_cp(iel)
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
    ( iflid  , rtp(1,ivar), smbrs  , rovsdt )
  else
    call uitsth &
    ( iflid  , rtp(1,ivar), smbrs  , rovsdt )
  endif
endif

call ustssc &
!==========
( nvar   , nscal  , ncepdp , ncesmp ,                            &
  iscal  ,                                                       &
  icepdc , icetsm , itypsm ,                                     &
  dt     ,                                                       &
  ckupdc , smacel , smbrs  , rovsdt )

! Atmospheric chemistry
! In case of a semi-coupled resolution, computation of the explicit
! chemical source term to be considered during dynamical resolution
! The first nespg user scalars are supposed to be chemical species
if ((ichemistry.ge.1).and.(isepchemistry.eq.2)                    &
     .and.(iscal.le.nespg).and.(ntcabs.gt.1)) then
  call chem_source_terms(iscal, rtpa, smbrs, rovsdt)
endif


! Si on extrapole les TS :
!   SMBRS recoit -theta PROPCE du pas de temps precedent
!     (on aurait pu le faire avant ustssc, mais avec le risque que
!      l'utilisateur l'ecrase)
!   SMBRS recoit la partie du terme source qui depend de la variable
!   A l'ordre 2, on suppose que le ROVSDT fourni par l'utilisateur est <0
!     on implicite le terme (donc ROVSDT*RTPA va dans SMBRS)
!   En std, on adapte le traitement au signe de ROVSDT, mais ROVSDT*RTPA va
!     quand meme dans SMBRS (pas d'autre choix)
if (st_prv_id .ge. 0) then
  do iel = 1, ncel
    ! Stockage temporaire pour economiser un tableau
    smbexp = c_st_scal(iel)
    ! Terme source utilisateur explicite
    c_st_scal(iel) = smbrs(iel)
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
  call cptssy(iscal, rtpa, smbrs, rovsdt)
  !==========
endif

! --> Physique particulieres
!     Ordre 2 non pris en compte

if (ippmod(iphpar).ge.1) then
  call pptssc &
  !==========
 ( iscal  ,                                                       &
   rtpa   , rtp    , propce ,                                     &
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
  endif

  !-> Charbon pulverise
  !   Ordre 2 non pris en compte
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

  if (iscal.eq.iscalt .and. (itherm.eq.1 .or. itherm.eq.2)) then

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
  if (st_prv_id .ge. 0) then
    do iel = 1, ncel
      c_st_scal(iel) = c_st_scal(iel) + w1(iel)
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
    allocate(grad(3,ncelet))

    ! Remarque : on a prevu la possibilite de scalaire associe non
    !  variable de calcul, mais des adaptations sont requises

    if (ivarsc.gt.0) then
      iii = ivarsc
    else
      write(nfecra,9000)ivarsc
      call csexit(1)
    endif

    iprev = 1
    inc = 1
    iccocg = 1

    call field_gradient_scalar(ivarfl(iii), iprev, imrgra, inc,   &
                               iccocg,                            &
                               grad)

    ! Traitement de la production
    ! On utilise MAX(PROPCE,ZERO) car en LES dynamique on fait un clipping
    ! tel que (mu + mu_t)>0, donc mu_t peut etre negatif et donc
    ! potentiellement (lambda/Cp + mu_t/sigma) aussi
    ! Ceci ne pose probleme que quand on resout une equation de variance
    ! de scalaire avec un modele LES ... ce qui serait curieux mais n'est
    ! pas interdit par le code.
    !   Si extrapolation : dans PROPCE
    if (st_prv_id .ge. 0) then
      ! On prend la viscosite a l'instant n, meme si elle est extrapolee
      ipcvso = ipproc(ivisct)
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
          c_st_scal(iel) = c_st_scal(iel) -2.d0*xcpp(iel)*volume(iel) &
                                          *(xut(1,iel)*grad(1,iel)    &
                                           +xut(2,iel)*grad(2,iel)    &
                                           +xut(3,iel)*grad(3,iel) )
        enddo
      ! SGDH model
      else
        do iel = 1, ncel
          c_st_scal(iel) = c_st_scal(iel)                                     &
               + 2.d0*xcpp(iel)*max(propce(iel,ipcvso),zero)                  &
               *volume(iel)/sigmas(iscal)                                     &
               *(grad(1,iel)**2 + grad(2,iel)**2 + grad(3,iel)**2)
        enddo
      endif
    ! Sinon : dans SMBRS
    else
      ipcvso = ipproc(ivisct)

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
                                  *(xut(1,iel)*grad(1,iel)      &
                                   +xut(2,iel)*grad(2,iel)      &
                                   +xut(3,iel)*grad(3,iel) )
        enddo

      ! SGDH model
      else
        do iel = 1, ncel
          smbrs(iel) = smbrs(iel)                                            &
                     + 2.d0*xcpp(iel)*max(propce(iel,ipcvso),zero)           &
                     * volume(iel)/sigmas(iscal)                             &
                     * (grad(1,iel)**2 + grad(2,iel)**2 + grad(3,iel)**2)
        enddo
      endif

      ! Production term for a variance  TODO compute ustdy when isso2t >0
      if (idilat.eq.4) then
        do iel = 1, ncel
          propce(iel,ipproc(iustdy(iscal))) =                     &
          propce(iel,ipproc(iustdy(iscal))) +                     &
               2.d0*xcpp(iel)*max(propce(iel,ipcvso),zero)        &
             *volume(iel)/sigmas(iscal)                           &
             *(grad(1,iel)**2 + grad(2,iel)**2 + grad(3,iel)**2)
        enddo
      endif
    endif

    ! Free memory
    deallocate(grad)

    ! Traitement de la dissipation
    if (st_prv_id .ge. 0) then
      thetap = thetv
    else
      thetap = 1.d0
    endif

    if (itytur.eq.2 .or. itytur.eq.5) then
      call field_get_val_prev_s(ivarfl(ik), cvara_k)
      call field_get_val_prev_s(ivarfl(iep), cvara_ep)
    elseif (itytur.eq.3) then
      call field_get_val_prev_s(ivarfl(iep), cvara_ep)
      call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
      call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
      call field_get_val_prev_s(ivarfl(ir33), cvara_r33)
    elseif(iturb.eq.60) then
      call field_get_val_prev_s(ivarfl(ik), cvara_k)
      call field_get_val_prev_s(ivarfl(iomg), cvara_omg)
    endif

    do iel = 1, ncel
      if (itytur.eq.2 .or. itytur.eq.5) then
        xk = cvara_k(iel)
        xe = cvara_ep(iel)
      elseif (itytur.eq.3) then
        xk = 0.5d0*(cvara_r11(iel)+cvara_r22(iel)+cvara_r33(iel))
        xe = cvara_ep(iel)
      elseif(iturb.eq.60) then
        xk = cvara_k(iel)
        xe = cmu*xk*cvara_omg(iel)
      endif
      rhovst = xcpp(iel)*crom(iel)*xe/(xk * rvarfl(iscal))       &
             *volume(iel)

      ! La diagonale recoit eps/Rk, (*theta eventuellement)
      rovsdt(iel) = rovsdt(iel) + rhovst*thetap
      ! SMBRS recoit la dissipation
      smbrs(iel) = smbrs(iel) - rhovst*rtpa(iel,ivar)
      ! Dissipation term for a variance
      if (idilat.eq.4) then
        propce(iel,ipproc(iustdy(iscal))) =                               &
          propce(iel,ipproc(iustdy(iscal))) - xcpp(iel)*rhovst*rtpa(iel,ivar)
      endif
    enddo

  endif

endif

if (st_prv_id .ge. 0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + thetp1 * c_st_scal(iel)
  enddo
endif

! Low Mach compressible algos (conservative in time). Same algo for cavitation.
if (idilat.gt.1 .or. icavit.ge.0) then
  call field_get_val_prev_s(icrom, pcrom)

! Standard algo
else
  call field_get_val_s(icrom, pcrom)
endif

! Get the the order of the diffusivity tensor of the variable
call field_get_key_int_by_name(ivarfl(ivar), "diffusivity_tensor", idftnp)

! "VITESSE" DE DIFFUSION FACETTE

! On prend le MAX(mu_t,0) car en LES dynamique mu_t peut etre negatif
! (clipping sur (mu + mu_t)). On aurait pu prendre
! MAX(K + K_t,0) mais cela autoriserait des K_t negatif, ce qui est
! considere ici comme non physique.
if (idiff(ivar).ge.1) then
  ! Scalar diffusivity
  if (idftnp.eq.1) then

    idifftp = idifft(ivar)
    if (ityturt(iscal).eq.3) then
      idifftp = 0
    endif
    if (ipcvsl.eq.0) then
      do iel = 1, ncel
        w1(iel) = visls0(iscal)                                     &
           + idifftp*xcpp(iel)*max(visct(iel),zero)/sigmas(iscal)
      enddo
    else
      do iel = 1, ncel
        w1(iel) = propce(iel,ipcvsl)                                &
           + idifftp*xcpp(iel)*max(visct(iel),zero)/sigmas(iscal)
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

    call field_get_val_v(ivsten, visten)

    if (ipcvsl.eq.0) then
      do iel = 1, ncel

        temp = idifft(ivar)*xcpp(iel)*ctheta(iscal)/csrij
        viscce(1,iel) = temp*visten(1,iel) + visls0(iscal)
        viscce(2,iel) = temp*visten(2,iel) + visls0(iscal)
        viscce(3,iel) = temp*visten(3,iel) + visls0(iscal)
        viscce(4,iel) = temp*visten(4,iel)
        viscce(5,iel) = temp*visten(5,iel)
        viscce(6,iel) = temp*visten(6,iel)

      enddo
    else
      do iel = 1, ncel

        temp = idifft(ivar)*xcpp(iel)*ctheta(iscal)/csrij
        viscce(1,iel) = temp*visten(1,iel) + propce(iel,ipcvsl)
        viscce(2,iel) = temp*visten(2,iel) + propce(iel,ipcvsl)
        viscce(3,iel) = temp*visten(3,iel) + propce(iel,ipcvsl)
        viscce(4,iel) = temp*visten(4,iel)
        viscce(5,iel) = temp*visten(5,iel)
        viscce(6,iel) = temp*visten(6,iel)

      enddo
    endif

    iwarnp = iwarni(ivar)

    call vitens &
    !==========
   ( viscce , iwarnp ,             &
     weighf , weighb ,             &
     viscf  , viscb  )

  endif

  ! AFM model or DFM models: add div(Cp*rho*T'u') to smbrs
  ! Compute T'u' for GGDH
  if (ityturt(iscal).ge.1) then

    call divrit &
    !==========
    ( nscal  ,                                                       &
      iscal  ,                                                       &
      dt     , rtp    , rtpa   ,                                     &
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

  ! --> Unsteady term and mass aggregation term
  do iel = 1, ncel
    rovsdt(iel) = rovsdt(iel)                                                 &
                + istat(ivar)*xcpp(iel)*pcrom(iel)*volume(iel)/dt(iel)
  enddo

! With porosity
else

  call field_get_val_s(ipori, porosi)

  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)*porosi(iel)
  enddo

  ! --> Unsteady term and mass aggregation term
  do iel = 1, ncel
    rovsdt(iel) = ( rovsdt(iel)                                             &
                  + istat(ivar)*xcpp(iel)*pcrom(iel)*volume(iel)/dt(iel)    &
                  ) * porosi(iel)
  enddo

endif

! Scalar with a Drift:
! compute the convective flux
!----------------------------

call field_get_key_int(iflid, keydri, iscdri)

if (btest(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)) then
 call driflu &
 !=========
 ( iflid  ,                                                       &
   dt     , rtpa   , propce ,                                     &
   imasfl , bmasfl ,                                              &
   rovsdt , smbrs  )
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
! all boundary convective flux with upwind
icvflb = 0

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

call codits &
!==========
 ( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscce , viscf  , viscb  , viscce ,          &
   weighf , weighb ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbrs  , rtp(1,ivar)     , dpvar  ,                   &
   xcpp   , rvoid  )

!===============================================================================
! 4. Writing and clipping
!===============================================================================

if (ivarsc.gt.0) then
  iii = ivarsc
else
! Valeur bidon
  iii = nflown
endif

call clpsca(ncelet, ncel, iscal, rtp(1,iii), rtp)
!==========

! BILAN EXPLICITE (VOIR CODITS : ON ENLEVE L'INCREMENT)
! Ceci devrait etre valable avec le theta schema sur les Termes source

if (iwarni(ivar).ge.2) then
  if (nswrsm(ivar).gt.1) then
    ibcl = 1
  else
    ibcl = 0
  endif
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                                 &
            - istat(ivar)*xcpp(iel)*(pcrom(iel)/dt(iel))*volume(iel)&
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
