!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
!> Please refer to the
!> <a href="../../theory.pdf#covofi"><b>covofi</b></a> section of the
!> theory guide for more informations.
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
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ncmast        number of cells with condensation source terms
!> \param[in]     iterns        Navier-Stokes iteration number
!> \param[in]     iscal         scalar number
!> \param[in]     itspdv        indicator to compute production/dissipation
!>                              terms for a variance:
!>                               - 0: no
!>                               - 1: yes
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     ifbpcd        index of faces with condensation source terms
!> \param[in]     ltmast        index of cells with condensation source terms
!> \param[in]     itypsm        type of mass source term for the variables
!> \param[in]     itypcd        type of surface condensation source term
!> \param[in]     itypst        type of volume  condensation source term
!> \param[in]     dt            time step (per cell)
!> \param[in]     tslagr        coupling term for the Lagrangian module
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!> \param[in]     spcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, spcond is the flow rate
!>                              \f$ \Gamma_{s, cond}^n \f$)
!> \param[in]     svcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, svcond is the flow rate
!>                              \f$ \Gamma_{v, cond}^n \f$)
!> \param[in]     flxmst        variable value associated to heat transfer flux
!>                              associated to the metal mass condensation
!> \param[in]     viscf         visc*surface/dist at internal faces
!> \param[in]     viscb         visc*surface/dist at boundary faces
!_______________________________________________________________________________

subroutine covofi &
 ( nvar   , nscal  , ncepdp , ncesmp , nfbpcd , ncmast ,          &
   iterns , iscal  , itspdv ,                                     &
   icepdc , icetsm , ifbpcd , ltmast ,                            &
   itypsm , itypcd , itypst ,                                     &
   dt     , tslagr ,                                              &
   ckupdc , smacel , spcond , svcond , flxmst ,                   &
   viscf  , viscb  )

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
use ppcpfu
use lagran
use radiat
use field
use field_operator
use mesh
use parall
use period
use cs_f_interfaces
use atchem
use darcy_module
use cs_c_bindings
use pointe, only: itypfb, pmapper_double_r1
use atincl, only: kopint

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp , nfbpcd ,  ncmast
integer          iterns , iscal  , itspdv

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ifbpcd(nfbpcd), itypcd(nfbpcd,nvar)
integer          ltmast(ncelet), itypst(ncelet,nvar)

double precision dt(ncelet)
double precision tslagr(ncelet,*)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision spcond(nfbpcd,nvar)
double precision svcond(ncelet,nvar), flxmst(ncelet)
double precision viscf(nfac), viscb(nfabor)

! Local variables

logical          lprev
character(len=80) :: chaine, fname
integer          ivar
integer          ii, ifac , iel, isou
integer          iprev , inc   , iccocg, iiun, ibcl
integer          ivarsc
integer          iiscav
integer          ifcvsl, iflmas, iflmab, f_oi_id
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imucpp, idftnp, iswdyp
integer          iflid , f_id, st_prv_id, st_id,  keydri, iscdri
integer          f_id_al
integer          icvflb, f_dim, iflwgr
integer          icla
integer          icrom_scal
integer          key_buoyant_id, is_buoyant_fld
integer          key_t_ext_id
integer          iviext

integer          ivoid(1)

double precision epsrgp, climgp, extrap, relaxp, blencp, epsilp
double precision epsrsp
double precision rhovst, xk    , xe    , sclnor
double precision thetv , thets , thetap, thetp1
double precision smbexp, dvar, cprovol, prod
double precision temp, idifftp
double precision turb_schmidt
double precision xR, prdtl, alpha_theta
double precision normp
double precision l2norm, l2errork

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1, smbrs, rovsdt
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: coefa_p, coefb_p
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:) :: xcpp
double precision, allocatable, dimension(:) :: srcmas
double precision, allocatable, dimension(:) :: srccond
double precision, allocatable, dimension(:) :: srcmst

double precision, dimension(:,:), pointer :: xut, visten
double precision, dimension(:,:), pointer :: vistet
double precision, dimension(:,:), pointer :: cpro_wgrec_v
double precision, dimension(:), pointer :: cpro_wgrec_s
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom, croma, pcrom
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_omg
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:,:), pointer :: cvara_rij
double precision, dimension(:), pointer :: cvar_al
double precision, dimension(:), pointer :: visct, viscl, cpro_cp, cproa_scal_st
double precision, dimension(:), pointer :: cpro_scal_st
double precision, dimension(:), pointer :: cpro_viscls, cpro_visct
double precision, dimension(:), pointer :: cpro_tsscal
double precision, dimension(:), pointer :: cpro_x2icla
double precision, dimension(:), pointer :: cvar_var, cvara_var, cvara_varsca
double precision, dimension(:), pointer :: cvark_var
double precision, allocatable, dimension(:), target :: wcvark_var
double precision, allocatable, dimension(:) :: errork
double precision, allocatable, dimension(:) :: divflu
! Darcy arrays
double precision, allocatable, dimension(:) :: diverg
double precision, dimension(:), pointer :: cpro_delay, cpro_sat
double precision, dimension(:), pointer :: cproa_delay, cproa_sat
! Radiat arrays
double precision, dimension(:), pointer :: cpro_tsre1, cpro_tsre, cpro_tsri1
character(len=80) :: f_name

type(var_cal_opt) :: vcopt, vcopt_varsc
type(gwf_soilwater_partition) :: sorption_scal

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! --- Variable number
ivar   = isca(iscal)

call field_get_val_s(ivarfl(ivar), cvar_var)
call field_get_val_prev_s(ivarfl(ivar), cvara_var)

! Index of the field
iflid = ivarfl(ivar)

! Key id for buoyant field (inside the Navier Stokes loop)
call field_get_key_id("is_buoyant", key_buoyant_id)
call field_get_key_int(iflid, key_buoyant_id, is_buoyant_fld)

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

! If the scalar is buoyant, it is inside the Navier Stokes loop, and so iterns >=1
! otherwise it is outside of the loop and iterns = -1.
if (  (is_buoyant_fld.eq. 1 .and. iterns.eq.-1) &
  .or.(is_buoyant_fld.eq. 0 .and. iterns.ne.-1)) return

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

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

if (vcopt%iwgrec.eq.1) then
  ! Id weighting field for gradient
  call field_get_key_int(iflid, kwgrec, iflwgr)
  call field_get_dim(iflwgr, f_dim)
  if (f_dim.gt.1) then
    call field_get_val_v(iflwgr, cpro_wgrec_v)
  else
    call field_get_val_s(iflwgr, cpro_wgrec_s)
  endif
endif

! Allocate temporary arrays
allocate(w1(ncelet))
allocate(dpvar(ncelet))
allocate(smbrs(ncelet), rovsdt(ncelet))

if (ippmod(idarcy).eq.1) then
  allocate(diverg(ncelet))
endif

! Initialize variables to avoid compiler warnings

xe = 0.d0
xk = 0.d0

! --- Numero du scalaire eventuel associe dans le cas fluctuation
!         et numero de variable de calcul
iiscav = iscavr(iscal)
if (iiscav.gt.0.and.iiscav.le.nscal) then
  ivarsc = isca(iiscav)
else
  ivarsc = 0
endif

! --- Numero des grandeurs physiques
call field_get_key_int (ivarfl(isca(iscal)), kromsl, icrom_scal)

if (icrom_scal.eq.-1) then
  icrom_scal = icrom
endif

call field_get_val_s(icrom_scal, crom)
call field_have_previous(icrom_scal, lprev)
if (lprev) then
  call field_get_val_prev_s(icrom_scal, croma)
endif
call field_get_val_s(ivisct, visct)
call field_get_val_s(iviscl, viscl)

call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_viscls)
endif

if (idilat.ge.4) then
  call field_get_val_s(iustdy(iscal), cpro_tsscal)
endif

! --- Numero de propriété du terme source si extrapolation
call field_get_key_int(iflid, kstprv, st_prv_id)
if (st_prv_id .ge.0) then
  call field_get_val_s(st_prv_id, cproa_scal_st)
else
  cproa_scal_st => null()
endif

! S pour Source, V pour Variable
thets  = thetss(iscal)
thetv  = vcopt%thetav

call field_get_name(ivarfl(ivar), chaine)

if(vcopt%iwarni.ge.1) then
  write(nfecra,1000) chaine(1:16)
endif

! When solving the Temperature, we solve:
!  rho*cp*Vol*dT/dt + ...

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
  if (icp.ge.0) then
    call field_get_val_s(icp, cpro_cp)
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

! Retrieve turbulent Schmidt value for current scalar
call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)

!===============================================================================
! 2. Source terms
!===============================================================================

! --> Initialization

do iel = 1, ncel
  rovsdt(iel) = 0.d0
  smbrs(iel) = 0.d0
enddo

if (iscal.ne.iscalt) then
  call uitssc(ippmod(idarcy), iflid, cvar_var, smbrs, rovsdt)
else
  call uitsth(iflid, cvar_var, smbrs, rovsdt)
endif

call ustssc &
( nvar   , nscal  , ncepdp , ncesmp ,                            &
  iscal  ,                                                       &
  icepdc , icetsm , itypsm ,                                     &
  dt     ,                                                       &
  ckupdc , smacel , smbrs  , rovsdt )

! C version
call user_source_terms(ivarfl(isca(iscal)), smbrs, rovsdt)

! Take into account radioactive decay rate (implicit source term)
if (ippmod(idarcy).eq.1) then
  call cs_gwf_decay_rate(ivarfl(ivar), rovsdt)
endif

! Store the source terms for convective limiter or time extrapolation for buoyant scalar
call field_get_key_int(iflid, kst, st_id)
if (st_id .ge.0) then
  call field_get_val_s(st_id, cpro_scal_st)

  do iel = 1, ncel
    !Fill the scalar source term field
    cpro_scal_st(iel) = smbrs(iel)
  end do
  ! Handle parallelism and periodicity
  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(cpro_scal_st)
  endif
end if

if (vcopt%ibdtso.gt.1.and.ntcabs.gt.ntinit &
    .and.(idtvar.eq.0.or.idtvar.eq.1)) then
  ! TODO: remove test on ntcabs and implement a "proper" condition for
  ! initialization.
  f_id = ivarfl(ivar)
  call cs_backward_differentiation_in_time(f_id, smbrs, rovsdt)
endif
! Skip first time step after restart if previous values have not been read.
if (vcopt%ibdtso.lt.0) vcopt%ibdtso = iabs(vcopt%ibdtso)

! Set ibdtso value
call field_set_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

! Nudging towards optimal interpolation for current scalar
if (ippmod(iatmos).ge.0) then
  call field_get_key_int(ivarfl(ivar), kopint, f_oi_id)
  if (f_oi_id.ge.0) then
    call cs_at_data_assim_source_term(ivarfl(ivar), smbrs, rovsdt)
  endif
endif

! Atmospheric chemistry
! In case of a semi-coupled resolution, computation of the explicit
! chemical source term to be considered during dynamical resolution
if ((ichemistry.ge.1) .and. (isepchemistry.eq.2) .and. (ntcabs.gt.1)) then
  if ((isca_chem(1).le.iscal).and.(iscal.le.isca_chem(nespg))) then
    call chem_source_terms(iscal, smbrs, rovsdt)
  endif
endif

! Precipitation/dissolution for lagrangian module
! Calculation of source terms du to precipitation and dissolution phenomena
if (ipreci.eq.1.and.iscal.eq.1) then
  call precst(dtref, crom, cvar_var, smbrs)
endif

! Si on extrapole les TS :
!   SMBRS recoit -theta TS du pas de temps precedent
!     (on aurait pu le faire avant ustssc, mais avec le risque que
!      l'utilisateur l'ecrase)
!   SMBRS recoit la partie du terme source qui depend de la variable
!   A l'ordre 2, on suppose que le ROVSDT fourni par l'utilisateur est <0
!     on implicite le terme (donc ROVSDT*RTPA va dans SMBRS)
!   En std, on adapte le traitement au signe de ROVSDT, mais ROVSDT*RTPA va
!     quand meme dans SMBRS (pas d'autre choix)
if (st_prv_id .ge. 0) then
  do iel = 1, ncel
    smbexp = cproa_scal_st(iel)
    ! If the scalar is not buoyant no need of saving the current source term,
    ! save directly the previous one
    if (iterns.eq.-1) then
      cproa_scal_st(iel) = smbrs(iel)
    endif
    ! Terme source du pas de temps precedent et
    ! On suppose -ROVSDT > 0 : on implicite
    !    le terme source utilisateur (le reste)
    smbrs(iel) = rovsdt(iel)*cvara_var(iel) - thets * smbexp
    ! Diagonale
    rovsdt(iel) = - thetv*rovsdt(iel)
  enddo

! Si on n'extrapole pas les TS :
else
  do iel = 1, ncel
    ! Terme source utilisateur
    smbrs(iel) = smbrs(iel) + rovsdt(iel)*cvara_var(iel)
    ! Diagonale
    rovsdt(iel) = max(-rovsdt(iel),zero)
  enddo
endif

! Add thermodynamic pressure variation for the low-Mach algorithm:
! NB: iscalt is the Enthalpy
if ((idilat.eq.3.or.ipthrm.eq.1).and.iscal.eq.iscalt) then
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + (pther - pthera)/dt(iel)*cell_f_vol(iel)
  enddo
endif

! --> Couplage volumique avec Syrthes
!     Ordre 2 non pris en compte

if (iscal.eq.iscalt) then
  call cptssy(iscal, smbrs, rovsdt)
endif

! --> Physique particulieres
!     Ordre 2 non pris en compte

if (ippmod(iphpar).ge.1) then
  call pptssc(iscal, smbrs, rovsdt, tslagr)
endif

! --> Rayonnement
!     Ordre 2 non pris en compte

if (iirayo.ge.1) then

  if (iscal.eq.iscalt) then
    call cs_rad_transfer_source_terms(smbrs, rovsdt)

    ! Store the explicit radiative source term
    if (idilat.ge.4) then
      call field_get_id("rad_st", f_id)
      call field_get_val_s(f_id,cpro_tsre1)
      do iel = 1, ncel
        cpro_tsscal(iel) = cpro_tsscal(iel)   &
        + cpro_tsre1(iel)*cell_f_vol(iel)
      enddo
    endif
  endif

  !-> Charbon pulverise
  !   Ordre 2 non pris en compte
  ! new model
  if (ippmod(iccoal) .ge. 0) then
    if (isca(iscal).ge.isca(ih2(1)) .and.       &
        isca(iscal).le.isca(ih2(nclacp))) then

      call cs_coal_radst &
      !=================
      ( ivar   , ncelet , ncel  ,               &
        cell_f_vol , smbrs , rovsdt )

    endif

    if (iscal .eq.ihgas) then
      call field_get_id("rad_st", f_id)
      call field_get_val_s(f_id, cpro_tsre1)
      do iel = 1, ncel
        smbrs(iel) = smbrs(iel)+volume(iel)*cpro_tsre1(iel)
      enddo

      do icla = 1, nclacp
        write(f_name,  '("rad_st_", i2.2)') icla+1
        call field_get_id(f_name, f_id)
        call field_get_val_s(f_id, cpro_tsre)
        call field_get_val_s(ix2(icla), cpro_x2icla)
        do iel = 1, ncel
          smbrs(iel) = smbrs(iel)-volume(iel)*cpro_tsre(iel) &
                                             *cpro_x2icla(iel)
        enddo
      enddo
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
     ( ivar   , ncelet , ncel  ,                &
       cell_f_vol , smbrs , rovsdt)

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

  ! On incremente SMBRS par -Gamma RTPA et ROVSDT par Gamma
  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   , isso2t(iscal) ,            &
   icetsm , itypsm(1,ivar) ,                                      &
   cell_f_vol , cvara_var    , smacel(1,ivar) , srcmas   ,        &
   smbrs  , rovsdt , w1)

  deallocate(srcmas)

  ! Si on extrapole les TS on met Gamma Pinj dans cproa_scal_st
  if (st_prv_id .ge. 0) then
    do iel = 1, ncel
      cproa_scal_st(iel) = cproa_scal_st(iel) + w1(iel)
    enddo
  ! Sinon on le met directement dans SMBRS
  else
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) + w1(iel)
    enddo
  endif

endif

! Condensation source terms for the scalars
! associated to a surface zone (icondb=0)
if (nfbpcd.gt.0) then

  allocate(srccond(nfbpcd))

  ! When treating the Temperature, the equation is multiplied by Cp
  do ii = 1, nfbpcd
    ifac= ifbpcd(ii)
    iel = ifabor(ifac)

    if (spcond(ii,ipr).lt.0.d0 .and.itypcd(ii,ivar).eq.1) then
      srccond(ii) = spcond(ii,ipr)*xcpp(iel)
    else
      srccond(ii) = 0.d0
    endif
  enddo

  call condensation_source_terms &
  !=============================
  (ncelet , ncel ,                                       &
   iscal  ,                                              &
   nfbpcd , ifbpcd  , itypcd(1,ivar) ,                   &
   0      , ivoid   , ivoid          ,                   &
   spcond(1,ivar)   , srccond        ,                   &
   rvoid  , rvoid   , rvoid          ,                   &
   cvara_var        ,                                    &
   smbrs            , rovsdt )

  deallocate(srccond)

endif

! Condensation source terms for the scalars
! associated to a volumic zone (icondv=0)
! taking into account the metal mass
! structures condensation modelling
if (icondv.eq.0) then
  allocate(srcmst(ncelet))

  ! When treating the Temperature, the equation is multiplied by Cp
  do ii = 1, ncmast
    iel = ltmast(ii)

    if (svcond(iel,ipr).lt.0.d0 .and.itypst(ii,ivar).eq.1) then
      srcmst(iel) = svcond(iel,ipr)*xcpp(iel)
    else
      srcmst(iel) = 0.d0
    endif
  enddo

  call condensation_source_terms &
  !=============================
  (ncelet , ncel ,                                       &
   iscal  ,                                              &
   0      , ivoid   , ivoid          ,                   &
   ncmast , ltmast  , itypst(1,ivar) ,                   &
   rvoid  , rvoid   ,                                    &
   svcond(1,ivar)   , srcmst         , flxmst  ,         &
   cvara_var        ,                                    &
   smbrs            , rovsdt )

  deallocate(srcmst)

endif


! If the current scalar is the variance of an other scalar,
! production and dissipation terms are added.
if (itspdv.eq.1) then

  if (itytur.eq.2 .or. itytur.eq.3 .or. itytur.eq.5 .or. iturb.eq.60) then

    ! Allocate a temporary array for the gradient reconstruction
    allocate(grad(3,ncelet))
    allocate(coefa_p(nfabor), coefb_p(nfabor))

    ! Remarque : on a prevu la possibilite de scalaire associe non
    !  variable de calcul, mais des adaptations sont requises

    if (ivarsc.le.0) then
      write(nfecra,9000)ivarsc
      call csexit(1)
    endif

    iprev = 1
    inc = 1
    iccocg = 1

    ! Homogeneous Neumann on convective inlet on the production term for the
    ! variance
    call field_get_val_prev_s(ivarfl(ivarsc), cvara_varsca)
    call field_get_coefa_s (ivarfl(ivarsc), coefap)
    call field_get_coefb_s (ivarfl(ivarsc), coefbp)

    ! pas de diffusion en entree
    do ifac = 1, nfabor
      coefa_p(ifac) = coefap(ifac)
      coefb_p(ifac) = coefbp(ifac)
      if (itypfb(ifac).eq.i_convective_inlet) then
        coefa_p(ifac) = 0.d0
        coefb_p(ifac) = 1.d0
      endif
    enddo

    call field_get_key_struct_var_cal_opt(ivarfl(ivarsc), vcopt_varsc)

    nswrgp = vcopt_varsc%nswrgr
    imligp = vcopt_varsc%imligr
    iwarnp = vcopt_varsc%iwarni
    epsrgp = vcopt_varsc%epsrgr
    climgp = vcopt_varsc%climgr
    extrap = vcopt_varsc%extrag

    call gradient_s                                                          &
     ( ivarfl(ivarsc)  , imrgra , inc    , iccocg , nswrgp , imligp ,        &
       iwarnp          , epsrgp , climgp , extrap ,                          &
       cvara_varsca    , coefa_p, coefb_p,                                   &
       grad )

    deallocate (coefa_p, coefb_p)

    ! Production Term

    ! NB: diffusivity is clipped to 0 because in LES, it might be negative. Problematic
    ! for the variance even if it is strange to use variance and LES ...

    ! Time extrapolation (2nd order)
    if (st_prv_id .ge. 0) then

      ! Not extapolated value of the viscosity
      call field_get_val_s(ivisct, cpro_visct)
      call field_get_key_int(ivisct, key_t_ext_id, iviext)
      if (iviext.gt.0) call field_get_val_prev_s(ivisct, cpro_visct)

      ! iscal is the variance of the scalar iiscav
      ! with modelized turbulent fluxes GGDH or AFM or DFM
      if (ityturt(iiscav).ge.1) then

        ! Name of the scalar iiscav associated to the variance iscal
        call field_get_name(ivarfl(ivarsc), fname)

        ! Index of the corresponding turbulent flux
        call field_get_id(trim(fname)//'_turbulent_flux', f_id)

        call field_get_val_v(f_id, xut)

        do iel = 1, ncel
          cproa_scal_st(iel) = cproa_scal_st(iel)                              &
                             - 2.d0*xcpp(iel)*cell_f_vol(iel)*crom(iel)        &
                                          *(xut(1,iel)*grad(1,iel)             &
                                           +xut(2,iel)*grad(2,iel)             &
                                           +xut(3,iel)*grad(3,iel) )
        enddo
      ! SGDH model
      else
        do iel = 1, ncel
          cproa_scal_st(iel) = cproa_scal_st(iel)                             &
               + 2.d0*xcpp(iel)*max(cpro_visct(iel),zero)                     &
               *cell_f_vol(iel)/turb_schmidt                                 &
               *(grad(1,iel)**2 + grad(2,iel)**2 + grad(3,iel)**2)
        enddo
      endif

    ! If not time extrapolation...
    else

      call field_get_val_s(ivisct, cpro_visct)

      ! iscal is the variance of the scalar iiscav
      ! with modelized turbulent fluxes GGDH or AFM or DFM
      if (ityturt(iiscav).ge.1) then


        ! Name of the scalar ivarsc associated to the variance iscal
        call field_get_name(ivarfl(ivarsc), fname)

        ! Index of the corresponding turbulent flux
        call field_get_id(trim(fname)//'_turbulent_flux', f_id)

        call field_get_val_v(f_id, xut)

        do iel = 1, ncel
          cprovol = xcpp(iel)*cell_f_vol(iel)*crom(iel)
          ! Special time stepping to ensure positivity of the variance
          prod = -2.d0 * (xut(1,iel)*grad(1,iel)      &
                         +xut(2,iel)*grad(2,iel)      &
                         +xut(3,iel)*grad(3,iel) )

          smbrs(iel) = smbrs(iel) + max(prod * cprovol, 0.d0)

          ! Implicit "production" term when negative, but check if the
          ! variance is non-zero.
          if (cvar_var(iel).gt. epzero * abs(prod * dt(iel)) &
            .and.prod*cprovol.lt.0.d0) then
            rovsdt(iel) = rovsdt(iel)  &
                        - prod * cprovol / cvar_var(iel)
          endif
        enddo

      ! SGDH model
      else
        do iel = 1, ncel
          smbrs(iel) = smbrs(iel)                                            &
                     + 2.d0*xcpp(iel)*max(cpro_visct(iel),zero)           &
                     * cell_f_vol(iel)/turb_schmidt                       &
                     * (grad(1,iel)**2 + grad(2,iel)**2 + grad(3,iel)**2)
        enddo
      endif

      ! Production term for a variance  TODO compute ustdy when isso2t >0
      if (idilat.ge.4) then
        do iel = 1, ncel
          cpro_tsscal(iel) = cpro_tsscal(iel) +                   &
               2.d0*xcpp(iel)*max(cpro_visct(iel),zero)        &
             *cell_f_vol(iel)/turb_schmidt                     &
             *(grad(1,iel)**2 + grad(2,iel)**2 + grad(3,iel)**2)
        enddo
      endif
    endif

    ! Free memory
    deallocate(grad)

    ! Dissipation term
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
      if (irijco.eq.1) then
        call field_get_val_prev_v(ivarfl(irij), cvara_rij)
      else
        call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
        call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
        call field_get_val_prev_s(ivarfl(ir33), cvara_r33)
      endif
      ! EB- AFM or EB-DFM or EB-GGDH
      if (iturt(iiscav).eq.11 .or. iturt(iiscav).eq.21 .or. iturt(iiscav).eq.31) then
        ! Name of the scalar corresponding to the current variance
        call field_get_name(ivarfl(ivarsc), fname)
        ! Index of the corresponding alpha
        call field_get_id(trim(fname)//'_alpha', f_id_al)
        call field_get_val_s(f_id_al, cvar_al)
      endif
    elseif(iturb.eq.60) then
      call field_get_val_prev_s(ivarfl(ik), cvara_k)
      call field_get_val_prev_s(ivarfl(iomg), cvara_omg)
    endif

    do iel = 1, ncel
      if (itytur.eq.2 .or. itytur.eq.5) then
        xk = cvara_k(iel)
        xe = cvara_ep(iel)
      elseif (itytur.eq.3) then
        if (irijco.eq.1) then
          xk = 0.5d0*(cvara_rij(1,iel)+cvara_rij(2,iel)+cvara_rij(3,iel))
          xe = cvara_ep(iel)
        else
          xk = 0.5d0*(cvara_r11(iel)+cvara_r22(iel)+cvara_r33(iel))
          xe = cvara_ep(iel)
        endif
      elseif(iturb.eq.60) then
        xk = cvara_k(iel)
        xe = cmu*xk*cvara_omg(iel)
      endif

      ! Implicit term -1/Rh * Eps/k * Variance
      !  with
      ! Rh = R = 0.8 for SGDH
      ! Rh = (1-alpha_T) * Pr + R * alpha_T
      ! with - R = 0.5
      !      - alpha_T = 1.0 for GGDH/DFM/AFM
      if(iturt(iiscav).eq.11.or.iturt(iiscav).eq.21.or.iturt(iiscav).eq.31) then
        alpha_theta = cvar_al(iel)
      else
        alpha_theta = 1.d0
      endif

      if (ifcvsl.ge.0) then
        prdtl = viscl(iel)*xcpp(iel)/cpro_viscls(iel)
      else
        prdtl = viscl(iel)*xcpp(iel)/visls0(iscal)
      endif
      xR = ( 1.d0 - alpha_theta ) * prdtl + alpha_theta * rvarfl(iscal)

      rhovst = xcpp(iel)*crom(iel)*xe/(xk * xR)       &
             *cell_f_vol(iel)

      ! La diagonale recoit eps/Rk, (*theta eventuellement)
      rovsdt(iel) = rovsdt(iel) + rhovst*thetap
      ! SMBRS recoit la dissipation
      smbrs(iel) = smbrs(iel) - rhovst*cvara_var(iel)
    enddo

  endif

endif

if (st_prv_id .ge. 0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + thetp1 * cproa_scal_st(iel)
  enddo
endif

! Compressible algorithm
! or Low Mach compressible algos with mass flux prediction
if (ippmod(icompf).ge.0.or.(idilat.gt.1.and.ipredfl.eq.1.and.irovar.eq.1)) then
  pcrom => croma

! Low Mach compressible algos (conservative in time).
! Same algo. for Volume of Fluid method
else if ((idilat.gt.1.or.ivofmt.gt.0).and.irovar.eq.1) then
  if (iterns.eq.1) then
    call field_get_val_prev2_s(icrom_scal, pcrom)
  else
    call field_get_val_prev_s(icrom_scal, pcrom)
  endif

! Deprecated algo or constant density
else
  call field_get_val_s(icrom_scal, pcrom)
endif

! "VITESSE" DE DIFFUSION FACETTE

! On prend le MAX(mu_t,0) car en LES dynamique mu_t peut etre negatif
! (clipping sur (mu + mu_t)). On aurait pu prendre
! MAX(K + K_t,0) mais cela autoriserait des K_t negatif, ce qui est
! considere ici comme non physique.
if (vcopt%idiff.ge.1) then

  allocate(vistet(6,ncelet))

  ! AFM model or DFM models: add div(Cp*rho*T'u') to smbrs
  ! Compute T'u' for GGDH
  if (ityturt(iscal).ge.1) then

    ! EB-GGDH/AFM/DFM: solving alpha for the scalar
    if (iturt(iscal).eq.11.or.iturt(iscal).eq.21.or.iturt(iscal).eq.31) then
      ! Name of the scalar
      call field_get_name(ivarfl(isca(iscal)), fname)

      ! Index of the corresponding turbulent flux
      call field_get_id(trim(fname)//'_alpha', f_id)

      call resalp(f_id, xclt)
    endif

    call divrit &
    !==========
    ( nscal  ,                                                       &
      iscal  ,                                                       &
      dt     ,                                                       &
      xcpp   ,                                                       &
      vistet ,                                                       &
      smbrs  )

  endif

  ! Scalar diffusivity
  if (iand(vcopt%idften, ISOTROPIC_DIFFUSION).ne.0) then

    idifftp = vcopt%idifft
    if (ityturt(iscal).eq.3) then
      idifftp = 0
    endif

    if (ifcvsl.lt.0) then
      do iel = 1, ncel
        w1(iel) = visls0(iscal)                                     &
           + idifftp*xcpp(iel)*max(visct(iel),zero)/turb_schmidt
      enddo
    else
      do iel = 1, ncel
        w1(iel) = cpro_viscls(iel)                                &
           + idifftp*xcpp(iel)*max(visct(iel),zero)/turb_schmidt
      enddo
    endif

    if (vcopt%iwgrec.eq.1) then
      ! Weighting for gradient
      do iel = 1, ncel
        cpro_wgrec_s(iel) = w1(iel)
      enddo
      call synsca(cpro_wgrec_s)
      if (irangp.ge.0.or.iperio.eq.1) then
         call synsca(cpro_wgrec_s)
      endif
    endif

    call viscfa &
    !==========
   ( imvisf ,                      &
     w1     ,                      &
     viscf  , viscb  )

  ! Symmetric tensor diffusivity (GGDH)
  else if (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then

    ! Allocate temporary arrays
    allocate(viscce(6,ncelet))
    allocate(weighf(2,nfac))
    allocate(weighb(nfabor))

    if (iturb.ne.32) then
      call field_get_val_v(ivsten, visten)
    else ! EBRSM and (GGDH or AFM)
      call field_get_val_v(ivstes, visten)
    endif

    if (iturt(iscal).eq.11.or.iturt(iscal).eq.20.or.iturt(iscal).eq.21) then
      if (ifcvsl.lt.0) then
        do iel = 1, ncel

          temp = vcopt%idifft*xcpp(iel)
          viscce(1,iel) = temp*vistet(1,iel) + visls0(iscal)
          viscce(2,iel) = temp*vistet(2,iel) + visls0(iscal)
          viscce(3,iel) = temp*vistet(3,iel) + visls0(iscal)
          viscce(4,iel) = temp*vistet(4,iel)
          viscce(5,iel) = temp*vistet(5,iel)
          viscce(6,iel) = temp*vistet(6,iel)

        enddo
      else
        do iel = 1, ncel

          temp = vcopt%idifft*xcpp(iel)
          viscce(1,iel) = temp*vistet(1,iel) + cpro_viscls(iel)
          viscce(2,iel) = temp*vistet(2,iel) + cpro_viscls(iel)
          viscce(3,iel) = temp*vistet(3,iel) + cpro_viscls(iel)
          viscce(4,iel) = temp*vistet(4,iel)
          viscce(5,iel) = temp*vistet(5,iel)
          viscce(6,iel) = temp*vistet(6,iel)

        enddo
      endif
    else
      if (ifcvsl.lt.0) then
        do iel = 1, ncel

          temp = vcopt%idifft*xcpp(iel)*ctheta(iscal)/csrij
          viscce(1,iel) = temp*visten(1,iel) + visls0(iscal)
          viscce(2,iel) = temp*visten(2,iel) + visls0(iscal)
          viscce(3,iel) = temp*visten(3,iel) + visls0(iscal)
          viscce(4,iel) = temp*visten(4,iel)
          viscce(5,iel) = temp*visten(5,iel)
          viscce(6,iel) = temp*visten(6,iel)

        enddo
      else
        do iel = 1, ncel

          temp = vcopt%idifft*xcpp(iel)*ctheta(iscal)/csrij
          viscce(1,iel) = temp*visten(1,iel) + cpro_viscls(iel)
          viscce(2,iel) = temp*visten(2,iel) + cpro_viscls(iel)
          viscce(3,iel) = temp*visten(3,iel) + cpro_viscls(iel)
          viscce(4,iel) = temp*visten(4,iel)
          viscce(5,iel) = temp*visten(5,iel)
          viscce(6,iel) = temp*visten(6,iel)

        enddo
      endif
    endif

    iwarnp = vcopt%iwarni

    if (vcopt%iwgrec.eq.1) then
      ! Weighting for gradient
      do iel = 1, ncel
        do isou = 1, 6
          cpro_wgrec_v(isou,iel) = viscce(isou,iel)
        enddo
      enddo
      call syntis(cpro_wgrec_v)
    endif

    call vitens &
    !==========
   ( viscce , iwarnp ,             &
     weighf , weighb ,             &
     viscf  , viscb  )

  endif

  deallocate(vistet)

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

! Not Darcy
if (ippmod(idarcy).eq.-1) then

  ! --> Unsteady term and mass aggregation term
  do iel = 1, ncel
    rovsdt(iel) = rovsdt(iel)                                                 &
                + vcopt%istat*xcpp(iel)*pcrom(iel)*cell_f_vol(iel)/dt(iel)
  enddo

! Darcy : we take into account the porosity and delay for underground transport
else
  ! Retrieve sorption options for current scalar for ground water flow module
  call field_get_key_struct_gwf_soilwater_partition(ivarfl(ivar), sorption_scal)
  call field_get_val_s(sorption_scal%idel, cpro_delay)
  call field_get_val_prev_s(sorption_scal%idel, cproa_delay)
  call field_get_val_prev_s_by_name('saturation', cproa_sat)
  call field_get_val_s_by_name('saturation', cpro_sat)

  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)*cpro_delay(iel)*cpro_sat(iel)
    rovsdt(iel) = (rovsdt(iel) + vcopt%istat*xcpp(iel)*pcrom(iel) &
                                *volume(iel)/dt(iel))             &
                                *cpro_delay(iel)*cpro_sat(iel)
  enddo

  ! treatment of kinetic sorption
  if (sorption_scal%kinetic.eq.1) then
    call cs_gwf_kinetic_reaction(ivarfl(ivar), rovsdt, smbrs)
  endif

endif

! Scalar with a Drift:
! compute the convective flux
!----------------------------

call field_get_key_int(iflid, keydri, iscdri)

if (iscdri.ge.1) then
  allocate(divflu(ncelet))

  call driflu &
  !=========
  ( iflid  ,                                                       &
    dt     ,                                                       &
    imasfl , bmasfl ,                                              &
    divflu )

  iconvp = vcopt%iconv
  thetap = vcopt%thetav

  ! NB: if the porosity module is swiched on, the the porosity is already
  ! taken into account in divflu

  ! --> mass aggregation term
  do iel = 1, ncel
    rovsdt(iel) = rovsdt(iel) + iconvp*thetap*divflu(iel)
    smbrs(iel) = smbrs(iel) - iconvp*divflu(iel)*cvara_var(iel)
  enddo

  deallocate(divflu)
endif

! Darcy:
! This step is necessary because the divergence of the
! hydraulic head (or pressure), which is taken into account as the
! 'aggregation term' via codits -> bilsca -> bilsc2,
! is not exactly equal to the loss of mass of water in the unsteady case
! (just as far as the precision of the Newton scheme is good),
! and does not take into account the sorption of the tracer
! (represented by the 'delay' coefficient).
! We choose to 'cancel' the aggregation term here, and to add to smbr
! (right hand side of the transport equation)
! the necessary correction to take sorption into account and get the exact
! conservation of the mass of tracer.

if (ippmod(idarcy).eq.1) then
  if (darcy_unsteady.eq.1) then
    call divmas(1, imasfl , bmasfl , diverg)
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) - diverg(iel)*cvar_var(iel)                   &
        + cell_f_vol(iel)/dt(iel)*cvar_var(iel)                             &
        *( cproa_delay(iel)*cproa_sat(iel) - cpro_delay(iel)*cpro_sat(iel))
    enddo
    do iel = 1, ncel
      rovsdt(iel) = rovsdt(iel) + thetv*diverg(iel)
    enddo
  endif
endif

!===============================================================================
! 3. Solving
!===============================================================================

if (iterns.ge.1) then
  allocate(wcvark_var(ncelet))
  cvark_var => wcvark_var
  do iel = 1, ncelet
    cvark_var(iel) = cvar_var(iel)
  enddo
else
  call field_get_val_s(ivarfl(ivar), cvark_var)
endif

iconvp = vcopt%iconv
idiffp = vcopt%idiff
idftnp = vcopt%idften
ndircp = vcopt%ndircl
nswrsp = vcopt%nswrsm
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
ircflp = vcopt%ircflu
ischcp = vcopt%ischcv
isstpp = vcopt%isstpc
iescap = 0
iswdyp = vcopt%iswdyn
iwarnp = vcopt%iwarni
blencp = vcopt%blencv
epsilp = vcopt%epsilo
epsrsp = vcopt%epsrsm
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
extrap = vcopt%extrag
relaxp = vcopt%relaxv
! all boundary convective flux with upwind
icvflb = 0
normp = -1.d0

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

call codits &
!==========
 ( idtvar , iterns , ivarfl(ivar)    , iconvp , idiffp , ndircp , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp , normp  ,                                              &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   cvara_var       , cvark_var       ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , viscce ,                   &
   weighf , weighb ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbrs  , cvar_var        , dpvar  ,                   &
   xcpp   , rvoid  )

!===============================================================================
! 4. clipping and log
!===============================================================================

call clpsca(iscal)

! Update solid phase concentration for kinetic or precipitation models
if (ippmod(idarcy).eq.1) then
  ! Update of sorbed concentration
  if (sorption_scal%kinetic.eq.1) then
    call cs_gwf_sorbed_concentration_update(ivarfl(ivar))
  endif
  ! Treatment of precipitation for groundwater flow module.
  if (sorption_scal%imxsol.ge.0) then
    call cs_gwf_precipitation(ivarfl(ivar))
  endif
endif

if (idilat.ge.4.and.itspdv.eq.1) then

  do iel = 1, ncel
    if (itytur.eq.2 .or. itytur.eq.5) then
      xk = cvara_k(iel)
      xe = cvara_ep(iel)
    elseif (itytur.eq.3) then
        if (irijco.eq.1) then
          xk = 0.5d0*(cvara_rij(1,iel)+cvara_rij(2,iel)+cvara_rij(3,iel))
          xe = cvara_ep(iel)
        else
          xk = 0.5d0*(cvara_r11(iel)+cvara_r22(iel)+cvara_r33(iel))
          xe = cvara_ep(iel)
        endif
    elseif(iturb.eq.60) then
      xk = cvara_k(iel)
      xe = cmu*xk*cvara_omg(iel)
    endif
    rhovst = xcpp(iel)*crom(iel)*xe/(xk * rvarfl(iscal))       &
           *cell_f_vol(iel)

    cpro_tsscal(iel) = cpro_tsscal(iel) - rhovst*cvar_var(iel)

  enddo

endif

! Store the implicit part of the radiative source term
if (idilat.ge.4.and.iirayo.ge.1.and.iscal.eq.iscalt) then
  call field_get_id("rad_st_implicit", f_id)
  call field_get_val_s(f_id,cpro_tsri1)
  do iel = 1, ncel
    dvar = cvar_var(iel)-cvara_var(iel)
    cpro_tsscal(iel) = cpro_tsscal(iel) &
                     - cpro_tsri1(iel)*dvar*cell_f_vol(iel)
  enddo
endif

! BILAN EXPLICITE (VOIR CODITS : ON ENLEVE L'INCREMENT)
! Ceci devrait etre valable avec le theta schema sur les Termes source

if (vcopt%iwarni.ge.2) then
  if (vcopt%nswrsm.gt.1) then
    ibcl = 1
  else
    ibcl = 0
  endif
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                                 &
            - vcopt%istat*xcpp(iel)*(pcrom(iel)/dt(iel))*cell_f_vol(iel)        &
                *(cvar_var(iel)-cvara_var(iel))*ibcl
  enddo
  sclnor = sqrt(cs_gdot(ncel,smbrs,smbrs))
  write(nfecra,1200)chaine(1:16) ,sclnor
endif

! Log in case of PISO-like sub iterations
if (iterns.ge.1.and.vcopt%iwarni.ge.1) then

  allocate(errork(ncelet))
  do iel = 1, ncel
    errork(iel) = cvar_var(iel) - cvark_var(iel)
  enddo

  l2errork = sqrt(cs_gres(ncel, cell_f_vol, errork, errork))
  deallocate(errork)

  l2norm = sqrt(cs_gres(ncel, cell_f_vol, cvara_var, cvara_var))

  write(nfecra,2601) ivarfl(ivar), iterns, l2errork, l2errork/l2norm, l2norm

endif

! Free memory
if (allocated(wcvark_var)) deallocate(wcvark_var)
deallocate(w1)
deallocate(smbrs, rovsdt)
if (allocated(viscce)) deallocate(viscce)
if (allocated(weighf)) deallocate(weighf, weighb)
deallocate(dpvar)
deallocate(xcpp)
if (allocated(diverg)) deallocate(diverg)

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

2601 format('PISO scalar',I10, 'iter=', I10, 'L2 error = ',E12.4,' L2 normalized error', E12.4, 'L2 nomr', E12.4 ,/)

!----
! End
!----

return

end subroutine
