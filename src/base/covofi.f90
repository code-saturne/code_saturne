!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ncmast        number of cells with condensation source terms
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
!> \param[in]     propce        physical properties at cell centers
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
   iscal  , itspdv ,                                              &
   icepdc , icetsm , ifbpcd , ltmast ,                            &
   itypsm , itypcd , itypst ,                                     &
   dt     , propce , tslagr ,                                     &
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
use darcy_module
use cs_c_bindings
use pointe, only: itypfb

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp , nfbpcd ,  ncmast
integer          iscal  , itspdv

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ifbpcd(nfbpcd), itypcd(nfbpcd,nvar)
integer          ltmast(ncelet), itypst(ncelet,nvar)

double precision dt(ncelet)
double precision propce(ncelet,*)
double precision tslagr(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision spcond(nfbpcd,nvar)
double precision svcond(ncelet,nvar), flxmst(ncelet)
double precision viscf(nfac), viscb(nfabor)

! Local variables

logical          lprev
character(len=80) :: chaine, fname
integer          ivar
integer          ii, ifac , iel, isou
integer          iprev , inc   , iccocg, isqrt, iii, iiun, ibcl
integer          ivarsc
integer          iiscav
integer          ifcvsl, iflmas, iflmab
integer          ipcvso
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imucpp, idftnp, iswdyp
integer          iflid , f_id, st_prv_id,  keydri, iscdri
integer          icvflb, f_dim, iflwgr
integer          delay_id, icla

integer          ivoid(1)

logical          interleaved

double precision epsrgp, climgp, extrap, relaxp, blencp, epsilp
double precision epsrsp
double precision rhovst, xk    , xe    , sclnor
double precision thetv , thets , thetap, thetp1
double precision smbexp, dvar, cprovol, prod
double precision temp, idifftp

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1, smbrs, rovsdt
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, allocatable, dimension(:,:) :: grad, grdni
double precision, allocatable, dimension(:) :: coefa_p, coefb_p
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:) :: xcpp
double precision, allocatable, dimension(:) :: srcmas
double precision, allocatable, dimension(:) :: srccond
double precision, allocatable, dimension(:) :: srcmst

double precision, dimension(:,:), pointer :: xut, visten
double precision, dimension(:,:), pointer :: cpro_wgrec_v
double precision, dimension(:), pointer :: cpro_wgrec_s
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom, croma, pcrom
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: porosi
double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_omg
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: visct, cpro_cp, c_st_scal
double precision, dimension(:), pointer :: cpro_viscls
! Darcy arrays
double precision, allocatable, dimension(:) :: diverg
double precision, dimension(:), pointer :: cpro_delay, cpro_sat
double precision, dimension(:), pointer :: cproa_delay, cproa_sat
double precision, dimension(:), pointer :: cvar_var, cvara_var, cvara_varsca

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

if (iwgrec(ivar).eq.1) then
  ! Id weighting field for gradient
  call field_get_key_int(iflid, kwgrec, iflwgr)
  call field_get_dim(iflid, f_dim, interleaved)
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
call field_get_val_s(icrom, crom)
call field_have_previous(icrom, lprev)
if (lprev) then
  call field_get_val_prev_s(icrom, croma)
endif
call field_get_val_s(iprpfl(ivisct), visct)

call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_viscls)
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
    ( iflid  , cvar_var , smbrs  , rovsdt )
  else
    call uitsth &
    ( iflid  , cvar_var , smbrs  , rovsdt )
  endif
endif

call ustssc &
!==========
( nvar   , nscal  , ncepdp , ncesmp ,                            &
  iscal  ,                                                       &
  icepdc , icetsm , itypsm ,                                     &
  dt     ,                                                       &
  ckupdc , smacel , smbrs  , rovsdt )

if (ibdtso(ivar).gt.ntinit.and.ntcabs.gt.1 &
    .and.(idtvar.eq.0.or.idtvar.eq.1)) then
  ! TODO: remove test on ntcabs and implemente a "proper" condition for
  ! initialization.
  f_id = ivarfl(ivar)
  call cs_backward_differentiation_in_time(f_id, smbrs, rovsdt)
endif
! Skip first time step after restart if previous values have not been read.
if (ibdtso(ivar).lt.0) ibdtso(ivar) = iabs(ibdtso(ivar))

! Atmospheric chemistry
! In case of a semi-coupled resolution, computation of the explicit
! chemical source term to be considered during dynamical resolution
! The first nespg user scalars are supposed to be chemical species
if ((ichemistry.ge.1).and.(isepchemistry.eq.2)                    &
     .and.(iscal.le.nespg).and.(ntcabs.gt.1)) then
  call chem_source_terms(iscal, smbrs, rovsdt)
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
    smbrs(iel) = rovsdt(iel)*cvara_var(iel) - thets*smbexp
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
  call cptssy(iscal, smbrs, rovsdt)
  !==========
endif

! --> Physique particulieres
!     Ordre 2 non pris en compte

if (ippmod(iphpar).ge.1) then
  call pptssc &
  !==========
 ( iscal  ,                                                       &
   propce ,                                                       &
   smbrs  , rovsdt , tslagr )
endif

! --> Rayonnement
!     Ordre 2 non pris en compte

if (iirayo.ge.1) then

  if (iscal.eq.iscalt) then
    call raysca &
    !==========
  ( iscalt,ncelet,ncel,     &
    smbrs, rovsdt,volume)

    ! Store the explicit radiative source term
    if (idilat.ge.4) then
      do iel = 1, ncel
        propce(iel,ipproc(iustdy(iscalt))) = &
        propce(iel,ipproc(iustdy(iscalt)))   &
        + propce(iel,ipproc(itsre(1)))*volume(iel)
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
        volume , propce , smbrs , rovsdt )

    endif

    if (iscal .eq.ihgas) then

      do iel = 1, ncel

        smbrs(iel) = smbrs(iel)+volume(iel)*propce(iel,ipproc(itsre(1)))
        do icla = 1, nclacp
          smbrs(iel) = smbrs(iel)-volume(iel)*propce(iel,ipproc(itsre(icla+1))) &
                                             *propce(iel,ipproc(ix2(icla)))
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
       volume , propce , smbrs , rovsdt)

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
   volume , cvara_var    , smacel(1,ivar) , srcmas   ,            &
   smbrs  , rovsdt , w1)

  deallocate(srcmas)

  ! Si on extrapole les TS on met Gamma Pinj dans c_st_scal
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
    allocate(grad(3,ncelet),grdni(ncelet,3))
    allocate(coefa_p(nfabor), coefb_p(nfabor))

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

    ! Homogeneous Neumann on convective inlet on the production term for the
    ! variance
    call field_get_val_prev_s(ivarfl(iii), cvara_varsca)
    call field_get_coefa_s (ivarfl(iii), coefap)
    call field_get_coefb_s (ivarfl(iii), coefbp)

    ! pas de diffusion en entree
    do ifac = 1, nfabor
      coefa_p(ifac) = coefap(ifac)
      coefb_p(ifac) = coefbp(ifac)
      if (itypfb(ifac).eq.i_convective_inlet) then
        coefa_p(ifac) = 0.d0
        coefb_p(ifac) = 1.d0
      endif
    enddo

    nswrgp = nswrgr(iii)
    imligp = imligr(iii)
    iwarnp = iwarni(iii)
    epsrgp = epsrgr(iii)
    climgp = climgr(iii)
    extrap = extrag(iii)

    call grdcel &
    !==========
     ( iii    , imrgra , inc    , iccocg , nswrgp , imligp ,          &
       iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
       cvara_varsca    , coefa_p, coefb_p,                            &
       grdni   )

    do iel = 1, ncel
      grad(1,iel) = grdni(iel,1)
      grad(2,iel) = grdni(iel,2)
      grad(3,iel) = grdni(iel,3)
    enddo

    deallocate (grdni, coefa_p, coefb_p)

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
          c_st_scal(iel) = c_st_scal(iel)                                  &
                         - 2.d0*xcpp(iel)*volume(iel)*crom(iel)            &
                                          *(xut(1,iel)*grad(1,iel)         &
                                           +xut(2,iel)*grad(2,iel)         &
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
          cprovol = xcpp(iel)*volume(iel)*crom(iel)
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
                     + 2.d0*xcpp(iel)*max(propce(iel,ipcvso),zero)           &
                     * volume(iel)/sigmas(iscal)                             &
                     * (grad(1,iel)**2 + grad(2,iel)**2 + grad(3,iel)**2)
        enddo
      endif

      ! Production term for a variance  TODO compute ustdy when isso2t >0
      if (idilat.ge.4) then
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
      smbrs(iel) = smbrs(iel) - rhovst*cvara_var(iel)
    enddo

  endif

endif

if (st_prv_id .ge. 0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + thetp1 * c_st_scal(iel)
  enddo
endif

! Low Mach compressible algos (conservative in time).
! Same algo for cavitation and compressible algorithm.
if (ippmod(icompf).ge.0.or.idilat.gt.1.or.icavit.ge.0) then
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
    if (ifcvsl.lt.0) then
      do iel = 1, ncel
        w1(iel) = visls0(iscal)                                     &
           + idifftp*xcpp(iel)*max(visct(iel),zero)/sigmas(iscal)
      enddo
    else
      do iel = 1, ncel
        w1(iel) = cpro_viscls(iel)                                &
           + idifftp*xcpp(iel)*max(visct(iel),zero)/sigmas(iscal)
      enddo
    endif

    if (iwgrec(ivar).eq.1) then
      ! Weighting for gradient
      do iel = 1, ncel
        cpro_wgrec_s(iel) = w1(iel)
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

    if (iturb.ne.32) then
      call field_get_val_v(ivsten, visten)
    else ! EBRSM and (GGDH or AFM)
      call field_get_val_v(ivstes, visten)
    endif

    if (ifcvsl.lt.0) then
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
        viscce(1,iel) = temp*visten(1,iel) + cpro_viscls(iel)
        viscce(2,iel) = temp*visten(2,iel) + cpro_viscls(iel)
        viscce(3,iel) = temp*visten(3,iel) + cpro_viscls(iel)
        viscce(4,iel) = temp*visten(4,iel)
        viscce(5,iel) = temp*visten(5,iel)
        viscce(6,iel) = temp*visten(6,iel)

      enddo
    endif

    iwarnp = iwarni(ivar)

    if (iwgrec(ivar).eq.1) then
      ! Weighting for gradient
      do iel = 1, ncel
        do isou = 1, 6
          cpro_wgrec_v(isou,iel) = viscce(isou,iel)
        enddo
      enddo
    endif

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
      dt     ,                                                       &
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

if (ippmod(idarcy).eq.1) then
  call field_get_name(ivarfl(isca(iscal)), fname)
  call field_get_id(trim(fname)//'_delay', delay_id)
  call field_get_val_s(delay_id, cpro_delay)
  call field_get_val_prev_s(delay_id, cproa_delay)
  call field_get_val_prev_s_by_name('saturation', cproa_sat)
  call field_get_val_s_by_name('saturation', cpro_sat)
endif

! Without porosity neither Darcy
if ((iporos.eq.0).and.(ippmod(idarcy).eq.-1)) then

  ! --> Unsteady term and mass aggregation term
  do iel = 1, ncel
    rovsdt(iel) = rovsdt(iel)                                                 &
                + istat(ivar)*xcpp(iel)*pcrom(iel)*volume(iel)/dt(iel)
  enddo
! With porosity but not Darcy
elseif (ippmod(idarcy).eq.-1) then
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
! Darcy : we take into account the porosity and delay for underground transport
else
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)*cpro_delay(iel)*cpro_sat(iel)
  enddo
  do iel = 1, ncel
    rovsdt(iel) = (rovsdt(iel)                                                &
                + istat(ivar)*xcpp(iel)*pcrom(iel)*volume(iel)/dt(iel) )      &
                * cpro_delay(iel)*cpro_sat(iel)
  enddo
endif

! Scalar with a Drift:
! compute the convective flux
!----------------------------

call field_get_key_int(iflid, keydri, iscdri)

if (iscdri.ge.1) then
 call driflu &
 !=========
 ( iflid  ,                                                       &
   dt     ,                                                       &
   imasfl , bmasfl ,                                              &
   rovsdt , smbrs  )
endif

! Darcy:
! This step is necessary because the divergence of the
! hydraulic head (or pressure), which is taken into account as the
! 'aggregation term' via codits -> bilsca -> bilsc2,
! is not exactly equal to the loss of mass of water in the unsteady case
! (just as far as the precision of the Newton scheme is good),
! and does not take into account the sorbption of the tracer
! (represented by the 'delay' coefficient).
! We choose to 'cancel' the aggregation term here, and to add to smbr
! (right hand side of the transport equation)
! the necessary correction to take sorption into account and get the exact
! conservation of the mass of tracer.

if (ippmod(idarcy).eq.1) then
  if (darcy_unsteady.eq.1) then
    call divmas(1, imasfl , bmasfl , diverg)
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) - diverg(iel)*cvar_var(iel)              &
        + volume(iel)/dt(iel)*cvar_var(iel)                            &
        *( cproa_delay(iel)*cproa_sat(iel)                             &
         - cpro_delay(iel)*cpro_sat(iel) )
    enddo
    do iel = 1, ncel
      rovsdt(iel) = rovsdt(iel) + thetv*diverg(iel)
    enddo
  endif
endif

!===============================================================================
! 3. Solving
!===============================================================================

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ndircp = ndircl(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
iescap = 0
iswdyp = iswdyn(ivar)
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
 ( idtvar , ivar   , iconvp , idiffp , ndircp ,                   &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   cvara_var       , cvara_var       ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscce , viscf  , viscb  , viscce ,          &
   weighf , weighb ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbrs  , cvar_var        , dpvar  ,                   &
   xcpp   , rvoid  )

!===============================================================================
! 4. Writing and clipping
!===============================================================================

call clpsca(iscal)

if (idilat.ge.4.and.itspdv.eq.1) then

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

    propce(iel,ipproc(iustdy(iscal))) =                               &
      propce(iel,ipproc(iustdy(iscal))) - rhovst*cvar_var(iel)

  enddo

endif

! Store the implicit part of the radiative source term
if (idilat.ge.4.and.iirayo.ge.1.and.iscal.eq.iscalt) then
  do iel = 1, ncel
    ivar = isca(iscalt)
    dvar = cvar_var(iel)-cvara_var(iel)
    propce(iel,ipproc(iustdy(iscalt))) = &
    propce(iel,ipproc(iustdy(iscalt)))   &
    - propce(iel,ipproc(itsri(1)))*dvar*volume(iel)
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
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                                 &
            - istat(ivar)*xcpp(iel)*(pcrom(iel)/dt(iel))*volume(iel)        &
                *(cvar_var(iel)-cvara_var(iel))*ibcl
  enddo
  isqrt = 1
  call prodsc(ncel,isqrt,smbrs,smbrs,sclnor)
  write(nfecra,1200)chaine(1:16) ,sclnor
endif

! Free memory
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

!----
! End
!----

return

end subroutine
