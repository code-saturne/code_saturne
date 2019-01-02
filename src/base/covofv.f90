!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file covofv.f90
!>
!> \brief This subroutine performs the solving the convection/diffusion
!> equation (with eventually source terms and/or drift) for a vectorial quantity
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
!> \param[in]     iterns        Navier-Stokes iteration number
!> \param[in]     iscal         scalar number
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        type of mass source term for the variables
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!> \param[in]     viscf         visc*surface/dist at internal faces
!> \param[in]     viscb         visc*surface/dist at boundary faces
!_______________________________________________________________________________

subroutine covofv &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iterns , iscal  ,                                              &
   icepdc , icetsm ,                                              &
   itypsm ,                                                       &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iterns , iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)

! Local variables

logical          lprev
character(len=80) :: chaine
integer          ivar
integer          ifac , iel, isou, jsou
integer          iiun, ibcl
integer          ivarsc
integer          iiscav
integer          ifcvsl, iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap, ivissv
integer          idftnp, iswdyp
integer          iflid , st_prv_id, st_id,  keydri, iscdri
integer          icvflb, f_dim, iflwgr
integer          key_buoyant_id, is_buoyant_fld

integer          ivoid(1)

double precision epsrgp, climgp, relaxp, blencp, epsilp
double precision epsrsp
double precision sclnor
double precision thetv , thets , thetap, thetp1
double precision smbexp(3)
double precision temp, idifftp
double precision turb_schmidt

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:,:) :: smbrv, gavinj
double precision, allocatable, dimension(:,:,:) :: fimp
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb

double precision, dimension(:,:), pointer :: visten
double precision, dimension(:,:), pointer :: cpro_wgrec_v
double precision, dimension(:), pointer :: cpro_wgrec_s
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom, croma, pcrom
double precision, dimension(:,:), pointer :: coefap, cofafp
double precision, dimension(:,:,:), pointer :: coefbp, cofbfp
double precision, dimension(:), pointer :: visct
double precision, dimension(:,:), pointer :: cpro_vect_st, cproa_vect_st
double precision, dimension(:), pointer :: cpro_viscls
! Darcy arrays
double precision, allocatable, dimension(:) :: diverg, divflu
double precision, dimension(:), pointer :: cpro_delay, cpro_sat
double precision, dimension(:), pointer :: cproa_delay, cproa_sat
double precision, dimension(:,:), pointer :: cvar_var, cvara_var

type(var_cal_opt) :: vcopt

type(gwf_soilwater_partition) :: sorption_scal

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! --- Variable number
ivar   = isca(iscal)
! Index of the field
iflid = ivarfl(ivar)

call field_get_val_v(iflid, cvar_var)
call field_get_val_prev_v(iflid, cvara_var)

! Key id for buoyant field (inside the Navier Stokes loop)
call field_get_key_id("is_buoyant", key_buoyant_id)
call field_get_key_int(iflid, key_buoyant_id, is_buoyant_fld)

! If the vector is buoyant, it is inside the Navier Stokes loop, and so iterns >=1
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

call field_get_key_struct_var_cal_opt(iflid, vcopt)

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
allocate(smbrv(3,ncelet), fimp(3,3,ncelet))
allocate(weighf(2,nfac))
allocate(weighb(nfabor))

if (ippmod(idarcy).eq.1) then
  allocate(diverg(ncelet))
endif

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
call field_get_val_s(ivisct, visct)

call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_viscls)
endif

! --- Numero de propriété du terme source si extrapolation
call field_get_key_int(iflid, kstprv, st_prv_id)
if (st_prv_id .ge.0) then
  call field_get_val_v(st_prv_id, cproa_vect_st)
else
  cproa_vect_st => null()
endif

! S pour Source, V pour Variable
thets  = thetss(iscal)
thetv  = vcopt%thetav

call field_get_name(iflid, chaine)

if(vcopt%iwarni.ge.1) then
  write(nfecra,1000) chaine(1:16)
endif

! Retrieve turbulent Schmidt value for current vector
call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)

!===============================================================================
! 2. Source terms
!===============================================================================

! --> Initialization

do iel = 1, ncel
  do isou = 1, 3
    do jsou = 1, 3
      fimp(isou,jsou,iel) = 0.d0
    enddo
    smbrv(isou,iel) = 0.d0
  enddo
enddo

call ustsvv &
!==========
( nvar   , nscal  , ncepdp , ncesmp ,                            &
  iscal  ,                                                       &
  icepdc , icetsm , itypsm ,                                     &
  dt     ,                                                       &
  ckupdc , smacel , smbrv  , fimp )

! Store the source terms for convective limiter
call field_get_key_int(iflid, kst, st_id)
if (st_id .ge.0) then
  call field_get_dim(st_id, f_dim)
  if (f_dim.ne.3) then
    call csexit(1)
  endif
  call field_get_val_v(st_id, cpro_vect_st)

  do iel = 1, ncel
    !Fill the scalar source term field
    do isou = 1, 3
      cpro_vect_st(isou,iel) = smbrv(isou,iel)
    enddo
  end do
  ! Handle parallelism and periodicity
  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(cpro_vect_st)
  endif
end if

! Si on extrapole les TS :
!   SMBRV recoit -theta TS du pas de temps precedent
!     (on aurait pu le faire avant ustssc, mais avec le risque que
!      l'utilisateur l'ecrase)
!   SMBRV recoit la partie du terme source qui depend de la variable
!   A l'ordre 2, on suppose que le ROVSDT fourni par l'utilisateur est <0
!     on implicite le terme (donc ROVSDT*RTPA va dans SMBRV)
!   En std, on adapte le traitement au signe de ROVSDT, mais ROVSDT*RTPA va
!     quand meme dans SMBRV (pas d'autre choix)
if (st_prv_id .ge. 0) then
  do iel = 1, ncel
    ! Stockage temporaire pour economiser un tableau
    do isou = 1, 3
      smbexp(isou) = cproa_vect_st(isou,iel)

      ! Terme source utilisateur explicite
      cproa_vect_st(isou,iel) = smbrv(isou,iel)
      ! Terme source du pas de temps precedent et
      ! On suppose -fimp(isou,isou) > 0 : on implicite
      !    le terme source utilisateur (le reste)
      smbrv(isou,iel) = fimp(isou,isou,iel)*cvara_var(isou,iel)               &
                      - thets*smbexp(isou)
      ! Diagonale
      fimp(isou,isou,iel) = - thetv*fimp(isou,isou,iel)
    enddo
  enddo

! Si on n'extrapole pas les TS :
else
  do iel = 1, ncel
    do isou = 1, 3
      ! Terme source utilisateur
      smbrv(isou,iel) = smbrv(isou,iel) &
                      + fimp(isou,isou,iel)*cvara_var(isou,iel)
      ! Diagonale
      fimp(isou,isou,iel) = max(-fimp(isou,isou,iel),zero)
    enddo
  enddo
endif

! --> Physique particulieres
!     Ordre 2 non pris en compte

if (ippmod(iphpar).ge.1) then
  call pptsvv(iscal, smbrv, fimp)
endif

! Mass source term
if (ncesmp.gt.0) then

  ! Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

  ! On incremente SMBRV par -Gamma RTPA et ROVSDT par Gamma
  allocate(gavinj(3,ncelet))
  call catsmv &
 ( ncelet , ncel   , ncesmp , iiun   , isso2t(iscal) ,            &
   icetsm , itypsm(1,ivar) ,                                      &
   cell_f_vol , cvara_var    , smacel(1,ivar) , smacel(1,ipr),    &
   smbrv  , fimp , gavinj)

  ! Si on extrapole les TS on met Gamma Pinj dans cproa_vect_st
  if (st_prv_id .ge. 0) then
    do iel = 1, ncel
      do isou = 1, 3
        cproa_vect_st(isou,iel) = cproa_vect_st(isou,iel) + gavinj(isou,iel)
      enddo
    enddo
  ! Sinon on le met directement dans SMBRV
  else
    do iel = 1, ncel
      do isou = 1, 3
        smbrv(isou,iel) = smbrv(isou,iel) + gavinj(isou,iel)
      enddo
    enddo
  endif

endif

if (st_prv_id .ge. 0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    do isou = 1, 3
      smbrv(isou,iel) = smbrv(isou,iel) + thetp1 * cproa_vect_st(isou,iel)
    enddo
  enddo
endif

! Compressible algorithm
! or Low Mach compressible algos with mass flux prediction
if (ippmod(icompf).ge.0.or.(idilat.gt.1.and.ipredfl.eq.1.and.irovar.eq.1)) then
  pcrom => croma

! Low Mach compressible algos (conservative in time).
! Same algo for Volume of Fluid method.
else if ((idilat.gt.1.or.ivofmt.ge.0).and.irovar.eq.1) then
  if (iterns.eq.1) then
    call field_get_val_prev2_s(icrom, pcrom)
  else
    call field_get_val_prev_s(icrom, pcrom)
  endif

! Deprecated algo or constant density
else
  call field_get_val_s(icrom, pcrom)
endif

! Get the the order of the diffusivity tensor of the variable
call field_get_key_int_by_name(iflid, "diffusivity_tensor", idftnp)

! "VITESSE" DE DIFFUSION FACETTE

! On prend le MAX(mu_t,0) car en LES dynamique mu_t peut etre negatif
! (clipping sur (mu + mu_t)). On aurait pu prendre
! MAX(K + K_t,0) mais cela autoriserait des K_t negatif, ce qui est
! considere ici comme non physique.
if (vcopt%idiff.ge.1) then
  ! Scalar diffusivity
  if (iand(vcopt%idften, ISOTROPIC_DIFFUSION).ne.0) then

    allocate(w1(ncelet))
    idifftp = vcopt%idifft
    if (ifcvsl.lt.0) then
      do iel = 1, ncel
        w1(iel) = visls0(iscal)                                     &
           + idifftp*max(visct(iel),zero)/turb_schmidt
      enddo
    else
      do iel = 1, ncel
        w1(iel) = cpro_viscls(iel)                                &
           + idifftp*max(visct(iel),zero)/turb_schmidt
      enddo
    endif

    if (vcopt%iwgrec.eq.1) then
      ! Weighting for gradient
      do iel = 1, ncel
        cpro_wgrec_s(iel) = w1(iel)
      enddo
      call synsca(cpro_wgrec_s)
    endif

    call viscfa &
   ( imvisf ,                      &
     w1     ,                      &
     viscf  , viscb  )

    deallocate(w1)

  ! Symmetric tensor diffusivity (GGDH)
  elseif (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then

    ! Allocate temporary arrays
    allocate(viscce(6,ncelet))

    if (iturb.ne.32) then
      call field_get_val_v(ivsten, visten)
    else ! EBRSM and (GGDH or AFM)
      call field_get_val_v(ivstes, visten)
    endif

    if (ifcvsl.lt.0) then
      do iel = 1, ncel

        temp = vcopt%idifft*ctheta(iscal)/csrij
        viscce(1,iel) = temp*visten(1,iel) + visls0(iscal)
        viscce(2,iel) = temp*visten(2,iel) + visls0(iscal)
        viscce(3,iel) = temp*visten(3,iel) + visls0(iscal)
        viscce(4,iel) = temp*visten(4,iel)
        viscce(5,iel) = temp*visten(5,iel)
        viscce(6,iel) = temp*visten(6,iel)

      enddo
    else
      do iel = 1, ncel

        temp = vcopt%idifft*ctheta(iscal)/csrij
        viscce(1,iel) = temp*visten(1,iel) + cpro_viscls(iel)
        viscce(2,iel) = temp*visten(2,iel) + cpro_viscls(iel)
        viscce(3,iel) = temp*visten(3,iel) + cpro_viscls(iel)
        viscce(4,iel) = temp*visten(4,iel)
        viscce(5,iel) = temp*visten(5,iel)
        viscce(6,iel) = temp*visten(6,iel)

      enddo
    endif

    if (vcopt%iwgrec.eq.1) then
      ! Weighting for gradient
      do iel = 1, ncel
        do isou = 1, 6
          cpro_wgrec_v(isou,iel) = viscce(isou,iel)
        enddo
      enddo
      call syntis(cpro_wgrec_v)
    endif

    iwarnp = vcopt%iwarni

    call vitens &
    ( viscce , iwarnp ,             &
      weighf , weighb ,             &
      viscf  , viscb  )

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
  call field_get_key_struct_gwf_soilwater_partition(iflid, sorption_scal)
  call field_get_val_s(sorption_scal%idel, cpro_delay)
  call field_get_val_prev_s(sorption_scal%idel, cproa_delay)
  call field_get_val_prev_s_by_name('saturation', cproa_sat)
  call field_get_val_s_by_name('saturation', cpro_sat)
endif

! Not Darcy
if (ippmod(idarcy).eq.-1) then
  ! --> Unsteady term and mass aggregation term
  do iel = 1, ncel
    do isou = 1, 3
      fimp(isou,isou,iel) = fimp(isou,isou,iel)                               &
                          + vcopt%istat*pcrom(iel)*cell_f_vol(iel)/dt(iel)
    enddo
  enddo
! Darcy : we take into account the porosity and delay for underground transport
else
  do iel = 1, ncel
    do isou = 1, 3
      smbrv(isou,iel) = smbrv(isou,iel)*cpro_delay(iel)*cpro_sat(iel)
      fimp(isou,isou,iel) = (fimp(isou,isou,iel)                              &
                          + vcopt%istat*pcrom(iel)*volume(iel)/dt(iel) )      &
                            * cpro_delay(iel)*cpro_sat(iel)
    enddo
  enddo
endif

! Scalar with a Drift:
! compute the convective flux
!----------------------------

call field_get_key_int(iflid, keydri, iscdri)

if (iscdri.ge.1) then
  allocate(divflu(ncelet))

  call driflu &
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
    do isou = 1, 3
      fimp(isou,isou,iel) = fimp(isou,isou,iel) + iconvp*thetap*divflu(iel)
      smbrv(isou,iel) = smbrv(isou,iel) - iconvp*divflu(iel)*cvara_var(isou,iel)
    enddo
  enddo

  deallocate(divflu)
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
      do isou = 1, 3
        smbrv(isou,iel) = smbrv(isou,iel) - diverg(iel)*cvar_var(isou,iel)     &
          + cell_f_vol(iel)/dt(iel)*cvar_var(isou,iel)                         &
          *( cproa_delay(iel)*cproa_sat(iel)                                   &
           - cpro_delay(iel)*cpro_sat(iel) )
      enddo
    enddo
    do iel = 1, ncel
      do isou = 1, 3
        fimp(isou,isou,iel) = fimp(isou,isou,iel) + thetv*diverg(iel)
      enddo
    enddo
  endif
endif

!===============================================================================
! 3. Solving
!===============================================================================

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
relaxp = vcopt%relaxv
! all boundary convective flux with upwind
icvflb = 0
! transposed gradient term only for NS
ivissv = 0

call field_get_coefa_v(iflid, coefap)
call field_get_coefb_v(iflid, coefbp)
call field_get_coefaf_v(iflid, cofafp)
call field_get_coefbf_v(iflid, cofbfp)

call coditv &
!==========
 ( idtvar , iterns , iflid  , iconvp , idiffp , ndircp ,          &
   imrgra , nswrsp , nswrgp , imligp , ircflp , ivissv ,          &
   ischcp , isstpp , iescap , idftnp , iswdyp ,                   &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp ,                   &
   relaxp , thetv  ,                                              &
   cvara_var       , cvara_var       ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  , rvoid  ,          &
   viscce , weighf , weighb ,                                     &
   icvflb , ivoid  ,                                              &
   fimp   , smbrv  , cvar_var        ,                            &
   rvoid  )

!===============================================================================
! 4. Writing
!===============================================================================

! BILAN EXPLICITE (VOIR CODITS : ON ENLEVE L'INCREMENT)
! Ceci devrait etre valable avec le theta schema sur les Termes source

if (vcopt%iwarni.ge.2) then
  if (vcopt%nswrsm.gt.1) then
    ibcl = 1
  else
    ibcl = 0
  endif
  do iel = 1, ncel
    do isou = 1, 3
      smbrv(isou,iel) = smbrv(isou,iel)                                    &
                      - vcopt%istat*(pcrom(iel)/dt(iel))*cell_f_vol(iel)   &
                        *(cvar_var(isou,iel)-cvara_var(isou,iel))*ibcl
    enddo
  enddo
  sclnor = sqrt(cs_gdot(3*ncel,smbrv,smbrv))
  write(nfecra,1200)chaine(1:16) ,sclnor
endif

! Free memory
deallocate(smbrv, fimp)
deallocate(weighf, weighb)
if (allocated(viscce)) deallocate(viscce)
if (allocated(diverg)) deallocate(diverg)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A16                       ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A16,' : BILAN EXPLICITE = ',E14.5)
! 9000 format( &
!'@                                                            ',/,&
!'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
!'@                                                            ',/,&
!'@ @@ ATTENTION : ERREUR DANS COVOFV                          ',/,&
!'@    =========                                               ',/,&
!'@    IVARSC DOIT ETRE UN ENTIER POSITIF STRICTEMENT          ',/,&
!'@    IL VAUT ICI ',I10                                        ,/,&
!'@                                                            ',/,&
!'@  Le calcul ne peut etre execute.                           ',/,&
!'@                                                            ',/,&
!'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
!'@                                                            ',/)

#else

 1000 format(/,                                                   &
'   ** SOLVING VARIABLE ',A16                                  ,/,&
'      ----------------'                                       ,/)
 1200 format(1X,A16,' : EXPLICIT BALANCE = ',E14.5)
! 9000 format( &
!'@'                                                            ,/,&
!'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
!'@'                                                            ,/,&
!'@ @@ WARNING: ERROR IN COVOFV'                                ,/,&
!'@    ========'                                                ,/,&
!'@    IVARSC MUST BE A STRICTLY POSITIVE INTEGER'              ,/,&
!'@    ITS VALUE IS ',I10                                       ,/,&
!'@'                                                            ,/,&
!'@  The calculation will not be run.'                          ,/,&
!'@'                                                            ,/,&
!'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
!'@                                                            ',/)

#endif

!----
! End
!----

return

end subroutine
