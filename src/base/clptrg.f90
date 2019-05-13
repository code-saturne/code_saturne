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
! Function :
! --------

!> \file clptrg.f90
!>
!> \brief Boundary conditions for rough walls (icodcl = 6).
!>
!> The wall functions may change the value of the diffusive flux.
!>
!> The values at a boundary face \f$ \fib \f$ stored in the face center
!> \f$ \centf \f$ of the variable \f$ P \f$ and its diffusive flux \f$ Q \f$
!> are written as:
!> \f[
!> P_\centf = A_P^g + B_P^g P_\centi
!> \f]
!> and
!> \f[
!> Q_\centf = A_P^f + B_P^f P_\centi
!> \f]
!> where \f$ P_\centi \f$ is the value of the variable \f$ P \f$ at the
!> neighboring cell.
!>
!> Warning:
!>
!> - for a vector field such as the velocity \f$ \vect{u} \f$ the boundary
!>   conditions may read:
!>   \f[
!>   \vect{u}_\centf = \vect{A}_u^g + \tens{B}_u^g \vect{u}_\centi
!>   \f]
!>   and
!>   \f[
!>   \vect{Q}_\centf = \vect{A}_u^f + \tens{B}_u^f \vect{u}_\centi
!>   \f]
!>   where \f$ \tens{B}_u^g \f$ and \f$ \tens{B}_u^f \f$ are 3x3 tensor matrix
!>   which coupled veclocity components next to a boundary.
!>
!> Please refer to the <a href="../../theory.pdf#cpltrg"><b>clptrg</b></a> section
!> of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nscal         total number of scalars
!> \param[in]     isvhb         indicator to save exchange coeffient
!> \param[in,out] icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!> \param[in]     velipb        value of the velocity at \f$ \centip \f$
!>                               of boundary cells
!> \param[in]     rijipb        value of \f$ R_{ij} \f$ at \f$ \centip \f$
!>                               of boundary cells
!> \param[out]    visvdr        viscosite dynamique ds les cellules
!>                               de bord apres amortisst de v driest
!> \param[out]    hbord         coefficients d'echange aux bords
!>
!> \param[in]     theipb        boundary temperature in \f$ \centip \f$
!>                               (more exaclty the energetic variable)
!_______________________________________________________________________________

subroutine clptrg &
 ( nscal  , isvhb  , icodcl ,                                     &
   rcodcl ,                                                       &
   velipb , rijipb , visvdr ,                                     &
   hbord  , theipb )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use pointe
use entsor
use albase
use parall
use ppppar
use ppthch
use ppincl
use radiat
use cplsat
use mesh
use field
use lagran
use turbomachinery
use cs_c_bindings
use atincl, only: itotwt, iliqwt, modsedi, moddep

!===============================================================================

implicit none

! Arguments

integer          nscal, isvhb

integer, pointer, dimension(:,:) :: icodcl

double precision, pointer, dimension(:,:,:) :: rcodcl
double precision, dimension(:,:) :: velipb
double precision, pointer, dimension(:,:) :: rijipb
double precision, pointer, dimension(:) :: visvdr, hbord, theipb

! Local variables

integer          ifac, iel, isou, ii, jj, kk
integer          iscal, clsyme
integer          modntl
integer          iuntur, iuiptn, f_id, iustar
integer          kdflim

double precision rnx, rny, rnz, rxnn
double precision tx, ty, tz, txn, txn0, t2x, t2y, t2z
double precision utau, upx, upy, upz, usn
double precision uiptn, uiptmn, uiptmx
double precision uetmax, uetmin, ukmax, ukmin, yplumx, yplumn
double precision tetmax, tetmin, tplumx, tplumn
double precision uk, uet, yplus, uplus, phit
double precision gredu, temp
double precision cfnnu, cfnns, cfnnk, cfnne
double precision sqrcmu, ek
double precision xmutlm
double precision rcprod, rcflux
double precision hflui, hint, pimp, qimp
double precision und0, deuxd0
double precision eloglo(3,3), alpha(6,6)
double precision rcodcx, rcodcy, rcodcz, rcodcn
double precision visclc, visctc, romc  , distbf, srfbnf
double precision cofimp
double precision distb0, rugd  , ydep
double precision dsa0
double precision rinfiv(3)
double precision visci(3,3), fikis, viscis, distfi
double precision fcoefa(6), fcoefb(6), fcofaf(6)
double precision fcofbf(6), fcofad(6), fcofbd(6)

double precision rxx, rxy, rxz, ryy, ryz, rzz, rnnb
double precision rttb, alpha_rnn, liqwt, totwt

double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: viscl, visct, cpro_cp, yplbr, ustar
double precision, dimension(:), allocatable :: byplus, buk
double precision, dimension(:), allocatable, target :: buet, bcfnns_loc

double precision, dimension(:), pointer :: cvar_k, bcfnns
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:), pointer :: cvara_nusa

double precision, dimension(:), pointer :: cvar_totwt, cvar_t, cpro_liqwt
double precision, dimension(:), pointer :: cpro_rugt
double precision, dimension(:), pointer :: cpro_diff_lim_k
double precision, dimension(:), pointer :: cpro_diff_lim_eps
double precision, dimension(:), pointer :: cpro_diff_lim_rij

double precision, dimension(:,:), pointer :: coefau, cofafu, visten
double precision, dimension(:,:,:), pointer :: coefbu, cofbfu
double precision, dimension(:), pointer :: coefa_k, coefb_k, coefaf_k, coefbf_k
double precision, dimension(:), pointer :: coefa_ep, coefaf_ep
double precision, dimension(:), pointer :: coefb_ep, coefbf_ep
double precision, dimension(:), pointer :: coefa_r11, coefaf_r11, coefad_r11
double precision, dimension(:), pointer :: coefb_r11, coefbf_r11, coefbd_r11
double precision, dimension(:), pointer :: coefa_r22, coefaf_r22, coefad_r22
double precision, dimension(:), pointer :: coefb_r22, coefbf_r22, coefbd_r22
double precision, dimension(:), pointer :: coefa_r33, coefaf_r33, coefad_r33
double precision, dimension(:), pointer :: coefb_r33, coefbf_r33, coefbd_r33
double precision, dimension(:), pointer :: coefa_r12, coefaf_r12, coefad_r12
double precision, dimension(:), pointer :: coefb_r12, coefbf_r12, coefbd_r12
double precision, dimension(:), pointer :: coefa_r13, coefaf_r13, coefad_r13
double precision, dimension(:), pointer :: coefb_r13, coefbf_r13, coefbd_r13
double precision, dimension(:), pointer :: coefa_r23, coefaf_r23, coefad_r23
double precision, dimension(:), pointer :: coefb_r23, coefbf_r23, coefbd_r23
double precision, dimension(:,:), pointer :: coefa_rij, coefaf_rij, coefad_rij
double precision, dimension(:,:,:), pointer :: coefb_rij, coefbf_rij, coefbd_rij
double precision, dimension(:), pointer :: coefa_omg, coefaf_omg
double precision, dimension(:), pointer :: coefb_omg, coefbf_omg
double precision, dimension(:), pointer :: coefa_al, coefaf_al
double precision, dimension(:), pointer :: coefb_al, coefbf_al
double precision, dimension(:), pointer :: coefa_phi, coefaf_phi
double precision, dimension(:), pointer :: coefb_phi, coefbf_phi
double precision, dimension(:), pointer :: coefa_fb, coefaf_fb
double precision, dimension(:), pointer :: coefb_fb, coefbf_fb
double precision, dimension(:), pointer :: coefa_nu, coefaf_nu
double precision, dimension(:), pointer :: coefb_nu, coefbf_nu

integer          ntlast , iaff
data             ntlast , iaff /-1 , 0/
save             ntlast , iaff

type(var_cal_opt) :: vcopt
type(var_cal_opt) :: vcopt_rij, vcopt_ep

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine clptrg_scalar(iscal, isvhb, icodcl, rcodcl,              &
                           byplus, buk, buet, bcfnns, hbord, theipb,  &
                           tetmax , tetmin , tplumx , tplumn)

    implicit none
    integer          iscal, isvhb
    integer, pointer, dimension(:,:) :: icodcl
    double precision, pointer, dimension(:,:,:) :: rcodcl
    double precision, dimension(:) :: byplus, buk, buet, bcfnns
    double precision, pointer, dimension(:) :: hbord, theipb
    double precision tetmax, tetmin, tplumx, tplumn

  end subroutine clptrg_scalar

 end interface

!===============================================================================
! 1. Initializations
!===============================================================================

! Initialize variables to avoid compiler warnings

ek = 0.d0
phit = 0.d0
uiptn = 0.d0

rinfiv(1) = rinfin
rinfiv(2) = rinfin
rinfiv(3) = rinfin

! --- Constants
uet = 1.d0
utau = 1.d0
sqrcmu = sqrt(cmu)

! --- Correction factors for stratification (used in atmospheric models)
cfnnu=1.d0
cfnns=1.d0
cfnnk=1.d0
cfnne=1.d0

und0   = 1.d0
deuxd0 = 2.d0

! Alpha constant for a realisable BC for R12 with the SSG model
alpha_rnn = 0.47d0

! pointers to y+ if saved

yplbr => null()

if (iyplbr.ge.0) then
  call field_get_val_s(iyplbr, yplbr)
endif
if (itytur.eq.3 .and. idirsm.eq.1) call field_get_val_v(ivsten, visten)

! Diffusion limiter
call field_get_key_id("diffusion_limiter_id", kdflim)

! --- Store wall friction velocity

call field_get_id_try('ustar', iustar)
if (iustar.ge.0) then !TODO remove, this information is in cofaf cofbf
  call field_get_val_s(iustar, ustar)
else
  allocate(buet(nfabor))
  ustar => buet
endif

! --- Gradient and flux Boundary Conditions

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)
call field_get_coefaf_v(ivarfl(iu), cofafu)
call field_get_coefbf_v(ivarfl(iu), cofbfu)

if (ik.gt.0) then
  call field_get_coefa_s(ivarfl(ik), coefa_k)
  call field_get_coefb_s(ivarfl(ik), coefb_k)
  call field_get_coefaf_s(ivarfl(ik), coefaf_k)
  call field_get_coefbf_s(ivarfl(ik), coefbf_k)

  ! Diffuion limiter
  call field_get_key_int(ivarfl(ik), kdflim, f_id)
  if (f_id.ge.0) then
    call field_get_val_s(f_id, cpro_diff_lim_k)
  else
    cpro_diff_lim_k => null()
  endif
else
  coefa_k => null()
  coefb_k => null()
  coefaf_k => null()
  coefbf_k => null()
endif

if (iep.gt.0) then
  call field_get_coefa_s(ivarfl(iep), coefa_ep)
  call field_get_coefb_s(ivarfl(iep), coefb_ep)
  call field_get_coefaf_s(ivarfl(iep), coefaf_ep)
  call field_get_coefbf_s(ivarfl(iep), coefbf_ep)
  ! Diffuion limiter
  call field_get_key_int(ivarfl(iep), kdflim, f_id)
  if (f_id.ge.0) then
    call field_get_val_s(f_id, cpro_diff_lim_eps)
  else
    cpro_diff_lim_eps => null()
  endif
else
  coefa_ep => null()
  coefb_ep => null()
  coefaf_ep => null()
  coefbf_ep => null()
endif

if (itytur.eq.3) then! Also have boundary conditions for the momentum equation
  call field_get_key_struct_var_cal_opt(ivarfl(irij), vcopt_rij)
  call field_get_key_struct_var_cal_opt(ivarfl(iep), vcopt_ep)

  if (irijco.eq.1) then
    call field_get_coefa_v(ivarfl(irij), coefa_rij)
    call field_get_coefb_v(ivarfl(irij), coefb_rij)
    call field_get_coefaf_v(ivarfl(irij), coefaf_rij)
    call field_get_coefbf_v(ivarfl(irij), coefbf_rij)
    call field_get_coefad_v(ivarfl(irij), coefad_rij)
    call field_get_coefbd_v(ivarfl(irij), coefbd_rij)

    coefb_r11 => null()
    coefaf_r11 => null()
    coefbf_r11 => null()
    coefad_r11 => null()
    coefbd_r11 => null()

    coefa_r22 => null()
    coefb_r22 => null()
    coefaf_r22 => null()
    coefbf_r22 => null()
    coefad_r22 => null()
    coefbd_r22 => null()

    coefa_r33 => null()
    coefb_r33 => null()
    coefaf_r33 => null()
    coefbf_r33 => null()
    coefad_r33 => null()
    coefbd_r33 => null()

    coefa_r12 => null()
    coefb_r12 => null()
    coefaf_r12 => null()
    coefbf_r12 => null()
    coefad_r12 => null()
    coefbd_r12 => null()

    coefa_r13 => null()
    coefb_r13 => null()
    coefaf_r13 => null()
    coefbf_r13 => null()
    coefad_r13 => null()
    coefbd_r13 => null()

    coefa_r23 => null()
    coefb_r23 => null()
    coefaf_r23 => null()
    coefbf_r23 => null()
    coefad_r23 => null()
    coefbd_r23 => null()
    ! Diffuion limiter
    call field_get_key_int(ivarfl(irij), kdflim, f_id)
    if (f_id.ge.0) then
      call field_get_val_s(f_id, cpro_diff_lim_rij)
    else
      cpro_diff_lim_rij => null()
    endif

  else

    call field_get_coefa_s(ivarfl(ir11), coefa_r11)
    call field_get_coefb_s(ivarfl(ir11), coefb_r11)
    call field_get_coefaf_s(ivarfl(ir11), coefaf_r11)
    call field_get_coefbf_s(ivarfl(ir11), coefbf_r11)
    call field_get_coefad_s(ivarfl(ir11), coefad_r11)
    call field_get_coefbd_s(ivarfl(ir11), coefbd_r11)

    call field_get_coefa_s(ivarfl(ir22), coefa_r22)
    call field_get_coefb_s(ivarfl(ir22), coefb_r22)
    call field_get_coefaf_s(ivarfl(ir22), coefaf_r22)
    call field_get_coefbf_s(ivarfl(ir22), coefbf_r22)
    call field_get_coefad_s(ivarfl(ir22), coefad_r22)
    call field_get_coefbd_s(ivarfl(ir22), coefbd_r22)

    call field_get_coefa_s(ivarfl(ir33), coefa_r33)
    call field_get_coefb_s(ivarfl(ir33), coefb_r33)
    call field_get_coefaf_s(ivarfl(ir33), coefaf_r33)
    call field_get_coefbf_s(ivarfl(ir33), coefbf_r33)
    call field_get_coefad_s(ivarfl(ir33), coefad_r33)
    call field_get_coefbd_s(ivarfl(ir33), coefbd_r33)

    call field_get_coefa_s(ivarfl(ir12), coefa_r12)
    call field_get_coefb_s(ivarfl(ir12), coefb_r12)
    call field_get_coefaf_s(ivarfl(ir12), coefaf_r12)
    call field_get_coefbf_s(ivarfl(ir12), coefbf_r12)
    call field_get_coefad_s(ivarfl(ir12), coefad_r12)
    call field_get_coefbd_s(ivarfl(ir12), coefbd_r12)

    call field_get_coefa_s(ivarfl(ir13), coefa_r13)
    call field_get_coefb_s(ivarfl(ir13), coefb_r13)
    call field_get_coefaf_s(ivarfl(ir13), coefaf_r13)
    call field_get_coefbf_s(ivarfl(ir13), coefbf_r13)
    call field_get_coefad_s(ivarfl(ir13), coefad_r13)
    call field_get_coefbd_s(ivarfl(ir13), coefbd_r13)

    call field_get_coefa_s(ivarfl(ir23), coefa_r23)
    call field_get_coefb_s(ivarfl(ir23), coefb_r23)
    call field_get_coefaf_s(ivarfl(ir23), coefaf_r23)
    call field_get_coefbf_s(ivarfl(ir23), coefbf_r23)
    call field_get_coefad_s(ivarfl(ir23), coefad_r23)
    call field_get_coefbd_s(ivarfl(ir23), coefbd_r23)

    coefa_rij => null()
    coefb_rij => null()
    coefaf_rij => null()
    coefbf_rij => null()
    coefad_rij => null()
    coefbd_rij => null()
  endif
endif

if (ial.gt.0) then
  call field_get_coefa_s(ivarfl(ial), coefa_al)
  call field_get_coefb_s(ivarfl(ial), coefb_al)
  call field_get_coefaf_s(ivarfl(ial), coefaf_al)
  call field_get_coefbf_s(ivarfl(ial), coefbf_al)
else
  coefa_al => null()
  coefb_al => null()
  coefaf_al => null()
  coefbf_al => null()
endif

if (iomg.gt.0) then
  call field_get_coefa_s(ivarfl(iomg), coefa_omg)
  call field_get_coefb_s(ivarfl(iomg), coefb_omg)
  call field_get_coefaf_s(ivarfl(iomg), coefaf_omg)
  call field_get_coefbf_s(ivarfl(iomg), coefbf_omg)
else
  coefa_omg => null()
  coefb_omg => null()
  coefaf_omg => null()
  coefbf_omg => null()
endif

if (iphi.gt.0) then
  call field_get_coefa_s(ivarfl(iphi), coefa_phi)
  call field_get_coefb_s(ivarfl(iphi), coefb_phi)
  call field_get_coefaf_s(ivarfl(iphi), coefaf_phi)
  call field_get_coefbf_s(ivarfl(iphi), coefbf_phi)
else
  coefa_phi => null()
  coefb_phi => null()
  coefaf_phi => null()
  coefbf_phi => null()
endif

if (ifb.gt.0) then
  call field_get_coefa_s(ivarfl(ifb), coefa_fb)
  call field_get_coefb_s(ivarfl(ifb), coefb_fb)
  call field_get_coefaf_s(ivarfl(ifb), coefaf_fb)
  call field_get_coefbf_s(ivarfl(ifb), coefbf_fb)
else
  coefa_fb => null()
  coefb_fb => null()
  coefaf_fb => null()
  coefbf_fb => null()
endif

if (inusa.gt.0) then
  call field_get_val_prev_s(ivarfl(inusa), cvara_nusa)
  call field_get_coefa_s(ivarfl(inusa), coefa_nu)
  call field_get_coefaf_s(ivarfl(inusa), coefaf_nu)
  call field_get_coefb_s(ivarfl(inusa), coefb_nu)
  call field_get_coefbf_s(ivarfl(inusa), coefbf_nu)
  call field_get_key_struct_var_cal_opt(ivarfl(inusa), vcopt)
else
  cvara_nusa => null()
  coefa_nu => null()
  coefb_nu => null()
  coefaf_nu => null()
  coefbf_nu => null()
endif

! --- Physical quantities
call field_get_val_s(icrom, crom)

if (itytur.eq.2 .or. itytur.eq.5                             &
    .or. iturb.eq.60 .or. iturb.eq.50 .or. iturb.eq.51) then
  call field_get_val_s(ivarfl(ik), cvar_k)
endif

if(itytur.eq.3) then
  if (irijco.eq.1) then
    call field_get_val_v(ivarfl(irij), cvar_rij)
  else
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
    call field_get_val_s(ivarfl(ir12), cvar_r12)
    call field_get_val_s(ivarfl(ir13), cvar_r13)
    call field_get_val_s(ivarfl(ir23), cvar_r23)
  endif
endif

call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)
if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
endif

! min. and max. of wall tangential velocity
uiptmx = -grand
uiptmn =  grand

! min. and max. of wall friction velocity
uetmax = -grand
uetmin =  grand
ukmax  = -grand
ukmin  =  grand

! min. and max. of y+
yplumx = -grand
yplumn =  grand

! min. and max. of wall friction of the thermal scalar
tetmax = -grand
tetmin =  grand

! min. and max. of T+
tplumx = -grand
tplumn =  grand

! Counter (reversal)
iuiptn = 0


! With v2f type model, (phi-fbar et BL-v2/k) u=0 is set directly, so
! uiptmx and uiptmn are necessarily 0
if (itytur.eq.5) then
  uiptmx = 0.d0
  uiptmn = 0.d0
endif

! Pointers to specific fields
allocate(byplus(nfabor))
allocate(buk(nfabor))

call field_get_id_try("non_neutral_scalar_correction", f_id)
if (f_id.ge.0) then
  call field_get_val_s(f_id, bcfnns)
else
  allocate(bcfnns_loc(nfabor))
  bcfnns => bcfnns_loc
endif

cvar_t => null()
cvar_totwt => null()
cpro_liqwt => null()
cpro_rugt => null()

if (ippmod(iatmos).ge.1) then
  call field_get_val_s(ivarfl(isca(iscalt)), cvar_t)
  if (ippmod(iatmos).eq.2) then
    call field_get_val_s(ivarfl(isca(itotwt)), cvar_totwt)
    call field_get_val_s(iliqwt, cpro_liqwt)

    if (modsedi.eq.1.and.moddep.gt.0) then
      call field_get_val_s_by_name('boundary_thermal_roughness', cpro_rugt)
    endif
  endif
endif

! --- Loop on boundary faces
do ifac = 1, nfabor

  ! Test on the presence of a rough wall
  if (icodcl(ifac,iu).eq.6) then

    iel = ifabor(ifac)

    ! Physical properties
    visclc = viscl(iel)
    visctc = visct(iel)
    romc   = crom(iel)

    ! Geometric quantities
    distbf = distb(ifac)
    srfbnf = surfbn(ifac)

    !===========================================================================
    ! 1. Local framework
    !===========================================================================

    ! Unit normal

    rnx = surfbo(1,ifac)/srfbnf
    rny = surfbo(2,ifac)/srfbnf
    rnz = surfbo(3,ifac)/srfbnf

    ! Handle displacement velocity

    rcodcx = rcodcl(ifac,iu,1)
    rcodcy = rcodcl(ifac,iv,1)
    rcodcz = rcodcl(ifac,iw,1)

    ! If we are not using ALE, force the displacement velocity for the face
    !  to be tangential (and update rcodcl for possible use)
    ! In frozen rotor (iturbo = 1), the velocity is neither tangential to the
    !  wall (absolute velocity solved in a relative frame of reference)
    if (iale.eq.0.and.iturbo.eq.0) then
      rcodcn = rcodcx*rnx+rcodcy*rny+rcodcz*rnz
      rcodcx = rcodcx -rcodcn*rnx
      rcodcy = rcodcy -rcodcn*rny
      rcodcz = rcodcz -rcodcn*rnz
      rcodcl(ifac,iu,1) = rcodcx
      rcodcl(ifac,iv,1) = rcodcy
      rcodcl(ifac,iw,1) = rcodcz
    endif

    ! Relative tangential velocity

    upx = velipb(ifac,1) - rcodcx
    upy = velipb(ifac,2) - rcodcy
    upz = velipb(ifac,3) - rcodcz

    usn = upx*rnx+upy*rny+upz*rnz
    tx  = upx -usn*rnx
    ty  = upy -usn*rny
    tz  = upz -usn*rnz
    txn = sqrt(tx**2 +ty**2 +tz**2)
    utau= txn

    ! Unit tangent

    if (txn.ge.epzero) then

      txn0 = 1.d0

      tx  = tx/txn
      ty  = ty/txn
      tz  = tz/txn

    else

      ! If the velocity is zero,
      !  Tx, Ty, Tz is not used (we cancel the velocity), so we assign any
      !  value (zero for example)

      txn0 = 0.d0

      tx  = 0.d0
      ty  = 0.d0
      tz  = 0.d0

    endif

    ! Complete if necessary for Rij-Epsilon

    if (itytur.eq.3) then

      ! --> T2 = RN X T (where X is the cross product)

      t2x = rny*tz - rnz*ty
      t2y = rnz*tx - rnx*tz
      t2z = rnx*ty - rny*tx

      ! --> Orthogonal matrix for change of reference frame ELOGLOij
      !     (from local to global reference frame)

      !                      |TX  -RNX  T2X|
      !             ELOGLO = |TY  -RNY  T2Y|
      !                      |TZ  -RNZ  T2Z|

      !    Its transpose ELOGLOt is its inverse

      eloglo(1,1) =  tx
      eloglo(1,2) = -rnx
      eloglo(1,3) =  t2x
      eloglo(2,1) =  ty
      eloglo(2,2) = -rny
      eloglo(2,3) =  t2y
      eloglo(3,1) =  tz
      eloglo(3,2) = -rnz
      eloglo(3,3) =  t2z

      ! Compute Reynolds stress transformation matrix

      clsyme = 0
      call turbulence_bc_rij_transform(clsyme, eloglo, alpha)

    endif

    !===========================================================================
    ! 2. Friction velocities
    !===========================================================================

    ! ---> Compute Uet depending if we are in the log zone or not
    !      in 1 or 2 velocity scales
    !      and uk based on ek


    if (abs(utau).le.epzero) utau = epzero

    ! rugd: rugosite de paroi pour les variables dynamiques
    !       seule la valeur stockee pour iu est utilisee
    rugd = rcodcl(ifac,iu,3)

    ! NB: for rough walls, yplus is computed from the roughness and not uk.
    yplus = distbf/rugd

    ! Pseudo shift of wall by rugd ((distbf+rugd)/rugd)
    uet = utau/log(yplus+1.d0)*xkappa

    ! Dimensionless velocity
    uplus = log(yplus+1.d0)/xkappa

    if (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then
      ek = cvar_k(iel)
    else if(itytur.eq.3) then
      if (irijco.eq.1) then
        ek = 0.5d0*(cvar_rij(1,iel)+cvar_rij(2,iel)+cvar_rij(3,iel))
        rxx = cvar_rij(1,iel)
        rxy = cvar_rij(4,iel)
        rxz = cvar_rij(6,iel)
        ryy = cvar_rij(2,iel)
        ryz = cvar_rij(5,iel)
        rzz = cvar_rij(3,iel)
      else
        ek = 0.5d0*(cvar_r11(iel)+cvar_r22(iel)+cvar_r33(iel))
        rxx = cvar_r11(iel)
        rxy = cvar_r12(iel)
        rxz = cvar_r13(iel)
        ryy = cvar_r22(iel)
        ryz = cvar_r23(iel)
        rzz = cvar_r33(iel)
      endif
      rnnb =   rnx * (rxx * rnx + rxy * rny + rxz * rnz) &
        + rny * (rxy * rnx + ryy * rny + ryz * rnz) &
        + rnz * (rxz * rnx + ryz * rny + rzz * rnz)

      rttb =   tx * (rxx * tx + rxy * ty + rxz * tz) &
        + ty * (rxy * tx + ryy * ty + ryz * tz) &
        + tz * (rxz * tx + ryz * ty + rzz * tz)
    endif

    if (iwallf.le.2) then
      uk = uet

    ! Si iwallf= 3, 4 ou 5 on calcule uk et uet
    else

      uk = cmu025*sqrt(ek)

    endif

    if (ippmod(iatmos).ge.1) then

      ! Compute reduced gravity for non horizontal walls :
      gredu = gx*rnx + gy*rny + gz*rnz

      temp = cvar_t(iel)
      totwt = 0.d0
      liqwt = 0.d0

      if (ippmod(iatmos).eq.2) then
        totwt = cvar_totwt(iel)
        liqwt = cpro_liqwt(iel)

        if (modsedi.eq.1.and.moddep.gt.0) then
          cpro_rugt(ifac) = rcodcl(ifac,iv,3)
        endif
      endif

      call atmcls &
      !==========
    ( ifac   ,                                                       &
      utau   , yplus  ,                                              &
      uet    ,                                                       &
      gredu  ,                                                       &
      cfnnu  , cfnns  , cfnnk  , cfnne  ,                            &
      temp   , totwt  , liqwt  ,                                     &
      icodcl , rcodcl )

    endif

    if (iwallf.le.2) then
      uk = uet
    endif

    uetmax = max(uet,uetmax)
    uetmin = min(uet,uetmin)
    ukmax  = max(uk,ukmax)
    ukmin  = min(uk,ukmin)
    yplumx = max(yplus,yplumx)
    yplumn = min(yplus,yplumn)

    ! save turbulent subgrid viscosity after van Driest damping in LES
    ! care is taken to not dampen it twice at boundary cells having more
    ! one boundary face
    if (itytur.eq.4.and.idries.eq.1) then
      if (visvdr(iel).lt.-900.d0) then
        ! FIXME amortissement de van Driest a revoir en rugueux :
        ! visct(iel) = visct(iel)*(1.d0-exp(-yplus/cdries))**2
        visvdr(iel) = visct(iel)
        visctc      = visct(iel)
      endif
    endif

    ! Save yplus if post-processed
    if (iyplbr.ge.0) then
      yplbr(ifac) = yplus
    endif

    !===========================================================================
    ! 3. Velocity boundary conditions
    !===========================================================================

    ! uiptn respecte la production de k
    !  de facon conditionnelle --> Coef RCPROD

    ! --> All turbulence models (except v2f and EBRSM)
    !-------------------------------------------------
    if (itytur.eq.2 .or. iturb.eq.60 .or.        &
         iturb.eq.0 .or. iturb.eq.10 .or.        &
         iturb.eq.30.or. iturb.eq.31 .or.        &
        itytur.eq.4 .or.                         &
         iturb.eq.70        ) then

      if (visctc.gt.epzero) then

        ! Pseudo decalage de la paroi de la distance rugd :
        distb0=distbf+rugd
        xmutlm = xkappa*uk*distb0*romc

        rcprod = distbf/distb0*max(1.d0,                                &
                     2.d0*sqrt(xmutlm/visctc) - distb0/distbf/(2.d0+rugd/distb0))

        rcflux = max(xmutlm,visctc)/(visclc+visctc)*distbf/distb0

        ! modified for non neutral boundary layer (cfnnu)
        ! FIXME: it should be uet/xkappa and not uk!
        uiptn  = min(utau,max(utau - uk/xkappa*rcprod*cfnnu,0.d0))
        iuntur = 1

        ! Coupled solving of the velocity components
        ! The boundary term for velocity gradient is implicit
        cofimp  = max(1.d0 - 1.d0/(xkappa*uplus)*rcprod*cfnnu, 0.d0)
        ! The term (rho*uet*uk) is implicit
        hflui = (visclc+visctc)/distbf/(xkappa*uplus)*rcflux*cfnnu

      ! In the viscous sub-layer
      else
        uiptn  = 0.d0
        iuntur = 0

        ! Coupled solving of the velocity components
        cofimp  = 0.d0
        hflui = visclc / distbf

      endif

      ! Clipping :
      ! On borne U_f,grad entre 0 et Utau (il y a surement mieux...)
      ! - 0    : on interdit le retournement en face de bord, qui est en
      !          contradiction avec l'hypoth\E8se de loi log.
      ! - Utau : la production turbulente ne peut etre nulle
      ! On empeche U_f,flux d'etre negatif

    ! --> v2f and EBRSM !FIXME EBRSM
    !------------------
    elseif (itytur.eq.5) then

      ! Avec ces conditions, pas besoin de calculer uiptmx, uiptmn
      ! et iuiptn qui sont nuls (valeur d'initialisation)
      iuntur = 0
      uiptn  = 0.d0

      ! Coupled solving of the velocity components
      hflui = (visclc + visctc) / distbf
      cofimp = 0.d0

    endif

    ! Min and Max and counter of reversal layer
    uiptmn = min(uiptn*iuntur,uiptmn)
    uiptmx = max(uiptn*iuntur,uiptmx)
    if (uiptn*iuntur.lt.-epzero) iuiptn = iuiptn + 1

    if (itytur.eq.3) then
      hint =  visclc          /distbf
    else
      hint = (visclc + visctc)/distbf
    endif

    ! Coupled solving of the velocity components

    ! Gradient boundary conditions
    !-----------------------------
    rcodcn = rcodcx*rnx+rcodcy*rny+rcodcz*rnz

    coefau(1,ifac) = (1.d0-cofimp)*(rcodcx - rcodcn*rnx) + rcodcn*rnx
    coefau(2,ifac) = (1.d0-cofimp)*(rcodcy - rcodcn*rny) + rcodcn*rny
    coefau(3,ifac) = (1.d0-cofimp)*(rcodcz - rcodcn*rnz) + rcodcn*rnz

    ! Projection in order to have the velocity parallel to the wall
    ! B = cofimp * ( IDENTITY - n x n )

    coefbu(1,1,ifac) = cofimp*(1.d0-rnx**2)
    coefbu(2,2,ifac) = cofimp*(1.d0-rny**2)
    coefbu(3,3,ifac) = cofimp*(1.d0-rnz**2)
    coefbu(1,2,ifac) = -cofimp*rnx*rny
    coefbu(1,3,ifac) = -cofimp*rnx*rnz
    coefbu(2,1,ifac) = -cofimp*rny*rnx
    coefbu(2,3,ifac) = -cofimp*rny*rnz
    coefbu(3,1,ifac) = -cofimp*rnz*rnx
    coefbu(3,2,ifac) = -cofimp*rnz*rny

    ! Flux boundary conditions
    !-------------------------

    cofafu(1,ifac)   = -hflui*(rcodcx - rcodcn*rnx) - hint*rcodcn*rnx
    cofafu(2,ifac)   = -hflui*(rcodcy - rcodcn*rny) - hint*rcodcn*rny
    cofafu(3,ifac)   = -hflui*(rcodcz - rcodcn*rnz) - hint*rcodcn*rnz

    ! Projection in order to have the shear stress parallel to the wall
    !  B = hflui*( IDENTITY - n x n )

    cofbfu(1,1,ifac) = hflui*(1.d0-rnx**2) + hint*rnx**2
    cofbfu(2,2,ifac) = hflui*(1.d0-rny**2) + hint*rny**2
    cofbfu(3,3,ifac) = hflui*(1.d0-rnz**2) + hint*rnz**2

    cofbfu(1,2,ifac) = (hint - hflui)*rnx*rny
    cofbfu(1,3,ifac) = (hint - hflui)*rnx*rnz
    cofbfu(2,1,ifac) = (hint - hflui)*rny*rnx
    cofbfu(2,3,ifac) = (hint - hflui)*rny*rnz
    cofbfu(3,1,ifac) = (hint - hflui)*rnz*rnx
    cofbfu(3,2,ifac) = (hint - hflui)*rnz*rny

    ! In case of transient turbomachinery computations, save the coefficents
    ! associated to rough wall velocity BC, in order to update the wall velocity
    ! after the geometry update (between prediction and correction step)
    if (iturbo.eq.2 .and. irotce(iel).ne.0) then
      coftur(ifac) = cofimp
      hfltur(ifac) = hflui
    endif

    !===========================================================================
    ! 4. Boundary conditions on k and epsilon
    !===========================================================================

    ydep = distbf*0.5d0+rugd

    if (itytur.eq.2) then

      ! Dirichlet Boundary Condition on k
      !----------------------------------

      pimp = uk**2/sqrcmu*cfnnk
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_k(ifac), coefaf_k(ifac),             &
           coefb_k(ifac), coefbf_k(ifac),             &
           pimp         , hint          , rinfin )


      ! No diffusion reconstruction when using wall functions
      if (associated(cpro_diff_lim_k)) cpro_diff_lim_k(ifac) = 0.d0

      ! Neumann Boundary Condition on epsilon
      !--------------------------------------

      hint = (visclc+visctc/sigmae)/distbf

      pimp = uk**3/(xkappa*ydep**2)*distbf*cfnne
      qimp = -pimp*hint !TODO transform it to use d eps / d y directly

      call set_neumann_scalar &
           !==================
         ( coefa_ep(ifac), coefaf_ep(ifac),             &
           coefb_ep(ifac), coefbf_ep(ifac),             &
           qimp          , hint )

      ! No diffusion reconstruction when using wall functions
      if (associated(cpro_diff_lim_eps)) cpro_diff_lim_eps(ifac) = 0.d0

    !===========================================================================
    ! 5. Boundary conditions on Rij-epsilon
    !===========================================================================

    elseif (itytur.eq.3) then

      ! Exchange coefficient

      ! Symmetric tensor diffusivity (Daly Harlow -- GGDH)
      if (iand(vcopt_rij%idften, ANISOTROPIC_RIGHT_DIFFUSION).ne.0) then

        visci(1,1) = visclc + visten(1,iel)
        visci(2,2) = visclc + visten(2,iel)
        visci(3,3) = visclc + visten(3,iel)
        visci(1,2) =          visten(4,iel)
        visci(2,1) =          visten(4,iel)
        visci(2,3) =          visten(5,iel)
        visci(3,2) =          visten(5,iel)
        visci(1,3) =          visten(6,iel)
        visci(3,1) =          visten(6,iel)

        ! ||Ki.S||^2
        viscis = ( visci(1,1)*surfbo(1,ifac)       &
                 + visci(1,2)*surfbo(2,ifac)       &
                 + visci(1,3)*surfbo(3,ifac))**2   &
               + ( visci(2,1)*surfbo(1,ifac)       &
                 + visci(2,2)*surfbo(2,ifac)       &
                 + visci(2,3)*surfbo(3,ifac))**2   &
               + ( visci(3,1)*surfbo(1,ifac)       &
                 + visci(3,2)*surfbo(2,ifac)       &
                 + visci(3,3)*surfbo(3,ifac))**2

        ! IF.Ki.S
        fikis = ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,1)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,1)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,1)   &
                )*surfbo(1,ifac)                              &
              + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,2)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,2)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,2)   &
                )*surfbo(2,ifac)                              &
              + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,3)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,3)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,3)   &
                )*surfbo(3,ifac)

        distfi = distb(ifac)

        ! Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
        ! NB: eps =1.d-1 must be consistent with vitens.f90
        fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

        hint = viscis/surfbn(ifac)/fikis

      ! Scalar diffusivity
      else
        hint = (visclc+visctc*csrij/cmu)/distbf
      endif

      ! ---> Tensor Rij (Partially implicited)

      do isou = 1, 6
        fcoefa(isou) = 0.0d0
        fcoefb(isou) = 0.0d0
        fcofad(isou) = 0.0d0
        fcofbd(isou) = 0.0d0
      enddo

      do isou = 1, 6

        if (isou.eq.1) then
          jj = 1
          kk = 1
        else if (isou.eq.2) then
          jj = 2
          kk = 2
        else if (isou.eq.3) then
          jj = 3
          kk = 3
        else if (isou.eq.4) then
          jj = 1
          kk = 2
        else if (isou.eq.5) then
          jj = 2
          kk = 3
        else if (isou.eq.6) then
          jj = 1
          kk = 3
        endif

        ! Coupled version: we implicit as much as possible
        if (irijco.eq.1) then
          coefa_rij(isou, ifac) = - (  eloglo(jj,1)*eloglo(kk,2)         &
                                     + eloglo(jj,2)*eloglo(kk,1))        &
                                    * alpha_rnn * sqrt(rnnb * rttb) * cfnnk
          coefaf_rij(isou, ifac)      = -hint * coefa_rij(isou, ifac)
          coefad_rij(isou, ifac)      = 0.d0
          do ii = 1, 6
            coefb_rij(isou,ii, ifac)  = alpha(ii,isou)
            if (ii.eq.isou) then
              coefbf_rij(isou,ii, ifac) = hint &
                                        * (1.d0 - coefb_rij(isou,ii, ifac))
            else
              coefbf_rij(isou,ii, ifac) = - hint &
                                        * coefb_rij(isou,ii, ifac)
            endif
            coefbd_rij(isou,ii, ifac) = coefb_rij(isou,ii, ifac)
          enddo

        else if (iclptr.eq.1) then

          do ii = 1, 6
            if (ii.ne.isou) then
              fcoefa(isou) = fcoefa(isou) + alpha(isou,ii) * rijipb(ifac,ii)
            endif
          enddo
          fcoefb(isou) = alpha(isou,isou)
        else
          do ii = 1, 6
            fcoefa(isou) = fcoefa(isou) + alpha(isou,ii) * rijipb(ifac,ii)
          enddo
          fcoefb(isou) = 0.d0
        endif

        ! Boundary conditions for the momentum equation
        fcofad(isou) = fcoefa(isou)
        fcofbd(isou) = fcoefb(isou)

        fcoefa(isou) = fcoefa(isou)                                 &
                       - (  eloglo(jj,1)*eloglo(kk,2)               &
                          + eloglo(jj,2)*eloglo(kk,1))*uet*uk*cfnnk

        ! Translate into Diffusive flux BCs
        fcofaf(isou) = -hint*fcoefa(isou)
        fcofbf(isou) = hint*(1.d0-fcoefb(isou))
      enddo

      if (irijco.ne.1) then
        do isou = 1, 6
          if (isou.eq.1) then
            coefa_r11(ifac) = fcoefa(isou)
            coefb_r11(ifac) = fcoefb(isou)
            coefaf_r11(ifac) = fcofaf(isou)
            coefbf_r11(ifac) = fcofbf(isou)
            coefad_r11(ifac) = fcofad(isou)
            coefbd_r11(ifac) = fcofbd(isou)
          else if (isou.eq.2) then
            coefa_r22(ifac) = fcoefa(isou)
            coefb_r22(ifac) = fcoefb(isou)
            coefaf_r22(ifac) = fcofaf(isou)
            coefbf_r22(ifac) = fcofbf(isou)
            coefad_r22(ifac) = fcofad(isou)
            coefbd_r22(ifac) = fcofbd(isou)
          else if (isou.eq.3) then
            coefa_r33(ifac) = fcoefa(isou)
            coefb_r33(ifac) = fcoefb(isou)
            coefaf_r33(ifac) = fcofaf(isou)
            coefbf_r33(ifac) = fcofbf(isou)
            coefad_r33(ifac) = fcofad(isou)
            coefbd_r33(ifac) = fcofbd(isou)
          else if (isou.eq.4) then
            coefa_r12(ifac) = fcoefa(isou)
            coefb_r12(ifac) = fcoefb(isou)
            coefaf_r12(ifac) = fcofaf(isou)
            coefbf_r12(ifac) = fcofbf(isou)
            coefad_r12(ifac) = fcofad(isou)
            coefbd_r12(ifac) = fcofbd(isou)
          else if (isou.eq.5) then
            coefa_r23(ifac) = fcoefa(isou)
            coefb_r23(ifac) = fcoefb(isou)
            coefaf_r23(ifac) = fcofaf(isou)
            coefbf_r23(ifac) = fcofbf(isou)
            coefad_r23(ifac) = fcofad(isou)
            coefbd_r23(ifac) = fcofbd(isou)
          else if (isou.eq.6) then
            coefa_r13(ifac) = fcoefa(isou)
            coefb_r13(ifac) = fcoefb(isou)
            coefaf_r13(ifac) = fcofaf(isou)
            coefbf_r13(ifac) = fcofbf(isou)
            coefad_r13(ifac) = fcofad(isou)
            coefbd_r13(ifac) = fcofbd(isou)
          endif
        enddo
      endif

      ! No diffusion reconstruction when using wall functions
      if (associated(cpro_diff_lim_rij)) cpro_diff_lim_rij(ifac) = 0.d0

      ! ---> Epsilon

      ! Symmetric tensor diffusivity (Daly Harlow -- GGDH)
      if (iand(vcopt_ep%idften, ANISOTROPIC_DIFFUSION).ne.0) then

        visci(1,1) = visclc + visten(1,iel)/sigmae
        visci(2,2) = visclc + visten(2,iel)/sigmae
        visci(3,3) = visclc + visten(3,iel)/sigmae
        visci(1,2) =          visten(4,iel)/sigmae
        visci(2,1) =          visten(4,iel)/sigmae
        visci(2,3) =          visten(5,iel)/sigmae
        visci(3,2) =          visten(5,iel)/sigmae
        visci(1,3) =          visten(6,iel)/sigmae
        visci(3,1) =          visten(6,iel)/sigmae

        ! ||Ki.S||^2
        viscis = ( visci(1,1)*surfbo(1,ifac)       &
                 + visci(1,2)*surfbo(2,ifac)       &
                 + visci(1,3)*surfbo(3,ifac))**2   &
               + ( visci(2,1)*surfbo(1,ifac)       &
                 + visci(2,2)*surfbo(2,ifac)       &
                 + visci(2,3)*surfbo(3,ifac))**2   &
               + ( visci(3,1)*surfbo(1,ifac)       &
                 + visci(3,2)*surfbo(2,ifac)       &
                 + visci(3,3)*surfbo(3,ifac))**2

        ! IF.Ki.S
        fikis = ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,1)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,1)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,1)   &
                )*surfbo(1,ifac)                              &
              + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,2)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,2)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,2)   &
                )*surfbo(2,ifac)                              &
              + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,3)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,3)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,3)   &
                )*surfbo(3,ifac)

        distfi = distb(ifac)

        ! Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
        ! NB: eps =1.d-1 must be consistent with vitens.f90
        fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

        hint = viscis/surfbn(ifac)/fikis

      ! Scalar diffusivity
      else
        hint = (visclc+visctc/sigmae)/distbf
      endif

      ! Neumann Boundary Condition on epsilon
      !--------------------------------------

      pimp = uk**3/(xkappa*ydep**2)*distbf*cfnne
      qimp = -pimp*hint !TODO transform it to use d eps / d y directly

      call set_neumann_scalar &
           !==================
         ( coefa_ep(ifac), coefaf_ep(ifac),             &
           coefb_ep(ifac), coefbf_ep(ifac),             &
           qimp          , hint )

      ! No diffusion reconstruction when using wall functions
      if (associated(cpro_diff_lim_eps)) cpro_diff_lim_eps(ifac) = 0.d0

    !===========================================================================
    ! 6a.Boundary conditions on k, epsilon, f_bar and phi in the phi_Fbar model
    !===========================================================================

    elseif (iturb.eq.50) then

      ! Dirichlet Boundary Condition on k
      !----------------------------------

      pimp = 0.d0
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_k(ifac), coefaf_k(ifac),             &
           coefb_k(ifac), coefbf_k(ifac),             &
           pimp         , hint          , rinfin )

      ! Dirichlet Boundary Condition on epsilon
      !----------------------------------------

      pimp = 2.0d0*visclc/romc*cvar_k(iel)/distbf**2
      hint = (visclc+visctc/sigmae)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_ep(ifac), coefaf_ep(ifac),             &
           coefb_ep(ifac), coefbf_ep(ifac),             &
           pimp          , hint          , rinfin )

      ! Dirichlet Boundary Condition on Phi
      !------------------------------------

      pimp = 0.d0
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_phi(ifac), coefaf_phi(ifac),             &
           coefb_phi(ifac), coefbf_phi(ifac),             &
           pimp           , hint          , rinfin )

      ! Dirichlet Boundary Condition on Fb
      !-----------------------------------

      pimp = 0.d0
      hint = 1.d0/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_fb(ifac), coefaf_fb(ifac),             &
           coefb_fb(ifac), coefbf_fb(ifac),             &
           pimp          , hint           , rinfin )

    !===========================================================================
    ! 6b.Boundary conditions on k, epsilon, phi and alpha in the Bl-v2/k model
    !===========================================================================

    elseif (iturb.eq.51) then

      ! Dirichlet Boundary Condition on k
      !----------------------------------

      pimp = 0.d0
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_k(ifac), coefaf_k(ifac),             &
           coefb_k(ifac), coefbf_k(ifac),             &
           pimp         , hint          , rinfin )

      ! Dirichlet Boundary Condition on epsilon
      !----------------------------------------

      pimp = visclc/romc*cvar_k(iel)/distbf**2
      hint = (visclc+visctc/sigmae)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_ep(ifac), coefaf_ep(ifac),             &
           coefb_ep(ifac), coefbf_ep(ifac),             &
           pimp          , hint           , rinfin )

      ! Dirichlet Boundary Condition on Phi
      !------------------------------------

      pimp = 0.d0
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_phi(ifac), coefaf_phi(ifac),             &
           coefb_phi(ifac), coefbf_phi(ifac),             &
           pimp           , hint            , rinfin )

      ! Dirichlet Boundary Condition on alpha
      !--------------------------------------

      pimp = 0.d0
      hint = 1.d0/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_al(ifac), coefaf_al(ifac),             &
           coefb_al(ifac), coefbf_al(ifac),             &
           pimp          , hint           , rinfin )


    !===========================================================================
    ! 7. Boundary conditions on k and omega
    !===========================================================================

    elseif (iturb.eq.60) then

      ! Always out of the viscous sub-layer

      ! Dirichlet Boundary Condition on k
      !----------------------------------

      pimp = uk**2/sqrcmu

      !FIXME it is wrong because sigma is computed within the model
      ! see turbkw.f90
      hint = (visclc+visctc/ckwsk2)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_k(ifac), coefaf_k(ifac),             &
           coefb_k(ifac), coefbf_k(ifac),             &
           pimp         , hint          , rinfin )

      ! Neumann Boundary Condition on omega
      !------------------------------------

      !FIXME it is wrong because sigma is computed within the model
      ! see turbkw.f90 (So the flux is not the one we impose!)
      hint = (visclc+visctc/ckwsw2)/distbf

      pimp = distbf*4.d0*uk**3*romc**2/           &
            (sqrcmu*xkappa*visclc**2*yplus**2)
      qimp = -pimp*hint !TODO transform it to use d eps / d y directly

      call set_neumann_scalar &
           !==================
         ( coefa_omg(ifac), coefaf_omg(ifac),             &
           coefb_omg(ifac), coefbf_omg(ifac),             &
           qimp           , hint )

    !===========================================================================
    ! 7.1 Boundary conditions on the Spalart Allmaras turbulence model
    !===========================================================================

    elseif (iturb.eq.70) then

      dsa0 = rugd
      hint = (visclc + vcopt%idifft*cvara_nusa(iel)*romc*dsa0/(distbf+dsa0) ) &
            / distbf / csasig

      ! If we have a rough wall then:
      ! nusa_wall*(1- I'F/d0)=nusa_I'
      ! which is a Robin type BC
      coefa_nu(ifac) = 0.d0
      coefb_nu(ifac) = dsa0/(dsa0+distbf)


      coefaf_nu(ifac) = 0.d0
      coefbf_nu(ifac) = hint*distbf/(dsa0+distbf)
    endif

    byplus(ifac) = yplus
    buk(ifac) = uk
    ustar(ifac) = uet
    bcfnns(ifac) = cfnns

  endif
  ! Test on the presence of a rough wall (End)

enddo
! --- End of loop over faces

!===========================================================================
! 8. Boundary conditions on the other scalars
!    (Specific treatment for the variances of the scalars next to walls:
!     see condli)
!===========================================================================

do iscal = 1, nscal

  if (iscavr(iscal).le.0) then

    call clptrg_scalar(iscal, isvhb, icodcl, rcodcl,              &
                       byplus, buk, ustar, bcfnns,                &
                       hbord, theipb,                             &
                       tetmax, tetmin, tplumx, tplumn)

  endif

enddo

if (irangp.ge.0) then
  call parmin (uiptmn)
  call parmax (uiptmx)
  call parmin (uetmin)
  call parmax (uetmax)
  call parmin (ukmin)
  call parmax (ukmax)
  call parmin (yplumn)
  call parmax (yplumx)
  call parcpt (iuiptn)
  if (iscalt.gt.0) then
    call parmin (tetmin)
    call parmax (tetmax)
    call parmin (tplumn)
    call parmax (tplumx)
  endif
endif

deallocate(byplus)
deallocate(buk)
if (allocated(buet)) deallocate(buet)
if (allocated(bcfnns_loc)) deallocate(bcfnns_loc)

!===============================================================================
! 9. Writings
!===============================================================================

!     Remarque : afin de ne pas surcharger les logs dans le cas ou
!       quelques yplus ne sont pas corrects, on ne produit le message
!       qu'aux deux premiers pas de temps ou le message apparait et
!       aux deux derniers pas de temps du calcul, ou si IWARNI est >= 2
!       On indique aussi le numero du dernier pas de temps auquel on
!       a rencontre des yplus hors bornes admissibles

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)

if (vcopt%iwarni.ge.0) then
  if (ntlist.gt.0) then
    modntl = mod(ntcabs,ntlist)
  elseif (ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
    modntl = 0
  else
    modntl = 1
  endif

  if ((modntl.eq.0 .or. vcopt%iwarni.ge.2).and.iscalt.gt.0) then
    write(nfecra,2011) &
         uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,   &
         tetmin, tetmax, tplumn, tplumx, iuiptn
  elseif (modntl.eq.0 .or. vcopt%iwarni.ge.2) then
    write(nfecra,2010) &
         uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,   &
         iuiptn
  endif

endif

!===============================================================================
! 10. Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(/,' LA NORMALE A LA FACE DE BORD DE PAROI ',I10,/,   &
         ' EST NULLE ; COORDONNEES : ',3E12.5)

 2010 format(/,                                                   &
 3X,'** CONDITIONS AUX LIMITES EN PAROI RUGUEUSE',/,        &
 '   -------------------------------------------',/,        &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Vitesse rel. en paroi    uiptn : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uet   : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uk    : ',2E12.5                 ,/,&
 '   Distance adim. rugueuse  yplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------   ',/,&
 '   Nbre de retournements de la vitesse en paroi : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)

 2011 format(/,                                                   &
 3X,'** CONDITIONS AUX LIMITES EN PAROI RUGUEUSE',/,        &
 '   -------------------------------------------',/,        &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Vitesse rel. en paroi    uiptn : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uet   : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uk    : ',2E12.5                 ,/,&
 '   Distance adim. rugueuse  yplus : ',2E12.5                 ,/,&
 '   Sca. thermal de frott.   tstar : ',2E12.5                 ,/,&
 '   Sca. thermal adim. rug.  tplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------'   ,/,&
 '   Nbre de retournements de la vitesse en paroi : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)


#else

 1000 format(/,' THE NORMAL TO THE WALL BOUNDARY FACE ',I10,/,    &
         ' IS NULL; COORDINATES: ',3E12.5)

 2010 format(/,                                                   &
 3X,'** BOUNDARY CONDITIONS FOR ROUGH WALLS',/,             &
 '   --------------------------------------',/,             &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Rel velocity at the wall uiptn : ',2E12.5                 ,/,&
 '   Friction velocity        uet   : ',2E12.5                 ,/,&
 '   Friction velocity        uk    : ',2E12.5                 ,/,&
 '   Rough dimensionless dist yplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------   ',/,&
 '   Nb of reversal of the velocity at the wall   : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)

 2011 format(/,                                                   &
 3X,'** BOUNDARY CONDITIONS FOR ROUGH WALLS',/,             &
 '   --------------------------------------',/,             &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Rel velocity at the wall uiptn : ',2E12.5                 ,/,&
 '   Friction velocity        uet   : ',2E12.5                 ,/,&
 '   Friction velocity        uk    : ',2E12.5                 ,/,&
 '   Rough dimensionless dist yplus : ',2E12.5                 ,/,&
 '   Friction thermal sca.    tstar : ',2E12.5                 ,/,&
 '   Rough dim-less th. sca.  tplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------   ',/,&
 '   Nb of reversal of the velocity at the wall   : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)


#endif

!----
! End
!----

return
end subroutine

!===============================================================================
! Local functions
!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iscal         scalar id
!> \param[in]     isvhb         indicator to save exchange coeffient
!> \param[in,out] icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!> \param[in]     byplus        dimensionless distance to the wall
!> \param[in]     buk           dimensionless velocity
!> \param[in]     buet          boundary ustar value
!> \param[in]     bcfnns        boundary correction factor
!> \param[in,out] hbord         exchange coefficient at boundary
!> \param[in]     theipb        boundary temperature in \f$ \centip \f$
!>                               (more exaclty the energetic variable)
!> \param[out]    tetmax        maximum local ustar value
!> \param[out]    tetmin        minimum local ustar value
!> \param[out]    tplumx        maximum local tplus value
!> \param[out]    tplumn        minimum local tplus value
!_______________________________________________________________________________

subroutine clptrg_scalar &
 ( iscal  , isvhb  , icodcl ,                                     &
   rcodcl ,                                                       &
   byplus , buk    , buet   , bcfnns ,                            &
   hbord  , theipb ,                                              &
   tetmax , tetmin , tplumx , tplumn )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use pointe
use entsor
use albase
use parall
use ppppar
use ppthch
use ppincl
use radiat
use cplsat
use mesh
use field
use lagran
use turbomachinery
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal, isvhb
integer, pointer, dimension(:,:) :: icodcl
double precision, pointer, dimension(:,:,:) :: rcodcl
double precision, dimension(:) :: byplus, buk, buet, bcfnns
double precision, pointer, dimension(:) :: hbord, theipb
double precision tetmax, tetmin, tplumx, tplumn

! Local variables

integer          ivar, f_id, b_f_id, isvhbl
integer          ifac, iel, isou, jsou
integer          ifcvsl, itplus, itstar

double precision cpp, rkl, prdtl, visclc, romc, tplus, cpscv
double precision distfi, distbf, fikis, hint, heq, hflui, hext
double precision yplus, phit, pimp, temp, tet, uk
double precision viscis, visctc, cofimp
double precision act, rugt
double precision rinfiv(3), pimpv(3)
double precision visci(3,3), hintt(6)
double precision turb_schmidt, exchange_coef

character(len=80) :: fname

double precision, dimension(:), pointer :: val_s, bval_s, crom, viscls
double precision, dimension(:), pointer :: viscl, visct, cpro_cp, cpro_cv

double precision, dimension(:), pointer :: bfconv, bhconv
double precision, dimension(:), pointer :: tplusp, tstarp
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:,:), pointer :: coefaut, cofafut, cofarut, visten
double precision, dimension(:,:,:), pointer :: coefbut, cofbfut, cofbrut

integer, save :: kbfid = -1

type(var_cal_opt) :: vcopt

!===============================================================================

ivar = isca(iscal)
f_id = ivarfl(ivar)

call field_get_val_s(ivarfl(ivar), val_s)

call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

call field_get_key_int (f_id, kivisl, ifcvsl)

if (ifcvsl .ge. 0) then
  call field_get_val_s(ifcvsl, viscls)
endif

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

if (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0.or.ityturt(iscal).eq.3) then
  if (iturb.ne.32.or.ityturt(iscal).eq.3) then
    call field_get_val_v(ivsten, visten)
  else ! EBRSM and (GGDH or AFM)
    call field_get_val_v(ivstes, visten)
  endif
endif

call field_get_coefa_s(f_id, coefap)
call field_get_coefb_s(f_id, coefbp)
call field_get_coefaf_s(f_id, cofafp)
call field_get_coefbf_s(f_id, cofbfp)

call field_get_val_s(icrom, crom)
if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
endif

if (ippmod(icompf) .ge. 0) then
  if (icv.ge.0) then
    call field_get_val_s(icv, cpro_cv)
  endif
endif

isvhbl = 0
if (iscal.eq.isvhb) then
  isvhbl = isvhb
endif

if (iscal.eq.iscalt) then
  ! min. and max. of wall friction of the thermal scalar
  tetmax = -grand
  tetmin =  grand
  ! min. and max. of T+
  tplumx = -grand
  tplumn =  grand
endif

rinfiv(1) = rinfin
rinfiv(2) = rinfin
rinfiv(3) = rinfin

if (ityturt(iscal).eq.3) then

  ! Name of the scalar ivar
  call field_get_name(ivarfl(ivar), fname)

  ! Index of the corresponding turbulent flux
  call field_get_id(trim(fname)//'_turbulent_flux', f_id)

  call field_get_coefa_v(f_id,coefaut)
  call field_get_coefb_v(f_id,coefbut)
  call field_get_coefaf_v(f_id,cofafut)
  call field_get_coefbf_v(f_id,cofbfut)
  call field_get_coefad_v(f_id,cofarut)
  call field_get_coefbd_v(f_id,cofbrut)

endif

! pointers to T+ and T* if saved

itplus = -1
itstar = -1

tplusp => null()
tstarp => null()

if (iscal.eq.iscalt) then
  call field_get_id_try('tplus', itplus)
  if (itplus.ge.0) then
    call field_get_val_s (itplus, tplusp)
  endif
  call field_get_id_try('tstar', itstar)
  if (itstar.ge.0) then
    call field_get_val_s (itstar, tstarp)
  endif
endif

! Pointers to specific fields

if (iirayo.ge.1 .and. iscal.eq.iscalt) then
  call field_get_val_s_by_name("rad_convective_flux", bfconv)
  call field_get_val_s_by_name("rad_exchange_coefficient", bhconv)
endif

if (kbfid.lt.0) call field_get_key_id("boundary_value_id", kbfid)

call field_get_key_int(f_id, kbfid, b_f_id)

if (b_f_id .ge. 0) then
  call field_get_val_s(b_f_id, bval_s)
else
  bval_s => null()
  ! if thermal variable has no boundary but temperature does, use it
  if (itherm.eq.2 .and. itempb.ge.0) then
    b_f_id = itempb
    call field_get_val_s(b_f_id, bval_s)
  endif
endif

! retrieve turbulent Schmidt value for current scalar
call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)

! --- Loop on boundary faces
do ifac = 1, nfabor

  yplus = byplus(ifac)
  uk = buk(ifac)

  ! Test on the presence of a rough wall condition (start)
  if (icodcl(ifac,iu).eq.6) then

    iel = ifabor(ifac)

    ! Physical quantities

    visclc = viscl(iel)
    visctc = visct(iel)
    romc   = crom(iel)

    ! Geometric quantities
    distbf = distb(ifac)

    cpp = 1.d0
    if (iscacp(iscal).eq.1) then
      if (icp.ge.0) then
        cpp = cpro_cp(iel)
      else
        cpp = cp0
      endif
    endif

    if (ifcvsl.lt.0) then
      rkl = visls0(iscal)
      prdtl = cpp*visclc/rkl
    else
      rkl = viscls(iel)
      prdtl = cpp*visclc/rkl
    endif

    ! --> Compressible module:
    ! On suppose que le nombre de Pr doit etre
    ! defini de la meme facon que l'on resolve
    ! en enthalpie ou en energie, soit Mu*Cp/Lambda.
    ! Si l'on resout en energie, on a calcule ci-dessus
    ! Mu*Cv/Lambda.

    if (iscal.eq.iscalt .and. itherm.eq.3) then
      if (icp.ge.0) then
        prdtl = prdtl*cpro_cp(iel)
      else
        prdtl = prdtl*cp0
      endif
      if (icv.ge.0) then
        prdtl = prdtl/cpro_cv(iel)
      else
        prdtl = prdtl/cv0
      endif
    endif

    ! Scalar diffusivity
    if (iand(vcopt%idften, ISOTROPIC_DIFFUSION).ne.0) then
      ! En compressible, pour l'energie LAMBDA/CV+CP/CV*(MUT/TURB_SCHMIDT)
      if (ippmod(icompf) .ge. 0) then
        if (icp.ge.0) then
          cpscv = cpro_cp(iel)
        else
          cpscv = cp0
        endif
        if (icv.ge.0) then
          cpscv = cpscv/cpro_cv(iel)
        else
          cpscv = cpscv/cv0
        endif
        hint = (rkl+vcopt%idifft*cpp*cpscv*visctc/turb_schmidt)/distbf
      else
        hint = (rkl+vcopt%idifft*cpp*visctc/turb_schmidt)/distbf
      endif

      ! Symmetric tensor diffusivity (GGDH or AFM)
    elseif (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then
      ! En compressible, pour l'energie LAMBDA/CV+CP/CV*(MUT/SIGMAS)
      if (ippmod(icompf) .ge. 0) then
        if (icp.ge.0) then
          cpscv = cpro_cp(iel)
        else
          cpscv = cp0
        endif
        if (icv.ge.0) then
          cpscv = cpscv/cpro_cv(iel)
        else
          cpscv = cpscv/cv0
        endif
        temp = vcopt%idifft*cpp*cpscv*ctheta(iscal)
      else
        temp = vcopt%idifft*cpp*ctheta(iscal)
      endif
      visci(1,1) = temp*visten(1,iel) + rkl
      visci(2,2) = temp*visten(2,iel) + rkl
      visci(3,3) = temp*visten(3,iel) + rkl
      visci(1,2) = temp*visten(4,iel)
      visci(2,1) = temp*visten(4,iel)
      visci(2,3) = temp*visten(5,iel)
      visci(3,2) = temp*visten(5,iel)
      visci(1,3) = temp*visten(6,iel)
      visci(3,1) = temp*visten(6,iel)

      ! ||Ki.S||^2
      viscis =   ( visci(1,1)*surfbo(1,ifac)       &
                 + visci(1,2)*surfbo(2,ifac)       &
                 + visci(1,3)*surfbo(3,ifac))**2   &
               + ( visci(2,1)*surfbo(1,ifac)       &
                 + visci(2,2)*surfbo(2,ifac)       &
                 + visci(2,3)*surfbo(3,ifac))**2   &
               + ( visci(3,1)*surfbo(1,ifac)       &
                 + visci(3,2)*surfbo(2,ifac)       &
                 + visci(3,3)*surfbo(3,ifac))**2

      ! IF.Ki.S
      fikis = ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,1)   &
              + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,1)   &
              + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,1)   &
              )*surfbo(1,ifac)                              &
            + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,2)   &
              + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,2)   &
              + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,2)   &
              )*surfbo(2,ifac)                              &
            + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,3)   &
              + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,3)   &
              + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,3)   &
              )*surfbo(3,ifac)

      distfi = distb(ifac)

      ! Take I" so that I"F= eps*||FI||*Ki.n when I" is not in cell i
      ! NB: eps =1.d-1 must be consistent with vitens.f90
      fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

      hint = viscis/surfbn(ifac)/fikis
    endif

    ! Dirichlet on the scalar, with wall function
    if (iturb.ne.0.and.icodcl(ifac,ivar).eq.6) then
      rugt = rcodcl(ifac,iv,3)
      act = xkappa/log((distbf+rugt)/rugt)
      ! T+ = (T_I - T_w) / Tet
      tplus = log((distbf+rugt)/rugt)/xkappa
      hflui = romc*cpp*buet(ifac)*act*bcfnns(ifac)
      ! Neumann on the scalar, with wall function (for post-processing)
    else
      hflui = hint
      ! Tplus is chosen for post-processing
      tplus = romc*cpp*max(buet(ifac)*bcfnns(ifac),epzero)/hint
    endif

    hext = rcodcl(ifac,ivar,2)
    pimp = rcodcl(ifac,ivar,1)

    if (abs(hext).gt.rinfin*0.5d0) then
      heq = hflui
    else
      heq = hflui*hext/(hflui+hext)
    endif

    ! --->  C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE

    ! Si on a deux types de conditions aux limites (ICLVAR, ICLVAF)
    !   il faut que le flux soit traduit par ICLVAF.
    ! Si on n'a qu'un type de condition, peu importe (ICLVAF=ICLVAR)
    ! Pour le moment, dans cette version compressible, on impose un
    !   flux nul pour ICLVAR, lorsqu'il est different de ICLVAF (cette
    !   condition ne sert qu'a la reconstruction des gradients et
    !   s'applique a l'energie totale qui inclut l'energie cinetique :

    if (icodcl(ifac,ivar).eq.6) then

      cofimp = 1.d0 - heq/hint

      ! To be coherent with a wall function, clip it to 0
      cofimp = max(cofimp, 0.d0)

      ! Gradient BCs
      coefap(ifac) = (1.d0 - cofimp)*pimp
      coefbp(ifac) = cofimp

      ! Flux BCs
      cofafp(ifac) = -heq*pimp
      cofbfp(ifac) =  heq

      !--> Turbulent heat flux
      if (ityturt(iscal).eq.3) then

        phit = (cofafp(ifac) + cofbfp(ifac)*val_s(iel))

        hintt(1) =   0.5d0*(visclc+rkl)/distbf                        &
                   + visten(1,iel)*ctheta(iscal)/distbf/csrij
        hintt(2) =   0.5d0*(visclc+rkl)/distbf                        &
                   + visten(2,iel)*ctheta(iscal)/distbf/csrij
        hintt(3) =   0.5d0*(visclc+rkl)/distbf                        &
                   + visten(3,iel)*ctheta(iscal)/distbf/csrij
        hintt(4) = visten(4,iel)*ctheta(iscal)/distbf/csrij
        hintt(5) = visten(5,iel)*ctheta(iscal)/distbf/csrij
        hintt(6) = visten(6,iel)*ctheta(iscal)/distbf/csrij

        ! Dirichlet Boundary Condition
        !-----------------------------

        ! Add rho*uk*Tet to T'v' in High Reynolds
        do isou = 1, 3
          pimpv(isou) = surfbo(isou,ifac)*phit/(surfbn(ifac)*cpp*romc)
        enddo

        call set_dirichlet_vector_aniso &
             !========================
           ( coefaut(:,ifac)  , cofafut(:,ifac)  ,           &
             coefbut(:,:,ifac), cofbfut(:,:,ifac),           &
             pimpv            , hintt            , rinfiv )

        ! Boundary conditions used in the temperature equation
        do isou = 1, 3
          cofarut(isou,ifac) = 0.d0
          do jsou = 1, 3
            cofbrut(isou,jsou,ifac) = 0.d0
          enddo
        enddo

      endif

      ! Storage of the thermal exchange coefficient
      ! (conversion in case of energy or enthaly)
      ! the exchange coefficient is in W/(m2 K)
      ! Usefull for thermal coupling or radiative transfer
      if (iirayo.ge.1 .and. iscal.eq.iscalt.or.isvhbl.gt.0) then
        ! Enthalpy
        if (itherm.eq.2) then
          ! If Cp is variable
          if (icp.ge.0) then
            exchange_coef = hflui*cpro_cp(iel)
          else
            exchange_coef = hflui*cp0
          endif

          ! Total energy (compressible module)
        elseif (itherm.eq.3) then
          ! If Cv is variable
          if (icv.ge.0) then
            exchange_coef = hflui*cpro_cv(iel)
          else
            exchange_coef = hflui*cv0
          endif

          ! Temperature
        elseif (iscacp(iscal).eq.1) then
          exchange_coef = hflui
        endif
      endif

      ! ---> Thermal coupling, store h = lambda/d
      if (isvhbl.gt.0) hbord(ifac) = exchange_coef

      ! ---> Radiative transfer
      if (iirayo.ge.1 .and. iscal.eq.iscalt) then
        bhconv(ifac) = exchange_coef

        ! The outgoing flux is stored (Q = h(Ti'-Tp): negative if
        !  gain for the fluid) in W/m2
        bfconv(ifac) = cofafp(ifac) + cofbfp(ifac)*theipb(ifac)
      endif

    endif ! End if icodcl=6

    ! Save the value of T^star and T^+ for post-processing

    if (b_f_id.ge.0 .or. iscal.eq.iscalt) then

      ! Rough wall function
      if (icodcl(ifac,ivar).eq.6) then
        phit = cofafp(ifac)+cofbfp(ifac)*theipb(ifac)
        ! Imposed flux with wall function for post-processing
      elseif (icodcl(ifac,ivar).eq.3) then
        phit = rcodcl(ifac,ivar,3)
      else
        phit = 0.d0
      endif

      tet = phit/(romc*cpp*max(buet(ifac)*bcfnns(ifac),epzero))

      if (b_f_id .ge. 0) bval_s(ifac) = bval_s(ifac) - tplus*tet

      if (itplus .ge. 0) tplusp(ifac) = tplus
      if (itstar .ge. 0) tstarp(ifac) = tet

      if (iscal.eq.iscalt) then
        tetmax = max(tet, tetmax)
        tetmin = min(tet, tetmin)
        tplumx = max(tplus,tplumx)
        tplumn = min(tplus,tplumn)
      endif

    endif

  endif ! rough wall condition

enddo

!----
! End
!----

return
end subroutine
