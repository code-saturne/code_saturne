!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file clptur.f90
!>
!> \brief Boundary conditions for smooth walls (icodcl = 5).
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
!> - for a vector field such as the veclocity \f$ \vect{u} \f$ the boundary
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
!> Please refer to the
!> <a href="../../theory.pdf#wallboundary"><b>wall boundary conditions</b></a>
!> section of the theory guide for more informations, as well as the
!> <a href="../../theory.pdf#clptur"><b>clptur</b></a> section.
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
!> \param[out]    visvdr        dynamic viscosity after V. Driest damping in
!>                               boundary cells
!> \param[out]    hbord         exchange coefficient at boundary
!> \param[in]     theipb        value of thermal scalar at \f$ \centip \f$
!>                               of boundary cells
!_______________________________________________________________________________

subroutine clptur &
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
use dimens, only: nvar
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

integer          nscal, isvhb

integer          icodcl(nfabor,nvar)

double precision rcodcl(nfabor,nvar,3)
double precision velipb(nfabor,ndim), rijipb(nfabor,6)
double precision visvdr(ncelet)
double precision hbord(nfabor),theipb(nfabor)

! Local variables

integer          ifac, iel, isou, ii, jj, kk
integer          iscal
integer          modntl
integer          iuntur, f_dim
integer          nlogla, nsubla, iuiptn
integer          f_id_rough, f_id

double precision rnx, rny, rnz, rxnn
double precision tx, ty, tz, txn, txn0, t2x, t2y, t2z
double precision utau, upx, upy, upz, usn
double precision uiptn, uiptmn, uiptmx
double precision uetmax, uetmin, ukmax, ukmin, yplumx, yplumn
double precision tetmax, tetmin, tplumx, tplumn
double precision uk, uet, nusury, yplus, dplus
double precision sqrcmu, clsyme, ek
double precision xnuii, xnuit, xmutlm
double precision rcprod
double precision hflui, hint, pimp, qimp
double precision eloglo(3,3), alpha(6,6)
double precision rcodcx, rcodcy, rcodcz, rcodcn
double precision visclc, visctc, romc  , distbf, srfbnf, efvisc
double precision cofimp, ypup
double precision bldr12
double precision xkip
double precision rinfiv(3)
double precision visci(3,3), fikis, viscis, distfi
double precision fcoefa(6), fcoefb(6), fcofaf(6), fcofbf(6), fcofad(6), fcofbd(6)
double precision rxx, rxy, rxz, ryy, ryz, rzz, rnnb
double precision rttb, alpha_rnn
double precision roughness

double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: viscl, visct, cpro_cp, yplbr, uetbor
double precision, dimension(:), pointer :: bfpro_roughness
double precision, dimension(:), allocatable :: byplus, bdplus, buk
double precision, dimension(:), pointer :: cvar_k, cvar_ep
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:,:), pointer :: cvar_rij

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

double precision, dimension(:,:), pointer :: coefa_rij, coefaf_rij, coefad_rij
double precision, dimension(:,:,:), pointer :: coefb_rij, coefbf_rij, coefbd_rij

double precision  pimp_lam, pimp_turb, gammap

integer          ntlast , iaff
data             ntlast , iaff /-1 , 0/
save             ntlast , iaff

type(var_cal_opt) :: vcopt_u, vcopt_rij, vcopt_ep

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

! Initialize variables to avoid compiler warnings

cofimp  = 0.d0
ek = 0.d0
rcprod = 0.d0
uiptn = 0.d0
rnnb = 0.d0

rinfiv(1) = rinfin
rinfiv(2) = rinfin
rinfiv(3) = rinfin

uet = 1.d0
utau = 1.d0

! --- Constants
sqrcmu = sqrt(cmu)

yplbr => null()
bfpro_roughness => null()

call field_get_id_try("boundary_roughness", f_id_rough)
if (f_id_rough.ge.0) call field_get_val_s(f_id_rough, bfpro_roughness)

if (iyplbr.ge.0) call field_get_val_s(iyplbr, yplbr)
if (itytur.eq.3 .and. idirsm.eq.1) call field_get_val_v(ivsten, visten)

uetbor => null()

if (     (itytur.eq.4 .and. idries.eq.1) &
    .or. (iilagr.ge.1 .and. idepst.gt.0) ) then
  call field_get_id_try('ustar', f_id)
  if (f_id.ge.0) then
    call field_get_val_s(f_id, uetbor)
  endif
endif

! --- Gradient and flux boundary conditions

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)
call field_get_coefaf_v(ivarfl(iu), cofafu)
call field_get_coefbf_v(ivarfl(iu), cofbfu)

if (ik.gt.0) then
  call field_get_coefa_s(ivarfl(ik), coefa_k)
  call field_get_coefb_s(ivarfl(ik), coefb_k)
  call field_get_coefaf_s(ivarfl(ik), coefaf_k)
  call field_get_coefbf_s(ivarfl(ik), coefbf_k)
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
else
  coefa_ep => null()
  coefb_ep => null()
  coefaf_ep => null()
  coefbf_ep => null()
endif

if (itytur.eq.3) then! Also have boundary conditions for the momentum equation

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
  call field_get_coefa_s(ivarfl(inusa), coefa_nu)
  call field_get_coefaf_s(ivarfl(inusa), coefaf_nu)
  call field_get_coefb_s(ivarfl(inusa), coefb_nu)
  call field_get_coefbf_s(ivarfl(inusa), coefbf_nu)
else
  coefa_nu => null()
  coefb_nu => null()
  coefaf_nu => null()
  coefbf_nu => null()
endif

! --- Physical quantities
call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)
if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
endif

if (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then
  call field_get_val_s(ivarfl(ik), cvar_k)
endif
if ((iturb.eq.30).or.(iturb.eq.31)) then
  call field_get_val_s(ivarfl(iep), cvar_ep)
endif

if (itytur.eq.3) then
  call field_get_key_struct_var_cal_opt(ivarfl(ir11), vcopt_rij)
  call field_get_key_struct_var_cal_opt(ivarfl(iep), vcopt_ep)
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

! Counters (turbulent, laminar, reversal, scale correction)
nlogla = 0
nsubla = 0
iuiptn = 0

! Alpha constant for a realisable BC for R12 with the SSG model
alpha_rnn = 0.47d0

! With v2f type model, (phi-fbar et BL-v2/k) u=0 is set directly, so
! uiptmx and uiptmn are necessarily 0
if (itytur.eq.5) then
  uiptmx = 0.d0
  uiptmn = 0.d0
endif

! Pointers to specific fields
allocate(byplus(nfabor))
allocate(bdplus(nfabor))
allocate(buk(nfabor))

! --- Loop on boundary faces
do ifac = 1, nfabor

  ! Test on the presence of a smooth wall condition (start)
  if (icodcl(ifac,iu).eq.5) then

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

    elseif (itytur.eq.3) then

      ! If the velocity is zero, vector T is normal and random;
      ! we need it for the reference change for Rij, and we cancel the velocity.

      txn0 = 0.d0

      if (abs(rny).ge.epzero.or.abs(rnz).ge.epzero)then
        rxnn = sqrt(rny**2+rnz**2)
        tx  =  0.d0
        ty  =  rnz/rxnn
        tz  = -rny/rxnn
      elseif (abs(rnx).ge.epzero.or.abs(rnz).ge.epzero)then
        rxnn = sqrt(rnx**2+rnz**2)
        tx  =  rnz/rxnn
        ty  =  0.d0
        tz  = -rnx/rxnn
      else
        write(nfecra,1000)ifac,rnx,rny,rnz
        call csexit (1)
      endif

    else

      ! If the velocity is zero, and we are not using Reynolds Stresses,
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

      ! --> Commpute alpha(6,6)

      ! Let f be the center of the boundary faces and
      !   I the center of the matching cell

      ! We noteE Rg (resp. Rl) indexed by f or by I
      !   the Reynolds Stress tensor in the global basis (resp. local)

      ! The alpha matrix applied to the global vector in I'
      !   (Rg11,I'|Rg22,I'|Rg33,I'|Rg12,I'|Rg13,I'|Rg23,I')t
      !    must provide the values to prescribe to the face
      !   (Rg11,f |Rg22,f |Rg33,f |Rg12,f |Rg13,f |Rg23,f )t
      !    except for the Dirichlet boundary conditions (added later)

      ! We define it by computing Rg,f as a function of Rg,I' as follows

      !   RG,f = ELOGLO.RL,f.ELOGLOt (matrix products)

      !                     | RL,I'(1,1)     B*U*.Uk     C*RL,I'(1,3) |
      !      with    RL,f = | B*U*.Uk       RL,I'(2,2)       0        |
      !                     | C*RL,I'(1,3)     0         RL,I'(3,3)   |

      !             with    RL,I = ELOGLOt.RG,I'.ELOGLO
      !                     B = 0
      !              and    C = 0 at the wall (1 with symmetry)

      ! We compute in fact  ELOGLO.projector.ELOGLOt

      clsyme=0.d0
      call clca66 (clsyme , eloglo , alpha)
      !==========

    endif

    !===========================================================================
    ! 2. Friction velocities
    !===========================================================================

    ! ---> Compute Uet depending if we are in the log zone or not
    !      in 1 or 2 velocity scales
    !      and uk based on ek

    nusury = visclc/(distbf*romc)
    xnuii = visclc/romc
    xnuit = visctc/romc

    if (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then
      ek = cvar_k(iel)
      ! TODO: we could add 2*nu_T dv/dy to rnnb
      rnnb = 2.d0 / 3.d0 * ek
    else if (itytur.eq.3) then
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

    if (f_id_rough.ge.0) then
      roughness = bfpro_roughness(ifac)
    else
      roughness = 0.d0
    endif

    call wallfunctions &
  ( iwallf, ifac  ,                                        &
    xnuii , xnuit , utau  , distbf, roughness, rnnb, ek,   &
    iuntur, nsubla, nlogla,                                &
    uet   , uk    , yplus , ypup  , cofimp, dplus )

    uetmax = max(uet,uetmax)
    uetmin = min(uet,uetmin)
    ukmax  = max(uk,ukmax)
    ukmin  = min(uk,ukmin)
    yplumx = max(yplus-dplus,yplumx)
    yplumn = min(yplus-dplus,yplumn)

    ! Sauvegarde de la vitesse de frottement et de la viscosite turbulente
    ! apres amortissement de van Driest pour la LES
    ! On n'amortit pas mu_t une seconde fois si on l'a deja fait
    ! (car une cellule peut avoir plusieurs faces de paroi)
    ! ou
    ! Sauvegarde de la vitesse de frottement et distance a la paroi yplus
    ! si le modele de depot de particules est active.

    if (itytur.eq.4.and.idries.eq.1) then
      uetbor(ifac) = uet !TODO remove, this information is in cofaf cofbf
      if (visvdr(iel).lt.-900.d0) then
        visct(iel) = visct(iel)*(1.d0-exp(-(yplus-dplus)/cdries))**2
        visvdr(iel) = visct(iel)
        visctc      = visct(iel)
      endif
    else if (iilagr.gt.0.and.idepst.gt.0) then
      uetbor(ifac) = uet
    endif

    ! Save yplus if post-processed or condensation modelling
    if (iyplbr.ge.0) then
      yplbr(ifac) = yplus-dplus
    endif

    !===========================================================================
    ! 3. Velocity boundary conditions
    !===========================================================================

    ! Deprecated power law (Werner & Wengle)
    if (iwallf.eq.1) then
      uiptn  = utau + uet*apow*bpow*yplus**bpow*(2.d0**(bpow-1.d0)-2.d0)

    ! Dependant on the turbulence Model
    else

      ! uiptn respecte la production de k
      !  de facon conditionnelle --> Coef RCPROD

      ! --> k-epsilon and k-omega
      !--------------------------
      if (itytur.eq.2.or.iturb.eq.60) then

        xmutlm = xkappa*visclc*yplus

        ! If yplus=0, uiptn is set to 0 to avoid division by 0.
        ! By the way, in this case: iuntur=0
        if (yplus.gt.epzero) then !TODO use iuntur.eq.1
          rcprod = min(xkappa , max(1.0d0,sqrt(xmutlm/visctc))/yplus)

          uiptn  = utau + distbf*uet*uk*romc/xkappa/visclc*(       &
               1.0d0/(2.0d0*yplus-dplus) - 2.0d0*rcprod )
        else
          uiptn = 0.d0
        endif

      ! --> No turbulence, mixing length or Rij-espilon
      !------------------------------------------------
      elseif (iturb.eq.0.or.iturb.eq.10.or.itytur.eq.3) then

        ! Dans le cadre de la ponderation elliptique, on ne doit pas
        ! tenir compte des lois de paroi. On fait donc un test sur le modele
        ! de turbulence :
        ! si on est en LRR ou SSG on laisse les lois de paroi, si on est en
        ! EBRSM, on impose l adherence.
        if (iturb.eq.32.or.iturb.eq.0) then
          uiptn = 0.d0
        else

          ! If yplus=0, uiptn is set to 0 to avoid division by 0.
          ! By the way, in this case: iuntur=0
          if (yplus.gt.epzero) then !FIXME use iuntur
            uiptn = utau - distbf*romc*uet*uk/xkappa/visclc                    &
                                 *(2.0d0/yplus - 1.0d0/(2.0d0*yplus-dplus))
          else
            uiptn = 0.d0
          endif

        endif

      ! --> LES and Spalart Allmaras
      !-----------------------------
      elseif (itytur.eq.4.or.iturb.eq.70) then

        uiptn  = utau - uet/xkappa*1.5d0

        ! If (mu+mut) becomes zero (dynamic models), an arbitrary value is set
        ! (nul flux) but without any problems because the flux
        ! is really zero at this face.
        if (visctc+visclc.le.0) then
          hflui = 0.d0 !FIXME

        endif

      ! --> v2f
      !--------
      elseif (itytur.eq.5) then

        ! Avec ces conditions, pas besoin de calculer uiptmx, uiptmn
        ! et iuiptn qui sont nuls (valeur d'initialisation)
        uiptn  = 0.d0

      endif
    endif

    ! Min and Max and counter of reversal layer
    uiptmn = min(uiptn*iuntur,uiptmn)
    uiptmx = max(uiptn*iuntur,uiptmx)
    if (uiptn*iuntur.lt.-epzero) iuiptn = iuiptn + 1

    ! To be coherent with a wall function, clip it to 0
    cofimp = max(cofimp, 0.d0)

    ! On implicite le terme (rho*uet*uk)
    hflui = visclc / distbf * ypup

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
    ! associated to turbulent wall velocity BC, in order to update the wall
    ! velocity after the geometry update(between prediction and correction step)
    if (iturbo.eq.2 .and. irotce(iel).ne.0) then
      coftur(ifac) = cofimp
      hfltur(ifac) = hflui
    endif

    !===========================================================================
    ! 4. Boundary conditions on k and espilon
    !===========================================================================

    if (itytur.eq.2) then

      ! Dirichlet Boundary Condition on k
      !----------------------------------

      if (iuntur.eq.1) then
        pimp = uk**2/sqrcmu
      else
        pimp = 0.d0
      endif
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa_k(ifac), coefaf_k(ifac),             &
           coefb_k(ifac), coefbf_k(ifac),             &
           pimp         , hint          , rinfin )


      ! Neumann Boundary Condition on epsilon
      !--------------------------------------

      hint = (visclc+visctc/sigmae)/distbf

      ! If yplus=0, uiptn is set to 0 to avoid division by 0.
      ! By the way, in this case: iuntur=0
      if (yplus.gt.epzero.and.iuntur.eq.1) then !FIXME use only iuntur
        efvisc = visclc/romc + exp(-xkappa*(8.5-5.2)) * roughness * uk
        pimp = distbf*4.d0*uk**5/           &
            (xkappa*efvisc**2*(yplus+dplus)**2)

        qimp = -pimp*hint !TODO transform it, it is only to be fully equivalent
      else
        qimp = 0.d0
      endif

      call set_neumann_scalar &
           !==================
         ( coefa_ep(ifac), coefaf_ep(ifac),             &
           coefb_ep(ifac), coefbf_ep(ifac),             &
           qimp          , hint )

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

      ! ---> Tensor Rij (Partially or totally implicited)

      do isou = 1, 6
        fcoefa(isou) = 0.0d0
        fcoefb(isou) = 0.0d0
        fcofad(isou) = 0.0d0
        fcofbd(isou) = 0.0d0
      enddo

      ! blending factor so that the component R(n,tau) have only
      ! -mu_T/(mu+mu_T)*uet*uk
      bldr12 = visctc/(visclc + visctc)

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

        ! LRR and the Standard SGG.
        if ((iturb.eq.30).or.(iturb.eq.31).and.iuntur.eq.1) then

          if (irijco.eq.1) then
            coefa_rij(isou, ifac) = - (  eloglo(jj,1)*eloglo(kk,2)         &
                                       + eloglo(jj,2)*eloglo(kk,1))        &
                                      * alpha_rnn * sqrt(rnnb * rttb)
            coefaf_rij(isou, ifac)      = -hint * coefa_rij(isou, ifac)
            coefad_rij(isou, ifac)      = 0.d0
            do ii = 1, 6
              coefb_rij(isou,ii, ifac)  = alpha(ii, isou)
              if (ii.eq.isou) then
                coefbf_rij(isou,ii, ifac) = hint * (1.d0 - coefb_rij(isou,ii, ifac))
              else
                coefbf_rij(isou,ii, ifac) = - hint * coefb_rij(isou,ii, ifac)
              endif
              coefbd_rij(isou,ii, ifac) = coefb_rij(isou,ii, ifac)
            enddo

          else if ((iclptr.eq.1)) then

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

          fcoefa(isou) = fcoefa(isou)                                   &
                         - (  eloglo(jj,1)*eloglo(kk,2)                 &
                            + eloglo(jj,2)*eloglo(kk,1))*bldr12*uet*uk

        ! In the viscous sublayer or for EBRSM: zero Reynolds' stresses
        else
          if (irijco.eq.1) then
            coefa_rij(isou, ifac) = 0.d0
            coefaf_rij(isou, ifac) = 0.d0
            coefad_rij(isou, ifac) = 0.d0
            do ii = 1, 6
              coefb_rij(isou,ii, ifac)  = 0.d0
              if (ii.eq.isou) then
                coefbf_rij(isou,ii, ifac) = hint
              else
                coefbf_rij(isou,ii, ifac) = 0.d0
              endif
              coefbd_rij(isou,ii, ifac) = 0.d0
            enddo

          else

            fcoefa(isou) = 0.d0
            fcofad(isou) = 0.d0
            fcoefb(isou) = 0.d0
            fcofbd(isou) = 0.d0
          endif

        endif

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

      ! ---> Epsilon
      !      NB: no reconstruction, possibility of partial implicitation

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

      if ((iturb.eq.30).or.(iturb.eq.31)) then

        ! Si yplus=0, on met coefa a 0 directement pour eviter une division
        ! par 0.
        if (yplus.gt.epzero.and.iuntur.eq.1) then
          efvisc = visclc/romc + exp(-xkappa*(8.5-5.2)) * roughness * uk
          pimp = distbf*4.d0*uk**5/           &
                (xkappa*efvisc**2*(yplus+dplus)**2)
        else
          pimp = 0.d0
        endif

        ! Neumann Boundary Condition
        !---------------------------

        if (iclptr.eq.1) then !TODO not available for k-eps

          qimp = -pimp*hint !TODO transform it, it is only to be fully equivalent

           call set_neumann_scalar &
           !==================
         ( coefa_ep(ifac), coefaf_ep(ifac),             &
           coefb_ep(ifac), coefbf_ep(ifac),             &
           qimp          , hint )

        ! Dirichlet Boundary Condition
        !-----------------------------

        else

          pimp = pimp + cvar_ep(iel)

          call set_dirichlet_scalar &
               !====================
             ( coefa_ep(ifac), coefaf_ep(ifac),             &
               coefb_ep(ifac), coefbf_ep(ifac),             &
               pimp          , hint           , rinfin )

        endif

      elseif (iturb.eq.32) then
        ! Use k at I'
        xkip = 0.5d0*(rijipb(ifac,1)+rijipb(ifac,2)+rijipb(ifac,3))

        ! Dirichlet Boundary Condition
        !-----------------------------

        pimp = 2.d0*visclc*xkip/(distbf**2*romc)

        call set_dirichlet_scalar &
             !====================
           ( coefa_ep(ifac), coefaf_ep(ifac),             &
             coefb_ep(ifac), coefbf_ep(ifac),             &
             pimp          , hint           , rinfin )

        ! ---> Alpha

        ! Dirichlet Boundary Condition
        !-----------------------------

        pimp = 0.d0

        hint = 1.d0/distbf

        call set_dirichlet_scalar &
             !====================
           ( coefa_al(ifac), coefaf_al(ifac),             &
             coefb_al(ifac), coefbf_al(ifac),             &
             pimp          , hint           , rinfin )

      endif

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

      ! Dirichlet Boundary Condition on k
      !----------------------------------

      ! Si on est hors de la sous-couche visqueuse (reellement ou via les
      ! scalable wall functions)
      if (iuntur.eq.1) then
        pimp = uk**2/sqrcmu

      ! Si on est en sous-couche visqueuse
      else
        pimp = 0.d0
      endif

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

      ! In viscous sub layer
      pimp_lam  = 120.d0*8.d0*visclc/(romc*ckwbt1*distbf**2)

      ! If we are outside the viscous sub-layer (either naturally, or
      ! artificialy using scalable wall functions)

      if (yplus > epzero) then
        pimp_turb = distbf*4.d0*uk**3*romc**2/           &
                   (sqrcmu*xkappa*visclc**2*(yplus+dplus)**2)

        ! Use gamma function of Kader to weight
        !between high and low Reynolds meshes

        gammap    = -0.01d0*(yplus+dplus)**4.d0/(1.d0+5.d0*(yplus+dplus))

        pimp      = pimp_lam*exp(gammap) + exp(1.d0/gammap)*pimp_turb
      else
        pimp      = pimp_lam
      endif

      qimp      = -pimp*hint !TODO transform it, it is only to be fully equivalent

      call set_neumann_scalar &
           !==================
         ( coefa_omg(ifac), coefaf_omg(ifac),             &
           coefb_omg(ifac), coefbf_omg(ifac),             &
           qimp           , hint )

    !===========================================================================
    ! 7.1 Boundary conditions on the Spalart Allmaras turbulence model
    !===========================================================================

    elseif (iturb.eq.70) then

      ! Dirichlet Boundary Condition on nusa
      !-------------------------------------

      pimp = 0.d0

      hint = visclc / distbf / csasig ! Note: nusa is zero at the wall

      call set_dirichlet_scalar &
           !====================
         ( coefa_nu(ifac), coefaf_nu(ifac),             &
           coefb_nu(ifac), coefbf_nu(ifac),             &
           pimp          , hint           , rinfin )

    endif

    byplus(ifac) = yplus
    bdplus(ifac) = dplus
    buk(ifac) = uk

  endif
  ! Test on the presence of a smooth wall (End)

enddo
! --- End of loop over faces

!===========================================================================
! 8. Boundary conditions on the other scalars
!    (Specific treatment for the variances of the scalars next to walls:
!     see condli)
!===========================================================================

do iscal = 1, nscal

  if (iscavr(iscal).le.0) then

    f_id = ivarfl(isca(iscal))

    call field_get_dim(f_id, f_dim)

    if (f_dim.eq.1) then

      call clptur_scalar &
   ( iscal  , isvhb  , icodcl ,                                     &
     rcodcl ,                                                       &
     byplus , bdplus , buk    ,                                     &
     hbord  , theipb ,                                              &
     tetmax , tetmin , tplumx , tplumn )


    ! Vector field
    else

      call clptur_vector &
   ( iscal  , isvhb  , icodcl ,                                     &
     rcodcl ,                                                       &
     byplus , bdplus , buk    ,                                     &
     hbord  )

    endif

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
  call parcpt (nlogla)
  call parcpt (nsubla)
  call parcpt (iuiptn)
  if (iscalt.gt.0) then
    call parmin (tetmin)
    call parmax (tetmax)
    call parmin (tplumn)
    call parmax (tplumx)
  endif
endif

!===============================================================================
! 9. Writings
!===============================================================================

!     Remarque : afin de ne pas surcharger les listings dans le cas ou
!       quelques yplus ne sont pas corrects, on ne produit le message
!       qu'aux deux premiers pas de temps ou le message apparait et
!       aux deux derniers pas de temps du calcul, ou si IWARNI est >= 2
!       On indique aussi le numero du dernier pas de temps auquel on
!       a rencontre des yplus hors bornes admissibles

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)

if (vcopt_u%iwarni.ge.0) then
  if (ntlist.gt.0) then
    modntl = mod(ntcabs,ntlist)
  elseif (ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
    modntl = 0
  else
    modntl = 1
  endif

  if ( (iturb.eq.0.and.nlogla.ne.0)                      .or.     &
       (itytur.eq.5.and.nlogla.ne.0)                     .or.     &
       ((itytur.eq.2.or.itytur.eq.3) .and. nsubla.gt.0)      )    &
       ntlast = ntcabs

  if ( (ntlast.eq.ntcabs.and.iaff.lt.2         ) .or.             &
       (ntlast.ge.0     .and.ntcabs.ge.ntmabs-1) .or.             &
       (ntlast.eq.ntcabs.and.vcopt_u%iwarni.ge.2) ) then
    iaff = iaff + 1

    if (iscalt.gt.0) then
      write(nfecra,2011) &
           uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,      &
           tetmin, tetmax, tplumn, tplumx, iuiptn,nsubla,nsubla+nlogla
    else
      write(nfecra,2010) &
           uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,   &
           iuiptn,nsubla,nsubla+nlogla
    endif

    if (iturb.eq. 0) write(nfecra,2020)  ntlast,ypluli
    if (itytur.eq.5) write(nfecra,2030)  ntlast,ypluli
    ! No warnings in EBRSM
    if (itytur.eq.2.or.iturb.eq.30.or.iturb.eq.31)                &
      write(nfecra,2040)  ntlast,ypluli
    if (vcopt_u%iwarni.lt.2.and.iturb.ne.32) then
      write(nfecra,2050)
    elseif (iturb.ne.32) then
      write(nfecra,2060)
    endif

  else if (modntl.eq.0 .or. vcopt_u%iwarni.ge.2) then

    if (iscalt.gt.0) then
      write(nfecra,2011) &
           uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,      &
           tetmin, tetmax, tplumn, tplumx, iuiptn,nsubla,nsubla+nlogla
    else
      write(nfecra,2010) &
           uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,   &
           iuiptn,nsubla,nsubla+nlogla
    endif
  endif

endif

!===============================================================================
! 10. Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(/,' LA NORMALE A LA FACE DE BORD DE PAROI ',I10,/,   &
         ' EST NULLE ; COORDONNEES : ',3E12.5)

 2010 format(/,                                                   &
 3X,'** CONDITIONS AUX LIMITES EN PAROI LISSE',/,                 &
 '   ----------------------------------------',/,                 &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Vitesse rel. en paroi    uiptn : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uet   : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uk    : ',2E12.5                 ,/,&
 '   Distance adimensionnelle yplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------   ',/,&
 '   Nbre de retournements de la vitesse en paroi : ',I10      ,/,&
 '   Nbre de faces en sous couche visqueuse       : ',I10      ,/,&
 '   Nbre de faces de paroi total                 : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)

 2011 format(/,                                                   &
 3X,'** CONDITIONS AUX LIMITES EN PAROI LISSE',/,                 &
 '   ----------------------------------------',/,                 &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Vitesse rel. en paroi    uiptn : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uet   : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uk    : ',2E12.5                 ,/,&
 '   Distance adimensionnelle yplus : ',2E12.5                 ,/,&
 '   Sca. thermal de frott.   tstar : ',2E12.5                 ,/,&
 '   Sca. thermal adim. rug.  tplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------'   ,/,&
 '   Nbre de retournements de la vitesse en paroi : ',I10      ,/,&
 '   Nbre de faces en sous couche visqueuse       : ',I10      ,/,&
 '   Nbre de faces de paroi total                 : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)

 2020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAFFINEMENT INSUFFISANT DU MAILLAGE EN PAROI',/,&
'@    =========                                               ',/,&
'@    Le maillage semble insuffisamment raffine en paroi      ',/,&
'@      pour pouvoir realiser un calcul laminaire.            ',/,&
'@                                                            ',/,&
'@    Le dernier pas de temps auquel ont ete observees de trop',/,&
'@      grandes valeurs de la distance adimensionnelle a la   ',/,&
'@      paroi (yplus) est le pas de temps ',I10                ,/,&
'@                                                            ',/,&
'@    La valeur minimale de yplus doit etre inferieure a la   ',/,&
'@      valeur limite YPLULI = ',E14.5                         ,/,&
'@                                                            ',/,&
'@    Observer la repartition de yplus en paroi (sous EnSight ',/,&
'@      ou ParaView par exemple) pour determiner dans quelle  ',/,&
'@      mesure la qualite des resultats est susceptible d etre',/,&
'@      affectee.')

 2030 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAFFINEMENT INSUFFISANT DU MAILLAGE EN PAROI',/,&
'@    =========                                               ',/,&
'@    Le maillage semble insuffisamment raffine en paroi      ',/,&
'@      pour pouvoir realiser un calcul type v2f              ',/,&
'@            (phi-fbar ou BL-v2/k)                           ',/,&
'@                                                            ',/,&
'@    Le dernier pas de temps auquel ont ete observees de trop',/,&
'@      grandes valeurs de la distance adimensionnelle a la   ',/,&
'@      paroi (yplus) est le pas de temps ',I10                ,/,&
'@                                                            ',/,&
'@    La valeur minimale de yplus doit etre inferieure a la   ',/,&
'@      valeur limite YPLULI = ',E14.5                         ,/,&
'@                                                            ',/,&
'@    Observer la repartition de yplus en paroi (sous EnSight ',/,&
'@      ou ParaView par exemple) pour determiner dans quelle  ',/,&
'@      mesure la qualite des resultats est susceptible d etre',/,&
'@      affectee.')

 2040 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MAILLAGE TROP FIN EN PAROI                  ',/,&
'@    =========                                               ',/,&
'@    Le maillage semble trop raffine en paroi pour utiliser  ',/,&
'@      un modele de turbulence haut Reynolds.                ',/,&
'@                                                            ',/,&
'@    Le dernier pas de temps auquel ont ete observees des    ',/,&
'@      valeurs trop faibles de la distance adimensionnelle a ',/,&
'@      la paroi (yplus) est le pas de temps ',I10             ,/,&
'@                                                            ',/,&
'@    La valeur minimale de yplus doit etre superieure a la   ',/,&
'@      valeur limite YPLULI = ',E14.5                         ,/,&
'@                                                            ',/,&
'@    Observer la repartition de yplus en paroi (sous EnSight ',/,&
'@      ou ParaView par exemple) pour determiner dans quelle  ',/,&
'@      mesure la qualite des resultats est susceptible d etre',/,&
'@      affectee.')
 2050 format(                                                     &
'@                                                            ',/,&
'@    Ce message ne s''affiche qu''aux deux premieres         ',/,&
'@      occurences du probleme et aux deux derniers pas de    ',/,&
'@      temps du calcul. La disparition du message ne signifie',/,&
'@      pas forcement la disparition du probleme.             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2060 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(/,' THE NORMAL TO THE WALL BOUNDARY FACE ',I10,/,    &
         ' IS NULL; COORDINATES: ',3E12.5)

 2010 format(/,                                                   &
 3X,'** BOUNDARY CONDITIONS FOR SMOOTH WALLS',/,                  &
 '   ---------------------------------------',/,                  &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Rel velocity at the wall uiptn : ',2E12.5                 ,/,&
 '   Friction velocity        uet   : ',2E12.5                 ,/,&
 '   Friction velocity        uk    : ',2E12.5                 ,/,&
 '   Dimensionless distance   yplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------'   ,/,&
 '   Nb of reversal of the velocity at the wall   : ',I10      ,/,&
 '   Nb of faces within the viscous sub-layer     : ',I10      ,/,&
 '   Total number of wall faces                   : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)

 2011 format(/,                                                   &
 3X,'** BOUNDARY CONDITIONS FOR SMOOTH WALLS',/,                  &
 '   ---------------------------------------',/,                  &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Rel velocity at the wall uiptn : ',2E12.5                 ,/,&
 '   Friction velocity        uet   : ',2E12.5                 ,/,&
 '   Friction velocity        uk    : ',2E12.5                 ,/,&
 '   Dimensionless distance   yplus : ',2E12.5                 ,/,&
 '   Friction thermal sca.    tstar : ',2E12.5                 ,/,&
 '   Rough dim-less th. sca.  tplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------'   ,/,&
 '   Nb of reversal of the velocity at the wall   : ',I10      ,/,&
 '   Nb of faces within the viscous sub-layer     : ',I10      ,/,&
 '   Total number of wall faces                   : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)


 2020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: MESH NOT ENOUGH REFINED AT THE WALL            ',/,&
'@    ========                                                ',/,&
'@    The mesh does not seem to be enough refined at the wall ',/,&
'@      to be able to run a laminar simulation.               ',/,&
'@                                                            ',/,&
'@    The last time step at which too large values for the    ',/,&
'@      dimensionless distance to the wall (yplus) have been  ',/,&
'@      observed is the time step ',I10                        ,/,&
'@                                                            ',/,&
'@    The minimum value for yplus must be lower than the      ',/,&
'@      limit value YPLULI = ',E14.5                           ,/,&
'@                                                            ',/,&
'@    Have a look at the distribution of yplus at the wall    ',/,&
'@      (with EnSight or ParaView for example) to conclude on ',/,&
'@      the way the results quality might be affected.        ')

 2030 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: MESH NOT ENOUGH REFINED AT THE WALL            ',/,&
'@    ========                                                ',/,&
'@    The mesh does not seem to be enough refined at the wall ',/,&
'@      to be able to run a v2f simulation                    ',/,&
'@      (phi-fbar or BL-v2/k)                                 ',/,&
'@                                                            ',/,&
'@    The last time step at which too large values for the    ',/,&
'@      dimensionless distance to the wall (yplus) have been  ',/,&
'@      observed is the time step ',I10                        ,/,&
'@                                                            ',/,&
'@    The minimum value for yplus must be lower than the      ',/,&
'@      limit value YPLULI = ',E14.5                           ,/,&
'@                                                            ',/,&
'@    Have a look at the distribution of yplus at the wall    ',/,&
'@      (with EnSight or ParaView for example) to conclude on ',/,&
'@      the way the results quality might be affected.        ')

 2040 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: MESH TOO REFINED AT THE WALL                   ',/,&
'@    ========                                                ',/,&
'@    The mesh seems to be too refined at the wall to use     ',/,&
'@      a high-Reynolds turbulence model.                     ',/,&
'@                                                            ',/,&
'@    The last time step at which too small values for the    ',/,&
'@      dimensionless distance to the wall (yplus) have been  ',/,&
'@      observed is the time step ',I10                        ,/,&
'@                                                            ',/,&
'@    The minimum value for yplus must be greater than the    ',/,&
'@      limit value YPLULI = ',E14.5                           ,/,&
'@                                                            ',/,&
'@    Have a look at the distribution of yplus at the wall    ',/,&
'@      (with EnSight or ParaView for example) to conclude on ',/,&
'@      the way the results quality might be affected.        ')
 2050 format(                                                     &
'@                                                            ',/,&
'@    This warning is only printed at the first two           ',/,&
'@      occurences of the problem and at the last time step   ',/,&
'@      of the calculation. The vanishing of the message does ',/,&
'@      not necessarily mean the vanishing of the problem.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

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
!> \param[in]     bdplus        dimensionless shift to the wall
!>                               for scalable wall functions
!> \param[in]     buk           dimensionless velocity
!> \param[in,out] hbord         exchange coefficient at boundary
!> \param[in]     theipb        value of thermal scalar at \f$ \centip \f$
!>                               of boundary cells
!> \param[out]    tetmax        maximum local ustar value
!> \param[out]    tetmin        minimum local ustar value
!> \param[out]    tplumx        maximum local tplus value
!> \param[out]    tplumn        minimum local tplus value
!_______________________________________________________________________________

subroutine clptur_scalar &
 ( iscal  , isvhb  , icodcl ,                                     &
   rcodcl ,                                                       &
   byplus , bdplus , buk    ,                                     &
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
use dimens, only: nvar
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

integer          icodcl(nfabor,nvar)

double precision rcodcl(nfabor,nvar,3)
double precision byplus(nfabor), bdplus(nfabor)
double precision hbord(nfabor), theipb(nfabor), buk(nfabor)
double precision tetmax, tetmin, tplumx, tplumn

! Local variables

integer          ivar, f_id, b_f_id, isvhbl
integer          ifac, iel, isou, jsou
integer          ifcvsl, itplus, itstar

double precision cpp, rkl, prdtl, visclc, romc, tplus, cofimp, cpscv
double precision distfi, distbf, fikis, hint_al, heq, hflui, hext
double precision yplus, dplus, phit, pimp, pimp_al, rcprod, temp, tet, uk
double precision viscis, visctc, xmutlm, ypth, xnuii
double precision rinfiv(3), pimpv(3)
double precision visci(3,3), hintt(6)
double precision turb_schmidt, exchange_coef

character(len=80) :: fname

double precision, dimension(:), pointer :: val_s, bval_s, crom, viscls
double precision, dimension(:), pointer :: viscl, visct, cpro_cp, cpro_cv

double precision, dimension(:), pointer :: bfconv, bhconv
double precision, dimension(:), pointer :: tplusp, tstarp, dist_theipb
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp, hextp
double precision, dimension(:), pointer :: a_al, b_al, af_al, bf_al
double precision, dimension(:,:), pointer :: coefaut, cofafut, cofarut, visten
double precision, dimension(:,:,:), pointer :: coefbut, cofbfut, cofbrut

double precision, allocatable, dimension(:) :: hbnd, hint, yptp

integer, save :: kbfid = -1

type(var_cal_opt) :: vcopt

logical(c_bool), dimension(:), pointer ::  cpl_faces

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

if (    iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0             &
    .or.ityturt(iscal).eq.3) then
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

ypth = 0.d0

if (ityturt(iscal).eq.3) then

  ! Turbulent diffusive flux of the scalar T
  ! (blending factor so that the component v'T' have only
  !  mu_T/(mu+mu_T)* Phi_T)

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

! EB-GGDH/AFM/DFM alpha boundary conditions
if (iturt(iscal).eq.11 .or. iturt(iscal).eq.21 .or. iturt(iscal).eq.31) then

  ! Name of the scalar ivar
  call field_get_name(ivarfl(ivar), fname)

  ! Index of the corresponding turbulent flux
  call field_get_id(trim(fname)//'_alpha', f_id)

  call field_get_coefa_s (f_id, a_al)
  call field_get_coefb_s (f_id, b_al)
  call field_get_coefaf_s(f_id, af_al)
  call field_get_coefbf_s(f_id, bf_al)
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
  if (vcopt%icoupl.gt.0) then
    allocate(dist_theipb(nfabor))
    call cs_ic_field_dist_data_by_face_id(f_id, 1, theipb, dist_theipb)
  endif
endif

if (vcopt%icoupl.gt.0) then
  call field_get_coupled_faces(f_id, cpl_faces)
endif

allocate(hbnd(nfabor),hint(nfabor),yptp(nfabor))

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
  dplus = bdplus(ifac)
  uk = buk(ifac)

  ! Test on the presence of a smooth wall condition (start)
  if (icodcl(ifac,iu).eq.5) then

    iel = ifabor(ifac)

    ! Physical quantities

    visclc = viscl(iel)
    visctc = visct(iel)
    romc   = crom(iel)

    xnuii = visclc / romc

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
      if (iscal.eq.iscalt .and. itherm.eq.3) then
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
        hint(ifac) = (rkl+vcopt%idifft*cpscv*visctc/turb_schmidt)/distbf
      else
        hint(ifac) = (rkl+vcopt%idifft*cpp*visctc/turb_schmidt)/distbf
      endif

    ! Symmetric tensor diffusivity (GGDH or AFM)
    else if (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then
      ! En compressible, pour l'energie LAMBDA/CV+CP/CV*(MUT/SIGMAS)
      if (iscal.eq.iscalt .and. itherm.eq.3) then
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
        temp = vcopt%idifft*cpscv*ctheta(iscal)/csrij
      else
        temp = vcopt%idifft*cpp*ctheta(iscal)/csrij
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

      hint(ifac) = viscis/surfbn(ifac)/fikis
    endif

    ! wall function and Dirichlet or Neumann on the scalar
    if (iturb.ne.0.and.(icodcl(ifac,ivar).eq.5.or.icodcl(ifac,ivar).eq.3)) then

      call hturbp(iwalfs,prdtl,turb_schmidt,yplus,dplus,hflui,ypth)
      ! Compute (y+-d+)/T+ *PrT
      yptp(ifac) = hflui/prdtl
      ! Compute lambda/y * (y+-d+)/T+
      hflui = rkl/distbf *hflui

    else

      ! y+/T+ *PrT
      yptp(ifac) = 1.d0/prdtl
      hflui = hint(ifac)

    endif

    hbnd(ifac) = hflui

  endif ! smooth wall condition

enddo

! internal coupling
if (vcopt%icoupl.gt.0) then
  ! Update exchange coef. in coupling entity of current scalar
  call cs_ic_field_set_exchcoeff(f_id, hbnd)
  ! Get external exchange coef.
  call field_get_hext(f_id, hextp)
endif

! --- Loop on boundary faces
do ifac = 1, nfabor

  yplus = byplus(ifac)
  dplus = bdplus(ifac)

  ! Test on the presence of a smooth wall condition (start)
  if (icodcl(ifac,iu).eq.5) then

    iel = ifabor(ifac)

    ! Physical quantities

    visclc = viscl(iel)
    visctc = visct(iel)
    romc   = crom(iel)

    xnuii = visclc / romc

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
    else
      rkl = viscls(iel)
    endif

    hext = rcodcl(ifac,ivar,2)
    pimp = rcodcl(ifac,ivar,1)

    if (vcopt%icoupl.gt.0) then
      if (cpl_faces(ifac)) then
        hext = hextp(ifac)/surfbn(ifac)
      endif
    endif

    hflui = hbnd(ifac)

    if (abs(hext).gt.rinfin*0.5d0) then
      heq = hflui
    else
      heq = hflui*hext/(hflui+hext)
    endif

    ! ---> Dirichlet Boundary condition with a wall function correction
    !      with or without an additional exchange coefficient hext

    if (icodcl(ifac,ivar).eq.5) then
      ! DFM: the gradient BCs are so that the production term
      !      of u'T' is correcty computed
      if (ityturt(iscal).ge.1) then
        ! In the log layer
        if (yplus.ge.ypth.and.iturb.ne.0) then
          xmutlm = xkappa*visclc*yplus
          rcprod = min(xkappa , max(1.0d0,sqrt(xmutlm/visctc))/yplus)

          cofimp = 1.d0 - yptp(ifac)*turb_schmidt/xkappa*                   &
                          (2.0d0*rcprod - 1.0d0/(2.0d0*yplus-dplus))

          ! In the viscous sub-layer
        else
          cofimp = 0.d0
        endif
      else
        cofimp = 1.d0 - heq/hint(ifac)
      endif

      ! To be coherent with a wall function, clip it to 0
      cofimp = max(cofimp, 0.d0)

      ! Gradient BCs
      coefap(ifac) = (1.d0 - cofimp)*pimp
      coefbp(ifac) = cofimp

      ! Flux BCs
      cofafp(ifac) = -heq*pimp
      cofbfp(ifac) =  heq

      ! Storage of the thermal exchange coefficient
      ! (conversion in case of energy or enthaly)
      ! the exchange coefficient is in W/(m2 K)
      ! Useful for thermal coupling or radiative transfer
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

    endif ! End if icodcl.eq.5

    !--> Turbulent heat flux
    if (ityturt(iscal).eq.3) then

      ! Turbulent diffusive flux of the scalar T
      ! (blending factor so that the component v'T' have only
      !  mu_T/(mu+mu_T)* Phi_T)
      if (icodcl(ifac,ivar).eq.5) then
        phit = cofafp(ifac)+cofbfp(ifac)*val_s(iel)
      elseif (icodcl(ifac,ivar).eq.3) then
        phit = rcodcl(ifac,ivar,3)
      elseif (icodcl(ifac,ivar).eq.1) then
        phit = heq *(val_s(iel) - pimp)
      else
        phit = 0.d0
      endif

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
      if (yplus.ge.ypth) then
        do isou = 1, 3
          pimpv(isou) = surfbo(isou,ifac)*phit/(surfbn(ifac)*cpp*romc)
        enddo
      else
        do isou = 1, 3
          pimpv(isou) = 0.d0
        enddo
      endif

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

    ! EB-GGDH/AFM/DFM alpha boundary conditions
    if (iturt(iscal).eq.11 .or. iturt(iscal).eq.21 .or. iturt(iscal).eq.31) then

      ! Dirichlet Boundary Condition
      !-----------------------------
      pimp_al = 0.d0
      hint_al = 1.d0/distbf

      call set_dirichlet_scalar &
           !====================
         ( a_al(ifac), af_al(ifac),             &
           b_al(ifac), bf_al(ifac),             &
           pimp_al   , hint_al    , rinfin )

    endif

    ! Save the values of T^star and T^+ for post-processing

    if (b_f_id.ge.0 .or. iscal.eq.iscalt) then

      ! Wall function
      if (icodcl(ifac,ivar).eq.5) then
        if (iscal.eq.iscalt) then
          phit = cofafp(ifac)+cofbfp(ifac)*theipb(ifac)
        else
          phit = cofafp(ifac)+cofbfp(ifac)*bval_s(ifac)
        endif
      ! Imposed flux with wall function for post-processing
      elseif (icodcl(ifac,ivar).eq.3) then
        phit = rcodcl(ifac,ivar,3) ! = 0 if current face is coupled
      elseif (icodcl(ifac,ivar).eq.1) then
        if (iscal.eq.iscalt) then
          phit = heq *(theipb(ifac) - pimp)
        else
          phit = heq *(bval_s(ifac) - pimp)
        endif
      else
        phit = 0.d0
      endif

      ! if face is coupled
      if (vcopt%icoupl.gt.0) then
        if (cpl_faces(ifac)) then
          phit = heq*(theipb(ifac)-dist_theipb(ifac))
        endif
      endif

      tet = phit/(romc*cpp*max(yplus-dplus, epzero)*xnuii/distbf)
      ! T+ = (T_I - T_w) / Tet
      tplus = max((yplus-dplus), epzero)/yptp(ifac)

      if (b_f_id.ge.0) bval_s(ifac) = bval_s(ifac) - tplus*tet

      if (itplus.ge.0) tplusp(ifac) = tplus
      if (itstar.ge.0) tstarp(ifac) = tet

      if (iscal.eq.iscalt) then
        tetmax = max(tet, tetmax)
        tetmin = min(tet, tetmin)
        tplumx = max(tplus,tplumx)
        tplumn = min(tplus,tplumn)
      endif

    endif

  endif ! smooth wall condition

enddo

deallocate(hbnd, hint, yptp)

if (iscal.eq.iscalt.and.vcopt%icoupl.gt.0) then
  deallocate(dist_theipb)
endif

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
!> \param[in]     bdplus        dimensionless shift to the wall
!>                               for scalable wall functions
!> \param[in]     buk           dimensionless velocity
!> \param[in,out] hbord         exchange coefficient at boundary
!_______________________________________________________________________________

subroutine clptur_vector &
 ( iscal  , isvhb  , icodcl ,                                     &
   rcodcl ,                                                       &
   byplus , bdplus , buk    ,                                     &
   hbord  )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only: nvar
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
use field_operator
use lagran
use turbomachinery
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal, isvhb

integer          icodcl(nfabor,nvar)

double precision rcodcl(nfabor,nvar,3)
double precision byplus(nfabor), bdplus(nfabor), buk(nfabor)
double precision hbord(nfabor)

! Local variables

integer          ivar, f_id, isvhbl
integer          ifac, iel
integer          ifcvsl

double precision cpp, rkl, prdtl, visclc, romc, cofimp
double precision distbf, heq, yptp, hflui, hext
double precision yplus, dplus, rcprod, uk
double precision visctc, xmutlm, ypth, xnuii, srfbnf
double precision rcodcx, rcodcy, rcodcz, rcodcn, rnx, rny, rnz
double precision turb_schmidt

double precision, dimension(:), pointer :: crom, viscls
double precision, dimension(:,:), pointer :: val_p_v
double precision, dimension(:), pointer :: viscl, visct, cpro_cp, hextp

double precision, dimension(:,:), pointer :: coefav, cofafv
double precision, dimension(:,:,:), pointer :: coefbv, cofbfv

double precision, allocatable, dimension(:) :: hbnd, hint

type(var_cal_opt) :: vcopt

logical(c_bool), dimension(:), pointer ::  cpl_faces

!===============================================================================

ivar = isca(iscal)
f_id = ivarfl(ivar)

call field_get_val_prev_v(f_id, val_p_v)

call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

call field_get_key_int (f_id, kivisl, ifcvsl)

if (ifcvsl .ge. 0) then
  call field_get_val_s(ifcvsl, viscls)
endif

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

call field_get_coefa_v(f_id, coefav)
call field_get_coefb_v(f_id, coefbv)
call field_get_coefaf_v(f_id, cofafv)
call field_get_coefbf_v(f_id, cofbfv)

call field_get_val_s(icrom, crom)
if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
endif

! retrieve turbulent Schmidt value for current vector
call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)

isvhbl = 0
if (iscal.eq.isvhb) then
  isvhbl = isvhb
endif

ypth = 0.d0

if (vcopt%icoupl.gt.0) then
  call field_get_coupled_faces(f_id, cpl_faces)
endif

allocate(hbnd(nfabor),hint(nfabor))

! --- Loop on boundary faces
do ifac = 1, nfabor

  yplus = byplus(ifac)
  dplus = bdplus(ifac)
  uk = buk(ifac)

  ! Geometric quantities
  distbf = distb(ifac)

  ! Test on the presence of a smooth wall condition (start)
  if (icodcl(ifac,iu).eq.5) then

    iel = ifabor(ifac)

    ! Physical quantities

    visclc = viscl(iel)
    visctc = visct(iel)
    romc   = crom(iel)

    xnuii = visclc / romc

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

    ! Scalar diffusivity
    if (iand(vcopt%idften, ISOTROPIC_DIFFUSION).ne.0) then
      hint(ifac) = (rkl+vcopt%idifft*cpp*visctc/turb_schmidt)/distbf
    else
      ! TODO if (vcopt%idften.eq.6)
      call csexit(1)
    endif

    ! wall function and Dirichlet or Neumann on the scalar
    if (iturb.ne.0.and.(icodcl(ifac,ivar).eq.5.or.icodcl(ifac,ivar).eq.3)) then

      call hturbp(iwalfs,prdtl,turb_schmidt,yplus,dplus,hflui,ypth)
      ! Compute (y+-d+)/T+ *PrT
      yptp = hflui/prdtl
      ! Compute lambda/y * (y+-d+)/T+
      hflui = rkl/distbf *hflui

    else

      ! y+/T+ *PrT
      yptp = 1.d0/prdtl
      hflui = hint(ifac)

    endif

    hbnd(ifac) = hflui

  endif ! smooth wall condition

enddo

! internal coupling
if (vcopt%icoupl.gt.0) then
  ! Update exchange coef. in coupling entity of current scalar
  call cs_ic_field_set_exchcoeff(f_id, hbnd)
  ! Get external exchange coef.
  call field_get_hext(f_id, hextp)
endif

! --- Loop on boundary faces
do ifac = 1, nfabor

  yplus = byplus(ifac)
  dplus = bdplus(ifac)

  ! Test on the presence of a smooth wall condition (start)
  if (icodcl(ifac,iu).eq.5) then

    srfbnf = surfbn(ifac)

    ! Local framework

    ! Unit normal

    rnx = surfbo(1,ifac)/srfbnf
    rny = surfbo(2,ifac)/srfbnf
    rnz = surfbo(3,ifac)/srfbnf

    ! Handle Dirichlet vector values

    rcodcx = rcodcl(ifac,ivar  ,1)
    rcodcy = rcodcl(ifac,ivar+1,1)
    rcodcz = rcodcl(ifac,ivar+2,1)

    !  Keep tangential part

    rcodcn = rcodcx*rnx+rcodcy*rny+rcodcz*rnz
    rcodcx = rcodcx -rcodcn*rnx
    rcodcy = rcodcy -rcodcn*rny
    rcodcz = rcodcz -rcodcn*rnz

    ! Physical quantities

    visclc = viscl(iel)
    visctc = visct(iel)

    hext = rcodcl(ifac,ivar,2)

    if (vcopt%icoupl.gt.0.and.cpl_faces(ifac)) then
      hext = hextp(ifac)/surfbn(ifac)
    endif

    hflui = hbnd(ifac)

    if (abs(hext).gt.rinfin*0.5d0) then
      heq = hflui
    else
      heq = hflui*hext/(hflui+hext)
    endif

    ! ---> Dirichlet Boundary condition with a wall function correction
    !      with or without an additional exchange coefficient hext

    if (icodcl(ifac,ivar).eq.5) then
      ! DFM: the gradient BCs are so that the production term
      !      of u'T' is correcty computed
      if (ityturt(iscal).ge.1) then
        ! In the log layer
        if (yplus.ge.ypth.and.iturb.ne.0) then
          xmutlm = xkappa*visclc*yplus
          rcprod = min(xkappa , max(1.0d0,sqrt(xmutlm/visctc))/yplus)

          cofimp = 1.d0 - yptp*turb_schmidt/xkappa*                        &
                  (2.0d0*rcprod - 1.0d0/(2.0d0*yplus-dplus))

          ! In the viscous sub-layer
        else
          cofimp = 0.d0
        endif
      else
        cofimp = 1.d0 - heq/hint(ifac)
      endif

      ! To be coherent with a wall function, clip it to 0
      cofimp = max(cofimp, 0.d0)

      ! Gradient boundary conditions
      !-----------------------------
      rcodcn = rcodcx*rnx+rcodcy*rny+rcodcz*rnz

      coefav(1,ifac) = (1.d0-cofimp)*(rcodcx - rcodcn*rnx) + rcodcn*rnx
      coefav(2,ifac) = (1.d0-cofimp)*(rcodcy - rcodcn*rny) + rcodcn*rny
      coefav(3,ifac) = (1.d0-cofimp)*(rcodcz - rcodcn*rnz) + rcodcn*rnz

      ! Projection in order to have the vector parallel to the wall
      ! B = cofimp * ( IDENTITY - n x n )

      coefbv(1,1,ifac) = cofimp*(1.d0-rnx**2)
      coefbv(2,2,ifac) = cofimp*(1.d0-rny**2)
      coefbv(3,3,ifac) = cofimp*(1.d0-rnz**2)
      coefbv(1,2,ifac) = -cofimp*rnx*rny
      coefbv(1,3,ifac) = -cofimp*rnx*rnz
      coefbv(2,1,ifac) = -cofimp*rny*rnx
      coefbv(2,3,ifac) = -cofimp*rny*rnz
      coefbv(3,1,ifac) = -cofimp*rnz*rnx
      coefbv(3,2,ifac) = -cofimp*rnz*rny

      ! Flux boundary conditions
      !-------------------------

      cofafv(1,ifac)   = -heq*(rcodcx - rcodcn*rnx) - hint(ifac)*rcodcn*rnx
      cofafv(2,ifac)   = -heq*(rcodcy - rcodcn*rny) - hint(ifac)*rcodcn*rny
      cofafv(3,ifac)   = -heq*(rcodcz - rcodcn*rnz) - hint(ifac)*rcodcn*rnz

      ! Projection
      !  B = heq*( IDENTITY - n x n )

      cofbfv(1,1,ifac) = heq*(1.d0-rnx**2) + hint(ifac)*rnx**2
      cofbfv(2,2,ifac) = heq*(1.d0-rny**2) + hint(ifac)*rny**2
      cofbfv(3,3,ifac) = heq*(1.d0-rnz**2) + hint(ifac)*rnz**2

      cofbfv(1,2,ifac) = (hint(ifac) - heq)*rnx*rny
      cofbfv(1,3,ifac) = (hint(ifac) - heq)*rnx*rnz
      cofbfv(2,1,ifac) = (hint(ifac) - heq)*rny*rnx
      cofbfv(2,3,ifac) = (hint(ifac) - heq)*rny*rnz
      cofbfv(3,1,ifac) = (hint(ifac) - heq)*rnz*rnx
      cofbfv(3,2,ifac) = (hint(ifac) - heq)*rnz*rny

    endif ! End if icodcl.eq.5

    ! TODO : postprocessing at the boundary

  endif ! smooth wall condition

enddo

deallocate(hbnd,hint)

!----
! End
!----

return
end subroutine
