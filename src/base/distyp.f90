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

!> \file distyp.f90
!>
!> \brief This subroutine computes the dimensionless distance to the wall
!> solving a steady transport equation.
!>
!> This function solves the following steady pure convection equation on
!> \f$ \varia \f$:
!> \f[
!> \divs \left( \varia \vect{V} \right)
!>     - \divs \left( \vect{V} \right) \varia = 0
!> \f]
!> where the vector field \f$ \vect{V} \f$ is defined by:
!> \f[
!>  \vect{V} = \dfrac{ \grad y }{\norm{\grad y} }
!> \f]
!> The boundary conditions on \f$ \varia \f$ read:
!> \f[
!>  \varia = \dfrac{u_\star}{\nu} \textrm{ on walls}
!> \f]
!> \f[
!>  \dfrac{\partial \varia}{\partial n} = 0 \textrm{ elsewhere}
!> \f]
!>
!> Then the dimensionless distance is deduced by:
!> \f[
!>  y^+ = y \varia
!> \f]
!>
!>
!> Then, Imposition of an amortization of Van Driest type for the LES.
!>        \f$ \nu_T \f$ is absorbed by \f$ (1-\exp(\dfrac{-y^+}{d^+}))^2 \f$
!>        where \f$ d^+ \f$ is set at 26.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     itypfb        boundary face types
!> \param[in]     visvdr        dynamic viscosity in edge cells after
!>                               driest velocity amortization
!_______________________________________________________________________________

subroutine distyp &
 ( itypfb , visvdr)

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
use parall
use period
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
double precision visvdr(ncelet)

! Local variables

integer          idtva0, f_id  , iconvp, idiffp
integer          f_id_yplus
integer          ndircp
integer          iphydp
integer          ifac  , iel   , init
integer          inc   , iccocg, isweep
integer          imucpp, idftnp
integer          nswrgp, nswrsp
integer          icvflb, iescap, imligp, ircflp, iswdyp, isstpp, ischcp, iwarnp
integer          iflmas, iflmab

integer          infpar
save             infpar

integer          ivoid(1)

double precision relaxp, blencp, climgp, epsilp, epsrgp, epsrsp, extrap
double precision thetap
double precision wall_surf
double precision xusnmx, xusnmn, xnorm0
double precision dismax, dismin, usna
double precision hint, pimp, qimp

double precision rvoid(1)

double precision, allocatable, dimension(:) :: dvarp, smbdp, rovsdp
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: viscap
double precision, dimension(:), pointer :: w_dist
double precision, pointer, dimension(:)   :: cvar_var
double precision, pointer, dimension(:)   :: cvara_var
double precision, dimension(:), pointer :: crom, uetbor
double precision, dimension(:), pointer :: viscl
double precision, dimension(:), pointer :: visct
double precision, pointer, dimension(:) :: coefap, coefbp
double precision, pointer, dimension(:) :: cofafp, cofbfp
double precision, pointer, dimension(:) :: a_y, b_y
double precision, pointer, dimension(:) :: af_y, bf_y
double precision, dimension(:), pointer :: imasfl, bmasfl
type(var_cal_opt) :: vcopt

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays for the distance resolution
allocate(dvarp(ncelet), smbdp(ncelet), rovsdp(ncelet))
allocate(dpvar(ncelet))
allocate(viscf(nfac), viscb(nfabor))
allocate(viscap(ncelet))

ipass  = ipass + 1

call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)

uetbor => null()

call field_get_id_try('ustar', f_id)
if (f_id.ge.0) then
  call field_get_val_s(f_id, uetbor)
endif

call field_get_id("wall_distance", f_id)
call field_get_val_s(f_id, w_dist)

call field_get_coefa_s( f_id, a_y)
call field_get_coefb_s( f_id, b_y)
call field_get_coefaf_s(f_id, af_y)
call field_get_coefbf_s(f_id, bf_y)


call field_get_id("wall_yplus", f_id_yplus)
call field_get_key_struct_var_cal_opt(f_id_yplus, vcopt)

call field_get_val_s(f_id_yplus, cvar_var)
call field_get_val_prev_s(f_id_yplus, cvara_var)

call field_get_coefa_s( f_id_yplus, coefap)
call field_get_coefb_s( f_id_yplus, coefbp)
call field_get_coefaf_s(f_id_yplus, cofafp)
call field_get_coefbf_s(f_id_yplus, cofbfp)

call field_get_key_int(f_id_yplus, kimasf, iflmas)
call field_get_key_int(f_id_yplus, kbmasf, iflmab)

! Get pointer to the convective mass flux
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

! Number of wall faces
if (ipass.eq.1) then
  infpar = 0
  do ifac = 1, nfabor
    if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
      infpar = infpar+1
    endif
  enddo
  if (irangp.ge.0) then
    call parcpt(infpar)
  endif
endif

! If no wall, no wall distance
if (infpar.eq.0) then
  do iel = 1, ncelet
    cvar_var(iel) = grand
  enddo

  return
endif

!===============================================================================
! 2. At the first time step
!===============================================================================

! Au premier pas de temps, on a en general u* = 0 (ou faux)
!   on ne calcule pas y+

! En effet ca prend du temps, d'autant plus que u* est petit, car il
!   alors calculer y+ jusqu'a une grande distance des parois

if (ntcabs.eq.1) then

  do iel = 1, ncel
    cvar_var(iel) = grand
  enddo

  if (vcopt%iwarni.ge.1) then
    write(nfecra,7000)
  endif

  return

endif

!===============================================================================
! 3. Boundary conditions
!===============================================================================

! Dirichlet u*/nu at walls, homogeneous Neumann elsewhere

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then
    iel = ifabor(ifac)

    ! Dirichlet Boundary Condition
    !-----------------------------

    hint = 1.d0/distb(ifac)
    pimp = uetbor(ifac)*crom(iel)/viscl(iel)

    call set_dirichlet_scalar &
         !====================
       ( coefap(ifac), cofafp(ifac),             &
         coefbp(ifac), cofbfp(ifac),             &
         pimp        , hint        , rinfin )

    ! Dirichlet Boundary Condition
    !-----------------------------

    hint = 1.d0/distb(ifac)
    pimp = 0.d0

    call set_dirichlet_scalar &
         !====================
       ( a_y(ifac), af_y(ifac),             &
         b_y(ifac), bf_y(ifac),             &
         pimp     , hint      , rinfin )


  else
    ! Neumann Boundary Conditions
    !----------------------------

    hint = 1.d0/distb(ifac)
    qimp = 0.d0

    call set_neumann_scalar &
         !==================
       ( coefap(ifac), cofafp(ifac),             &
         coefbp(ifac), cofbfp(ifac),             &
         qimp        , hint )

    ! Neumann Boundary Conditions
    !----------------------------

    hint = 1.d0/distb(ifac)
    qimp = 0.d0

    call set_neumann_scalar &
         !==================
       ( a_y(ifac), af_y(ifac),             &
         b_y(ifac), bf_y(ifac),             &
         qimp     , hint )

  endif
enddo

!===============================================================================
! 4. Compute the mass flux due to V = Grad(y)
!===============================================================================

call field_get_id("wall_distance", f_id)

! Default initilization at 0
init   = 1
! Take Dirichlet into account
inc    = 1
iccocg = 1
iphydp = 0

epsrgp = vcopt%epsrgr
extrap = vcopt%extrag
climgp = vcopt%climgr
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
iwarnp = vcopt%iwarni

! Pseudo viscosity, to compute the convective flux "1 grad(y). Sij"
do iel = 1, ncelet
  viscap(iel) = 1.d0
enddo

call viscfa &
  ( imvisf ,          &
  viscap ,            &
  viscf  , viscb  )

! If the equation on the wall distance has no flux-reconstruction (ircflu=0)
! then no reconstruction on the mass-flux (nswrgr)
if (vcopt%ircflu.eq.0) then
  nswrgp = 0
endif

! Compute convective mass flux
! here -div(1 grad(y))
call itrmas &
 ( f_id   , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydp ,     &
   0      , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   rvoid  ,                                                                    &
   w_dist ,                                                                    &
   a_y , b_y , af_y , bf_y ,                                                   &
   viscf  , viscb  ,                                                           &
   viscap ,                                                                    &
   imasfl , bmasfl )

! Now take the opposite
do ifac = 1, nfac
  imasfl(ifac) = - imasfl(ifac)
enddo
do ifac = 1, nfabor
  bmasfl(ifac) = - bmasfl(ifac)
enddo

!===============================================================================
! 5. Diagonal part of the matrix
!===============================================================================

do iel = 1, ncelet
  rovsdp(iel) = 0.d0
enddo

! Reinforce diagonal
do ifac = 1, nfac
  rovsdp(ifacel(1, ifac)) =  rovsdp(ifacel(1, ifac)) + imasfl(ifac)
  rovsdp(ifacel(2, ifac)) =  rovsdp(ifacel(2, ifac)) - imasfl(ifac)
enddo

do ifac = 1, nfabor
  rovsdp(ifabor(ifac)) =  rovsdp(ifabor(ifac)) + bmasfl(ifac)
enddo

do iel = 1, ncel
  rovsdp(iel) = 1.d-6 * abs(rovsdp(iel))
enddo

call synsca(rovsdp)

!===============================================================================
! 6. Time loop
!===============================================================================

! Initializations
!=================

! Inconnue
!   Au cas ou on n'atteint pas tout a fait l'etat stationnaire,
!   il faut que le yplus ne soit pas nul dans la zone ou les
!   conditions aux limites n'ont pas ete convectees. On voudrait
!   plutot que yplus y soit maximum.
!   Si on utilise zero ou une valeur negative comme initialisation,
!   on risque de se retrouver avec des valeurs proches de
!   zero issues de la diffusion due au schema upwind au voisinage
!   du front convecte et donc avec des yplus proches de zero
!   n'importe ou.
!   On va donc utiliser la valeur max de u*/nu.

!   A partir du second pas de temps, on a egalement le yplus du pas
!     de temps precedent

! On calcule le min et le max
xusnmx = -grand
xusnmn =  grand
do ifac = 1, nfabor
  if(itypfb(ifac).eq.iparoi .or.                            &
     itypfb(ifac).eq.iparug) then
    xusnmx = max(xusnmx,coefap(ifac))
    xusnmn = min(xusnmn,coefap(ifac))
  endif
enddo
if(irangp.ge.0) then
  call parmax (xusnmx)
  call parmin (xusnmn)
endif

if (ntcabs.eq.1) then
  do iel = 1, ncelet
    dvarp(iel) = xusnmx
  enddo
else
  do iel = 1, ncel
    usna = cvar_var(iel)/max(w_dist(iel),epzero)
    usna = max(usna,xusnmn)
    usna = min(usna,xusnmx)
    dvarp(iel) = usna
  enddo
endif

! L2 norm of (u*/nu) over wall boundary faces
xnorm0 = 0.d0
wall_surf = 0.d0
do ifac = 1, nfabor
  if(itypfb(ifac).eq.iparoi .or.                            &
    itypfb(ifac).eq.iparug) then
    wall_surf = wall_surf + surfbn(ifac)
    xnorm0 = xnorm0 + coefap(ifac)**2 * surfbn(ifac)
  endif
enddo
if(irangp.ge.0) then
  call parsom (wall_surf)
  call parsom (xnorm0)
endif
xnorm0 = sqrt(xnorm0 / wall_surf) * voltot

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(dvarp)
endif

! Right hand side
!=================
do iel = 1, ncel
  smbdp(iel) = 0.d0
enddo

! Solving
!=========
iconvp = vcopt%iconv
idiffp = vcopt%idiff
idftnp = vcopt%idften
nswrsp = vcopt%nswrsm
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
ircflp = vcopt%ircflu
ischcp = vcopt%ischcv
isstpp = vcopt%isstpc
! No error estimate
iescap = 0
imucpp = 0
iswdyp = vcopt%iswdyn
iwarnp = vcopt%iwarni
blencp = vcopt%blencv
epsilp = vcopt%epsilo
epsrsp = vcopt%epsrsm
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
extrap = vcopt%extrag
relaxp = vcopt%relaxv
thetap = vcopt%thetav
! all boundary convective flux with upwind
icvflb = 0
init   = 1

! There are som Dirichlet BCs
ndircp = 1
! No steady state algo
idtva0 = 0
! no over loops
isweep = -1

! Warning: no diffusion so no need of other diffusive Boundary coefficient

call codits &
 ( idtva0 , isweep , f_id_yplus, iconvp , idiffp , ndircp ,       &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp , xnorm0 ,                                              &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   dvarp  , dvarp  ,                                              &
   coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
   imasfl , bmasfl ,                                              &
   imasfl , bmasfl , imasfl , bmasfl , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdp , smbdp  , dvarp  , dpvar  ,                            &
   rvoid  , rvoid  )

! Clipping (indispensable si on initialise par u*/nu du pas de
!==========                                   temps precedent)

do iel = 1, ncel
  dvarp(iel) = max(dvarp(iel),xusnmn)
  dvarp(iel) = min(dvarp(iel),xusnmx)
enddo

!===============================================================================
! 7. Finalization and printing
!===============================================================================

do iel = 1, ncel
  cvar_var(iel) = dvarp(iel)*w_dist(iel)
enddo

dismax = -grand
dismin =  grand

do iel = 1, ncel
  dismin = min(cvar_var(iel),dismin)
  dismax = max(cvar_var(iel),dismax)
enddo

if (irangp.ge.0) then
  call parmin(dismin)
  call parmax(dismax)
endif

if (vcopt%iwarni.ge.1) then
  write(nfecra,1000) dismin, dismax
endif

!===============================================================================
! 8. Van Driest amortization
!===============================================================================

call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

do iel = 1, ncel
  visct(iel) = visct(iel)*(1.0d0-exp(-cvar_var(iel)/cdries))**2
enddo

! For the wall cells we add the turbulent viscosity which was absorbed
! in clptur and which has served to calculate the boundary conditions
do iel = 1, ncel
  if (visvdr(iel).gt.-900.d0) visct(iel) = visvdr(iel)
enddo

! Free memory
deallocate(dvarp, smbdp, rovsdp)
deallocate(dpvar)
deallocate(viscf, viscb)
deallocate(viscap)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format( &
''                                                             ,/,&
' ** DISTANCE A LA PAROI ADIMENSIONNELLE'                      ,/,&
'    ------------------------------------'                     ,/,&
''                                                             ,/,&
'  Distance+ min = ',E14.5    ,' Distance+ max = ',E14.5       ,/)
 7000 format( &
''                                                             ,/,&
' ** DISTANCE A LA PAROI ADIMENSIONNELLE'                      ,/,&
'    ------------------------------------'                     ,/,&
''                                                             ,/,&
'  Elle n''est pas calculee au premier pas de temps'           ,/)

#else

 1000 format( &
'                                                             ',/,&
' ** DIMENSIONLESS WALL DISTANCE                              ',/,&
'    ---------------------------                              ',/,&
'                                                             ',/,&
'  Min distance+ = ',E14.5    ,' Max distance+ = ',E14.5       ,/,&
'                                                             ',/,&
'     (Distance calculation done in ',I10   ,' iterations)'    ,/)
 7000 format( &
''                                                             ,/,&
' ** DIMENSIONLESS WALL DISTANCE'                              ,/,&
'    ---------------------------'                              ,/,&
''                                                             ,/,&
'  It is not computed at the first time step'                  ,/)

#endif

!----
! End
!----

return
end subroutine
