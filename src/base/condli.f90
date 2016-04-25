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
! Function :
! --------

!> \file condli.f90
!>
!> \brief Translation of the boundary conditions given by
!> cs_user_boundary_conditions
!> in a form that fits to the solver.
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
!> - if we consider an increment of a variable, the boundary conditions
!>   read:
!>   \f[
!>   \delta P_\centf = B_P^g \delta P_\centi
!>   \f]
!>   and
!>   \f[
!>   \delta Q_\centf = B_P^f \delta P_\centi
!>   \f]
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
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     isvhb         indicator to save exchange coeffient
!>                               at the walls
!> \param[in]     iterns        iteration number on Navier-Stokes equations
!> \param[in,out] icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[in,out] isostd        indicator for standard outlet
!>                               and reference face index
!> \param[in]     dt            time step (per cell)
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
!> \param[out]    visvdr        viscosite dynamique ds les cellules
!>                               de bord apres amortisst de v driest
!> \param[out]    hbord         coefficients d'echange aux bords
!> \param[out]    theipb        boundary temperature in \f$ \centip \f$
!>                               (more exaclty the energetic variable)
!> \param[in]     frcxt         external force responsible for the hydrostatic
!>                               pressure
!_______________________________________________________________________________

subroutine condli &
 ( nvar   , nscal  , iterns ,                                     &
   isvhb  ,                                                       &
   icodcl , isostd ,                                              &
   dt     , rcodcl ,                                              &
   visvdr , hbord  , theipb , frcxt  )

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
use field_operator
use turbomachinery
use darcy_module

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal , iterns
integer          isvhb

integer          icodcl(nfabor,nvarcl)
integer          isostd(nfabor+1)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvarcl,3)
double precision frcxt(3,ncelet)
double precision visvdr(ncelet)
double precision hbord(nfabor),theipb(nfabor)

! Local variables

integer          ifac  , iel   , ivar
integer          isou  , jsou  , ii
integer          ihcp  , iscal
integer          inc   , iprev , iccocg
integer          iok   , iok1
integer          icodcu
integer          isoent, isorti, ncpt,   isocpt(2)
integer          iclsym, ipatur, ipatrg, isvhbl
integer          ifcvsl, ipccv
integer          itplus, itstar
integer          f_id  ,  iut  , ivt   , iwt, iflmab

double precision sigma , cpp   , rkl
double precision hint  , hext  , pimp  , xdis, qimp, cfl
double precision hintt(6)
double precision flumbf, visclc, visctc, distbf, srfbn2
double precision xxp0, xyp0, xzp0
double precision srfbnf, normal(3)
double precision vistot
double precision rinfiv(3), pimpv(3), qimpv(3), hextv(3), cflv(3)
double precision visci(3,3), fikis, viscis, distfi
double precision temp

character(len=80) :: fname

double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:,:) :: velipb, rijipb
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:,:,:) :: gradv
double precision, pointer, dimension(:,:) :: dttens, visten
double precision, dimension(:), pointer :: tplusp, tstarp, yplbr
double precision, pointer, dimension(:,:) :: forbr
double precision, dimension(:,:), pointer :: coefaut, cofafut, cofarut
double precision, dimension(:,:,:), pointer :: coefbut, cofbfut, cofbrut
double precision, dimension(:), pointer :: bmasfl
double precision, dimension(:), pointer :: bfconv, bhconv
double precision, dimension(:,:), pointer :: coefau, cofafu, cfaale, claale
double precision, dimension(:,:,:), pointer :: coefbu, cofbfu, cfbale, clbale
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: cofadp, cofbdp
double precision, dimension(:,:), pointer :: vel, vela
double precision, dimension(:), pointer :: crom

double precision, dimension(:), pointer :: viscl, visct, viscls
double precision, dimension(:), pointer :: cpro_cp, cpro_cv, cvar_s, cvara_s
double precision, dimension(:), pointer :: cpro_visma1, cpro_visma2, cpro_visma3

! Darcy arrays
double precision, dimension(:), pointer :: permeability
double precision, dimension(:,:), pointer :: tensor_permeability
integer fid

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

if (ippmod(idarcy).eq.1) then
  if (darcy_anisotropic_permeability.eq.0) then
    call field_get_val_s_by_name('permeability', permeability)
  else
    call field_get_id('permeability', fid)
    call field_get_val_v(fid, tensor_permeability)
  endif
endif

! Allocate temporary arrays
allocate(velipb(nfabor,3))

! coefa and coefb are required to compute the cell gradients for the wall
!  turbulent boundary conditions.
! So, their initial values are kept (note that at the first time step, they are
!  initialized to zero flux in inivar.f90)

! velipb stores the velocity in I' of boundary cells

! Initialize variables to avoid compiler warnings

ipccv = 0

rinfiv(1) = rinfin
rinfiv(2) = rinfin
rinfiv(3) = rinfin

! pointers to y+, T+ and T* if saved

yplbr => null()
tplusp => null()
tstarp => null()

! Initialization of the array storing yplus
!  which is computed in clptur.f90 and/or clptrg.f90
call field_get_id_try('yplus', f_id)
if (f_id.ge.0) then
  call field_get_val_s(f_id, yplbr)
  do ifac = 1, nfabor
    yplbr(ifac) = 0.d0
  enddo
endif

call field_get_id_try('tplus', itplus)
if (itplus.ge.0) then
  call field_get_val_s (itplus, tplusp)
  do ifac = 1, nfabor
    tplusp(ifac) = 0.d0
  enddo
endif

call field_get_id_try('tstar', itstar)
if (itstar.ge.0) then
  call field_get_val_s (itstar, tstarp)
  do ifac = 1, nfabor
    tstarp(ifac) = 0.d0
  enddo
endif

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)
call field_get_val_prev_v(ivarfl(iu), vela)

! Pointers to the mass fluxes
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

! Pointers to specific fields
if (ifconv.ge.0) call field_get_val_s(ifconv, bfconv)
if (ihconv.ge.0) call field_get_val_s(ihconv, bhconv)

if (idtten.ge.0) call field_get_val_v(idtten, dttens)

if (ineedf.eq.1 .and. iterns.eq.1) call field_get_val_v(iforbr, forbr)

! Pointers to velcocity BC coefficients
call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)
call field_get_coefaf_v(ivarfl(iu), cofafu)
call field_get_coefbf_v(ivarfl(iu), cofbfu)

!===============================================================================
! 2. Treatment of types of BCs given by itypfb
!===============================================================================

if (ippmod(iphpar).ge.1.and.ippmod(igmix).eq.-1) then
  call pptycl &
  !==========
 ( nvar   ,                                                       &
   icodcl , itypfb , izfppp ,                                     &
   dt     ,                                                       &
   rcodcl )
endif

if (iale.eq.1) then
  call altycl &
  !==========
 ( itypfb , ialtyb , icodcl , impale ,                            &
   dt     ,                                                       &
   rcodcl , xyzno0 )
endif

if (imobil.eq.1 .or. iturbo.eq.1 .or. iturbo.eq.2) then
  call mmtycl(itypfb, rcodcl)
  !==========
endif

call typecl &
!==========
 ( nvar   , nscal  ,                                              &
   itypfb , itrifb , icodcl , isostd ,                            &
   rcodcl , frcxt  )

!===============================================================================
! 3. Check the consistency of the BCs
!===============================================================================

call vericl                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   itypfb , icodcl ,                                              &
   rcodcl )

!===============================================================================
! 4. Deprecated model to compute wall distance
!===============================================================================
! attention, si on s'en sert pour autre chose, disons le module X
!   bien faire attention dans verini avec quelles options le module
!   X doit alors etre incompatible (perio, parall).

iok1 = 0

if(ineedy.eq.1.and.abs(icdpar).eq.2) then

  ! Allocate a temporary array
  allocate(w1(ncelet))

  ! ON FERA ATTENTION EN PARALLELISME OU PERIODICITE
  !    (UNE PAROI PEUT ETRE PLUS PROCHE EN TRAVERSANT UN BORD ...)

  do iel = 1, ncel
    w1(iel) = grand
  enddo

  do ifac = 1, nfabor
    icodcu = icodcl(ifac,iu)
    if( icodcu.eq.5 .or. icodcu.eq.6 ) then
      do iel = 1, ncel
        xdis =                                                &
              (cdgfbo(1,ifac)-xyzcen(1,iel))**2                   &
             +(cdgfbo(2,ifac)-xyzcen(2,iel))**2                   &
             +(cdgfbo(3,ifac)-xyzcen(3,iel))**2
        if(w1(iel).gt.xdis) then
          w1(iel) = xdis
          ifapat(iel) = ifac
        endif
      enddo
    endif
  enddo

  ! Free memory
  deallocate(w1)

  iok = 0
  do iel = 1, ncel
    if(ifapat(iel).le.0)then
      iok = iok + 1
    endif
  enddo
  if(iok.gt.0) then
    write(nfecra,1000) irijec, idries
    iok1 = 1
  endif

endif

! Normalement, on ne passe pas en parallele ici,
!   mais au cas ou ...
if(irangp.ge.0) then
  call parcpt(iok1)
endif

if(iok1.ne.0) then
  call csexit (1)
  !==========
endif

!===============================================================================
! 5. Variables
!===============================================================================

! --- Variables
xxp0   = xyzp0(1)
xyp0   = xyzp0(2)
xzp0   = xyzp0(3)

! --- Physical quantities
call field_get_val_s(iprpfl(iviscl), viscl)
call field_get_val_s(iprpfl(ivisct), visct)

! --- Compressible
if (ippmod(icompf).ge.0) then
  if(icv.gt.0) then
    ipccv  = ipproc(icv)
  else
    ipccv = 0
  endif
endif

!===============================================================================
! 6. Compute the temperature or the enthalpy in I' for boundary cells
!     (thanks to the formula: Fi + GRAD(Fi).II')

!    For the coupling with SYRTHES
!     theipb is used by coupbo after condli
!    For the coupling with the 1D wall thermal module
!     theipb is used by cou1do after condli
!    For the radiation module
!     theipb is used to compute the required flux in raypar

!        CECI POURRAIT EN PRATIQUE ETRE HORS DE LA BOUCLE.

!===============================================================================

! Allocate a temporary array for the gradient reconstruction
allocate(grad(3,ncelet))

!  Pour le couplage SYRTHES ou module thermique 1D
!  -----------------------------------------------
!  Ici, on fait une boucle "inutile"  (on ne fait quelque chose
!    que pour icpsyr(iscal) = 1). C'est pour preparer le traitement
!    eventuel de plusieurs temperatures (ie plusieurs couplages
!    SYRTHES a la fois ; noter cependant que meme dans ce cas,
!    une seule temperature sera recue de chaque couplage. En polyph,
!    il faudrait ensuite reconstruire les enthalpies ...
!    plus tard si necessaire).
!  Ici, il ne peut y avoir qu'un seul scalaire avec icpsyr = 1 et
!    ce uniquement s'il y a effectivement couplage avec SYRTHES
!    (sinon, on s'est arrete dans verini)
!  Dans le cas du couplage avec le module 1D, on utilise iscalt.

!  Pour le rayonnement
!  -------------------
!  On calcule la valeur en I' s'il y a une variable
!    thermique


!  On recherche l'unique scalaire qui convient
!     (ce peut etre T, H, ou E (en compressible))

! Compute the boundary value of the thermal scalar in I' if required
if (iscalt.gt.0) then

  ivar = isca(iscalt)

  if (ntcabs.gt.1 .and. itbrrb.eq.1 .and. ircflu(ivar).eq.1) then

    call field_get_val_s(ivarfl(ivar), cvar_s)

    inc = 1
    iprev = 1
    iccocg = 1

    call field_gradient_scalar(ivarfl(ivar), iprev, imrgra, inc,  &
                               iccocg,                            &
                               grad)

    do ifac = 1 , nfabor
      iel = ifabor(ifac)
      theipb(ifac) = cvar_s(iel) &
                   + grad(1,iel)*diipb(1,ifac) &
                   + grad(2,iel)*diipb(2,ifac) &
                   + grad(3,iel)*diipb(3,ifac)
    enddo

  else

    call field_get_val_prev_s(ivarfl(ivar), cvara_s)

    do ifac = 1 , nfabor
      iel = ifabor(ifac)
      theipb(ifac) = cvara_s(iel)
    enddo

  endif

endif

!===============================================================================
! 6.bis Compute the velocity and Renolds stesses tensor in I' for boundary cells
!        (thanks to the formula: Fi + GRAD(Fi).II') if there are symmetry or
!         wall with wall functions boundary conditions
!===============================================================================

! ---> Indicator for symmetries or wall with wall functions

iclsym = 0
ipatur = 0
ipatrg = 0
do ifac = 1, nfabor
  if (icodcl(ifac,iu).eq.4) then
    iclsym = 1
  elseif (icodcl(ifac,iu).eq.5) then
    ipatur = 1
  elseif (icodcl(ifac,iu).eq.6) then
    ipatrg = 1
  endif
  if (iclsym.ne.0.and.ipatur.ne.0.and.ipatrg.ne.0) goto 100
enddo
100 continue

if (irangp.ge.0) then
  call parcmx(iclsym)
  call parcmx(ipatur)
  call parcmx(ipatrg)
endif

! ---> Compute the velocity in I' for boundary cells

if (iclsym.ne.0.or.ipatur.ne.0.or.ipatrg.ne.0.or.ineedf.eq.1) then

  if (ntcabs.gt.1) then

    ! Allocate a temporary array
    allocate(gradv(3,3,ncelet))

    inc = 1
    iprev = 1

    call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                               gradv)

    do isou = 1, 3
      do ifac = 1, nfabor
        iel = ifabor(ifac)
        velipb(ifac,isou) = gradv(1,isou,iel)*diipb(1,ifac)    &
                          + gradv(2,isou,iel)*diipb(2,ifac)    &
                          + gradv(3,isou,iel)*diipb(3,ifac)    &
                          + vela(isou,iel)
      enddo
    enddo

    deallocate(gradv)

  ! Nb: at the first time step, coefa and coefb are unknown, so the walue
  !     in I is stored instead of the value in I'
  else

    do isou = 1, 3

      do ifac = 1, nfabor
        iel = ifabor(ifac)
        velipb(ifac,isou) = vela(isou,iel)
      enddo

    enddo

  endif

endif

! ---> Compute Rij in I' for boundary cells

if ((iclsym.ne.0.or.ipatur.ne.0.or.ipatrg.ne.0).and.itytur.eq.3) then

  ! Allocate a work array to store Rij values at boundary faces
  allocate(rijipb(nfabor,6))

  do isou = 1 , 6

    if (isou.eq.1) ivar = ir11
    if (isou.eq.2) ivar = ir22
    if (isou.eq.3) ivar = ir33
    if (isou.eq.4) ivar = ir12
    if (isou.eq.5) ivar = ir13
    if (isou.eq.6) ivar = ir23


    if (ntcabs.gt.1.and.irijrb.eq.1) then

      call field_get_val_s(ivarfl(ivar), cvar_s)

      inc = 1
      iprev = 1
      iccocg = 1

      call field_gradient_scalar(ivarfl(ivar), iprev, imrgra, inc,  &
                                 iccocg,                            &
                                 grad)

      do ifac = 1 , nfabor
        iel = ifabor(ifac)
        rijipb(ifac,isou) = cvar_s(iel)               &
                          + grad(1,iel)*diipb(1,ifac) &
                          + grad(2,iel)*diipb(2,ifac) &
                          + grad(3,iel)*diipb(3,ifac)
      enddo

    ! Nb: at the first time step, coefa and coefb are unknown, so the walue
    !     in I is stored instead of the value in I'
    else

      call field_get_val_prev_s(ivarfl(ivar), cvara_s)

      do ifac = 1 , nfabor
        iel = ifabor(ifac)
        rijipb(ifac,isou) = cvara_s(iel)
      enddo

    endif

  enddo

endif

! Free memory
deallocate(grad)

!===============================================================================
! 7. Turbulence at walls:
!       (u,v,w,k,epsilon,Rij,temperature)
!===============================================================================
! --- On a besoin de velipb et de rijipb (et theipb pour le rayonnement)

!     On initialise visvdr a -999.d0.
!     Dans clptur, on amortit la viscosite turbulente sur les cellules
!     de paroi si on a active van Driest. La valeur finale est aussi
!     stockee dans visvdr.
!     Plus loin, dans vandri, la viscosite sur les cellules
!     de paroi sera amortie une seconde fois. On se sert alors de
!     visvdr pour lui redonner une valeur correcte.
if(itytur.eq.4.and.idries.eq.1) then
  do iel=1,ncel
    visvdr(iel) = -999.d0
  enddo
endif

if (ipatur.ne.0) then

  ! Smooth wall laws
  call clptur &
  !==========
 ( nscal  , isvhb  , icodcl ,                                     &
   rcodcl ,                                                       &
   velipb , rijipb , visvdr ,                                     &
   hbord  , theipb )

endif

if (ipatrg.ne.0) then

  ! Rough wall laws
  call clptrg &
  !==========
 ( nscal  , isvhb  , icodcl ,                                     &
   rcodcl ,                                                       &
   velipb , rijipb , visvdr ,                                     &
   hbord  , theipb )

endif

!===============================================================================
! 8. Symmetry for vectors and tensors
!       (u,v,w,Rij)
!===============================================================================
!   On a besoin de velipb et de rijipb

do ifac = 1, nfabor
  isympa(ifac) = 1
enddo

if (iclsym.ne.0) then

  call clsyvt &
  !==========
 ( nscal  , icodcl ,                                              &
   rcodcl ,                                                       &
   velipb , rijipb )

endif

!===============================================================================
! 9. Velocity: Outlet, Dirichlet and Neumann and convective outlet
!===============================================================================

! ---> Outlet: in case of incomming mass flux, the mass flux is set to zero.

isoent = 0
isorti = 0
do ifac = 1, nfabor

  if (icodcl(ifac,iu).eq.9) then

    flumbf = bmasfl(ifac)

    ! --- Physical Properties
    iel = ifabor(ifac)
    visclc = viscl(iel)
    visctc = visct(iel)

    ! --- Geometrical quantities
    distbf = distb(ifac)

    if (itytur.eq.3) then
      hint =  visclc          /distbf
    else
      hint = (visclc + visctc)/distbf
    endif

    isorti = isorti + 1

    if (flumbf.lt.-epzero) then

      ! Dirichlet Boundary Condition
      !-----------------------------

      ! Coupled solving of the velocity components

      pimpv(1) = 0.d0
      pimpv(2) = 0.d0
      pimpv(3) = 0.d0

      call set_dirichlet_vector &
           !====================
         ( coefau(:,ifac)  , cofafu(:,ifac)  ,             &
           coefbu(:,:,ifac), cofbfu(:,:,ifac),             &
           pimpv           , hint            , rinfiv )


      isoent = isoent + 1

    else

      ! Neumann Boundary Conditions
      !----------------------------

      qimp = 0.d0

      ! Coupled solving of the velocity components

      qimpv(1) = 0.d0
      qimpv(2) = 0.d0
      qimpv(3) = 0.d0

      call set_neumann_vector &
           !==================
         ( coefau(:,ifac)  , cofafu(:,ifac)  ,             &
           coefbu(:,:,ifac), cofbfu(:,:,ifac),             &
           qimpv           , hint )

    endif

  endif

enddo

if (mod(ntcabs,ntlist).eq.0 .or. iwarni(iu).ge. 0) then
  isocpt(1) = isoent
  isocpt(2) = isorti
  if (irangp.ge.0) then
    ncpt = 2
    call parism(ncpt, isocpt)
  endif
  if (isocpt(2).gt.0 .and. (iwarni(iu).ge.2.or.isocpt(1).gt.0)) then
    write(nfecra,3010) isocpt(1), isocpt(2)
  endif
endif

! ---> Dirichlet and Neumann

do ifac = 1, nfabor

  iel = ifabor(ifac)

  ! --- Physical Propreties
  visclc = viscl(iel)
  visctc = visct(iel)

  ! --- Geometrical quantities
  distbf = distb(ifac)

  if (itytur.eq.3) then
    hint =   visclc         /distbf
  else
    hint = ( visclc+visctc )/distbf
  endif

  ! Dirichlet Boundary Conditions
  !------------------------------

  if (icodcl(ifac,iu).eq.1) then


    pimpv(1) = rcodcl(ifac,iu,1)
    pimpv(2) = rcodcl(ifac,iv,1)
    pimpv(3) = rcodcl(ifac,iw,1)
    hextv(1) = rcodcl(ifac,iu,2)
    hextv(2) = rcodcl(ifac,iv,2)
    hextv(3) = rcodcl(ifac,iw,2)

    call set_dirichlet_vector &
         !====================
       ( coefau(:,ifac)  , cofafu(:,ifac)  ,             &
         coefbu(:,:,ifac), cofbfu(:,:,ifac),             &
         pimpv           , hint            , hextv )


  ! Neumann Boundary Conditions
  !----------------------------

  elseif (icodcl(ifac,iu).eq.3) then

    ! Coupled solving of the velocity components

    qimpv(1) = rcodcl(ifac,iu,3)
    qimpv(2) = rcodcl(ifac,iv,3)
    qimpv(3) = rcodcl(ifac,iw,3)

    call set_neumann_vector &
         !==================
       ( coefau(:,ifac)  , cofafu(:,ifac)  ,             &
         coefbu(:,:,ifac), cofbfu(:,:,ifac),             &
         qimpv           , hint )

  ! Convective Boundary Conditions
  !-------------------------------

  elseif (icodcl(ifac,iu).eq.2) then

    ! Coupled solving of the velocity components

    pimpv(1) = rcodcl(ifac,iu,1)
    cflv(1) = rcodcl(ifac,iu,2)
    pimpv(2) = rcodcl(ifac,iv,1)
    cflv(2) = rcodcl(ifac,iv,2)
    pimpv(3) = rcodcl(ifac,iw,1)
    cflv(3) = rcodcl(ifac,iw,2)

    call set_convective_outlet_vector &
         !==================
       ( coefau(:,ifac)  , cofafu(:,ifac)  ,             &
         coefbu(:,:,ifac), cofbfu(:,:,ifac),             &
         pimpv           , cflv            , hint )

  ! Imposed value for the convection operator, imposed flux for diffusion
  !----------------------------------------------------------------------

  elseif (icodcl(ifac,iu).eq.13) then

    pimpv(1) = rcodcl(ifac,iu,1)
    pimpv(2) = rcodcl(ifac,iv,1)
    pimpv(3) = rcodcl(ifac,iw,1)

    qimpv(1) = rcodcl(ifac,iu,3)
    qimpv(2) = rcodcl(ifac,iv,3)
    qimpv(3) = rcodcl(ifac,iw,3)

    call set_dirichlet_conv_neumann_diff_vector &
         !=====================================
       ( coefau(:,ifac)  , cofafu(:,ifac)  ,             &
         coefbu(:,:,ifac), cofbfu(:,:,ifac),             &
         pimpv           , qimpv           )

  ! Convective Boundary For Marangoni Effects (generalized symmetry condition)
  !---------------------------------------------------------------------------

  elseif (icodcl(ifac,iu).eq.14) then

    pimpv(1) = rcodcl(ifac,iu,1)
    pimpv(2) = rcodcl(ifac,iv,1)
    pimpv(3) = rcodcl(ifac,iw,1)

    qimpv(1) = rcodcl(ifac,iu,3)
    qimpv(2) = rcodcl(ifac,iv,3)
    qimpv(3) = rcodcl(ifac,iw,3)

    normal(1) = surfbo(1,ifac)/surfbn(ifac)
    normal(2) = surfbo(2,ifac)/surfbn(ifac)
    normal(3) = surfbo(3,ifac)/surfbn(ifac)


    ! Coupled solving of the velocity components

    call set_generalized_sym_vector &
         !=========================
       ( coefau(:,ifac)  , cofafu(:,ifac)  ,             &
         coefbu(:,:,ifac), cofbfu(:,:,ifac),             &
         pimpv           , qimpv            , hint , normal )

  endif

enddo

!===============================================================================
! 10. Pressure: Dirichlet and Neumann and convective outlet
!===============================================================================

call field_get_coefa_s(ivarfl(ipr), coefap)
call field_get_coefb_s(ivarfl(ipr), coefbp)
call field_get_coefaf_s(ivarfl(ipr), cofafp)
call field_get_coefbf_s(ivarfl(ipr), cofbfp)
if (icavit.ge.0)  call field_get_val_s(icrom, crom)

do ifac = 1, nfabor

  iel = ifabor(ifac)

  ! --- Geometrical quantities
  distbf = distb(ifac)

  ! If a flux dt.grad P (W/m2) is set in cs_user_boundary
  if (idften(ipr).eq.1) then
    hint = dt(iel)/distbf
    if (icavit.ge.0)  hint = hint/crom(iel)
    if (ippmod(idarcy).eq.1) hint = permeability(iel)/distbf
  else if (idften(ipr).eq.3) then
    hint = ( dttens(1, iel)*surfbo(1,ifac)**2              &
           + dttens(2, iel)*surfbo(2,ifac)**2              &
           + dttens(3, iel)*surfbo(3,ifac)**2              &
           ) / (surfbn(ifac)**2 * distbf)
    if (icavit.ge.0)  hint = hint/crom(iel)
  ! Symmetric tensor diffusivity
  elseif (idften(ipr).eq.6) then
    if (ippmod(idarcy).eq.-1) then
      visci(1,1) = dttens(1,iel)
      visci(2,2) = dttens(2,iel)
      visci(3,3) = dttens(3,iel)
      visci(1,2) = dttens(4,iel)
      visci(2,1) = dttens(4,iel)
      visci(2,3) = dttens(5,iel)
      visci(3,2) = dttens(5,iel)
      visci(1,3) = dttens(6,iel)
      visci(3,1) = dttens(6,iel)
    else
      visci(1,1) = tensor_permeability(1,iel)
      visci(2,2) = tensor_permeability(2,iel)
      visci(3,3) = tensor_permeability(3,iel)
      visci(1,2) = tensor_permeability(4,iel)
      visci(2,1) = tensor_permeability(4,iel)
      visci(2,3) = tensor_permeability(5,iel)
      visci(3,2) = tensor_permeability(5,iel)
      visci(1,3) = tensor_permeability(6,iel)
      visci(3,1) = tensor_permeability(6,iel)
    endif

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
    if (icavit.ge.0)  hint = hint/crom(iel)

  endif

  ! On doit remodifier la valeur du  Dirichlet de pression de maniere
  !  a retrouver P*. Car dans typecl.f90 on a travaille avec la pression
  !  totale fournie par l'utilisateur :  Ptotale= P*+ rho.g.r
  ! En compressible, on laisse rcodcl tel quel

  ! Dirichlet Boundary Condition
  !-----------------------------

  if (icodcl(ifac,ipr).eq.1) then

    hext = rcodcl(ifac,ipr,2)
    if (ippmod(icompf).ge.0) then
      pimp = rcodcl(ifac,ipr,1)
    else
      pimp = rcodcl(ifac,ipr,1)                              &
           - ro0*( gx*(cdgfbo(1,ifac)-xxp0)                  &
           + gy*(cdgfbo(2,ifac)-xyp0)                  &
           + gz*(cdgfbo(3,ifac)-xzp0) )                &
           + pred0 - p0
    endif

    call set_dirichlet_scalar &
         !===================
       ( coefap(ifac), cofafp(ifac),                         &
         coefbp(ifac), cofbfp(ifac),                         &
         pimp              , hint             , hext )

  endif

  ! Neumann Boundary Conditions
  !----------------------------

  if (icodcl(ifac,ipr).eq.3) then

    qimp = rcodcl(ifac,ipr,3)

    call set_neumann_scalar &
         !=================
       ( coefap(ifac), cofafp(ifac),                         &
         coefbp(ifac), cofbfp(ifac),                         &
         qimp              , hint )

  ! Convective Boundary Conditions
  !-------------------------------

  elseif (icodcl(ifac,ipr).eq.2) then

    pimp = rcodcl(ifac,ipr,1)
    cfl = rcodcl(ifac,ipr,2)

    call set_convective_outlet_scalar &
         !===========================
       ( coefap(ifac), cofafp(ifac),                         &
         coefbp(ifac), cofbfp(ifac),                         &
         pimp              , cfl               , hint )

  ! Imposed value for the convection operator, imposed flux for diffusion
  !----------------------------------------------------------------------

  elseif (icodcl(ifac,ipr).eq.13) then

    pimp = rcodcl(ifac,ipr,1)
    qimp = rcodcl(ifac,ipr,3)

    call set_dirichlet_conv_neumann_diff_scalar &
         !=====================================
       ( coefap(ifac), cofafp(ifac),                         &
         coefbp(ifac), cofbfp(ifac),                         &
         pimp              , qimp )

  endif

enddo

!===============================================================================
! 11. Void fraction (cavitation): Dirichlet and Neumann and convective outlet
!===============================================================================

if (icavit.ge.0) then

  call field_get_coefa_s(ivarfl(ivoidf), coefap)
  call field_get_coefb_s(ivarfl(ivoidf), coefbp)
  call field_get_coefaf_s(ivarfl(ivoidf), cofafp)
  call field_get_coefbf_s(ivarfl(ivoidf), cofbfp)

  do ifac = 1, nfabor

    ! hint in unused since there is no diffusion for the void fraction
    hint = 1.d0

    ! Dirichlet Boundary Condition
    !-----------------------------

    if (icodcl(ifac,ivoidf).eq.1) then

      pimp = rcodcl(ifac,ivoidf,1)
      hext = rcodcl(ifac,ivoidf,2)

      call set_dirichlet_scalar &
      !====================
    ( coefap(ifac), cofafp(ifac),                        &
      coefbp(ifac), cofbfp(ifac),                        &
      pimp              , hint              , hext )

    endif

    ! Neumann Boundary Conditions
    !----------------------------

    if (icodcl(ifac,ivoidf).eq.3) then

      qimp = rcodcl(ifac,ivoidf,3)

      call set_neumann_scalar &
      !==================
    ( coefap(ifac), cofafp(ifac),                        &
      coefbp(ifac), cofbfp(ifac),                        &
      qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

    elseif (icodcl(ifac,ivoidf).eq.2) then

      pimp = rcodcl(ifac,ivoidf,1)
      cfl = rcodcl(ifac,ivoidf,2)

      call set_convective_outlet_scalar &
      !==================
    ( coefap(ifac), cofafp(ifac),                        &
      coefbp(ifac), cofbfp(ifac),                        &
      pimp              , cfl               , hint )

    endif

  enddo

endif

!===============================================================================
! 12. Turbulent quantities: Dirichlet and Neumann and convective outlet
!===============================================================================

! ---> k-epsilon and k-omega

if (itytur.eq.2.or.iturb.eq.60) then

  do ii = 1, 2

    !     Pour le k-omega, on met les valeurs sigma_k2 et sigma_w2 car ce terme
    !     ne concerne en pratique que les entrees (pas de pb en paroi ou en flux
    !     nul)
    if (ii.eq.1.and.itytur.eq.2) then
      ivar   = ik
      sigma  = sigmak
    elseif (ii.eq.1.and.iturb.eq.60) then
      ivar   = ik
      sigma  = ckwsk2 !FIXME it is not consistent with the model
    elseif (itytur.eq.2) then
      ivar   = iep
      sigma  = sigmae
    else
      ivar   = iomg
      sigma  = ckwsw2 !FIXME it is not consistent with the model
    endif

    call field_get_coefa_s(ivarfl(ivar), coefap)
    call field_get_coefb_s(ivarfl(ivar), coefbp)
    call field_get_coefaf_s(ivarfl(ivar), cofafp)
    call field_get_coefbf_s(ivarfl(ivar), cofbfp)

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Physical Propreties
      visclc = viscl(iel)
      visctc = visct(iel)

      ! --- Geometrical quantities
      distbf = distb(ifac)

      hint = (visclc+visctc/sigma)/distbf

      ! Dirichlet Boundary Condition
      !-----------------------------

      if (icodcl(ifac,ivar).eq.1) then

        pimp = rcodcl(ifac,ivar,1)
        hext = rcodcl(ifac,ivar,2)

        call set_dirichlet_scalar &
             !====================
           ( coefap(ifac), cofafp(ifac),                        &
             coefbp(ifac), cofbfp(ifac),                        &
             pimp              , hint              , hext )

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                        &
             coefbp(ifac), cofbfp(ifac),                        &
             qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                        &
             coefbp(ifac), cofbfp(ifac),                        &
             pimp              , cfl               , hint )

      ! Imposed value for the convection operator, imposed flux for diffusion
      !----------------------------------------------------------------------

      elseif (icodcl(ifac,ivar).eq.13) then

        pimp = rcodcl(ifac,ivar,1)
        qimp = rcodcl(ifac,ivar,3)

        call set_dirichlet_conv_neumann_diff_scalar &
             !=====================================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , qimp )


      endif

    enddo

  enddo

! ---> Rij-epsilon
elseif (itytur.eq.3) then

  ! --> Rij

  do isou = 1, 6

    if(isou.eq.1) then
      ivar   = ir11
    elseif(isou.eq.2) then
      ivar   = ir22
    elseif(isou.eq.3) then
      ivar   = ir33
    elseif(isou.eq.4) then
      ivar   = ir12
    elseif(isou.eq.5) then
      ivar   = ir13
    elseif(isou.eq.6) then
      ivar   = ir23
    endif

    call field_get_coefa_s(ivarfl(ivar), coefap)
    call field_get_coefb_s(ivarfl(ivar), coefbp)
    call field_get_coefaf_s(ivarfl(ivar), cofafp)
    call field_get_coefbf_s(ivarfl(ivar), cofbfp)
    call field_get_coefad_s(ivarfl(ivar), cofadp)
    call field_get_coefbd_s(ivarfl(ivar), cofbdp)

    if (idften(ivar).eq.6) call field_get_val_v(ivsten, visten)

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Physical Propreties
      visclc = viscl(iel)

      ! --- Geometrical quantities
      distbf = distb(ifac)

      ! Symmetric tensor diffusivity (Daly Harlow -- GGDH)
      if (idften(ivar).eq.6) then

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
        visctc = visct(iel)
        hint = (visclc+visctc*csrij/cmu)/distbf
      endif

      ! Dirichlet Boundary Condition
      !-----------------------------

      if (icodcl(ifac,ivar).eq.1) then
        pimp = rcodcl(ifac,ivar,1)
        hext = rcodcl(ifac,ivar,2)

        call set_dirichlet_scalar &
             !====================
           ( coefap(ifac), cofafp(ifac),                        &
             coefbp(ifac), cofbfp(ifac),                        &
             pimp              , hint              , hext )

        ! Boundary conditions for the momentum equation
        cofadp(ifac) = coefap(ifac)
        cofbdp(ifac) = coefbp(ifac)
      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then
        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                        &
             coefbp(ifac), cofbfp(ifac),                        &
             qimp              , hint )

        ! Boundary conditions for the momentum equation
        cofadp(ifac) = coefap(ifac)
        cofbdp(ifac) = coefbp(ifac)

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then
        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                        &
             coefbp(ifac), cofbfp(ifac),                        &
             pimp              , cfl               , hint )
        ! Boundary conditions for the momentum equation
        cofadp(ifac) = coefap(ifac)
        cofbdp(ifac) = coefbp(ifac)

      ! Imposed value for the convection operator, imposed flux for diffusion
      !----------------------------------------------------------------------

      elseif (icodcl(ifac,ivar).eq.13) then

        pimp = rcodcl(ifac,ivar,1)
        qimp = rcodcl(ifac,ivar,3)

        call set_dirichlet_conv_neumann_diff_scalar &
             !=====================================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , qimp )

        ! Boundary conditions for the momentum equation
        cofadp(ifac) = coefap(ifac)
        cofbdp(ifac) = coefbp(ifac)

      endif

    enddo

  enddo

  ! --> epsilon

  ivar   = iep

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)
  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  if (idften(ivar).eq.6) call field_get_val_v(ivsten, visten)

  do ifac = 1, nfabor

    iel = ifabor(ifac)

    ! --- Physical Propreties
    visclc = viscl(iel)
    visctc = visct(iel)

    ! --- Geometrical quantities
    distbf = distb(ifac)

    ! Symmetric tensor diffusivity (Daly Harlow -- GGDH)
    if (idften(ivar).eq.6) then

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

    ! Dirichlet Boundary Condition
    !-----------------------------

    if (icodcl(ifac,ivar).eq.1) then

      pimp = rcodcl(ifac,ivar,1)
      hext = rcodcl(ifac,ivar,2)

      call set_dirichlet_scalar &
           !====================
         ( coefap(ifac), cofafp(ifac),                        &
           coefbp(ifac), cofbfp(ifac),                        &
           pimp              , hint              , hext )

    endif

    ! Neumann Boundary Conditions
    !----------------------------

    if (icodcl(ifac,ivar).eq.3) then

      qimp = rcodcl(ifac,ivar,3)

      call set_neumann_scalar &
           !==================
         ( coefap(ifac), cofafp(ifac),                        &
           coefbp(ifac), cofbfp(ifac),                        &
           qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                      &
             coefbp(ifac), cofbfp(ifac),                      &
             pimp              , cfl               , hint )


      ! Imposed value for the convection operator, imposed flux for diffusion
      !----------------------------------------------------------------------

      elseif (icodcl(ifac,ivar).eq.13) then

        pimp = rcodcl(ifac,ivar,1)
        qimp = rcodcl(ifac,ivar,3)

        call set_dirichlet_conv_neumann_diff_scalar &
             !=====================================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , qimp )

    endif

  enddo

  ! --> alpha for the EBRSM

  if (iturb.eq.32) then
    ivar   = ial

    call field_get_coefa_s(ivarfl(ivar), coefap)
    call field_get_coefb_s(ivarfl(ivar), coefbp)
    call field_get_coefaf_s(ivarfl(ivar), cofafp)
    call field_get_coefbf_s(ivarfl(ivar), cofbfp)

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      distbf = distb(ifac)

      hint = 1.d0/distbf

      ! Dirichlet Boundary Condition
      !-----------------------------

      if (icodcl(ifac,ivar).eq.1) then

        pimp = rcodcl(ifac,ivar,1)
        hext = rcodcl(ifac,ivar,2)

        call set_dirichlet_scalar &
             !====================
           ( coefap(ifac), cofafp(ifac),                        &
             coefbp(ifac), cofbfp(ifac),                        &
             pimp              , hint              , hext )

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                        &
             coefbp(ifac), cofbfp(ifac),                        &
             qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                        &
             coefbp(ifac), cofbfp(ifac),                        &
             pimp              , cfl               , hint )

      ! Imposed value for the convection operator, imposed flux for diffusion
      !----------------------------------------------------------------------

      elseif (icodcl(ifac,ivar).eq.13) then

        pimp = rcodcl(ifac,ivar,1)
        qimp = rcodcl(ifac,ivar,3)

        call set_dirichlet_conv_neumann_diff_scalar &
             !=====================================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , qimp )


      endif

    enddo
  endif

! ---> v2f type models (phi_bar and Bl-v2/k)

elseif (itytur.eq.5) then

  !   --> k, epsilon  and phi
  do ii = 1, 3

    if (ii.eq.1) then
      ivar   = ik
      sigma  = sigmak
    elseif (ii.eq.2) then
      ivar   = iep
      sigma  = sigmae
    else
      ivar   = iphi
      sigma  = sigmak
    endif

    call field_get_coefa_s(ivarfl(ivar), coefap)
    call field_get_coefb_s(ivarfl(ivar), coefbp)
    call field_get_coefaf_s(ivarfl(ivar), cofafp)
    call field_get_coefbf_s(ivarfl(ivar), cofbfp)

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Physical Propreties
      visclc = viscl(iel)
      visctc = visct(iel)

      ! --- Geometrical quantities
      distbf = distb(ifac)

      hint = (visclc+visctc/sigma)/distbf

      ! Dirichlet Boundary Condition
      !-----------------------------

      if (icodcl(ifac,ivar).eq.1) then

        pimp = rcodcl(ifac,ivar,1)
        hext = rcodcl(ifac,ivar,2)

        call set_dirichlet_scalar &
             !====================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , hint              , hext )

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , cfl               , hint )

      ! Imposed value for the convection operator, imposed flux for diffusion
      !----------------------------------------------------------------------

      elseif (icodcl(ifac,ivar).eq.13) then

        pimp = rcodcl(ifac,ivar,1)
        qimp = rcodcl(ifac,ivar,3)

        call set_dirichlet_conv_neumann_diff_scalar &
             !=====================================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , qimp )


      endif

    enddo

  enddo

  if (iturb.eq.50) then

    ! --> FB

    ivar   = ifb

    call field_get_coefa_s(ivarfl(ivar), coefap)
    call field_get_coefb_s(ivarfl(ivar), coefbp)
    call field_get_coefaf_s(ivarfl(ivar), cofafp)
    call field_get_coefbf_s(ivarfl(ivar), cofbfp)

    do ifac = 1, nfabor

      ! --- Physical Propreties
      visclc = 1.d0

      ! --- Geometrical quantities
      distbf = distb(ifac)

      hint = visclc/distbf

      ! Dirichlet Boundary Condition
      !-----------------------------

      if (icodcl(ifac,ivar).eq.1) then

        pimp = rcodcl(ifac,ivar,1)
        hext = rcodcl(ifac,ivar,2)

        call set_dirichlet_scalar &
             !====================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , hint              , hext )

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , cfl               , hint )

      ! Imposed value for the convection operator, imposed flux for diffusion
      !----------------------------------------------------------------------

      elseif (icodcl(ifac,ivar).eq.13) then

        pimp = rcodcl(ifac,ivar,1)
        qimp = rcodcl(ifac,ivar,3)

        call set_dirichlet_conv_neumann_diff_scalar &
             !=====================================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , qimp )


      endif

    enddo

  elseif (iturb.eq.51) then

    ! --> alpha

    ivar   = ial

    call field_get_coefa_s(ivarfl(ivar), coefap)
    call field_get_coefb_s(ivarfl(ivar), coefbp)
    call field_get_coefaf_s(ivarfl(ivar), cofafp)
    call field_get_coefbf_s(ivarfl(ivar), cofbfp)

    do ifac = 1, nfabor

      ! --- Physical Propreties
      visclc = 1.d0

      ! --- Geometrical quantities
      distbf = distb(ifac)

      hint = visclc/distbf

      ! Dirichlet Boundary Condition
      !-----------------------------

      if (icodcl(ifac,ivar).eq.1) then

        pimp = rcodcl(ifac,ivar,1)
        hext = rcodcl(ifac,ivar,2)

        call set_dirichlet_scalar &
             !====================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , hint              , hext )

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , cfl               , hint )

      ! Imposed value for the convection operator, imposed flux for diffusion
      !----------------------------------------------------------------------

      elseif (icodcl(ifac,ivar).eq.13) then

        pimp = rcodcl(ifac,ivar,1)
        qimp = rcodcl(ifac,ivar,3)

        call set_dirichlet_conv_neumann_diff_scalar &
             !=====================================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , qimp )


      endif

    enddo

  endif

! ---> Spalart Allmaras

elseif (iturb.eq.70) then

  ivar   = inusa

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)
  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  do ifac = 1, nfabor

    iel = ifabor(ifac)

    ! --- Physical Propreties
    visclc = viscl(iel)

    ! --- Geometrical quantities
    distbf = distb(ifac)

    hint = visclc/distbf

    ! Dirichlet Boundary Condition
    !-----------------------------

    if (icodcl(ifac,ivar).eq.1) then

      pimp = rcodcl(ifac,ivar,1)
      hext = rcodcl(ifac,ivar,2)

      call set_dirichlet_scalar &
           !====================
         ( coefap(ifac), cofafp(ifac),                         &
           coefbp(ifac), cofbfp(ifac),                         &
           pimp              , hint              , hext )

    endif

    ! Neumann Boundary Conditions
    !----------------------------

    if (icodcl(ifac,ivar).eq.3) then

      qimp = rcodcl(ifac,ivar,3)

      call set_neumann_scalar &
           !==================
         ( coefap(ifac), cofafp(ifac),                         &
           coefbp(ifac), cofbfp(ifac),                         &
           qimp              , hint )

    ! Convective Boundary Conditions
    !-------------------------------

    elseif (icodcl(ifac,ivar).eq.2) then

      pimp = rcodcl(ifac,ivar,1)
      cfl = rcodcl(ifac,ivar,2)

      call set_convective_outlet_scalar &
           !==================
         ( coefap(ifac), cofafp(ifac),                         &
           coefbp(ifac), cofbfp(ifac),                         &
           pimp              , cfl               , hint )

    ! Imposed value for the convection operator, imposed flux for diffusion
    !----------------------------------------------------------------------

    elseif (icodcl(ifac,ivar).eq.13) then

      pimp = rcodcl(ifac,ivar,1)
      qimp = rcodcl(ifac,ivar,3)

      call set_dirichlet_conv_neumann_diff_scalar &
           !=====================================
         ( coefap(ifac), cofafp(ifac),                         &
           coefbp(ifac), cofbfp(ifac),                         &
           pimp              , qimp )


    endif

  enddo

endif

!===============================================================================
! 13. Other scalars (except variances):
!     Dirichlet and Neumann and convectiv outlet
!===============================================================================

if (nscal.ge.1) then

  if(icp.gt.0) then
    call field_get_val_s(iprpfl(icp), cpro_cp)
  endif

  if (ipccv.gt.0) then
    call field_get_val_s(iprpfl(ipccv), cpro_cv)
  endif

  do ii = 1, nscal

    ivar   = isca(ii)

    isvhbl = 0
    if (ii.eq.isvhb) then
      isvhbl = isvhb
    endif

    call field_get_key_int (ivarfl(ivar), kivisl, ifcvsl)
    if (ifcvsl .ge. 0) then
      call field_get_val_s(ifcvsl, viscls)
    endif

    ! --- Indicateur de prise en compte de Cp ou non
    !       (selon si le scalaire (scalaire associe pour une fluctuation)
    !        doit etre ou non traite comme une temperature)
    !      Si le scalaire est une variance et que le
    !        scalaire associe n'est pas resolu, on suppose alors qu'il
    !        doit etre traite comme un scalaire passif (defaut IHCP = 0)
    ihcp = 0

    iscal = ii
    if (iscavr(ii).gt.0) then
      iscal = iscavr(ii)
    endif

    if (iscacp(iscal).eq.1) then
      if(icp.gt.0) then
        ihcp = 2
      else
        ihcp = 1
      endif
    endif

    call field_get_coefa_s(ivarfl(ivar), coefap)
    call field_get_coefb_s(ivarfl(ivar), coefbp)
    call field_get_coefaf_s(ivarfl(ivar), cofafp)
    call field_get_coefbf_s(ivarfl(ivar), cofbfp)

    if (idften(ivar).eq.6.or.ityturt(ii).eq.3) then
      if (iturb.ne.32.or.ityturt(ii).eq.3) then
        call field_get_val_v(ivsten, visten)
      else ! EBRSM and (GGDH or AFM)
        call field_get_val_v(ivstes, visten)
      endif
    endif

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Physical Properties
      visctc = visct(iel)

      ! --- Geometrical quantities
      distbf = distb(ifac)

      ! --- Prise en compte de Cp ou CV
      !      (dans le Cas compressible ihcp=0)

      cpp = 1.d0
      if (ihcp.eq.0) then
        cpp = 1.d0
      elseif (ihcp.eq.2) then
        cpp = cpro_cp(iel)
      elseif (ihcp.eq.1) then
        cpp = cp0
      endif

      ! --- Viscosite variable ou non
      if (ifcvsl.lt.0) then
        rkl = visls0(ii)
      else
        rkl = viscls(iel)
      endif

      ! Scalar diffusivity
      if (idften(ivar).eq.1) then
        if (ippmod(idarcy).eq.-1) then !FIXME
          hint = (rkl+idifft(ivar)*cpp*visctc/sigmas(ii))/distbf
        else ! idarcy = 1
          hint = rkl/distbf
        endif

      ! Symmetric tensor diffusivity
      elseif (idften(ivar).eq.6) then

        temp = idifft(ivar)*cpp*ctheta(ii)/csrij
        visci(1,1) = rkl + temp*visten(1,iel)
        visci(2,2) = rkl + temp*visten(2,iel)
        visci(3,3) = rkl + temp*visten(3,iel)
        visci(1,2) =       temp*visten(4,iel)
        visci(2,1) =       temp*visten(4,iel)
        visci(2,3) =       temp*visten(5,iel)
        visci(3,2) =       temp*visten(5,iel)
        visci(1,3) =       temp*visten(6,iel)
        visci(3,1) =       temp*visten(6,iel)

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

      endif

      ! Dirichlet Boundary Condition
      !-----------------------------

      if (icodcl(ifac,ivar).eq.1) then

        pimp = rcodcl(ifac,ivar,1)
        hext = rcodcl(ifac,ivar,2)

        call set_dirichlet_scalar &
             !====================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , hint              , hext )

        ! ---> COUPLAGE : on stocke le hint (lambda/d      en temperature,
        !                                    lambda/(cp d) en enthalpie,
        !                                    lambda/(cv d) en energie)
        !FIXME useless
        if (isvhbl.gt.0) then
          hbord(ifac) = hint
        endif

        !--> Rayonnement :

        !      On stocke le coefficient d'echange lambda/distance
        !      (ou son equivalent en turbulent) quelle que soit la
        !      variable thermique transportee (temperature ou enthalpie)
        !      car on l'utilise pour realiser des bilans aux parois qui
        !      sont faits en temperature (on cherche la temperature de
        !      paroi quelle que soit la variable thermique transportee pour
        !      ecrire des eps sigma T4).

        !     donc :

        !       lorsque la variable transportee est la temperature
        !         ABS(ISCSTH(II)).EQ.1 : RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) = HINT
        !         puisque HINT = VISLS * CP / DISTBR
        !                      = lambda/distance en W/(m2 K)

        !       lorsque la variable transportee est l'enthalpie
        !         ISCSTH(II).EQ.2 : RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) = HINT*CPR
        !         avec
        !            IF(ICP.GT.0) THEN
        !              CPR = CPRO_CP(IEL)
        !            ELSE
        !              CPR = CP0
        !            ENDIF
        !         puisque HINT = VISLS / DISTBR
        !                      = lambda/(CP * distance)

        !       lorsque la variable transportee est l'energie
        !         ISCSTH(II).EQ.3 :
        !         on procede comme pour l'enthalpie avec CV au lieu de CP
        !         (rq : il n'y a pas d'hypothese, sf en non orthogonal :
        !               le flux est le bon et le coef d'echange aussi)

        !      De meme plus bas et de meme dans clptur.

        !               Si on rayonne et que
        !                  le scalaire est la variable energetique

        if (iirayo.ge.1 .and. ii.eq.iscalt) then

          ! We compute the exchange coefficient in W/(m2 K)

          ! Enthalpy
          if (itherm.eq.2) then
            ! If Cp is variable
            if (icp.gt.0) then
              bhconv(ifac) = hint*cpro_cp(iel)
            else
              bhconv(ifac) = hint*cp0
            endif

          ! Energie (compressible module)
          elseif (itherm.eq.3) then
            ! If Cv is variable
            if (ipccv.gt.0) then
              bhconv(ifac) = hint*cpro_cv(iel)
            else
              bhconv(ifac) = hint*cv0
            endif

          ! Temperature
          elseif(iscacp(ii).eq.1) then
            bhconv(ifac) = hint
          endif

          ! The outgoing flux is stored (Q = h(Ti'-Tp): negative if
          !  gain for the fluid) in W/m2
          bfconv(ifac) = cofafp(ifac) + cofbfp(ifac)*theipb(ifac)

        endif

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             qimp              , hint )

        if (isvhbl.gt.0) hbord(ifac) = hint

        !--> Rayonnement :

        if (iirayo.ge.1 .and. ii.eq.iscalt) then

          ! We compute the exchange coefficient in W/(m2 K)

          ! Enthalpy
          if (itherm.eq.2) then
            ! If Cp is variable
            if (icp.gt.0) then
              bhconv(ifac) = hint*cpro_cp(iel)
            else
              bhconv(ifac) = hint*cp0
            endif

          ! Energy (compressible module)
          elseif (itherm.eq.3) then
            ! If Cv is variable
            if (ipccv.gt.0) then
              bhconv(ifac) = hint*cpro_cv(iel)
            else
              bhconv(ifac) = hint*cv0
            endif

          ! Temperature
          elseif (itherm.eq.1) then
            bhconv(ifac) = hint
          endif

          ! The outgoing flux is stored (Q = h(Ti'-Tp): negative if
          !  gain for the fluid) in W/m2
          bfconv(ifac) = rcodcl(ifac,ivar,3)
        endif

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , cfl               , hint )

      ! Imposed value for the convection operator, imposed flux for diffusion
      !----------------------------------------------------------------------

      elseif (icodcl(ifac,ivar).eq.13) then

        pimp = rcodcl(ifac,ivar,1)
        qimp = rcodcl(ifac,ivar,3)

        call set_dirichlet_conv_neumann_diff_scalar &
             !=====================================
           ( coefap(ifac), cofafp(ifac),                         &
             coefbp(ifac), cofbfp(ifac),                         &
             pimp              , qimp )


      endif

      ! Thermal heat flux boundary conditions
      if (ityturt(ii).eq.3) then

        ! Name of the scalar ivar !TODO move outside of the loop
        call field_get_name(ivarfl(ivar), fname)

        ! Index of the corresponding turbulent flux
        call field_get_id(trim(fname)//'_turbulent_flux', f_id)

        call field_get_coefa_v(f_id,coefaut)
        call field_get_coefb_v(f_id,coefbut)
        call field_get_coefaf_v(f_id,cofafut)
        call field_get_coefbf_v(f_id,cofbfut)
        call field_get_coefad_v(f_id,cofarut)
        call field_get_coefbd_v(f_id,cofbrut)

        ! --- Physical Propreties
        visclc = viscl(iel)

        ! --- Geometrical quantities
        distbf = distb(ifac)

        if (ifcvsl.lt.0) then
          rkl = visls0(iscal)/cpp
        else
          rkl = viscls(iel)/cpp
        endif
        hintt(1) = 0.5d0*(visclc+rkl)/distbf                        &
                 + visten(1,iel)*ctheta(iscal)/distbf/csrij !FIXME ctheta (iscal)
        hintt(2) = 0.5d0*(visclc+rkl)/distbf                        &
                 + visten(2,iel)*ctheta(iscal)/distbf/csrij
        hintt(3) = 0.5d0*(visclc+rkl)/distbf                        &
                 + visten(3,iel)*ctheta(iscal)/distbf/csrij
        hintt(4) = visten(4,iel)*ctheta(iscal)/distbf/csrij
        hintt(5) = visten(5,iel)*ctheta(iscal)/distbf/csrij
        hintt(6) = visten(6,iel)*ctheta(iscal)/distbf/csrij

        ! Set pointer values of turbulent fluxes in icodcl
        iut = nvar + 3*(ifltur(ii) - 1) + 1
        ivt = nvar + 3*(ifltur(ii) - 1) + 2
        iwt = nvar + 3*(ifltur(ii) - 1) + 3

        ! Dirichlet Boundary Condition
        !-----------------------------

        if (icodcl(ifac,iut).eq.1) then

          pimpv(1) = rcodcl(ifac,iut,1)
          pimpv(2) = rcodcl(ifac,ivt,1)
          pimpv(3) = rcodcl(ifac,iwt,1)
          hextv(1) = rcodcl(ifac,iut,2)
          hextv(2) = rcodcl(ifac,ivt,2)
          hextv(3) = rcodcl(ifac,iwt,2)

          call set_dirichlet_vector_ggdh &
               !========================
             ( coefaut(:,ifac)  , cofafut(:,ifac)  ,           &
               coefbut(:,:,ifac), cofbfut(:,:,ifac),           &
               pimpv            , hintt            , hextv )

          ! Boundary conditions for thermal transport equation
          do isou = 1, 3
            cofarut(isou,ifac) = coefaut(isou,ifac)
            do jsou =1, 3
              cofbrut(isou,jsou,ifac) = coefbut(isou,jsou,ifac)
            enddo
          enddo

        ! Neumann Boundary Conditions
        !----------------------------

        elseif (icodcl(ifac,iut).eq.3) then

          qimpv(1) = rcodcl(ifac,iut,3)
          qimpv(2) = rcodcl(ifac,ivt,3)
          qimpv(3) = rcodcl(ifac,iwt,3)

          call set_neumann_vector_ggdh &
          !===========================
             ( coefaut(:,ifac)  , cofafut(:,ifac)  ,           &
               coefbut(:,:,ifac), cofbfut(:,:,ifac),           &
               qimpv            , hintt )

          ! Boundary conditions for thermal transport equation
          do isou = 1, 3
            cofarut(isou,ifac) = coefaut(isou,ifac)
            do jsou =1, 3
              cofbrut(isou,jsou,ifac) = coefbut(isou,jsou,ifac)
            enddo
          enddo

        ! Convective Boundary Conditions
        !-------------------------------

        elseif (icodcl(ifac,iut).eq.2) then

          pimpv(1) = rcodcl(ifac,iut,1)
          cflv(1) = rcodcl(ifac,iut,2)
          pimpv(2) = rcodcl(ifac,ivt,1)
          cflv(2) = rcodcl(ifac,ivt,2)
          pimpv(3) = rcodcl(ifac,iwt,1)
          cflv(3) = rcodcl(ifac,iwt,2)

          call set_convective_outlet_vector_ggdh &
          !=====================================
             ( coefaut(:,ifac)  , cofafut(:,ifac)  ,           &
               coefbut(:,:,ifac), cofbfut(:,:,ifac),           &
               pimpv            , cflv             , hintt )

          ! Boundary conditions for thermal transport equation
          do isou = 1, 3
            cofarut(isou,ifac) = coefaut(isou,ifac)
            do jsou =1, 3
              cofbrut(isou,jsou,ifac) = coefbut(isou,jsou,ifac)
            enddo
          enddo

        endif

      endif

    enddo

  enddo

endif

!===============================================================================
! 14. Mesh velocity (ALE module): Dirichlet and Neumann and convective outlet
!===============================================================================

if (iale.eq.1) then

  call field_get_coefa_v(ivarfl(iuma), claale)
  call field_get_coefb_v(ivarfl(iuma), clbale)
  call field_get_coefaf_v(ivarfl(iuma), cfaale)
  call field_get_coefbf_v(ivarfl(iuma), cfbale)

  call field_get_val_s(iprpfl(ivisma(1)), cpro_visma1)
  call field_get_val_s(iprpfl(ivisma(2)), cpro_visma2)
  call field_get_val_s(iprpfl(ivisma(3)), cpro_visma3)

  do ifac = 1, nfabor

    iel = ifabor(ifac)
    distbf = distb(ifac)
    srfbn2 = surfbn(ifac)**2
    if (iortvm.eq.0) then
      hint = cpro_visma1(iel)/distbf
    else
      hint = ( cpro_visma1(iel)*surfbo(1,ifac)**2    &
             + cpro_visma2(iel)*surfbo(2,ifac)**2    &
             + cpro_visma3(iel)*surfbo(3,ifac)**2 )  &
           /distbf/srfbn2
    endif

    ! Dirichlet Boundary Conditions
    !------------------------------

    if (icodcl(ifac,iuma).eq.1) then

      pimpv(1) = rcodcl(ifac,iuma,1)
      pimpv(2) = rcodcl(ifac,ivma,1)
      pimpv(3) = rcodcl(ifac,iwma,1)
      hextv(1) = rcodcl(ifac,iuma,2)
      hextv(2) = rcodcl(ifac,ivma,2)
      hextv(3) = rcodcl(ifac,iwma,2)

      call set_dirichlet_vector &
           !====================
         ( claale(:,ifac)  , cfaale(:,ifac)  ,             &
           clbale(:,:,ifac), cfbale(:,:,ifac),             &
           pimpv           , hint            , hextv )

    ! Neumann Boundary Conditions
    !----------------------------

    elseif (icodcl(ifac,iuma).eq.3) then

      ! Coupled solving of the velocity components

      qimpv(1) = rcodcl(ifac,iuma,3)
      qimpv(2) = rcodcl(ifac,ivma,3)
      qimpv(3) = rcodcl(ifac,iwma,3)

      call set_neumann_vector &
           !==================
         ( claale(:,ifac)  , cfaale(:,ifac)  ,             &
           clbale(:,:,ifac), cfbale(:,:,ifac),             &
           qimpv           , hint )


    ! Convective Boundary Conditions
    !-------------------------------

    elseif (icodcl(ifac,iuma).eq.2) then

      ! Coupled solving of the velocity components

      pimpv(1) = rcodcl(ifac,iuma,1)
      cflv(1) = rcodcl(ifac,iuma,2)
      pimpv(2) = rcodcl(ifac,ivma,1)
      cflv(2) = rcodcl(ifac,ivma,2)
      pimpv(3) = rcodcl(ifac,iwma,1)
      cflv(3) = rcodcl(ifac,iwma,2)

      call set_convective_outlet_vector &
           !==================
         ( claale(:,ifac)  , cfaale(:,ifac)  ,             &
           clbale(:,:,ifac), cfbale(:,:,ifac),             &
           pimpv           , cflv            , hint )

    endif

  enddo

endif

!===============================================================================
! 15. Compute stresses at boundary (step 1 over 5)
!===============================================================================

if (ineedf.eq.1 .and. iterns.eq.1) then

  ! Coupled solving of the velocity components
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    visclc = viscl(iel)
    visctc = visct(iel)
    if (itytur.eq.3) then
      vistot = visclc
    else
      vistot = visclc + visctc
    endif
    distbf = distb(ifac)
    srfbnf = surfbn(ifac)

    ! The implicit term is added after having updated the velocity
    forbr(1,ifac) = ( cofafu(1,ifac)                              &
                    + cofbfu(1,1,ifac) * velipb(ifac,1)           &
                    + cofbfu(1,2,ifac) * velipb(ifac,2)           &
                    + cofbfu(1,3,ifac) * velipb(ifac,3) )*srfbnf
    forbr(2,ifac) = ( cofafu(2,ifac)                              &
                    + cofbfu(2,1,ifac) * velipb(ifac,1)           &
                    + cofbfu(2,2,ifac) * velipb(ifac,2)           &
                    + cofbfu(2,3,ifac) * velipb(ifac,3) )*srfbnf
    forbr(3,ifac) = ( cofafu(3,ifac)                              &
                    + cofbfu(3,1,ifac) * velipb(ifac,1)           &
                    + cofbfu(3,2,ifac) * velipb(ifac,2)           &
                    + cofbfu(3,3,ifac) * velipb(ifac,3) )*srfbnf
  enddo
endif

! Free memory
deallocate(velipb)
if (allocated(rijipb)) deallocate(rijipb)

!===============================================================================
! 16. Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,&
'@    =========                                               ',/,&
'@      INCOHERENCE ENTRE OPTIONS DE CALCUL ET COND. LIM.     ',/,&
'@                                                            ',/,&
'@      la prise en compte des termes d''echo de paroi        ',/,&
'@      du modele de turbulence Rij-epsilon est activee       ',/,&
'@      IRIJEC = ',I10,'                                      ',/,&
'@      Ou bien l amortissement de la viscosite turbulente    ',/,&
'@      est active IDRIES = ',I10,'en LES                     ',/,&
'@    mais aucune face de bord de type paroi n''est detectee. ',/,&
'@    L''incoherence indiquee ci-dessus n''est pas bloquante  ',/,&
'@      mais peut resulter d''une erreur lors de la           ',/,&
'@      specification des conditions aux limites.             ',/,&
'@                                                            ',/,&
'@    Par securite, le calcul ne sera pas execute.            ',/,&
'@                                                            ',/,&
'@    Verifier les conditions aux limites dans                ',/,&
'@      cs_user_boundary si le domaine comporte des parois.   ',/,&
'@    Eliminer l''option IRIJEC de usipsu si le domaine ne    ',/,&
'@      comporte pas de paroi (ou conditions ICODCL = 5 en    ',/,&
'@      vitesse).                                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                           &
 'Debit entrant retenu en ',I10   ,                  &
                                      ' faces de sortie sur ',I10)
#else

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE BOUNDARY CONDITIONS SPECIFICATION ',/,&
'@    =========                                               ',/,&
'@      INCOHERENCY BETWEEN CALCULATION OPTIONS AND BOUND COND',/,&
'@                                                            ',/,&
'@      The wall-echo terms of the Rij-epsilon turbulence     ',/,&
'@      model are taken into account                          ',/,&
'@      IRIJEC = ',I10,'                                      ',/,&
'@      Or the Van Driest damping of the turbulent viscosity  ',/,&
'@      is active IDRIES = ',I10,'in LES                      ',/,&
'@    but no wall boundary face is detected.                  ',/,&
'@    This incoherency is not blocking but may result from    ',/,&
'@      an error during the boundary conditions               ',/,&
'@      specification.                                        ',/,&
'@                                                            ',/,&
'@    By safety, the calculation will not be run.             ',/,&
'@                                                            ',/,&
'@    Verify the boundary conditions in cs_user_boundary in   ',/,&
'@      if the domain has any walls.                          ',/,&
'@    Remove the option IRIJEC from usipsu if the domain does ',/,&
'@      not have any wall (or conditions ICODCL = 5 for the   ',/,&
'@      velocity).                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                           &
 'Incoming flow detained for ', I10   ,              &
                                          ' outlet faces on ',I10)

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
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimp          Dirichlet value to impose
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     hext          External exchange coefficient (10^30 by default)
!_______________________________________________________________________________

subroutine set_dirichlet_scalar &
 ( coefa , cofaf, coefb , cofbf, pimp  , hint, hext)

!===============================================================================
! Module files
!===============================================================================

use cstnum, only: rinfin

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, pimp, hint, hext

! Local variables

double precision heq

!===============================================================================

if (abs(hext).gt.rinfin*0.5d0) then

  ! Gradient BCs
  coefa = pimp
  coefb = 0.d0

  ! Flux BCs
  cofaf = -hint*pimp
  cofbf =  hint

else

  ! Gradient BCs
  coefa = hext*pimp/(hint + hext)
  coefb = hint     /(hint + hext)

  ! Flux BCs
  heq = hint*hext/(hint + hext)
  cofaf = -heq*pimp
  cofbf =  heq

endif

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     hextv         External exchange coefficient (10^30 by default)
!_______________________________________________________________________________

subroutine set_dirichlet_vector &
 ( coefa , cofaf, coefb , cofbf, pimpv  , hint , hextv)

!===============================================================================
! Module files
!===============================================================================

use cstnum, only: rinfin

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision pimpv(3)
double precision hint
double precision hextv(3)

! Local variables

integer          isou  , jsou
double precision heq

!===============================================================================

do isou = 1, 3

  if (abs(hextv(isou)).gt.rinfin*0.5d0) then
    ! Gradient BCs
    coefa(isou) = pimpv(isou)
    do jsou = 1, 3
      coefb(isou,jsou) = 0.d0
    enddo

    ! Flux BCs
    cofaf(isou) = -hint*pimpv(isou)
    do jsou = 1, 3
      if (jsou.eq.isou) then
        cofbf(isou,jsou) = hint
      else
        cofbf(isou,jsou) = 0.d0
      endif
    enddo

  else

    heq = hint*hextv(isou)/(hint + hextv(isou))

    ! Gradient BCs
    coefa(isou) = hextv(isou)*pimpv(isou)/(hint + hextv(isou))
    do jsou = 1, 3
      if (jsou.eq.isou) then
        coefb(isou,jsou) = hint/(hint + hextv(isou))
      else
        coefb(isou,jsou) = 0.d0
      endif
    enddo

    ! Flux BCs
    cofaf(isou) = -heq*pimpv(isou)
    do jsou = 1, 3
      if (jsou.eq.isou) then
        cofbf(isou,jsou) = heq
      else
        cofbf(isou,jsou) = 0.d0
      endif
    enddo

  endif

enddo

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     hextv         External exchange coefficient (10^30 by default)
!_______________________________________________________________________________

subroutine set_dirichlet_vector_ggdh &
 ( coefa , cofaf, coefb , cofbf, pimpv  , hint , hextv)

!===============================================================================
! Module files
!===============================================================================

use cstnum, only: rinfin

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision pimpv(3)
double precision hint(6)
double precision hextv(3)

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  if (abs(hextv(isou)).gt.rinfin*0.5d0) then
    ! Gradient BCs
    coefa(isou) = pimpv(isou)
    do jsou = 1, 3
      coefb(isou,jsou) = 0.d0
    enddo

  else

    call csexit(1)

  endif

enddo

! Flux BCs
cofaf(1) = -(hint(1)*pimpv(1) + hint(4)*pimpv(2) + hint(6)*pimpv(3))
cofaf(2) = -(hint(4)*pimpv(1) + hint(2)*pimpv(2) + hint(5)*pimpv(3))
cofaf(3) = -(hint(6)*pimpv(1) + hint(5)*pimpv(2) + hint(3)*pimpv(3))
cofbf(1,1) = hint(1)
cofbf(2,2) = hint(2)
cofbf(3,3) = hint(3)
cofbf(1,2) = hint(4)
cofbf(2,1) = hint(4)
cofbf(2,3) = hint(5)
cofbf(3,2) = hint(5)
cofbf(1,3) = hint(6)
cofbf(3,1) = hint(6)

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     qimp          Flux value to impose
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_neumann_scalar &
 ( coefa , cofaf, coefb , cofbf, qimp  , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, qimp, hint

! Local variables

!===============================================================================

! Gradient BCs
coefa = -qimp/hint
coefb = 1.d0

! Flux BCs
cofaf = qimp
cofbf = 0.d0

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     qimpv         Flux value to impose
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_neumann_vector &
 ( coefa , cofaf, coefb , cofbf, qimpv  , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision qimpv(3)
double precision hint

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  coefa(isou) = -qimpv(isou)/hint
  do jsou = 1, 3
    if (jsou.eq.isou) then
      coefb(isou,jsou) = 1.d0
    else
      coefb(isou,jsou) = 0.d0
    endif
  enddo

  ! Flux BCs
  cofaf(isou) = qimpv(isou)
  do jsou = 1, 3
    cofbf(isou,jsou) = 0.d0
  enddo

enddo

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     qimpv         Flux value to impose
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_neumann_vector_ggdh &
 ( coefa , cofaf, coefb , cofbf, qimpv  , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision qimpv(3)
double precision hint(6)

! Local variables

integer          isou  , jsou
double precision invh(6), invdet, m(6)

!===============================================================================
m(1) = hint(2)*hint(3) - hint(5)*hint(5)
m(2) = hint(1)*hint(3) - hint(6)*hint(6)
m(3) = hint(1)*hint(2) - hint(4)*hint(4)
m(4) = hint(5)*hint(6) - hint(4)*hint(3)
m(5) = hint(4)*hint(6) - hint(1)*hint(5)
m(6) = hint(4)*hint(5) - hint(2)*hint(6)

invdet = 1.d0/(hint(1)*m(1) + hint(4)*m(4) + hint(6)*m(6))

invh(1) = m(1) * invdet
invh(2) = m(2) * invdet
invh(3) = m(3) * invdet
invh(4) = m(4) * invdet
invh(5) = m(5) * invdet
invh(6) = m(6) * invdet

! Gradient BCs
coefa(1) = -(invh(1)*qimpv(1) + invh(4)*qimpv(2) + invh(6)*qimpv(3))
coefa(2) = -(invh(4)*qimpv(1) + invh(2)*qimpv(2) + invh(5)*qimpv(3))
coefa(3) = -(invh(6)*qimpv(1) + invh(5)*qimpv(2) + invh(3)*qimpv(3))
coefb(1,1) = invh(1)
coefb(2,2) = invh(2)
coefb(3,3) = invh(3)
coefb(1,2) = invh(4)
coefb(2,1) = invh(4)
coefb(2,3) = invh(5)
coefb(3,2) = invh(5)
coefb(1,3) = invh(6)
coefb(3,1) = invh(6)

do isou = 1, 3

  ! Flux BCs
  cofaf(isou) = qimpv(isou)
  do jsou = 1, 3
    cofbf(isou,jsou) = 0.d0
  enddo

enddo

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefau        explicit BC coefficient for gradients
!> \param[out]    cofafu        explicit BC coefficient for diffusive flux
!> \param[out]    coefav        explicit BC coefficient for gradients
!> \param[out]    cofafv        explicit BC coefficient for diffusive flux
!> \param[out]    coefaw        explicit BC coefficient for gradients
!> \param[out]    cofafw        explicit BC coefficient for diffusive flux
!> \param[out]    coefbu        implicit BC coefficient for gradients
!> \param[out]    cofbfu        implicit BC coefficient for diffusive flux
!> \param[out]    coefbv        implicit BC coefficient for gradients
!> \param[out]    cofbfv        implicit BC coefficient for diffusive flux
!> \param[out]    coefbw        implicit BC coefficient for gradients
!> \param[out]    cofbfw        implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose on the normal
!>                              component
!> \param[in]     qimpv         Flux value to impose on the
!>                              tangential components
!> \param[in]     vect          value of the vector at time n
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     normal        normal
!_______________________________________________________________________________

subroutine set_generalized_sym_scalar &
 ( coefau, cofafu,                                                &
   coefav, cofafv,                                                &
   coefaw, cofafw,                                                &
   coefbu, cofbfu,                                                &
   coefbv, cofbfv,                                                &
   coefbw, cofbfw,                                                &
   pimpv , qimpv , vect  , hint, normal)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefau, coefav, coefaw
double precision cofafu, cofafv, cofafw
double precision coefbu, coefbv, coefbw
double precision cofbfu, cofbfv, cofbfw
double precision hint
double precision normal(3)
double precision vect(3)
double precision pimpv(3), qimpv(3)

! Local variables

integer          isou  , jsou
double precision coefa(3), cofaf(3)
double precision coefb(3), cofbf(3)

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  coefa(isou) = pimpv(isou)*normal(isou)                                 &
              - (1.d0-normal(isou)*normal(isou))*qimpv(isou)/hint
  coefb(isou) = 1.d0 - normal(isou)*normal(isou)

  ! Part which cannot be implicited
  do jsou = 1, 3
    if (jsou.ne.isou) then
      coefa(isou) = coefa(isou) - normal(isou)*normal(jsou)*vect(jsou)
    endif
  enddo

  ! Flux BCs
  cofaf(isou) = -hint*pimpv(isou)*normal(isou)                           &
              + (1.d0-normal(isou)*normal(isou))*qimpv(isou)

  cofbf(isou) = hint*normal(isou)*normal(isou)

  ! Part which cannot be implicited
  do jsou = 1, 3
    if (jsou.ne.isou) then
      cofaf(isou) = cofaf(isou) + hint*normal(isou)*normal(jsou)*vect(jsou)
    endif
  enddo

enddo

coefau = coefa(1)
coefav = coefa(2)
coefaw = coefa(3)

coefbu = coefb(1)
coefbv = coefb(2)
coefbw = coefb(3)

cofafu = cofaf(1)
cofafv = cofaf(2)
cofafw = cofaf(3)

cofbfu = cofbf(1)
cofbfv = cofbf(2)
cofbfw = cofbf(3)

return
end subroutine set_generalized_sym_scalar


!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose on the normal
!>                              component
!> \param[in]     qimpv         Flux value to impose on the
!>                              tangential components
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     normal        normal
!_______________________________________________________________________________

subroutine set_generalized_sym_vector &
 ( coefa , cofaf, coefb , cofbf, pimpv, qimpv, hint, normal)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision hint
double precision normal(3)
double precision pimpv(3), qimpv(3)

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  coefa(isou) = pimpv(isou)*normal(isou)                                 &
              - (1.d0-normal(isou)*normal(isou))*qimpv(isou)/hint
  do jsou = 1, 3
    if (jsou.eq.isou) then
      coefb(isou,jsou) = 1.d0 - normal(isou)*normal(jsou)
    else
      coefb(isou,jsou) = - normal(isou)*normal(jsou)
    endif
  enddo

  ! Flux BCs
  cofaf(isou) = -hint*pimpv(isou)*normal(isou)                           &
              + (1.d0-normal(isou)*normal(isou))*qimpv(isou)
  do jsou = 1, 3
    cofbf(isou,jsou) = hint*normal(isou)*normal(jsou)
  enddo

enddo

return
end subroutine set_generalized_sym_vector

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimp          Flux value to impose
!> \param[in]     cfl           Local Courant number used to convect
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_convective_outlet_scalar &
 ( coefa , cofaf, coefb , cofbf, pimp  , cfl   , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, pimp, cfl, hint

! Local variables

!===============================================================================

! Gradient BCs
coefb = cfl/(1.d0+cfl)
coefa = (1.d0-coefb)*pimp

! Flux BCs
cofaf = -hint*coefa
cofbf =  hint*(1.d0 - coefb)

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose
!> \param[in]     cflv          Local Courant number used to convect
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_convective_outlet_vector &
 ( coefa , cofaf, coefb , cofbf, pimpv  , cflv  , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision pimpv(3), cflv(3)
double precision hint

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  do jsou = 1, 3
    if (jsou.eq.isou) then
      coefb(isou,jsou) = cflv(isou)*(1.d0+cflv(isou))
    else
      coefb(isou,jsou) = 0.d0
    endif
  enddo
  coefa(isou) = (1.d0-coefb(isou,isou))*pimpv(isou)

  ! Flux BCs
  cofaf(isou) = -hint*coefa(isou)
  do jsou = 1, 3
    if (jsou.eq.isou) then
      cofbf(isou,jsou) = hint*(1.d0 - coefb(isou,jsou))
    else
      cofbf(isou,jsou) = 0.d0
    endif
  enddo

enddo

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose
!> \param[in]     cflv          Local Courant number used to convect
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_convective_outlet_vector_ggdh &
 ( coefa , cofaf, coefb , cofbf, pimpv  , cflv  , hint )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision pimpv(3), cflv(3)
double precision hint(6)

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  do jsou = 1, 3
    if (jsou.eq.isou) then
      coefb(isou,jsou) = cflv(isou)*(1.d0+cflv(isou))
    else
      coefb(isou,jsou) = 0.d0
    endif
  enddo
  coefa(isou) = (1.d0-coefb(isou,isou))*pimpv(isou)

enddo

! Flux BCs
cofaf(1) = -(hint(1)*coefa(1) + hint(4)*coefa(2) + hint(6)*coefa(3))
cofaf(2) = -(hint(4)*coefa(1) + hint(2)*coefa(2) + hint(5)*coefa(3))
cofaf(3) = -(hint(6)*coefa(1) + hint(5)*coefa(2) + hint(3)*coefa(3))
cofbf(1,1) = hint(1)*(1.d0 - coefb(1,1))
cofbf(2,2) = hint(2)*(1.d0 - coefb(2,2))
cofbf(3,3) = hint(3)*(1.d0 - coefb(3,3))
cofbf(1,2) = hint(4)*(1.d0 - coefb(1,1))
cofbf(2,1) = hint(4)*(1.d0 - coefb(1,1))
cofbf(2,3) = hint(5)*(1.d0 - coefb(2,2))
cofbf(3,2) = hint(5)*(1.d0 - coefb(2,2))
cofbf(1,3) = hint(6)*(1.d0 - coefb(3,3))
cofbf(3,1) = hint(6)*(1.d0 - coefb(3,3))

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimp          Dirichlet value to impose
!> \param[in]     qimp          Flux value to impose
!_______________________________________________________________________________

subroutine set_dirichlet_conv_neumann_diff_scalar &
 ( coefa, cofaf, coefb, cofbf, pimp, qimp )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, pimp, qimp

! BS test sur hext ? if (abs(hext).gt.rinfin*0.5d0) then

! Gradients BCs
coefa = pimp
coefb = 0.d0

! Flux BCs
cofaf = qimp
cofbf = 0.d0

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose
!> \param[in]     qimpv         Flux value to impose
!_______________________________________________________________________________

subroutine set_dirichlet_conv_neumann_diff_vector &
 ( coefa, cofaf, coefb, cofbf, pimpv, qimpv )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision pimpv(3), qimpv(3)

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! BS test sur hextv ? if (abs(hextv(isou)).gt.rinfin*0.5d0) then

  ! Gradient BCs
  coefa(isou) = pimpv(isou)
  do jsou = 1, 3
    coefb(isou,jsou) = 0.d0
  enddo

  ! Flux BCs
  cofaf(isou) = qimpv(isou)
  do jsou = 1, 3
    cofbf(isou,jsou) = 0.d0
  enddo

enddo

return
end subroutine

!===============================================================================
