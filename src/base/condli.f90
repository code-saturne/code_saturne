!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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
!> \brief Translation of the boundary conditions given by cs_user_boundary_conditions
!> in a form that fits to the solver.
!>
!> The values at a border face \f$ \fib \f$ stored in the face center
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
!> neighbooring cell.
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
!>   which coupled veclocity components next to a boundary. This is only
!>   available when the option ivelco is set to 1.
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
!> \param[in]     isvtb         indicator to save the temperature at
!>                               the walls
!> \param[in,out] icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rought wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[in,out] isostd        indicator for standard outlet
!>                               and reference face index
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughtness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!> \param[out]    coefa         explicit boundary condition coefficient
!> \param[out]    coefb         implicit boundary condition coefficient
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
   isvhb  , isvtb  ,                                              &
   icodcl , isostd ,                                              &
   dt     , rtp    , rtpa   , propce , propfb , rcodcl ,          &
   coefa  , coefb  , visvdr , hbord  , theipb , frcxt  )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal , iterns
integer          isvhb  , isvtb

integer          icodcl(nfabor,nvarcl)
integer          isostd(nfabor+1)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(ndimfb,*)
double precision rcodcl(nfabor,nvarcl,3)
double precision frcxt(3,ncelet)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision visvdr(ncelet)
double precision hbord(nfabor),theipb(nfabor)

! Local variables

integer          ifac  , iel   , ivar
integer          isou  , jsou  , ii
integer          ihcp  , iscal , iscat
integer          inc   , iccocg
integer          iok   , iok1
integer          icodcu
integer          isoent, isorti, ncpt,   isocpt(2)
integer          iclsym, ipatur, ipatrg, isvhbl
integer          ipcvis, ipcvst, ipccp , ipcvsl, ipccv
integer          iclpr , iclu  , iclv  , iclw  , iclk  , iclep
integer          iclnu
integer          icl11 , icl22 , icl33 , icl12 , icl13 , icl23
integer          icl11r, icl22r, icl33r, icl12r, icl13r, icl23r
integer          iclvrr
integer          iclphi, iclfb , iclal , iclomg
integer          iclvar, icluma, iclvma, iclwma
integer          iclprf, icluf , iclvf , iclwf , iclkf , iclepf
integer          iclnuf
integer          icl11f, icl22f, icl33f, icl12f, icl13f, icl23f
integer          iclphf, iclfbf, iclalf, iclomf
integer          iclvaf, iclumf, iclvmf, iclwmf
integer          nswrgp, imligp, iwarnp
integer          itplus, itstar
integer          f_id  ,  iut  , ivt   , iwt, iflmab

double precision sigma , cpp   , rkl
double precision hint  , hext  , pimp  , xdis, qimp, cfl
double precision hintt(6)
double precision flumbf, visclc, visctc, distbf, srfbn2
double precision epsrgp, climgp, extrap
double precision xxp0, xyp0, xzp0
double precision srfbnf, normal(3)
double precision vistot
double precision rinfiv(3), pimpv(3), qimpv(3), hextv(3), cflv(3), vect(3)
double precision visci(3,3), fikis, viscis, distfi
double precision temp

logical          ilved
character*80     fname

double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:,:) :: velipb, rijipb
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:), pointer :: tplusp, tstarp
double precision, dimension(:,:), pointer :: coefaut, cofafut, cofarut
double precision, dimension(:,:,:), pointer :: coefbut, cofbfut, cofbrut
double precision, dimension(:), pointer :: bmasfl

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

! Allocate temporary arrays
allocate(velipb(nfabor,3))

! coefa and coefb are required to compute the cell gradients for the wall
!  turbulent boundary conditions.
! So, their initial values are kept (note that at the first time step, they are
!  initialized to zero flux in inivar.f90)

! velipb stores the velocity in I' of boundary cells

! Initialize variables to avoid compiler warnings

iclep = 0
iclk = 0
iclomg = 0
iclfb = 0
iclphi = 0
iclnu = 0
icl11 = 0
icl22 = 0
icl33 = 0
icl12 = 0
icl13 = 0
icl23 = 0
icl11r = 0
icl12r = 0
icl13r = 0
icl22r = 0
icl23r = 0
icl33r = 0
icl11f = 0
icl12f = 0
icl13f = 0
icl22f = 0
icl23f = 0
icl33f = 0
iclal = 0
iclalf= 0
iclvar = 0
iclvaf = 0
icluf = 0
iclvf = 0
iclwf = 0
ipccv = 0
iclepf = 0
iclfbf = 0
iclkf = 0
iclnuf = 0
iclomf = 0
iclphf = 0
iclvrr = 0

rinfiv(1) = rinfin
rinfiv(2) = rinfin
rinfiv(3) = rinfin

! pointers to T+ and T* if saved

tplusp => null()
tstarp => null()

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

! Pointers to the mass fluxes
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

!===============================================================================
! 2. Treatment of types of BCs given by itypfb
!===============================================================================

if (ippmod(iphpar).ge.1) then
  call pptycl &
  !==========
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  , rcodcl )
endif

if (iale.eq.1) then
  call altycl &
  !==========
 ( itypfb , ialtyb , icodcl , impale ,                            &
   dt     ,                                                       &
   rcodcl , xyzno0 , depale )
endif

if (imobil.eq.1) then
  call mmtycl(itypfb, rcodcl)
  !==========
endif

call typecl &
!==========
 ( nvar   , nscal  ,                                              &
   itypfb , itrifb , icodcl , isostd ,                            &
   rtpa   , propce , rcodcl , frcxt  )

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

! --- Gradient Boundary Conditions
iclpr = iclrtp(ipr,icoef)
iclu  = iclrtp(iu, icoef)
iclv  = iclrtp(iv, icoef)
iclw  = iclrtp(iw, icoef)
if (itytur.eq.2) then
  iclk  = iclrtp(ik ,icoef)
  iclep = iclrtp(iep,icoef)
elseif (itytur.eq.3) then
  icl11 = iclrtp(ir11,icoef)
  icl22 = iclrtp(ir22,icoef)
  icl33 = iclrtp(ir33,icoef)
  icl12 = iclrtp(ir12,icoef)
  icl13 = iclrtp(ir13,icoef)
  icl23 = iclrtp(ir23,icoef)
  ! Boundary conditions for the momentum equation
  icl11r = iclrtp(ir11,icoefr)
  icl22r = iclrtp(ir22,icoefr)
  icl33r = iclrtp(ir33,icoefr)
  icl12r = iclrtp(ir12,icoefr)
  icl13r = iclrtp(ir13,icoefr)
  icl23r = iclrtp(ir23,icoefr)
  iclep = iclrtp(iep,icoef)
  if (iturb.eq.32) iclal = iclrtp(ial,icoef)
elseif (itytur.eq.5) then
  iclk   = iclrtp(ik ,icoef)
  iclep  = iclrtp(iep,icoef)
  iclphi = iclrtp(iphi,icoef)
  if (iturb.eq.50) then
    iclfb = iclrtp(ifb,icoef)
  elseif (iturb.eq.51) then
    iclal = iclrtp(ial,icoef)
  endif
elseif (iturb.eq.60) then
  iclk   = iclrtp(ik ,icoef)
  iclomg = iclrtp(iomg,icoef)
elseif (iturb.eq.70) then
  iclnu = iclrtp(inusa,icoef)
endif

! --- Flux Boundary Conditions
iclprf = iclrtp(ipr,icoeff)
icluf  = iclrtp(iu, icoeff)
iclvf  = iclrtp(iv, icoeff)
iclwf  = iclrtp(iw, icoeff)
if (itytur.eq.2) then
  iclkf  = iclrtp(ik, icoeff)
  iclepf = iclrtp(iep,icoeff)
elseif (itytur.eq.3) then
  icl11f = iclrtp(ir11,icoeff)
  icl22f = iclrtp(ir22,icoeff)
  icl33f = iclrtp(ir33,icoeff)
  icl12f = iclrtp(ir12,icoeff)
  icl13f = iclrtp(ir13,icoeff)
  icl23f = iclrtp(ir23,icoeff)
  iclepf = iclrtp(iep,icoeff)
  if (iturb.eq.32) iclalf = iclrtp(ial,icoeff)
elseif (itytur.eq.5) then
  iclkf  = iclrtp(ik ,icoeff)
  iclepf = iclrtp(iep,icoeff)
  iclphf = iclrtp(iphi,icoeff)
  if (iturb.eq.50) then
    iclfbf = iclrtp(ifb,icoeff)
  elseif (iturb.eq.51) then
    iclalf = iclrtp(ial,icoeff)
  endif
elseif (iturb.eq.60) then
  iclkf  = iclrtp(ik, icoeff)
  iclomf = iclrtp(iomg,icoeff)
elseif (iturb.eq.70) then
  iclnuf = iclrtp(inusa,icoeff)
endif

! --- Physical quantities
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
if (icp.gt.0) then
  ipccp  = ipproc(icp)
else
  ipccp = 0
endif

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
allocate(grad(ncelet,3))

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
!  Dans le cas du couplage avec le module 1D, on utilise le scalaire
!    couple avec Syrthes s'il y a couplage, sinon iscalt.
!  La valeur de isvtb a ete initialisee dans tridim
!    au numero du scalaire couple.

!  Pour le rayonnement
!  -------------------
!  On calcule la valeur en I' s'il y a une variable
!    thermique


!  On recherche l'unique scalaire qui convient
!     (ce peut etre T, H, ou E (en compressible))

iscat = 0

! Si un scalaire est couple a SYRTHES ou au module 1D
if (isvtb.ne.0) then
  ! si ce n'est pas la variable thermique, ca ne va pas.
  if(isvtb.ne.iscalt) then
    write(nfecra,8000)isvtb,iscalt
    call csexit (1)
    !==========
    !         sinon, on calcule le gradient.
  else
    iscat = isvtb
  endif
else
  iscat = iscalt
endif

! Compute the boundary value of the thermal scalar in I' if required
if (iscat.gt.0) then

  ivar = isca(iscat)

  if (ntcabs.gt.1 .and. itbrrb.eq.1 .and. ircflu(ivar).eq.1) then

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
    iclvar = iclrtp(ivar,icoef)

    call grdcel &
    !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtpa(1,ivar)    , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   grad   )

    do ifac = 1 , nfabor
      iel = ifabor(ifac)
      theipb(ifac) = rtpa(iel,ivar) &
                   + grad(iel,1)*diipb(1,ifac) &
                   + grad(iel,2)*diipb(2,ifac) &
                   + grad(iel,3)*diipb(3,ifac)
    enddo

  else

    do ifac = 1 , nfabor
      iel = ifabor(ifac)
      theipb(ifac) = rtpa(iel,ivar)
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

if (iclsym.ne.0.or.ipatur.ne.0.or.ipatrg.ne.0) then

  if(ntcabs.gt.1) then

    ! Allocate a temporary array
    allocate(gradv(ncelet,3,3))

    iccocg = 1
    inc    = 1
    nswrgp = nswrgr(iu)
    imligp = imligr(iu)
    iwarnp = iwarni(iu)
    epsrgp = epsrgr(iu)
    climgp = climgr(iu)
    extrap = extrag(iu)
    iclvar = iclrtp(iu,icoef)

    if (ivelco.eq.1) then

      ilved = .false.

      call grdvec &
      !==========
    ( iu     , imrgra , inc    , nswrgp , imligp ,                   &
      iwarnp , nfecra ,                                              &
      epsrgp , climgp , extrap ,                                     &
      ilved ,                                                        &
      rtpa(1,iu) ,  coefau , coefbu,                                 &
      gradv  )

    else

      call grdvni &
      !==========
    ( iu     , imrgra , inc    , iccocg , nswrgp , imligp ,          &
      iwarnp , nfecra ,                                              &
      epsrgp , climgp , extrap ,                                     &
      rtpa(1,iu)      , coefa(1,iclvar) , coefb(1,iclvar) ,          &
      gradv  )

    endif

    do isou = 1, 3
      if(isou.eq.1) ivar = iu
      if(isou.eq.2) ivar = iv
      if(isou.eq.3) ivar = iw

      do ifac = 1, nfabor
        iel = ifabor(ifac)
        velipb(ifac,isou) = gradv(iel,1,isou)*diipb(1,ifac)    &
                          + gradv(iel,2,isou)*diipb(2,ifac)    &
                          + gradv(iel,3,isou)*diipb(3,ifac)    &
                          + rtpa(iel,ivar)
      enddo
    enddo

    deallocate(gradv)

  ! Nb: at the first time step, coefa and coefb are unknown, so the walue
  !     in I is stored instead of the value in I'
  else

    do isou = 1, 3
      if(isou.eq.1) ivar = iu
      if(isou.eq.2) ivar = iv
      if(isou.eq.3) ivar = iw

      do ifac = 1, nfabor
        iel = ifabor(ifac)
        velipb(ifac,isou) = rtpa(iel,ivar)
      enddo

    enddo

  endif

endif

! ---> Compute Rij in I' for boundary cells

if ((iclsym.ne.0.or.ipatur.ne.0.or.ipatrg.ne.0).and.itytur.eq.3) then

  ! Allocate a work array to store Rij values at boundary faces
  allocate(rijipb(nfabor,6))

  do isou = 1 , 6

    if(isou.eq.1) ivar = ir11
    if(isou.eq.2) ivar = ir22
    if(isou.eq.3) ivar = ir33
    if(isou.eq.4) ivar = ir12
    if(isou.eq.5) ivar = ir13
    if(isou.eq.6) ivar = ir23


    if(ntcabs.gt.1.and.irijrb.eq.1) then

      inc = 1
      iccocg = 1
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)
      iclvar = iclrtp(ivar,icoef)

      call grdcel &
      !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtpa(1,ivar)    , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   grad   )

      do ifac = 1 , nfabor
        iel = ifabor(ifac)
        rijipb(ifac,isou) = rtpa(iel,ivar)            &
                          + grad(iel,1)*diipb(1,ifac) &
                          + grad(iel,2)*diipb(2,ifac) &
                          + grad(iel,3)*diipb(3,ifac)
      enddo

    ! Nb: at the first time step, coefa and coefb are unknown, so the walue
    !     in I is stored instead of the value in I'
    else

      do ifac = 1 , nfabor
        iel = ifabor(ifac)
        rijipb(ifac,isou) = rtpa(iel,ivar)
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

! Initialization of the array storing yplus
!  which is computed in clptur.f90 and/or clptrg.f90

if (ipstdv(ipstyp).ne.0) then
  do ifac = 1, nfabor
    yplbr(ifac) = 0.d0
  enddo
endif

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
   rtp    , rtpa   , propce , propfb , rcodcl ,                   &
   velipb , rijipb , coefa  , coefb  , visvdr ,                   &
   hbord  , theipb )

endif

if (ipatrg.ne.0) then

  ! Rough wall laws
  call clptrg &
  !==========
 ( nscal  , isvhb  ,                                              &
   icodcl ,                                                       &
   dt     , rtp    , rtpa   , propce , propfb , rcodcl ,          &
   velipb , rijipb , coefa  , coefb  , visvdr ,                   &
   hbord  , theipb )

endif

!===============================================================================
! 7. Symmetry for vectors and tensors
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
   rtp    , rtpa   , propce , rcodcl ,                            &
   velipb , rijipb , coefa  , coefb  )

endif

!===============================================================================
! 8. Velocity: Outlet, Dirichlet and Neumann and convectiv outlet
!===============================================================================

! ---> Outlet: in case of incomming mass flux, the mass flux is set to zero.

isoent = 0
isorti = 0
do ifac = 1, nfabor

  if (icodcl(ifac,iu).eq.9) then

    flumbf = bmasfl(ifac)

    ! --- Physical Properties
    iel = ifabor(ifac)
    visclc = propce(iel,ipcvis)
    visctc = propce(iel,ipcvst)

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

      pimp = 0.d0

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclu), coefa(ifac,icluf),             &
           coefb(ifac,iclu), coefb(ifac,icluf),             &
           pimp            , hint             , rinfin )

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclv), coefa(ifac,iclvf),             &
           coefb(ifac,iclv), coefb(ifac,iclvf),             &
           pimp            , hint             , rinfin )

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclw), coefa(ifac,iclwf),             &
           coefb(ifac,iclw), coefb(ifac,iclwf),             &
           pimp            , hint             , rinfin )


      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then

        pimpv(1) = 0.d0
        pimpv(2) = 0.d0
        pimpv(3) = 0.d0

        call set_dirichlet_vector &
             !====================
           ( coefau(1,ifac)  , cofafu(1,ifac)  ,             &
             coefbu(1,1,ifac), cofbfu(1,1,ifac),             &
             pimpv           , hint            , rinfiv )

      endif

      isoent = isoent + 1

    else

      ! Neumann Boundary Conditions
      !----------------------------

      qimp = 0.d0

      call set_neumann_scalar &
           !==================
         ( coefa(ifac,iclu), coefa(ifac,icluf),             &
           coefb(ifac,iclu), coefb(ifac,icluf),             &
           qimp            , hint )

      call set_neumann_scalar &
           !==================
         ( coefa(ifac,iclv), coefa(ifac,iclvf),             &
           coefb(ifac,iclv), coefb(ifac,iclvf),             &
           qimp            , hint )

      call set_neumann_scalar &
           !==================
         ( coefa(ifac,iclw), coefa(ifac,iclwf),             &
           coefb(ifac,iclw), coefb(ifac,iclwf),             &
           qimp            , hint )


      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then

        qimpv(1) = 0.d0
        qimpv(2) = 0.d0
        qimpv(3) = 0.d0

        call set_neumann_vector &
             !==================
           ( coefau(1,ifac)  , cofafu(1,ifac)  ,             &
             coefbu(1,1,ifac), cofbfu(1,1,ifac),             &
             qimpv           , hint )

      endif

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
  visclc = propce(iel,ipcvis)
  visctc = propce(iel,ipcvst)

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

    pimp = rcodcl(ifac,iu,1)
    hext = rcodcl(ifac,iu,2)

    call set_dirichlet_scalar &
         !====================
       ( coefa(ifac,iclu), coefa(ifac,icluf),             &
         coefb(ifac,iclu), coefb(ifac,icluf),             &
         pimp            , hint             , hext )

    pimp = rcodcl(ifac,iv,1)
    hext = rcodcl(ifac,iv,2)

    call set_dirichlet_scalar &
         !====================
       ( coefa(ifac,iclv), coefa(ifac,iclvf),             &
         coefb(ifac,iclv), coefb(ifac,iclvf),             &
         pimp            , hint             , hext )

    pimp = rcodcl(ifac,iw,1)
    hext = rcodcl(ifac,iw,2)

    call set_dirichlet_scalar &
         !====================
       ( coefa(ifac,iclw), coefa(ifac,iclwf),             &
         coefb(ifac,iclw), coefb(ifac,iclwf),             &
         pimp            , hint             , hext )


    ! Coupled solving of the velocity components
    if (ivelco.eq.1) then

      pimpv(1) = rcodcl(ifac,iu,1)
      pimpv(2) = rcodcl(ifac,iv,1)
      pimpv(3) = rcodcl(ifac,iw,1)
      hextv(1) = rcodcl(ifac,iu,2)
      hextv(2) = rcodcl(ifac,iv,2)
      hextv(3) = rcodcl(ifac,iw,2)

      call set_dirichlet_vector &
           !====================
         ( coefau(1,ifac)  , cofafu(1,ifac)  ,             &
           coefbu(1,1,ifac), cofbfu(1,1,ifac),             &
           pimpv           , hint            , hextv )

    endif

  ! Neumann Boundary Conditions
  !----------------------------

  elseif (icodcl(ifac,iu).eq.3) then

    qimp = rcodcl(ifac,iu,3)

    call set_neumann_scalar &
         !==================
       ( coefa(ifac,iclu), coefa(ifac,icluf),             &
         coefb(ifac,iclu), coefb(ifac,icluf),             &
         qimp            , hint )

    qimp = rcodcl(ifac,iv,3)

    call set_neumann_scalar &
         !==================
       ( coefa(ifac,iclv), coefa(ifac,iclvf),             &
         coefb(ifac,iclv), coefb(ifac,iclvf),             &
         qimp            , hint )

    qimp = rcodcl(ifac,iw,3)

    call set_neumann_scalar &
         !==================
       ( coefa(ifac,iclw), coefa(ifac,iclwf),             &
         coefb(ifac,iclw), coefb(ifac,iclwf),             &
         qimp            , hint )


    ! Coupled solving of the velocity components
    if (ivelco.eq.1) then

      qimpv(1) = rcodcl(ifac,iu,3)
      qimpv(2) = rcodcl(ifac,iv,3)
      qimpv(3) = rcodcl(ifac,iw,3)

      call set_neumann_vector &
           !==================
         ( coefau(1,ifac)  , cofafu(1,ifac)  ,             &
           coefbu(1,1,ifac), cofbfu(1,1,ifac),             &
           qimpv           , hint )

    endif

  ! Convective Boundary Conditions
  !-------------------------------

  elseif (icodcl(ifac,iu).eq.2) then

    pimp = rcodcl(ifac,iu,1)
    cfl = rcodcl(ifac,iu,2)

    call set_convective_outlet_scalar &
         !==================
       ( coefa(ifac,iclu), coefa(ifac,icluf),             &
         coefb(ifac,iclu), coefb(ifac,icluf),             &
         pimp            , cfl              , hint )

    pimp = rcodcl(ifac,iv,1)
    cfl = rcodcl(ifac,iv,2)

    call set_convective_outlet_scalar &
         !==================
       ( coefa(ifac,iclv), coefa(ifac,iclvf),             &
         coefb(ifac,iclv), coefb(ifac,iclvf),             &
         pimp            , cfl              , hint )

    pimp = rcodcl(ifac,iw,1)
    cfl = rcodcl(ifac,iw,2)

    call set_convective_outlet_scalar &
         !==================
       ( coefa(ifac,iclw), coefa(ifac,iclwf),             &
         coefb(ifac,iclw), coefb(ifac,iclwf),             &
         pimp            , cfl              , hint )


    ! Coupled solving of the velocity components
    if (ivelco.eq.1) then

      pimpv(1) = rcodcl(ifac,iu,1)
      cflv(1) = rcodcl(ifac,iu,2)
      pimpv(2) = rcodcl(ifac,iv,1)
      cflv(2) = rcodcl(ifac,iv,2)
      pimpv(3) = rcodcl(ifac,iw,1)
      cflv(3) = rcodcl(ifac,iw,2)

      call set_convective_outlet_vector &
           !==================
         ( coefau(1,ifac)  , cofafu(1,ifac)  ,             &
           coefbu(1,1,ifac), cofbfu(1,1,ifac),             &
           pimpv           , cflv            , hint )

    endif

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
    if (ivelco.eq.1) then

      call set_generalized_sym_vector &
           !=========================
         ( coefau(1,ifac)  , cofafu(1,ifac)  ,             &
           coefbu(1,1,ifac), cofbfu(1,1,ifac),             &
           pimpv           , qimpv            , hint , normal )

    else

      vect(1) = velipb(ifac, 1)
      vect(2) = velipb(ifac, 2)
      vect(3) = velipb(ifac, 3)

       call set_generalized_sym_scalar &
           !==========================
         ( coefa(ifac,iclu), coefa(ifac,icluf),                       &
           coefa(ifac,iclv), coefa(ifac,iclvf),                       &
           coefa(ifac,iclw), coefa(ifac,iclwf),                       &
           coefb(ifac,iclu), coefb(ifac,icluf),                       &
           coefb(ifac,iclv), coefb(ifac,iclvf),                       &
           coefb(ifac,iclw), coefb(ifac,iclwf),                       &
           pimpv           , qimpv            , vect , hint , normal )

    endif


  endif

enddo

!===============================================================================
! 9. Pressure: Dirichlet and Neumann and convectiv outlet
!===============================================================================

iclvar = iclpr
iclvaf = iclprf

do ifac = 1, nfabor

  iel = ifabor(ifac)

  ! --- Geometrical quantities
  distbf = distb(ifac)

  ! If a flux dt.grad P (W/m2) is set in cs_user_boundary
  if (idften(ipr).eq.1) then
    hint = dt(iel)/distbf
  else if (idften(ipr).eq.3) then
    hint = ( dttens(1, iel)*surfbo(1,ifac)**2              &
           + dttens(2, iel)*surfbo(2,ifac)**2              &
           + dttens(3, iel)*surfbo(3,ifac)**2              &
           ) / (surfbn(ifac)**2 * distbf)
  ! Symmetric tensor diffusivity
  elseif (idften(ipr).eq.6) then

    visci(1,1) = dttens(1,iel)
    visci(2,2) = dttens(2,iel)
    visci(3,3) = dttens(3,iel)
    visci(1,2) = dttens(4,iel)
    visci(2,1) = dttens(4,iel)
    visci(2,3) = dttens(5,iel)
    visci(3,2) = dttens(5,iel)
    visci(1,3) = dttens(6,iel)
    visci(3,1) = dttens(6,iel)

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

  ! On doit remodifier la valeur du  Dirichlet de pression de manière
  !  à retrouver P*. Car dans typecl.f90 on a travaillé avec la pression
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
       ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
         coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
         pimp              , hint             , hext )

  endif

  ! Neumann Boundary Conditions
  !----------------------------

  if (icodcl(ifac,ipr).eq.3) then

    qimp = rcodcl(ifac,ipr,3)

    call set_neumann_scalar &
         !=================
       ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
         coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
         qimp              , hint )

  ! Convective Boundary Conditions
  !-------------------------------

  elseif (icodcl(ifac,ipr).eq.2) then

    pimp = rcodcl(ifac,ipr,1)
    cfl = rcodcl(ifac,ipr,2)

    call set_convective_outlet_scalar &
         !===========================
       ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
         coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
         pimp              , cfl               , hint )

  elseif (icodcl(ifac,ipr).eq.13) then

    pimp = rcodcl(ifac,ipr,1)
    qimp = rcodcl(ifac,ipr,3)

    call set_dirichlet_conv_neumann_diff_scalar &
         !===========================
       ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
         coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
         pimp              , qimp )

  endif

enddo


!===============================================================================
! 10. Turbulent quantities: Dirichlet and Neumann and convectiv outlet
!===============================================================================

! ---> k-epsilon and k-omega

if (itytur.eq.2.or.iturb.eq.60) then

  do ii = 1, 2

    !     Pour le k-omega, on met les valeurs sigma_k2 et sigma_w2 car ce terme
    !     ne concerne en pratique que les entrees (pas de pb en paroi ou en flux
    !     nul)
    if (ii.eq.1.and.itytur.eq.2) then
      ivar   = ik
      iclvar = iclk
      iclvaf = iclkf
      sigma  = sigmak
    elseif (ii.eq.1.and.iturb.eq.60) then
      ivar   = ik
      iclvar = iclk
      iclvaf = iclkf
      sigma  = ckwsk2 !FIXME it is not consistent with the model
    elseif (itytur.eq.2) then
      ivar   = iep
      iclvar = iclep
      iclvaf = iclepf
      sigma  = sigmae
    else
      ivar   = iomg
      iclvar = iclomg
      iclvaf = iclomf
      sigma  = ckwsw2 !FIXME it is not consistent with the model
    endif

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Physical Propreties
      visclc = propce(iel,ipcvis)
      visctc = propce(iel,ipcvst)

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
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , hint              , hext )

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , cfl               , hint )

      endif

    enddo

  enddo

! ---> Rij-epsilon
elseif (itytur.eq.3) then

  ! --> Rij

  do isou = 1, 6

    if(isou.eq.1) then
      ivar   = ir11
      iclvar = icl11
      iclvrr = icl11r
      iclvaf = icl11f
    elseif(isou.eq.2) then
      ivar   = ir22
      iclvar = icl22
      iclvrr = icl22r
      iclvaf = icl22f
    elseif(isou.eq.3) then
      ivar   = ir33
      iclvar = icl33
      iclvrr = icl33r
      iclvaf = icl22f
    elseif(isou.eq.4) then
      ivar   = ir12
      iclvar = icl12
      iclvrr = icl12r
      iclvaf = icl12f
    elseif(isou.eq.5) then
      ivar   = ir13
      iclvar = icl13
      iclvrr = icl13r
      iclvaf = icl13f
    elseif(isou.eq.6) then
      ivar   = ir23
      iclvar = icl23
      iclvrr = icl23r
      iclvaf = icl23f
    endif

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Physical Propreties
      visclc = propce(iel,ipcvis)

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

      else
        call csexit(1)
      endif

      ! Dirichlet Boundary Condition
      !-----------------------------

      if (icodcl(ifac,ivar).eq.1) then
        pimp = rcodcl(ifac,ivar,1)
        hext = rcodcl(ifac,ivar,2)

        call set_dirichlet_scalar &
             !====================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , hint              , hext )

        ! Boundary conditions for the momentum equation
        coefa(ifac,iclvrr) = coefa(ifac,iclvar)
        coefb(ifac,iclvrr) = coefb(ifac,iclvar)
      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then
        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             qimp              , hint )

        ! Boundary conditions for the momentum equation
        coefa(ifac,iclvrr) = coefa(ifac,iclvar)
        coefb(ifac,iclvrr) = coefb(ifac,iclvar)
      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then
        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , cfl               , hint )
        ! Boundary conditions for the momentum equation
        coefa(ifac,iclvrr) = coefa(ifac,iclvar)
        coefb(ifac,iclvrr) = coefb(ifac,iclvar)

      endif

    enddo

  enddo

  ! --> epsilon

  ivar   = iep
  iclvar = iclep
  iclvaf = iclepf

  do ifac = 1, nfabor

    iel = ifabor(ifac)

    ! --- Physical Propreties
    visclc = propce(iel,ipcvis)
    visctc = propce(iel,ipcvst)

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

    else
      call csexit(1)
    endif

    ! Dirichlet Boundary Condition
    !-----------------------------

    if (icodcl(ifac,ivar).eq.1) then

      pimp = rcodcl(ifac,ivar,1)
      hext = rcodcl(ifac,ivar,2)

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , hext )

    endif

    ! Neumann Boundary Conditions
    !----------------------------

    if (icodcl(ifac,ivar).eq.3) then

      qimp = rcodcl(ifac,ivar,3)

      call set_neumann_scalar &
           !==================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , cfl               , hint )

    endif

  enddo

  ! --> alpha for the EBRSM

  if (iturb.eq.32) then
    ivar   = ial
    iclvar = iclal
    iclvaf = iclalf


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
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , hint              , hext )

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , cfl               , hint )

      endif

    enddo
  endif

! ---> v2f type models (phi_bar and Bl-v2/k)

elseif (itytur.eq.5) then

  !   --> k, epsilon  and phi
  do ii = 1, 3

    if (ii.eq.1) then
      ivar   = ik
      iclvar = iclk
      iclvaf = iclkf
      sigma  = sigmak
    elseif (ii.eq.2) then
      ivar   = iep
      iclvar = iclep
      iclvaf = iclepf
      sigma  = sigmae
    else
      ivar   = iphi
      iclvar = iclphi
      iclvaf = iclphf
      sigma  = sigmak
    endif

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Physical Propreties
      visclc = propce(iel,ipcvis)
      visctc = propce(iel,ipcvst)

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
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , hint              , hext )

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , cfl               , hint )

      endif

    enddo

  enddo

  if (iturb.eq.50) then

    ! --> FB

    ivar   = ifb
    iclvar = iclfb
    iclvaf = iclfbf

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
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , hint              , hext )

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , cfl               , hint )

      endif

    enddo

  elseif (iturb.eq.51) then

    ! --> alpha

    ivar   = ial
    iclvar = iclal
    iclvaf = iclalf

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
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , hint              , hext )

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             qimp              , hint )

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , cfl               , hint )

      endif

    enddo

  endif

! ---> Spalart Allmaras

elseif (iturb.eq.70) then

  ivar   = inusa
  iclvar = iclnu
  iclvaf = iclnuf

  do ifac = 1, nfabor

    iel = ifabor(ifac)

    ! --- Physical Propreties
    visclc = propce(iel,ipcvis)

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
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , hext )

    endif

    ! Neumann Boundary Conditions
    !----------------------------

    if (icodcl(ifac,ivar).eq.3) then

      qimp = rcodcl(ifac,ivar,3)

      call set_neumann_scalar &
           !==================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           qimp              , hint )

    ! Convective Boundary Conditions
    !-------------------------------

    elseif (icodcl(ifac,ivar).eq.2) then

      pimp = rcodcl(ifac,ivar,1)
      cfl = rcodcl(ifac,ivar,2)

      call set_convective_outlet_scalar &
           !==================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , cfl               , hint )

    endif

  enddo

endif

!===============================================================================
! 11. Other scalars (except variances):
!     Dirichlet and Neumann and convectiv outlet
!===============================================================================

if (nscal.ge.1) then

  do ii = 1, nscal

    ivar   = isca(ii)
    iclvar = iclrtp(ivar,icoef)
    iclvaf = iclrtp(ivar,icoeff)

    isvhbl = 0
    if(ii.eq.isvhb) then
      isvhbl = isvhb
    endif

    if(ivisls(ii).gt.0) then
      ipcvsl = ipproc(ivisls(ii))
    else
      ipcvsl = 0
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

    if (iscsth(iscal).eq.0.or.             &
        iscsth(iscal).eq.2.or.             &
        iscsth(iscal).eq.3) then
      ihcp = 0
    elseif(abs(iscsth(iscal)).eq.1) then
      if(ipccp.gt.0) then
        ihcp = 2
      else
        ihcp = 1
      endif
    endif

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Physical Properties
      visctc = propce(iel,ipcvst)

      ! --- Geometrical quantities
      distbf = distb(ifac)

      ! --- Prise en compte de Cp ou CV
      !      (dans le Cas compressible ihcp=0)

      cpp = 1.d0
      if (ihcp.eq.0) then
        cpp = 1.d0
      elseif (ihcp.eq.2) then
        cpp = propce(iel,ipccp)
      elseif (ihcp.eq.1) then
        cpp = cp0
      endif

      ! --- Viscosite variable ou non
      if (ipcvsl.le.0) then
        rkl = visls0(ii)
      else
        rkl = propce(iel,ipcvsl)
      endif

      ! Scalar diffusivity
      if (idften(ivar).eq.1) then
        hint = (rkl+idifft(ivar)*cpp*visctc/sigmas(ii))/distbf

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
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
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
        !            IF(IPCCP.GT.0) THEN
        !              CPR = PROPCE(IEL,IPCCP )
        !            ELSE
        !              CPR = CP0
        !            ENDIF
        !         puisque HINT = VISLS / DISTBR
        !                      = lambda/(CP * distance)

        !       lorsque la variable transportee est l'energie
        !         ISCSTH(II).EQ.3 :
        !         on procede comme pour l'enthalpie avec CV au lieu de CP
        !         (rq : il n'y a pas d'hypothèse, sf en non orthogonal :
        !               le flux est le bon et le coef d'echange aussi)

        !      De meme plus bas et de meme dans clptur.

        !               Si on rayonne et que
        !                  le scalaire est la variable energetique

        if (iirayo.ge.1 .and. ii.eq.iscalt) then

          ! We compute the exchange coefficient in W/(m2 K)

          ! Enthalpy
          if (iscsth(ii).eq.2) then
            ! If Cp is variable
            if (ipccp.gt.0) then
              propfb(ifac,ipprob(ihconv)) = hint*propce(iel,ipccp )
            else
              propfb(ifac,ipprob(ihconv)) = hint*cp0
            endif

          ! Energie (compressible module)
          elseif (iscsth(ii).eq.3) then
            ! If Cv is variable
            if (ipccv.gt.0) then
              propfb(ifac,ipprob(ihconv)) = hint*propce(iel,ipccv)
            else
              propfb(ifac,ipprob(ihconv)) = hint*cv0
            endif

          ! Temperature
          elseif(abs(iscsth(ii)).eq.1) then
            propfb(ifac,ipprob(ihconv)) = hint
          endif

          ! The outgoing flux is stored (Q = h(Ti'-Tp): negative if
          !  gain for the fluid) in W/m2
          propfb(ifac,ipprob(ifconv)) = coefa(ifac,iclvaf)              &
                                      + coefb(ifac,iclvaf)*theipb(ifac)

        endif

      endif

      ! Neumann Boundary Conditions
      !----------------------------

      if (icodcl(ifac,ivar).eq.3) then

        qimp = rcodcl(ifac,ivar,3)

        call set_neumann_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             qimp              , hint )

        if (isvhbl.gt.0) hbord(ifac) = hint

        !--> Rayonnement :

        if (iirayo.ge.1 .and. ii.eq.iscalt) then

          ! We compute the exchange coefficient in W/(m2 K)

          ! Enthalpy
          if (iscsth(ii).eq.2) then
            ! If Cp is variable
            if (ipccp.gt.0) then
              propfb(ifac,ipprob(ihconv)) = hint*propce(iel, ipccp)
            else
              propfb(ifac,ipprob(ihconv)) = hint*cp0
            endif

          ! Energie (compressible module)
          elseif (iscsth(ii).eq.3) then
            ! If Cv is variable
            if (ipccv.gt.0) then
              propfb(ifac,ipprob(ihconv)) = hint*propce(iel,ipccv)
            else
              propfb(ifac,ipprob(ihconv)) = hint*cv0
            endif

          ! Temperature
          elseif (abs(iscsth(ii)).eq.1) then
            propfb(ifac,ipprob(ihconv)) = hint
          endif

          ! The outgoing flux is stored (Q = h(Ti'-Tp): negative if
          !  gain for the fluid) in W/m2
          propfb(ifac,ipprob(ifconv)) = rcodcl(ifac,ivar,3)
        endif

      ! Convective Boundary Conditions
      !-------------------------------

      elseif (icodcl(ifac,ivar).eq.2) then

        pimp = rcodcl(ifac,ivar,1)
        cfl = rcodcl(ifac,ivar,2)

        call set_convective_outlet_scalar &
             !==================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , cfl               , hint )

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
        visclc = propce(iel,ipcvis)

        ! --- Geometrical quantities
        distbf = distb(ifac)

        if (ivisls(iscal).le.0) then
          rkl = visls0(iscal)/cpp
        else
          rkl = propce(iel,ipproc(ivisls(iscal)))/cpp
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
! 13. Mesh velocity (ALE module): Dirichlet and Neumann and convectiv outlet
!===============================================================================

if (iale.eq.1) then

  icluma = iclrtp(iuma,icoef)
  iclvma = iclrtp(ivma,icoef)
  iclwma = iclrtp(iwma,icoef)
  iclumf = iclrtp(iuma,icoeff)
  iclvmf = iclrtp(ivma,icoeff)
  iclwmf = iclrtp(iwma,icoeff)

  do ifac = 1, nfabor

    iel = ifabor(ifac)
    distbf = distb(ifac)
    srfbn2 = surfbn(ifac)**2
    if (iortvm.eq.0) then
      hint = propce(iel,ipproc(ivisma(1)))/distbf
    else
      hint = ( propce(iel,ipproc(ivisma(1)))*surfbo(1,ifac)**2    &
             + propce(iel,ipproc(ivisma(2)))*surfbo(2,ifac)**2    &
             + propce(iel,ipproc(ivisma(3)))*surfbo(3,ifac)**2 )  &
           /distbf/srfbn2
    endif

    ! Dirichlet Boundary Conditions
    !------------------------------

    if (icodcl(ifac,iuma).eq.1) then

      pimp = rcodcl(ifac,iuma,1)
      hext = rcodcl(ifac,iuma,2)

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,icluma), coefa(ifac,iclumf),             &
           coefb(ifac,icluma), coefb(ifac,iclumf),             &
           pimp              , hint              , hext )

      pimp = rcodcl(ifac,ivma,1)
      hext = rcodcl(ifac,ivma,2)

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvma), coefa(ifac,iclvmf),             &
           coefb(ifac,iclvma), coefb(ifac,iclvmf),             &
           pimp              , hint              , hext )

      pimp = rcodcl(ifac,iwma,1)
      hext = rcodcl(ifac,iwma,2)

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclwma), coefa(ifac,iclwmf),             &
           coefb(ifac,iclwma), coefb(ifac,iclwmf),             &
           pimp              , hint              , hext )


      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then

        pimpv(1) = rcodcl(ifac,iuma,1)
        pimpv(2) = rcodcl(ifac,ivma,1)
        pimpv(3) = rcodcl(ifac,iwma,1)
        hextv(1) = rcodcl(ifac,iuma,2)
        hextv(2) = rcodcl(ifac,ivma,2)
        hextv(3) = rcodcl(ifac,iwma,2)

        call set_dirichlet_vector &
             !====================
           ( claale(1,ifac)  , cfaale(1,ifac)  ,             &
             clbale(1,1,ifac), cfbale(1,1,ifac),             &
             pimpv           , hint            , hextv )

      endif

    ! Neumann Boundary Conditions
    !----------------------------

    elseif (icodcl(ifac,iuma).eq.3) then

      qimp = rcodcl(ifac,iuma,3)

      call set_neumann_scalar &
           !==================
         ( coefa(ifac,icluma), coefa(ifac,iclumf),             &
           coefb(ifac,icluma), coefb(ifac,iclumf),             &
           qimp              , hint )

      qimp = rcodcl(ifac,ivma,3)

      call set_neumann_scalar &
           !==================
         ( coefa(ifac,iclvma), coefa(ifac,iclvmf),             &
           coefb(ifac,iclvma), coefb(ifac,iclvmf),             &
           qimp              , hint )

      qimp = rcodcl(ifac,iwma,3)

      call set_neumann_scalar &
           !==================
         ( coefa(ifac,iclwma), coefa(ifac,iclwmf),             &
           coefb(ifac,iclwma), coefb(ifac,iclwmf),             &
           qimp              , hint )


      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then

        qimpv(1) = rcodcl(ifac,iuma,3)
        qimpv(2) = rcodcl(ifac,ivma,3)
        qimpv(3) = rcodcl(ifac,iwma,3)

        call set_neumann_vector &
             !==================
           ( claale(1,ifac)  , cfaale(1,ifac)  ,             &
             clbale(1,1,ifac), cfbale(1,1,ifac),             &
             qimpv           , hint )

      endif

    ! Convective Boundary Conditions
    !-------------------------------

    elseif (icodcl(ifac,iuma).eq.2) then

      pimp = rcodcl(ifac,iuma,1)
      cfl = rcodcl(ifac,iuma,2)

      call set_convective_outlet_scalar &
           !==================
         ( coefa(ifac,icluma), coefa(ifac,iclumf),             &
           coefb(ifac,icluma), coefb(ifac,iclumf),             &
           pimp              , cfl               , hint )

      pimp = rcodcl(ifac,ivma,1)
      cfl = rcodcl(ifac,ivma,2)

      call set_convective_outlet_scalar &
           !==================
         ( coefa(ifac,iclvma), coefa(ifac,iclvmf),             &
           coefb(ifac,iclvma), coefb(ifac,iclvmf),             &
           pimp              , cfl               , hint )

      pimp = rcodcl(ifac,iwma,1)
      cfl = rcodcl(ifac,iwma,2)

      call set_convective_outlet_scalar &
           !==================
         ( coefa(ifac,iclwma), coefa(ifac,iclwmf),             &
           coefb(ifac,iclwma), coefb(ifac,iclwmf),             &
           pimp              , cfl               , hint )


      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then

        pimpv(1) = rcodcl(ifac,iuma,1)
        cflv(1) = rcodcl(ifac,iuma,2)
        pimpv(2) = rcodcl(ifac,ivma,1)
        cflv(2) = rcodcl(ifac,ivma,2)
        pimpv(3) = rcodcl(ifac,iwma,1)
        cflv(3) = rcodcl(ifac,iwma,2)

        call set_convective_outlet_vector &
             !==================
           ( claale(1,ifac)  , cfaale(1,ifac)  ,             &
             clbale(1,1,ifac), cfbale(1,1,ifac),             &
             pimpv           , cflv            , hint )

      endif

    endif

  enddo

endif

!===============================================================================
! 14. Compute stresses at boundary (step 1 over 5)
!===============================================================================

if (ineedf.eq.1 .and. iterns.eq.1) then

  if (ivelco.eq.0) then
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      visclc = propce(iel,ipproc(iviscl))
      visctc = propce(iel,ipproc(ivisct))
      if (itytur.eq.3) then
        vistot = visclc
      else
        vistot = visclc + visctc
      endif
      distbf = distb(ifac)
      srfbnf = surfbn(ifac)
      forbr(1,ifac) = ( coefa(ifac,icluf)                          &
                      + coefb(ifac,icluf)*velipb(ifac,1) )*srfbnf
      forbr(2,ifac) = ( coefa(ifac,iclvf)                          &
                      + coefb(ifac,iclvf)*velipb(ifac,2) )*srfbnf
      forbr(3,ifac) = ( coefa(ifac,iclwf)                          &
                      + coefb(ifac,iclwf)*velipb(ifac,3) )*srfbnf
    enddo

  ! Coupled solving of the velocity components
  else
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      visclc = propce(iel,ipproc(iviscl))
      visctc = propce(iel,ipproc(ivisct))
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
      forbr(3,ifac) = ( coefau(3,ifac)                              &
                      + cofbfu(3,1,ifac) * velipb(ifac,1)           &
                      + cofbfu(3,2,ifac) * velipb(ifac,2)           &
                      + cofbfu(3,3,ifac) * velipb(ifac,3) )*srfbnf
    enddo
  endif
endif

! Free memory
deallocate(velipb)
if (allocated(rijipb)) deallocate(rijipb)

!===============================================================================
! 15. Formats
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
 8000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : COND. LIM.                                  ',/,&
'@    =========                                               ',/,&
'@     Le scalaire ',I10   ,' est couple a SYRTHES            ',/,&
'@      mais n''est pas la variable energetique               ',/,&
'@         ISCALT = ',I10                               ,/,&
'@                                                            ',/,&
'@     Le calcul ne sera pas execute.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

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
 8000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: BOUNDARY CONDITIONS                            ',/,&
'@    ========                                                ',/,&
'@     The scalar ',I10   ,' is coupled with SYRTHES          ',/,&
'@      but is not the energy variable                        ',/,&
'@         ISCALT = ',I10                               ,/,&
'@                                                            ',/,&
'@     The calculation will not be run.                       ',/,&
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

! Gradients BCs
coefa = pimp
coefb = 0.d0

! Flux BCs
cofaf = qimp
cofbf = 0.d0

return
end subroutine

!===============================================================================
