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

subroutine cs_coal_varini &
!========================

 ( nvar   , nscal  ,                                            &
   dt     )

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL
!    POUR LA PHYSIQUE PARTICULIERE : COMBUSTION CP
!    PENDANT DE USINIV.F

! Cette routine est appelee en debut de calcul (suite ou non)
!     avant le debut de la boucle en temps

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul,
!     les valeurs du pas de temps


! On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
!     ou relues d'un fichier suite
! On ne dispose des variables VISCLS, CP (quand elles sont
!     definies) que si elles ont pu etre relues dans un fichier
!     suite de calcul

! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
!     ET PAS ICI

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
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
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use mesh
use field

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

integer          iel, ige, mode, icla, icha, ifac

double precision t1init, h1init, coefe(ngazem)
double precision t2init
double precision f1mc(ncharm), f2mc(ncharm)
double precision xkent, xeent, d2s3

integer ioxy
double precision wmh2o,wmco2,wmn2,wmo2,dmas

double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi
double precision, dimension(:), pointer :: cvar_fb, cvar_omg
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:), pointer :: cvar_xchcl, cvar_xckcl, cvar_xnpcl
double precision, dimension(:), pointer :: cvar_xwtcl
double precision, dimension(:), pointer :: cvar_h2cl
double precision, dimension(:), pointer :: cvar_scalt, cvar_xch
double precision, dimension(:), pointer :: cvar_f1m, cvar_f2m
double precision, dimension(:), pointer :: cvar_f4m, cvar_f5m, cvar_f6m
double precision, dimension(:), pointer :: cvar_f7m, cvar_f8m, cvar_f9m
double precision, dimension(:), pointer :: cvar_fvp2m
double precision, dimension(:), pointer :: cvar_yco2, cvar_yhcn, cvar_ynh3
double precision, dimension(:), pointer :: cvar_yno, cvar_hox
double precision, dimension(:), pointer :: cpro_x1, bpro_x1

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
! 1. Initializations
!===============================================================================

! Massic fraction of gas
call field_get_val_s_by_name("x_c", cpro_x1)
call field_get_val_s_by_name("b_x_c", bpro_x1)

ipass = ipass + 1

d2s3 = 2.d0/3.d0

if (itytur.eq.2) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (itytur.eq.3) then
  call field_get_val_s(ivarfl(ir11), cvar_r11)
  call field_get_val_s(ivarfl(ir22), cvar_r22)
  call field_get_val_s(ivarfl(ir33), cvar_r33)
  call field_get_val_s(ivarfl(ir12), cvar_r12)
  call field_get_val_s(ivarfl(ir13), cvar_r13)
  call field_get_val_s(ivarfl(ir23), cvar_r23)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (iturb.eq.50) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
  call field_get_val_s(ivarfl(iphi), cvar_phi)
  call field_get_val_s(ivarfl(ifb), cvar_fb)
elseif (iturb.eq.60) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iomg), cvar_omg)
endif

call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
call field_get_val_s_by_name("x_c_h",cvar_xch)

if (noxyd.ge.2) then
  call field_get_val_s(ivarfl(isca(if4m)), cvar_f4m)
endif
if (noxyd.ge.3) then
  call field_get_val_s(ivarfl(isca(if5m)), cvar_f5m)
endif
if (ippmod(iccoal).ge.1) then
  call field_get_val_s(ivarfl(isca(if6m)), cvar_f6m)
endif
call field_get_val_s(ivarfl(isca(if7m)), cvar_f7m)
if (ihtco2.eq.1) then
  call field_get_val_s(ivarfl(isca(if8m)), cvar_f8m)
endif
if (ihth2o.eq.1) then
  call field_get_val_s(ivarfl(isca(if9m)), cvar_f9m)
endif
call field_get_val_s(ivarfl(isca(ifvp2m)), cvar_fvp2m)
if (ieqco2.ge.1) then
  call field_get_val_s(ivarfl(isca(iyco2)), cvar_yco2)
endif
if (ieqnox.eq.1) then
  call field_get_val_s(ivarfl(isca(iyhcn)), cvar_yhcn)
  call field_get_val_s(ivarfl(isca(iynh3)), cvar_ynh3)
  call field_get_val_s(ivarfl(isca(iyno)), cvar_yno)
  call field_get_val_s(ivarfl(isca(ihox)), cvar_hox)
endif

!===============================================================================
! 2. Variable initialization
!===============================================================================

! RQ IMPORTANTE : pour la combustion CP, 1 seul passage suffit

if ( isuite.eq.0 .and. ipass.eq.1 ) then

! --> Initialisation de k et epsilon

  xkent = 1.d-10
  xeent = 1.d-10

! ---- TURBULENCE

  if (itytur.eq.2) then

    do iel = 1, ncel
      cvar_k(iel)  = xkent
      cvar_ep(iel) = xeent
    enddo

  elseif (itytur.eq.3) then

    do iel = 1, ncel
      cvar_r11(iel) = d2s3*xkent
      cvar_r22(iel) = d2s3*xkent
      cvar_r33(iel) = d2s3*xkent
      cvar_r12(iel) = 0.d0
      cvar_r13(iel) = 0.d0
      cvar_r23(iel) = 0.d0
      cvar_ep(iel)  = xeent
    enddo

  elseif (iturb.eq.50) then

    do iel = 1, ncel
      cvar_k(iel)   = xkent
      cvar_ep(iel)  = xeent
      cvar_phi(iel) = d2s3
      cvar_fb(iel)  = 0.d0
    enddo

  elseif (iturb.eq.60) then

    do iel = 1, ncel
      cvar_k(iel)   = xkent
      cvar_omg(iel) = xeent/cmu/xkent
    enddo

  endif

  ! --> All the domain is filled with the first oxidizer at TINITK
  !                   ============================================

  ! ---- Computation of H1INIT and H2INIT

  t1init = t0
  t2init = t0

! ------ Variables de transport relatives a la phase solide :
!        calcul de H

  do icla = 1, nclacp
    icha = ichcor(icla)

    call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xchcl)
    call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xckcl)
    call field_get_val_s(ivarfl(isca(inp(icla))), cvar_xnpcl)
    call field_get_val_s(ivarfl(isca(ih2(icla))), cvar_h2cl)
    if ( ippmod(iccoal) .eq. 1 ) then
      call field_get_val_s(ivarfl(isca(ixwt(icla))), cvar_xwtcl)
    endif

    do iel = 1, ncel

      cvar_xchcl(iel) = zero
      cvar_xckcl(iel) = zero
      cvar_xnpcl(iel) = zero
      cvar_h2cl(iel) = zero
      if ( ippmod(iccoal) .eq. 1 ) then
        cvar_xwtcl(iel) = zero
      endif

    enddo
  enddo

  ! ------ Transported variables for the mix (solid+carrying gas)^2

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo

  ! Oxidizer are mix of O2, N2 (air), CO2 and H2O (recycled exhaust)
  ! the composition of the fisrt oxidiser is taken in account
  coefe(io2) = wmole(io2)*oxyo2(1)                                &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ih2o) = wmole(ih2o)*oxyh2o(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ico2) = wmole(ico2)*oxyco2(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(in2) = 1.d0-coefe(io2)-coefe(ih2o)-coefe(ico2)

  do icha = 1, ncharm
    f1mc(icha) = zero
    f2mc(icha) = zero
  enddo
!
  mode = -1
  call cs_coal_htconvers1(mode,h1init,coefe,f1mc,f2mc,t1init)
 !============================

  do iel = 1, ncel
    cvar_scalt(iel) = h1init
    cvar_xch(iel)   = h1init
  enddo

  ! ------ Transported variables for the mix (passive scalars, variance)
!        (scalaires passifs et variances associees)

  do icha = 1, ncharb
    call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m)
    call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m)
    do iel = 1, ncel
      cvar_f1m(iel) = zero
      cvar_f2m(iel) = zero
    enddo
  enddo

  do iel = 1, ncel

    if ( noxyd .ge. 2 ) then
      cvar_f4m(iel) = zero
    endif
    if ( noxyd .ge. 3 ) then
      cvar_f5m(iel) = zero
    endif
    if ( ippmod(iccoal) .ge. 1 ) then
      cvar_f6m(iel) = zero
    endif
!
    cvar_f7m(iel) = zero
!
    if ( ihtco2 .eq. 1 ) then
      cvar_f8m(iel) = zero
    endif
    if ( ihth2o .eq. 1 ) then
      cvar_f9m(iel) = zero
    endif
!
    cvar_fvp2m(iel) = zero
!
    if ( ieqco2.ge.1 ) then

      ioxy   =  1
      wmo2   = wmole(io2)
      wmco2  = wmole(ico2)
      wmh2o  = wmole(ih2o)
      wmn2   = wmole(in2)
      dmas = ( oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2               &
              +oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2 )
      xco2 = oxyco2(ioxy)*wmco2/dmas
      cvar_yco2(iel) = oxyco2(ioxy)*wmco2/dmas
    endif
    if ( ieqnox.eq.1 ) then
      cvar_yhcn(iel) = zero
! Initialisation du NH3
      cvar_ynh3(iel) = zero
      cvar_yno(iel) = zero
      cvar_hox(iel) = h1init
    endif

  enddo

endif

if (ipass.eq.1) then

  do iel = 1, ncel
    ! Initialization of the continuous mass fraction
    cpro_x1(iel) = 1.d0
  enddo

  ! Initialization of the continuous mass fraction AT the BCs
  do ifac = 1, nfabor
    bpro_x1(ifac) = 1.d0
  enddo

!===============================================================================
! 3. User initialization
!===============================================================================

  call cs_user_f_initialization &
  !==========================
( nvar   , nscal  ,                                            &
  dt     )

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
