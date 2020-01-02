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

subroutine cpltsv &
!================

 ( iscal  , iscala ,                                              &
   itypfb ,                                                       &
   smbrs  , rovsdt )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!      TERMES SOURCES DE PRODUCTION ET DE DISSIPATION POUR
!      LA VARIANCE (BILANS EXPLICITE ET IMPLICITE)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iscal            ! i  ! <-- ! scalar number                                  !
! iscala           ! e  ! <-- ! numero du scalaire associe                     !
! itypfb           ! ia ! <-- ! boundary face types                            !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
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
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal  , iscala

integer          itypfb(nfabor)

double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

integer          ivar   , ivarsc , ivarut, f_id0
integer          iel, ifac
integer          icha
integer          inc , iccocg , nswrgp , imligp , iwarnp

double precision xk , xe , rhovst
double precision epsrgp , climgp , extrap
double precision turb_schmidt

double precision, allocatable, dimension(:) :: coefap, coefbp
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1, w2
double precision, allocatable, dimension(:) :: w7
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: visct
double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_omg
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: cvar_f1m, cvar_f2m, cvar_f3m
double precision, dimension(:), pointer :: cvar_varsc, cvara_varsc

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Initialize variables to avoid compiler warnings

xe = 0.d0
xk = 0.d0

! Memoire


! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar   = isca(iscal)

! --- Numero du scalaire eventuel associe dans le cas fluctuation
!         ISCALA et numero de variable de calcul
if (iscala.gt.0) then
  ivarsc = isca(iscala)
else
  ivarsc = 0
endif

! --- Numero des grandeurs physiques
call field_get_val_s(icrom, crom)
call field_get_val_s(ivisct, visct)

if ( itytur.eq.2 .or. iturb.eq.50 ) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
elseif ( itytur.eq.3 ) then
  call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
  call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
  call field_get_val_prev_s(ivarfl(ir33), cvara_r33)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
elseif ( iturb.eq.60 ) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iomg), cvara_omg)
endif

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES DE PRODUCTION PAR LES GRADIENTS
!    ET DE DISSIPATION
!===============================================================================

if ( itytur.eq.2 .or. itytur.eq.3                   &
     .or. iturb.eq.50 .or. iturb.eq.60 ) then

  inc = 1
  iccocg = 1
  if (ivarsc.gt.0) then
    ivarut = ivarsc
  else
! A defaut de savoir pour F4M on prend comme pour F3M
    ivarut = isca(if3m)
  endif
  call field_get_key_struct_var_cal_opt(ivarfl(ivarut), vcopt)
  nswrgp = vcopt%nswrgr
  imligp = vcopt%imligr
  iwarnp = vcopt%iwarni
  epsrgp = vcopt%epsrgr
  climgp = vcopt%climgr
  extrap = vcopt%extrag

  ! Allocate work arrays
  allocate(w7(ncelet))

  do iel = 1, ncel
    w7(iel) = zero
  enddo

! ---- W7 = FJM (kg/kg du melange gazeux)

  if (ivarsc.eq.0) then
    ! Allocate work arrays
    allocate(w1(ncelet), w2(ncelet))
    do iel = 1, ncel
      w1(iel) = zero
      w2(iel) = zero
    enddo
    do icha = 1, ncharb
      call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m)
      call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m)
      do iel = 1, ncel
        w1(iel) =  w1(iel) + cvar_f1m(iel)
        w2(iel) =  w2(iel) + cvar_f2m(iel)
      enddo
    enddo
    call field_get_val_s(ivarfl(isca(if3m)), cvar_f3m)
    do iel = 1, ncel
      w7(iel) = 1.d0 - (w1(iel) + w2(iel) + cvar_f3m(iel))
    enddo
    ! Free some work arrays
    deallocate(w1, w2)
  else
    call field_get_val_s(ivarfl(ivarsc), cvar_varsc)
    do iel = 1, ncel
      w7(iel) = cvar_varsc(iel)
    enddo
  endif

  ! Allocate temporary arrays
  allocate(coefap(nfabor), coefbp(nfabor))

  do ifac = 1, nfabor
    coefap(ifac) = zero
    coefbp(ifac) = 1.d0
    if ( itypfb(ifac).eq.ientre ) then
      coefap(ifac) = zero
      coefbp(ifac) = zero
      if (ivarsc.eq.0) coefap(ifac) = 1.d0
    endif
  enddo

  ! Allocate a temporary array for gradient computation
  allocate(grad(3,ncelet))

!  f_id0 = -1 (indique pour la periodicite de rotation que la variable
!              n'est pas Rij)
  f_id0 = -1

  call gradient_s                                                 &
 ( f_id0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , epsrgp , climgp , extrap ,                            &
   w7     , coefap , coefbp ,                                     &
   grad   )

  ! Free memory
  deallocate(coefap, coefbp)

  call field_get_val_prev_s(ivarfl(ivarsc), cvara_varsc)
  call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)

  do iel = 1, ncel
    if ( itytur.eq.2 .or. iturb.eq.50 ) then
      xk = cvara_k(iel)
      xe = cvara_ep(iel)
    elseif ( itytur.eq.3 ) then
      xk = 0.5d0*(cvara_r11(iel)+cvara_r22(iel)+cvara_r33(iel))
      xe = cvara_ep(iel)
    elseif ( iturb.eq.60 ) then
      xk = cvara_k(iel)
      xe = cmu*xk*cvara_omg(iel)
    endif

    rhovst = crom(iel)*xe/                                           &
             (xk * rvarfl(iscal))*volume(iel)
    rovsdt(iel) = rovsdt(iel) + max(zero,rhovst)
    smbrs(iel) = smbrs(iel) +                                        &
                2.d0*visct(iel)*volume(iel)/turb_schmidt             &
                * (grad(1,iel)**2 + grad(2,iel)**2 + grad(3,iel)**2) &
                - rhovst*cvara_varsc(iel)
  enddo

  ! Free memory
  deallocate(grad)
  deallocate(w7)

endif

!--------
! FORMATS
!--------



!----
! FIN
!----

return

end subroutine
