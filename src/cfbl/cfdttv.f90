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

subroutine cfdttv &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   wcf    ,                                                       &
   wflmas , wflmab , viscb  )

!===============================================================================
! FUNCTION :
! ----------

! Computation of the constraint for the CFL (compressible algorithm)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. cs_user_mass_source_terms)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,nvar)    !    !     !  source de masse                               !
!                  !    !     ! pour ivar=ipr, smacel=flux de masse            !
! wcf(ncelet)      ! tr ! --> ! contrainte compressible                        !
! wflmas(nfac)     ! tr ! --- ! tab de trav aux faces internes                 !
! wflmab(nfabor    ! tr ! --- ! tab de trav aux faces de bord                  !
! viscb(nfabor     ! tr ! --- ! tab de trav aux faces de bord                  !
!__________________!____!_____!________________________________________________!


!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe, only:rvoid1
use numvar
use cstnum
use cstphy
use optcal
use entsor
use parall
use ppppar
use ppthch
use ppincl
use mesh
use field
use cs_cf_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision wcf(ncelet)
double precision wflmas(nfac), wflmab(nfabor), viscb(nfabor)

! Local variables

integer          ifac, iel, iterns, idtcfl
integer          iconvp, idiffp, isym

double precision, allocatable, dimension(:) :: viscf
double precision, allocatable, dimension(:) :: coefbt, cofbft
double precision, allocatable, dimension(:) :: w1, c2

double precision, dimension(:,:), pointer :: vela
double precision, dimension(:), pointer :: crom, cpro_cp, cpro_cv
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:), pointer :: cvar_fracv, cvar_fracm, cvar_frace

!===============================================================================

! Map field arrays
call field_get_val_prev_v(ivarfl(iu), vela)

call field_get_val_s(icrom, crom)

call field_get_val_s(ivarfl(ipr), cvar_pr)

if (icfhgn.gt.0) then
   call field_get_val_s(ivarfl(isca(ifracv)), cvar_fracv)
   call field_get_val_s(ivarfl(isca(ifracm)), cvar_fracm)
   call field_get_val_s(ivarfl(isca(ifrace)), cvar_frace)
else
  cvar_fracv => null()
  cvar_fracm => null()
  cvar_frace => null()
endif

!===============================================================================
! 0.  INITIALIZATION
!===============================================================================

! Allocate temporary arrays
allocate(viscf(nfac))
allocate(coefbt(nfabor),cofbft(nfabor))

! Allocate work arrays
allocate(w1(ncelet))

!===============================================================================
! 1. COMPUTATION OF THE CFL CONDITION ASSOCIATED TO THE PRESSURE EQUATION
!===============================================================================

! Map specific heats fields for sound celerity computation

if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
else
  cpro_cp => rvoid1
endif

if (icv.ge.0) then
  call field_get_val_s(icv, cpro_cv)
else
  cpro_cv => rvoid1
endif

! Computation of the convective flux associated to the density

do ifac = 1, nfac
  wflmas(ifac) = 0.d0
enddo
do ifac = 1, nfabor
  wflmab(ifac) = 0.d0
enddo

iterns = 1
idtcfl = 1
call cfmsfp                                                       &
!==========
 ( nvar   , nscal  , idtcfl , iterns , ncepdp , ncesmp ,          &
   icepdc , icetsm , itypsm ,                                     &
   dt     , vela   ,                                              &
   ckupdc , smacel ,                                              &
   wflmas , wflmab )


! Summation at each cell taking only outward flux

iconvp = 1
idiffp = 0
isym   = 2

do ifac = 1, nfac
  viscf(ifac) = 0.d0
enddo
do ifac = 1, nfabor
  coefbt(ifac) = 0.d0
  cofbft(ifac) = 0.d0
  viscb(ifac)  = 0.d0
enddo

call matrdt &
 ( iconvp  , idiffp  , isym   , coefbt , cofbft , wflmas ,          &
   wflmab , viscf  , viscb  , w1 )

! Compute the square of the sound celerity

allocate(c2(ncelet))

call cs_cf_thermo_c_square(cpro_cp, cpro_cv, cvar_pr, crom,         &
                           cvar_fracv, cvar_fracm, cvar_frace, c2, ncel)

! Compute the coefficient CFL/dt

if (iporos.ge.1) then
  do iel = 1, ncel
    if (isolid_0(iel).eq.1) then
      wcf(iel) = epzero
    else
      wcf(iel) =  w1(iel) * c2(iel) * crom(iel)                      &
                / ((cvar_pr(iel) + psginf)*cell_f_vol(iel))
    endif
  enddo
else
  do iel = 1, ncel
    wcf(iel) =  w1(iel) * c2(iel) * crom(iel)                        &
              / ((cvar_pr(iel) + psginf)*cell_f_vol(iel))
  enddo
endif

! Free memory
deallocate(viscf)
deallocate(w1)
deallocate(c2)
deallocate(coefbt,cofbft)

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
