!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     iscal         index number of the current scalar   
!> \param[in]     dt            time step (per cell)
!> \param[out]    crvexp        explicit part of the source term
!_______________________________________________________________________________



!===============================================================================


subroutine precst &
!================

 ( nvar   , nscal  ,                                              &
   iscal  ,                                                       &
   dt     ,                                            &
   crvexp  )

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iscal            ! i  ! <-- ! index number of the current scalar             !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! crvexp           ! ra ! --> ! explicit part of the source term               !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use optcal
use cstnum
use mesh
use field
use lagran
use lagdim
use pointe,  only:  solub, nbprec, mp_diss

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          iscal

double precision dt(ncelet)
double precision crvexp(ncelet)

! Local variables

integer          ivar , ivar2, iel ,  k

double precision, dimension(:), pointer ::  cpro_rom

integer           npt

double precision , dimension(:), allocatable :: mp_prec
integer, dimension(:), allocatable :: part_tot

double precision, dimension(:), pointer :: cvar_scal
double precision, dimension(:), pointer :: cvar_scalt

!===============================================================================
! 1. Initialization
!===============================================================================
! --- Index number of the variable associated to scalar iscal
ivar = isca(iscal)
ivar2= isca(iscalt)

! --- Density
call field_get_val_s(icrom, cpro_rom)

call field_get_val_s(ivarfl(ivar), cvar_scal)
call field_get_val_s(ivarfl(ivar2), cvar_scalt)

allocate(mp_prec(ncelet))
allocate(part_tot(ncelet))
!===============================================================================
! 2. Calculation of the mass source terms due to
!    precipitation and dissolution phenomena
!===============================================================================

if (iscal.eq.1 .and. nbrclas .gt. 0 .and. iprec .eq. 1  ) then

   part_tot = 0
   if (associated(nbpart)) then
      do iel = 1 , ncel
         do npt  = 1 , nbpart
            if (ipepa(jisor,npt) .eq. iel  .and. &
                 eptp(jmp,npt) .eq. rho_prec * pi /6.d0 * eptp(jdp,npt)* & 
                 eptp(jdp,npt) * eptp(jdp,npt)) then
               ! number of magnetite particles in the cell iel
               part_tot(iel) = part_tot(iel) + 1 
            endif
         enddo
      enddo
   endif

   !Source term applied to second scalar

   mp_diss = 0.d0
   mp_prec = 0.d0

   crvexp = 0.d0

   do iel = 1 , ncel

      nbprec(iel) = 0

      !PRECIPITATION
      if (cvar_scal(iel) .ge. solub(iel) ) then
         nbprec(iel) = ((cvar_scal(iel) - solub(iel))* volume(iel)) & 
              / (pi/6.d0 * dprec**3 * rho_prec)
         crvexp(iel) = - cpro_rom(iel) * (nbprec(iel) * (pi/6.d0 * &
              dprec**3 * rho_prec)) / dtref
         mp_prec(iel) = nbprec(iel) * (pi/6.d0 * dprec**3 * rho_prec)
      endif
      !to do:  impose a limit on  nbprec

      !DISSOLUTION
      if (cvar_scal(iel) .lt. solub(iel) .and. part_tot(iel ) .ge. 1) then
         if (associated(nbpart)) then
            do npt  = 1 , nbpart
               do k = 1, nbrclas
                  if (ipepa(jisor,npt) .eq. iel .and. eptp(jdp,npt) .eq. ruslag(k,1,idpt) .and. &
                       eptp(jmp,npt) .eq. rho_prec * pi /6.d0 * eptp(jdp,npt) & 
                       * eptp(jdp,npt) * eptp(jdp,npt)) then
                     if ( (((solub(iel) - cvar_scal(iel)) * volume(iel)) -            &
                          (mp_diss(iel,k) + pepa(jrpoi,npt)* & 
                          ( pi / 6.d0 * (eptp(jdp,npt))**3 * rho_prec ))) .ge. 0) then
                        mp_diss(iel,k) = mp_diss(iel,k) + pepa(jrpoi,npt) * (pi/6.d0 & 
                             * (eptp(jdp,npt))**3 * rho_prec)
                     endif
                  endif
               enddo
            enddo
         endif

         do k = 1, nbrclas
            crvexp(iel) =  crvexp(iel) + (cpro_rom(iel) * mp_diss(iel,k)/ dtref)
         enddo
      endif

   enddo

endif

deallocate(mp_prec)
deallocate(part_tot)
!----
! End
!----

return
end subroutine precst


!===============================================================================


