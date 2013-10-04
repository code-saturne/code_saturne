!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine pthrbm &
!================

 ( nvar   , ncesmp ,                                              &
   dt     , smacel )

!===============================================================================
! FONCTION :
! ----------

! Update the density rho^(n+1) with the rho^(n-1/2) density from the state law
! and the thermodynamic pressure pther^(n+1)
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nvar             ! e  ! <-- ! nombre total de variables                      !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate                     !
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
use pointe, only: itypfb, icetsm
use entsor
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , ncesmp

double precision dt(ncelet)
double precision smacel(ncesmp,nvar)

! VARIABLES LOCALES

logical          lromo

integer          iel , ifac, ieltsm

integer          iflmab

double precision ro0amoy,ro0moy
double precision debt

double precision debin, debout, debtot

double precision, dimension(:), pointer :: bmasfl
double precision, dimension(:), pointer :: brom, crom, cromo

!===============================================================================



!===============================================================================
! 0. Initialization
!===============================================================================

!   pointers for the different variables

call field_get_val_s(icrom, crom)

call field_have_previous(icrom, lromo)
if (lromo) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
endif

call field_get_val_s(ibrom, brom)

call field_get_key_int(ivarfl(ipr), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

!===============================================================================
! 1. Flow rate computation for the inlet and oulet conditions
!===============================================================================

!-----------------
!- Initialization
!-----------------

debin  = 0.d0
debout = 0.d0
debtot = 0.d0
!--------------

! Computation of mass flux imposed on the boundary faces
do ifac = 1, nfabor
  if (itypfb(ifac).eq.ientre) then
    ! the inlet mass flux integrated in space
    debin = debin + bmasfl(ifac)
  else if (itypfb(ifac).eq.isolib) then
    ! the outlet mass flux integrated in space:
    debout = debout + bmasfl(ifac)
  endif
enddo

debtot = debin + debout

! Computation of the inlet mass flux imposed on the cells volume
if (ncesmp.gt.0) then
  do ieltsm = 1, ncesmp
    iel = icetsm(ieltsm)
    debtot = debtot + smacel(ieltsm,ipr)*volume(iel)
  enddo
endif

if (irangp.ge.0) then
  call parsom (debtot)
endif

!===============================================================================
! 2. Thermodynamic pressure and density computation
!===============================================================================

! for the first time step : rho^(n-1) = rho^(n)
if (isuite.eq.0 .and. ntcabs.eq.1) then
  do iel = 1, ncel
    cromo(iel) = crom(iel)
  enddo
endif

!initialization
ro0amoy = 0.d0
ro0moy  = 0.d0

do iel = 1, ncel
 ro0amoy = ro0amoy + cromo(iel)*volume(iel)
 ro0moy  = ro0moy  + crom(iel)*volume(iel)
enddo

if (irangp.ge.0) then
  call parsom(ro0amoy)
  call parsom(ro0moy)
endif

! Compute the thermodynamic pressure p_ther^(n+1)
pther = pthera*(ro0amoy + dt(1)*debtot)/ro0moy

! Actualize the density rho[p_ther^(n+1)] = rho^(n+1)
do iel = 1, ncel
  crom(iel) = pther/pthera*crom(iel)
enddo

! Update the density at the boundary face
do ifac = 1, nfabor
  iel = ifabor(ifac)
  brom(ifac) = crom(iel)
enddo

!===============================================================================
! 3. Change the reference variable rho0
!===============================================================================

!initialization
ro0moy = 0.d0

do iel=1,ncel
  ro0moy = ro0moy + crom(iel)*volume(iel)
enddo

if (irangp.ge.0) then
  call parsom(ro0moy)
endif

ro0 = ro0moy/voltot

!===============================================================================
! 4.Printing
!===============================================================================

if (mod(ntcabs,ntlist).eq.0 .or. ntcabs.eq.1) then

  ! Compute the mass flux at the boundary faces
  debt = 0.d0
  do ifac = 1, nfabor
    debt = debt + bmasfl(ifac)
  enddo
  if (irangp.ge.0) then
    call parsom(debt)
  endif

  write (nfecra, 2002) ttcabs, pther, (pther-pthera)/dt(1), ro0, -debt

endif

!===============================================================================
! FORMATS
!----

  !================================================================
  2002 format                                                      &
  (/,                                                              &
   3X,'** LOW-MACH ALGORITHM: AVERAGED QUANTITIES '            , /,&
   3X,'   --------------------------------------- '            , /,&
   '---',                                                          &
   '-------------------------------------------------------',      &
   '-------------'                                             , /,&
   3X,'    Time      pther^(n+1)    Dp/Dt     ',                   &
   ' ro0_moy     mass_flux   '                                 , /,&
   '---',                                                          &
   '-------------------------------------------------------',      &
   '-------------'                                             , /,&
   3X, 5e12.4, /,                                                  &

   '---',                                                          &
   '-------------------------------------------------------',      &
   '-------------' )
  !================================================================


!----
! FIN
!----

return

end subroutine
