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

subroutine elprop
!================

!===============================================================================
!  FONCTION  :
!  ---------

!     INIT DES POSITIONS DES VARIABLES D'ETAT POUR
!                LE MODULE ELECTRIQUE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ipropp           ! e  ! ->  ! numero de la derniere case utlisee             !
!                  !    !     ! dans ipproc, ipprob, ipprof                    !
! ipppst           ! e  ! <-- ! pointeur indiquant le rang de la               !
!                  !    !     !  derniere grandeur definie aux                 !
!                  !    !     !  cellules (rtp,propce...) pour le              !
!                  !    !     !  post traitement                               !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppincl
use elincl
use ihmpre
use field

!===============================================================================

implicit none


! Local variables

character(len=80) :: f_name, f_label
integer       idimve
integer       ipropp

ipropp = nproce
!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================

!     Pointeurs dans propce (ca n'implique pas qu'on ne calcule pas
!     les variables non definies ici)

! ---> Temperature en K

call add_property_field('temperature', 'Temper', itemp)

! ---> Puissance volumique dissipee par effet Joule W/m3

call add_property_field('joule_power', 'PuisJoul', iefjou)

! ---> Densite de courant electrique reelle A/m2

do idimve = 1, ndimve
  write(f_name,  '(a11,i1)')  'current_re_', idimve
  write(f_label, '(a7,i1.1)') 'Cour_re', idimve
  call add_property_field(f_name, f_label, idjr(idimve))
enddo

! Variables specifiques Effet Joule
! =================================

if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then

! ---> Densite de courant electrique imaginaire A/m2

  do idimve = 1, ndimve
    write(f_name,  '(a11,i1)')  'current_im_', idimve
    write(f_label, '(a7,i1.1)') 'CouImag', idimve
    call add_property_field(f_name, f_label, idji(idimve))
  enddo

endif


! Variables specifiques Arc Electrique
! ====================================

if (ippmod(ielarc).ge.1) then

! ---> Forces electromagnetiques de Laplace en N/m3

  do idimve = 1, ndimve
    write(f_name,  '(a14,i1)')  'laplace_force_', idimve
    write(f_label, '(a7,i1.1)') 'For_Lap', idimve
    call add_property_field(f_name, f_label, ilapla(idimve))
  enddo

! ---> Puissance volumique rayonnee W/m3
!      ou coefficient d'absorption

  if (ixkabe .eq.1) then
    call add_property_field('absorption_coeff', 'Coef_Abso', idrad)
  else if (ixkabe .eq.2) then
    call add_property_field('radiation_source', 'TS_radia', idrad)
  endif

endif

! Variables specifiques Conduction Ionique
! ========================================

if (ippmod(ielion).ge.1) then

! ---> Charge electrique volumique C/m3

  call add_property_field('elec_charge', 'Charge', iqelec)

endif

! ----  Nb de variables algebriques (ou d'etat)
!         propre a la physique particuliere NSALPP
!         total NSALTO

nsalpp = nproce - ipropp
nsalto = nproce

return
end subroutine
