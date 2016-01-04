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

subroutine rayout

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  1) ECRITURE FICHIER SUITE,

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use pointe
use ppppar
use ppthch
use cpincl
use ppincl
use radiat
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

character        rubriq*64
integer          itysup, nbval, ifac
integer          ival(1)
double precision rval(1)

double precision, dimension(:), pointer :: btemp_s
double precision, allocatable, dimension(:) :: tb_save

type(c_ptr) :: rp

!===============================================================================
! 1. ECRITURE DU FICHIER SUITE DU MODULE DE RAYONNEMENT
!===============================================================================

! Open output

write(nfecra,6010)

call restart_create('radiative_transfer', '', 1, rp)

write(nfecra,6011)

! Entete et Dimensions ou on saute si erreur
!     On inclut une rubrique destinee a distinguer ce fichier
!       d'un autre fichier suite
!     Pour le moment, IVERS n'est pas utilise

itysup = 0
nbval  = 1
rubriq = 'version_fichier_suite_rayonnement'
ival(1) = 400000
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

write(nfecra,6012)

! Temps (par securite)

rubriq = 'nbre_pas_de_temps'
itysup = 0
nbval  = 1
ival(1) = ntcabs
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

rubriq = 'instant_precedent'
itysup = 0
nbval  = 1
rval(1) = ttcabs
call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

! Boundary values

if (itpscl.eq.1) then
  call restart_write_field_vals(rp, itempb, 0)
else
  allocate(tb_save(nfabor))
  call field_get_val_s(itempb, btemp_s)
  do ifac = 1, nfabor
    tb_save(ifac) = btemp_s(ifac) + tkelvi
  enddo
  rubriq = 'boundary_temperature'
  itysup = 3
  nbval = 1
  call restart_write_section_real_t(rp,rubriq,itysup,nbval,tb_save)
  deallocate(tb_save)
endif

call restart_write_field_vals(rp, iqinci, 0)
call restart_write_field_vals(rp, ihconv, 0)
call restart_write_field_vals(rp, ifconv, 0)

! Cell values

call restart_write_field_vals(rp, iprpfl(itsri(1)), 0)
call restart_write_field_vals(rp, iprpfl(itsre(1)), 0)
call restart_write_field_vals(rp, iprpfl(ilumin), 0)

write(nfecra,6013)

! Close file

call restart_destroy(rp)

write(nfecra,6014)

return

!--------
! Formats
!--------

 6010 format (/, &
           3X,'** Information on the radiative module',/,  &
           3X,'   -----------------------------------',/,  &
           3X,' Writing a restart file', /)

 6011 format (   3x,'   Write start',                   /)
 6012 format (   3x,'   End of output for dimensions',  /)
 6013 format (   3x,'   End of output for data',        /)
 6014 format (   3x,' End of output to restart file',   /)

!----
! End
!----

end subroutine rayout
