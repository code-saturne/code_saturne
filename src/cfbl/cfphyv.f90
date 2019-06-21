!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file cfphyv.f90
!> \brief Computation of variable physical properties for the specific physics
!> compressible.
!>
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine cfphyv

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
use ppincl
use mesh
use field
use cs_cf_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer :: iel, ifcven, ifclam
double precision, dimension(:), pointer :: cpro_venerg, cpro_lambda
double precision, dimension(:), pointer :: cpro_cp, cpro_cv
double precision, dimension(:), pointer :: mix_mol_mas

!===============================================================================

!===============================================================================
! 1. Update Lambda/Cv
!===============================================================================

! It has been checked before this subroutine that cv0 was non zero.
! If Cv is variable and zero, it is an error due to the user.
! Here a test is performed at each call (not optimal).
! If the diffusivity of the total energy is constant, then the thermal
! conductivity and the isochoric specific heat should be constant.

call field_get_key_int (ivarfl(isca(ienerg)), kivisl, ifcven)
if (ifcven.ge.0) then

  call field_get_val_s(ifcven, cpro_venerg)

  call field_get_key_int (ivarfl(isca(itempk)), kivisl, ifclam)
  if (ifclam.ge.0) then
    call field_get_val_s(ifclam, cpro_lambda)
    do iel = 1, ncel
      cpro_venerg(iel) = cpro_lambda(iel)
    enddo
  else
    do iel = 1, ncel
      cpro_venerg(iel) = visls0(itempk)
    enddo
  endif

  if (icv.ge.0) then
    call field_get_val_s(icp, cpro_cp)
    call field_get_val_s(icv, cpro_cv)
    call field_get_val_s_by_name("mix_mol_mas", mix_mol_mas)

    call cs_cf_thermo_cv(cpro_cp, mix_mol_mas, cpro_cv, ncel)

    do iel = 1, ncel
      if(cpro_cv(iel).le.0.d0) then
        write(nfecra,2000)iel,cpro_cv(iel)
        call csexit (1)
        !==========
      endif
    enddo

    do iel = 1, ncel
      cpro_venerg(iel) = cpro_venerg(iel) / cpro_cv(iel)
    enddo
  else

    do iel = 1, ncel
      cpro_venerg(iel) = cpro_venerg(iel) / cv0
    enddo

  endif

else

  visls0(ienerg) = visls0(itempk)/cv0

endif


!--------
! Formats
!--------

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP DURING EXECUTION (COMPRESSIBLE MODULE)   ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@  The isochoric specific heat has at leat one value         ',/,&
'@    negative or zero:                                       ',/,&
'@    cell    ',I10,   '  Cv = ',E18.9                         ,/,&
'@                                                            ',/,&
'@  The computation will not run further.                     ',/,&
'@                                                            ',/,&
'@  Check cs_user_physical_properties.'                        ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!----
! End
!----

return
end subroutine
