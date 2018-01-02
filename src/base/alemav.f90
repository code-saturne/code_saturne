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

!===============================================================================
! Function:
! ---------

!> \file alemav.f90
!>
!> \brief This subroutine updates the mesh in the ALE framework.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     itrale        number of the current ALE iteration
!> \param[in,out] xyzno0        nodes coordinates of the initial mesh
!_______________________________________________________________________________

subroutine alemav &
 ( itrale ,                                                       &
   xyzno0 )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use dimens
use entsor
use cstphy
use cstnum
use pointe
use parall
use period
use mesh
use albase, only:fdiale
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itrale

double precision xyzno0(3,nnod)

! Local variables

integer          inod
integer          iel
integer          idim

double precision, dimension(:,:), pointer :: mshvel, mshvela
double precision, dimension(:,:), pointer :: disale, disala

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1. Initialisation
!===============================================================================

call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt)

if (vcopt%iwarni.ge.1) then
  write(nfecra,1000)
endif

call field_get_val_v(ivarfl(iuma), mshvel)
call field_get_val_prev_v(ivarfl(iuma), mshvela)

call field_get_val_v(fdiale, disale)
call field_get_val_prev_v(fdiale, disala)

!===============================================================================
! 2. Update geometry
!===============================================================================

do inod = 1, nnod
  do idim = 1, ndim
    xyznod(idim,inod) = xyzno0(idim,inod) + disale(idim,inod)
    disala(idim,inod) = xyznod(idim,inod) - xyzno0(idim,inod)
  enddo
enddo

call algrma(volmin, volmax, voltot)

! Abort at the end of the current time-step if there is a negative volume
if (volmin.le.0.d0) ntmabs = ntcabs

! Si on est a l'iteration d'initialisation, on remet les vitesses de maillage
!   a leur valeur initiale
if (itrale.eq.0) then
  do iel = 1, ncelet
    do idim = 1, ndim
      mshvel(idim, iel) = mshvela(idim,iel)
    enddo
  enddo
endif

!--------
! Formats
!--------

 1000 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  Update the mesh (ALE)'                                      ,/,&
'  ====================='                                      ,/)

!----
! End
!----

end subroutine
