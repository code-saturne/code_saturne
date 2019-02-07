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

!-------------------------------------------------------------------------------

!> \function clprij
!> \brief Clipping of the turbulent Reynods stress tensor and the turbulent
!> dissipation (segregated version).
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in]     iclip         indicator = 0 if viscl0 is used
!>                              otherwise viscl is used.
!______________________________________________________________________________

subroutine clprij &
 ( ncelet , ncel   ,                                              &
   iclip  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use numvar
use cstnum
use parall
use cs_c_bindings
use field

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel
integer          iclip

! Local variables

integer          iel, ivar, ivar1, ivar2, isou
integer          iclrij(7),icl_max(1)
double precision vmin(7), vmax(7), var, rijmin, varrel, und0, epz2

double precision, dimension(:), pointer :: cvar_ep, cvara_ep
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: cvar_var, cvar_var1, cvar_var2
double precision, dimension(:), pointer :: cvara_var

!===============================================================================

! Initialization to avoid compiler warnings

ivar = 0
ivar1 = 0
ivar2 = 0
icl_max(1) = 0
! Une petite valeur pour eviter des valeurs exactement nulles.

epz2 = epzero**2

call field_get_val_s(ivarfl(iep), cvar_ep)
call field_get_val_s(ivarfl(ir11), cvar_r11)
call field_get_val_s(ivarfl(ir22), cvar_r22)
call field_get_val_s(ivarfl(ir33), cvar_r33)
call field_get_val_s(ivarfl(ir12), cvar_r12)
call field_get_val_s(ivarfl(ir13), cvar_r13)
call field_get_val_s(ivarfl(ir23), cvar_r23)

!===============================================================================
!  ---> Stockage Min et Max pour log
!===============================================================================

do isou = 1, 7
  if    (isou.eq.1) then
    cvar_var => cvar_r11
  elseif(isou.eq.2) then
    cvar_var => cvar_r22
  elseif(isou.eq.3) then
    cvar_var => cvar_r33
  elseif(isou.eq.4) then
    cvar_var => cvar_r12
  elseif(isou.eq.5) then
    cvar_var => cvar_r23
  elseif(isou.eq.6) then
    cvar_var => cvar_r13
  elseif(isou.eq.7) then
    cvar_var => cvar_ep
  endif

  iclrij(isou) = 0
  vmin(isou) =  grand
  vmax(isou) = -grand
  do iel = 1, ncel
    var = cvar_var(iel)
    vmin(isou) = min(vmin(isou),var)
    vmax(isou) = max(vmax(isou),var)
  enddo
enddo

! ---> Clipping (modif pour eviter les valeurs exactement nulles)

if (iclip.eq.1) then

  do isou = 1, 3

    if(isou.eq.1) cvar_var => cvar_r11
    if(isou.eq.2) cvar_var => cvar_r22
    if(isou.eq.3) cvar_var => cvar_r33

    do iel = 1, ncel
      if (cvar_var(iel).le.epz2) then
        iclrij(isou) = iclrij(isou) + 1
        cvar_var(iel) = epz2
      endif
    enddo

  enddo

  do iel = 1, ncel
    if (abs(cvar_ep(iel)).le.epz2) then
      iclrij(7) = iclrij(7) + 1
      cvar_ep(iel) = max(cvar_ep(iel),epz2)
    elseif(cvar_ep(iel).le.0.d0) then
      iclrij(7) = iclrij(7) + 1
      cvar_ep(iel) = abs(cvar_ep(iel))
    endif
  enddo

else

  varrel = 1.1d0

  call field_get_val_prev_s(ivarfl(iep), cvara_ep)

  call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
  call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
  call field_get_val_prev_s(ivarfl(ir33), cvara_r33)

  do isou = 1, 3

    if    (isou.eq.1) then
      cvar_var => cvar_r11
      cvara_var => cvara_r11
    elseif(isou.eq.2) then
      cvar_var => cvar_r22
      cvara_var => cvara_r22
    elseif(isou.eq.3) then
      cvar_var => cvar_r33
      cvara_var => cvara_r33
    endif

    do iel = 1, ncel
      if (abs(cvar_var(iel)).le.epz2) then
        iclrij(isou) = iclrij(isou) + 1
        cvar_var(iel) = max(cvar_var(iel),epz2)
      elseif(cvar_var(iel).le.0.d0) then
        iclrij(isou) = iclrij(isou) + 1
        cvar_var(iel) = min(abs(cvar_var(iel)), varrel*abs(cvara_var(iel)))
      endif
    enddo

  enddo

  iclrij(7) = 0
  do iel = 1, ncel
    if (abs(cvar_ep(iel)).lt.epz2) then
      iclrij(7) = iclrij(7) + 1
      cvar_ep(iel) = max(cvar_ep(iel),epz2)
    elseif(cvar_ep(iel).le.0.d0) then
      iclrij(7) = iclrij(7) + 1
      cvar_ep(iel) = min(abs(cvar_ep(iel)), varrel*abs(cvara_ep(iel)))
    endif
  enddo

endif

! On force l'inegalite de Cauchy Schwarz

do isou = 4, 6

  if(isou.eq.4) then
    cvar_var  => cvar_r12
    cvar_var1 => cvar_r11
    cvar_var2 => cvar_r22
  elseif(isou.eq.5) then
    cvar_var  => cvar_r23
    cvar_var1 => cvar_r22
    cvar_var2 => cvar_r33
  elseif(isou.eq.6) then
    cvar_var  => cvar_r13
    cvar_var1 => cvar_r11
    cvar_var2 => cvar_r33
  endif
  und0 = 1.d0
  do iel = 1, ncel
    rijmin = sqrt(cvar_var1(iel)*cvar_var2(iel))
    if (rijmin.lt.abs(cvar_var(iel))) then
      cvar_var(iel) = sign(und0,cvar_var(iel)) * rijmin
      iclrij(isou) = iclrij(isou) + 1
    endif
  enddo

enddo

! ---> Stockage nb de clippings pour log

do isou = 1, 7
  if    (isou.eq.1) then
    ivar = ir11
  elseif(isou.eq.2) then
    ivar = ir22
  elseif(isou.eq.3) then
    ivar = ir33
  elseif(isou.eq.4) then
    ivar = ir12
  elseif(isou.eq.5) then
    ivar = ir23
  elseif(isou.eq.6) then
    ivar = ir13
  elseif(isou.eq.7) then
    ivar = iep
  endif
  call log_iteration_clipping_field(ivarfl(ivar), iclrij(isou), 0,  &
                                    vmin(isou:isou), vmax(isou:isou),iclrij(isou),&
                                    icl_max)

enddo
return

end subroutine clprij

!-------------------------------------------------------------------------------

!> \function clprij2
!> \brief Clipping of the turbulent Reynods stress tensor and the turbulent
!> dissipation (coupled version).
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in]     iclip         indicator = 0 if viscl0 is used
!>                              otherwise viscl is used.
!______________________________________________________________________________

subroutine clprij2 &
 ( ncelet , ncel   ,                                              &
   iclip  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use numvar
use cstnum
use parall
use cs_c_bindings
use field

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel
integer          iclip

! Local variables
integer          iel, ivar, ivar1, ivar2, isou, icltot, iclpep(1)
integer          is_clipped
integer          iclrij(6),iclrij_max(6), iclep_max(1)
integer          kclipp, clip_e_id, clip_r_id

double precision vmin(7), vmax(7), rijmin, varrel, und0, epz2,cvar_var1, cvar_var2
double precision eigen_min, eigen_max, trrij, eigen_offset
double precision eigen_vals(3)
double precision tensor(6)

double precision, parameter :: eigen_tol = 1.0d-4

double precision :: rijmax(3), rijref, trref

double precision, dimension(:), pointer :: cvar_ep, cvara_ep
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:,:), pointer :: cpro_rij_clipped
double precision, dimension(:), pointer :: cpro_eps_clipped

!===============================================================================

! Initialization to avoid compiler warnings

ivar = 0
ivar1 = 0
ivar2 = 0

! Une petite valeur pour eviter des valeurs exactement nulles.

epz2 = epzero**2

call field_get_val_s(ivarfl(iep), cvar_ep)
call field_get_val_v(ivarfl(irij), cvar_rij)

call field_get_key_id("clipping_id", kclipp)

! Postprocess clippings?
call field_get_key_int(ivarfl(iep), kclipp, clip_e_id)
if (clip_e_id.ge.0) then
  call field_get_val_s(clip_e_id, cpro_eps_clipped)
endif

call field_get_key_int(ivarfl(irij), kclipp, clip_r_id)
if (clip_r_id.ge.0) then
  call field_get_val_v(clip_r_id, cpro_rij_clipped)
endif

!===============================================================================
!  Compute and store Min Max values for the log
!===============================================================================

do isou = 1, 7
  vmin(isou) =  grand
  vmax(isou) = -grand
  do iel = 1, ncel
    if (isou.lt.7) then
      iclrij_max(isou) = 0
      vmin(isou) = min(vmin(isou),cvar_rij(isou,iel))
      vmax(isou) = max(vmax(isou),cvar_rij(isou,iel))
      if (clip_r_id.ge.0) &
        cpro_rij_clipped(isou, iel) = 0.d0
    else
      iclep_max(1) = 0
      vmin(isou) = min(vmin(isou),cvar_ep(iel))
      vmax(isou) = max(vmax(isou),cvar_ep(iel))
      if (clip_e_id.ge.0) &
        cpro_eps_clipped(iel) = 0.d0
    endif
  enddo
enddo

! ---> Clipping (modif pour eviter les valeurs exactement nulles)

varrel = 1.1d0

call field_get_val_prev_s(ivarfl(iep), cvara_ep)
icltot = 0
iclpep(1) = 0
do isou = 1, 6
  iclrij(isou) = 0
enddo

! Computing the maximal value of each of the diagonal components over the
! entire domain. A reference value "rijref", used to determine if a value is
! small is then calculated as: rijref = (r11max + r22max + r33max)/3.
! New test is rii < epzero*rijref
rijmax(:) = 0.0d0
do iel = 1, ncel
  do isou = 1, 3
    if (cvar_rij(isou,iel).gt.rijmax(isou)) then
      rijmax(isou) = cvar_rij(isou,iel)
    end if
  end do
end do

call parrmx(3, rijmax)

trref  = sum(rijmax)
rijref = max(trref/3.0d0,epzero)

do iel = 1, ncel

  is_clipped = 0

  ! Check if R is positive and ill-conditionned (since the former
  ! will induce the latter after clipping ...

  trrij = cvar_rij(1,iel) + cvar_rij(2,iel) + cvar_rij(3,iel)

  if (trrij.le.epzero*trref) then
    do isou = 1, 3
      if (clip_r_id.ge.0) then
        cpro_rij_clipped(isou, iel)  = cvar_rij(isou,iel) - epzero*rijref
        cpro_rij_clipped(isou+3,iel) = cvar_rij(isou+3,iel)
      end if

      cvar_rij(isou,iel)   = epzero*rijref
      cvar_rij(isou+3,iel) = 0.0d0

      iclrij(isou) = iclrij(isou) + 1
      iclrij(isou+3) = iclrij(isou+3) + 1
    end do

    is_clipped = 1

  else
    do isou = 1, 6
      tensor(isou) = cvar_rij(isou,iel)/trrij
    enddo

    call calc_symtens_eigvals(tensor,eigen_vals)

    eigen_min = minval(eigen_vals(1:3))
    eigen_max = maxval(eigen_vals(1:3))

    if ( (eigen_min .le. (eigen_tol*eigen_max)) .or. &
      (eigen_min .le. epzero) ) then

      is_clipped = 1

      eigen_offset = (eigen_tol - eigen_min) * trrij

      eigen_offset = max(eigen_offset, epzero*rijref)

      do isou = 1, 3
        cvar_rij(isou,iel) = cvar_rij(isou,iel) + eigen_offset

        if (clip_r_id.ge.0) &
          cpro_rij_clipped(isou, iel) = eigen_offset

        iclrij(isou) = iclrij(isou) + 1
      end do
    end if
  end if

  ! Epsilon
  if (abs(cvar_ep(iel)).lt.epz2) then
    iclpep(1) = iclpep(1) + 1
    if (clip_e_id.ge.0) &
      cpro_eps_clipped(iel) = abs(cvar_ep(iel) - epz2)
    cvar_ep(iel) = max(cvar_ep(iel),epz2)
  elseif(cvar_ep(iel).le.0.d0) then
    iclpep(1) = iclpep(1) + 1
    if (clip_e_id.ge.0) &
      cpro_eps_clipped(iel) = 2.d0*abs(cvar_ep(iel))
    cvar_ep(iel) = min(abs(cvar_ep(iel)), varrel*abs(cvara_ep(iel)))
  endif

  ! Enforced Cauchy Schwarz inequality (only for x, y, z direction)
  do isou = 4, 6
    if(isou.eq.4) then
      cvar_var1 = cvar_rij(1,iel)
      cvar_var2 = cvar_rij(2,iel)
    elseif(isou.eq.6) then
      cvar_var1 = cvar_rij(1,iel)
      cvar_var2 = cvar_rij(3,iel)
    elseif(isou.eq.5) then
      cvar_var1 = cvar_rij(2,iel)
      cvar_var2 = cvar_rij(3,iel)
    endif
    und0 = 1.d0

    rijmin = sqrt(cvar_var1*cvar_var2)
    if (rijmin.lt.abs(cvar_rij(isou,iel))) then
      is_clipped = 1
      if (clip_r_id.ge.0) &
        cpro_rij_clipped(isou, iel) = cvar_rij(isou,iel)
      cvar_rij(isou,iel) = sign(und0,cvar_rij(isou,iel)) &
        * rijmin / (1.d0 + epzero)
      iclrij(isou) = iclrij(isou) + 1
    endif
  enddo

  icltot = icltot + is_clipped
enddo

call log_iteration_clipping_field(ivarfl(irij), icltot, 0,      &
                                  vmin, vmax,iclrij,iclrij_max)
call log_iteration_clipping_field(ivarfl(iep), iclpep(1), 0,    &
                                  vmin(7), vmax(7),iclpep, iclep_max)

return

end subroutine clprij2
