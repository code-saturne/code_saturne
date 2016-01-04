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

subroutine clprij &
!================

 ( ncelet , ncel   , nvar   ,                                     &
   iclip  )

!===============================================================================
! FONCTION :
! ----------

! CLIPPING DE Rij ET EPSILON

!-------------------------------------------------------------------------------
! Arguments
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! e  ! <-- ! nombre de cellules                             !
! nvar             ! e  ! <-- ! nombre de variables                            !
! iclip            ! e  ! <-- ! indicateur = 1 on n'utilise pas les champs au  !
!                  !    !     !     pas de temps precedent (inivar)            !
!                  !    !     !            sinon on peut (turrij)              !
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
use entsor
use numvar
use cstnum
use parall
use cs_c_bindings
use field

!===============================================================================

implicit none

! Arguments

integer          nvar, ncelet, ncel
integer          iclip

! Local variables

integer          iel, ivar, ivar1, ivar2, isou
integer          iclrij(7)
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
!  ---> Stockage Min et Max pour listing
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
    cvar_var => cvar_r13
  elseif(isou.eq.6) then
    cvar_var => cvar_r23
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
    cvar_var  => cvar_r13
    cvar_var1 => cvar_r11
    cvar_var2 => cvar_r33
  elseif(isou.eq.6) then
    cvar_var  => cvar_r23
    cvar_var1 => cvar_r22
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

! ---> Stockage nb de clippings pour listing

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
    ivar = ir13
  elseif(isou.eq.6) then
    ivar = ir23
  elseif(isou.eq.7) then
    ivar = iep
  endif
  call log_iteration_clipping_field(ivarfl(ivar), iclrij(isou), 0,  &
                                    vmin(isou:isou), vmax(isou:isou))

enddo

return

end subroutine
