!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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
   iclip  ,                                                       &
   rtpa   , rtp    )

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
! iclip            ! e  ! <-- ! indicateur = 1 on n'utilise pas rtpa           !
!                  !    !     !  (inivar)                                      !
!                  !    !     !            sinon on peut (turrij)              !
! rtpa             ! tr ! <-- ! tableaux des variables au pdt precedt          !
! (ncelet,nvar)    !    !     !                                                !
! rtp              ! tr ! <-- ! tableaux des variables au pdt courant          !
! (ncelet,nvar)    !    !     !                                                !
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

!===============================================================================

implicit none

! Arguments

integer          nvar, ncelet, ncel
integer          iclip
double precision rtpa(ncelet,nflown:nvar)
double precision rtp(ncelet,nflown:nvar)

! Local variables

integer          iel, ivar, ivar1, ivar2, isou
integer          iclrij(7)
double precision vmin(7), vmax(7), var, rijmin, varrel, und0, epz2

!===============================================================================

! Initialization to avoid compiler warnings

ivar = 0
ivar1 = 0
ivar2 = 0

! Une petite valeur pour eviter des valeurs exactement nulles.

epz2 = epzero**2

!===============================================================================
!  ---> Stockage Min et Max pour listing
!===============================================================================

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

  iclrij(isou) = 0
  vmin(isou) =  grand
  vmax(isou) = -grand
  do iel = 1, ncel
    var = rtp(iel,ivar)
    vmin(isou) = min(vmin(isou),var)
    vmax(isou) = max(vmax(isou),var)
  enddo
enddo

! ---> Clipping (modif pour eviter les valeurs exactement nulles)

if (iclip.eq.1) then

  do isou = 1, 3

    if(isou.eq.1) ivar=ir11
    if(isou.eq.2) ivar=ir22
    if(isou.eq.3) ivar=ir33

    do iel = 1, ncel
      if (rtp(iel,ivar).le.epz2) then
        iclrij(isou) = iclrij(isou) + 1
        rtp(iel,ivar) = epz2
      endif
    enddo

  enddo

  do iel = 1, ncel
    if (abs(rtp(iel,iep)).le.epz2) then
      iclrij(7) = iclrij(7) + 1
      rtp(iel,iep) = max(rtp(iel,iep),epz2)
    elseif(rtp(iel,iep).le.0.d0) then
      iclrij(7) = iclrij(7) + 1
      rtp(iel,iep) = abs(rtp(iel,iep))
    endif
  enddo

else

  varrel = 1.1d0

  do isou = 1, 3

    if(isou.eq.1) ivar=ir11
    if(isou.eq.2) ivar=ir22
    if(isou.eq.3) ivar=ir33

    do iel = 1, ncel
      if (abs(rtp(iel,ivar)).le.epz2) then
        iclrij(isou) = iclrij(isou) + 1
        rtp(iel,ivar) = max(rtp(iel,ivar),epz2)
      elseif(rtp(iel,ivar).le.0.d0) then
        iclrij(isou) = iclrij(isou) + 1
        rtp(iel,ivar) = min(abs(rtp(iel,ivar)), varrel*abs(rtpa(iel,ivar)))
      endif
    enddo

  enddo

  iclrij(7) = 0
  do iel = 1, ncel
    if (abs(rtp(iel,iep)).lt.epz2) then
      iclrij(7) = iclrij(7) + 1
      rtp(iel,iep) = max(rtp(iel,iep),epz2)
    elseif(rtp(iel,iep).le.0.d0) then
      iclrij(7) = iclrij(7) + 1
      rtp(iel,iep) = min(abs(rtp(iel,iep )), varrel*abs(rtpa(iel,iep )))
    endif
  enddo

endif

! On force l'inegalite de Cauchy Schwarz

do isou = 4, 6

  if(isou.eq.4) then
    ivar  = ir12
    ivar1 = ir11
    ivar2 = ir22
  elseif(isou.eq.5) then
    ivar  = ir13
    ivar1 = ir11
    ivar2 = ir33
  elseif(isou.eq.6) then
    ivar  = ir23
    ivar1 = ir22
    ivar2 = ir33
  endif
  und0 = 1.d0
  do iel = 1, ncel
    rijmin = sqrt(rtp(iel,ivar1)*rtp(iel,ivar2))
    if (rijmin.lt.abs(rtp(iel,ivar))) then
      rtp(iel,ivar) = sign(und0,rtp(iel,ivar)) * rijmin
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
