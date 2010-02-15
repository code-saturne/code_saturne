!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine clprij &
!================

 ( ncelet , ncel   , nvar   , nphas  ,                            &
   iphas  , iclip  ,                                              &
   propce , rtpa   , rtp    )

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
! nphas            ! i  ! <-- ! number of phases                               !
! iphas            ! e  ! <-- ! numero de la phase a traiter                   !
! iclip            ! e  ! <-- ! indicateur = 1 on n'utilise pas rtpa           !
!                  !    !     !  (inivar)                                      !
!                  !    !     !            sinon on peut (turrij)              !
! propce           ! tr ! <-- ! tableaux des variables au pdt courant          !
!(ncelet,*         !    !     !                                                !
! rtpa             ! tr ! <-- ! tableaux des variables au pdt precedt          !
! (ncelet,nvar)    !    !     !                                                !
! rtp              ! tr ! <-- ! tableaux des variables au pdt courant          !
! (ncelet,nvar)    !    !     !                                                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "entsor.h"
include "numvar.h"
include "cstnum.h"
include "parall.h"

!===============================================================================

! Arguments

integer          nvar, ncelet, ncel, nphas
integer          iphas, iclip
double precision propce(ncelet,*)
double precision rtpa(ncelet,nvar)
double precision rtp(ncelet,nvar)

! Local variables

integer          icleps, iel, ivar, ivar1, ivar2, isou, ipp
integer          ir11ip, ir22ip, ir33ip, ir12ip, ir13ip, ir23ip
integer          iepip
integer          iclrij(6)
double precision vmin, vmax, var, rijmin, varrel, und0, epz2

!===============================================================================


ir11ip = ir11(iphas)
ir22ip = ir22(iphas)
ir33ip = ir33(iphas)
ir12ip = ir12(iphas)
ir13ip = ir13(iphas)
ir23ip = ir23(iphas)
iepip  = iep(iphas)

! Une petite valeur pour eviter des valeurs exactement nulles.

epz2 = epzero**2

!===============================================================================
!  ---> Stockage Min et Max pour listing
!===============================================================================

do isou = 1, 7
  if    (isou.eq.1) then
    ivar = ir11ip
  elseif(isou.eq.2) then
    ivar = ir22ip
  elseif(isou.eq.3) then
    ivar = ir33ip
  elseif(isou.eq.4) then
    ivar = ir12ip
  elseif(isou.eq.5) then
    ivar = ir13ip
  elseif(isou.eq.6) then
    ivar = ir23ip
  elseif(isou.eq.7) then
    ivar = iepip
  endif
  ipp = ipprtp(ivar)

  vmin =  grand
  vmax = -grand
  do iel = 1, ncel
    var = rtp(iel,ivar)
    vmin = min(vmin,var)
    vmax = max(vmax,var)
  enddo
  if (irangp.ge.0) then
    call parmin(vmin)
    !==========
    call parmax(vmax)
    !==========
  endif
  varmna(ipp) = vmin
  varmxa(ipp) = vmax

enddo

! ---> Clipping (modif pour eviter les valeurs exactement nulles)

if(iclip.eq.1) then

  do isou = 1, 3

    iclrij(isou) = 0

    if(isou.eq.1) ivar=ir11ip
    if(isou.eq.2) ivar=ir22ip
    if(isou.eq.3) ivar=ir33ip

    do iel = 1, ncel
      if (abs(rtp(iel,ivar)).le.epz2) then
        iclrij(isou) = iclrij(isou) + 1
        rtp(iel,ivar) = max(rtp(iel,ivar),epz2)
      elseif(rtp(iel,ivar).le.0.d0) then
        iclrij(isou) = iclrij(isou) + 1
        rtp(iel,ivar) = abs(rtp(iel,ivar))
      endif
    enddo

  enddo

  icleps = 0
  do iel = 1, ncel
    if (abs(rtp(iel,iepip)).le.epz2) then
      icleps = icleps + 1
      rtp(iel,iepip) = max(rtp(iel,iepip),epz2)
    elseif(rtp(iel,iepip).le.0.d0) then
      icleps = icleps + 1
      rtp(iel,iepip) = abs(rtp(iel,iepip))
    endif
  enddo

else

  varrel = 1.1d0

  do isou = 1, 3

    iclrij(isou) = 0

    if(isou.eq.1) ivar=ir11ip
    if(isou.eq.2) ivar=ir22ip
    if(isou.eq.3) ivar=ir33ip

    do iel = 1, ncel
      if (abs(rtp(iel,ivar)).le.epz2) then
        iclrij(isou) = iclrij(isou) + 1
        rtp(iel,ivar) = max(rtp(iel,ivar),epz2)
      elseif(rtp(iel,ivar).le.0.d0) then
        iclrij(isou) = iclrij(isou) + 1
        rtp(iel,ivar) =                                           &
          min(abs(rtp(iel,ivar)), varrel*abs(rtpa(iel,ivar)))
      endif
    enddo

  enddo

  icleps = 0
  do iel = 1, ncel
    if (abs(rtp(iel,iepip)).lt.epz2) then
      icleps = icleps + 1
      rtp(iel,iepip) = max(rtp(iel,iepip),epz2)
    elseif(rtp(iel,iepip).le.0.d0) then
      icleps = icleps + 1
      rtp(iel,iepip) =                                            &
          min(abs(rtp(iel,iepip )), varrel*abs(rtpa(iel,iepip )))
    endif
  enddo

endif

! On force l'inegalite de Cauchy Schwarz

do isou = 4, 6

  iclrij(isou) = 0

  if(isou.eq.4) then
    ivar  = ir12ip
    ivar1 = ir11ip
    ivar2 = ir22ip
  elseif(isou.eq.5) then
    ivar  = ir13ip
    ivar1 = ir11ip
    ivar2 = ir33ip
  elseif(isou.eq.6) then
    ivar  = ir23ip
    ivar1 = ir22ip
    ivar2 = ir33ip
  endif
  und0 = 1.d0
  do iel = 1 , ncel
    rijmin = sqrt(rtp(iel,ivar1)*rtp(iel,ivar2))
    if(rijmin.lt.abs(rtp(iel,ivar))) then
      rtp(iel,ivar) = sign(und0,rtp(iel,ivar)) * rijmin
      iclrij(isou) = iclrij(isou) + 1
    endif
  enddo

enddo


! ---> Stockage nb de clippings pour listing

if (irangp.ge.0) then
  call parcpt (iclrij(1))
  !==========
  call parcpt (iclrij(2))
  !==========
  call parcpt (iclrij(3))
  !==========
  call parcpt (iclrij(4))
  !==========
  call parcpt (iclrij(5))
  !==========
  call parcpt (iclrij(6))
  !==========
  call parcpt (icleps)
  !==========
endif

iclpmn(ipprtp(ir11ip)) = iclrij(1)
iclpmn(ipprtp(ir22ip)) = iclrij(2)
iclpmn(ipprtp(ir33ip)) = iclrij(3)
iclpmn(ipprtp(ir12ip)) = iclrij(4)
iclpmn(ipprtp(ir13ip)) = iclrij(5)
iclpmn(ipprtp(ir23ip)) = iclrij(6)
iclpmn(ipprtp(iepip )) = icleps

return

end subroutine
