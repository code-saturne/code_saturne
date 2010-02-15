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

subroutine cou1di &
!================

 ( idbia0 , idbra0 ,                                              &
   nfabor ,                                                       &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   isvtb  , icodcl , idevel , ituser , ia     ,                   &
   rcodcl , rdevel , rtuser , ra     )

!===============================================================================

! FONCTION :
! ---------

! LECTURE DE DONNEES RELATIVES A UN COUPLAGE AVEC SYRTHES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! isvtb            ! e  ! <-- ! numero du scalaire couple                      !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail



!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstnum.h"
include "cstphy.h"
include "entsor.h"
include "pointe.h"

!===============================================================================

! Arguments

integer          idbia0, idbra0
integer          nfabor
integer          nvar , nscal , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          isvtb  , icodcl(nfabor,nvar)
integer          idevel(nideve) , ituser(nituse)
integer          ia(*)
double precision rcodcl(nfabor,nvar,3)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! Local variables


integer          ii , ivar
integer          ifac
integer          idebia, idebra
integer          icldef
integer          mode
double precision temper, enthal

!===============================================================================

idebia = idbia0
idebra = idbra0

!     Sans specification, une face couplee est une face de type paroi
icldef = 5

ivar = isca(isvtb)

do ii = 1, nfpt1d

   ifac = ia(iifpt1+ii-1)

   if ((icodcl(ifac,ivar) .ne. 1) .and.                           &
       (icodcl(ifac,ivar) .ne. 5) .and.                           &
       (icodcl(ifac,ivar) .ne. 6)) icodcl(ifac,ivar) = icldef

   rcodcl(ifac,ivar,1) = ra(itppt1+ii-1)
   rcodcl(ifac,ivar,2) = rinfin
   rcodcl(ifac,ivar,3) = 0.d0

enddo

! Conversion eventuelle temperature -> enthalpie

if (iscsth(isvtb).eq.2) then

  do ii = 1, nfpt1d

    ifac = ia(iifpt1+ii-1)

    temper = rcodcl(ifac,ivar,1)
    mode   = -1
    call usthht(mode,enthal,temper)
    !==========
    rcodcl(ifac,ivar,1) = enthal

  enddo

endif

end subroutine


