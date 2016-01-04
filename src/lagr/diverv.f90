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

subroutine diverv &
!================

 ( div    , u      , coefa  , coefb  )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!        CALCUL DE LA DIVERGENCE D'UN VECTEUR

!   (On ne s'embete pas, on appelle 3 fois le gradient)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! div(ncelet)      ! tr ! --> ! divergence du vecteur                          !
! ux,uy,uz         ! tr ! --> ! composante du vecteur                          !
! (ncelet)         !    !     !                                                !
! coefax,...       ! tr ! ->  ! conditions aux limites pour les                !
! coefbz           !    !     ! faces de bord                                  !
! (nfabor)         !    !     !                                                !
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
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use lagpar
use lagran
use mesh
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

double precision, dimension(ncelet)     :: div
double precision, dimension(3,ncelet)   :: u
double precision, dimension(3,nfabor)   :: coefa
double precision, dimension(3,3,nfabor) :: coefb

! Local variables

integer          f_id0
integer          iel
integer          inc
integer          nswrgp, imligp, iwarnp

double precision epsrgp, climgp

double precision, allocatable, dimension(:,:,:) :: grad

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(grad(3,3,ncelet))

f_id0 = -1

inc = 1
nswrgp = 100
imligp = -1
iwarnp = 2
epsrgp = 1.d-8
climgp = 1.5d0

!===============================================================================
! Calcul du gradient de U
!===============================================================================

call cgdvec                                                        &
( f_id0  , imrgra , inc    , nswrgp , iwarnp , imligp ,            &
  epsrgp , climgp , coefa  , coefb  , u      , grad   )

!===============================================================================
! Calcul de la divergence du vecteur
!===============================================================================

do iel = 1,ncel
  div(iel) = grad(1,1,iel) + grad(2,2,iel) + grad(3,3,iel)
enddo

! Free memory
deallocate(grad)

!----
! End
!----

end subroutine
