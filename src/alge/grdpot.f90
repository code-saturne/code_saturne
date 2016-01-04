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

subroutine grdpot &
!================

 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , epsrgp , climgp , extrap ,                            &
   frcxt  ,                                                       &
   pvar   , coefap , coefbp ,                                     &
   grad   )

!===============================================================================
! FONCTION :
! ----------

! APPEL DES DIFFERENTES ROUTINES DE CALCUL DE GRADIENT CELLULE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ivar             ! e  ! <-- ! numero de la variable                          !
!                  !    !     !   destine a etre utilise pour la               !
!                  !    !     !   periodicite uniquement (pering)              !
!                  !    !     !   on pourra donner ivar=0 si la                !
!                  !    !     !   variable n'est ni une composante de          !
!                  !    !     !   la vitesse, ni une composante du             !
!                  !    !     !   tenseur des contraintes rij                  !
! imrgra           ! e  ! <-- ! methode de reconstruction du gradient          !
!                  !    !     !  0 reconstruction 97                           !
!                  !    !     !  1 moindres carres                             !
!                  !    !     !  2 moindres carres support etendu              !
!                  !    !     !    complet                                     !
!                  !    !     !  3 moindres carres avec selection du           !
!                  !    !     !    support etendu                              !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! iccocg           ! e  ! <-- ! indicateur = 1 pour recalcul de cocg           !
!                  !    !     !              0 sinon                           !
! nswrgp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligp           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! iphydp           ! e  ! <-- ! indicateur de prise en compte de la            !
!                  !    !     ! pression hydrostatique                         !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! pvar  (ncelet    ! tr ! <-- ! variable (pression)                            !
! coefap,coefbp    ! tr ! <-- ! tableaux des cond lim pour pvar                !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! ppond(ncelet)    ! tr ! <-- ! ponderation "physique"                         !
! frcxt            ! tr ! <-- ! force exterieure generant la pression          !
!                  !    !     !  hydrostatique                                 !
! grad(ncelet,3)   ! tr ! --> ! gradient de pvar                               !
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
use albase
use cplsat
use pointe
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          ivar   , imrgra , inc    , iccocg , nswrgp
integer          imligp , iwarnp , iphydp
double precision epsrgp , climgp , extrap

double precision frcxt(3,ncelet)
double precision pvar(ncelet), coefap(nfabor), coefbp(nfabor)
double precision grad(ncelet,3)

! Local variables

integer          imrgrp, ilved
integer          idimtr, ipond

double precision rvoid(1)

!===============================================================================

!===============================================================================
! 0. Initialization
!===============================================================================

! Use iterative gradient

if (imrgra.lt.0) then
  imrgrp = 0
else
  imrgrp = imrgra
endif

! The gradient of a potential (pressure, ...) is a vector

ilved = 0
idimtr = 0
ipond = 0

!===============================================================================
! 1. Compute gradient
!===============================================================================

call cgdcel &
!==========
 ( ivar   , imrgrp , ilved  , inc    , iccocg , nswrgp ,          &
   idimtr , iphydp , ipond  , iwarnp , imligp , epsrgp , extrap , &
   climgp , frcxt  , coefap , coefbp ,                            &
   pvar   , rvoid  , grad   )

return
end subroutine
