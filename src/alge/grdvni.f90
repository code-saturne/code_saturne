!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine grdvni &
!================

 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   vel    , coefav , coefbv ,                                     &
   gradv  )

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
!                  !    !     !   on pourra donner ivar=0 si le                !
!                  !    !     !   numero de variable n'est pas connu           !
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
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! vel(3,ncelet)    ! tr ! <-- ! variable (vitesse)                             !
! coefav,coefbv    ! tr ! <-- ! tableaux des cond lim pour pvar                !
!   (3,nfabor)     !    !     !  sur la normale a la face de bord              !
! gradv            ! tr ! --> ! gradient d'un vecteur                          !
!   (ncelet,3,3)   !    !     !                                                !
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
use pointe
use parall
use period
use mesh
use cstphy
use cstnum
use albase
use cplsat
use dimens, only: ndimfb

!===============================================================================

implicit none

! Arguments

integer          ivar   , imrgra , inc    , iccocg , nswrgp
integer          imligp ,iwarnp  , nfecra
double precision epsrgp , climgp , extrap

double precision vel(3*ncelet)
double precision coefav(*), coefbv(*)
double precision gradv(3*3*ncelet)

! Local variables

integer          iel, isou, ivarloc
integer          iphydp, ipond
integer          idimtr, irpvar
integer          iiu,iiv,iiw
integer          imlini

double precision rvoid(1)

!===============================================================================

! Note: this function is now incompatible with periodicity of rotation, which
!       requires coupled velocity components.

idimtr = 0

!===============================================================================
! 1. COMPUTATION OF THE GRADIENT
!===============================================================================

! This subroutine is never used to compute the pressure gradient
iphydp = 0
ipond = 0

ivarloc = ivar

call cgdcel &
!==========
 ( ivarloc, imrgra , inc    , iccocg , nswrgp ,                   &
   idimtr , iphydp , ipond  , iwarnp , imligp , epsrgp , extrap , &
   climgp , isympa , rvoid  ,                                     &
   coefav(1)       , coefbv(1)       , vel(1) , rvoid  ,          &
   gradv(1)     )

ivarloc = ivarloc+1

call cgdcel &
!==========
 ( ivarloc, imrgra , inc    , iccocg , nswrgp ,                   &
   idimtr , iphydp , ipond  , iwarnp , imligp , epsrgp , extrap , &
   climgp , isympa , rvoid  ,                                     &
   coefav(1+ndimfb), coefbv(1+ndimfb), vel(1+ncelet)   , rvoid  , &
   gradv(1+3*ncelet)     )

ivarloc = ivarloc+1

call cgdcel &
!==========
 ( ivarloc, imrgra , inc    , iccocg , nswrgp ,                   &
   idimtr , iphydp , ipond  , iwarnp , imligp , epsrgp , extrap , &
   climgp , isympa , rvoid  ,                                     &
   coefav(1+2*ndimfb), coefbv(1+2*ndimfb), vel(1+2*ncelet),       &
   rvoid , gradv(1+6*ncelet)  )

return
end subroutine
