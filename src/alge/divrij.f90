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

subroutine divrij &
!================

 ( idim   , ivar   ,                                              &
   rtpa   ,                                                       &
   viscf  , viscb  )

!===============================================================================
! FONCTION :
! ---------

! DISPOSANT DU TENSEUR Rij
!  ON CALCULE LE TERME EN DIV INTERVENANT DANS L'EQUATION
!    DE LA VITESSE
!  ON PRODUIT DONC SOMME (Rij)kl Skl nkl
!    (Rij)kl EST LA VALEUR A LA FACE kl
!       Skl  EST LA SURFACE DE LA FACE kl
!       nkl  EST LE VECTEUR NORMAL A kl NORME
!       ON SOMME SUR TROIS COMPOSANTES DU TENSEUR
!  ON OBTIENT DONC UNE VALEUR PAR FACE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idim             ! e  ! <-- ! composante traitee                             !
! ivar             ! e  ! <-- ! numero de variable courante                    !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant prec)                     !
! viscf(nfac)      ! tr ! --> ! resultat du calcul                             !
! viscb(nfabor)    ! tr ! --> ! resultat du calcul                             !
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
use cstphy
use optcal
use pointe
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          idim   , ivar

double precision rtpa(ncelet,*)
double precision viscf(nfac), viscb(nfabor)

! Local variables

integer          ifac, ivar1, ivar2, ivar3, init, inc
integer          iccocg,iflmb0
integer          nswrgp, imligp, iwarnp
integer          imaspe, itypfl
double precision epsrgp, climgp, extrap
double precision, dimension(:), pointer :: crom, brom
double precision, dimension(:), pointer :: cofad1, cofad2, cofad3
double precision, dimension(:), pointer :: cofbd1, cofbd2, cofbd3

!===============================================================================

!===============================================================================
! 1. Initialisation
!===============================================================================

! Initialization to avoid compiler warnings
ivar1 = 0
ivar2 = 0
ivar2 = 0

! --- Density
call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

! --- Variables locales (Rij)
if(ivar.eq.iu) then
   ivar1 = ir11
   ivar2 = ir12
   ivar3 = ir13
elseif(ivar.eq.iv) then
   ivar1 = ir12
   ivar2 = ir22
   ivar3 = ir23
elseif(ivar.eq.iw) then
   ivar1 = ir13
   ivar2 = ir23
   ivar3 = ir33
endif

! --- Boundary conditions on the component Rij for the momentum equation

call field_get_coefad_s(ivarfl(ivar1), cofad1)
call field_get_coefad_s(ivarfl(ivar2), cofad2)
call field_get_coefad_s(ivarfl(ivar3), cofad3)

call field_get_coefbd_s(ivarfl(ivar1), cofbd1)
call field_get_coefbd_s(ivarfl(ivar2), cofbd2)
call field_get_coefbd_s(ivarfl(ivar3), cofbd3)

!===============================================================================
! 2. Calcul de la divergence
!===============================================================================

! --- Options de calcul
init = 1
inc  = 1
iccocg = 1
iflmb0 = 0
nswrgp = nswrgr(ir11)
imligp = imligr(ir11)
iwarnp = iwarni(ir11)
epsrgp = epsrgr(ir11)
climgp = climgr(ir11)
extrap = extrag(ir11)

imaspe = 2
itypfl = 1

call inimas                                                       &
!==========
 ( ivar1  , ivar2  , ivar3  , imaspe , itypfl ,                   &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   crom   , brom   ,                                              &
   rtpa(1,ivar1)   , rtpa(1,ivar2)   , rtpa(1,ivar3)   ,          &
   cofad1 , cofad2 , cofad3 , cofbd1 , cofbd2 , cofbd3 ,          &
   viscf  , viscb  )

!     Calcul des efforts aux bords (partie 5/5), si necessaire

if (ineedf.eq.1) then
  do ifac = 1, nfabor
    forbr(idim,ifac) = forbr(idim,ifac) + viscb(ifac)
  enddo
endif

return
end subroutine
