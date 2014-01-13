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

subroutine eliniv &
!================

 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce )

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL
!    POUR LA PHYSIQUE PARTICULIERE : VERSIONS ELECTRIQUES
!    PENDANT DE USINIV.F

! Cette routine est appelee en debut de calcul (suite ou non)
!     avant le debut de la boucle en temps

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul,
!     les valeurs du pas de temps


! On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
!     ou relues d'un fichier suite
! On ne dispose des variables VISCLS, CP (quand elles sont
!     definies) que si elles ont pu etre relues dans un fichier
!     suite de calcul

! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
!     ET PAS ICI

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)

! Local variables

integer          iel, mode, idimve , iesp

double precision tinit, hinit, coefe(ngazem)
double precision xkent, xeent, d2s3

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

ipass = ipass + 1


d2s3 = 2.d0/3.d0

!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

! RQ IMPORTANTE : pour module electrique, 1 seul passage suffit

if ( isuite.eq.0 .and. ipass.eq.1 ) then

! --> Initialisation de k et epsilon pas standard


! ---- Initialisation de k et epsilon

  xkent = 1.d-10
  xeent = 1.d-10

  do iel = 1, ncel

! ---- TURBULENCE

    if (itytur.eq.2) then

      rtp(iel,ik)  = xkent
      rtp(iel,iep) = xeent

    elseif (itytur.eq.3) then

      rtp(iel,ir11) = d2s3*xkent
      rtp(iel,ir22) = d2s3*xkent
      rtp(iel,ir33) = d2s3*xkent
      rtp(iel,ir12) = 0.d0
      rtp(iel,ir13) = 0.d0
      rtp(iel,ir23) = 0.d0
      rtp(iel,iep)  = xeent

    elseif (iturb.eq.50) then

      rtp(iel,ik)   = xkent
      rtp(iel,iep)  = xeent
      rtp(iel,iphi) = d2s3
      rtp(iel,ifb)  = 0.d0

    elseif (iturb.eq.60) then

      rtp(iel,ik)   = xkent
      rtp(iel,iomg) = xeent/cmu/xkent

    elseif (iturb.eq.70) then

      rtp(iel,inusa) = cmu*xkent**2/xeent

    endif

  enddo

! --> Enthalpie = H(T0) ou 0

!     En arc electrique, on initialise tout le domaine de calcul
!       a T0 avec la 1ere espece
!     En Joule, on initialise l'enthalpie a zero, et il faudra
!       que l'utilisateur intervienne, avec sa loi T->H ou une tabulation.

!   -- Calculs de HINIT

  if ( ippmod(ielarc).ge.1 ) then
    mode = -1
    tinit = t0
    coefe(1) = 1.d0
    if ( ngazg .gt. 1 ) then
      do iesp = 2, ngazg
        coefe(iesp) = 0.d0
      enddo
    endif
    call elthht(mode,ngazg,coefe,hinit,tinit)
  else
    hinit = 0.d0
  endif

!    -- Valeurs de l'enthalpie

  do iel = 1, ncel
    rtp(iel,isca(iscalt)) = hinit
  enddo


! --> Fractions massiques = 1 ou 0

  if ( ngazg .gt. 1 ) then
    do iel = 1, ncel
      rtp(iel,isca(iycoel(1))) = 1.d0
    enddo
    do iesp = 2, ngazg-1
      do iel = 1, ncel
        rtp(iel,isca(iycoel(iesp))) = 0.d0
      enddo
    enddo
  endif


! --> Potentiels Electrique = 0

!     -- Potentiel Reel
  do iel = 1, ncel
    rtp(iel,isca(ipotr)) = 0.d0
  enddo

!     -- Potentiel Imaginaire (Joule)
  if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
    do iel = 1, ncel
      rtp(iel,isca(ipoti)) = 0.d0
    enddo
  endif

!     -- Potentiel vecteur (arc elec. 3D)
  if ( ippmod(ielarc).ge.2 ) then
    do idimve = 1, ndimve
      do iel = 1, ncel
        rtp(iel,isca(ipotva(idimve))) = 0.d0
      enddo
    enddo
  endif


! --> Termes sources = 0

!     -- Puissance volumique dissipee par effet Joule W/m3
  do iel = 1, ncel
    propce(iel,ipproc(iefjou)) = 0.d0
  enddo

!     -- Force de Laplace (arc elec.)
  if ( ippmod(ielarc).ge.1 ) then
    do idimve = 1, ndimve
      do iel = 1, ncel
        propce(iel,ipproc(ilapla(idimve))) = 0.d0
      enddo
    enddo
  endif

endif

!===============================================================================
! 2.  ON PASSE LA MAIN A L'UTILISATEUR
!===============================================================================

if (ipass.eq.1) then

  call cs_user_initialization &
  !==========================
( nvar   , nscal  ,                                            &
  dt     , rtp    , propce )

endif


!----
! FORMATS
!----

!----
! FIN
!----

return

end subroutine
