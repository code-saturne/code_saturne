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

subroutine viscfa &
!================

 ( imvisf ,                                                       &
   vistot ,                                                       &
   viscf  , viscb  )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DE LA VITESSE DE DIFFUSION SUR LES FACETTES
! VISCF,B = VISCOSITE*SURFACE/DISTANCE, HOMOGENE A UN DEBIT EN KG/S

! RQE : A PRIORI, PAS BESOIN DE TECHNIQUE DE RECONSTRUCTION
!  ( A AMELIORER SI NECESSAIRE )

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! imvisf           ! e  ! <-- ! methode de calcul de la visc face              !
!                  !    !     !  = 0 arithmetique                              !
!                  !    !     !  = 1 harmonique                                !
! vistot(ncelet    ! tr ! <-- ! valeur de la viscosite                         !
! viscf(nfac)      ! tr ! --> ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --> ! surface aux faces de bord                      !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use field
use numvar, only: ipori
use paramx
use optcal, only: iporos
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          imvisf


double precision vistot(ncelet)
double precision viscf(nfac), viscb(nfabor)

! Local variables

integer          ifac, ii, jj
double precision visci, viscj, surfn, pnd
double precision, dimension(:), pointer :: porosi

!===============================================================================


! ---> Periodicity and parallelism treatment

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(vistot)
endif

! Without porosity
if (iporos.ne.1) then
  if (imvisf.eq.0) then

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      visci = vistot(ii)
      viscj = vistot(jj)

      viscf(ifac) = 0.5d0*(visci+viscj)*surfan(ifac)/dist(ifac)

    enddo

  else

    do ifac = 1,nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      visci = vistot(ii)
      viscj = vistot(jj)
      pnd  = pond(ifac)

      viscf(ifac) = visci*viscj                                          &
                  / (pnd*visci+(1.d0-pnd)*viscj)*surfan(ifac)/dist(ifac)

    enddo

  endif

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    viscb(ifac) = surfbn(ifac)

  enddo

! With porosity
else

  call field_get_val_s(ipori, porosi)

  if (imvisf.eq.0) then

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      visci = vistot(ii) * porosi(ii)
      viscj = vistot(jj) * porosi(jj)

      viscf(ifac) = 0.5d0*(visci+viscj)*surfan(ifac)/dist(ifac)

    enddo

  else

    do ifac = 1,nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      visci = vistot(ii) * porosi(ii)
      viscj = vistot(jj) * porosi(jj)
      surfn = surfan(ifac)
      pnd  = pond(ifac)
      viscf(ifac) = visci*viscj &
                  / (pnd*visci+(1.d0-pnd)*viscj)*surfan(ifac)/dist(ifac)

    enddo

  endif

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    viscb(ifac) = surfbn(ifac)*porosi(ii)

  enddo

endif

!----
! End
!----

return

end subroutine
