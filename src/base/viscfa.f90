!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine viscfa &
!================

 ( idbia0 , idbra0 ,                                              &
   nideve , nrdeve , nituse , nrtuse , imvisf ,                   &
   idevel , ituser , ia     ,                                     &
   vistot ,                                                       &
   viscf  , viscb  ,                                              &
   rdevel , rtuser , ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! imvisf           ! e  ! <-- ! methode de calcul de la visc face              !
!                  !    !     !  = 0 arithmetique                              !
!                  !    !     !  = 1 harmonique                                !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! vistot(ncelet    ! tr ! <-- ! valeur de la viscosite                         !
! viscf(nfac)      ! tr ! --> ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --> ! visc*surface/dist aux faces de bord            !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
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

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nideve , nrdeve , nituse , nrtuse , imvisf

integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision vistot(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          ifac, ii, jj
double precision visci, viscj, surfn, dist, pond

!===============================================================================

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(vistot)
  !==========
endif


if( imvisf.eq.0 ) then

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    visci = vistot(ii)
    viscj = vistot(jj)
    surfn = ra(isrfan-1+ifac)
    dist  = ra(idist-1+ifac)

    viscf(ifac) = 0.5d0*( visci +viscj )*surfn/dist

  enddo

else

  do ifac = 1,nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    visci = vistot(ii)
    viscj = vistot(jj)
    surfn = ra(isrfan-1+ifac)
    dist  = ra(idist-1+ifac)
    pond  = ra(ipond-1+ifac)

    viscf(ifac) =                                                 &
      visci*viscj / (pond*visci+(1.d0-pond)*viscj)  * surfn/dist

  enddo

endif

do ifac = 1, nfabor

  ii = ifabor(ifac)
  surfn = ra(isrfbn-1+ifac)
  dist  = ra(idistb-1+ifac)

  viscb(ifac) = vistot(ii)*surfn/dist

enddo

!----
! FIN
!----

return

end subroutine
