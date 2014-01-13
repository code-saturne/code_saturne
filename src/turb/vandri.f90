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

subroutine vandri &
!================

 (  ndim   , ncelet , ncel   , nfabor ,                           &
    itypfb , ifabor , ifapat ,                                    &
    xyzcen , cdgfbo , visvdr , yplusc , propce )

!===============================================================================
! FONCTION :
! ----------

! IMPOSITION D'UN AMORTISSEMENT DE TYPE VAN DRIEST POUR LA LES
! nut est amortie par (1-exp(-y+/d+))**2 ou d+ est mis par defaut a 26

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! itypfb           ! ia ! <-- ! boundary face types                            !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifapat           ! te ! <-- ! no de face de brd code 5 la + proche           !
! (ncelet)         !    !     !    (rij et echo de paroi      )                !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! visvdr(ncelet)   ! tr ! <-- ! viscosite dynamique ds les cellules            !
!                  !    !     !  de bord apres amortisst de v driest           !
! yplusc           ! tr ! <-- ! valeur de yplus aux cellules                   !
! (ncelet  )       !    !     !    dans le cas abs(icdpar).eq.1                !
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
use entsor
use cstphy
use parall
use pointe, only: uetbor
use field
!===============================================================================

implicit none

! Arguments

integer          ndim, ncelet , ncel   , nfabor
integer          itypfb(nfabor),ifabor(nfabor)
integer          ifapat(ncelet)

double precision xyzcen(ndim,ncelet),cdgfbo(ndim,nfabor)
double precision visvdr(ncelet)
double precision yplusc(ncelet)
double precision propce(ncelet,*)

! Local variables

integer          iel   , ifac  , ipcvis, ipcvst
double precision yplus , yminpa, viscos
double precision, dimension(:), pointer :: crom
!===============================================================================

call field_get_val_s(icrom, crom)
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)

!     Calcul direct de la distance a la paroi (non compatible parall/perio)
if(abs(icdpar).eq.2) then

!     En sequentiel, RAS
  if(irangp.lt.0) then
    do iel = 1, ncel
      ifac = ifapat(iel)
      viscos = propce(iel,ipcvis)/crom(iel)
      yminpa = sqrt((cdgfbo(1,ifac)-xyzcen(1,iel))**2             &
           +        (cdgfbo(2,ifac)-xyzcen(2,iel))**2             &
           +        (cdgfbo(3,ifac)-xyzcen(3,iel))**2)
      yplus = uetbor(ifac) * yminpa/ viscos
      propce(iel,ipcvst) = propce(iel,ipcvst)*                    &
           (1.0d0-exp(-yplus/cdries))**2
    enddo
!     En parallele, on n'amortit que la premiere maille de paroi :
!     dangereux mais a priori inutile (car l'utilisation de
!     ICDPAR=+/-2 en parallele est bloque dans verini)
  else
    write(nfecra,1000)
    do ifac = 1, nfabor
      if(itypfb(ifac).eq.iparoi .or.                        &
         itypfb(ifac).eq.iparug ) then
        iel = ifabor(ifac)
        viscos = propce(iel,ipcvis)/crom(iel)
        yminpa = sqrt((cdgfbo(1,ifac)-xyzcen(1,iel))**2           &
             +        (cdgfbo(2,ifac)-xyzcen(2,iel))**2           &
             +        (cdgfbo(3,ifac)-xyzcen(3,iel))**2)
        yplus = uetbor(ifac) * yminpa/ viscos
        propce(iel,ipcvst) = propce(iel,ipcvst)*                  &
             (1.0d0-exp(-yplus/cdries))**2
      endif
    enddo
endif

!     Nouveau mode de calcul : c'est plus simple
elseif(abs(icdpar).eq.1) then
  do iel = 1, ncel
    yplus = yplusc(iel)
    propce(iel,ipcvst) = propce(iel,ipcvst)*                      &
           (1.0d0-exp(-yplus/cdries))**2
  enddo
endif

!     Pour les cellules de paroi on remet la viscosite turbulente
!     qui avait ete amortie dans clptur et qui a servi a calculer
!     les conditions aux limites
do iel = 1, ncel
  if (visvdr(iel).gt.-900.d0)                               &
       propce(iel,ipcvst) = visvdr(iel)
enddo

!--------
! FORMATS
!--------
#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : DANS LE CAS DE LA LES AVEC AMORTISSEMENT    ',/,&
'@    =========                                               ',/,&
'@    L''AMORTISSEMENT DE VAN DRIEST N''EST FAIT QUE SUR LA   ',/,&
'@    PREMIERE CELLULE A LA PAROI EN CAS DE PARALLELISME      ',/,&
'@                                                            ',/,&
'@  Le calcul se poursuit.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                           &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: IN CASE OF LES WITH DAMPING'                    ,/,&
'@    ========'                                                ,/,&
'@    VAN DRIEST DAMPING IS ONLY EFFECTIVE ON THE FIRST CELL'  ,/,&
'@    OFF-WALL IN CASE OF PARALLELISM'                         ,/,&
'@'                                                            ,/,&
'@  The calculation will be run.'                              ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif
!----
! FIN
!----


return
end subroutine
