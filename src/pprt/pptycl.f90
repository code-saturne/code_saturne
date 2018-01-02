!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!===============================================================================
! Function :
! --------

!> \file pptycl.f90
!>
!> \brief Boundary conditions for specific physics modules.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in,out] icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[in,out] itypfb        boundary face types
!> \param[in,out] izfppp        index of the zone for the boundary faces
!>                               (for the specific physics)
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!_______________________________________________________________________________


subroutine pptycl &
 ( nvar   ,                                                       &
   icodcl , itypfb , izfppp ,                                     &
   dt     ,                                                       &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

! Arguments

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar

integer          icodcl(nfabor,nvar)
integer          itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, iok, ifvu, ii, izone, izonem

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================



!===============================================================================
! 2.  LISTE DES ZONES (pour n'importe quel modele)
!===============================================================================

! --> Faces appartiennent toutes a une zone frontiere

iok = 0

do ifac = 1, nfabor
  if(izfppp(ifac).le.0.or.izfppp(ifac).gt.nozppm) then
    iok = iok + 1
    write(nfecra,1000)ifac,nozppm,izfppp(ifac)
  endif
enddo

if(iok.gt.0) then
  call csexit (1)
  !==========
endif

! --> On construit une liste des numeros des zones frontieres.
!           (liste locale au processeur, en parallele)
nzfppp = 0
do ifac = 1, nfabor
  ifvu = 0
  do ii = 1, nzfppp
    if (ilzppp(ii).eq.izfppp(ifac)) then
      ifvu = 1
    endif
  enddo
  if(ifvu.eq.0) then
    nzfppp = nzfppp + 1
    if(nzfppp.le.nbzppm) then
      ilzppp(nzfppp) = izfppp(ifac)
    else
      write(nfecra,1001) nbzppm
      write(nfecra,1002)(ilzppp(ii),ii=1,nbzppm)
      call csexit (1)
      !==========
    endif
  endif
enddo

! ---> Plus grand numero de zone

izonem = 0
do ii = 1, nzfppp
  izone = ilzppp(ii)
  izonem = max(izonem,izone)
enddo
if(irangp.ge.0) then
  call parcmx(izonem)
  !==========
endif
nozapm = izonem

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PHYSIQUE PARTICULIERE                       ',/,&
'@    =========                                               ',/,&
'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES ',/,&
'@                                                            ',/,&
'@  Le numero de zone associee a la face ',I10   ,' doit etre ',/,&
'@    un entier strictement positif et inferieur ou egal a    ',/,&
'@    NOZPPM = ',I10                                           ,/,&
'@  Ce numero (IZFPPP(IFAC)) vaut ici ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites.                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PHYSIQUE PARTICULIERE                       ',/,&
'@    =========                                               ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@  Le nombre maximal de zones frontieres qui peuvent etre    ',/,&
'@    definies par l''utilisateur est NBZPPM = ',I10           ,/,&
'@    Il a ete depasse.                                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites.                      ',/,&
'@                                                            ',/,&
'@  Les NBZPPM premieres zones frontieres                     ',/,&
'@    portent ici les numeros suivants :                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(i10)



!===============================================================================
! 5.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LES CONDITIONS AUX LIMITES
!           POUR LES SCALAIRES
!    (selon le modele)
!===============================================================================


! ---> Chimie 3 points : USD3PC

if (ippmod(icod3p).ge.0) then

  call d3ptcl(itypfb, izfppp, icodcl, rcodcl)
  !==========

! ---> Combustion gaz USEBUC
!      Flamme de premelange modele EBU

elseif (ippmod(icoebu).ge.0) then

  call ebutcl(itypfb, izfppp, rcodcl)
  !==========

! ---> Combustion gaz USLWCC
!      Flamme de premelange modele LWC

elseif (ippmod(icolwc).ge.0) then

  call lwctcl(itypfb, izfppp, rcodcl)
  !==========

! ---> Combustion charbon pulverise

elseif ( ippmod(iccoal).ge.0 ) then

  call cs_coal_bcond(itypfb, izfppp, icodcl, rcodcl)
  !=================

! ---> Combustion charbon pulverise couple Lagrangien USCPLC

elseif (ippmod(icpl3c).ge.0) then

  call cpltcl(itypfb, izfppp, rcodcl)
  !==========

! ---> Combustion fuel

elseif (ippmod(icfuel).ge.0) then

  call cs_fuel_bcond(itypfb, izfppp, icodcl, rcodcl)
  !=================

! ---> Compressible

elseif (ippmod(icompf).ge.0) then

  call cfxtcl                                                     &
  !==========
 ( nvar   ,                                                       &
   icodcl , itypfb ,                                              &
   dt     ,                                                       &
   rcodcl )

! ---> Ecoulements atmospheriques

elseif (ippmod(iatmos).ge.0) then

  call attycl(itypfb, izfppp, rcodcl)
  !==========

! ---> Ecoulements electrique
! nothing

! ---> Cooling towers

elseif (ippmod(iaeros).ge.0) then

  call cs_ctwr_bcond(itypfb, izfppp, icodcl, rcodcl)

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
