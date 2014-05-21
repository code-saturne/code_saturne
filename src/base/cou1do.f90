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

subroutine cou1do &
!================

 ( nvar   , nscal  , nfpt1d ,                                     &
   ientha , ifpt1d , iclt1d ,                                     &
   tppt1d , tept1d , hept1d , fept1d ,                            &
   xlmbt1 , rcpt1d , dtpt1d , dt     ,                            &
   hbord  , tbord  )

!===============================================================================
! FONCTION :
! ---------

! ECRITURE DE DONNEES RELATIVES A UN COUPLAGE AVEC SYRTHES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nfpt1d           ! e  ! <-- ! nombre de faces avec module therm 1d           !
! ientha           ! e  ! <-- ! 1 si tparoi est une enthalpie                  !
! ifpt1d           ! te ! <-- ! numero de la face en traitement                !
!                  !    !     ! thermique en paroi                             !
! iclt1d           ! te ! <-- ! type de condition limite                       !
! hbord            ! tr ! <-- ! coefficients d'echange aux bords               !
! (nfabor)         !    !     !                                                !
! tbord            ! tr ! <-- ! temperatures aux bords                         !
! (nfabor)         !    !     !                                                !
! tppt1d           ! tr ! <-- ! temperature de paroi                           !
! tept1d           ! tr ! <-- ! temperature exterieure                         !
! hept1d           ! tr ! <-- ! coefficient d'echange exterieur                !
! fept1d           ! tr ! <-- ! flux exterieur                                 !
! xlmbt1           ! tr ! <-- ! diffusivite thermique                          !
! rcpt1d           ! tr ! <-- ! rocp                                           !
! dtpt1d           ! tr ! <-- ! pas de temps                                   !
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
use field
use parall
use period
use pointe, only: izft1d
use mesh

!===============================================================================

implicit none

! Arguments
integer          nfpt1d
integer          nvar   , nscal

integer          ifpt1d(nfpt1d), iclt1d(nfpt1d)
integer          ientha

double precision dt(ncelet)
double precision hbord(nfabor),tbord(nfabor)
double precision tppt1d(nfpt1d)
double precision tept1d(nfpt1d), hept1d(nfpt1d), fept1d(nfpt1d)
double precision xlmbt1(nfpt1d), rcpt1d(nfpt1d), dtpt1d(nfpt1d)

!     VARIABLES LOCALES

integer          mode
integer          iappel
integer          ifac, iel , ii

integer          ivoid(1)

double precision enthal, temper

double precision rvoid(1)
double precision, dimension(:), pointer :: cpro_cp

!===============================================================================

if(icp.gt.0) then
  call field_get_val_s(iprpfl(icp), cpro_cp)
endif

!     SI ENTHALPIE, ON TRANSFORME EN TEMPERATURE
!     Il est necessaire de transmettre a SYRTHES des Temperatures
!     Afin de conserver le flux Phi = (lambda/d     ) Delta T
!     ou Phi = (lambda/(d Cp)) Delta H
!     on multiplie HBORD = lambda/(d Cp) par Cp pris dans la
!     cellule adjacente.
!     Le resultat n'est pas garanti (conservation en particulier),
!     on ajoute donc un avertissement.

!     On ne change les TBORD et HBORD que sur les faces couplees. Du coup ces
!     tableaux contiennent des choses differentes suivant les faces.
!     C'est dangereux mais pas trop grave car on les jette juste apres
!     (COUPBO passe avant).

if(ientha.eq.1) then
   write(nfecra,1000)
   mode = 1
   do ii = 1, nfpt1d
      ifac = ifpt1d(ii)
      iel  = ifabor(ifac)
      enthal = tbord(ifac)
      call usthht (mode   , enthal , temper  )
      !==========
      tbord(ifac) = temper
      if(icp.gt.0) then
         hbord(ifac) = hbord(ifac)*cpro_cp(iel)
      else
         hbord(ifac) = hbord(ifac)*cp0
      endif
   enddo
endif

!     Pour l'instant on bloque le couplage si la variable est l'energie
!     -> on pourra le calquer sur coupbo si necessaire.
if (ientha.eq.2) then
  write(nfecra,2000)
  call csexit(1)
endif

!     Mise a jour des conditions aux limites externes du module 1D
iappel = 3

call  uspt1d &
!===========
 ( nvar   , nscal  , nfpt1d , iappel ,                            &
   ifpt1d , izft1d , ivoid  , iclt1d ,                            &
   tppt1d , rvoid  , rvoid  ,                                     &
   tept1d , hept1d , fept1d ,                                     &
   xlmbt1 , rcpt1d , dtpt1d ,                                     &
   dt     )

iappel = 3
call vert1d &
!==========
( nfabor , nfpt1d , iappel ,                                      &
  ifpt1d , ivoid  , iclt1d ,                                      &
  rvoid  , rvoid  ,                                               &
  xlmbt1 , rcpt1d , dtpt1d )

do ii = 1, nfpt1d

   ifac = ifpt1d(ii)

   call tpar1d                                                    &
   !==========
 ( ii-1      , iclt1d(ii), tbord(ifac), hbord(ifac),              &
   tept1d(ii), hept1d(ii), fept1d(ii) , xlmbt1(ii) ,              &
   rcpt1d(ii), dtpt1d(ii), tppt1d(ii) )


enddo

! --> FORMATS

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@ @@ ATTENTION : COUPLAGE SYRTHES/MODULE 1D AVEC CALCUL EN   ',/,&
'@                ENTHALPIE                                   ',/,&
'@    =========                                               ',/,&
'@      OPTION NON VALIDEE - CONTACTER L''EQUIPE DE DVPT      ',/,&
'@  ')
 2000 format(                                                           &
'@                                                            ',/,&
'@ @@ ATTENTION : COUPLAGE SYRTHES/MODULE 1D AVEC CALCUL EN   ',/,&
'@                ENERGIE                                   ',/,  &
'@    =========                                               ',/,&
'@      OPTION NON PERMISE - CONTACTER L''EQUIPE DE DVPT      ',/,&
'@                                                            ',/,&
'@      Le calcul s''arrete                                   ',/,&
'@  ')

#else

 1000 format(                                                           &
'@                                                            ',/,&
'@ @@ WARNING: 1D MODULE COUPLING WITH ENTHALPY CALCULATION   ',/,&
'@    ========                                                ',/,&
'@      OPTION NOT VALIDATED - CONTACT THE SUPPORT            ',/,&
'@                                                            ')
 2000 format(                                                           &
'@                                                            ',/,&
'@ @@ WARNING: 1D MODULE COUPLING WITH ENERGY CALCULATION     ',/,&
'@    ========                                                ',/,&
'@      OPTION NOT AVAILABLE - CONTACT THE SUPPORT            ',/,&
'@                                                            ',/,&
'@      The calculation will not be run                       ',/,&
'@  ')

#endif

return
end subroutine
