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

subroutine cou1do &
!================

 ( nvar   , nscal  , nfpt1d ,                                     &
   ifpt1d , iclt1d ,                                              &
   tppt1d , tept1d , hept1d , fept1d , eppt1d ,                   &
   xlmbt1 , rcpt1d , dtpt1d , dt     , cvcst  ,                   &
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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nfpt1d           ! e  ! <-- ! nombre de faces avec module therm 1d           !
! ifpt1d           ! te ! <-- ! numero de la face en traitement                !
!                  !    !     ! thermique en paroi                             !
! iclt1d           ! te ! <-- ! type de condition limite                       !
! tppt1d           ! tr ! <-- ! temperature de paroi                           !
! tept1d           ! tr ! <-- ! temperature exterieure                         !
! hept1d           ! tr ! <-- ! coefficient d'echange exterieur                !
! fept1d           ! tr ! <-- ! flux exterieur                                 !
! eppt1d           ! tr ! <-- ! epaisseur de paroi                             !
! xlmbt1           ! tr ! <-- ! diffusivite thermique                          !
! rcpt1d           ! tr ! <-- ! rocp                                           !
! dtpt1d           ! tr ! <-- ! pas de temps                                   !
! cvcst            ! r  ! <-- ! chaleur specifique si constante                !
! hbord(nfabor)    ! ra ! <-> ! coefficients d'echange aux bords               !
! tbord(nfabor)    ! ra ! <-> ! temperatures aux bords                         !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use pointe, only: izft1d, itypfb
use mesh
use cs_cf_bindings
use radiat

!===============================================================================

implicit none

! Arguments
integer          nfpt1d
integer          nvar   , nscal

integer          ifpt1d(nfpt1d), iclt1d(nfpt1d)

double precision dt(ncelet)
double precision hbord(nfabor),tbord(nfabor)
double precision tppt1d(nfpt1d)
double precision tept1d(nfpt1d), hept1d(nfpt1d), fept1d(nfpt1d), eppt1d(nfpt1d)
double precision xlmbt1(nfpt1d), rcpt1d(nfpt1d), dtpt1d(nfpt1d)
double precision cvcst

!     VARIABLES LOCALES

integer          iappel
integer          ifac, iel , ii

integer          ivoid(1)
double precision rvoid(1)

double precision energ, cvt

double precision, dimension(:), allocatable :: wa
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cpro_cp, cpro_cv, cpro_rho
double precision, dimension(:), pointer :: bqinci, beps

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine b_h_to_t(h_b, t_b)

    use mesh, only: nfabor
    implicit none

    double precision, dimension(nfabor), intent(in) :: h_b
    double precision, dimension(nfabor), intent(out), target :: t_b

  end subroutine b_h_to_t

 end interface

!===============================================================================

! Conversion to temperature for enthalpy or energy
! (check for surface couplings to make sure it is needed)

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

if (itherm.eq.2) then

  if (icp.gt.0) then
    call field_get_val_s(iprpfl(icp), cpro_cp)
  endif

  ! Temperature near boundary faces
  allocate(wa(nfabor))
  call b_h_to_t(tbord, wa)

  do ii = 1, nfpt1d
    ifac = ifpt1d(ii)
    iel  = ifabor(ifac)
    tbord(ifac) = wa(ifac)
    if (icp.gt.0) then
      hbord(ifac) = hbord(ifac)*cpro_cp(iel)
    else
      hbord(ifac) = hbord(ifac)*cp0
    endif
  enddo

else if (itherm.eq.3) then

  call field_get_val_v(ivarfl(iu), vel)
  call field_get_val_s(icrom, cpro_rho)
  if (icv.gt.0) then
    call field_get_val_s(iprpfl(icv), cpro_cv)
  endif

  ! Epsilon sup for perfect gas at cells
  allocate(wa(ncelet))
  call cs_cf_thermo_eps_sup(cpro_rho, wa, ncel)

  do ii = 1, nfpt1d
    ifac  = ifpt1d(ii)
    iel   = ifabor(ifac)
    energ = tbord(ifac)
    cvt   = energ                                               &
                  -(0.5d0*(  vel(1,iel)**2                      &
                           + vel(2,iel)**2                      &
                           + vel(3,iel)**2)                     &
                    + wa(iel) )
    if (icv.gt.0) then
      tbord(ifac) = cvt/cpro_cv(iel)
      hbord(ifac) = hbord(ifac)*cpro_cv(iel)
    else
      tbord(ifac) = cvt/cvcst
      hbord(ifac) = hbord(ifac)*cvcst
    endif
  enddo

endif

!     Mise a jour des conditions aux limites externes du module 1D
iappel = 3

call  uspt1d &
 ( nvar   , nscal  , nfpt1d , iappel ,                            &
   ifpt1d , izft1d , ivoid  , iclt1d ,                            &
   tppt1d , rvoid  , rvoid  ,                                     &
   tept1d , hept1d , fept1d ,                                     &
   xlmbt1 , rcpt1d , dtpt1d ,                                     &
   dt     )

iappel = 3
call vert1d &
( nfabor , nfpt1d , iappel ,                                      &
  ifpt1d , ivoid  , iclt1d ,                                      &
  rvoid  , rvoid  ,                                               &
  xlmbt1 , rcpt1d , dtpt1d )

! coupling with radiative module
if (iirayo.ge.1) then

  call field_get_val_s_by_name('rad_incident_flux', bqinci)
  call field_get_val_s_by_name('emissivity', beps)

  do ii = 1, nfpt1d

    ifac = ifpt1d(ii)

    ! FIXME pour gerer les faces qui ne sont pas des parois ou beps n'est pas
    ! renseigne, par exemple si une meme couleur est utilisee pour designer
    ! plusieurs faces, paroi + autre
    if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then

      call tpar1d                                                   &
      ( ii-1      , iclt1d(ii), tbord(ifac), hbord(ifac),           &
      bqinci(ifac), beps(ifac),                                     &
      tept1d(ii), hept1d(ii), fept1d(ii) ,                          &
      xlmbt1(ii), rcpt1d(ii), dtpt1d(ii) , tppt1d(ii) )

    endif

  enddo

! without coupling with radiative module
else

  do ii = 1, nfpt1d

    ifac = ifpt1d(ii)

    call tpar1d                                                 &
    ( ii-1    , iclt1d(ii), tbord(ifac), hbord(ifac),           &
    zero      , zero      ,                                     &
    tept1d(ii), hept1d(ii), fept1d(ii) ,                        &
    xlmbt1(ii), rcpt1d(ii), dtpt1d(ii) , tppt1d(ii) )

  enddo

endif

if (itherm .gt. 1) deallocate(wa)

return
end subroutine
