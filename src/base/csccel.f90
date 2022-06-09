!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file csccel.f90
!> \brief Exchange of coupling variables between tow instances of code_saturne
!> thanks to cells.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     ivar          variable number
!> \param[out]    crvimp        working table for implicit part
!> \param[out]    crvexp        working table for explicit part
!______________________________________________________________________________

subroutine csccel &
 ( ivar   ,                                                       &
   crvexp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use cplsat
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          ivar

double precision crvexp(3,ncelet)

! Local variables

integer          f_id, f_dim
integer          numcpl
integer          ncesup , nfbsup
integer          ncecpl , nfbcpl , ncencp , nfbncp
integer          ncedis , nfbdis
integer          ncecpg , ncedig
integer          ityloc , ityvar

integer, allocatable, dimension(:) :: lcecpl , lfbcpl
integer, allocatable, dimension(:) :: locpts

double precision, allocatable, dimension(:,:) :: coopts , djppts , dofpts
double precision, allocatable, dimension(:) :: pndpts
double precision, allocatable, dimension(:,:) :: rvdis, rvcel

!===============================================================================

! cannot be called for non variables ...
if (ivar.le.0) then
  f_id = -1
  call csexit(1)
else
  f_id = ivarfl(ivar)
endif

!get the dimension of the variable
call field_get_dim(f_id, f_dim)

do numcpl = 1, nbrcpl

!===============================================================================
! 1.  DEFINITION DE CHAQUE COUPLAGE
!===============================================================================

  call nbecpl                                                     &
  !==========
 ( numcpl ,                                                       &
   ncesup , nfbsup ,                                              &
   ncecpl , nfbcpl , ncencp , nfbncp )

  ! Allocate temporary arrays for coupling information
  allocate(lcecpl(ncecpl))
  allocate(lfbcpl(nfbcpl))

  call lelcpl                                                     &
  !==========
 ( numcpl ,                                                       &
   ncecpl , nfbcpl ,                                              &
   lcecpl , lfbcpl )

  deallocate(lfbcpl)

!===============================================================================
! 2.  PREPARATION DES VARIABLES A ENVOYER SUR LES CELLULES
!===============================================================================

  ityvar = 1

! --- Informations géométriques de localisation

  call npdcpl(numcpl, ncedis, nfbdis)
  !==========

  ! Allocate temporary arrays for geometric quantities
  allocate(locpts(ncedis))
  allocate(coopts(3,ncedis), djppts(3,ncedis), dofpts(3,ncedis))
  allocate(pndpts(ncedis))

  ! Allocate temporary arrays for variables exchange
  allocate(rvdis(f_dim, ncedis))
  allocate(rvcel(f_dim, ncecpl))

  call coocpl &
  !==========
( numcpl , ncedis , ityvar , &
  ityloc , locpts , coopts , &
  djppts , dofpts , pndpts )

  if (ityloc.eq.2) then
    write(nfecra,1000)
    call csexit(1)
    !==========
  endif

  ! We check that there exists some entities to be coupled at least on one
  ! rank
  ! otherwise there is no need to call gradients for instance
  ncecpg = ncecpl
  ncedig = ncedis
  if (irangp.ge.0) then
    call parcpt(ncecpg)
    call parcpt(ncedig)
  endif

  ! --- Prepare variables to be transfer

  if (ncedig.gt.0) then

    call cscpce &
  ( ncedis , f_id   , f_dim ,                                     &
    locpts ,                                                      &
    coopts , rvdis  )

  endif

  ! Free memory
  deallocate(locpts)
  deallocate(coopts, djppts, dofpts)
  deallocate(pndpts)

  ! Cet appel est symétrique, donc on teste sur NCEDIG et NCECPG
  ! (rien a envoyer, rien a recevoir)
  if (ncedig.gt.0.or.ncecpg.gt.0) then

    ! for vectorial exchange
    call varcpl &
  ( numcpl , ncedis , ncecpl , ityvar , f_dim  ,                  &
    rvdis  ,                                                      &
    rvcel  )

  endif

  ! Free memory
  deallocate(rvdis)

!===============================================================================
! 3.  TRADUCTION DU COUPLAGE EN TERME DE TERMES SOURCES
!===============================================================================

  if (ncecpg.gt.0) then

    call csc2ts(ncecpl, lcecpl, f_id, f_dim, rvcel, crvexp)

  endif

  ! Free memory
  deallocate(rvcel)
  deallocate(lcecpl)

enddo
!     Fin de la boucle sur les couplages


!--------
! FORMATS
!--------
 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :                                             ',/,&
'@    =========                                               ',/,&
'@    LE COUPLAGE VIA LES FACES EN TANT QU''ELEMENTS          ',/,&
'@    SUPPORTS N''EST PAS ENCORE GERE PAR LE NOYAU.           ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
!----
! FIN
!----

return
end subroutine
