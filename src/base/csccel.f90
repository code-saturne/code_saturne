!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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
!> \brief Exchange of coupling variables between two instances of code_saturne
!> thanks to cells.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     f_id          field id
!> \param[in,out] crvexp        working table for explicit part
!> \param[in,out] crvimp        working table for implicit part
!______________________________________________________________________________

subroutine csccel &
 ( f_id   ,                                                       &
   crvexp , crvimp )                                              &
   bind(C, name='cs_sat_coupling_exchange_at_cells')

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

integer(c_int), value :: f_id
real(c_double), dimension(*) :: crvexp
real(c_double), dimension(*) :: crvimp

! Local variables

integer          f_dim, f_dim2
integer          reverse
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
double precision, allocatable, dimension(:,:) :: cw1_dis, cw1_loc
double precision, allocatable, dimension(:,:,:) :: cw2_dis
double precision, allocatable, dimension(:,:,:) :: cw2_loc

procedure() :: csexit, nbecpl, lelcpl, npdcpl, coocpl, cscpce, varcpl, csc2ts

!===============================================================================

! cannot be called for non variables ...
if (f_id.lt.0) then
  call csexit(1)
endif

!get the dimension of the variable
call field_get_dim(f_id, f_dim)

do numcpl = 1, nbrcpl

!===============================================================================
! 1.  DEFINITION DE CHAQUE COUPLAGE
!===============================================================================

  call nbecpl &
 ( numcpl , reverse,                                              &
   ncesup , nfbsup ,                                              &
   ncecpl , nfbcpl , ncencp , nfbncp )

  ! Allocate temporary arrays for coupling information
  allocate(lcecpl(ncecpl))
  allocate(lfbcpl(nfbcpl))

  call lelcpl                                                     &
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

  ! Allocate temporary arrays for geometric quantities
  allocate(locpts(ncedis))
  allocate(coopts(3,ncedis), djppts(3,ncedis), dofpts(3,ncedis))
  allocate(pndpts(ncedis))

  allocate(cw1_dis(f_dim, ncedis))
  allocate(cw2_dis(f_dim,f_dim, ncedis))
  allocate(cw1_loc(f_dim, ncecpl))
  allocate(cw2_loc(f_dim, f_dim, ncecpl))

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

  ! Prepare variables to be transfer

  if (ncedig.gt.0.and.reverse.eq.0) then

    call cscpce &
  ( ncedis , f_id   , f_dim , &
    reverse, locpts ,         &
    coopts , cw1_dis  , cw2_dis )

  endif

  if (ncecpg.gt.0.and.reverse.eq.1) then

    call cscpce &
  ( ncecpl , f_id   , f_dim , &
    reverse, lcecpl ,         &
    coopts , cw1_loc  , cw2_loc)

  endif


  ! Free memory
  deallocate(coopts, djppts, dofpts)
  deallocate(pndpts)

  ! Symmetric call so test on both ncedig and ncecpg (otherwise nothing to do)
  if (ncedig.gt.0.or.ncecpg.gt.0) then
    call varcpl &
  ( numcpl , ncedis , ncecpl , ityvar , f_dim  , &
    reverse, cw1_dis  ,                          &
    cw1_loc )

    ! Second call for the implicit part in reverse mode
    if (reverse.eq.1) then
      f_dim2 = f_dim**2
      call varcpl &
    ( numcpl , ncedis , ncecpl , ityvar , f_dim2  , &
      reverse, cw2_dis,                             &
      cw2_loc )

    endif

  endif

  ! Free memory

!===============================================================================
! 3.  TRADUCTION DU COUPLAGE EN TERME DE TERMES SOURCES
!===============================================================================
  ! Compute source term (n_cells_loc contribution)
  if (ncecpg.gt.0.and.reverse.eq.0) then

    call csc2ts(ncecpl , lcecpl , f_id , f_dim , reverse,  &
                cw1_loc, cw2_loc,                          &
                crvexp , crvimp )

  endif

  ! Compute source term in reverse mode (n_cells_dist contributions)
  if (ncedig.gt.0.and.reverse.eq.1) then

    call csc2ts(ncedis , locpts , f_id , f_dim , reverse ,  &
                cw1_dis, cw2_dis,                          &
                crvexp , crvimp )

  endif

  ! Free memory
  deallocate(cw1_dis)
  deallocate(cw1_loc)
  deallocate(cw2_loc)
  deallocate(cw2_dis)
  deallocate(lcecpl)
  deallocate(locpts)
enddo
! End of loop over the coupling

!--------
! Formats
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
