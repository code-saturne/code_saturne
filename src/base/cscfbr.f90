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

!> \file cscfbr.f90
!> \brief Exchange of variables for coupling two code_saturne intances
!> with boundary faces.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     itypfb        boundary face types
!> \param[in]     dt            time step (per cell)
!______________________________________________________________________________

subroutine cscfbr(itypfb, dt)  &
  bind(C, name='cs_f_cscfbr')

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use cplsat
use mesh
use field
use cs_c_bindings
use dimens, only: nscal

!===============================================================================

implicit none

! Arguments

integer(c_int) ::  itypfb(nfabor)
real(c_double) ::  dt(ncelet)

! Local variables

integer          numcpl , ivarcp
integer          ncesup , nfbsup
integer          ncecpl , nfbcpl , ncencp , nfbncp
integer          ncedis , nfbdis
integer          nfbcpg , nfbdig
integer          ityloc , ityvar
integer          stride
integer          reverse
integer, allocatable, dimension(:) :: lcecpl , lfbcpl , lcencp , lfbncp
integer, allocatable, dimension(:) :: locpts

double precision, allocatable, dimension(:,:) :: coopts , djppts , dofpts
double precision, allocatable, dimension(:,:) :: dofcpl
double precision, allocatable, dimension(:)   :: pndpts
double precision, allocatable, dimension(:)   :: pndcpl
double precision, allocatable, dimension(:,:) :: rvdis , rvfbr

integer, pointer, dimension(:,:) :: icodcl
double precision, pointer, dimension(:,:,:) :: rcodcl

procedure() :: nbecpl, lelcpl, lencpl, npdcpl, coocpl, cscpfb, varcpl
procedure() :: pondcp, csc2cl

!===============================================================================

do numcpl = 1, nbrcpl

!===============================================================================
! 1.  DEFINITION DE CHAQUE COUPLAGE
!===============================================================================

  call nbecpl                                                     &
 ( numcpl , reverse,                                              &
   ncesup , nfbsup ,                                              &
   ncecpl , nfbcpl , ncencp , nfbncp )

  ! Allocate temporary arrays for coupling information
  allocate(lcecpl(ncecpl), lcencp(ncencp))
  allocate(lfbcpl(nfbcpl), lfbncp(nfbncp))

  ! Liste des cellules et faces de bord localisées
  call lelcpl                                                     &
 ( numcpl ,                                                       &
   ncecpl , nfbcpl ,                                              &
   lcecpl , lfbcpl )

  ! Liste des cellules et faces de bord non localisées
  call lencpl                                                     &
 ( numcpl ,                                                       &
   ncencp , nfbncp ,                                              &
   lcencp , lfbncp )

  ! Free memory
  deallocate(lcecpl, lcencp)

!===============================================================================
! 2.  PREPARATION DES VARIABLES A ENVOYER SUR LES FACES DE BORD
!===============================================================================

  ityvar = 2

! --- Informations géométriques de localisation

  call npdcpl(numcpl, ncedis, nfbdis)

  ! Allocate temporary arrays for geometric quantities
  allocate(locpts(nfbdis))
  allocate(coopts(3,nfbdis), djppts(3,nfbdis), dofpts(3,nfbdis))
  allocate(pndpts(nfbdis))
  ! Allocate temporary arrays for variables exchange
  if (nfbdis.gt.0) then
    allocate(rvdis(nfbdis,nvarto(numcpl)))
  else
    allocate(rvdis(1,nvarto(numcpl)))
  endif
  if (nfbcpl.gt.0) then
    allocate(rvfbr(nfbcpl,nvarto(numcpl)))
  else
    allocate(rvfbr(1,nvarto(numcpl)))
  endif

  call coocpl &
( numcpl , nfbdis , ityvar , &
  ityloc , locpts , coopts , &
  djppts , dofpts , pndpts )

  if (ityloc.eq.2) then
    write(nfecra,1000)
    call csexit(1)
  endif

!       On vérifie qu'il faut bien échanger quelque chose
!       de manière globale (à cause des appels à GRDCEL notamment)
  nfbcpg = nfbcpl
  nfbdig = nfbdis
  if (irangp.ge.0) then
    call parcpt(nfbcpg)
    call parcpt(nfbdig)
  endif

  ! Transfert des variables proprement dit.
  if (nfbdig.gt.0) then

    call cscpfb &
  ( nscal  ,                                                      &
    nfbdis , numcpl , nvarto(numcpl) ,                            &
    locpts ,                                                      &
    coopts , djppts , pndpts ,                                    &
    rvdis  , dofpts )

  endif

  ! Free memory
  deallocate(locpts)
  deallocate(coopts, djppts, dofpts)
  deallocate(pndpts)

!       Cet appel est symétrique, donc on teste sur NFBDIG et NFBCPG
!       (rien a envoyer, rien a recevoir)
  if (nfbdig.gt.0.or.nfbcpg.gt.0) then

    do ivarcp = 1, nvarto(numcpl)

      stride = 1

      call varcpl &
    ( numcpl , nfbdis , nfbcpl , ityvar , stride , &
      reverse,                                     &
      rvdis(1, ivarcp) ,                           &
      rvfbr(1, ivarcp) )

    enddo

  endif

  ! Free memory
  deallocate(rvdis)

!===============================================================================
! 3.  TRADUCTION DU COUPLAGE EN TERME DE CONDITIONS AUX LIMITES
!===============================================================================

  if (nfbcpg.gt.0) then

    ! Allocate temporary arrays for geometric quantities
    allocate(dofcpl(3,nfbcpl))
    allocate(pndcpl(nfbcpl))

    call pondcp &
  ( numcpl , nfbcpl , ityvar , pndcpl , dofcpl )

    call field_build_bc_codes_all(icodcl, rcodcl) ! Get map

    call csc2cl &
  ( nvarcp(numcpl), nvarto(numcpl) , nfbcpl , nfbncp ,            &
    icodcl , itypfb ,                                             &
    lfbcpl , lfbncp ,                                             &
    dt     ,                                                      &
    rcodcl ,                                                      &
    rvfbr  , pndcpl , dofcpl )

    ! Free memory
    deallocate(dofcpl, pndcpl)

  endif

  ! Free memory
  deallocate(rvfbr)
  deallocate(lfbcpl, lfbncp)

enddo

!     Fin de la boucle sur les couplages

!--------
! FORMATS
!--------
 1000 format(                                                           &
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

!===============================================================================

!> \brief Initialization for coupling two code_saturne intances
!>        with boundary faces.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     itypfb        boundary face types
!______________________________________________________________________________

subroutine cscfbr_init(itypfb)      &
  bind(C, name='cs_f_cscfbr_init')

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use cplsat
use mesh
use field

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

procedure() :: nbecpl, lelcpl, lencpl, csc2cl, csc2cl_init

! Arguments

integer(c_int) :: itypfb(nfabor)

! Local variables

integer          numcpl
integer          reverse
integer          ncesup , nfbsup
integer          ncecpl , nfbcpl , ncencp , nfbncp

integer, allocatable, dimension(:) :: lcecpl , lfbcpl , lcencp , lfbncp

integer, pointer, dimension(:,:) :: icodcl
double precision, pointer, dimension(:,:,:) :: rcodcl

!===============================================================================


do numcpl = 1, nbrcpl

  ! For each coupling

  call nbecpl(numcpl, reverse,                    &
              ncesup, nfbsup,                     &
              ncecpl, nfbcpl, ncencp, nfbncp)

  ! Allocate temporary arrays for coupling information
  allocate(lcecpl(ncecpl), lcencp(ncencp))
  allocate(lfbcpl(nfbcpl), lfbncp(nfbncp))

  ! Located cells and boundary faces
  call lelcpl(numcpl, ncecpl, nfbcpl, lcecpl, lfbcpl)

  ! Unlocated cells and boundary faces
  call lencpl(numcpl, ncencp, nfbncp, lcencp, lfbncp)

  ! Free memory
  deallocate(lcecpl, lcencp)

!===============================================================================
! 2.  TRADUCTION DU COUPLAGE EN TERME DE CONDITIONS AUX LIMITES
!===============================================================================

  if (nfbcpl.gt.0) then
    call field_build_bc_codes_all(icodcl, rcodcl) ! Get map
    call csc2cl_init(nvarcp(numcpl), nfbcpl, nfbncp,              &
                     icodcl, itypfb, lfbcpl, lfbncp)
  endif

  ! Free memory
  deallocate(lfbcpl, lfbncp)

enddo

return
end subroutine cscfbr_init
