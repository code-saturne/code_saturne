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

subroutine cscfbr &
!================

 ( nscal  ,                                                       &
   icodcl , itypfb ,                                              &
   dt     , rtp    , propce ,                                     &
   coefa  , coefb  , rcodcl )

!===============================================================================
! FONCTION :
! --------

! ECHANGE DES VARIABLES POUR UN COUPLAGE
!   ENTRE DEUX INSTANCES DE CODE_SATURNE VIA LES FACES DE BORD

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nscal            ! i  ! <-- ! total number of scalars                        !
! ivar             ! i  ! <-- ! variable number                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! crvexp(ncelet    ! tr ! --> ! tableau de travail pour part explicit          !
! crvimp(ncelet    ! tr ! --> ! tableau de travail pour part implicit          !
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
use parall
use period
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nscal

integer          icodcl(nfabor,nvarcl)
integer          itypfb(nfabor)

double precision dt(ncelet), rtp(ncelet,*)
double precision propce(ncelet,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvarcl,3)

! Local variables

integer          numcpl , ivarcp
integer          ncesup , nfbsup
integer          ncecpl , nfbcpl , ncencp , nfbncp
integer          ncedis , nfbdis
integer          nfbcpg , nfbdig
integer          ityloc , ityvar
integer          stride

integer, allocatable, dimension(:) :: lcecpl , lfbcpl , lcencp , lfbncp
integer, allocatable, dimension(:) :: locpts

double precision, allocatable, dimension(:,:) :: coopts , djppts , dofpts
double precision, allocatable, dimension(:,:) :: dofcpl
double precision, allocatable, dimension(:) :: pndpts
double precision, allocatable, dimension(:) :: pndcpl
double precision, allocatable, dimension(:,:) :: rvdis , rvfbr

!===============================================================================


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
  allocate(lcecpl(ncecpl), lcencp(ncencp))
  allocate(lfbcpl(nfbcpl), lfbncp(nfbncp))

!       Liste des cellules et faces de bord localisées
  call lelcpl                                                     &
  !==========
 ( numcpl ,                                                       &
   ncecpl , nfbcpl ,                                              &
   lcecpl , lfbcpl )

!       Liste des cellules et faces de bord non localisées
  call lencpl                                                     &
  !==========
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
  !==========

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
  !==========
( numcpl , nfbdis , ityvar , &
  ityloc , locpts , coopts , &
  djppts , dofpts , pndpts )

  if (ityloc.eq.2) then
    write(nfecra,1000)
    call csexit(1)
    !==========
  endif

!       On vérifie qu'il faut bien échanger quelque chose
!       de manière globale (à cause des appels à GRDCEL notamment)
  nfbcpg = nfbcpl
  nfbdig = nfbdis
  if (irangp.ge.0) then
    call parcpt(nfbcpg)
    !==========
    call parcpt(nfbdig)
    !==========
  endif


! --- Transfert des variables proprement dit.

  if (nfbdig.gt.0) then

    call cscpfb                                                   &
    !==========
  ( nscal  ,                                                      &
    nfbdis , numcpl , nvarto(numcpl) ,                            &
    locpts ,                                                      &
    rtp    , propce ,                                             &
    coefa  , coefb  ,                                             &
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
      !==========
    ( numcpl , nfbdis , nfbcpl , ityvar , stride , &
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
    !==========
  ( numcpl , nfbcpl , ityvar , pndcpl , dofcpl )

    call csc2cl &
    !==========
  ( nvarcp(numcpl), nvarto(numcpl) , nfbcpl , nfbncp ,            &
    icodcl , itypfb ,                                             &
    lfbcpl , lfbncp ,                                             &
    dt     , rtp    ,                                             &
    coefa  , coefb  , rcodcl ,                                    &
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
