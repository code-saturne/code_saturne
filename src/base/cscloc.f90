!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

subroutine cscloc &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!   LOCALISATION DES INTERFACES DE COUPLAGE (VIA FVM)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "period.h"
include "cplsat.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia , idebra , ifinia , ifinra
integer          maxelt , ils    , ifnia1
integer          numcpl
integer          ncesup , nfbsup
integer          ncecpl , nfbcpl , ncencp , nfbncp
integer          ilfbsu , ilcesu
integer          ilcecp , ilfbcp , ilcenc , ilfbnc

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

ipass  = ipass + 1

idebia = idbia0
idebra = idbra0

!     On surdimensionne les tableaux à NCEL et NFABOR
ncesup = ncel
nfbsup = nfabor
ncecpl = ncel
nfbcpl = nfabor
ncencp = 0
nfbncp = 0

call memcs1                                                       &
!==========
 ( idebia , idebra ,                                              &
   ncesup , nfbsup , ncecpl , nfbcpl , ncencp , nfbncp ,          &
   ilcesu , ilfbsu , ilcecp , ilfbcp , ilcenc , ilfbnc ,          &
   ifinia , ifinra )

do numcpl = 1, nbrcpl

!       On localise au premier passage ou si l'un des maillages
!         du couplage est mobile ou avec ALE (cf CSCINI).
  if (ipass.eq.1.or.imajcp(numcpl).eq.1) then

    maxelt = max(ncelet, nfac, nfabor)

    ils    = ifinia
    ifnia1 = ils + maxelt
    CALL IASIZE('CSCLOC',IFNIA1)

    call ussatc                                                   &
    !==========
  ( ifnia1 , ifinra , numcpl ,                                    &
    ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml ,&
    nnod   , lndfac , lndfbr , ncelbr ,                           &
    nituse , nrtuse ,                                             &
    ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils),&
    ipnfac , nodfac , ipnfbr , nodfbr ,                           &
    ia(ilcesu) , ia(ilfbsu) ,                                     &
    ia(ilcecp) , ia(ilfbcp) ,                                     &
    ncesup , nfbsup , ncecpl , nfbcpl ,                           &
    xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod )

    if (nfbsup.gt.0) then
      write(nfecra,1000)
      call csexit(1)
      !==========
    endif

!         Localisation proprement dite
    call  defcpl                                                  &
!         ============
  ( numcpl ,                                                      &
    ncesup , nfbsup , ncecpl , nfbcpl ,                           &
    ia(ilcesu)      , ia(ilfbsu)      ,                           &
    ia(ilcecp)      , ia(ilfbcp)      )

  endif

enddo


!--------
! FORMAT
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
