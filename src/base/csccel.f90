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

subroutine csccel &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                                                                 &
   ivar   ,                                                       &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   crvexp , crvimp ,                                              &
   dam    , xam    ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ivar             ! i  ! <-- ! variable number                                !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant            prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! crvexp(ncelet    ! tr ! --> ! tableau de travail pour part explicit          !
! crvimp(ncelet    ! tr ! --> ! tableau de travail pour part implicit          !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
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

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ivar

integer          ia(*)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision crvexp(ncelet), crvimp(ncelet)
double precision dam(ncelet ),xam(nfac ,2)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision ra(*)

! Local variables

integer          idebia , idebra , ifinia , ifinra
integer          numcpl , iun
integer          ncesup , nfbsup
integer          ncecpl , nfbcpl , ncencp , nfbncp
integer          ncedis , nfbdis
integer          ncecpg , ncedig
integer          ilfbsu , ilcesu
integer          ilcecp , ilfbcp , ilcenc , ilfbnc
integer          ilocpt , icoopt , idjppt , idofpt , ipndpt
integer          irvdis , irvcel , ipndcp , idofcp
integer          ityloc , ityvar

!====================================================================================

idebia = idbia0
idebra = idbra0

do numcpl = 1, nbrcpl

!======================================================================================
! 1.  DEFINITION DE CHAQUE COUPLAGE
!======================================================================================

  call nbecpl                                                     &
  !==========
 ( numcpl ,                                                       &
   ncesup , nfbsup ,                                              &
   ncecpl , nfbcpl , ncencp , nfbncp )

  call memcs1                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ncesup , nfbsup , ncecpl , nfbcpl , ncencp , nfbncp ,          &
   ilcesu , ilfbsu , ilcecp , ilfbcp , ilcenc , ilfbnc ,          &
   ifinia , ifinra )

  call lelcpl                                                     &
  !==========
 ( numcpl ,                                                       &
   ncecpl , nfbcpl ,                                              &
   ia(ilcecp) , ia(ilfbcp) )


!====================================================================================
! 2.  PREPARATION DES VARIABLES A ENVOYER SUR LES CELLULES
!====================================================================================

  ityvar = 1

! --- Informations géométriques de localisation

  call npdcpl(numcpl, ncedis, nfbdis)
  !==========

!       On n'échangera qu'une composante de la vitesse à la fois
  iun = 1

  call memcs2                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ncecpl , ncedis , iun    ,                                     &
   irvcel , ipndcp , idofcp ,                                     &
   irvdis , ilocpt , icoopt , idjppt , idofpt , ipndpt ,          &
   ifinia , ifinra )

  call coocpl(numcpl, ncedis, ityvar,                             &
  !==========
              ityloc, ia(ilocpt), ra(icoopt),                     &
              ra(idjppt), ra(idofpt), ra(ipndpt) )

  if (ityloc.eq.2) then
    write(nfecra,1000)
    call csexit(1)
    !==========
  endif

!       On vérifie qu'il faut bien échanger quelque chose
!       de manière globale (à cause des appels à GRDCEL notamment)
  ncecpg = ncecpl
  ncedig = ncedis
  if (irangp.ge.0) then
    call parcpt(ncecpg)
    !==========
    call parcpt(ncedig)
    !==========
  endif


! --- Transfert des variables proprement dit.

  if (ncedig.gt.0) then

    call cscpce                                                   &
    !==========
  ( ifinia , ifinra ,                                             &
    nvar   , nscal  ,                                             &
    ncedis , ityloc ,                                             &
    ivar   ,                                                      &
    ia(ilocpt) ,                                                  &
    ia     ,                                    &
    dt     , rtpa   , propce , propfa , propfb ,                  &
    coefa  , coefb  ,                                             &
    w1     , w2     , w3     , w4     , w5     , w6     ,         &
    ra(icoopt)      , ra(irvdis)      ,                           &
    ra     )

  endif

!       Cet appel est symétrique, donc on teste sur NCEDIG et NCECPG
!       (rien a envoyer, rien a recevoir)
  if (ncedig.gt.0.or.ncecpg.gt.0) then

    call varcpl                                                   &
    !==========
  ( numcpl , ncedis , ncecpl , ityvar ,                           &
    ra(irvdis) ,                                                  &
    ra(irvcel) )

  endif


!====================================================================================
! 3.  TRADUCTION DU COUPLAGE EN TERME DE TERMES SOURCES
!====================================================================================

  if (ncecpg.gt.0) then

    call csc2ts                                                   &
    !==========
  ( ifinia , ifinra ,                                             &
    nvar   , nscal  ,                                             &
    ncecpl ,                                                      &
    ivar   ,                                                      &
    ia(ilcecp)      ,                                             &
    ia     ,                                    &
    dt     , rtpa   , propce , propfa , propfb ,                  &
    coefa  , coefb  ,                                             &
    crvexp , crvimp ,                                             &
!         ------   ------
    dam    , xam    ,                                             &
    ra(irvcel)      ,                                             &
    w1     , w2     , w3     , w4     , w5     , w6     ,         &
    ra     )

  endif

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
