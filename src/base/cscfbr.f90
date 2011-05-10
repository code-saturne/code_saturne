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

subroutine cscfbr &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                                                                 &
   icodcl , itrifb , itypfb ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
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

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,ndim)
double precision ra(*)

! Local variables

integer          idebia , idebra , ifinia , ifinra
integer          numcpl , ivarcp
integer          ncesup , nfbsup
integer          ncecpl , nfbcpl , ncencp , nfbncp
integer          ncedis , nfbdis
integer          nfbcpg , nfbdig
integer          ilcesu , ilfbsu
integer          ilcecp , ilfbcp , ilcenc , ilfbnc
integer          ilocpt , icoopt , idjppt , ipndpt , idofpt
integer          irvdis , irvfbr , ipndcp , idofcp
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

!       Liste des cellules et faces de bord localisées
  call lelcpl                                                     &
  !==========
 ( numcpl ,                                                       &
   ncecpl , nfbcpl ,                                              &
   ia(ilcecp) , ia(ilfbcp) )

!       Liste des cellules et faces de bord non localisées
  call lencpl                                                     &
  !==========
 ( numcpl ,                                                       &
   ncencp , nfbncp ,                                              &
   ia(ilcenc) , ia(ilfbnc) )


!====================================================================================
! 2.  PREPARATION DES VARIABLES A ENVOYER SUR LES FACES DE BORD
!====================================================================================

  ityvar = 2

! --- Informations géométriques de localisation

  call npdcpl(numcpl, ncedis, nfbdis)
  !==========

  call memcs2                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   nfbcpl , nfbdis , nvarto(numcpl) ,                             &
   irvfbr , ipndcp , idofcp ,                                     &
   irvdis , ilocpt , icoopt , idjppt , idofpt , ipndpt ,          &
   ifinia , ifinra )

  call coocpl(numcpl, nfbdis, ityvar,                             &
  !==========
              ityloc, ia(ilocpt), ra(icoopt),                     &
              ra(idjppt), ra(idofpt), ra(ipndpt))

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
  ( ifinia , ifinra ,                                             &
    nvar   , nscal  ,                                             &
    nfbdis , ityloc , nvarcp(numcpl) , numcpl ,                   &
    nvarto(numcpl) ,                                              &
    ia(ilocpt) ,                                                  &
    ia     ,                                                      &
    dt     , rtp    , rtpa   , propce , propfa , propfb ,         &
    coefa  , coefb  ,                                             &
    w1     , w2     , w3     , w4     , w5     , w6     ,         &
    ra(icoopt)      , ra(idjppt)      , ra(ipndpt)      ,         &
    ra(irvdis)      , ra(idofpt)      ,                           &
    ra     )

  endif

!       Cet appel est symétrique, donc on teste sur NFBDIG et NFBCPG
!       (rien a envoyer, rien a recevoir)
  if (nfbdig.gt.0.or.nfbcpg.gt.0) then

    do ivarcp = 1, nvarto(numcpl)

      call varcpl                                                 &
      !==========
    ( numcpl , nfbdis , nfbcpl , ityvar ,                         &
      ra(irvdis + (ivarcp-1)*nfbdis) ,                            &
      ra(irvfbr + (ivarcp-1)*nfbcpl) )

    enddo

  endif


!====================================================================================
! 3.  TRADUCTION DU COUPLAGE EN TERME DE CONDITIONS AUX LIMITES
!====================================================================================

  if (nfbcpg.gt.0) then

    call pndcpl                                                   &
    !==========
  ( numcpl , nfbcpl , ityvar , ra(ipndcp) , ra(idofcp) )

    call csc2cl                                                   &
    !==========
  ( ifinia , ifinra ,                                             &
    nvar   , nscal  ,                                             &
    nvarcp(numcpl), nvarto(numcpl) , nfbcpl , nfbncp ,            &
    icodcl , itrifb , itypfb ,                                    &
    ia(ilfbcp) , ia(ilfbnc) ,                                     &
    ia     ,                                                      &
    dt     , rtp    , rtpa   , propce , propfa , propfb ,         &
    coefa  , coefb  , rcodcl ,                                    &
    w1     , w2     , w3     , w4     , w5     , w6     , coefu  ,&
    ra(irvfbr)      , ra(ipndcp)      , ra(idofcp)      ,         &
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
