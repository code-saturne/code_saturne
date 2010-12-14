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

subroutine cou1do &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncp    , nfpt1d ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ientha , ifpt1d , iclt1d ,                                     &
   idevel , ituser , ia     ,                                     &
   tppt1d , tept1d , hept1d , fept1d ,                            &
   xlmbt1 , rcpt1d , dtpt1d , dt     , rtpa   ,                   &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  ,                                              &
   cpcst  , cp     , hbord  , tbord  ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ---------

! ECRITURE DE DONNEES RELATIVES A UN COUPLAGE AVEC SYRTHES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncp              ! e  ! <-- ! dimension de cp (ncelet ou 1)                  !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nfpt1d           ! e  ! <-- ! nombre de faces avec module therm 1d           !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ientha           ! e  ! <-- ! 1 si tparoi est une enthalpie                  !
! ifpt1d           ! te ! <-- ! numero de la face en traitement                !
!                  !    !     ! thermique en paroi                             !
! iclt1d           ! te ! <-- ! type de condition limite                       !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! cpcst            ! r  ! <-- ! chaleur specifique si constante                !
! cp(ncp)          ! tr ! <-- ! chaleur specifique si variable                 !
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
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
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
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments
integer          idbia0 , idbra0
integer          nfpt1d
integer          nvar   , nscal  , nphas  , ncp
integer          nideve , nrdeve , nituse , nrtuse

integer          ifpt1d(nfpt1d), iclt1d(nfpt1d)
integer          ientha
integer          idevel(nideve), ituser(nituse), ia(*)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision hbord(nfabor),tbord(nfabor)
double precision cpcst, cp(ncp)
double precision tppt1d(nfpt1d)
double precision tept1d(nfpt1d), hept1d(nfpt1d), fept1d(nfpt1d)
double precision xlmbt1(nfpt1d), rcpt1d(nfpt1d), dtpt1d(nfpt1d)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

!     VARIABLES LOCALES

integer          idebia, idebra, mode
integer          iphas  , iappel
integer          ifac, iel , ii
integer          maxelt, ils, idbia1
double precision enthal, temper

!===============================================================================

idebia = idbia0
idebra = idbra0

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
      if(ncp.eq.ncelet) then
         hbord(ifac) = hbord(ifac)*cp(iel)
      else
         hbord(ifac) = hbord(ifac)*cpcst
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
iphas = 1
iappel = 3

maxelt = max(ncelet,nfac,nfabor)
ils    = idebia
idbia1 = ils + maxelt
call iasize('cou1do',idbia1)

call  uspt1d                                                      &
!===========
 ( idbia1 , idebra ,                                              &
   nvar   , nscal  , nphas  , nfpt1d , iphas  , iappel ,          &
   nideve , nrdeve , nituse , nrtuse ,                            &
   maxelt , ia(ils),                                              &
   ifpt1d , ia(idebia), iclt1d ,                                  &
   idevel , ituser , ia     ,                                     &
   tppt1d , ra(idebra), ra(idebra),                               &
   tept1d , hept1d , fept1d ,                                     &
   xlmbt1 , rcpt1d , dtpt1d ,                                     &
   dt     , rtpa   ,                                              &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  ,                                              &
   rdevel , rtuser , ra     )

iappel = 3
call vert1d                                                       &
!==========
 (idebia     , idebra     ,                                       &
  nfabor     , nfpt1d     , iappel    ,                           &
  ifpt1d     , ia(idebia) , iclt1d    , ia     ,                  &
  ra(idebra) , ra(idebra) ,                                       &
  xlmbt1     , rcpt1d     , dtpt1d    , ra      )

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
