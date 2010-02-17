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

subroutine memt1d &
!================

 ( idbia0 , idbra0 ,                                              &
   nfabor , ifnia1 , ifnra1 , ifnia2 , ifnra2 ,                   &
   ifinia , ifinra , ia     , ra     )

!===============================================================================

!  FONCTION
!  --------

!  GESTION MEMOIRE POUR LE MODULE THERMIQUE 1D EN PAROI

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0/idbra0    ! e  ! <-- ! pointeur de la premiere cas libre des          !
!                  !    !     !  tableaux ia/ra                                !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! ifnia1           ! e  ! --> ! pointeur de la premiere cas libre              !
!                  !    !     !  dans ia apres liberation tout sauf            !
!                  !    !     !                               ifpt1d           !
! ifnra1         e  ! --> ! pointeur de la premiere cas libre              !
!                  !    !     !  dans ra apres liberation tout sauf            !
!                  !    !     !                               tppt1d           !
! ifnia2           ! e  ! --> ! pointeur de la premiere cas libre              !
!                  !    !     !  dans ia apres liberation de nppt1d            !
! ifnra2         e  ! --> ! pointeur de la premiere cas libre              !
!                  !    !     !  dans ra apres liberation de rgpt1d            !
!                  !    !     !  et eppt1d                                     !
! ifinia           ! i  ! --> ! number of first free position in ia (at exit)  !
! ifinra           ! i  ! --> ! number of first free position in ra (at exit)  !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "cstnum.h"
include "optcal.h"
include "numvar.h"
include "entsor.h"
include "pointe.h"
include "parall.h"
!===============================================================================

! Arguments
integer          idbia0 ,idbra0
integer          nfabor
integer          ifinia, ifinra, ifnia1, ifnra1, ifnia2, ifnra2
integer          ia(*)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          iok1, ifac


!===============================================================================

!---> INITIALISATION

idebia = idbia0
idebra = idbra0

!---> VERIFICATION DES DIMENSIONS

iok1 = 0
if(nfpt1d.gt.nfabor .or. nfpt1d.lt.0) then
  write(nfecra,1000) nfpt1d
  iok1 = 1
endif
if(iok1.ne.0) then
  call csexit (1)
endif

!---> CALCUL DU NOMBRE DE FACES DE PAROI AVEC MODULE 1D TOTAL

nfpt1t = nfpt1d

if (irangp.ge.0) then
  call parcpt(nfpt1t)
endif

!---> QUELQUES MESSAGES

if(nfpt1t.eq.0) then
  write(nfecra,2000) nfpt1t
else
  write(nfecra,2001) nfpt1t, nfpt1d
endif
write(nfecra,3000)

!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA

ifinia = idebia
ifinra = idebra


iifpt1 = ifinia
ifnia1 = iifpt1 + nfpt1d
iiclt1 = ifnia1
ifnia2 = iiclt1 + nfpt1d
inppt1 = ifnia2
ifinia = inppt1 + nfpt1d

itppt1 = ifinra
ifnra1 = itppt1 + nfpt1d
itept1 = ifnra1
ihept1 = itept1 + nfpt1d
ifept1 = ihept1 + nfpt1d
ixlmt1 = ifept1 + nfpt1d
ircpt1 = ixlmt1 + nfpt1d
idtpt1 = ircpt1 + nfpt1d
ifnra2 = idtpt1 + nfpt1d
ieppt1 = ifnra2
irgpt1 = ieppt1 + nfpt1d
ifinra = irgpt1 + nfpt1d

!---> VERIFICATION

CALL IASIZE('MEMT1D',IFINIA)
!     ==========

CALL RASIZE('MEMT1D',IFINRA)
!     ==========

!---> INITIALISATION DES TABLEAUX
!     a des valeurs sortant en erreur dans vert1d
!     sauf pour les variables de conditions aux limites
!     qui sont initialisees a des valeurs standard
!     (flux nul, Timpose=0, coef d'echange infini)
do ifac = 1, nfpt1d
  ia(iifpt1+ifac-1) = -999
  ia(inppt1+ifac-1) = -999
  ia(iiclt1+ifac-1) = 3
  ra(ieppt1+ifac-1) = -999.d0
  ra(irgpt1+ifac-1) = -999.d0
  ra(itppt1+ifac-1) = 0.d0
  ra(itept1+ifac-1) = 0.d0
  ra(ihept1+ifac-1) = rinfin
  ra(ifept1+ifac-1) = 0.d0
  ra(ixlmt1+ifac-1) = -999.d0
  ra(ircpt1+ifac-1) = -999.d0
  ra(idtpt1+ifac-1) = -999.d0
enddo

!---> FORMATS

#if defined(_CS_LANG_FR)

 1000 format(/,' SORTIE DANS MEMT1D CAR LA DIMENSIONS DU TABLEAU ',/,   &
         '   RELATIF AUX FACES DE BORD COUPLEES AU MODULE ',/,    &
         '   THERMIQUE 1D DE PAROI EST INCORRECTE ',/,      &
         '   NFPT1D = ',I10)

 2000 format(                                                           &
    /,'TTES PHASES  : MODULE THERMIQUE 1D DE PAROI NON ACTIVE ',/,&
      '                 NFPT1D = ',I10,/)
 2001 format(                                                           &
    /,'TTES PHASES  : MODULE THERMIQUE 1D EN PAROI ACTIVE     ',/,&
      '   SUR UN TOTAL DE ',I10,' FACES DE BORD',/,         &
      '   (',I10,' FACES DE BORD EN LOCAL)',/)

 3000 format(                                                           &
'-------------------------------------------------------------',/)

#else

 1000 format(/,' ABORT IN MEMT1D BECAUSE THE DIMENSION OF THE ARRAY ',/,&
         '   RELATIVE TO THE COUPLED FACES OF THE 1D-WALL ',/,    &
         '   THERMAL MODULE IS INCORRECT ',/,               &
         '   NFPT1D = ',I10)

 2000 format(                                                           &
 /,'ALL PHASES  : 1D-WALL THERMAL MODULE NOT ACTIVATED ',/, &
   '                 NFPT1D = ',I10,/)
 2001 format(                                                           &
 /,'ALL PHASES  : 1D-WALL THERMAL MODULE ACTIVATED ',/,     &
   '   ON A TOTAL OF ',I10,' BOUNDARY FACES',/,             &
   '   (',I10,' LOCAL BOUNDARY FACES)',/)

 3000 format(                                                           &
'-------------------------------------------------------------',/)

#endif

return
end subroutine
