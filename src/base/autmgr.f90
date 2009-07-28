!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine autmgr &
!================

    ( igr    , isym   , iagmax , nagmax ,                         &
      ncelf  , ncelfe , nfacf  , iwarnp ,                         &
      ifacef ,                                                    &
      daf    , xaf    , surfaf , volumf , xyzfin ,                &
      irscel ,                                                    &
      indic  , inombr , irsfac , indicf , w1     , w2 )

!===============================================================================
! FONCTION :
! ----------

!  MULTIGRILLE ALGEBRIQUE :
!  CONSTRUCTION D'UN NIVEAU DE MAILLAGE GROSSIER A PARTIR
!  DU NIVEAU SUPERIEUR SUIVANT CRITERE AUTOMATIQUE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nfecra           ! e  ! <-- ! numero du fichier d'impressions                !
! isym             ! e  ! <-- ! indicateur = 1 matrice sym                     !
!                  !    !     !            = 2 matrice non sym                 !
! igr              ! e  ! <-- ! niveau du maillage grossier                    !
! ncelf            ! e  ! <-- ! nombre d'elements maillage fin                 !
! ncelfe           ! e  ! <-- ! nombre d'elements etendus fin                  !
! nfacf            ! e  ! <-- ! nombre de faces internes maill. fin            !
! iwarnp           ! e  ! <-- ! niveau d'impression                            !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! ifacef           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfacf)       !    !     !  du maillage fin                               !
! daf(ncelfe)      ! tr ! <-- ! diagonale matrice maillage fin                 !
! xaf              ! tr ! <-- ! extradiagonale matrice maillage fin            !
! (nfacf, isym)    !    !     !                                                !
! surfaf           ! tr ! <-- ! surfaces faces internes maillage fin           !
! (3, nfacf)       !    !     !                                                !
! volumf           ! tr ! <-- ! volumes des cellules du maillage fin           !
! (ncelfe)         !    !     !                                                !
! xyzfin           ! tr ! <-- ! centres des cellules du maillage fin           !
! (3, ncelfe)      !    !     !                                                !
! irscel           ! te ! --> ! cellule fine -> cellule grossiere              !
!  (ncelfe)        !    !     !                                                !
! indic(ncelfe)    ! te ! --- ! tableau de travail                             !
! inombr           ! te ! --- ! tableau de travail                             !
!  (ncelfe)        !    !     !                                                !
! irsfac           ! te ! --- ! face fine -> face grossiere                    !
!  (nfacf)         !    !     !  (tableau de travail)                          !
! indicf(nfacf)    ! te ! --- ! indicateur de regroupement des faces           !
! icelfa           ! te ! --- ! connectivite cellules->faces mailla-           !
!  (2*nfacf)       !    !     ! ge fin                                         !
! icelce           ! te ! --- ! connectivite cellules->cellules                !
!  (2*nfacf)       !    !     ! voisines du maillage fin                       !
! rw(ncelf)        ! tr ! --- ! tableau de travail                             !
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
include "entsor.h"
include "optcal.h"
include "parall.h"

!===============================================================================


!===============================================================================

! Arguments

integer          igr, isym, iagmax, nagmax
integer          ncelf, ncelfe, nfacf
integer          iwarnp

integer          ifacef(2, nfacf)
integer          irscel(ncelfe)
integer          indic(ncelfe), inombr(ncelfe)
integer          indicf(nfacf), irsfac(nfacf)

double precision daf(ncelfe), xaf(nfacf,isym)
double precision surfaf(3,nfacf), volumf(ncelfe)
double precision xyzfin(3,ncelfe)
double precision w1(ncelfe), w2(ncelfe)

! Local variables

double precision critr, epslon

integer          ncelg, icel, ifac , ifac1, icelg
integer          nfacn,nfacnr,npass,npasmx
integer          inditt, noaglo, ngros, incvoi
integer          ihist(10)
integer          i, j, imin, imax

!===============================================================================

! Default parameters

epslon = +1.d-6
ngros  = 8
npasmx = 10
incvoi = 1

!===============================================================================
! 1. Initialization
!===============================================================================

do icel = 1, ncelfe
  indic(icel) = -1
  irscel(icel) = 0
  inombr(icel) = 1
enddo
do ifac = 1, nfacf
  indicf(ifac) = ifac
  irsfac(ifac) = ifac
enddo

! Compute cardinality (number of neighbors for each cell -1)

do ifac = 1, nfacf
  i = ifacef(1,ifac)
  j = ifacef(2,ifac)
  indic(i) = indic(i) + 1
  indic(j) = indic(j) + 1
enddo

ncelg  = 0
nfacnr = nfacf
npass  = 0
noaglo = ncelf

! Passes

 100  continue

npass = npass+1
nfacn = nfacnr
iagmax= iagmax +1
iagmax= min(iagmax, nagmax)

do ifac=1,nfacn
  irsfac(ifac) = indicf(ifac)
  indicf(ifac) = 0
enddo
if (nfacn .lt. nfacf) then
  do ifac = nfacn+1, nfacf
    indicf(ifac) = 0
    irsfac(ifac) = 0
  enddo
endif

if (iwarnp .gt. 3) then
    write(nfecra,2001) npass, nfacnr, noaglo
endif

  ! Increment number of neighbors

do icel = 1, ncelf
  indic(icel) = indic(icel) + incvoi
enddo

  ! Initialize non-eliminated faces

nfacnr = 0

  ! Loop on non-eliminated faces

do ifac1 = 1, nfacn

  ifac = irsfac(ifac1)
  i = ifacef(1,ifac)
  j = ifacef(2,ifac)

    ! Exclude faces on parallel or periodic boundary, so as not to
    ! coarsen the grid across those boundaries (which would change
    ! the communication pattern and require a more complex algorithm).

  if (i.le.ncelf .and. j.le.ncelf) then

    inditt = 0
    critr  = (daf(i)/indic(i))*(daf(j)/indic(j))                  &
             /( xaf(ifac,1)*xaf(ifac,isym))

    if (       critr.lt.(1.d0-epslon)                             &
         .and. irscel(i)*irscel(j).le.0) then

      if (irscel(i).gt.0 .and. irscel(j).le.0) then
        if(inombr(irscel(i)) .le. iagmax) then
          irscel(j) = irscel(i)
          inombr(irscel(i)) = inombr(irscel(i)) +1
          inditt = inditt +1
        endif
      else if (irscel(i).le.0 .and. irscel(j).gt.0) then
        if (inombr(irscel(j)).le.iagmax) then
          irscel(i) = irscel(j)
          inombr(irscel(j)) = inombr(irscel(j)) + 1
          inditt = inditt +1
        endif
      else if (irscel(i).le.0 .and. irscel(j).le.0) then
        ncelg = ncelg+1
        irscel(i) = ncelg
        irscel(j) = ncelg
        inombr(ncelg) = inombr(ncelg) +1
        inditt = inditt +1
      endif

    endif

    if (inditt.ne.0 .and.inditt.ne.1) then
      write(nfecra,*) ' Bug in autmgr, contact support.'
      call csexit(1)
    endif

    if (inditt.eq.0 .and. irscel(i)*irscel(j).le.0) then
      nfacnr = nfacnr +1
      indicf(nfacnr) = ifac
    endif

  endif

enddo

  ! Check the number of coarse cells created

noaglo = 0
do icel=1,ncelf
  if (irscel(icel).le.0) noaglo = noaglo+1
enddo

! Loop on passes

if (noaglo.gt.0) then
  if ((ncelg+noaglo)*ngros .ge. ncelf) then
    if (npass.lt.npasmx .and. nfacnr.gt.0) then
      goto 100
    endif
  endif
endif

! Finish assembly

do icel = 1, ncelf
  if (irscel(icel).le.0) then
    ncelg = ncelg+1
    irscel(icel) = ncelg
  endif
enddo

! Various checks

imax = 0
imin = 2*ncelf
do icelg =1,ncelg
  imax = max(imax, inombr(icelg))
  imin = min(imin, inombr(icelg))
enddo

if (irangp .ge. 0) then
  call parcmn(imin)
  call parcmx(imax)
endif

if (iwarnp.gt.3) then

  write(nfecra,2002) imin, imax, nagmax
  write(nfecra,2003)
  noaglo=imax-imin+1
  if (noaglo.gt.0) then
    if (noaglo.gt.10) then
      write(nfecra,*) ' ihist badly dimensioned in autmgr'
      call csexit(1)
    endif
    do i = 1, noaglo
      ihist(i) = 0
    enddo
    do icelg = 1, ncelg
      do i = 1, noaglo
        if (inombr(icelg).eq.(imin+i-1))then
          ihist(i)=ihist(i)+1
        endif
      enddo
    enddo
    if (irangp .ge. 0) then
      call parism(noaglo, ihist)
    endif
    do i = 1, noaglo
      epslon = 100.d0*ihist(i)/ncelg
      write(nfecra,2004) imin+i-1, epslon
    enddo
  endif

endif

do icel = 1, ncelf
  indic(icel) = 0
enddo
do icel = 1, ncelf
  icelg = irscel(icel)
  indic(icelg) = indic(icelg) +1
enddo

i=0
j=2*ncelf
noaglo = 0
do icelg = 1, ncelg
  i = max(i, indic(icelg))
  j = min(j, indic(icelg))
  noaglo = noaglo + indic(icelg)
enddo

if (irangp .ge. 0) then
  call parcmn(j)
  call parcmx(j)
endif

if (iwarnp.gt.3) then
  write(nfecra,2005) j, i
endif

if (noaglo .ne. ncelf) then
  write(nfecra,*) ' Bug in autmgr, contact support.'
  call csexit(1)
endif

!--------
! Formats
!--------

 2001 format(&
  '    autmgr: pass ', i3, 'nfacnr = ', i10, ' noaglo = ', i10)
 2002 format(&
  '    autmgr: inombr min = ', i10, ' max = ', i10, ' target = ', i10)
 2003 format(&
  '      histogram ')
 2004 format(&
  '        regroupment ', i10,' = ', e12.5,' %')
 2005 format(&
  '    autmgr: agglomeration min = ', i10, ' max= ', i10)
!==============================================================================

return
end

