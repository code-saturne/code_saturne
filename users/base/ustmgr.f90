!-------------------------------------------------------------------------------

!VERS


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

subroutine ustmgr &
!================


    ( iappel , igr    , isym   ,                                  &
      ncelf  , ncelfe , nfacf  , iwarnp ,                         &
      iusmgr , niw    , nrw    ,                                  &
      ifacef ,                                                    &
      daf    , xaf    , surfaf , volumf , xyzfin ,                &
      irscel ,                                                    &
      iw     , rw     )

!==========================================================================
! FONCTION :
! ----------

!  MULTIGRILLE ALGEBRIQUE :

!  DETERMINATION DE LA CONNECTIVITE
!  CELLULE FINE -> CELLULE GROSSIERE (IRSCEL)
!  POUR LA CONSTRUCTION D'UN NIVEAU DE MAILLAGE GROSSIER
!  A PARTIR DU NIVEAU SUPERIEUR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! iappel           ! e  ! <-- ! 1 : selection/dimensionnement                  !
!                  !    !     ! 2 : determination de irscel                    !
! igr              ! e  ! <-- ! niveau du maillage grossier                    !
! isym             ! e  ! <-- ! indicateur = 1 matrice sym                     !
!                  !    !     !            = 2 matrice non sym                 !
! ncelf            ! e  ! <-- ! nombre d'elements maillage fin                 !
! ncelfe           ! e  ! <-- ! nombre d'elements etendus fin                  !
! nfacf            ! e  ! <-- ! nombre de faces internes maill. fin            !
! iwarnp           ! e  ! <-- ! niveau d'impression                            !
! iusmgr           ! e  ! --> ! 0 : agglomeration automatique                  !
!                  !    !     ! 1 : on utilise ce sous-programme               !
! niw / nrw        ! e  ! --> ! tailles des tableaux iw / rw                   !
!                  !    !     ! pour iappel = 2                                !
!                  !    !     ! (determinees lorsque iappel = 1)               !
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
! iw(niw)          ! te ! --- ! tableau de travail                             !
! rw(nrw)          ! tr ! --- ! tableau de travail                             !
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

!===============================================================================

! Arguments

integer          iappel, igr, isym
integer          ncelf, ncelfe, nfacf
integer          iwarnp
integer          iusmgr, niw, nrw

integer          ifacef(2, nfacf)
integer          irscel(ncelfe)

double precision daf(ncelfe), xaf(nfacf, isym)
double precision surfaf(3, nfacf), volumf(ncelfe)
double precision xyzfin(3, ncelfe)
integer          iw(niw)
double precision rw(nrw)


! VARIABLES LOCALES

integer          ncelg
integer          ifac
integer          indic, irspr2, ipaimp, iw1
integer          iedir

integer          imp, jmp, kmp, i, j, k, icel
integer          ntest, ip, jp, kp, ijk, iglob, i1, ig, jg
double precision dx, dy, dz, xyzmin(3), xyzmax(3)
double precision xmn, xmx, ymn, ymx, zmn, zmx
double precision xg, yg, zg, rap, rep, epslon
double precision ngros


!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

ncelg = 0
epslon = +1.d-6

!===============================================================================
! 0.  INITIALISATION ET DIMENSIONNEMENT
!===============================================================================

iusmgr = 0

!     DIMENSIONNEMENT MEMOIRE

indic  = 1
irspr2 = indic + ncelfe
ipaimp = irspr2 + ncelfe
niw    = ipaimp + ncelfe

iw1    = 1
nrw    = iw1 + ncelfe

if (iappel.eq.1) return

!===============================================================================
! 1.  CREATION DES CELLULES GROSSIERES : IRSCEL
!===============================================================================

!===============================================================================
!  MAILLAGE DU CARRE CAS 2D

!     BORNE MIN/MAX DU MAILLAGE

do i = 1, 3
  xyzmin(i) = +1.d12
  xyzmax(i) = -1.d12
  do icel = 1, ncelf
    xyzmin(i) = min(xyzfin(i, icel), xyzmin(i))
    xyzmax(i) = max(xyzfin(i, icel), xyzmax(i))
  enddo
  xyzmin(i) = xyzmin(i) -epslon
  xyzmax(i) = xyzmax(i) +epslon
enddo

!     EVALUATION DE IMP,JMP,KMP AVEC LA REGLE MAILLAGE REGULIER
!      IMP/(Xmax-Xmin) ~ JMP/(Ymax-Ymin) ~ KMP/(Zmax-Zmin)
!                       NCELG ~ IMP.JMP.KMP ~ NCEL/NGROS


!     CAS 2D

ngros = 4

rap    = ncelf/ngros
dx = xyzmax(1) -xyzmin(1)
dy = xyzmax(2) -xyzmin(2)
rep = sqrt(rap*dx/dy)
kmp   = 1
imp = int(rep)
jmp = int(rep*dy/dx)
imp = max(imp, 1)
jmp = max(jmp, 1)
ncelg = imp*jmp*kmp

!     INCLUSION DANS LE MAILLAGE STRUCTURÉ ORTHOGONAL

dx = (xyzmax(1) -xyzmin(1))/imp
dy = (xyzmax(2) -xyzmin(2))/jmp
dz = (xyzmax(3) -xyzmin(3))/kmp

do icel = 1, ncelf
  iw(indic+icel-1) = 0
  iw(irspr2+icel-1) = 0
enddo

do icel = 1, ncelf

!       INITIALISATION
  xg = xyzfin(1, icel)
  yg = xyzfin(2, icel)
  zg = xyzfin(3, icel)
  ip = 0
  jp = 0
  kp = 0

  do i = 1, imp
    xmn = xyzmin(1) +(i-1)*dx
    xmx = xyzmin(1) +   i *dx
    if (ip.eq.0) then
      if (xg.ge.xmn .and. xg.le.xmx) then
        ip = i
      endif
    endif
  enddo

  do j = 1, jmp
    ymn = xyzmin(2) +(j-1)*dy
    ymx = xyzmin(2) +   j *dy
    if (jp.eq.0) then
      if (yg.ge.ymn .and. yg.le.ymx) then
        jp = j
      endif
    endif
  enddo

  do k = 1, kmp
    zmn = xyzmin(3) +(k-1)*dz
    zmx = xyzmin(3) +   k *dz
    if (kp.eq.0) then
      if (zg.ge.zmn .and. zg.le.zmx) then
        kp = k
      endif
    endif
  enddo

  ijk = ip*jp*kp

  if (ijk.gt.0) then
    iglob = imp*jmp*(kp-1) +imp*(jp-1)+ip
    iw(irspr2+icel-1) = iglob
    iw(indic+iglob-1) = +1
  else
    WRITE(NFECRA,*) ' pb de programmation dans ustmgr.F  '
    WRITE(NFECRA,*) ' arret dans ustmgr.F '
    WRITE(NFECRA,*) ' IP, JP, KP = ', IP, JP, KP
    call csexit(1)
  endif

enddo


!     DEFINITION DU MAILLAGE PAIR-IMPAIR POUR VISUALISATION

do iglob =1,ncelf
  iw(ipaimp+iglob-1) = 0
enddo

do k = 1, kmp
  do j = 1, jmp
    i1 = mod(k+j,2) +1
    do i = i1, imp, 2
      iglob = imp*jmp*(k-1) +imp*(j-1) +i
      iw(ipaimp+iglob-1) = -1
    enddo
    i1 = mod(k+j+1,2)+1
    do i = i1, imp, 2
      iglob = imp*jmp*(k-1) +imp*(j-1) +i
      iw(ipaimp+iglob-1) = +1
    enddo
  enddo
enddo

!     Compactage du maillage structure si, au moins, une maille
!     grossiere ne contient aucune maille du maillage fin non structure

ntest = 0
do iglob = 1, ncelg
  if (iw(indic+iglob-1).le.0) ntest=ntest+1
enddo

if (ntest.le.0) then
  do icel=1,ncelf
    irscel(icel) = iw(irspr2+icel-1)
  enddo

else

  i = 0
  do iglob = 1, ncelg
    if (iw(indic+iglob-1).gt.0) then
      i = i +1
      iw(indic+iglob-1) = i
    endif
  enddo

!       NOUVEAU NCELG
  ncelg = i

  do icel=1,ncelf
    irscel(icel) = iw(indic+iw(irspr2+icel-1)-1)
  enddo

endif

i = 0
j = ncelf
do icel = 1, ncelf
  i = min(irscel(icel), i)
  j = max(irscel(icel), j)
enddo

iedir = 0
!     IEDIR = 1

if (iedir.gt.0) then

!       EXTRACTION DES POINTS DIRICHLETS (DIAGONALE FORTEMENT DOMINANTE)

  do icel = 1, ncelf
    iw(indic+icel-1) = 0
    iw(irspr2+icel-1) = 0
    iw(iw1+icel-1) = daf(icel)
  enddo

  do ifac = 1, nfacf
    i = ifacef(1,ifac)
    j = ifacef(2,ifac)
    iw(iw1+i-1) = iw(iw1+i-1) +xaf(ifac,1)
    iw(iw1+j-1) = iw(iw1+j-1) +xaf(ifac,isym)
    ig = irscel(i)
    jg = irscel(j)
    iw(indic+ig-1) = iw(indic+ig-1) +1
    iw(indic+jg-1) = iw(indic+jg-1) +1
  enddo

  rep = 0.1d0
  do icel = 1, ncelf
    xmx = iw(iw1+icel-1)/daf(icel)
    ig = irscel(icel)
    if (xmx.gt.rep) then
      if (iw(indic+ig-1).ge.2) then
        ncelg = ncelg+1
        irscel(icel) = ncelg
        iw(indic+ig-1) = iw(indic+ig-1) -1
        iw(indic+ncelg-1) = iw(indic+ncelg-1) -1
      endif
    endif
  enddo

endif


!==============================================================================

!--------
! FORMATS
!--------


!----
! FIN
!----

return
end
