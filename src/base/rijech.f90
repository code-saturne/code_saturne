!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine rijech &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iphas  , ivar   , isou   , ipp    ,                            &
   ia     ,                                                       &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , produc , smbr   ,                            &
   coefax , coefbx ,                                              &
   produk , w2     , w3     , w4     , epsk   , w6     ,          &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! TERMES D'ECHO DE PAROI
!   POUR Rij
! VAR  = R11 R22 R33 R12 R13 R23
! ISOU =  1   2   3   4   5   6

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! iphas            ! i  ! <-- ! phase number                                   !
! ivar             ! i  ! <-- ! variable number                                !
! isou             ! e  ! <-- ! numero de passage                              !
! ipp              ! e  ! <-- ! numero de variable pour sorties post           !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! produc           ! tr ! <-- ! production                                     !
!  (6,ncelet)      !    !     !                                                !
! smbr(ncelet      ! tr ! <-- ! tableau de travail pour sec mem                !
! coefax,coefbx    ! tr ! --- ! tableau de travail pour les cond.              !
!  (nfabor)        !    !     !    aux limites de la dist. paroi               !
! produk(ncelet    ! tr ! --- ! tableau de travail production                  !
! epsk  (ncelet    ! tr ! --- ! tableau de travail epsilon/k                   !
! w2...6(ncelet    ! tr ! --- ! tableau de travail production                  !
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
use pointe
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          iphas  , ivar   , isou   , ipp

integer          ia(*)

double precision rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision produc(6,ncelet)
double precision smbr(ncelet)
double precision coefax(nfabor), coefbx(nfabor)
double precision produk(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),epsk(ncelet),w6(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifacpt, iel   , ii    , jj    , kk    , mm
integer          irkm  , irki  , irkj  , iskm  , iski  , iskj
integer          ir11ip, ir22ip, ir33ip, ir12ip, ir13ip, ir23ip
integer          ieiph , ipcrom, ipcroo
integer          ifac
integer          inc   , iccocg, iphydp, ivar0 , iityph

double precision cmu075, distxn, d2s3  , trrij , xk
double precision unssur, vnk   , vnm   , vni   , vnj
double precision deltki, deltkj, deltkm, deltij
double precision aa    , bb    , xnorme


!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Initialize variables to avoid compiler warnings

ii = 0
jj = 0
irki = 0
irkj = 0
irkm = 0
iski = 0
iskj = 0
iskm = 0

vni = 0.d0
vnj = 0.d0
vnk = 0.d0
vnm = 0.d0

! Memoire

idebia = idbia0
idebra = idbra0

ir11ip = ir11
ir22ip = ir22
ir33ip = ir33
ir12ip = ir12
ir13ip = ir13
ir23ip = ir23
ieiph  = iep

ipcrom = ipproc(irom  )
ipcroo = ipcrom
if(isto2t.gt.0.and.iroext.gt.0) then
  ipcroo = ipproc(iroma)
endif

deltij = 1.0d0
if(isou.gt.3) then
  deltij = 0.0d0
endif

cmu075 = cmu**0.75d0
d2s3   = 2.d0/3.d0

!===============================================================================
! 2. CALCUL AUX CELLULES DES NORMALES EN PAROI CORRESPONDANTES
!===============================================================================

!     On stocke les composantes dans les tableaux de travail W2, W3 et W4

if(abs(icdpar).eq.2) then

!     On connait la face de paroi correspondante

    do iel = 1, ncel
      ifacpt = ia(iifapa-1+iel)
      unssur = 1.d0/surfbn(ifacpt)
      w2(iel)= surfbo(1,ifacpt)*unssur
      w3(iel)= surfbo(2,ifacpt)*unssur
      w4(iel)= surfbo(3,ifacpt)*unssur
    enddo

elseif(abs(icdpar).eq.1) then

!     La normale est definie comme - gradient de la distance
!       a la paroi

!       La distance a la paroi vaut 0 en paroi
!         par definition et obeit a un flux nul ailleurs

  iityph = iitypf+(iphas-1)*nfabor
  do ifac = 1, nfabor
    if(ia(iityph+ifac-1).eq.iparoi .or.                           &
       ia(iityph+ifac-1).eq.iparug) then
      coefax(ifac) = 0.0d0
      coefbx(ifac) = 0.0d0
    else
      coefax(ifac) = 0.0d0
      coefbx(ifac) = 1.0d0
    endif
  enddo

!       Calcul du gradient

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(ra(idipar))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iphydp = 0
  ivar0  = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   ivar0  , imrgra , inc    , iccocg , nswrgy , imligy , iphydp , &
   iwarny , nfecra ,                                              &
   epsrgy , climgy , extray ,                                     &
   ia     ,                                                       &
   ra     , ra     , ra     ,                                     &
   ra(idipar) , coefax , coefbx ,                                 &
   w2     , w3     , w4     ,                                     &
!        ------   ------   ------
   produk , epsk   , w6     ,                                     &
   ra     )


!     Normalisation (attention, le gradient peut etre nul, parfois)

  do iel = 1 ,ncel
    xnorme = max(sqrt(w2(iel)**2+w3(iel)**2+w4(iel)**2),epzero)
    w2(iel) = -w2(iel)/xnorme
    w3(iel) = -w3(iel)/xnorme
    w4(iel) = -w4(iel)/xnorme
  enddo

endif

!===============================================================================
! 2. CALCUL DE VARIABLES DE TRAVAIL
!===============================================================================



! ---> Production et k

do iel = 1 , ncel
  produk(iel) = 0.5d0 * (produc(1,iel)  + produc(2,iel)  +        &
                         produc(3,iel))
  xk          = 0.5d0 * (rtpa(iel,ir11ip) + rtpa(iel,ir22ip) +    &
                         rtpa(iel,ir33ip))
  epsk(iel)   = rtpa(iel,ieiph)/xk
enddo



! ---> Indices des tensions

if     ((isou.eq.1).or.(isou.eq.4).or.(isou.eq.5)) then
  ii = 1
elseif ((isou.eq.2).or.(isou.eq.6)) then
  ii = 2
elseif  (isou.eq.3) then
  ii = 3
endif

if     ((isou.eq.3).or.(isou.eq.5).or.(isou.eq.6)) then
  jj = 3
elseif ((isou.eq.2).or.(isou.eq.4)) then
  jj = 2
elseif  (isou.eq.1) then
  jj = 1
endif

! ---> Boucle pour construction des termes sources

do iel = 1, ncel
  w6(iel) = 0.d0
enddo

do kk = 1, 3

! ---> Sommes sur m

  do mm = 1, 3

!   --> Delta km

    if(kk.eq.mm) then
      deltkm = 1.0d0
    else
      deltkm = 0.0d0
    endif

!  --> R km

    if     ((kk*mm).eq.1) then
      irkm = ir11ip
      iskm = 1
    elseif ((kk*mm).eq.4) then
      irkm = ir22ip
      iskm = 2
    elseif ((kk*mm).eq.9) then
      irkm = ir33ip
      iskm = 3
    elseif ((kk*mm).eq.2) then
      irkm = ir12ip
      iskm = 4
    elseif ((kk*mm).eq.3) then
      irkm = ir13ip
      iskm = 5
    elseif ((kk*mm).eq.6) then
      irkm = ir23ip
      iskm = 6
    endif

!  --> Termes en R km et Phi km

    do iel = 1, ncel

      if    (kk.eq.1) then
        vnk    = w2(iel)
      elseif(kk.eq.2) then
        vnk    = w3(iel)
      elseif(kk.eq.3) then
        vnk    = w4(iel)
      endif

      if    (mm.eq.1) then
        vnm    = w2(iel)
      elseif(mm.eq.2) then
        vnm    = w3(iel)
      elseif(mm.eq.3) then
        vnm    = w4(iel)
      endif

      w6(iel) = w6(iel) + vnk*vnm*deltij*(                        &
             crijp1*rtpa(iel,irkm)*epsk(iel)                      &
            -crijp2                                               &
             *crij2*(produc(iskm,iel)-d2s3*produk(iel)*deltkm) )
    enddo

  enddo

! ---> Fin des sommes sur m


!  --> R ki

  if     ((kk*ii).eq.1) then
    irki = ir11ip
    iski = 1
  elseif ((kk*ii).eq.4) then
    irki = ir22ip
    iski = 2
  elseif ((kk*ii).eq.9) then
    irki = ir33ip
    iski = 3
  elseif ((kk*ii).eq.2) then
    irki = ir12ip
    iski = 4
  elseif ((kk*ii).eq.3) then
    irki = ir13ip
    iski = 5
  elseif ((kk*ii).eq.6) then
    irki = ir23ip
    iski = 6
  endif

!  --> R kj

  if     ((kk*jj).eq.1) then
    irkj = ir11ip
    iskj = 1
  elseif ((kk*jj).eq.4) then
    irkj = ir22ip
    iskj = 2
  elseif ((kk*jj).eq.9) then
    irkj = ir33ip
    iskj = 3
  elseif ((kk*jj).eq.2) then
    irkj = ir12ip
    iskj = 4
  elseif ((kk*jj).eq.3) then
    irkj = ir13ip
    iskj = 5
  elseif ((kk*jj).eq.6) then
    irkj = ir23ip
    iskj = 6
  endif

!   --> Delta ki

  if (kk.eq.ii) then
    deltki = 1.d0
  else
    deltki = 0.d0
  endif

!   --> Delta kj

  if (kk.eq.jj) then
    deltkj = 1.d0
  else
    deltkj = 0.d0
  endif

  do iel = 1, ncel

      if    (kk.eq.1) then
        vnk    = w2(iel)
      elseif(kk.eq.2) then
        vnk    = w3(iel)
      elseif(kk.eq.3) then
        vnk    = w4(iel)
      endif

      if    (ii.eq.1) then
        vni    = w2(iel)
      elseif(ii.eq.2) then
        vni    = w3(iel)
      elseif(ii.eq.3) then
        vni    = w4(iel)
      endif

      if    (jj.eq.1) then
        vnj    = w2(iel)
      elseif(jj.eq.2) then
        vnj    = w3(iel)
      elseif(jj.eq.3) then
        vnj    = w4(iel)
      endif

    w6(iel) = w6(iel) + 1.5d0*vnk*(                               &
    -crijp1*(rtpa(iel,irki)*vnj+rtpa(iel,irkj)*vni)*epsk(iel)     &
    +crijp2                                                       &
     *crij2*((produc(iski,iel)-d2s3*produk(iel)*deltki)*vnj       &
            +(produc(iskj,iel)-d2s3*produk(iel)*deltkj)*vni) )

  enddo

enddo


! ---> Distance a la paroi et fonction d'amortissement : W3
!     Pour chaque mode de calcul : meme code, test
!       en dehors de la boucle

if(abs(icdpar).eq.2) then
  do iel = 1 , ncel
    ifacpt = ia(iifapa-1+iel)
    distxn =                                                      &
          (cdgfbo(1,ifacpt)-xyzcen(1,iel))**2                     &
         +(cdgfbo(2,ifacpt)-xyzcen(2,iel))**2                     &
         +(cdgfbo(3,ifacpt)-xyzcen(3,iel))**2
    distxn = sqrt(distxn)
    trrij  = 0.5d0 * (rtpa(iel,ir11ip) + rtpa(iel,ir22ip) +       &
                         rtpa(iel,ir33ip))
    aa = 1.d0
    bb = cmu075*trrij**1.5d0/(xkappa*rtpa(iel,ieiph)*distxn)
    w3(iel) = min(aa, bb)
  enddo
else
  do iel = 1 , ncel
    distxn =  max(ra(idipar+iel-1),epzero)
    trrij  = 0.5d0 * (rtpa(iel,ir11ip) + rtpa(iel,ir22ip) +       &
                         rtpa(iel,ir33ip))
    aa = 1.d0
    bb = cmu075*trrij**1.5d0/(xkappa*rtpa(iel,ieiph)*distxn)
    w3(iel) = min(aa, bb)
  enddo
endif


! ---> Increment du terme source

do iel = 1, ncel
  smbr(iel) = smbr(iel)                                           &
            + propce(iel,ipcroo)*volume(iel)*w6(iel)*w3(iel)
enddo


return

end subroutine
