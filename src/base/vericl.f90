!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine vericl &
!================

 ( nvar   , nscal  ,                                              &
   icodcl ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl )

!===============================================================================
! FONCTION :
! --------

! VERIFICATION DE ICODCL

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
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
use optcal
use cstnum
use cstphy
use entsor
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)

! Local variables

character*80     chaine
integer          ifac, ivar, icode
integer          nstoni        , nstvit, nstopp
integer          nstoke, nstosc, nstovf
integer          nstuvw, nstoup, nstuke
integer          nstrij, nsurij, nstov2
integer          nstuv2, nstokw, nstukw
integer          nstunu, nstonu
integer          nstusc
integer          iis, icodcu, icodcv, icodcw, icodck, icodce
integer          icodcn
integer          icodcp, icodcf, icodca, icodom
integer          icor11, icor22, icor33, icor12, icor13, icor23
integer          ipp, iokcod, iok
integer          ippprp, ippuip, ippvip, ippwip, ippepp, ippkip
integer          ipp11p, ipp22p, ipp33p, ipp12p, ipp13p, ipp23p
integer          ippphp, ippfbp, ippalp, ippomg
integer          ippnup

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

ippepp = 0
ippkip = 0
ippfbp = 0
ippomg = 0
ippphp = 0
ippnup = 0
ipp11p = 0
ipp22p = 0
ipp33p = 0
ipp12p = 0
ipp13p = 0
ipp23p = 0

! Memoire


!===============================================================================
! 2.  VERIFICATION DE LA CONSISTANCE DES CL
!===============================================================================

! 2.1  INITIALISATION
! ====================

! Dans USCLIM, on se donne une grande liberte pour la specif. des c.l
!  sur les variables. neanmoins pour limiter la plage des tests, on se
!  donne, pour l'instant, les contraintes suivantes :

!   - meme type de c.l pour les 3 composantes de vitesse
!   - pas de conditions de frottemt sur la pression
!   - coherence entre les c.l vitesses et pression
!   - coherence entre les c.l vitesses et turbulence

nstoni = 0
nstosc = 0
nstovf = 0
nstusc = 0
nstvit = 0
nstopp = 0
nstoke = 0
nstrij = 0
nstov2 = 0
nstokw = 0
nstuvw = 0
nstoup = 0
nstuke = 0
nsurij = 0
nstuv2 = 0
nstukw = 0
nstunu = 0
nstonu = 0


! 2.2 VERIFICATIONS QUE TOUTES LES CL SONT INITIALISEES
! ======================================================

! --- Premiere boucle rapide
iokcod = 0
do ivar = 1, nvar
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if(icode.eq. 0) then
      iokcod = 1
    endif
  enddo
enddo

! --- Seconde boucle lente si pb plus haut
if(iokcod.ne.0) then
  do ipp = 2, nvppmx
    if (itrsvr(ipp).ge.1) then
      ivar = itrsvr(ipp)
      do ifac = 1, nfabor
        icode = icodcl(ifac,ivar)
        if(icode.eq. 0) then
          chaine=nomvar(ipp)
          write(nfecra,1000)ifac,iprfml(ifmfbr(ifac),1),          &
                            chaine(1:8),icodcl(ifac,ivar)
          nstoni = nstoni + 1
        endif
      enddo
    endif
  enddo
endif


! 2.3 VERIFICATIONS DE L'ADMISSIBILITE DES CONDITIONS
! ====================================================

ippprp = ipprtp(ipr)
ippuip = ipprtp(iu )
ippvip = ipprtp(iv )
ippwip = ipprtp(iw )
if(itytur.eq.2) then
  ippkip = ipprtp(ik )
  ippepp = ipprtp(iep)
elseif(itytur.eq.3) then
  ipp11p = ipprtp(ir11)
  ipp22p = ipprtp(ir22)
  ipp33p = ipprtp(ir33)
  ipp12p = ipprtp(ir12)
  ipp13p = ipprtp(ir13)
  ipp23p = ipprtp(ir23)
  ippepp = ipprtp(iep)
elseif(itytur.eq.5) then
  ippkip = ipprtp(ik )
  ippepp = ipprtp(iep)
  ippphp = ipprtp(iphi)
  if(iturb.eq.50) then
    ippfbp = ipprtp(ifb)
  elseif(iturb.eq.51) then
    ippalp = ipprtp(ial)
  endif
elseif(iturb.eq.60) then
  ippkip = ipprtp(ik )
  ippomg = ipprtp(iomg)
elseif(iturb.eq.70) then
  ippnup = ipprtp(inusa )
endif

! --- Conditions admissibles pour les composantes de vitesse
do ifac = 1, nfabor

  icodcu = icodcl(ifac,iu)
  icodcv = icodcl(ifac,iv)
  icodcw = icodcl(ifac,iw)

  if(icodcu.ne. 1.and.                 icodcu.ne. 3.and.        &
       icodcu.ne. 4.and.icodcu.ne. 5.and.icodcu.ne. 6.and.        &
       icodcu.ne. 9) then
    chaine=nomvar(ippuip)
    write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
         icodcl(ifac,iu)
    nstvit = nstvit + 1
  endif
  if(icodcv.ne. 1.and.                 icodcv.ne. 3.and.        &
       icodcv.ne. 4.and.icodcv.ne. 5.and.icodcv.ne. 6.and.        &
       icodcv.ne. 9) then
    chaine=nomvar(ippvip )
    write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
         icodcl(ifac,iv)
    nstvit = nstvit + 1
  endif
  if(icodcw.ne. 1.and.                 icodcw.ne. 3.and.        &
       icodcw.ne. 4.and.icodcw.ne. 5.and.icodcv.ne. 6.and.        &
       icodcw.ne. 9) then
    chaine=nomvar(ippwip)
    write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
         icodcl(ifac,iw)
    nstvit = nstvit + 1
  endif

  ! --- verification que la rugosite est initialisee si icodl=6
  if(icodcu.eq.6 .and. rcodcl(ifac,iu,3).lt.epzero)then
    CHAINE='RUGOSITV'
    write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
         icodcl(ifac,iu)
    nstvit = nstvit + 1
  endif

  ! --- on interdit les parois rugueuses en compressible
  if (icodcu.eq.6 .and. ippmod(icompf).gt.0) then
    chaine=nomvar(ippuip)
    write(nfecra,1015)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),  &
         icodcl(ifac,iu),ippmod(icompf)
    nstvit = nstvit + 1
  endif

enddo

! --- Conditions admissibles pour la pression
do ifac = 1, nfabor

  if(icodcl(ifac,ipr).ne. 1.and.                             &
       icodcl(ifac,ipr).ne. 3) then
    chaine=nomvar(ippprp)
    write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),  &
         icodcl(ifac,ipr)
    nstopp = nstopp + 1
  endif

enddo

! --- Conditions admissibles pour k et epsilon
if (itytur.eq.2) then

  do ifac = 1, nfabor

    if((icodcl(ifac,ik ).ne. 1.and.                          &
         icodcl(ifac,ik ).ne. 3.and.                          &
         icodcl(ifac,ik ).ne. 5.and.                          &
         icodcl(ifac,ik ).ne. 6     ).or.                     &
         (icodcl(ifac,iep).ne. 1.and.                          &
         icodcl(ifac,iep).ne. 3.and.                          &
         icodcl(ifac,iep).ne. 5.and.                          &
         icodcl(ifac,iep).ne. 6     ) )then
      chaine=nomvar(ippkip)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
           icodcl(ifac,ik )
      chaine=nomvar(ippepp)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
           icodcl(ifac,iep)
      nstoke = nstoke + 1
    endif

  enddo

  ! --- Conditions admissibles pour Rij et epsilon
elseif(itytur.eq.3) then

  ivar = ir11
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if(icode.ne. 1.and.                icode.ne. 3.and.         &
         icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6     ) then
      chaine=nomvar(ipp11p)
      write(nfecra,1010)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),icode
      nstrij = nstrij + 1
    endif
  enddo

  ivar = ir22
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if(icode.ne. 1.and.                icode.ne. 3.and.         &
         icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6     ) then
      chaine=nomvar(ipp22p)
      write(nfecra,1010)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),icode
      nstrij = nstrij + 1
    endif
  enddo

  ivar = ir33
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if(icode.ne. 1.and.                icode.ne. 3.and.         &
         icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6     ) then
      chaine=nomvar(ipp33p)
      write(nfecra,1010)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),icode
      nstrij = nstrij + 1
    endif
  enddo

  ivar = ir12
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if(icode.ne. 1.and.                icode.ne. 3.and.         &
         icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6     ) then
      chaine=nomvar(ipp12p)
      write(nfecra,1010)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),icode
      nstrij = nstrij + 1
    endif
  enddo

  ivar = ir13
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if(icode.ne. 1.and.                icode.ne. 3.and.         &
         icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6     ) then
      chaine=nomvar(ipp13p)
      write(nfecra,1010)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),icode
      nstrij = nstrij + 1
    endif
  enddo

  ivar = ir23
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if(icode.ne. 1.and.                icode.ne. 3.and.         &
         icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6     ) then
      chaine=nomvar(ipp23p)
      write(nfecra,1010)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),icode
      nstrij = nstrij + 1
    endif
  enddo

  do ifac = 1, nfabor
    icode = icodcl(ifac,iep)
    if(icode.ne. 1.and.                icode.ne. 3.and.         &
         icode.ne. 5.and.icode.ne. 6     ) then
      chaine=nomvar(ippepp)
      write(nfecra,1010)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),icode
      nstrij = nstrij + 1
    endif
  enddo

  ! --- Conditions admissibles pour k, epsilon, phi et f_barre
elseif (iturb.eq.50) then

  do ifac = 1, nfabor

    if((icodcl(ifac,ik ).ne. 1.and.                          &
         icodcl(ifac,ik ).ne. 3.and.                          &
         icodcl(ifac,ik ).ne. 5.and.                          &
         icodcl(ifac,ik ).ne. 6     ).or.                     &
         (icodcl(ifac,iep).ne. 1.and.                          &
         icodcl(ifac,iep).ne. 3.and.                          &
         icodcl(ifac,iep).ne. 5.and.                          &
         icodcl(ifac,iep).ne. 6     ).or.                     &
         (icodcl(ifac,iphi).ne. 1.and.                          &
         icodcl(ifac,iphi).ne. 3.and.                          &
         icodcl(ifac,iphi).ne. 5.and.                          &
         icodcl(ifac,iphi).ne. 6     ).or.                     &
         (icodcl(ifac,ifb).ne. 1.and.                          &
         icodcl(ifac,ifb).ne. 3.and.                          &
         icodcl(ifac,ifb).ne. 5.and.                          &
         icodcl(ifac,ifb).ne. 6     ) )then
      chaine=nomvar(ippkip)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
           icodcl(ifac,ik )
      chaine=nomvar(ippepp)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
           icodcl(ifac,iep)
      chaine=nomvar(ippphp)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
           icodcl(ifac,iphi )
      chaine=nomvar(ippfbp)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
           icodcl(ifac,ifb)
      nstov2 = nstov2 + 1

    endif

  enddo

! --- Conditions admissibles pour k, epsilon, phi et alpha
elseif (iturb.eq.51) then

  do ifac = 1, nfabor

    if((icodcl(ifac,ik ).ne. 1.and.                          &
        icodcl(ifac,ik ).ne. 3.and.                          &
        icodcl(ifac,ik ).ne. 5.and.                          &
        icodcl(ifac,ik ).ne. 6     ).or.                     &
       (icodcl(ifac,iep).ne. 1.and.                          &
        icodcl(ifac,iep).ne. 3.and.                          &
        icodcl(ifac,iep).ne. 5.and.                          &
        icodcl(ifac,iep).ne. 6     ).or.                     &
       (icodcl(ifac,iphi).ne. 1.and.                          &
        icodcl(ifac,iphi).ne. 3.and.                          &
        icodcl(ifac,iphi).ne. 5.and.                          &
        icodcl(ifac,iphi).ne. 6     ).or.                     &
        (icodcl(ifac,ial).ne. 1.and.                          &
        icodcl(ifac,ial).ne. 3.and.                          &
        icodcl(ifac,ial).ne. 5.and.                          &
        icodcl(ifac,ial).ne. 6     ) )then
      chaine=nomvar(ippkip)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
                        icodcl(ifac,ik )
      chaine=nomvar(ippepp)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
                        icodcl(ifac,iep)
      chaine=nomvar(ippphp)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
                        icodcl(ifac,iphi )
      chaine=nomvar(ippalp)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
                        icodcl(ifac,ial)
      nstov2 = nstov2 + 1

    endif

  enddo

! --- Conditions admissibles pour k et omega
elseif (iturb.eq.60) then

  do ifac = 1, nfabor

    if((icodcl(ifac,ik ).ne. 1.and.                          &
         icodcl(ifac,ik ).ne. 3.and.                          &
         icodcl(ifac,ik ).ne. 5.and.                          &
         icodcl(ifac,ik ).ne. 6     ).or.                     &
         (icodcl(ifac,iomg).ne. 1.and.                          &
         icodcl(ifac,iomg).ne. 3.and.                          &
         icodcl(ifac,iomg).ne. 5.and.                          &
         icodcl(ifac,iomg).ne. 6     ) )then
      chaine=nomvar(ippkip)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
           icodcl(ifac,ik )
      chaine=nomvar(ippomg)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
           icodcl(ifac,iomg)
      nstokw = nstokw + 1
    endif

  enddo

  ! --- Conditions admissibles pour Spalart-Allmaras
elseif (iturb.eq.70) then

  do ifac = 1, nfabor

    if(icodcl(ifac,inusa ).ne. 1.and.                           &
         icodcl(ifac,inusa ).ne. 3.and.                          &
         icodcl(ifac,inusa ).ne. 5.and.                          &
         icodcl(ifac,inusa ).ne. 6           )then
      chaine=nomvar(ippnup)
      write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
           icodcl(ifac,inusa )
      nstonu = nstonu + 1
    endif

  enddo

endif

! --- Conditions admissibles pour les scalaires
if(nscal.ge.1) then
  do iis = 1,nscal
    ivar = isca(iis)
    do ifac = 1, nfabor
      if(icodcl(ifac,ivar).ne. 1.and.                             &
         icodcl(ifac,ivar).ne. 3.and.                             &
         icodcl(ifac,ivar).ne. 5.and.                             &
         icodcl(ifac,ivar).ne. 6 ) then
        chaine=nomvar(ipprtp(ivar))
        write(nfecra,1010)                                        &
          ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
          icodcl(ifac,ivar)
        nstosc = nstosc + 1
      endif
      if(icodcl(ifac,ivar).eq. 5.and.                             &
         iscavr(iis).gt.0        ) then
        chaine=nomvar(ipprtp(ivar))
        write(nfecra,1010)                                        &
          ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
          icodcl(ifac,ivar)
        nstovf = nstovf + 1
      endif
      if(icodcl(ifac,ivar).eq. 6.and.                             &
         iscavr(iis).gt.0        ) then
        chaine=nomvar(ipprtp(ivar))
        write(nfecra,1010)                                        &
          ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
          icodcl(ifac,ivar)
        nstovf = nstovf + 1
      endif
! --- verification que la rugosite scalaire est initialisee si icodl=6
      if(icodcl(ifac,ivar).eq.6.and.                              &
         rcodcl(ifac,iv,3).lt.epzero)then
        CHAINE='RUGOSITS'
        write(nfecra,1010)ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),&
                          icodcl(ifac,ivar)
        nstosc = nstosc + 1
      endif
    enddo
  enddo
endif

! 2.4 VERIFICATIONS DES COHERENCES INTER VARIABLES INTRA PHASE
! =============================================================

ippprp = ipprtp(ipr)
ippuip = ipprtp(iu )
ippvip = ipprtp(iv )
ippwip = ipprtp(iw )
if(itytur.eq.2) then
  ippkip = ipprtp(ik )
  ippepp = ipprtp(iep)
elseif(itytur.eq.3) then
  ipp11p = ipprtp(ir11)
  ipp22p = ipprtp(ir22)
  ipp33p = ipprtp(ir33)
  ipp12p = ipprtp(ir12)
  ipp13p = ipprtp(ir13)
  ipp23p = ipprtp(ir23)
  ippepp = ipprtp(iep)
elseif(itytur.eq.5) then
  ippkip = ipprtp(ik )
  ippepp = ipprtp(iep)
  ippphp = ipprtp(iphi)
  if(iturb.eq.50) then
    ippfbp = ipprtp(ifb)
  elseif(iturb.eq.51) then
    ippalp = ipprtp(ial)
  endif
elseif(iturb.eq.60) then
  ippkip = ipprtp(ik )
  ippomg = ipprtp(iomg)
elseif(iturb.eq.70) then
  ippnup = ipprtp(inusa)
endif

! --- Coherence pour les composantes de vitesse
do ifac = 1, nfabor

  icodcu = icodcl(ifac,iu)
  icodcv = icodcl(ifac,iv)
  icodcw = icodcl(ifac,iw)

  if(icodcu.eq.4.or.icodcu.eq.5.or.icodcu.eq.6.or.              &
       icodcu.eq.9.or.                             &
       icodcv.eq.4.or.icodcv.eq.5.or.icodcv.eq.6.or.              &
       icodcv.eq.9.or.                             &
       icodcw.eq.4.or.icodcw.eq.5.or.icodcw.eq.6.or.              &
       icodcw.eq.9                )then

    if( icodcu.ne.icodcv .or. icodcu.ne.icodcw .or.             &
         icodcv.ne.icodcw ) then
      write(nfecra,1020)ifac,iprfml(ifmfbr(ifac),1),            &
           icodcu,icodcv,icodcw
      nstuvw = nstuvw + 1
    endif
  endif

  ! --- Coherence vitesse pression

  !      Remarques :
  !        Pas de regle stricte de coherence vitesse/pression.
  !        Avant on imposait un Dirichlet sur la pression pour en
  !        entree/sortie, mais cela ne semble pas imperatif. Le test
  !        est laisse en commentaire pour etre recupere si necessaire.

  !          IF( ICODCU.EQ.9 .OR. ICODCV.EQ.9 .OR. ICODCW.EQ.9 ) THEN
  !            IF( ICODCL(IFAC,IPRIPH).NE.1                    ) THEN
  !              CHAINE=NOMVAR(IPPPRP)
  !              WRITE(NFECRA,1030)
  !     &          IFAC,IPRFML(IFMFBR(IFAC),1),CHAINE(1:8)
  !     &          ICODCL(IFAC,IPRIPH),ICODCU,ICODCV,ICODCW
  !              NSTOUP = NSTOUP + 1
  !            ENDIF
  !          ENDIF

enddo

! --- Coherence vitesse turbulence

if(itytur.eq.2) then

  do ifac = 1, nfabor

    icodcu = icodcl(ifac,iu)
    icodcv = icodcl(ifac,iv)
    icodcw = icodcl(ifac,iw)
    icodck = icodcl(ifac,ik)
    icodce = icodcl(ifac,iep)

    if( (icodcu.eq.5 .or. icodcv.eq.5 .or. icodcw.eq.5 .or.     &
         icodck.eq.5 .or. icodce.eq.5) .and.                    &
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.     &
         icodck.ne.5 .or. icodce.ne.5)                    ) then
      chaine=nomvar(ippkip)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,ik),icodcu,icodcv,icodcw
      chaine=nomvar(ippepp)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,iep),icodcu,icodcv,icodcw
      nstuke = nstuke + 1
    endif

    if( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.     &
         icodck.eq.6 .or. icodce.eq.6) .and.                    &
         (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.     &
         icodck.ne.6 .or. icodce.ne.6)                    ) then
      chaine=nomvar(ippkip)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,ik),icodcu,icodcv,icodcw
      chaine=nomvar(ippepp)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,iep),icodcu,icodcv,icodcw
      nstuke = nstuke + 1
    endif

  enddo

elseif(itytur.eq.3) then

  do ifac = 1, nfabor

    icodcu = icodcl(ifac,iu)
    icodcv = icodcl(ifac,iv)
    icodcw = icodcl(ifac,iw)
    icor11 = icodcl(ifac,ir11)
    icor22 = icodcl(ifac,ir22)
    icor33 = icodcl(ifac,ir33)
    icor12 = icodcl(ifac,ir12)
    icor13 = icodcl(ifac,ir13)
    icor23 = icodcl(ifac,ir23)
    icodce = icodcl(ifac,iep)

    if(  (icodcu.eq.5 .or. icodcv.eq.5 .or. icodcw.eq.5 .or.    &
         icor11.eq.5 .or. icor22.eq.5 .or.                     &
         icor33.eq.5 .or. icor12.eq.5 .or.                     &
         icor13.eq.5 .or. icor23.eq.5 .or.                     &
         icodce.eq.5                      ) .and.              &
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.    &
         icor11.ne.5 .or. icor22.ne.5 .or.                     &
         icor33.ne.5 .or. icor12.ne.5 .or.                     &
         icor13.ne.5 .or. icor23.ne.5 .or.                     &
         icodce.ne.5                      )      ) then
      write(nfecra,1040)                                        &
           ifac,iprfml(ifmfbr(ifac),1),                            &
           icor11,icor22,icor33,                                   &
           icor12,icor13,icor23,                                   &
           icodce,icodcu,icodcv,icodcw
      nsurij = nsurij + 1
    endif

    if(  (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.    &
         icor11.eq.6 .or. icor22.eq.6 .or.                     &
         icor33.eq.6 .or. icor12.eq.6 .or.                     &
         icor13.eq.6 .or. icor23.eq.6 .or.                     &
         icodce.eq.6                      ) .and.              &
         (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.    &
         icor11.ne.6 .or. icor22.ne.6 .or.                     &
         icor33.ne.6 .or. icor12.ne.6 .or.                     &
         icor13.ne.6 .or. icor23.ne.6 .or.                     &
         icodce.ne.6                      )      ) then
      write(nfecra,1040)                                        &
           ifac,iprfml(ifmfbr(ifac),1),                            &
           icor11,icor22,icor33,                                   &
           icor12,icor13,icor23,                                   &
           icodce,icodcu,icodcv,icodcw
      nsurij = nsurij + 1
    endif

    if(  (icodcu.eq.4 .or. icodcv.eq.4 .or. icodcw.eq.4 .or.    &
         icor11.eq.4 .or. icor22.eq.4 .or.                     &
         icor33.eq.4 .or. icor12.eq.4 .or.                     &
         icor13.eq.4 .or. icor23.eq.4                          &
         ) .and.                         &
         (icodcu.ne.4 .or. icodcv.ne.4 .or. icodcw.ne.4 .or.    &
         icor11.ne.4 .or. icor22.ne.4 .or.                     &
         icor33.ne.4 .or. icor12.ne.4 .or.                     &
         icor13.ne.4 .or. icor23.ne.4 .or.                     &
         icodce.ne.3) ) then
      write(nfecra,1040)                                        &
           ifac,iprfml(ifmfbr(ifac),1),                            &
           icor11,icor22,icor33,                                   &
           icor12,icor13,icor23,                                   &
           icodce,icodcu,icodcv,icodcw
      nsurij = nsurij + 1
    endif

  enddo

elseif(iturb.eq.50 ) then

  do ifac = 1, nfabor

    icodcu = icodcl(ifac,iu)
    icodcv = icodcl(ifac,iv)
    icodcw = icodcl(ifac,iw)
    icodck = icodcl(ifac,ik)
    icodce = icodcl(ifac,iep)
    icodcp = icodcl(ifac,iphi)
    icodcf = icodcl(ifac,ifb)

    if( (icodcu.eq.5 .or. icodcv.eq.5 .or. icodcw.eq.5 .or.     &
         icodck.eq.5 .or. icodce.eq.5 .or. icodcp.eq.5 .or.     &
         icodcf.eq.5 ) .and.                                    &
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.     &
         icodck.ne.5 .or. icodce.ne.5 .or. icodcp.ne.5 .or.     &
         icodcf.ne.5 )                    ) then
      chaine=nomvar(ippkip)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,ik),icodcu,icodcv,icodcw
      chaine=nomvar(ippepp)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,iep),icodcu,icodcv,icodcw
      chaine=nomvar(ippphp)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,iphi),icodcu,icodcv,icodcw
      chaine=nomvar(ippfbp)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,ifb),icodcu,icodcv,icodcw
      nstuv2 = nstuv2 + 1

      if( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.     &
           icodck.eq.6 .or. icodce.eq.6 .or. icodcp.eq.6 .or.     &
           icodcf.eq.6 ) .and.                                    &
           (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.     &
           icodck.ne.6 .or. icodce.ne.6 .or. icodcp.ne.6 .or.     &
           icodcf.ne.6 )                    ) then
        chaine=nomvar(ippkip)
        write(nfecra,1030)                                        &
             ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
             icodcl(ifac,ik),icodcu,icodcv,icodcw
        chaine=nomvar(ippepp)
        write(nfecra,1030)                                        &
             ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
             icodcl(ifac,iep),icodcu,icodcv,icodcw
        chaine=nomvar(ippphp)
        write(nfecra,1030)                                        &
             ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
             icodcl(ifac,iphi),icodcu,icodcv,icodcw
        chaine=nomvar(ippfbp)
        write(nfecra,1030)                                        &
             ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
             icodcl(ifac,ifb),icodcu,icodcv,icodcw
        nstuv2 = nstuv2 + 1

      endif

    endif

  enddo

elseif(iturb.eq.51 ) then

  do ifac = 1, nfabor

    icodcu = icodcl(ifac,iu)
    icodcv = icodcl(ifac,iv)
    icodcw = icodcl(ifac,iw)
    icodck = icodcl(ifac,ik)
    icodce = icodcl(ifac,iep)
    icodcp = icodcl(ifac,iphi)
    icodca = icodcl(ifac,ial)

    if( (icodcu.eq.5 .or. icodcv.eq.5 .or. icodcw.eq.5 .or.     &
         icodck.eq.5 .or. icodce.eq.5 .or. icodcp.eq.5 .or.     &
         icodca.eq.5 ) .and.                                    &
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.     &
         icodck.ne.5 .or. icodce.ne.5 .or. icodcp.ne.5 .or.     &
         icodca.ne.5 )                    ) then
      chaine=nomvar(ippkip)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),          &
           icodcl(ifac,ik),icodcu,icodcv,icodcw
      chaine=nomvar(ippepp)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),          &
           icodcl(ifac,iep),icodcu,icodcv,icodcw
      chaine=nomvar(ippphp)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),          &
           icodcl(ifac,iphi),icodcu,icodcv,icodcw
      chaine=nomvar(ippalp)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),          &
           icodcl(ifac,ial),icodcu,icodcv,icodcw
      nstuv2 = nstuv2 + 1

      if( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.     &
           icodck.eq.6 .or. icodce.eq.6 .or. icodcp.eq.6 .or.     &
           icodca.eq.6 ) .and.                                    &
           (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.     &
           icodck.ne.6 .or. icodce.ne.6 .or. icodcp.ne.6 .or.     &
           icodca.ne.6 )                    ) then
        chaine=nomvar(ippkip)
        write(nfecra,1030)                                        &
             ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),          &
             icodcl(ifac,ik),icodcu,icodcv,icodcw
        chaine=nomvar(ippepp)
        write(nfecra,1030)                                        &
             ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),          &
             icodcl(ifac,iep),icodcu,icodcv,icodcw
        chaine=nomvar(ippphp)
        write(nfecra,1030)                                        &
             ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),          &
             icodcl(ifac,iphi),icodcu,icodcv,icodcw
        chaine=nomvar(ippalp)
        write(nfecra,1030)                                        &
             ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),          &
             icodcl(ifac,ial),icodcu,icodcv,icodcw
        nstuv2 = nstuv2 + 1

      endif

    endif

  enddo

elseif(iturb.eq.60 ) then

  do ifac = 1, nfabor

    icodcu = icodcl(ifac,iu)
    icodcv = icodcl(ifac,iv)
    icodcw = icodcl(ifac,iw)
    icodck = icodcl(ifac,ik)
    icodom = icodcl(ifac,iomg)

    if( (icodcu.eq.5 .or. icodcv.eq.5 .or. icodcw.eq.5 .or.     &
         icodck.eq.5 .or. icodom.eq.5 ) .and.                   &
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.     &
         icodck.ne.5 .or. icodom.ne.5 ) ) then
      chaine=nomvar(ippkip)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,ik),icodcu,icodcv,icodcw
      chaine=nomvar(ippomg)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,iomg),icodcu,icodcv,icodcw
      nstukw = nstukw + 1
    endif

    if( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.     &
         icodck.eq.6 .or. icodom.eq.6 ) .and.                   &
         (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.     &
         icodck.ne.6 .or. icodom.ne.6 ) ) then
      chaine=nomvar(ippkip)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,ik),icodcu,icodcv,icodcw
      chaine=nomvar(ippomg)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,iomg),icodcu,icodcv,icodcw
      nstukw = nstukw + 1
    endif

  enddo

elseif(iturb.eq.70 ) then

  do ifac = 1, nfabor

    icodcu = icodcl(ifac,iu)
    icodcv = icodcl(ifac,iv)
    icodcw = icodcl(ifac,iw)
    icodcn = icodcl(ifac,inusa)

    if( (icodcu.eq.5 .or. icodcv.eq.5 .or. icodcw.eq.5 .or.     &
         icodcn.eq.5 ) .and.                                    &
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.     &
         icodcn.ne.5 ) ) then
      chaine=nomvar(ippnup)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,inusa),icodcu,icodcv,icodcw
      nstunu = nstunu + 1
    endif

    if( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.     &
         icodcn.eq.6 ) .and.                                    &
         (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.     &
         icodcn.ne.6 ) ) then
      chaine=nomvar(ippnup)
      write(nfecra,1030)                                        &
           ifac,iprfml(ifmfbr(ifac),1),chaine(1:8),                &
           icodcl(ifac,inusa),icodcu,icodcv,icodcw
    endif

  enddo
endif

! 2.5 VERIFICATIONS DES COHERENCES INTER VARIABLES INTRA PHASE
! =============================================================

! --- Coherence vitesse scalaires

if( nscal.ge.1 ) then
  do iis = 1, nscal
    if(itytur.eq.2.or.itytur.eq.3) then
      ivar  = isca(iis)
      do ifac = 1, nfabor
        icodcu = icodcl(ifac,iu)
        if(icodcl(ifac,ivar).eq.5.and.icodcu.ne.5) then
          chaine=nomvar(ipprtp(ivar))
          write(nfecra,1050) ifac,iprfml(ifmfbr(ifac),1),         &
                        chaine(1:8), iis,                         &
                        icodcl(ifac,ivar), icodcu
          nstusc = nstusc + 1
        endif
      enddo
    endif
  enddo
endif

!===============================================================================
! 3.  IMPRESSIONS RECAPITULATIVES
!===============================================================================

iok = 0

if( nstoni.gt.0 .or. nstosc.gt.0 .or. nstovf.gt.0 .or.            &
                                      nstusc.gt.0 ) then
  write (nfecra,1901) nstoni, nstosc, nstovf, nstusc
  iok = 1
endif

if( nstvit.gt.0 .or. nstopp.gt.0 .or.             &
     nstoke.gt.0 .or. nstrij.gt.0 .or.             &
     nstov2.gt.0 .or. nstonu.gt.0 .or.             &
     nstuvw.gt.0 .or. nstoup.gt.0 .or.             &
     nstuke.gt.0 .or. nsurij.gt.0 .or.             &
     nstuv2.gt.0 .or. nstunu.gt.0       ) then
  write (nfecra,1902)        nstvit,nstopp,       &
       nstoke,nstrij,       &
       nstov2,nstonu,       &
       nstuvw,nstoup,       &
       nstuke,nsurij,       &
       nstuv2,nstunu
  iok = 1
endif

if(iok.ne.0) then
  call csexit (1)
  !==========
endif

!===============================================================================
! 3.  FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@ COND. LIM. NON INITIALISEES                                ',/,&
'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; VARIABLE ',A8    ,/,&
'@     ICODCL VARIABLE ', I10                                  ,/,&
'@                                                            '  )
 1010 format(                                                           &
'@                                                            ',/,&
'@ COND. LIM. NON PREVUES                                     ',/,&
'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; VARIABLE ',A8    ,/,&
'@     ICODCL VARIABLE ', I10                                  ,/,&
'@                                                            '  )
 1015 format(                                                           &
'@                                                            ',/,&
'@ CONDITIONS AUX LIMITES DE PAROI RUGUEUSE INCOMPATIBLES     ',/,&
'@ AVEC LE MODULE COMPRESSIBLE                                ',/,&
'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; VARIABLE ',A8    ,/,&
'@     ICODCL VARIABLE =',I10                                  ,/,&
'@     IPPMOD(ICOMPF)  =',I10                                  ,/,&
'@                                                            '  )
 1020 format(                                                           &
'@                                                            ',/,&
'@ INCOHERENCE COND. LIM. COMPOSANTES DE LA VITESSE           ',/,&
'@   FACE ',I10   ,'; PROPRIETE 1:',I10                        ,/,&
'@     ICODCL VITESSE  ',3I10                                  ,/,&
'@                                                            '  )
 1030 format(                                                           &
'@                                                            ',/,&
'@ INCOHERENCE COND. LIM. VITESSE-VARIABLE                    ',/,&
'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; VARIABLE ',A8    ,/,&
'@     ICODCL VARIABLE ', I10                                  ,/,&
'@     ICODCL VITESSE  ',3I10                                  ,/,&
'@                                                            '  )
 1040 format(                                                           &
'@                                                            ',/,&
'@ INCOHERENCE COND. LIM. VITESSE-RIJ-EPSILON                 ',/,&
'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; RIJ-EPSILON     ',/,&
'@     ICODCL RIJ-EPS ',7I5                                    ,/,&
'@     ICODCL VITESSE ',3I5                                    ,/,&
'@                                                            '  )
 1050 format(                                                           &
'@                                                            ',/,&
'@ INCOHERENCE COND. LIM. VITESSE-SCALAIRE                    ',/,&
'@   FACE ',I10   ,'; PROPRIETE 1:',I10   ,'; VARIABLE ',A8    ,/,&
'@     SCALAIRE NUMERO ',I10                                   ,/,&
'@     ICODCL SCALAIRE ',I10   ,'; ICODCL VITESSE ',I10        ,/,&
'@                                                            '  )
 1901 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@         Conditions aux limites non initialisees  : ',I10    ,/,&
'@         Conditions aux limites non prevues :               ',/,&
'@             sur les scalaires                    : ',I10    ,/,&
'@             sur les scalaires representant                 ',/,&
'@                                    une variance  : ',I10    ,/,&
'@         Incoherences :                                     ',/,&
'@             entre vitesse et scalaires           : ',I10    ,/,&
'@                                                            ',/,&
'@         Le calcul ne sera pas execute.                     ',/,&
'@                                                            ',/,&
'@         Verifier les parametres donnes via l''interface    ',/,&
'@           ou usclim.                                       ',/,&
'@                                                            ',/)
 1902 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@         Conditions aux limites non prevues :               ',/,&
'@             sur la vitesse                       : ',I10    ,/,&
'@             sur la pression                      : ',I10    ,/,&
'@             sur k et epsilon                     : ',I10    ,/,&
'@             sur Rij et epsilon                   : ',I10    ,/,&
'@             sur k, epsilon, phi et f_barre       : ',I10    ,/,&
'@             sur nu de Spalart Allmaras           : ',I10    ,/,&
'@         Incoherences :                                     ',/,&
'@             entre les composantes de la vitesse  : ',I10    ,/,&
'@             entre vitesse et pression            : ',I10    ,/,&
'@             entre vitesse et k-epsilon           : ',I10    ,/,&
'@             entre vitesse et Rij-epsilon         : ',I10    ,/,&
'@             entre vitesse et v2f                 : ',I10    ,/,&
'@             entre vitesse et nu                  : ',I10    ,/,&
'@                                                            ',/,&
'@         Le calcul ne sera pas execute.                     ',/,&
'@                                                            ',/,&
'@         Verifier les parametres donnes via l''interface    ',/,&
'@           ou usclim.                                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                           &
'@                                                            ',/,&
'@ UNINITIALIZED BOUNDARY CONDITIONS                          ',/,&
'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; VARIABLE ',A8     ,/,&
'@     ICODCL VARIABLE ', I10                                  ,/,&
'@                                                            '  )
 1010 format(                                                           &
'@                                                            ',/,&
'@ UNEXPECTED BOUNDARY CONDITIONS                             ',/,&
'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; VARIABLE ',A8     ,/,&
'@     ICODCL VARIABLE ', I10                                  ,/,&
'@                                                            '  )
 1015 format(                                                           &
'@                                                            ',/,&
'@ ROUGH WALL BOUNDARY CONDITIONS INCOMPATIBLE WITH THE       ',/,&
'@ COMPRESSIBLE MODULE                                        ',/,&
'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; VARIABLE ',A8     ,/,&
'@     ICODCL VARIABLE =',I10                                  ,/,&
'@     IPPMOD(ICOMPF)  =',I10                                  ,/,&
'@                                                            '  )
 1020 format(                                                           &
'@                                                            ',/,&
'@ INCOHERENCY BOUNDARY CONDITIONS VELOCITY COMPONENT         ',/,&
'@   FACE ',I10   ,'; PROPERTY 1:',I10                         ,/,&
'@     ICODCL VELOCITY ',3I10                                  ,/,&
'@                                                            '  )
 1030 format(                                                           &
'@                                                            ',/,&
'@ INCOHERENCY BOUNDARY CONDITIONS VELOCITY-VARIABLE          ',/,&
'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; VARIABLE ',A8     ,/,&
'@     ICODCL VARIABLE ', I10                                  ,/,&
'@     ICODCL VELOCITY ',3I10                                  ,/,&
'@                                                            '  )
 1040 format(                                                           &
'@                                                            ',/,&
'@ INCOHERENCY BOUNDARY CONDITIONS VELOCITY-RIJ-EPSILON       ',/,&
'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; RIJ-EPSILON      ',/,&
'@     ICODCL RIJ-EPS  ',7I5                                   ,/,&
'@     ICODCL VELOCITY ',3I5                                   ,/,&
'@                                                            '  )
 1050 format(                                                           &
'@                                                            ',/,&
'@ INCOHERENCY BOUNDARY CONDITIONS VELOCITY-SCALAR            ',/,&
'@   FACE ',I10   ,'; PROPERTY 1:',I10   ,'; VARIABLE ',A8     ,/,&
'@     SCALAR NUMBER ',I10                                     ,/,&
'@     ICODCL SCALAR ',I10   ,'; ICODCL VELOCITY ',I10         ,/,&
'@                                                            '  )
 1901 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE BOUNDARY CONDITIONS VERIF.    ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@         Uninitialized boundary conditions        : ',I10    ,/,&
'@         Unexpected  boundary conditions:                   ',/,&
'@             on the scalars                       : ',I10    ,/,&
'@             on the scalars representing                    ',/,&
'@                                      a variance  : ',I10    ,/,&
'@         Incoherencies:                                     ',/,&
'@             between velocity and scalars         : ',I10    ,/,&
'@                                                            ',/,&
'@         The calculation will not be run.                   ',/,&
'@                                                            ',/,&
'@         Verify the parameters given via the interface or   ',/,&
'@           usclim.                                          ',/,&
'@                                                            ',/)
 1902 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE BOUNDARY CONDITIONS VERIF.    ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@         Unexpeted boundary conditions:                     ',/,&
'@             on the velocity                      : ',I10    ,/,&
'@             on the pressure                      : ',I10    ,/,&
'@             on k and epsilon                     : ',I10    ,/,&
'@             on Rij and epsilon                   : ',I10    ,/,&
'@             on k, epsilon, phi and f_barre       : ',I10    ,/,&
'@             on nu of Spalart Allmaras model      : ',I10    ,/,&
'@         Incoherencies:                                     ',/,&
'@             between the velocity components      : ',I10    ,/,&
'@             between velocity and pressure        : ',I10    ,/,&
'@             between velocity and k-epsilon       : ',I10    ,/,&
'@             between velocity and Rij-epsilon     : ',I10    ,/,&
'@             between velocity and v2f             : ',I10    ,/,&
'@             between velocity and nu              : ',I10    ,/,&
'@                                                            ',/,&
'@         The calculation will not be run.                   ',/,&
'@                                                            ',/,&
'@         Verify the parameters given via the interface or   ',/,&
'@           usclim.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

return
end subroutine
