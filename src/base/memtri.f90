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

subroutine memtri &
!================

 ( idebra ,                                                       &
   nvar   , nscal  ,                                              &
   nproce ,                                                       &
   idt    , itpuco , irtp   , irtpa  , ipropc ,                   &
   ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE VARIABLES NON GEOMETRIQUES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nproce           ! e  ! <-- ! nombre de prop phy aux centres                 !
! idt              ! e  ! --> ! "pointeur" sur dt                              !
! itpuco           ! e  ! --> ! "pointeur" sur tpucou                          !
! irtp, irtpa      ! e  ! --> ! "pointeur" sur rtp, rtpa                       !
! ipropc           ! e  ! --> ! "pointeur" sur propce                          !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use cstphy
use numvar
use entsor
use pointe
use albase
use period
use ppppar
use ppthch
use ppincl
use cfpoin
use lagpar
use lagdim
use lagran
use ihmpre
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          idebra
integer          nvar   , nscal
integer          nproce
integer          idt    , itpuco
integer          irtp   , irtpa
integer          ipropc
integer          ifinra

! Local variables

integer          iis, ippu, ippv, ippw, ivar, iprop
integer          imom, idtnm
integer          iipuco

!===============================================================================


!===============================================================================
! 1. INITIALISATION
!===============================================================================


!===============================================================================
! 2. PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINRA
!===============================================================================

! Work arrays "tpucou"
iipuco = 0
if (ipucou.eq.1 .or. ncpdct.gt.0) then
  iipuco = 1
endif

! Allocate main real arrays
irtp   = idebra
irtpa  = irtp   + ncelet *nvar
ipropc = irtpa  + ncelet *nvar
idt    = ipropc + ncelet *nproce
itpuco = idt    + ncelet
ifinra = itpuco + ncelet *ndim*iipuco

!===============================================================================
! 3. CORRESPONDANCE POUR POST-TRAITEMENT
!===============================================================================

! --> Correspondance IPP2RA pour post-process
!       Variables de calcul et proprietes physiques

do iis = 1 , nvppmx
  ipp2ra(iis) = 1
enddo

!     IPPROC a ete complete au prealable dans VARPOS

ivar = ipr
ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
ivar = iu
ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
ivar = iv
ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
ivar = iw
ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet

if    (itytur.eq.2) then
  ivar = ik
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = iep
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
elseif(itytur.eq.3) then
  ivar = ir11
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = ir22
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = ir33
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = ir12
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = ir13
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = ir23
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = iep
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ! Alpha is a new variable for the EBRSM
  if (iturb.eq.32) then
    ivar = ial
    ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  endif
elseif(itytur.eq.5) then
  ivar = ik
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = iep
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = iphi
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  if(iturb.eq.50) then
    ivar = ifb
    ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  elseif(iturb.eq.51) then
    ivar = ial
    ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  endif
elseif(iturb.eq.60) then
  ivar = ik
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = iomg
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
elseif(iturb.eq.70) then
  ivar = inusa
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
endif

if (iale.eq.1) then
  ivar = iuma
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = ivma
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = iwma
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
endif

!     Le choix fait dans VARPOS indique qu'on ne s'interessera
!       qu'aux proprietes au centre des cellules (pas au flux
!       de masse en particulier, ni a la masse volumique au bord)

do iprop = 1, nproce
  ipp2ra(ipppro(iprop)) = ipropc+(iprop-1)*ncelet
enddo

!     Pour les moments, on repere dans IPPMOM le mode de division par le temps
!       = 0 : pas de division
!       > 0 : IPPMOM donne le pointeur dans RA sur le DT cumule
!                                                   (tableau NCEL dans PROPCE)
!       < 0 : IPPMOM donne le rang dans DTCMOM du DT cumule (uniforme)
do iprop = 1, nvppmx
  ippmom(iprop) = 0
enddo
do imom = 1, nbmomt
!       Pointeur iprop des moments pour IPP2RA(IPPPRO(IPROP)) et IPPMOM(IPPPRO(IPROP))
  iprop = ipproc(icmome(imom))
!       Type de DT cumule et numero
  idtnm = idtmom(imom)
  if(idtnm.gt.0) then
    ippmom(ipppro(iprop)) = ipropc+(ipproc(icdtmo(idtnm))-1)*ncelet
  elseif(idtnm.lt.0) then
    ippmom(ipppro(iprop)) = idtnm
  endif
enddo

do iis = 1 , nscal
  ivar = isca  (iis  )
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
enddo

if (idtvar.le.0) then
  ipp2ra(ippdt ) = 1
else
  ipp2ra(ippdt ) = idt
endif

!     Couplage instationnaire vitesse/pression
if (ipucou.eq.0) then
  ipp2ra(ipptx)= 1
  ipp2ra(ippty)= 1
  ipp2ra(ipptz)= 1
else
  ipp2ra(ipptx)= itpuco
  ipp2ra(ippty)= itpuco+ncelet
  ipp2ra(ipptz)= itpuco+2*ncelet
endif

!     Vecteur vitesse chrono
ippu = ipprtp(iu)
ippv = ipprtp(iv)
ippw = ipprtp(iw)
if(ichrvr(ippu).eq.1.and.ichrvr(ippv).eq.1.and.                 &
     ichrvr(ippw).eq.1) then
  ichrvr(ippv) = 0
  ichrvr(ippw) = 0
  ipp2ra(ippu) = - ipp2ra(ippu)
endif
!     Vecteur vitesse de maillage chrono
if (iale.eq.1) then
  ippu = ipprtp(iuma)
  ippv = ipprtp(ivma)
  ippw = ipprtp(iwma)
  if(ichrvr(ippu).eq.1.and.ichrvr(ippv).eq.1.and.ichrvr(ippw).eq.1) then
    ichrvr(ippv) = 0
    ichrvr(ippw) = 0
    ipp2ra(ippu) = - ipp2ra(ippu)
  endif
endif
!     Potentiel vecteur chrono
if(ippmod(ielarc).ge.2) then
  ippu = ipprtp(isca(ipotva(1)))
  ippv = ipprtp(isca(ipotva(2)))
  ippw = ipprtp(isca(ipotva(3)))
  if(ichrvr(ippu).eq.1.and.ichrvr(ippv).eq.1.and.ichrvr(ippw).eq.1) then
    ichrvr(ippv) = 0
    ichrvr(ippw) = 0
    ipp2ra(ippu) = - ipp2ra(ippu)
  endif
endif
!     Laplace vecteur chrono
if(ippmod(ielarc).ge.1) then
  ippu = ipppro(ipproc(ilapla(1)))
  ippv = ipppro(ipproc(ilapla(2)))
  ippw = ipppro(ipproc(ilapla(3)))
  if(ichrvr(ippu).eq.1.and.ichrvr(ippv).eq.1.and.ichrvr(ippw).eq.1) then
    ichrvr(ippv) = 0
    ichrvr(ippw) = 0
    ipp2ra(ippu) = - ipp2ra(ippu)
  endif
endif

return
end subroutine
