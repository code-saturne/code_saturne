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

subroutine memtri &
!================

 ( idbia0 , idbra0 , iverif ,                                     &
   nvar   , nscal  ,                                                                                 &
   ncofab , nproce , nprofa , nprofb ,                            &
   iisstd , ifrcx  ,                                              &
   idt    , irtp   , irtpa  , ipropc , ipropf , ipropb ,          &
   icoefa , icoefb ,                                              &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE VARIABLES NON GEOMETRIQUES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! iverif           ! e  ! <-- ! indicateur des tests elementaires              !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncofab           ! e  ! <-- ! nombre de couple de cl a prevoir               !
! nproce           ! e  ! <-- ! nombre de prop phy aux centres                 !
! nprofa           ! e  ! <-- ! nombre de prop phy aux faces internes          !
! nprofb           ! e  ! <-- ! nombre de prop phy aux faces de bord           !
! iisstd           ! e  ! --> ! "pointeur" sur isostd(reperage sortie          !
! idt              ! e  ! --> ! "pointeur" sur dt                              !
! irtp, irtpa      ! e  ! --> ! "pointeur" sur rtp, rtpa                       !
! ipropc           ! e  ! --> ! "pointeur" sur propce                          !
! ipropf           ! e  ! --> ! "pointeur" sur propfa                          !
! ipropb           ! e  ! --> ! "pointeur" sur propfb                          !
! icoefa, b        ! e  ! --> ! "pointeur" sur coefa, coefb                    !
! ifrcx            ! e  ! --> ! "pointeur" sur frcxt                           !
! ifinia           ! i  ! --> ! number of first free position in ia (at exit)  !
! ifinra           ! i  ! --> ! number of first free position in ra (at exit)  !
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
use dimens, only: ndimfb
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

integer          idbia0 , idbra0
integer          iverif
integer          nvar   , nscal
integer          ncofab , nproce , nprofa , nprofb
integer          iisstd , ifrcx
integer          idt
integer          irtp   , irtpa
integer          ipropc , ipropf , ipropb
integer          icoefa , icoefb
integer          ifinia , ifinra

! Local variables

integer          idebia , idebra
integer          iis, ippu, ippv, ippw, ivar, iprop
integer          iipero , iiirij , imom, idtnm
integer          idpar1 , idpar2 , iypar1, iiyplb, iiforb, iicoci

!===============================================================================


!===============================================================================
! 1. INITIALISATION
!===============================================================================


idebia = idbia0
idebra = idbra0





!===============================================================================
! 2. PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
!===============================================================================

! --> Remarques :

!     IPUCOU = 1 ne depend pas de la phase

!     NCOFAB, NPROCE, NPROFA et NPROFB ont ete determines dans VARPOS
!         et ne servent en tant que dimensions que dans le present
!         sous programme. On pourrait les passer en common dans numvar.h

!     ITYPFB, ITRIFB et ISYMPA peuvent passer en entier dans certains
!         sous-pgm, il convient donc qu'ils soient en un seul bloc.

!     Le tableau des zones frontieres des faces de bord pour les
!         physiques particulieres est de declare ci-dessous (voir PPCLIM)


! --> Preparations :

!     On regarde s'il existe au moins une periodicite de rotation
!                                 une phase avec Rij

iipero = 0
if (iperot.gt.0) then
  iipero = 1
endif
iiirij = 0
if(itytur.eq.3) then
  iiirij = 1
endif


!     Distance a la paroi

!       On reserve ici idipar (distance a la paroi) : c'est oblige
!       On reserve aussi iyppar (yplus) : on pourrait s'en passer
!         et en faire un tableau local a reserver dans memdyp
!         mais ca facilite la visualisation et l'initialisation
!         de yplus pour son calcul a tous les pas de temps

!     Distance a la paroi, tableau selon le mode de calcul
idpar1 = 0
idpar2 = 0
!     Calcul par eq de diffusion
if(ineedy.eq.1.and.abs(icdpar).eq.1) then
  idpar1 = 1
endif
!     Calcul direct
if(ineedy.eq.1.and.abs(icdpar).eq.2) then
  idpar2 = 1
endif
!     Yplus associe (calcul par mode de diffusion et LES+VanDriest
!       ou lagrangien+IROULE=2)
iypar1 = 0
if(ineedy.eq.1.and.abs(icdpar).eq.1) then
  if(itytur.eq.4) then
    if(idries.eq.1) then
      iypar1 = 1
    endif
  endif
  if (iilagr.ge.1 .and. iroule.eq.2) iypar1 = 1
endif

!     Stockage supplémentaire si on initialise le gradient
!       par moindre carre

iicoci = 0
if(imrgra.eq.4 .or. iverif.eq.1) then
  iicoci = 1
endif


!     Post-traitements particuliers  (faces de bord)


iiyplb = 0
!     Yplus au bord
if(mod(ipstdv,ipstyp).eq.0) then
  iiyplb = 1
endif
!     Efforts aux bords
iiforb = 0
if(ineedf.eq.1) then
  iiforb = 1
endif


! --> Reservation de memoire entiere


iitypf = idebia
iitrif = iitypf + nfabor
iisymp = iitrif + nfabor
ifinia = iisymp + nfabor
if(idpar2.eq.1) then
  iifapa = ifinia
  ifinia        = iifapa + ncelet
else
  !         cette valeur nulle est utilisee dans les tests
  iifapa = 0
endif

!  Zones de face de bord : on utilise provisoirement les zones des physiques
!    particulieres, meme sans physique particuliere
!    -> sera modifie lors de la restructuration des zones de bord
iizfpp = ifinia
if(ippmod(iphpar).ge.1 .or. iihmpr.eq.1) then
  ifinia = iizfpp + nfabor
else
  ifinia = iizfpp
endif

iisstd = ifinia
ifinia = iisstd + (nfabor+1)*iphydr

if(ippmod(icompf).ge.0) then
  iifbet = ifinia
  iifbru = iifbet + nfabor
  ifinia = iifbru + nfabor
else
  iifbet = 0
  iifbru = 0
endif

! --> Reservation de memoire reelle

icoefa = idebra
icoefb = icoefa + ndimfb *ncofab
irtp   = icoefb + ndimfb *ncofab
irtpa  = irtp   + ncelet *nvar
ipropc = irtpa  + ncelet *nvar
ipropf = ipropc + ncelet *nproce
ipropb = ipropf + nfac   *nprofa
idt    = ipropb + ndimfb *nprofb
icocg  = idt    + ncelet
icocgb = icocg  + ncelet *9
icoci  = icocgb + ncelbr *9
icocib = icoci  + ncelet *9 * iicoci
itpuco = icocib + ncelbr *9 * iicoci
idipar = itpuco + ncelet *ndim*ipucou
iyppar = idipar + ncelet *idpar1
idudxy = iyppar + ncelet *iypar1
iwdudx = idudxy + (ncelet-ncel) * 3 * 3 * iipero
idrdxy = iwdudx + (ncelet-ncel) * 3 * 3 * iipero
iwdrdx = idrdxy + (ncelet-ncel) * 6 * 3 * iipero*iiirij
ifrcx  = iwdrdx + (ncelet-ncel) * 6 * 3 * iipero*iiirij
iyplbr = ifrcx  + ncelet*ndim*iphydr
iforbr = iyplbr + nfabor*iiyplb
ifinra = iforbr + nfabor*ndim*iiforb

!     On rajoute des tableaux pour le k-omega SST si necessaire
!     En k-omega, on a besoin de calculer 2 Sij.Sij pour etre utilise
!     dans PHYVAR et dans TURBKW. On reserve un tableau pour divU en meme temps.

is2kw  = ifinra
idvukw = ifinra
if (iturb.eq.60) then
  idvukw = is2kw  + ncelet
  ifinra        = idvukw + ncelet
endif

! En ALE ou maillage mobile, on reserve des tableaux supplementaires
! de position initiale
if (iale.eq.1.or.imobil.eq.1) then
  ixyzn0 = ifinra
  ifinra = ixyzn0 + ndim*nnod
else
  ixyzn0 = 0
endif

! En ALE, on reserve des tableaux supplementaires
! de deplacement et de type de faces de bord
if (iale.eq.1) then
  iimpal = ifinia
  iialty = iimpal + nnod
  ifinia = iialty + nfabor

  idepal = ifinra
  ifinra = idepal + ndim*nnod
else
  iimpal = 0
  iialty = 0
  idepal = 0
endif

! --> Verification

call iasize('memtri',ifinia)
!==========

call rasize('memtri',ifinra)
!==========


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
elseif(iturb.eq.50) then
  ivar = ik
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = iep
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = iphi
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
  ivar = ifb
  ipp2ra(ipprtp(ivar)) = irtp  +(ivar-1)*ncelet
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
    ippmom(ipppro(iprop)) =                                       &
         ipropc+(ipproc(icdtmo(idtnm))-1)*ncelet
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
  if(ichrvr(ippu).eq.1.and.ichrvr(ippv).eq.1.and.                 &
    ichrvr(ippw).eq.1) then
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
  if(ichrvr(ippu).eq.1.and.ichrvr(ippv).eq.1.and.                 &
                           ichrvr(ippw).eq.1) then
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
  if(ichrvr(ippu).eq.1.and.ichrvr(ippv).eq.1.and.                 &
                           ichrvr(ippw).eq.1) then
    ichrvr(ippv) = 0
    ichrvr(ippw) = 0
    ipp2ra(ippu) = - ipp2ra(ippu)
  endif
endif

return
end subroutine
