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

subroutine cscpfb &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nptdis , ityloc , nvcp   , numcpl , nvcpto,                    &
   locpts ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   coopts , djppts , pndpts ,                                     &
   rvdis  , dofpts ,                                              &
   ra     )

!===============================================================================
! FONCTION :
! --------

! PREPARATION DE L'ENVOI DES VARIABLES POUR UN COUPLAGE
!   ENTRE DEUX INSTANCES DE CODE_SATURNE VIA LES FACES DE BORD

! L'INFORMATION RECUE SERA TRANSFORMEE EN CONDITION LIMITE DANS
!   LA SUBROUTINE CSC2CL

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nptdis , nvcp   , numcpl , nvcpto , ityloc

integer          locpts(nptdis)
integer          ia(*)


double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision coopts(3,nptdis), djppts(3,nptdis)
double precision pndpts(nptdis), dofpts(3,nptdis)
double precision rvdis(nptdis,nvcpto)
double precision ra(*)

! Local variables

integer          idebia , idebra , ifinia , ifinra

integer          itrav1 , itrav2 , itrav3 , itrav4 , itrav5
integer          itrav6 , itrav7 , itrav8 , itrav

integer          ipt    , ifac   , iel    , isou
integer          ivar   , iscal  , ipcrom
integer          ipriph , iuiph  , iviph  , iwiph
integer          ikiph  , iepiph , ifbiph , iomiph , iphiph
integer          ir11ip , ir22ip , ir33ip , ir12ip , ir13ip , ir23ip
integer          itravx , itravy , itravz
integer          igradx , igrady , igradz
integer          inc    , iccocg , iphydp , iclvar, nswrgp
integer          iwarnp , imligp
integer          ipos
integer          itytu0

double precision epsrgp , climgp , extrap
double precision xjjp   , yjjp   , zjjp
double precision d2s3
double precision xjpf,yjpf,zjpf,jpf
double precision xx, yy, zz
double precision omegal(3), omegad(3), omegar(3), omgnrl, omgnrd, omgnrr
double precision vitent, daxis2

!===============================================================================

!=========================================================================
! 1.  INITIALISATIONS
!=========================================================================

! Initialize variables to avoid compiler warnings

itrav = 0

vitent = 0.d0

! Memoire

idebia = idbia0
idebra = idbra0

itrav1 = idebra
itrav2 = itrav1 + nptdis
itrav3 = itrav2 + nptdis
itrav4 = itrav3 + nptdis
itrav5 = itrav4 + nptdis
itrav6 = itrav5 + nptdis
itrav7 = itrav6 + nptdis
itrav8 = itrav7 + nptdis
ifinra = itrav8 + nptdis

call rasize('cscpfb',ifinra)


d2s3 = 2.d0/3.d0

ipcrom = ipproc(irom)

if (icormx(numcpl).eq.1) then

  ! On récupère dans tous les cas le vecteur rotation de l'autre instance
  omegal(1) = omegax
  omegal(2) = omegay
  omegal(3) = omegaz
  call tbrcpl(numcpl,3,3,omegal,omegad)

  ! Vecteur vitesse relatif d'une instance a l'autre
  omegar(1) = omegal(1) - omegad(1)
  omegar(2) = omegal(2) - omegad(2)
  omegar(3) = omegal(3) - omegad(3)

  omgnrl = sqrt(omegal(1)**2 + omegal(2)**2 + omegal(3)**2)
  omgnrd = sqrt(omegad(1)**2 + omegad(2)**2 + omegad(3)**2)
  omgnrr = sqrt(omegar(1)**2 + omegar(2)**2 + omegar(3)**2)

else

  omegal(1) = 0.d0
  omegal(2) = 0.d0
  omegal(3) = 0.d0

  omegar(1) = 0.d0
  omegar(2) = 0.d0
  omegar(3) = 0.d0

  omgnrl = 0.d0
  omgnrd = 0.d0
  omgnrr = 0.d0

endif

! On part du principe que l'on envoie les bonnes variables à
! l'instance distante et uniquement celles-là.

! De plus, les variables sont envoyées dans l'ordre de VARPOS :

!     - pression (unique pour toute les phases)
!     - vitesse
!     - grandeurs turbulentes (selon le modèle)
!   Et ensuite :
!     - scalaires physique particulière (pas encore traités)
!     - scalaires utilisateur
!     - vitesse de maillage (non couplée, donc non envoyée)


ipos = 1

!=========================================================================
! 1.  PREPARATION DE LA PRESSION
!=========================================================================

! La pression est unique pour chaque phase, inutile donc de
! le faire pour toutes les phases.

ipriph = ipr

! --- Calcul du gradient de la pression pour interpolation

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtp(1,ipriph))
  !==========
endif

inc    = 1
iccocg = 1
iphydp = 0
iclvar = iclrtp(ipriph,icoef)
nswrgp = nswrgr(ipriph)
imligp = imligr(ipriph)
iwarnp = iwarni(ipriph)
epsrgp = epsrgr(ipriph)
climgp = climgr(ipriph)
extrap = extrag(ipriph)

call grdcel                                                   &
!==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    ipriph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,ipriph) , coefa(1,iclvar) , coefb(1,iclvar) ,           &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )

! For a specific face to face coupling, geometric assumptions are made

if (ifaccp.eq.1) then

  do ipt = 1, nptdis

    iel = locpts(ipt)

    ! --- Pour la pression on veut imposer un dirichlet tel que le gradient
    !     de pression se conserve entre les deux domaines couplés Pour cela
    !     on impose une interpolation centrée

    xjjp = djppts(1,ipt)
    yjjp = djppts(2,ipt)
    zjjp = djppts(3,ipt)

    rvdis(ipt,ipos) = rtp(iel,ipriph) &
         + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    ! On prend en compte le potentiel centrifuge en repère relatif
    if (icormx(numcpl).eq.1) then

      ! Calcul de la distance a l'axe de rotation
      ! On suppose que les axes sont confondus...

      xx = xyzcen(1,iel) + xjjp
      yy = xyzcen(2,iel) + yjjp
      zz = xyzcen(3,iel) + zjjp

      daxis2 =   (omegar(2)*zz - omegar(3)*yy)**2 &
           + (omegar(3)*xx - omegar(1)*zz)**2 &
           + (omegar(1)*yy - omegar(2)*xx)**2

      daxis2 = daxis2 / omgnrr**2

      rvdis(ipt,ipos) = rvdis(ipt,ipos)                         &
           + 0.5d0*propce(iel,ipcrom)*(omgnrl**2 - omgnrd**2)*daxis2

    endif

  enddo

  ! For a generic coupling, no assumption can be made

else

  do ipt = 1, nptdis

    iel = locpts(ipt)

    xjpf = coopts(1,ipt) - xyzcen(1,iel)- djppts(1,ipt)
    yjpf = coopts(2,ipt) - xyzcen(2,iel)- djppts(2,ipt)
    zjpf = coopts(3,ipt) - xyzcen(3,iel)- djppts(3,ipt)

    if(pndpts(ipt).ge.0.d0.and.pndpts(ipt).le.1.d0) then
      jpf = -1.d0*sqrt(xjpf**2+yjpf**2+zjpf**2)
    else
      jpf =       sqrt(xjpf**2+yjpf**2+zjpf**2)
    endif

    rvdis(ipt,ipos) = (xjpf*w1(iel)+yjpf*w2(iel)+zjpf*w3(iel))  &
         /jpf

  enddo

endif
!       FIn pour la pression


!=========================================================================
! 2.  PREPARATION DE LA VITESSE
!=========================================================================

iuiph = iu
iviph = iv
iwiph = iw

! --- Calcul du gradient de la vitesse pour interpolation

if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(rtp(1,iuiph), rtp(1,iviph), rtp(1,iwiph))
  !==========
endif

do isou = 1, 3

  ipos = ipos + 1

  if(isou.eq.1) ivar = iuiph
  if(isou.eq.2) ivar = iviph
  if(isou.eq.3) ivar = iwiph

  inc    = 1
  iccocg = 1
  iphydp = 0
  iclvar = iclrtp(ivar,icoef)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  call grdcel                                                   &
  !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,             &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

! --- Pour la vitesse on veut imposer un dirichlet de vitesse qui "imite"
!     ce qui se passe pour une face interne. On se donne le choix entre
!     UPWIND, SOLU et CENTRE (parties commentées selon le choix retenu).
!     Pour l'instant seul le CENTRE respecte ce qui se passerait pour la
!     diffusion si on avait un seul domaine

      ! -- UPWIND

      !        xjjp = djppts(1,ipt)
      !        yjjp = djppts(2,ipt)
      !        zjjp = djppts(3,ipt)

      !        rvdis(ipt,ipos) = rtp(iel,ivar)

      ! -- SOLU

      !        xjf = coopts(1,ipt) - xyzcen(1,iel)
      !        yjf = coopts(2,ipt) - xyzcen(2,iel)
      !        zjf = coopts(3,ipt) - xyzcen(3,iel)

      !        rvdis(ipt,ipos) = rtp(iel,ivar) &
      !          + xjf*w1(iel) + yjf*w2(iel) + zjf*W3(iel)

      ! -- CENTRE

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      rvdis(ipt,ipos) = rtp(iel,ivar) &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

      ! On prend en compte la vitesse d'entrainement en repère relatif
      if (icormx(numcpl).eq.1) then

        if (isou.eq.1) then
          vitent =   omegar(2)*(xyzcen(3,iel)+zjjp) &
               - omegar(3)*(xyzcen(2,iel)+yjjp)
        elseif (isou.eq.2) then
          vitent =   omegar(3)*(xyzcen(1,iel)+xjjp) &
               - omegar(1)*(xyzcen(3,iel)+zjjp)
        elseif (isou.eq.3) then
          vitent =   omegar(1)*(xyzcen(2,iel)+yjjp) &
               - omegar(2)*(xyzcen(1,iel)+xjjp)
        endif

        rvdis(ipt,ipos) = rvdis(ipt,ipos) + vitent

      endif

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = dofpts(1,ipt) + djppts(1,ipt)
      yjjp = dofpts(2,ipt) + djppts(2,ipt)
      zjjp = dofpts(3,ipt) + djppts(3,ipt)


      rvdis(ipt,ipos) = rtp(iel,ivar)                             &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

  endif

enddo
!       Fin de la boucle sur les composantes de la vitesse


!=========================================================================
! 3.  PREPARATION DES GRANDEURS TURBULENTES
!=========================================================================

itytu0 = iturcp(numcpl)/10


!=========================================================================
!       3.1 Turbulence dans l'instance locale : modèles k-epsilon
!=========================================================================

if (itytur.eq.2) then

!=======================================================================
!          3.1.1. INTERPOLATION EN J'
!=======================================================================

!         Préparation des données: interpolation de k en J'

  ikiph  = ik
  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ikiph))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iphydp = 0
  iclvar = iclrtp(ikiph,icoef)
  nswrgp = nswrgr(ikiph)
  imligp = imligr(ikiph)
  iwarnp = iwarni(ikiph)
  epsrgp = epsrgr(ikiph)
  climgp = climgr(ikiph)
  extrap = extrag(ikiph)

  call grdcel                                                   &
  !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    ikiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp,  &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,ikiph) , coefa(1,iclvar) , coefb(1,iclvar) ,            &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav1 + ipt-1) = rtp(iel,ikiph)                         &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav1 + ipt-1) = rtp(iel,ikiph)                         &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

  endif

  !         Préparation des données: interpolation de epsilon en J'

  iepiph  = iep

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,iepiph))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iphydp = 0
  iclvar = iclrtp(iepiph,icoef)
  nswrgp = nswrgr(iepiph)
  imligp = imligr(iepiph)
  iwarnp = iwarni(iepiph)
  epsrgp = epsrgr(iepiph)
  climgp = climgr(iepiph)
  extrap = extrag(iepiph)

  call grdcel                                                   &
  !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    iepiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,iepiph) , coefa(1,iclvar) , coefb(1,iclvar) ,           &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav2 + ipt-1) = rtp(iel,iepiph)                        &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav2 + ipt-1) = rtp(iel,iepiph)                        &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

  endif


!=======================================================================
!          3.1.2.   Transfert de variable à "iso-modèle"
!=======================================================================

  if (itytu0.eq.2) then

    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1
    ikiph  = ik

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav1 + ipt-1)
    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1
    iepiph = iep

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav2 + ipt-1)
    enddo

    !=======================================================================
    !          3.1.3.   Transfert de k-eps vers Rij-eps
    !=======================================================================

  elseif (itytu0.eq.3) then

    !           Tenseur Rij
    !           ------------
    !           Termes de la diagonal R11,R22,R33

    do isou =1, 3

      ipos = ipos + 1

      do ipt = 1, nptdis
        rvdis(ipt,ipos) = d2s3*ra(itrav1 + ipt-1)
      enddo

    enddo

    !           Termes R12,R13,R23

    iuiph = iu
    iviph = iv
    iwiph = iw

    !           La synchronisation des halos a deja ete faite plus haut

    do isou = 1, 3

      if(isou.eq.1) ivar = iuiph
      if(isou.eq.2) ivar = iviph
      if(isou.eq.3) ivar = iwiph

      inc    = 1
      iccocg = 1
      iphydp = 0
      iclvar = iclrtp(ivar,icoef)
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)

      call grdcel                                               &
      !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,             &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


      do ipt = 1, nptdis

        iel = locpts(ipt)

        if(isou.eq.1) then
          ra(itrav3 + ipt-1) = w2(iel)
          ra(itrav4 + ipt-1) = w3(iel)
        elseif(isou.eq.2) then
          ra(itrav5 + ipt-1) = w1(iel)
          ra(itrav6 + ipt-1) = w3(iel)
        elseif(isou.eq.3) then
          ra(itrav7 + ipt-1) = w1(iel)
          ra(itrav8 + ipt-1) = w2(iel)
        endif

      enddo

    enddo
    !           Fin de la boucle sur les composantes de la vitesse

    !           R12
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*ra(itrav1 + ipt-1)**2*cmu        &
           /max(1.0d-10, ra(itrav2 + ipt-1))                       &
           *0.5d0*(ra(itrav3 + ipt-1) + ra(itrav5 + ipt-1))
    enddo

    !           R13
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*ra(itrav1 + ipt-1)**2*cmu        &
           /max(1.0d-10, ra(itrav2 + ipt-1))                       &
           *0.5d0*(ra(itrav4 + ipt-1) + ra(itrav7 + ipt-1))
    enddo

    !           R23
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*ra(itrav1 + ipt-1)**2*cmu        &
           /max(1.0d-10,ra(itrav2 + ipt-1))                        &
           *0.5d0*(ra(itrav6 + ipt-1) + ra(itrav8 + ipt-1))
    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1
    iepiph = iep

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav2 + ipt -1)
    enddo

    !=======================================================================
    !          3.1.4.   Transfert de k-eps vers v2f
    !=======================================================================

  elseif (iturcp(numcpl).eq.50) then

    !   ATTENTION: CAS NON PRIS EN COMPTE (ARRET DU CALCUL DANS CSCINI.F)

    !=======================================================================
    !          3.1.5.   Transfert de k-eps vers k-omega
    !=======================================================================

  elseif (iturcp(numcpl).eq.60) then

    !           Energie turbulente
    !           -----------------
    ipos = ipos + 1
    ikiph  = ik

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav1 + ipt -1)
    enddo

    !           Omega
    !           -----
    ipos = ipos + 1
    iomiph = iomg

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav2 + ipt-1)/cmu                  &
           /max(1.0d-10, ra(itrav1 + ipt-1))
    enddo


  endif

  !=========================================================================
  !       3.2 Turbulence dans l'instance locale : modèle Rij-epsilon
  !=========================================================================

elseif (itytur.eq.3) then

  !=======================================================================
  !          3.2.1. INTERPOLATION EN J'
  !=======================================================================

  !         Préparation des données: interpolation des Rij en J'

  if (irangp.ge.0.or.iperio.eq.1) then
    call synten &
    !==========
  ( rtp(1,ir11), rtp(1,ir12), rtp(1,ir13),  &
    rtp(1,ir12), rtp(1,ir22), rtp(1,ir23),  &
    rtp(1,ir13), rtp(1,ir23), rtp(1,ir33) )
  endif

  do isou = 1, 6

    if (isou.eq.1) ivar = ir11
    if (isou.eq.2) ivar = ir22
    if (isou.eq.3) ivar = ir33
    if (isou.eq.4) ivar = ir12
    if (isou.eq.5) ivar = ir13
    if (isou.eq.6) ivar = ir23

    inc    = 1
    iccocg = 1
    iphydp = 0
    iclvar = iclrtp(ivar,icoef)
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    call grdcel                                                 &
    !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,             &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )

    if (isou.eq.1) itrav = itrav1
    if (isou.eq.2) itrav = itrav2
    if (isou.eq.3) itrav = itrav3
    if (isou.eq.4) itrav = itrav4
    if (isou.eq.5) itrav = itrav5
    if (isou.eq.6) itrav = itrav6

    ! For a specific face to face coupling, geometric assumptions are made

    if (ifaccp.eq.1) then

      do ipt = 1, nptdis

        iel = locpts(ipt)

        xjjp = djppts(1,ipt)
        yjjp = djppts(2,ipt)
        zjjp = djppts(3,ipt)

        ra(itrav + ipt-1) = rtp(iel,ivar)                         &
             + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

      enddo

      ! For a generic coupling, no assumption can be made

    else

      do ipt = 1, nptdis

        iel = locpts(ipt)

        xjjp = djppts(1,ipt)
        yjjp = djppts(2,ipt)
        zjjp = djppts(3,ipt)

        ra(itrav + ipt-1) = rtp(iel,ivar)                         &
             + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

      enddo

    endif

  enddo

  !         Préparation des données: interpolation de epsilon en J'

  iepiph  = iep

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,iepiph))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iphydp = 0
  iclvar = iclrtp(iepiph,icoef)
  nswrgp = nswrgr(iepiph)
  imligp = imligr(iepiph)
  iwarnp = iwarni(iepiph)
  epsrgp = epsrgr(iepiph)
  climgp = climgr(iepiph)
  extrap = extrag(iepiph)

  call grdcel                                                   &
  !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    iepiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,iepiph) , coefa(1,iclvar) , coefb(1,iclvar) ,           &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav7 + ipt-1) = rtp(iel,iepiph)                        &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav7 + ipt-1) = rtp(iel,iepiph)                        &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

  endif

!=======================================================================
!          3.2.2. Transfert de variable à "iso-modèle"
!=======================================================================

  if (itytu0.eq.3) then

    !           Tensions de Reynolds
    !           --------------------
    do isou = 1, 6

      ipos = ipos + 1

      if (isou.eq.1) itrav = itrav1
      if (isou.eq.2) itrav = itrav2
      if (isou.eq.3) itrav = itrav3
      if (isou.eq.4) itrav = itrav4
      if (isou.eq.5) itrav = itrav5
      if (isou.eq.6) itrav = itrav6

      do ipt = 1, nptdis
        rvdis(ipt,ipos) = ra(itrav + ipt-1)
      enddo

    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1
    iepiph = iep

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav7 + ipt-1)
    enddo

    !=======================================================================
    !          3.2.3. Transfert de Rij-epsilon vers k-epsilon
    !=======================================================================

  elseif (itytu0.eq.2) then

    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1
    ikiph = ik

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = 0.5d0*(ra(itrav1 + ipt-1)               &
           + ra(itrav2 + ipt-1) + ra(itrav3 + ipt-1))
    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1
    iepiph = iep

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav7 + ipt-1)
    enddo

    !=======================================================================
    !          3.2.4. Transfert de Rij-epsilon vers v2f
    !=======================================================================

  elseif (iturcp(numcpl).eq.50) then

    !    ATTENTION: CAS NON PRIS EN COMPTE (ARRET DU CALCUL DANS CSCINI.F)

    !=======================================================================
    !          3.2.5. Transfert de Rij-epsilon vers k-omega
    !=======================================================================

  elseif (iturcp(numcpl).eq.60) then

    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1
    ikiph = ik

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = 0.5d0*(ra(itrav1 + ipt-1)               &
           + ra(itrav2 + ipt-1) + ra(itrav3 + ipt-1))
    enddo

    !           Omega
    !           -----
    ipos = ipos + 1
    iomiph = iomg

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav7 + ipt-1)/cmu                  &
           /max(1.0d-10, rvdis(ipt,ipos-1))
    enddo

  endif

  !==============================================================================
  !       3.3 Turbulence dans l'instance locale : modèle v2f (phi-model)
  !==============================================================================

elseif (iturb.eq.50) then

  !=======================================================================
  !          3.3.1. INTERPOLATION EN J'
  !=======================================================================

  !         Préparation des données: interpolation de k en J'

  ikiph  = ik

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ikiph))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iphydp = 0
  iclvar = iclrtp(ikiph,icoef)
  nswrgp = nswrgr(ikiph)
  imligp = imligr(ikiph)
  iwarnp = iwarni(ikiph)
  epsrgp = epsrgr(ikiph)
  climgp = climgr(ikiph)
  extrap = extrag(ikiph)

  call grdcel                                                   &
  !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    ikiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp,  &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,ikiph) , coefa(1,iclvar) , coefb(1,iclvar) ,            &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav1 + ipt-1) = rtp(iel,ikiph)                         &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav1 + ipt-1) = rtp(iel,ikiph)                         &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

  endif

  !         Préparation des données: interpolation de epsilon en J'

  iepiph  = iep

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,iepiph))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iphydp = 0
  iclvar = iclrtp(iepiph,icoef)
  nswrgp = nswrgr(iepiph)
  imligp = imligr(iepiph)
  iwarnp = iwarni(iepiph)
  epsrgp = epsrgr(iepiph)
  climgp = climgr(iepiph)
  extrap = extrag(iepiph)

  call grdcel                                                   &
  !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    iepiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,iepiph) , coefa(1,iclvar) , coefb(1,iclvar) ,           &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav2 + ipt-1) = rtp(iel,iepiph)                        &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav2 + ipt-1) = rtp(iel,iepiph)                        &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

  endif

  !         Préparation des données: interpolation de Phi en J'

  iphiph = iphi

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,iphiph))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iphydp = 0
  iclvar = iclrtp(iphiph,icoef)
  nswrgp = nswrgr(iphiph)
  imligp = imligr(iphiph)
  iwarnp = iwarni(iphiph)
  epsrgp = epsrgr(iphiph)
  climgp = climgr(iphiph)
  extrap = extrag(iphiph)

  call grdcel                                                   &
  !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    iphiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,iphiph) , coefa(1,iclvar) , coefb(1,iclvar) ,           &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


  do ipt = 1, nptdis

    iel = locpts(ipt)

    xjjp = djppts(1,ipt)
    yjjp = djppts(2,ipt)
    zjjp = djppts(3,ipt)

    ra(itrav3 + ipt-1) = rtp(iel,iphiph)                        &
         + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

  enddo

  !         Préparation des données: interpolation de F-barre en J'

  ifbiph = ifb

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ifbiph))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iphydp = 0
  iclvar = iclrtp(ifbiph,icoef)
  nswrgp = nswrgr(ifbiph)
  imligp = imligr(ifbiph)
  iwarnp = iwarni(ifbiph)
  epsrgp = epsrgr(ifbiph)
  climgp = climgr(ifbiph)
  extrap = extrag(ifbiph)

  call grdcel                                                   &
  !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    ifbiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,ifbiph) , coefa(1,iclvar) , coefb(1,iclvar) ,           &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav4 + ipt-1) = rtp(iel,ifbiph)                        &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav4 + ipt-1) = rtp(iel,ifbiph)                        &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

  endif

!=======================================================================
!          3.3.2. Transfert de variable à "iso-modèle"
!=======================================================================

  if (iturcp(numcpl).eq.50) then

    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1
    ikiph  = ik

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav1 + ipt-1)
    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1
    iepiph = iep

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav2 + ipt-1)
    enddo

    !           Phi
    !           ---
    ipos = ipos + 1
    iphiph = iphi

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav3 + ipt-1)
    enddo

    !           F-barre
    !           -------
    ipos = ipos + 1
    ifbiph = ifb

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav4 + ipt-1)
    enddo


    !         ATTENTION: LE COUPLAGE ENTRE UN MODELE V2F ET UN MODELE DE
    !         TURBULENCE DIFFERENT N'EST PAS PRIS EN COMPTE

  elseif (itytu0.eq.2) then
  elseif (itytu0.eq.3) then
  elseif (iturcp(numcpl).eq.60) then
  endif

!==============================================================================
!       3.4 Turbulence dans l'instance locale : modèle omega SST
!==============================================================================

elseif (iturb.eq.60) then

  !=======================================================================
  !          3.4.1. INTERPOLATION EN J'
  !=======================================================================

  !         Préparation des données: interpolation de k en J'

  ikiph  = ik

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ikiph))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iphydp = 0
  iclvar = iclrtp(ikiph,icoef)
  nswrgp = nswrgr(ikiph)
  imligp = imligr(ikiph)
  iwarnp = iwarni(ikiph)
  epsrgp = epsrgr(ikiph)
  climgp = climgr(ikiph)
  extrap = extrag(ikiph)

  call grdcel                                                   &
  !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    ikiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp,  &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,ikiph) , coefa(1,iclvar) , coefb(1,iclvar) ,            &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav1 + ipt-1) = rtp(iel,ikiph)                         &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav1 + ipt-1) = rtp(iel,ikiph)                         &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

  endif

  !         Préparation des données: interpolation de omega en J'

  iomiph  = iomg

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,iomiph))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iphydp = 0
  iclvar = iclrtp(iomiph,icoef)
  nswrgp = nswrgr(iomiph)
  imligp = imligr(iomiph)
  iwarnp = iwarni(iomiph)
  epsrgp = epsrgr(iomiph)
  climgp = climgr(iomiph)
  extrap = extrag(iomiph)

  call grdcel                                                   &
  !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    iomiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,iomiph) , coefa(1,iclvar) , coefb(1,iclvar) ,           &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )

  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav2 + ipt-1) = rtp(iel,iomiph)                        &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      ra(itrav2 + ipt-1) = rtp(iel,iomiph)                        &
           + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

    enddo

  endif

  !=======================================================================
  !          3.4.2. Transfert de variable à "iso-modèle"
  !=======================================================================

  if (iturcp(numcpl).eq.60) then

    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1
    ikiph  = ik

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav1 + ipt-1)
    enddo

    !           Omega
    !           -----
    ipos = ipos + 1
    iomiph = iomg

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav2 + ipt-1)
    enddo

  elseif (itytu0.eq.2) then

    !========================================================================
    !          3.4.3. Transfert de k-omega vers k-epsilon
    !========================================================================
    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1
    ikiph  = ik

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav1 + ipt-1)
    enddo

    !           Omega
    !           -----
    ipos = ipos + 1
    iomiph = iomg

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav2 + ipt-1)*cmu                  &
           *ra(itrav1 + ipt-1)
    enddo

    !========================================================================
    !          3.4.3. Transfert de k-omega vers Rij-epsilon
    !========================================================================

  elseif (itytu0.eq.3) then

    !           Tenseur Rij
    !            ----------
    !           Termes de la diagonal R11,R22,R33

    do isou =1, 3

      ipos = ipos + 1

      do ipt = 1, nptdis
        rvdis(ipt,ipos) = d2s3*ra(itrav1 + ipt-1)
      enddo

    enddo

    !           Termes R12,R13,R23

    iuiph = iu
    iviph = iv
    iwiph = iw

    do isou = 1, 3

      if(isou.eq.1) ivar = iuiph
      if(isou.eq.2) ivar = iviph
      if(isou.eq.3) ivar = iwiph

      inc    = 1
      iccocg = 1
      iphydp = 0
      iclvar = iclrtp(ivar,icoef)
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)

      call grdcel                                               &
      !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,             &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )


      do ipt = 1, nptdis

        iel = locpts(ipt)

        if (isou.eq.1) then
          ra(itrav3 + ipt-1) = w2(iel)
          ra(itrav4 + ipt-1) = w3(iel)
        elseif (isou.eq.2) then
          ra(itrav5 + ipt-1) = w1(iel)
          ra(itrav6 + ipt-1) = w3(iel)
        elseif (isou.eq.3) then
          ra(itrav7 + ipt-1) = w1(iel)
          ra(itrav8 + ipt-1) = w2(iel)
        endif

      enddo

    enddo
    !           Fin de la boucle sur les composantes de la vitesse

    !           R12
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*ra(itrav1 + ipt-1)               &
           /max(1.0d-10, ra(itrav2 + ipt-1))                       &
           *0.5d0*(ra(itrav3 + ipt-1) + ra(itrav5 + ipt-1))
    enddo

    !           R13
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*ra(itrav1 + ipt-1)               &
           /max(1.0d-10, ra(itrav2 + ipt-1))                       &
           *0.5d0*(ra(itrav4 + ipt-1) + ra(itrav7 + ipt-1))
    enddo

    !           R23
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*ra(itrav1 + ipt-1)               &
           /max(1.0d-10, ra(itrav2 + ipt-1))                       &
           *0.5d0*(ra(itrav6 + ipt-1) + ra(itrav8 + ipt-1))
    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1
    iepiph = iep

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = ra(itrav2 + ipt-1)*cmu                  &
           *ra(itrav1 + ipt-1)
    enddo


    !=======================================================================
    !          3.3.4. Transfert de k-omega vers v2f
    !=======================================================================

  elseif (iturcp(numcpl).eq.50) then

    !  ATTENTION: CAS NON PRIS EN COMPTE. ARRET DU CALCUL DANS CSCINI.F

  endif

endif

!=========================================================================
! 4.  PREPARATION DES SCALAIRES
!=========================================================================

if (nscal.gt.0) then

  do iscal = 1, nscal

    ipos = ipos + 1

    ivar = isca(iscal)

! --- Calcul du gradient du scalaire pour interpolation

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,ivar))
      !==========
    endif

    inc    = 1
    iccocg = 1
    iphydp = 0
    iclvar = iclrtp(ivar,icoef)
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    call grdcel                                                   &
    !==========
  ( ifinia , ifinra ,                                             &
    nphas  ,                                                      &
    ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp, &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    ia     ,                                    &
    w4     , w4     , w4     ,                                    &
    rtp(1,ivar)     , coefa(1,iclvar) , coefb(1,iclvar) ,         &
    w1     , w2     , w3     ,                                    &
!         ------   ------   ------
    w4     , w5     , w6     ,                                    &
    ra     )

    ! For a specific face to face coupling, geometric assumptions are made

    if (ifaccp.eq.1) then

      do ipt = 1, nptdis

        iel = locpts(ipt)

! --- Pour les scalaires on veut imposer un dirichlet. On se laisse
!     le choix entre UPWIND, SOLU ou CENTRE. Seul le centré respecte
!     la diffusion si il n'y avait qu'un seul domaine

! -- UPWIND

!        rvdis(ipt,ipos) = rtp(iel,ivar)

! -- SOLU

!        xjf = coopts(1,ipt) - xyzcen(1,iel)
!        yjf = coopts(2,ipt) - xyzcen(2,iel)
!        zjf = coopts(3,ipt) - xyzcen(3,iel)

!        rvdis(ipt,ipos) = rtp(iel,ivar) &
!          + xjf*w1(iel) + yjf*w2(iel) + zjf*w3(iel)

! -- CENTRE

        xjjp = djppts(1,ipt)
        yjjp = djppts(2,ipt)
        zjjp = djppts(3,ipt)

        rvdis(ipt,ipos) = rtp(iel,ivar) &
          + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

      enddo

    ! For a generic coupling, no assumption can be made

    else

      do ipt = 1, nptdis

        iel = locpts(ipt)

        xjjp = djppts(1,ipt)
        yjjp = djppts(2,ipt)
        zjjp = djppts(3,ipt)

        rvdis(ipt,ipos) = rtp(iel,ivar)                             &
          + xjjp*w1(iel) + yjjp*w2(iel) + zjjp*w3(iel)

      enddo

    endif

  enddo

endif


return
end subroutine
