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

!> \file cscpfb.f90
!> \brief Preparation of sending variables for coupling between two instances
!> of Code_Saturne via boundary faces.
!> Received indformation will be transformed into boundary condition
!> in subroutine \ref csc2cl.

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nscal         total number of scalars
!> \param[in]     nptdis
!> \param[in]     numcpl
!> \param[in]     nvcpto
!> \param[in]     locpts
!> \param[in]     rtp           calculated variables at cell centers
!>                              (at current and previous time steps)
!> \param[in]     coefa         boundary conditions for boundary faces
!> \param[in]     coefb         boundary conditions for boundary faces
!> \param[in]     coopts
!> \param[in]     djppts
!> \param[in]     pndpts
!> \param[in,out] rvdis
!> \param[in]     dofpts
!______________________________________________________________________________

subroutine cscpfb &
 ( nscal  ,                                                       &
   nptdis , numcpl , nvcpto,                                      &
   locpts ,                                                       &
   rtp    ,                                                       &
   coefa  , coefb  ,                                              &
   coopts , djppts , pndpts ,                                     &
   rvdis  , dofpts )

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
use field
!===============================================================================

implicit none

! Arguments

integer          nscal
integer          nptdis , numcpl , nvcpto

integer          locpts(nptdis)

double precision rtp(ncelet,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision coopts(3,nptdis), djppts(3,nptdis)
double precision pndpts(nptdis), dofpts(3,nptdis)
double precision rvdis(nptdis,nvcpto)

! Local variables


integer          ipt    , iel    , isou
integer          ivar   , iscal
integer          inc    , iccocg , iclvar, nswrgp
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

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: trav1, trav2, trav3, trav4
double precision, allocatable, dimension(:) :: trav5, trav6, trav7, trav8
double precision, dimension(:), pointer :: crom
!===============================================================================

!=========================================================================
! 1.  INITIALISATIONS
!=========================================================================

! Allocate a temporary array
allocate(grad(ncelet,3))

! Initialize variables to avoid compiler warnings

vitent = 0.d0

! Memoire


! Allocate temporary arrays

allocate(trav1(nptdis))
allocate(trav2(nptdis))
allocate(trav3(nptdis))
allocate(trav4(nptdis))
allocate(trav5(nptdis))
allocate(trav6(nptdis))
allocate(trav7(nptdis))
allocate(trav8(nptdis))

d2s3 = 2.d0/3.d0

call field_get_val_s(icrom, crom)

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

! --- Calcul du gradient de la pression pour interpolation

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtp(1,ipr))
  !==========
endif

inc    = 1
iccocg = 1
iclvar = iclrtp(ipr,icoef)
nswrgp = nswrgr(ipr)
imligp = imligr(ipr)
iwarnp = iwarni(ipr)
epsrgp = epsrgr(ipr)
climgp = climgr(ipr)
extrap = extrag(ipr)

call grdcel &
!==========
  ( ipr , imrgra , inc    , iccocg , nswrgp , imligp ,            &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,ipr) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
    grad   )

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

    rvdis(ipt,ipos) = rtp(iel,ipr) &
         + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

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
           + 0.5d0*crom(iel)*(omgnrl**2 - omgnrd**2)*daxis2

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

    rvdis(ipt,ipos) = (xjpf*grad(iel,1)+yjpf*grad(iel,2)+zjpf*grad(iel,3))  &
         /jpf

  enddo

endif
!       FIn pour la pression


!=========================================================================
! 2.  PREPARATION DE LA VITESSE
!=========================================================================

! --- Calcul du gradient de la vitesse pour interpolation

if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(rtp(1,iu), rtp(1,iv), rtp(1,iw))
  !==========
endif

do isou = 1, 3

  ipos = ipos + 1

  if(isou.eq.1) ivar = iu
  if(isou.eq.2) ivar = iv
  if(isou.eq.3) ivar = iw

  inc    = 1
  iccocg = 1
  iclvar = iclrtp(ivar,icoef)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  call grdcel                                                   &
  !==========
  ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,         &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,             &
    grad   )


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
      !          + xjf*grad(iel,1) + yjf*grad(iel,2) + zjf*grad(iel,3)

      ! -- CENTRE

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      rvdis(ipt,ipos) = rtp(iel,ivar) &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

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
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

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

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ik))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iclvar = iclrtp(ik,icoef)
  nswrgp = nswrgr(ik)
  imligp = imligr(ik)
  iwarnp = iwarni(ik)
  epsrgp = epsrgr(ik)
  climgp = climgr(ik)
  extrap = extrag(ik)

  call grdcel &
  !==========
  ( ik , imrgra , inc    , iccocg , nswrgp , imligp ,             &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,ik) , coefa(1,iclvar) , coefb(1,iclvar) ,               &
    grad   )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav1(ipt) = rtp(iel,ik)                         &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav1(ipt) = rtp(iel,ik)                         &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

  endif

  !         Préparation des données: interpolation de epsilon en J'

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,iep))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iclvar = iclrtp(iep,icoef)
  nswrgp = nswrgr(iep)
  imligp = imligr(iep)
  iwarnp = iwarni(iep)
  epsrgp = epsrgr(iep)
  climgp = climgr(iep)
  extrap = extrag(iep)

  call grdcel &
  !==========
  ( iep , imrgra , inc    , iccocg , nswrgp , imligp ,            &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,iep) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
    grad   )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav2(ipt) = rtp(iel,iep)                        &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav2(ipt) = rtp(iel,iep)                        &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

  endif


!=======================================================================
!          3.1.2.   Transfert de variable à "iso-modèle"
!=======================================================================

  if (itytu0.eq.2) then

    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav1(ipt)
    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav2(ipt)
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
        rvdis(ipt,ipos) = d2s3*trav1(ipt)
      enddo

    enddo

    !           Termes R12,R13,R23

    !           La synchronisation des halos a deja ete faite plus haut

    do isou = 1, 3

      if(isou.eq.1) ivar = iu
      if(isou.eq.2) ivar = iv
      if(isou.eq.3) ivar = iw

      inc    = 1
      iccocg = 1
      iclvar = iclrtp(ivar,icoef)
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)

      call grdcel                                               &
      !==========
  ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,         &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,             &
    grad   )


      do ipt = 1, nptdis

        iel = locpts(ipt)

        if(isou.eq.1) then
          trav3(ipt) = grad(iel,2)
          trav4(ipt) = grad(iel,3)
        elseif(isou.eq.2) then
          trav5(ipt) = grad(iel,1)
          trav6(ipt) = grad(iel,3)
        elseif(isou.eq.3) then
          trav7(ipt) = grad(iel,1)
          trav8(ipt) = grad(iel,2)
        endif

      enddo

    enddo
    !           Fin de la boucle sur les composantes de la vitesse

    !           R12
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*trav1(ipt)**2*cmu / max(1.0d-10, trav2(ipt)) &
           *0.5d0*(trav3(ipt) + trav5(ipt))
    enddo

    !           R13
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*trav1(ipt)**2*cmu / max(1.0d-10, trav2(ipt)) &
           *0.5d0*(trav4(ipt) + trav7(ipt))
    enddo

    !           R23
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*trav1(ipt)**2*cmu / max(1.0d-10,trav2(ipt)) &
           *0.5d0*(trav6(ipt) + trav8(ipt))
    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav2(ipt)
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

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav1(ipt)
    enddo

    !           Omega
    !           -----
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav2(ipt) / cmu / max(1.0d-10, trav1(ipt))
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
    iclvar = iclrtp(ivar,icoef)
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    call grdcel                                                 &
    !==========
  ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,         &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,             &
    grad   )

    ! For a specific face to face coupling, geometric assumptions are made

    if (ifaccp.eq.1) then

      do ipt = 1, nptdis

        iel = locpts(ipt)

        xjjp = djppts(1,ipt)
        yjjp = djppts(2,ipt)
        zjjp = djppts(3,ipt)

        if (isou.eq.1) then
          trav1(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        else if (isou.eq.2) then
          trav2(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        else if (isou.eq.3) then
          trav3(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        else if (isou.eq.4) then
          trav4(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        else if (isou.eq.5) then
          trav5(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        else if (isou.eq.6) then
          trav6(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        endif

      enddo

      ! For a generic coupling, no assumption can be made

    else

      do ipt = 1, nptdis

        iel = locpts(ipt)

        xjjp = djppts(1,ipt)
        yjjp = djppts(2,ipt)
        zjjp = djppts(3,ipt)

        if (isou.eq.1) then
          trav1(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        else if (isou.eq.2) then
          trav2(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        else if (isou.eq.3) then
          trav3(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        else if (isou.eq.4) then
          trav4(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        else if (isou.eq.5) then
          trav5(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        else if (isou.eq.6) then
          trav6(ipt) = rtp(iel,ivar) &
             + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)
        endif

      enddo

    endif

  enddo

  !         Préparation des données: interpolation de epsilon en J'

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,iep))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iclvar = iclrtp(iep,icoef)
  nswrgp = nswrgr(iep)
  imligp = imligr(iep)
  iwarnp = iwarni(iep)
  epsrgp = epsrgr(iep)
  climgp = climgr(iep)
  extrap = extrag(iep)

  call grdcel &
  !==========
  ( iep , imrgra , inc    , iccocg , nswrgp , imligp ,            &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,iep) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
    grad   )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav7(ipt) = rtp(iel,iep)                        &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav7(ipt) = rtp(iel,iep)                        &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

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

      if (isou.eq.1) then
        do ipt = 1, nptdis
          rvdis(ipt,ipos) = trav1(ipt)
        enddo
      else if (isou.eq.2) then
        do ipt = 1, nptdis
          rvdis(ipt,ipos) = trav2(ipt)
        enddo
      else if (isou.eq.3) then
        do ipt = 1, nptdis
          rvdis(ipt,ipos) = trav3(ipt)
        enddo
      else if (isou.eq.4) then
        do ipt = 1, nptdis
          rvdis(ipt,ipos) = trav4(ipt)
        enddo
      else if (isou.eq.5) then
        do ipt = 1, nptdis
          rvdis(ipt,ipos) = trav5(ipt)
        enddo
      else if (isou.eq.6) then
        do ipt = 1, nptdis
          rvdis(ipt,ipos) = trav6(ipt)
        enddo
      endif

    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav7(ipt)
    enddo

    !=======================================================================
    !          3.2.3. Transfert de Rij-epsilon vers k-epsilon
    !=======================================================================

  elseif (itytu0.eq.2) then

    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = 0.5d0*(trav1(ipt) + trav2(ipt) + trav3(ipt))
    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav7(ipt)
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

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = 0.5d0*(trav1(ipt) + trav2(ipt) + trav3(ipt))
    enddo

    !           Omega
    !           -----
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav7(ipt) / cmu / max(1.0d-10, rvdis(ipt,ipos-1))
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

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ik))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iclvar = iclrtp(ik,icoef)
  nswrgp = nswrgr(ik)
  imligp = imligr(ik)
  iwarnp = iwarni(ik)
  epsrgp = epsrgr(ik)
  climgp = climgr(ik)
  extrap = extrag(ik)

  call grdcel &
  !==========
  ( ik , imrgra , inc    , iccocg , nswrgp , imligp ,             &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,ik) , coefa(1,iclvar) , coefb(1,iclvar) ,               &
    grad   )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav1(ipt) = rtp(iel,ik)                         &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav1(ipt) = rtp(iel,ik)                         &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

  endif

  !         Préparation des données: interpolation de epsilon en J'

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,iep))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iclvar = iclrtp(iep,icoef)
  nswrgp = nswrgr(iep)
  imligp = imligr(iep)
  iwarnp = iwarni(iep)
  epsrgp = epsrgr(iep)
  climgp = climgr(iep)
  extrap = extrag(iep)

  call grdcel &
  !==========
  ( iep , imrgra , inc    , iccocg , nswrgp , imligp ,            &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,iep) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
    grad   )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav2(ipt) = rtp(iel,iep)                        &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav2(ipt) = rtp(iel,iep)                        &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

  endif

  !         Préparation des données: interpolation de Phi en J'

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,iphi))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iclvar = iclrtp(iphi,icoef)
  nswrgp = nswrgr(iphi)
  imligp = imligr(iphi)
  iwarnp = iwarni(iphi)
  epsrgp = epsrgr(iphi)
  climgp = climgr(iphi)
  extrap = extrag(iphi)

  call grdcel &
  !==========
  ( iphi , imrgra , inc    , iccocg , nswrgp , imligp ,           &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,iphi) , coefa(1,iclvar) , coefb(1,iclvar) ,             &
    grad   )


  do ipt = 1, nptdis

    iel = locpts(ipt)

    xjjp = djppts(1,ipt)
    yjjp = djppts(2,ipt)
    zjjp = djppts(3,ipt)

    trav3(ipt) = rtp(iel,iphi)                        &
         + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

  enddo

  !         Préparation des données: interpolation de F-barre en J'

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ifb))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iclvar = iclrtp(ifb,icoef)
  nswrgp = nswrgr(ifb)
  imligp = imligr(ifb)
  iwarnp = iwarni(ifb)
  epsrgp = epsrgr(ifb)
  climgp = climgr(ifb)
  extrap = extrag(ifb)

  call grdcel &
  !==========
  ( ifb , imrgra , inc    , iccocg , nswrgp , imligp ,            &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,ifb) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
    grad   )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav4(ipt) = rtp(iel,ifb)                        &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav4(ipt) = rtp(iel,ifb)                        &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

  endif

!=======================================================================
!          3.3.2. Transfert de variable à "iso-modèle"
!=======================================================================

  if (iturcp(numcpl).eq.50) then

    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav1(ipt)
    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav2(ipt)
    enddo

    !           Phi
    !           ---
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav3(ipt)
    enddo

    !           F-barre
    !           -------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav4(ipt)
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

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ik))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iclvar = iclrtp(ik,icoef)
  nswrgp = nswrgr(ik)
  imligp = imligr(ik)
  iwarnp = iwarni(ik)
  epsrgp = epsrgr(ik)
  climgp = climgr(ik)
  extrap = extrag(ik)

  call grdcel &
  !==========
  ( ik , imrgra , inc    , iccocg , nswrgp , imligp ,             &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,ik) , coefa(1,iclvar) , coefb(1,iclvar) ,               &
    grad   )


  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav1(ipt) = rtp(iel,ik)                         &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav1(ipt) = rtp(iel,ik)                         &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

  endif

  !         Préparation des données: interpolation de omega en J'

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,iomg))
    !==========
  endif

  inc    = 1
  iccocg = 1
  iclvar = iclrtp(iomg,icoef)
  nswrgp = nswrgr(iomg)
  imligp = imligr(iomg)
  iwarnp = iwarni(iomg)
  epsrgp = epsrgr(iomg)
  climgp = climgr(iomg)
  extrap = extrag(iomg)

  call grdcel &
  !==========
  ( iomg , imrgra , inc    , iccocg , nswrgp , imligp ,           &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,iomg) , coefa(1,iclvar) , coefb(1,iclvar) ,             &
    grad   )

  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav2(ipt) = rtp(iel,iomg)                        &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

    ! For a generic coupling, no assumption can be made

  else

    do ipt = 1, nptdis

      iel = locpts(ipt)

      xjjp = djppts(1,ipt)
      yjjp = djppts(2,ipt)
      zjjp = djppts(3,ipt)

      trav2(ipt) = rtp(iel,iomg)                        &
           + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

    enddo

  endif

  !=======================================================================
  !          3.4.2. Transfert de variable à "iso-modèle"
  !=======================================================================

  if (iturcp(numcpl).eq.60) then

    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav1(ipt)
    enddo

    !           Omega
    !           -----
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav2(ipt)
    enddo

  elseif (itytu0.eq.2) then

    !========================================================================
    !          3.4.3. Transfert de k-omega vers k-epsilon
    !========================================================================
    !           Energie turbulente
    !           ------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav1(ipt)
    enddo

    !           Omega
    !           -----
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav2(ipt)*cmu*trav1(ipt)
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
        rvdis(ipt,ipos) = d2s3*trav1(ipt)
      enddo

    enddo

    !           Termes R12,R13,R23

    do isou = 1, 3

      if(isou.eq.1) ivar = iu
      if(isou.eq.2) ivar = iv
      if(isou.eq.3) ivar = iw

      inc    = 1
      iccocg = 1
      iclvar = iclrtp(ivar,icoef)
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)

      call grdcel                                               &
      !==========
  ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,         &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,             &
    grad   )


      do ipt = 1, nptdis

        iel = locpts(ipt)

        if (isou.eq.1) then
          trav3(ipt) = grad(iel,2)
          trav4(ipt) = grad(iel,3)
        elseif (isou.eq.2) then
          trav5(ipt) = grad(iel,1)
          trav6(ipt) = grad(iel,3)
        elseif (isou.eq.3) then
          trav7(ipt) = grad(iel,1)
          trav8(ipt) = grad(iel,2)
        endif

      enddo

    enddo
    !           Fin de la boucle sur les composantes de la vitesse

    !           R12
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*trav1(ipt) / max(1.0d-10, trav2(ipt))  &
           *0.5d0*(trav3(ipt) + trav5(ipt))
    enddo

    !           R13
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*trav1(ipt) / max(1.0d-10, trav2(ipt))  &
           *0.5d0*(trav4(ipt) + trav7(ipt))
    enddo

    !           R23
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = -2.0d0*trav1(ipt) / max(1.0d-10, trav2(ipt))  &
           *0.5d0*(trav6(ipt) + trav8(ipt))
    enddo

    !           Dissipation turbulente
    !           ----------------------
    ipos = ipos + 1

    do ipt = 1, nptdis
      rvdis(ipt,ipos) = trav2(ipt)*cmu*trav1(ipt)
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
    iclvar = iclrtp(ivar,icoef)
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    call grdcel &
    !==========
  ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,         &
    iwarnp , nfecra ,                                             &
    epsrgp , climgp , extrap ,                                    &
    rtp(1,ivar)     , coefa(1,iclvar) , coefb(1,iclvar) ,         &
    grad   )

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
!          + xjf*grad(iel,1) + yjf*grad(iel,2) + zjf*grad(iel,3)

! -- CENTRE

        xjjp = djppts(1,ipt)
        yjjp = djppts(2,ipt)
        zjjp = djppts(3,ipt)

        rvdis(ipt,ipos) = rtp(iel,ivar) &
          + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

      enddo

    ! For a generic coupling, no assumption can be made

    else

      do ipt = 1, nptdis

        iel = locpts(ipt)

        xjjp = djppts(1,ipt)
        yjjp = djppts(2,ipt)
        zjjp = djppts(3,ipt)

        rvdis(ipt,ipos) = rtp(iel,ivar)                             &
          + xjjp*grad(iel,1) + yjjp*grad(iel,2) + zjjp*grad(iel,3)

      enddo

    endif

  enddo

endif

! Free memory
deallocate(grad)
deallocate(trav1, trav2, trav3, trav4)
deallocate(trav5, trav6, trav7, trav8)

return
end subroutine
