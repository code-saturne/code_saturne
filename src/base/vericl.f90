!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!===============================================================================
! Function :
! --------

!> \file vericl.f90
!>
!> \brief Check boundary condition code.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in,out] itypfb        face boundary condition type
!> \param[in,out] icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!_______________________________________________________________________________

subroutine vericl &
 ( nvar   , nscal  ,                                              &
   itypfb , icodcl ,                                              &
   rcodcl )

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
use albase
use ppppar
use ppthch
use ppincl
use parall
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          itypfb(nfabor)
integer          icodcl(nfabor,nvar)

double precision rcodcl(nfabor,nvar,3)

! Local variables

character(len=80) :: chaine

double precision grav2
integer          ifac, ivar, icode, f_dim
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
integer          icodni(2), icodvi(3), icodpp(2), icodtb(8), icodsc(2)
integer          icodvf(2), icoduv(3), icodct(11), icodus(4)

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

do ipp = 1, 2
  icodni(ipp) = -1
  icodpp(ipp) = -1
  icodsc(ipp) = -1
  icodvf(ipp) = -1
enddo
do ipp = 1, 3
  icodvi(ipp) = -1
  icoduv(ipp) = -1
enddo
do ipp = 1, 8
  icodtb(ipp) = -1
enddo
do ipp = 1, 11
  icodct(ipp) = -1
enddo
do ipp = 1, 4
  icodus(ipp) = -1
enddo


!===============================================================================
! 2.  VERIFICATION DE LA CONSISTANCE DES CL
!===============================================================================

! 2.1  INITIALISATION
! ====================

! Dans cs_user_boundary_conditions, on se donne une grande liberte
!  pour la specif. des c.l sur les variables.
!  Neanmoins pour limiter la plage des tests, on se
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
if (iokcod.ne.0) then
  do ivar = 1, nvar
    do ifac = 1, nfabor
      icode = icodcl(ifac,ivar)
      if(icode.eq. 0) then
        if (itypfb(ifac).gt.0) then
          itypfb(ifac) = -itypfb(ifac)
        endif
        icodni(1) = ivar
        icodni(2) = ifac
        nstoni = nstoni + 1
      endif
    enddo
  enddo
endif


! 2.3 VERIFICATIONS DE L'ADMISSIBILITE DES CONDITIONS
! ====================================================

! --- Conditions admissibles pour les composantes de vitesse
do ifac = 1, nfabor

  icodcu = icodcl(ifac,iu)
  icodcv = icodcl(ifac,iv)
  icodcw = icodcl(ifac,iw)

  if(icodcu.ne. 1.and.icodcu.ne. 2.and.icodcu.ne. 3.and.        &
     icodcu.ne. 4.and.icodcu.ne. 5.and.icodcu.ne. 6.and.        &
     icodcu.ne.11.and.                                          &
     icodcu.ne. 9.and.icodcu.ne.13.and.icodcu.ne.14) then
    if (itypfb(ifac).gt.0) then
      itypfb(ifac) = -itypfb(ifac)
    endif
    icodvi(1) = iu
    icodvi(2) = icodcl(ifac,iu)
    icodvi(3) = -1
    nstvit = nstvit + 1
  endif
  if(icodcv.ne. 1.and.icodcv.ne. 2.and.icodcv.ne. 3.and.        &
     icodcv.ne. 4.and.icodcv.ne. 5.and.icodcv.ne. 6.and.        &
     icodcv.ne.11.and.                                          &
     icodcv.ne. 9.and.icodcv.ne.13.and.icodcv.ne.14) then
    if (itypfb(ifac).gt.0) then
      itypfb(ifac) = -itypfb(ifac)
    endif
    icodvi(1) = iv
    icodvi(2) = icodcl(ifac,iv)
    icodvi(3) = -1
    nstvit = nstvit + 1
  endif
  if(icodcw.ne. 1.and.icodcw.ne. 2.and.icodcw.ne. 3.and.        &
     icodcw.ne. 4.and.icodcw.ne. 5.and.icodcw.ne. 6.and.        &
     icodcw.ne.11.and.                                          &
     icodcw.ne. 9.and.icodcw.ne.13.and.icodcw.ne.14) then
    if (itypfb(ifac).gt.0) then
      itypfb(ifac) = -itypfb(ifac)
    endif
    icodvi(1) = iw
    icodvi(2) = icodcl(ifac,iw)
    icodvi(3) = -1
    nstvit = nstvit + 1
  endif

  ! --- verification que la rugosite est initialisee si icodl=6
  if(icodcu.eq.6 .and. rcodcl(ifac,iu,3).lt.epzero)then
    if (itypfb(ifac).gt.0) then
      itypfb(ifac) = -itypfb(ifac)
    endif
    icodvi(1) = -6
    icodvi(2) = icodcl(ifac,iw)
    icodvi(3) = -1
    nstvit = nstvit + 1
  endif

  ! --- on interdit les parois rugueuses en compressible
  if (icodcu.eq.6 .and. ippmod(icompf).gt.0) then
    icodvi(1) = iu
    icodvi(2) = icodcl(ifac,iu)
    icodvi(3) = 1
    nstvit = nstvit + 1
  endif

enddo

! --- Conditions admissibles pour la pression
do ifac = 1, nfabor

  if (icodcl(ifac,ipr).ne.1 .and. icodcl(ifac,ipr).ne.2 .and.     &
      icodcl(ifac,ipr).ne.3 .and. icodcl(ifac,ipr).ne.11.and.     &
      icodcl(ifac,ipr).ne.13) then
    if (itypfb(ifac).gt.0) then
      itypfb(ifac) = -itypfb(ifac)
    endif
    icodpp(1) = ipr
    icodpp(2) = icodcl(ifac,ipr)
    nstopp = nstopp + 1
  endif

enddo

! --- Conditions admissibles pour k et epsilon
if (itytur.eq.2) then

  do ifac = 1, nfabor

    if((icodcl(ifac,ik ).ne. 1.and.                          &
        icodcl(ifac,ik ).ne. 2.and.                          &
        icodcl(ifac,ik ).ne. 3.and.                          &
        icodcl(ifac,ik ).ne. 5.and.                          &
        icodcl(ifac,ik ).ne.13.and.                          &
        icodcl(ifac,ik ).ne. 6     ).or.                     &
       (icodcl(ifac,iep).ne. 1.and.                          &
        icodcl(ifac,iep).ne. 2.and.                          &
        icodcl(ifac,iep).ne. 3.and.                          &
        icodcl(ifac,iep).ne. 5.and.                          &
        icodcl(ifac,iep).ne.13.and.                          &
        icodcl(ifac,iep).ne. 6     ) )then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = ik
      icodtb(2) = icodcl(ifac,ik)
      icodtb(3) = iep
      icodtb(4) = icodcl(ifac,iep)
      nstoke = nstoke + 1
    endif

  enddo

  ! --- Conditions admissibles pour Rij et epsilon
elseif(itytur.eq.3) then

  ivar = ir11
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if (icode.ne. 1.and.icode.ne. 2.and. icode.ne. 3.and.         &
        icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6 .and.icode.ne.13) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = ivar
      icodtb(2) = icode
      nstrij = nstrij + 1
    endif
  enddo

  ivar = ir22
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if (icode.ne. 1.and.icode.ne. 2.and. icode.ne. 3.and.         &
        icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6 .and.icode.ne.13) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = ivar
      icodtb(2) = icode
      nstrij = nstrij + 1
    endif
  enddo

  ivar = ir33
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if (icode.ne. 1.and.icode.ne. 2.and. icode.ne. 3.and.         &
        icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6 .and.icode.ne.13) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = ivar
      icodtb(2) = icode
      nstrij = nstrij + 1
    endif
  enddo

  ivar = ir12
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if (icode.ne. 1.and.icode.ne. 2.and. icode.ne. 3.and.         &
        icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6 .and.icode.ne.13) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = ivar
      icodtb(2) = icode
      nstrij = nstrij + 1
    endif
  enddo

  ivar = ir23
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if (icode.ne. 1.and.icode.ne. 2.and. icode.ne. 3.and.         &
        icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6 .and.icode.ne.13) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = ivar
      icodtb(2) = icode
      nstrij = nstrij + 1
    endif
  enddo

  ivar = ir13
  do ifac = 1, nfabor
    icode = icodcl(ifac,ivar)
    if (icode.ne. 1.and.icode.ne. 2.and. icode.ne. 3.and.         &
        icode.ne. 4.and.icode.ne. 5.and.icode.ne. 6 .and.icode.ne.13) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = ivar
      icodtb(2) = icode
      nstrij = nstrij + 1
    endif
  enddo

  do ifac = 1, nfabor
    icode = icodcl(ifac,iep)
    if (icode.ne. 1.and.icode.ne. 2.and. icode.ne. 3.and.         &
        icode.ne. 5.and.icode.ne. 6 .and.icode.ne.13) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = iep
      icodtb(2) = icode
      icodtb(3) = -1
      nstrij = nstrij + 1
    endif
  enddo

  if (iturb.eq.32) then
    do ifac = 1, nfabor
      icode = icodcl(ifac,ial)
      if (icode.ne. 1.and.icode.ne. 2.and. icode.ne. 3.and.         &
          icode.ne. 5.and.icode.ne. 6 .and.icode.ne.13) then
        if (itypfb(ifac).gt.0) then
          itypfb(ifac) = -itypfb(ifac)
        endif
        icodtb(1) = ial
        icodtb(2) = icode
        icodtb(3) = -1
        nstrij = nstrij + 1
      endif
      ! No rough wall with EBRSM
      do ivar = 1, nvar
        if (icodcl(ifac,ivar).eq.6) then
          if (itypfb(ifac).gt.0) then
            itypfb(ifac) = -itypfb(ifac)
          endif
          icodtb(1) = ial
          icodtb(2) = icode
          icodtb(3) = ivar
          nstrij = nstrij + 1
        endif
      enddo
    enddo
  endif

  ! --- Conditions admissibles pour k, epsilon, phi et f_barre
elseif (iturb.eq.50) then

  do ifac = 1, nfabor

    if((icodcl(ifac,ik ).ne. 1.and.                          &
        icodcl(ifac,ik ).ne. 2.and.                          &
        icodcl(ifac,ik ).ne. 3.and.                          &
        icodcl(ifac,ik ).ne. 5.and.                          &
        icodcl(ifac,ik ).ne.13.and.                          &
        icodcl(ifac,ik ).ne. 6     ).or.                     &
       (icodcl(ifac,iep).ne. 1.and.                          &
        icodcl(ifac,iep).ne. 2.and.                          &
        icodcl(ifac,iep).ne. 3.and.                          &
        icodcl(ifac,iep).ne. 5.and.                          &
        icodcl(ifac,iep).ne. 6     ).or.                     &
       (icodcl(ifac,iphi).ne. 1.and.                         &
        icodcl(ifac,iphi).ne. 2.and.                         &
        icodcl(ifac,iphi).ne. 3.and.                         &
        icodcl(ifac,iphi).ne. 5.and.                         &
        icodcl(ifac,iphi).ne.13.and.                         &
        icodcl(ifac,iphi).ne. 6     ).or.                    &
       (icodcl(ifac,ifb).ne. 1.and.                          &
        icodcl(ifac,ifb).ne. 2.and.                          &
        icodcl(ifac,ifb).ne. 3.and.                          &
        icodcl(ifac,ifb).ne. 5.and.                          &
        icodcl(ifac,ifb).ne.13.and.                          &
        icodcl(ifac,ifb).ne. 6     ) )then

      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = ik
      icodtb(2) = icodcl(ifac,ik)
      icodtb(3) = iep
      icodtb(4) = icodcl(ifac,iep)
      icodtb(5) = iphi
      icodtb(6) = icodcl(ifac,iphi)
      icodtb(7) = ifb
      icodtb(8) = icodcl(ifac,ifb)
      nstov2 = nstov2 + 1

    endif

  enddo

! --- Conditions admissibles pour k, epsilon, phi et alpha
elseif (iturb.eq.51) then

  do ifac = 1, nfabor

    if((icodcl(ifac,ik ).ne. 1.and.                          &
        icodcl(ifac,ik ).ne. 2.and.                          &
        icodcl(ifac,ik ).ne. 3.and.                          &
        icodcl(ifac,ik ).ne. 5.and.                          &
        icodcl(ifac,ik ).ne.13.and.                          &
        icodcl(ifac,ik ).ne. 6     ).or.                     &
       (icodcl(ifac,iep).ne. 1.and.                          &
        icodcl(ifac,iep).ne. 2.and.                          &
        icodcl(ifac,iep).ne. 3.and.                          &
        icodcl(ifac,iep).ne. 5.and.                          &
        icodcl(ifac,iep).ne.13.and.                          &
        icodcl(ifac,iep).ne. 6     ).or.                     &
       (icodcl(ifac,iphi).ne. 1.and.                         &
        icodcl(ifac,iphi).ne. 2.and.                         &
        icodcl(ifac,iphi).ne. 3.and.                         &
        icodcl(ifac,iphi).ne. 5.and.                         &
        icodcl(ifac,iphi).ne.13.and.                         &
        icodcl(ifac,iphi).ne. 6     ).or.                    &
       (icodcl(ifac,ial).ne. 1.and.                          &
        icodcl(ifac,ial).ne. 2.and.                          &
        icodcl(ifac,ial).ne. 3.and.                          &
        icodcl(ifac,ial).ne. 5.and.                          &
        icodcl(ifac,ial).ne.13.and.                          &
        icodcl(ifac,ial).ne. 6     ) )then

      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = ik
      icodtb(2) = icodcl(ifac,ik)
      icodtb(3) = iep
      icodtb(4) = icodcl(ifac,iep)
      icodtb(5) = iphi
      icodtb(6) = icodcl(ifac,iphi)
      icodtb(7) = ial
      icodtb(8) = icodcl(ifac,ial)
      nstov2 = nstov2 + 1

    endif

  enddo

! --- Conditions admissibles pour k et omega
elseif (iturb.eq.60) then

  do ifac = 1, nfabor

    if((icodcl(ifac,ik ).ne. 1.and.                          &
        icodcl(ifac,ik ).ne. 2.and.                          &
        icodcl(ifac,ik ).ne. 3.and.                          &
        icodcl(ifac,ik ).ne. 5.and.                          &
        icodcl(ifac,ik ).ne.13.and.                          &
        icodcl(ifac,ik ).ne. 6     ).or.                     &
       (icodcl(ifac,iomg).ne. 1.and.                         &
        icodcl(ifac,iomg).ne. 2.and.                         &
        icodcl(ifac,iomg).ne. 3.and.                         &
        icodcl(ifac,iomg).ne. 5.and.                         &
        icodcl(ifac,iomg).ne.13.and.                         &
        icodcl(ifac,iomg).ne. 6     ) )then

      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = ik
      icodtb(2) = icodcl(ifac,ik)
      icodtb(3) = iomg
      icodtb(4) = icodcl(ifac,iomg)
      nstokw = nstokw + 1

    endif

  enddo

  ! --- Conditions admissibles pour Spalart-Allmaras
elseif (iturb.eq.70) then

  do ifac = 1, nfabor

    if(icodcl(ifac,inusa ).ne. 1.and.                          &
       icodcl(ifac,inusa ).ne. 2.and.                          &
       icodcl(ifac,inusa ).ne. 3.and.                          &
       icodcl(ifac,inusa ).ne. 5.and.                          &
       icodcl(ifac,inusa ).ne.13.and.                          &
       icodcl(ifac,inusa ).ne. 6           )then

      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodtb(1) = inusa
      icodtb(2) = icodcl(ifac,inusa)
      nstonu = nstonu + 1

    endif

  enddo

endif

! --- Conditions admissibles pour les scalaires
if(nscal.ge.1) then
  do iis = 1,nscal
    ivar = isca(iis)
    call field_get_dim(ivarfl(ivar), f_dim)
    if (f_dim.eq.1) then
      do ifac = 1, nfabor
        if(icodcl(ifac,ivar).ne. 1.and.                             &
           icodcl(ifac,ivar).ne. 2.and.                             &
           icodcl(ifac,ivar).ne. 3.and.                             &
           icodcl(ifac,ivar).ne. 5.and.                             &
           icodcl(ifac,ivar).ne.12.and.                             &
           icodcl(ifac,ivar).ne.13.and.                             &
           icodcl(ifac,ivar).ne. 6 ) then
          if (itypfb(ifac).gt.0) then
            itypfb(ifac) = -itypfb(ifac)
          endif
          icodsc(1) = ivar
          icodsc(2) = icodcl(ifac,ivar)
          nstosc = nstosc + 1
        endif
        if(icodcl(ifac,ivar).eq. 5.and.                             &
           iscavr(iis).gt.0        ) then
          if (itypfb(ifac).gt.0) then
            itypfb(ifac) = -itypfb(ifac)
          endif
          icodvf(1) = ivar
          icodvf(2) = icodcl(ifac,ivar)
          nstovf = nstovf + 1
        endif
        if(icodcl(ifac,ivar).eq. 6.and.                             &
           iscavr(iis).gt.0        ) then
          if (itypfb(ifac).gt.0) then
            itypfb(ifac) = -itypfb(ifac)
          endif
          icodvf(1) = ivar
          icodvf(2) = icodcl(ifac,ivar)
          nstovf = nstovf + 1
        endif
  ! --- verification que la rugosite scalaire est initialisee si icodl=6
        if(icodcl(ifac,ivar).eq.6.and.                              &
           rcodcl(ifac,iv,3).lt.epzero)then
          if (itypfb(ifac).gt.0) then
            itypfb(ifac) = -itypfb(ifac)
          endif
          icodsc(1) = -6
          icodsc(2) = icodcl(ifac,ivar)
          nstosc = nstosc + 1
        endif
      enddo
    else ! additional vector variables
      do ifac = 1, nfabor
        if(icodcl(ifac,ivar).ne. 1.and.                             &
           icodcl(ifac,ivar).ne. 2.and.                             &
           icodcl(ifac,ivar).ne. 3.and.                             &
           icodcl(ifac,ivar).ne. 4.and.                             &
           icodcl(ifac,ivar).ne. 5.and.                             &
           icodcl(ifac,ivar).ne.11.and.                             &
           icodcl(ifac,ivar).ne.13.and.                             &
           icodcl(ifac,ivar).ne.14.and.                             &
           icodcl(ifac,ivar).ne. 6 ) then
          if (itypfb(ifac).gt.0) then
            itypfb(ifac) = -itypfb(ifac)
          endif
          icodsc(1) = ivar
          icodsc(2) = icodcl(ifac,ivar)
          nstosc = nstosc + 1
        endif
  ! --- verification que la rugosite scalaire est initialisee si icodl=6
        if(icodcl(ifac,ivar).eq.6.and.                              &
           rcodcl(ifac,iv,3).lt.epzero)then
          if (itypfb(ifac).gt.0) then
            itypfb(ifac) = -itypfb(ifac)
          endif
          icodsc(1) = -6
          icodsc(2) = icodcl(ifac,ivar)
          nstosc = nstosc + 1
        endif
      enddo
    endif ! dim
  enddo
endif

! --- Check if the gravity is non zero in case of free-surface
iok = 0
if (iale.ge.1) then
  grav2 = gx**2 + gy**2 + gz**2
  do ifac = 1, nfabor
    if (ialtyb(ifac).eq.ifresf.and.grav2.le.epzero**2) then
      iok = 1
    endif
  enddo
endif

if (irangp.ge.0) call parcmx(iok)

if (iok.ne.0) then
  write(nfecra,2001)
  call csexit (1)
  !==========
endif

! 2.4 VERIFICATIONS DES COHERENCES INTER VARIABLES
! ================================================

! --- Coherence pour les composantes de vitesse
do ifac = 1, nfabor

  icodcu = icodcl(ifac,iu)
  icodcv = icodcl(ifac,iv)
  icodcw = icodcl(ifac,iw)

  if(icodcu.eq.4.or.icodcu.eq.5.or.icodcu.eq.6.or.                &
       icodcu.eq.9.or.                                            &
       icodcv.eq.4.or.icodcv.eq.5.or.icodcv.eq.6.or.              &
       icodcv.eq.9.or.                                            &
       icodcw.eq.4.or.icodcw.eq.5.or.icodcw.eq.6.or.              &
       icodcw.eq.9                )then

    if (icodcu.ne.icodcv .or. icodcu.ne.icodcw .or.               &
        icodcv.ne.icodcw ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icoduv(1) = icodcu
      icoduv(2) = icodcv
      icoduv(3) = icodcw
      nstuvw = nstuvw + 1
    endif
  endif

  ! --- Coherence vitesse pression

  !      Remarques :
  !        Pas de regle stricte de coherence vitesse/pression.
  !        Avant on imposait un Dirichlet sur la pression pour en
  !        entree/sortie, mais cela ne semble pas imperatif. Le test
  !        est laisse en commentaire pour etre recupere si necessaire.

  !        if (icodcu.eq.9 .or. icodcv.eq.9 .or. icodcw.eq.9) THEN
  !          if (icodcl(ifac,ipriph).ne.1) then
  !            nstoup = nstoup + 1
  !          endif
  !        endif

enddo

! --- Coherence vitesse turbulence

if(itytur.eq.2) then

  do ifac = 1, nfabor

    icodcu = icodcl(ifac,iu)
    icodcv = icodcl(ifac,iv)
    icodcw = icodcl(ifac,iw)
    icodck = icodcl(ifac,ik)
    icodce = icodcl(ifac,iep)

    if( (icodcu.eq.5 .or. icodcv.eq.5 .or. icodcw.eq.5 .or.      &
         icodck.eq.5 .or. icodce.eq.5) .and.                     &
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.     &
         icodck.ne.5 .or. icodce.ne.5)                    ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = ik
      icodct(2) = icodcl(ifac,ik)
      icodct(3) = iep
      icodct(4) = icodcl(ifac,iep)
      icodct(5) = icodcu
      icodct(6) = icodcv
      icodct(7) = icodcw
      nstuke = nstuke + 1
    endif

    if( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.      &
         icodck.eq.6 .or. icodce.eq.6) .and.                     &
         (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.     &
         icodck.ne.6 .or. icodce.ne.6)                    ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = ik
      icodct(2) = icodcl(ifac,ik)
      icodct(3) = iep
      icodct(4) = icodcl(ifac,iep)
      icodct(5) = icodcu
      icodct(6) = icodcv
      icodct(7) = icodcw
      nstuke = nstuke + 1
    endif

  enddo

elseif(iturb.eq.30.or.iturb.eq.31) then

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

    if(  (icodcu.eq.5 .or. icodcv.eq.5 .or. icodcw.eq.5 .or.       &
         icor11.eq.5 .or. icor22.eq.5 .or.                         &
         icor33.eq.5 .or. icor12.eq.5 .or.                         &
         icor13.eq.5 .or. icor23.eq.5 .or.                         &
         icodce.eq.5                      ) .and.                  &
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.       &
         icor11.ne.5 .or. icor22.ne.5 .or.                         &
         icor33.ne.5 .or. icor12.ne.5 .or.                         &
         icor13.ne.5 .or. icor23.ne.5 .or.                         &
         icodce.ne.5                      )      ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = icor11
      icodct(2) = icor22
      icodct(3) = icor33
      icodct(4) = icor12
      icodct(5) = icor13
      icodct(6) = icor23
      icodct(7) = icodce
      icodct(8) = icodcu
      icodct(9) = icodcv
      icodct(10) = icodcw
      nsurij = nsurij + 1
    endif

    if(  (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.       &
         icor11.eq.6 .or. icor22.eq.6 .or.                         &
         icor33.eq.6 .or. icor12.eq.6 .or.                         &
         icor13.eq.6 .or. icor23.eq.6 .or.                         &
         icodce.eq.6                      ) .and.                  &
         (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.       &
         icor11.ne.6 .or. icor22.ne.6 .or.                         &
         icor33.ne.6 .or. icor12.ne.6 .or.                         &
         icor13.ne.6 .or. icor23.ne.6 .or.                         &
         icodce.ne.6                      )      ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = icor11
      icodct(2) = icor22
      icodct(3) = icor33
      icodct(4) = icor12
      icodct(5) = icor13
      icodct(6) = icor23
      icodct(7) = icodce
      icodct(8) = icodcu
      icodct(9) = icodcv
      icodct(10) = icodcw
      nsurij = nsurij + 1
    endif

    if(  (icodcu.eq.4 .or. icodcv.eq.4 .or. icodcw.eq.4 .or.       &
         icor11.eq.4 .or. icor22.eq.4 .or.                         &
         icor33.eq.4 .or. icor12.eq.4 .or.                         &
         icor13.eq.4 .or. icor23.eq.4                              &
         ) .and.                         &
         (icodcu.ne.4 .or. icodcv.ne.4 .or. icodcw.ne.4 .or.       &
         icor11.ne.4 .or. icor22.ne.4 .or.                         &
         icor33.ne.4 .or. icor12.ne.4 .or.                         &
         icor13.ne.4 .or. icor23.ne.4 .or.                         &
         icodce.ne.3) ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = icor11
      icodct(2) = icor22
      icodct(3) = icor33
      icodct(4) = icor12
      icodct(5) = icor13
      icodct(6) = icor23
      icodct(7) = icodce
      icodct(8) = icodcu
      icodct(9) = icodcv
      icodct(10) = icodcw
      nsurij = nsurij + 1
    endif

  enddo

elseif (iturb.eq.32) then
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
    icodca = icodcl(ifac,ial)

    if ( (icodcu.eq.5 .or. icodcv.eq.5 .or. icodcw.eq.5 .or.    &
          icor11.eq.5 .or. icor22.eq.5 .or.                     &
          icor33.eq.5 .or. icor12.eq.5 .or.                     &
          icor13.eq.5 .or. icor23.eq.5 .or.                     &
          icodce.eq.5 .or. icodca.eq.5            ) .and.       &
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.    &
          icor11.ne.5 .or. icor22.ne.5 .or.                     &
          icor33.ne.5 .or. icor12.ne.5 .or.                     &
          icor13.ne.5 .or. icor23.ne.5 .or.                     &
          icodce.ne.5 .or. icodca.ne.5     )      ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = icor11
      icodct(2) = icor22
      icodct(3) = icor33
      icodct(4) = icor12
      icodct(5) = icor13
      icodct(6) = icor23
      icodct(7) = icodce
      icodct(8) = icodcu
      icodct(9) = icodca
      icodct(10) = icodcv
      icodct(11) = icodcw
      nsurij = nsurij + 1
    endif

    if ( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.    &
          icor11.eq.6 .or. icor22.eq.6 .or.                     &
          icor33.eq.6 .or. icor12.eq.6 .or.                     &
          icor13.eq.6 .or. icor23.eq.6 .or.                     &
          icodce.eq.6 .or. icodca.eq.6     ) .and.              &
         (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.    &
          icor11.ne.6 .or. icor22.ne.6 .or.                     &
          icor33.ne.6 .or. icor12.ne.6 .or.                     &
          icor13.ne.6 .or. icor23.ne.6 .or.                     &
          icodce.ne.6 .or. icodca.ne.6     )      ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = icor11
      icodct(2) = icor22
      icodct(3) = icor33
      icodct(4) = icor12
      icodct(5) = icor13
      icodct(6) = icor23
      icodct(7) = icodce
      icodct(8) = icodcu
      icodct(9) = icodca
      icodct(10) = icodcv
      icodct(11) = icodcw
      nsurij = nsurij + 1
    endif

    if ( (icodcu.eq.4 .or. icodcv.eq.4 .or. icodcw.eq.4 .or.    &
          icor11.eq.4 .or. icor22.eq.4 .or.                     &
          icor33.eq.4 .or. icor12.eq.4 .or.                     &
          icor13.eq.4 .or. icor23.eq.4                          &
                                ) .and.                         &
         (icodcu.ne.4 .or. icodcv.ne.4 .or. icodcw.ne.4 .or.    &
          icor11.ne.4 .or. icor22.ne.4 .or.                     &
          icor33.ne.4 .or. icor12.ne.4 .or.                     &
          icor13.ne.4 .or. icor23.ne.4 .or.                     &
          icodce.ne.3                                           &
          .or.icodca.ne.3      ) ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = icor11
      icodct(2) = icor22
      icodct(3) = icor33
      icodct(4) = icor12
      icodct(5) = icor13
      icodct(6) = icor23
      icodct(7) = icodce
      icodct(8) = icodcu
      icodct(9) = icodca
      icodct(10) = icodcv
      icodct(11) = icodcw
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

    if( (icodcu.eq.5 .or. icodcv.eq.5 .or. icodcw.eq.5 .or.       &
         icodck.eq.5 .or. icodce.eq.5 .or. icodcp.eq.5 .or.       &
         icodcf.eq.5 ) .and.                                      &
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.      &
         icodck.ne.5 .or. icodce.ne.5 .or. icodcp.ne.5 .or.       &
         icodcf.ne.5 )                    ) then

      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = ik
      icodct(2) = icodcl(ifac,ik)
      icodct(3) = iep
      icodct(4) = icodcl(ifac,iep)
      icodct(5) = iphi
      icodct(6) = icodcl(ifac,iphi)
      icodct(7) = ifb
      icodct(8) = icodcl(ifac,ifb)
      icodct(9) = icodcu
      icodct(10) = icodcv
      icodct(11) = icodcw
      nstuv2 = nstuv2 + 1

      if( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.     &
           icodck.eq.6 .or. icodce.eq.6 .or. icodcp.eq.6 .or.     &
           icodcf.eq.6 ) .and.                                    &
           (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.    &
           icodck.ne.6 .or. icodce.ne.6 .or. icodcp.ne.6 .or.     &
           icodcf.ne.6 )                    ) then

        if (itypfb(ifac).gt.0) then
          itypfb(ifac) = -itypfb(ifac)
        endif
        icodct(1) = ik
        icodct(2) = icodcl(ifac,ik)
        icodct(3) = iep
        icodct(4) = icodcl(ifac,iep)
        icodct(5) = iphi
        icodct(6) = icodcl(ifac,iphi)
        icodct(7) = ifb
        icodct(8) = icodcl(ifac,ifb)
        icodct(9) = icodcu
        icodct(10) = icodcv
        icodct(11) = icodcw
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
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.    &
         icodck.ne.5 .or. icodce.ne.5 .or. icodcp.ne.5 .or.     &
         icodca.ne.5 )                    ) then

      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = ik
      icodct(2) = icodcl(ifac,ik)
      icodct(3) = iep
      icodct(4) = icodcl(ifac,iep)
      icodct(5) = iphi
      icodct(6) = icodcl(ifac,iphi)
      icodct(7) = ial
      icodct(8) = icodcl(ifac,ial)
      icodct(9) = icodcu
      icodct(10) = icodcv
      icodct(11) = icodcw
      nstuv2 = nstuv2 + 1

      if( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.     &
           icodck.eq.6 .or. icodce.eq.6 .or. icodcp.eq.6 .or.     &
           icodca.eq.6 ) .and.                                    &
           (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.    &
           icodck.ne.6 .or. icodce.ne.6 .or. icodcp.ne.6 .or.     &
           icodca.ne.6 )                    ) then

        if (itypfb(ifac).gt.0) then
          itypfb(ifac) = -itypfb(ifac)
        endif
        icodct(1) = ik
        icodct(2) = icodcl(ifac,ik)
        icodct(3) = iep
        icodct(4) = icodcl(ifac,iep)
        icodct(5) = iphi
        icodct(6) = icodcl(ifac,iphi)
        icodct(7) = ial
        icodct(8) = icodcl(ifac,ial)
        icodct(9) = icodcu
        icodct(10) = icodcv
        icodct(11) = icodcw
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
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.    &
         icodck.ne.5 .or. icodom.ne.5 ) ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = ik
      icodct(2) = icodcl(ifac,ik)
      icodct(3) = iomg
      icodct(4) = icodcl(ifac,iomg)
      icodct(5) = icodcu
      icodct(6) = icodcv
      icodct(7) = icodcw
      nstukw = nstukw + 1
    endif

    if( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.     &
         icodck.eq.6 .or. icodom.eq.6 ) .and.                   &
         (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.    &
         icodck.ne.6 .or. icodom.ne.6 ) ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = ik
      icodct(2) = icodcl(ifac,ik)
      icodct(3) = iomg
      icodct(4) = icodcl(ifac,iomg)
      icodct(5) = icodcu
      icodct(6) = icodcv
      icodct(7) = icodcw
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
         (icodcu.ne.5 .or. icodcv.ne.5 .or. icodcw.ne.5 .or.    &
         icodcn.ne.5 ) ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = inusa
      icodct(2) = icodcl(ifac,inusa)
      icodct(3) = icodcu
      icodct(4) = icodcv
      icodct(5) = icodcw
      nstunu = nstunu + 1
    endif

    if( (icodcu.eq.6 .or. icodcv.eq.6 .or. icodcw.eq.6 .or.     &
         icodcn.eq.6 ) .and.                                    &
         (icodcu.ne.6 .or. icodcv.ne.6 .or. icodcw.ne.6 .or.    &
         icodcn.ne.6 ) ) then
      if (itypfb(ifac).gt.0) then
        itypfb(ifac) = -itypfb(ifac)
      endif
      icodct(1) = inusa
      icodct(2) = icodcl(ifac,inusa)
      icodct(3) = icodcu
      icodct(4) = icodcv
      icodct(5) = icodcw
      nstunu = nstunu + 1
    endif

  enddo
endif

! 2.5 VERIFICATIONS DES COHERENCES INTER VARIABLES
! ================================================

! --- Coherence vitesse scalaires

if( nscal.ge.1 ) then
  do iis = 1, nscal
    if(itytur.eq.2.or.itytur.eq.3) then
      ivar  = isca(iis)
      do ifac = 1, nfabor
        icodcu = icodcl(ifac,iu)
        if(icodcl(ifac,ivar).eq.5.and.icodcu.ne.5) then
          if (itypfb(ifac).gt.0) then
            itypfb(ifac) = -itypfb(ifac)
          endif
          icodus(1) = ivar
          icodus(2) = iis
          icodus(3) = icodcl(ifac,ivar)
          icodus(4) = icodcu
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

if (nstoni.gt.0 .or. nstosc.gt.0 .or. nstovf.gt.0 .or. nstusc.gt.0 ) then
  iok = 1
endif

if (nstvit.gt.0 .or. nstopp.gt.0 .or. nstoke.gt.0 .or. nstrij.gt.0 .or.        &
    nstov2.gt.0 .or. nstonu.gt.0 .or. nstuvw.gt.0 .or. nstoup.gt.0 .or.        &
    nstuke.gt.0 .or. nsurij.gt.0 .or. nstuv2.gt.0 .or. nstunu.gt.0     ) then
  iok = 1
endif

if (irangp.ge.0) call parcmx(iok)

if(iok.ne.0) then

  call sync_bc_err(nstoni, 2, icodni)
  if (nstoni.ne.0) then
    call field_get_label(ivarfl(icodni(1)), chaine)
    write(nfecra,1000) nstoni, chaine, icodni(2)
  endif

  call sync_bc_err(nstvit, 3, icodvi)
  if (nstvit.ne.0) then
    if (icodvi(1) .eq. -6) then
      chaine = 'rugosity'
    else
      call field_get_label(ivarfl(icodvi(1)), chaine)
    endif
    if (icodvi(3).eq.-1) then
      write(nfecra,1010) nstvit, chaine, icodvi(2)
    elseif (icodvi(3).eq.1) then
      write(nfecra,1015) nstvit, chaine, icodvi(2), ippmod(icompf)
    endif
  endif

  if (itytur.eq.2) then

    call sync_bc_err(nstoke, 4, icodtb)
    if (nstoke.ne.0) then
      call field_get_label(ivarfl (icodtb(1)), chaine)
      write(nfecra,1010) nstoke, chaine, icodtb(2)
      call field_get_label(ivarfl (icodtb(3)), chaine)
      write(nfecra,1010) nstoke, chaine, icodtb(4)
    endif

  elseif (itytur.eq.3) then

    call sync_bc_err(nstrij, 3, icodtb)
    if (nstrij.ne.0) then
      call field_get_label(ivarfl (icodtb(1)), chaine)
      write(nfecra,1010) nstrij, chaine, icodtb(2)
      if (iturb.eq.32 .and. icodtb(3) .ge. 0) then
        write(nfecra,2000)
      endif
    endif

  elseif (itytur.eq.5) then

    call sync_bc_err(nstov2, 8, icodtb)
    if (nstov2.ne.0) then
      call field_get_label(ivarfl (icodtb(1)), chaine)
      write(nfecra,1010) nstov2, chaine, icodtb(2)
      call field_get_label(ivarfl (icodtb(3)), chaine)
      write(nfecra,1010) nstov2, chaine, icodtb(4)
      call field_get_label(ivarfl (icodtb(5)), chaine)
      write(nfecra,1010) nstov2, chaine, icodtb(6)
      call field_get_label(ivarfl (icodtb(7)), chaine)
      write(nfecra,1010) nstov2, chaine, icodtb(8)
    endif

  elseif (itytur.eq.6) then

    call sync_bc_err(nstokw, 4, icodtb)
    if (nstokw.ne.0) then
      call field_get_label(ivarfl (icodtb(1)), chaine)
      write(nfecra,1010) nstokw, chaine, icodtb(2)
      call field_get_label(ivarfl (icodtb(3)), chaine)
      write(nfecra,1010) nstokw, chaine, icodtb(4)
    endif

  elseif (itytur.eq.7) then

    call sync_bc_err(nstonu, 2, icodtb)
    if (nstonu.ne.0) then
      call field_get_label(ivarfl (icodtb(1)), chaine)
      write(nfecra,1000) nstonu, chaine, icodtb(2)
    endif

  endif

  call sync_bc_err(nstosc, 2, icodsc)
  if (nstosc.ne.0) then
    if (icodsc(1) .eq. -6) then
      chaine = 'rugosity'
    else
      call field_get_label(ivarfl (icodsc(1)), chaine)
    endif
    write(nfecra,1010) nstosc, chaine, icodsc(2)
  endif

  call sync_bc_err(nstovf, 2, icodvf)
  if (nstovf.ne.0) then
    call field_get_label(ivarfl (icodvf(1)), chaine)
    write(nfecra,1010) nstovf, chaine, icodvf(2)
  endif

  call sync_bc_err(nstuvw, 3, icoduv)
  if (nstuvw.ne.0) then
    write(nfecra,1020) nstuvw, icoduv(1), icoduv(2), icoduv(3)
  endif

  if (itytur.eq.2) then

    call sync_bc_err(nstuke, 7, icodct)
    if (nstuke.ne.0) then
      call field_get_label(ivarfl (icodct(1)), chaine)
      write(nfecra,1030) nstuke,  chaine, icodct(2),             &
                         icodct(5), icodct(6), icodct(7)
      call field_get_label(ivarfl (icodct(3)), chaine)
      write(nfecra,1030) nstuke, chaine, icodct(4),              &
                         icodct(5), icodct(6), icodct(7)
    endif

  elseif (iturb.eq.30 .or. iturb.eq.31) then

    call sync_bc_err(nsurij, 10, icodct)
    if (nsurij.ne.0) then
      write(nfecra,1040) nsurij, icodct(1), icodct(2),           &
                         icodct(3), icodct(4), icodct(5),        &
                         icodct(6), icodct(7), icodct(8),        &
                         icodct(9), icodct(10)
    endif

  elseif (iturb.eq.32) then

    call sync_bc_err(nsurij, 11, icodct)
    if (nsurij.ne.0) then
      write(nfecra,1041) nsurij, icodct(1), icodct(2),           &
                         icodct(3), icodct(4), icodct(5),        &
                         icodct(6), icodct(7), icodct(8),        &
                         icodct(9), icodct(10), icodct(11)
    endif

  elseif (itytur.eq.5) then

    call sync_bc_err(nstuv2, 11, icodct)
    if (nstuv2.ne.0) then
      call field_get_label(ivarfl (icodct(1)), chaine)
      write(nfecra,1030) nstuv2, chaine, icodct(2),              &
                         icodct(9), icodct(10), icodct(11)
      call field_get_label(ivarfl (icodct(3)), chaine)
      write(nfecra,1030) nstuv2, chaine, icodct(4),              &
                         icodct(9), icodct(10), icodct(11)
      call field_get_label(ivarfl (icodct(5)), chaine)
      write(nfecra,1030) nstuv2, chaine, icodct(6),              &
                         icodct(9), icodct(10), icodct(11)
      call field_get_label(ivarfl (icodct(7)), chaine)
      write(nfecra,1030) nstuv2, chaine, icodct(8),              &
                         icodct(9), icodct(10), icodct(11)
    endif

  elseif (itytur.eq.6) then

    call sync_bc_err(nstukw, 7, icodct)
    if (nstukw.ne.0) then
      call field_get_label(ivarfl (icodct(1)), chaine)
      write(nfecra,1030) nstukw, chaine, icodct(2),              &
                         icodct(5), icodct(6), icodct(7)
      call field_get_label(ivarfl (icodct(3)), chaine)
      write(nfecra,1030) nstukw, chaine, icodct(4),              &
                         icodct(5), icodct(6), icodct(7)
    endif

  elseif (itytur.eq.7) then

    call sync_bc_err(nstunu, 5, icodct)
    if (nstunu.ne.0) then
      call field_get_label(ivarfl (icodct(1)), chaine)
      write(nfecra,1030) nstunu, chaine, icodct(2),              &
                         icodct(3), icodct(4), icodct(5)
    endif
  endif

  call sync_bc_err(nstusc, 4, icodus)
  if (nstusc.ne.0) then
    call field_get_label(ivarfl (icodus(1)), chaine)
    write(nfecra,1050) nstusc, chaine,                           &
         icodus(2), icodus(3), icodus(4)
  endif

  if (nstoni.gt.0 .or. nstosc.gt.0 .or. nstovf.gt.0 .or. nstusc.gt.0 ) then
    write (nfecra,1901) nstoni, nstosc, nstovf, nstusc
  endif

  if (nstvit.gt.0 .or. nstopp.gt.0 .or. nstoke.gt.0 .or. nstrij.gt.0 .or.    &
      nstov2.gt.0 .or. nstonu.gt.0 .or. nstuvw.gt.0 .or. nstoup.gt.0 .or.    &
      nstuke.gt.0 .or. nsurij.gt.0 .or. nstuv2.gt.0 .or. nstunu.gt.0     ) then
    write (nfecra,1902)  nstvit,nstopp, nstoke,nstrij,                       &
                         nstov2,nstonu, nstuvw,nstoup,                       &
                         nstuke,nsurij, nstuv2,nstunu
  endif

  call boundary_conditions_error(itypfb)

endif

!===============================================================================
! 3.  FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@                                                            ',/,&
'@ COND. LIM. NON INITIALISEES                                ',/,&
'@   Nombre de faces de bord ',i10   ,'; variable ',a16        ,/,&
'@     icodcl variable derniere face ', i10                    ,/,&
'@                                                            '  )
 1010 format(                                                     &
'@                                                            ',/,&
'@ COND. LIM. NON PREVUES                                     ',/,&
'@   Nombre de faces de bord ',i10   ,'; variable ',a16        ,/,&
'@     icodcl variable derniere face ', i10                    ,/,&
'@                                                            '  )
 1015 format(                                                     &
'@                                                            ',/,&
'@ CONDITIONS AUX LIMITES DE PAROI RUGUEUSE INCOMPATIBLES     ',/,&
'@ AVEC LE MODULE COMPRESSIBLE                                ',/,&
'@   Nombre de faces de bord ',i10   ,'; variable ',a16        ,/,&
'@     icodcl variable derniere face ', i10                    ,/,&
'@     ippmod(icompf)', i10                                    ,/,&
'@                                                            '  )
 1020 format(                                                     &
'@                                                            ',/,&
'@ INCOHERENCE COND. LIM. COMPOSANTES DE LA VITESSE           ',/,&
'@   Nombre de faces de bord ',i10                             ,/,&
'@     icodcl derniere face ', 3i10                            ,/,&
'@                                                            '  )
 1030 format(                                                     &
'@                                                            ',/,&
'@ INCOHERENCE COND. LIM. VITESSE-VARIABLE                    ',/,&
'@   Nombre de faces de bord ',i10   ,'; variable ',a16        ,/,&
'@     icodcl variable ', i10                                  ,/,&
'@     icodcl vitesse  ',3i10                                  ,/,&
'@                                                            '  )
 1040 format(                                                     &
'@                                                            ',/,&
'@ INCOHERENCE COND. LIM. VITESSE-RIJ-EPSILON                 ',/,&
'@   Nombre de faces de bord ',i10                             ,/,&
'@     icodcl derniere face Rij-eps ', 7i5                     ,/,&
'@     icodcl               vitesse ', 3i5                     ,/,&
'@                                                            '  )
 1041 format(                                                     &
'@                                                            ',/,&
'@ INCOHERENCE COND. LIM. VITESSE-RIJ-EPSILON EBRSM           ',/,&
'@   Nombre de faces de bord ',i10                             ,/,&
'@     icodcl derniere face Rij-eps ', 8i5                     ,/,&
'@     icodcl               vitesse ', 3i5                     ,/,&
'@                                                            '  )
 1050 format(                                                     &
'@                                                            ',/,&
'@ INCOHERENCE COND. LIM. VITESSE-SCALAIRE                    ',/,&
'@   Nombre de faces de bord ',i10   ,'; variable ',a16        ,/,&
'@     derniere face : '                                       ,/,&
'@       scalaire numero ',i10                                 ,/,&
'@       icodcl scalaire ',i10   ,'; icodcl vitesse ',i10      ,/,&
'@                                                            '  )
 1901 format(                                                     &
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
'@           ou cs_user_boundary_conditions.                  ',/,&
'@                                                            ',/)
 1902 format(                                                     &
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
'@           ou cs_user_boundary_conditions.                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@                CONDITION DE PAROI RUGUEUSE CHOISIE         ',/,&
'@   CONDITION LIMITE INCOMPATIBLE AVEC LE MODELE EBRSM       ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute                             ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou        ',/,&
'@    cs_user_boundary_conditions :                           ',/,&
'@      modifier la condition de paroi rugueuse en            ',/,&
'@               paroi lisse                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2001 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@                CONDITION DE SURFACE LIBRE EN ALE.          ',/,&
'@   CONDITION LIMITE INCOMPATIBLE SANS GRAVITE               ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute                             ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@      modifier la gravite                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                     &
'@                                                            ',/,&
'@ UNINITIALIZED BOUNDARY CONDITIONS                          ',/,&
'@   Number of boundary faces ',i10   ,'; variable ',a16       ,/,&
'@     icodcl variable last face ', i10                        ,/,&
'@                                                            '  )
 1010 format(                                                     &
'@                                                            ',/,&
'@ UNEXPECTED BOUNDARY CONDITIONS                             ',/,&
'@   Number of boundary faces ',i10   ,'; variable ',a16       ,/,&
'@     icodcl variable last face ', i10                        ,/,&
'@                                                            '  )
 1015 format(                                                     &
'@                                                            ',/,&
'@ ROUGH WALL BOUNDARY CONDITIONS INCOMPATIBLE WITH THE       ',/,&
'@ COMPRESSIBLE MODULE                                        ',/,&
'@   Number of boundary faces ',i10   ,'; variable ',a16       ,/,&
'@     icodcl variable last face ', i10                        ,/,&
'@     ippmod(icompf)', i10                                    ,/,&
'@                                                            '  )
 1020 format(                                                     &
'@                                                            ',/,&
'@ INCOHERENT BOUNDARY CONDITIONS VELOCITY COMPONENT          ',/,&
'@   Number of boundary faces ',i10                            ,/,&
'@     icodcl last face ', 3i10                                ,/,&
'@                                                            '  )
 1030 format(                                                     &
'@                                                            ',/,&
'@ INCOHERENCY BOUNDARY CONDITIONS VELOCITY-VARIABLE          ',/,&
'@   Number of boundary faces ',i10   ,'; variable ',a16       ,/,&
'@     icodcl variable last face ', i10                        ,/,&
'@     icodcl velocity ',3i10                                  ,/,&
'@                                                            '  )
 1040 format(                                                     &
'@                                                            ',/,&
'@ INCOHERENT BOUNDARY CONDITIONS VELOCITY-RIJ-EPSILON        ',/,&
'@   Number of boundary faces ',i10                            ,/,&
'@     icodcl last face Rij-eps  ', 7i5                        ,/,&
'@     icodcl           velocity ', 3i5                        ,/,&
'@                                                            '  )
 1041 format(                                                     &
'@                                                            ',/,&
'@ INCOHERENT BOUNDARY CONDITIONS VELOCITY-RIJ-EPSILON EBRSM  ',/,&
'@   Number of boundary faces ',i10                            ,/,&
'@     icodcl last face Rij-eps  ', 8i5                        ,/,&
'@     icodcl           velocity ', 3i5                        ,/,&
'@                                                            '  )
 1050 format(                                                     &
'@                                                            ',/,&
'@ INCOHERENT BOUNDARY CONDITIONS VELOCITY-SCALAR             ',/,&
'@   Number of boundary faces ',i10   ,'; variable ',a16       ,/,&
'@     last face:'                                             ,/,&
'@       scalar number ',i10                                   ,/,&
'@       icodcl scalar ',i10   ,'; icodcl velocity ',i10       ,/,&
'@                                                            '  )
 1901 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE BOUNDARY CONDITIONS CHECK     ',/,&
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
'@           cs_user_boundary_conditions.                     ',/,&
'@                                                            ',/)
 1902 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE BOUNDARY CONDITIONS CHECK     ',/,&
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
'@           cs_user_boundary_conditions.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE BOUNDARY CONDITIONS CHECK     ',/,&
'@    =========                                               ',/,&
'@             ROUGH WALL BOUNDARY CONDITIONS INCOMPATIBLE    ',/,&
'@                   WITH EBRSM MODEL                         ',/,&
'@                                                            ',/,&
'@         The calculation will not be run.                   ',/,&
'@                                                            ',/,&
'@         Verify the parameters given via the interface or   ',/,&
'@           cs_user_boundary_conditions.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
2001 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE BOUNDARY CONDITIONS CHECK     ',/,&
'@    =======                                                 ',/,&
'@             FREE-SURFACE CONDITION IN ALE                  ',/,&
'@   MUST BE COMBINED WITH A NON ZERO GRAVITY TERM!           ',/,&
'@                                                            ',/,&
'@         The calculation will not be run.                   ',/,&
'@                                                            ',/,&
'@      Verify the parameters or modify the gravity.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


#endif

return
end subroutine vericl

!===============================================================================
! Local functions
!===============================================================================

!> \brief synchronize boundary condition error logging across MPI ranks.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in, out] nerloc       number of errors (local rank in, global out)
!> \param[in]      nerrcd       number of codes saved at error faces
!> \param[in, out] errcod       codes saved at one error face (local in,
!                               broadcast out)
!_______________________________________________________________________________

subroutine sync_bc_err &
 ( nerloc , nerrcd , errcod )

!===============================================================================
! Module files
!===============================================================================

use parall

!===============================================================================

implicit none

! Arguments

integer nerloc, nerrcd
integer errcod(nerrcd)

! Local variables

integer irkerr

!===============================================================================

if (irangp.ge.0) then
  irkerr = -1
  if (nerloc.gt.0) irkerr = irangp
  call parcpt(nerloc)
  if (nerloc .ne. 0) then
    call parcmx(irkerr)
    call parbci(irkerr, nerrcd, errcod)
  endif
endif

return
end subroutine sync_bc_err
