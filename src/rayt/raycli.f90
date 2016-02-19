!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

!> \file raycli.f90
!>
!> \brief Compute wall temperature for radiative transfer, and update BCs
!>
!> 1) Compute wall temperature for radiative transfer
!>
!> 2) Update BCs for the energy computation
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
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
!> \param[in]     itypfb        face boundary condition type
!> \param[in]     izfrad        zone index for boundary faces
!>                               and reference face index
!> \param[in]     dt            time step (per cell)
!> \param[in]     propce        physical properties at cell centers
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


subroutine raycli &
 ( nvar   , nscal  ,                                              &
   icodcl , itypfb ,                                              &
   izfrad ,                                                       &
   dt     , propce , rcodcl )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ihmpre
use ppppar
use ppthch
use ppincl
use ppcpfu
use radiat
use dimens, only: ndimfb
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(ndimfb,nvarcl)
integer          itypfb(ndimfb)
integer          izfrad(ndimfb)

double precision dt(ncelet)
double precision propce(ncelet,*)
double precision rcodcl(ndimfb,nvarcl,3)

! Local variables

integer          iok, ifac, iel, ideb, ivart
integer          mode, ifvu, ii, izonem, izone
integer          nrferr(15), icoerr(16)
integer          ivahg, nlst, icodw

double precision tmin , tmax   , tx
double precision xmtk
double precision rvferr(27)

integer, allocatable, dimension(:) :: isothm, lstfac

double precision, allocatable, dimension(:) :: tempk, thwall, tparo
double precision, allocatable, dimension(:) :: text, tint
double precision, dimension(:), pointer :: bhconv, bfconv
double precision, dimension(:), pointer :: b_temp, bqinci
double precision, dimension(:), pointer :: bxlam, bepa, beps, bfnet

double precision, dimension(:), pointer :: cvara_scalt

integer    ipacli
data       ipacli /0/
save       ipacli

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine b_t_to_h(nlst, lstfac, t_b, h_b)

    use mesh, only: nfabor
    implicit none

    integer :: nlst
    integer, dimension(nlst) :: lstfac
    double precision, dimension(nfabor), intent(in) :: t_b
    double precision, dimension(nfabor), intent(out), target :: h_b

  end subroutine b_t_to_h

 end interface

!===============================================================================
! 1. Initializations
!===============================================================================

! Allocate temporary arrays
allocate(isothm(nfabor), lstfac(nfabor))
allocate(tempk(ncelet), thwall(ndimfb))
allocate(text(nfabor), tint(nfabor))
allocate(tparo(nfabor))

! Map field arrays
call field_get_val_s(itempb, b_temp)
call field_get_val_s(iqinci, bqinci)
call field_get_val_s(ixlam, bxlam)
call field_get_val_s(iepa, bepa)
call field_get_val_s(ieps, beps)
call field_get_val_s(ifnet, bfnet)

!---> NUMERO DE PASSAGE RELATIF
ipacli = ipacli + 1
ideb = 0

!---> Min and Max values of the temperature (in Kelvin)
tmin = 0.d0
tmax = grand + tkelvi

!---> COEFF DE RELAX

!     tx est strictement superieur a 0 et inferieur ou egal a 1

!     Pour calculer la temperature de paroi, on calcule un increment
!     de temperature DeltaT entre l'etape courante n et l'etape
!     precedente n-1, puis on calcule :
!          n    n-1                                 n-1
!         T  = T    + DeltaT si le rapport DeltaT/T    =< tx, sinon

!          n    n-1                      n-1             n-1
!         T  = T    * (1 + tx *((DeltaT/T   ) / |DeltaT/T   |))

tx = 0.1d0

! Temperature scale

if (itpscl.eq.2) then
  xmtk = -tkelvi
else if (itpscl.eq.1) then
  xmtk = 0.d0
endif

! Wall temperature

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
    tparo(ifac) = b_temp(ifac) - xmtk
  else
    tparo(ifac) = 0.d0
  endif
enddo

!---> Default initialization
do ifac = 1, nfabor
  izfrad(ifac) = -1
  isothm(ifac) = -1
  bxlam(ifac) = -grand
  bepa(ifac)  = -grand
  beps(ifac)  = -grand
  text  (ifac) = -grand
  tint  (ifac) = -grand
enddo

! Index of the thermal variable
ivart = isca(iscalt)
if (ihgas.ge.0) ivahg = isca(ihgas)

! Pointers to specific fields
if (ifconv.ge.0) call field_get_val_s(ifconv, bfconv)
if (ihconv.ge.0) call field_get_val_s(ihconv, bhconv)

! Error checking

do ii = 1, 15
  nrferr(ii) = 0
enddo

!===============================================================================
! 2. SI PAS DE FICHIER SUITE ALORS INITIALISATION AU PREMIER PASSAGE
!    DE TPAROI ET QINCID :
!      LECTURE DE L'INITIALISATION DE TPAROI A TINT
!      QINCID EST INITIALISE A STEPHN*TINT**4 (SI ON INITIALISE QINCID
!      A ZERO, ON AURA UN DEFICIT SUR LA CONDITION LIMITE DE LUMINANCE
!      AUX PAROIS AU 1er PAS DE TEMPS EN DOM)
!===============================================================================

if (ipacli.eq.1 .and. isuird.eq.0) then

  ! Indicateur : si non suite et premier pas de temps.
  ideb = 1

  do iel = 1,ncelet
    propce(iel,ipproc(itsri(1))) = zero
    propce(iel,ipproc(itsre(1))) = zero
  enddo

  do ifac = 1,nfabor
    bhconv(ifac) = zero
    bfconv(ifac) = zero
  enddo

  !     On utilise TBORD comme auxiliaire pour l'appel a USRAY2
  !       pour etre sur que TPAROI ne sera pas modifie
  !       (puisqu'on a TBORD libre)
  !     On utilise FLUNET comme auxiliaire pour l'appel a USRAY2
  !       pour etre sur que QINCID ne sera pas modifie
  !       (puisqu'on a FLUNET libre)

  do ifac = 1, nfabor
    thwall(ifac) = zero
    bfnet(ifac) = zero
  enddo

  ! - Interface Code_Saturne
  !   ======================

  if (iihmpr.eq.1) then

    !---> NUMERO DU SCALAIRE ET DE LA VARIABLE THERMIQUE

    call uiray2 &
   ( itypfb, iparoi, iparug, ivart , izfrad,                  &
     isothm, itpimp, ipgrno, iprefl, ifgrno, ifrefl,          &
     nozppm, nfabor, nvar,                                    &
     beps, bepa,                                              &
     tint, text,                                              &
     bxlam, rcodcl)

  endif

  call cs_user_radiative_transfer_bcs &
( nvar   , nscal  ,                                              &
  itypfb ,                                                       &
  icodcl , izfrad , isothm ,                                     &
  tmin   , tmax   , tx     ,                                     &
  dt     , rcodcl ,                                              &
  thwall , bfnet  , bhconv , bfconv ,                            &
  bxlam  , bepa   , beps   ,                                     &
  text   , tint   )

  write(nfecra,1000)

  ! Tparoi en Kelvin et QINCID en W/m2
  do ifac = 1, nfabor
    if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then
      tparo(ifac) = tint(ifac)
      bqinci(ifac) = stephn*tint(ifac)**4
    else
      tparo(ifac) = 0.d0
      bqinci(ifac) = 0.d0
    endif
  enddo

endif

!===============================================================================
! 3.1 DONNEES SUR LES FACES FRONTIERES
!===============================================================================

!     On utilise flunet comme auxiliaire pour l'appel a USRAY2
!       pour etre sur que QINCID ne sera pas modifie
!       (puisqu'on a flunet libre)

do ifac = 1, nfabor
  thwall (ifac) = tparo(ifac)
  bfnet(ifac) = bqinci(ifac)
enddo

!     - Interface Code_Saturne
!       ======================

if (iihmpr.eq.1) then

  call uiray2 &
( itypfb, iparoi, iparug, ivart , izfrad,                       &
  isothm, itpimp, ipgrno, iprefl, ifgrno, ifrefl,               &
  nozppm, nfabor, nvar,                                         &
  beps, bepa, tint, text,                                       &
  bxlam, rcodcl)

endif

call cs_user_radiative_transfer_bcs &
( nvar   , nscal  ,                                              &
  itypfb ,                                                       &
  icodcl , izfrad , isothm ,                                     &
  tmin   , tmax   , tx     ,                                     &
  dt     , rcodcl ,                                              &
  thwall , bfnet ,  bhconv ,                                     &
  bfconv , bxlam , bepa   , beps ,                               &
  text   , tint   )

!===============================================================================
! 3.2 CONTROLE DES DONNEES UTILISATEUR
!===============================================================================

!--> Arret si le numero de zone est non renseigne ou mal renseigne

do ifac = 1, nfabor
  if (izfrad(ifac).le.0.or.izfrad(ifac).gt.nozrdm) then
    nrferr(1) = nrferr(1) + 1
    icoerr(1) = izfrad(ifac)
    itypfb(ifac) = - iabs(itypfb(ifac))
  endif
enddo

! --> On construit une liste des numeros des zones frontieres.
!           (liste locale au processeur, en parallele)
!     Erreur si depassement.

nzfrad = 0
do ifac = 1, nfabor
  ifvu = 0
  do ii = 1, nzfrad
    if (ilzrad(ii).eq.izfrad(ifac)) then
      ifvu = 1
    endif
  enddo
  if (ifvu.eq.0) then
    nzfrad = nzfrad + 1
    if (nzfrad.le.nbzrdm) then
      ilzrad(nzfrad) = izfrad(ifac)
    else
      nrferr(2) = nrferr(2) + 1
    endif
  endif
enddo

! ---> Plus grand numero de zone atteint

izonem = 0
do ii = 1, nzfrad
  izone = ilzrad(ii)
  izonem = max(izonem,izone)
enddo
if (irangp.ge.0) then
  call parcmx(izonem)
endif
nozarm = izonem

! On verra si ca coute cher ou non.
!   Pour le moment on le fait tout le temps.
!        IF(IWARNI(IVART).GE.-1.OR.IPACLI.LE.3) THEN
if (1.eq.1) then

  !--> Si en paroi ISOTHM non renseignee : stop
  do ifac = 1, nfabor
    if ((itypfb(ifac).eq.iparoi  .or.                            &
         itypfb(ifac).eq.iparug) .and.                           &
         isothm(ifac).eq.-1) then
      nrferr(3) = nrferr(3) + 1
      icoerr(3) = izfrad(ifac)
      itypfb(ifac) = - iabs(itypfb(ifac))
    endif
  enddo

  !--> Si ISOTHM renseignee en non paroi : stop
  do ifac = 1, nfabor
    if (itypfb(ifac).ne.iparoi .and.                             &
        itypfb(ifac).ne.iparug .and.                             &
        isothm(ifac).ne.-1         ) then
      nrferr(4) = nrferr(4) + 1
      icoerr(4) = izfrad(ifac)
      itypfb(ifac) = - iabs(itypfb(ifac))
    endif
  enddo

  !--> Si valeur physique erronee : stop
  do ifac = 1, nfabor
    if (isothm(ifac).eq.itpimp) then
      if (beps(ifac).lt.0.d0.or.                                 &
          beps(ifac).gt.1.d0.or.                                 &
          tint(ifac).le.0.d0                      ) then
        nrferr(5) = nrferr(5) + 1
        icoerr(5) = izfrad(ifac)
        rvferr(1) = beps(ifac)
        rvferr(2) = tint(ifac)
        itypfb(ifac) = - iabs(itypfb(ifac))
      endif
    elseif (isothm(ifac).eq.ipgrno) then
      if (beps(ifac) .lt.0.d0.or.                                &
          beps(ifac) .gt.1.d0.or.                                &
          bxlam(ifac).le.0.d0.or.                                &
          bepa(ifac) .le.0.d0.or.                                &
          text(ifac).le.0.d0.or.                                 &
          tint(ifac).le.0.d0                      ) then
        nrferr(6) = nrferr(6) + 1
        icoerr(6) = izfrad(ifac)
        rvferr(3) = beps(ifac)
        rvferr(4) = bxlam(ifac)
        rvferr(5) = bepa(ifac)
        rvferr(6) = text(ifac)
        rvferr(7) = tint(ifac)
        itypfb(ifac) = - iabs(itypfb(ifac))
      endif
    elseif (isothm(ifac).eq.iprefl) then
      if (bxlam(ifac).le.0.d0.or.                                &
          bepa(ifac) .le.0.d0.or.                                &
          text(ifac).le.0.d0.or.                                 &
          tint(ifac).le.0.d0                      ) then
        nrferr(7) = nrferr(7) + 1
        icoerr(7) = izfrad(ifac)
        rvferr(8) = bxlam(ifac)
        rvferr(9) = bepa(ifac)
        rvferr(10) = text(ifac)
        rvferr(11) = tint(ifac)
        itypfb(ifac) = - iabs(itypfb(ifac))
      endif
    elseif (isothm(ifac).eq.ifgrno) then
      if (beps(ifac).lt.0.d0.or.                                 &
          beps(ifac).gt.1.d0.or.                                 &
          tint(ifac).le.0.d0                      ) then
        nrferr(8) = nrferr(8) + 1
        icoerr(8) = izfrad(ifac)
        rvferr(12) = beps(ifac)
        rvferr(13) = tint(ifac)
        itypfb(ifac) = - iabs(itypfb(ifac))
      endif
    elseif (isothm(ifac).eq.ifrefl) then
      if (tint(ifac).le.0.d0) then
        nrferr(9) = nrferr(9) + 1
        icoerr(9) = izfrad(ifac)
        rvferr(14) = tint(ifac)
        itypfb(ifac) = - iabs(itypfb(ifac))
      endif
    elseif (isothm(ifac).eq.itpt1d) then
      if (beps(ifac).lt.0.d0.or.                                 &
          beps(ifac).gt.1.d0.or.                                 &
          tint(ifac).le.0.d0) then
        nrferr(15) = nrferr(15) + 1
        icoerr(16) = izfrad(ifac)
        rvferr(26) = beps(ifac)
        rvferr(27) = tint(ifac)
        itypfb(ifac) = - iabs(itypfb(ifac))
      endif
    elseif (isothm(ifac).ne.-1) then
      nrferr(10) = nrferr(10) + 1
      icoerr(10) = izfrad(ifac)
      icoerr(11) = isothm(ifac)
      itypfb(ifac) = - iabs(itypfb(ifac))
    endif
  enddo

  !--> Si valeur renseignee sans raison : stop
  do ifac = 1, nfabor
   if (isothm(ifac).eq.itpimp) then
      if (bxlam(ifac).gt.0.d0.or.                &
          bepa(ifac) .gt.0.d0.or.                &
          text(ifac).gt.0.d0                      ) then
        nrferr(11) = nrferr(11) + 1
        icoerr(12) = izfrad(ifac)
        rvferr(15) = bxlam(ifac)
        rvferr(16) = bepa(ifac)
        rvferr(17) = text(ifac)
        itypfb(ifac) = - iabs(itypfb(ifac))
      endif
    elseif (isothm(ifac).eq.iprefl) then
      if (beps(ifac).ge.0.d0) then
        nrferr(12) = nrferr(12) + 1
        icoerr(13) = izfrad(ifac)
        rvferr(18) = beps(ifac)
        itypfb(ifac) = - iabs(itypfb(ifac))
      endif
    elseif (isothm(ifac).eq.ifgrno) then
      if (bxlam(ifac).gt.0.d0.or.                &
          bepa(ifac) .gt.0.d0.or.                &
          text(ifac).gt.0.d0                      ) then
        nrferr(13) = nrferr(13) + 1
        icoerr(14) = izfrad(ifac)
        rvferr(19) = bxlam(ifac)
        rvferr(20) = bepa(ifac)
        rvferr(21) = text(ifac)
        itypfb(ifac) = - iabs(itypfb(ifac))
      endif
    elseif (isothm(ifac).eq.ifrefl) then
      if (beps(ifac) .ge.0.d0.or.                &
          bxlam(ifac).gt.0.d0.or.                &
          bepa(ifac) .gt.0.d0.or.                &
          text(ifac).gt.0.d0                      ) then
        nrferr(14) = nrferr(14) + 1
        icoerr(15) = izfrad(ifac)
        rvferr(22) = beps(ifac)
        rvferr(23) = bxlam(ifac)
        rvferr(24) = bepa(ifac)
        rvferr(25) = text(ifac)
        itypfb(ifac) = - iabs(itypfb(ifac))
      endif
    endif
  enddo

endif

!===============================================================================
! Error logging
!===============================================================================

iok = 0

do ii = 1, 15
  if (nrferr(ii).gt.0) iok = 1
enddo

if (irangp.ge.0) call parcmx(iok)

if (iok.ne.0) then

  call sync_bc_err(nrferr(1), 1, icoerr(1:1))
  if (nrferr(1).gt.0) then
    write(nfecra,2000) nozrdm, nrferr(1), icoerr(1)
  endif

  call sync_bc_err(nrferr(2), nbzrdm, ilzrad)
  if (nrferr(2).gt.0) then
    write(nfecra,2001) nbzrdm
    write(nfecra,2002)(ilzrad(ii),ii=1,nbzrdm)
    call csexit (1)
  endif

  call sync_bc_err(nrferr(3), 1, icoerr(3:3))
  if (nrferr(3).gt.0) then
    write(nfecra,2110) nrferr(3), icoerr(3)
  endif

  call sync_bc_err(nrferr(4), 1, icoerr(4:4))
  if (nrferr(4).gt.0) then
    write(nfecra,2111) nrferr(4), icoerr(4)
  endif

  call sync_rad_bc_err(nrferr(5), 2, icoerr(5), rvferr(1:2))
  if (nrferr(5).gt.0) then
    write(nfecra,2120) nrferr(5), icoerr(5), rvferr(1), rvferr(2)
  endif

  call sync_rad_bc_err(nrferr(6), 5, icoerr(6), rvferr(3:7))
  if (nrferr(6).gt.0) then
    write(nfecra,2130) nrferr(6), icoerr(6), rvferr(3), rvferr(4), &
                       rvferr(5), rvferr(6), rvferr(7)
  endif

  call sync_rad_bc_err(nrferr(7), 4, icoerr(7), rvferr(8:11))
  if (nrferr(7).gt.0) then
    write(nfecra,2140) nrferr(7), icoerr(7), rvferr(8), rvferr(9), &
                       rvferr(10), rvferr(11)
  endif

  call sync_rad_bc_err(nrferr(8), 2, icoerr(8), rvferr(12:13))
  if (nrferr(8).gt.0) then
    write(nfecra,2150) nrferr(8), icoerr(8), rvferr(12), rvferr(13)
  endif

  call sync_rad_bc_err(nrferr(9), 1, icoerr(9), rvferr(14:14))
  if (nrferr(9).gt.0) then
    write(nfecra,2160) nrferr(9), icoerr(9), rvferr(14)
  endif

  call sync_rad_bc_err(nrferr(15), 2, icoerr(15), rvferr(26:27))
  if (nrferr(15).gt.0) then
    write(nfecra,2170) nrferr(15), icoerr(16), rvferr(26), rvferr(27)
  endif

  call sync_bc_err(nrferr(10), 2, icoerr(10:11))
  if (nrferr(10).gt.0) then
    write(nfecra,2180) nrferr(10), icoerr(10), icoerr(11)
  endif

  call sync_rad_bc_err(nrferr(11), 3, icoerr(12), rvferr(15:17))
  if (nrferr(11).gt.0) then
    write(nfecra,2220) nrferr(11), icoerr(12), rvferr(15), rvferr(16), &
                       rvferr(17)
  endif

  call sync_rad_bc_err(nrferr(12), 1, icoerr(13), rvferr(18:18))
  if (nrferr(12).gt.0) then
    write(nfecra,2240) nrferr(12), icoerr(13), rvferr(18)
  endif

  call sync_rad_bc_err(nrferr(13), 3, icoerr(14), rvferr(19:21))
  if (nrferr(13).gt.0) then
    write(nfecra,2250) nrferr(13), icoerr(14), rvferr(19), rvferr(20), &
                       rvferr(21)
  endif

  call sync_rad_bc_err(nrferr(14), 4, icoerr(15), rvferr(22:25))
  if (nrferr(14).gt.0) then
    write(nfecra,2260) nrferr(14), icoerr(15), rvferr(22), rvferr(23), &
                       rvferr(24), rvferr(25)
  endif

  call boundary_conditions_error(itypfb)

endif

!===============================================================================
! 3.2 COMPLETION DES DONNEES UTILISATEUR
!===============================================================================

! ICODCL et EPS (quand il est nul)

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparug) then
    icodw = 6
  else
    icodw = 5
  endif
  if (isothm(ifac).eq.itpimp) then
    icodcl(ifac,ivart) = icodw
    if (ihgas.ge.0) icodcl(ifac,ivahg) = icodw
  elseif (isothm(ifac).eq.ipgrno) then
    icodcl(ifac,ivart) = icodw
    if (ihgas.ge.0) icodcl(ifac,ivahg) = icodw
  elseif (isothm(ifac).eq.iprefl) then
    icodcl(ifac,ivart) = icodw
    if (ihgas.ge.0) icodcl(ifac,ivahg) = icodw
    beps(ifac) = 0.d0
  elseif (isothm(ifac).eq.ifgrno) then
    icodcl(ifac,ivart) = icodw
    if (ihgas.ge.0) icodcl(ifac,ivahg) = icodw
  elseif (isothm(ifac).eq.ifrefl) then
    icodcl(ifac,ivart) = 3
    if (ihgas.ge.0) icodcl(ifac,ivahg) = 3
    beps(ifac) = 0.d0
  elseif (isothm(ifac).eq.itpt1d) then
    icodcl(ifac,ivart) = icodw
    if (ihgas.ge.0) icodcl(ifac,ivahg) = icodw
  endif
enddo

!===============================================================================
! 4. STOCKAGE DE LA TEMPERATURE (en Kelvin) dans TEMPK(IEL)
!===============================================================================

call field_get_val_prev_s(ivarfl(isca(iscalt)), cvara_scalt)

if (itherm.eq.1) then

  if (itpscl.eq.2) then
    do iel = 1, ncel
      tempk(iel) = cvara_scalt(iel) + tkelvi
    enddo
  else if (itpscl.eq.1) then
    do iel = 1, ncel
      tempk(iel) = cvara_scalt(iel)
    enddo
  endif

elseif (itherm.eq.2) then

  call c_h_to_t(cvara_scalt, tempk) ! Resultat : T en K

endif

!===============================================================================
! 5. CALCUL DES TEMPERATURES DE PAROIS
!===============================================================================

! DANS TOUS LES CAS HFCONV CONTIENT Lambda * Hturb / distance
!   (HFCONV : W/(m2 K) ; Hturb est sans dimension)
!  (au premier passage, il est nul)

!--> CALCUL DU FLUX CONVECTIF
!      Par flux convectif, on entend bien sur
!        flux convectif parallele a la paroi,
!        on suppose que la paroi est etanche...
!      Le flux est calcule dans condli clptur, sauf au premier
!        passage sans suite de calcul, puisque raycli est appele avant.

if (ideb.eq.1) then

  do ifac = 1, nfabor
    if (isothm(ifac).ne.-1) then
      bfconv(ifac) = bhconv(ifac)*(tempk(ifabor(ifac))-tparo(ifac))
    endif
  enddo

endif

!--> Les cas ou il faut calculer TPAROI sont, au premier passage sans suite
!      des cas a temperature imposee TPAROI = TINT

if (ideb.eq.1) then

  do ifac = 1,nfabor
    if (isothm(ifac).eq.ipgrno .or.                             &
        isothm(ifac).eq.iprefl .or.                             &
        isothm(ifac).eq.ifgrno .or.                             &
        isothm(ifac).eq.itpt1d) then
      isothm(ifac) = itpimp
    endif
  enddo

endif

if (ideb.eq.0) then

  call raypar &
 ( isothm , izfrad ,                                              &
   tmin   , tmax   , tx     ,                                     &
   rcodcl ,                                                       &
   tparo  , bqinci , text   , tint   ,                            &
   bxlam  , bepa   , beps   , bhconv ,                            &
   bfconv , tempk  )

endif

!===============================================================================
! 6.  CHANGEMENT DES CONDITIONS LIMITES UTILISATEUR
!===============================================================================

!===============================================================================
! 6.1  LA VARIABLE TRANSPORTEE EST LA TEMPERATURE
!===============================================================================

if (itherm.eq.1) then

  do ifac = 1, nfabor

    if (isothm(ifac).eq.itpimp .or.                             &
        isothm(ifac).eq.ipgrno .or.                             &
        isothm(ifac).eq.ifgrno) then
      rcodcl(ifac,ivart,1) = tparo(ifac)+xmtk
      rcodcl(ifac,ivart,2) = rinfin
      rcodcl(ifac,ivart,3) = 0.d0

    elseif (isothm(ifac).eq.iprefl) then
      rcodcl(ifac,ivart,1) = text(ifac)+xmtk
      rcodcl(ifac,ivart,2) = bxlam(ifac)/        &
                              bepa(ifac)
      rcodcl(ifac,ivart,3) = 0.d0

    elseif (isothm(ifac).eq.ifrefl) then
      icodcl(ifac,ivart) = 3
      rcodcl(ifac,ivart,1) = 0.d0
      rcodcl(ifac,ivart,2) = rinfin
    endif

  enddo

!===============================================================================
! 6.2  LA VARIABLE TRANSPORTEE EST L'ENTHALPIE
!===============================================================================

elseif (itherm.eq.2) then

  !---> LECTURES DES DONNEES UTILISATEURS
  !     ON CONVERTIT TPAROI EN ENTHALPIE DE BORD, STOCKEE DANS FLUNET,
  !     QUI EST UTILISE COMME AUXILIAIRE

  mode = 0

  do ifac = 1, nfabor
    if (isothm(ifac).eq.itpimp.or.                              &
        isothm(ifac).eq.ipgrno.or.                              &
        isothm(ifac).eq.ifgrno) then
      mode = -1
    endif
  enddo

  if (mode.eq.-1) then

    nlst = 0
    do ifac = 1, nfabor
      if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
        nlst = nlst+1
        lstfac(nlst) = ifac
      endif
    enddo

    call b_t_to_h(nlst, lstfac, tparo, bfnet)

  endif

  mode = 0

  do ifac = 1, nfabor
    if (isothm(ifac).eq.iprefl) then
      mode = -1
    endif
  enddo

  if (mode.eq.-1) then

    nlst = 0
    do ifac = 1, nfabor
      if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
        nlst = nlst+1
        lstfac(nlst) = ifac
      endif
    enddo

    call b_t_to_h(nlst, lstfac, text, thwall)

  endif

  do ifac = 1, nfabor

    if (isothm(ifac).eq.itpimp.or.                              &
        isothm(ifac).eq.ipgrno.or.                              &
        isothm(ifac).eq.ifgrno) then
      rcodcl(ifac,ivart,1) = bfnet(ifac)
      rcodcl(ifac,ivart,2) = rinfin
      rcodcl(ifac,ivart,3) = 0.d0

      if (ihgas.ge.0) then
        rcodcl(ifac,ivahg,1) = bfnet(ifac)
        rcodcl(ifac,ivahg,2) = rinfin
        rcodcl(ifac,ivahg,3) = 0.d0
      endif
    elseif (isothm(ifac).eq.iprefl) then

      rcodcl(ifac,ivart,1) = thwall(ifac)
      ! hext
      rcodcl(ifac,ivart,2) =  bxlam(ifac) / (bepa(ifac))
      rcodcl(ifac,ivart,3) = 0.d0

      if (ihgas.ge.0) then
        rcodcl(ifac,ivahg,1) = thwall(ifac)
        rcodcl(ifac,ivahg,2) = bxlam(ifac) / bepa(ifac)
        rcodcl(ifac,ivahg,3) = 0.d0
      endif
    elseif (isothm(ifac).eq.ifrefl) then
      icodcl(ifac,ivart) = 3
      rcodcl(ifac,ivart,1) = 0.d0
      rcodcl(ifac,ivart,2) = rinfin

      if (ihgas.ge.0) then
        icodcl(ifac,ivahg) = 3
        rcodcl(ifac,ivahg,1) = 0.d0
        rcodcl(ifac,ivahg,2) = rinfin
      endif
    endif

  enddo

endif

! Update boundary temperature field

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
    if (itpscl.eq.2) then
      b_temp(ifac) = tparo(ifac) - tkelvi
    else
      b_temp(ifac) = tparo(ifac)
    endif
  endif
enddo

! Free memory
deallocate(isothm,lstfac)
deallocate(tempk,thwall)
deallocate(text, tint)
deallocate(tparo)

!--------
! Formats
!--------

 1000 format (/, &
 3X,'** INFORMATIONS SUR LE MODULE DE RAYONNEMENT'             ,/,&
3X,'   ------------------------------------------'             ,/,&
3X,' Initialisation de la temperature de paroi'                ,/,&
3X,' (TPAROI) avec le profil utilisateur (TINTP)'              ,/,&
3X,' et du flux incident aux parois (QINCID).'                 ,/)

#if defined(_CS_LANG_FR)

 2000 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : RAYONNEMENT'                                 ,/,&
'@    ========='                                               ,/,&
'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES' ,/,&
'@'                                                            ,/,&
'@  Le numero de zone associee a une face doit etre'           ,/,&
'@    un entier strictement positif et inferieur ou egal a'    ,/,&
'@    nozrdm = ',i10                                           ,/,&
'@  Ce numero (IZFRDP(IFAC)) est hors de ces bornes'           ,/,&
'@    pour ',i10, ' faces'                                     ,/,&
'@    derniere face :'                                         ,/,&
'@      zone         ', i10                                    ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@'                                                            ,/,&
'@  Verifier les conditions aux limites dans'                  ,/,&
'@   cs_user_radiative_transfer_bcs.f90.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : CONDITIONS AUX LIMITES DE RAYONNEMENT'       ,/,&
'@    ========='                                               ,/,&
'@'                                                            ,/,&
'@  Le nombre maximal de zones frontieres qui peuvent etre'    ,/,&
'@    definies par l''utilisateur est NBZRDM = ',i10           ,/,&
'@    Il a ete depasse'                                        ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute'                            ,/,&
'@'                                                            ,/,&
'@  Verifier les conditions aux limites de rayonnement'        ,/,&
'@'                                                            ,/,&
'@  Les NBZRDM premieres zones frontieres'                     ,/,&
'@    portent ici les numeros suivants :'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2002 format(i10)
 2110 format(                                                     &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    ISOTHP DOIT ETRE RENSEIGNE SUR TOUTES LES FACES DE PAROI',/,&
'@                                                            ',/,&
'@  Il ne l''a pas ete pour ', i10, ' faces'                   ,/,&
'@    derniere face :'                                         ,/,&
'@      zone         ', i10                                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
2111 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    ISOTHP A ETE RENSEIGNE SUR UNE FACE NON PAROI           ',/,&
'@                                                            ',/,&
'@  Il l''a ete pour ', i10, ' faces'                         ,/,&
'@    derniere face :'                                         ,/,&
'@      zone         ', i10                                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2120 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@    Avec ISOTHP = ITPIMP,'                                   ,/,&
'@      EPSP  doit etre un reel inclus dans [0.; 1.]'          ,/,&
'@      TINTP doit etre un reel strictement positif'           ,/,&
'@'                                                            ,/,&
'@  Ceci n''est pas le cas pour ', i10, ' faces'               ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2130 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@    Avec ISOTHP = IPGRNO,'                                   ,/,&
'@      EPSP  doit etre un reel inclus dans [0.; 1.]'          ,/,&
'@      XLAMP, EPAP, TEXTP, et TINTP doivent etre des reels'   ,/,&
'@                                   strictement positifs'     ,/,&
'@'                                                            ,/,&
'@  Ceci n''est pas le cas pour ', i10, ' faces'               ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@      XLAMP =      ', e12.4                                  ,/,&
'@      EPAP  =      ', e12.4                                  ,/,&
'@      TEXTP =      ', e12.4                                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2140 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@    Avec ISOTHP = IPREFL,'                                   ,/,&
'@      XLAMP, EPAP, TEXTP, et TINTP doivent etre des reels'   ,/,&
'@                                   strictement positifs'     ,/,&
'@'                                                            ,/,&
'@  Ceci n''est pas le cas pour ', i10, ' faces'               ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      XLAMP =      ', e12.4                                  ,/,&
'@      EPAP  =      ', e12.4                                  ,/,&
'@      TEXTP =      ', e12.4                                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2150 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@    Avec ISOTHP = IFGRNO,'                                   ,/,&
'@      EPSP  doit etre un reel inclus dans [0.; 1.]'          ,/,&
'@      TINTP doit etre un reel strictement positif'           ,/,&
'@'                                                            ,/,&
'@  Ceci n''est pas le cas pour ', i10, ' faces'               ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2160 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@    Avec ISOTHP = IFREFL,'                                   ,/,&
'@      TINTP doit etre un reel strictement positif'           ,/,&
'@'                                                            ,/,&
'@  Ceci n''est pas le cas pour ', i10, ' faces'               ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2170 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@    Avec ISOTHP = ITPT1D,'                                   ,/,&
'@      EPSP  doit etre un reel inclus dans [0.; 1.]'          ,/,&
'@      TINTP doit etre un reel strictement positif'           ,/,&
'@'                                                            ,/,&
'@  Ceci n''est pas le cas pour ', i10, ' faces'               ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2180 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@  Valeur interdite de ISOTHM pour ', i10, ' faces'           ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      ISOTHM =     ', i10                                    ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2220 format(                                                     &
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@    Avec ISOTHP = ITPIMP,'                                   ,/,&
'@    XLAMP, EPAP et TEXTP ne doivent pas etre renseignes'     ,/,&
'@                                                            ',/,&
'@  Ceci n''est pas le cas pour ', i10, ' faces'               ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      XLAMP =      ', e12.4                                  ,/,&
'@      EPAP  =      ', e12.4                                  ,/,&
'@      TEXTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2240 format(                                                     &
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@    Avec ISOTHP = IPREFL,'                                   ,/,&
'@    EPSP ne doit pas etre renseigne'                         ,/,&
'@                                                            ',/,&
'@  Ceci n''est pas le cas pour ', i10, ' faces'               ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2250 format(                                                     &
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@    Avec ISOTHP = IFGRNO,'                                   ,/,&
'@    XLAMP, EPAP et TEXTP ne doivent pas etre renseignes'     ,/,&
'@                                                            ',/,&
'@  Ceci n''est pas le cas pour ', i10, ' faces'               ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      XLAMP =      ', e12.4                                  ,/,&
'@      EPAP  =      ', e12.4                                  ,/,&
'@      TEXTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2260 format(                                                     &
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT'   ,/,&
'@    ========='                                               ,/,&
'@    Avec ISOTHP = IFREFL,'                                   ,/,&
'@    EPSP, XLAMP, EPAP et TEXTP ne doivent pas etre renseignes',/,&
'@                                                            ',/,&
'@  Ceci n''est pas le cas pour ', i10, ' faces'               ,/,&
'@    derniere face avec erreur (zone = :',i10,')'             ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@      XLAMP =      ', e12.4                                  ,/,&
'@      EPAP  =      ', e12.4                                  ,/,&
'@      TEXTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites de rayonnement.'       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#else

 2000 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@'                                                            ,/,&
'@  The zone number associated with a face must be'            ,/,&
'@    a strictly positive integer equal to or less than'       ,/,&
'@    nozrdm = ',i10                                           ,/,&
'@  This number (IZFRDP(IFAC)) is out of these bounds'         ,/,&
'@    for ',i10, ' faces'                                      ,/,&
'@    last face:'                                              ,/,&
'@      zone         ', i10                                    ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@'                                                            ,/,&
'@  The maximum number of boundary zones which may be'         ,/,&
'@    defined by the user ist NBZRDM = ',i10                   ,/,&
'@    It has been exceeded.                                   ',/,&
'@                                                            ',/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@  The nbzrdm first boundary zones'                           ,/,&
'@    have the following numbers:'                             ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2002 format(i10)
 2110 format(                                                     &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    ISOTHP MUST BE DEFINED ON ALL WALL FACES'                ,/,&
'@'                                                            ,/,&
'@  It was not defined for ', i10, ' faces'                    ,/,&
'@    last face:'                                              ,/,&
'@      zone         ', i10                                    ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run.'                            ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2111 format(                                                     &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    ISOTHP HAS BEEN DEFINED ON NON-WALL FACES'               ,/,&
'@'                                                            ,/,&
'@  It was defined for ', i10, ' faces'                        ,/,&
'@    last face:'                                              ,/,&
'@      zone         ', i10                                    ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2120 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    With ISOTHP = ITPIMP,'                                   ,/,&
'@      EPSP  must be a real in the range [0.; 1.]'            ,/,&
'@      TINTP must be a strictly positive real'                ,/,&
'@'                                                            ,/,&
'@  This is not the case for ', i10, ' faces'                  ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2130 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    With ISOTHP = ITGRNO,'                                   ,/,&
'@      EPSP  must be a real in the range [0.; 1.]'            ,/,&
'@      XLAMP, EPAP, TINTP, TEXTP must be strictly'            ,/,&
'@                                        positive reals'      ,/,&
'@'                                                            ,/,&
'@  This is not the case for ', i10, ' faces'                  ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@      XLAMP =      ', e12.4                                  ,/,&
'@      EPAP  =      ', e12.4                                  ,/,&
'@      TEXTP =      ', e12.4                                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2140 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    With ISOTHP = IPREFL,'                                   ,/,&
'@      XLAMP, EPAP, TINTP, TEXTP must be strictly'            ,/,&
'@                                        positive reals'      ,/,&
'@'                                                            ,/,&
'@  This is not the case for ', i10, ' faces'                  ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      XLAMP =      ', e12.4                                  ,/,&
'@      EPAP  =      ', e12.4                                  ,/,&
'@      TEXTP =      ', e12.4                                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2150 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    With ISOTHP = IFGRNO,'                                   ,/,&
'@      EPSP  must be a real in the range [0.; 1.]'            ,/,&
'@      TINTP must be a strictly positive real'                ,/,&
'@'                                                            ,/,&
'@  This is not the case for ', i10, ' faces'                  ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2160 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    With ISOTHP = IFREFL,'                                   ,/,&
'@      TINTP must be a strictly positive real'                ,/,&
'@'                                                            ,/,&
'@  This is not the case for ', i10, ' faces'                  ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2170 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    With ISOTHP = ITPT1D,'                                   ,/,&
'@      EPSP  must be a real in the range [0.; 1.]'            ,/,&
'@      TINTP must be a strictly positive real'                ,/,&
'@'                                                            ,/,&
'@  This is not the case for ', i10, ' faces'                  ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@      TINTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2180 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@  Forbidden value of ISOTHM for ', i10, ' faces'             ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      ISOTHM =     ', i10                                    ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2220 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    With ISOTHP = ITPIMP,'                                   ,/,&
'@    XLAMP, EPAP and TEXTP must not be defined'               ,/,&
'@'                                                            ,/,&
'@  This is not the case for ', i10, ' faces'                  ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      XLAMP =      ', e12.4                                  ,/,&
'@      EPAP  =      ', e12.4                                  ,/,&
'@      TEXTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2240 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    With ISOTHP = IPREFL,'                                   ,/,&
'@    EPSP must not be defined'                                ,/,&
'@'                                                            ,/,&
'@  This is not the case for ', i10, ' faces'                  ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2250 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    With ISOTHP = IFGRNO,'                                   ,/,&
'@    XLAMP, EPAP and TEXTP must not be defined'               ,/,&
'@'                                                            ,/,&
'@  This is not the case for ', i10, ' faces'                  ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      XLAMP =      ', e12.4                                  ,/,&
'@      EPAP  =      ', e12.4                                  ,/,&
'@      TEXTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2260 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN RADIATIVE BOUNDARY CONDITIONS CHECK'   ,/,&
'@    ======='                                                 ,/,&
'@    With ISOTHP = IFREFL,'                                   ,/,&
'@    EPSP, XLAMP, EPAP and TEXTP must not be defined'         ,/,&
'@'                                                            ,/,&
'@  This is not the case for ', i10, ' faces'                  ,/,&
'@    last face with error (zone = :',i10,')'                  ,/,&
'@      EPSP  =      ', e12.4                                  ,/,&
'@      XLAMP =      ', e12.4                                  ,/,&
'@      EPAP  =      ', e12.4                                  ,/,&
'@      TEXTP =      ', e12.4                                  ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be run'                             ,/,&
'@'                                                            ,/,&
'@  Check radiative boundary conditions.'                      ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

!----
! End
!----

return

end subroutine raycli

!===============================================================================
! Local functions
!===============================================================================

!> \brief synchronize radiative boundary condition error logging across
!         MPI ranks.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in, out] nerloc       number of errors (local rank in, global out)
!> \param[in]      nerrcd       number of codes saved at error faces
!> \param[in, out] znferr       zone number for one error face (local in,
!                               broadcast out)
!> \param[in, out] rvferr       values saved at one error face (local in,
!                               broadcast out)
!_______________________________________________________________________________

subroutine sync_rad_bc_err &
 ( nerloc , nerrcd , znferr, rvferr )

!===============================================================================
! Module files
!===============================================================================

use parall

!===============================================================================

implicit none

! Arguments

integer nerloc, nerrcd
integer znferr
double precision rvferr(nerrcd)

! Local variables

integer irkerr

!===============================================================================

if (irangp.ge.0) then
  irkerr = -1
  if (nerloc.gt.0) irkerr = irangp
  call parcpt(nerloc)
  if (nerloc .ne. 0) then
    call parcmx(irkerr)
    call parall_bcast_i(irkerr, znferr)
    call parbcr(irkerr, nerrcd, rvferr)
  endif
endif

return
end subroutine sync_rad_bc_err
