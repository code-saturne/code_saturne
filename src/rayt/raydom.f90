!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine raydom &
!================

 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   izfrad ,                                                       &
   dt     , rtp    , rtpa   , propce )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  Enveloppe principale du module de resolution de l'equation
!  des transferts radiatifs

!  Deux methodes sont disponibles :

!    1) La methode : "Discretes Ordinates Methods" (DOM)
!    2) L'approximation P-1 (recommande uniquement pour le CP)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfrad(nfabor    ! te ! <-- ! numero de zone des faces de bord               !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use parall
use period
use ppppar
use ppthch
use cs_fuel_incl
use ppincl
use cpincl
use radiat
use ihmpre
use dimens, only: ndimfb
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          itypfb(ndimfb)
integer          izfrad(ndimfb)

double precision dt(ncelet), rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)

! Local variables

integer          iappel
integer          ifac   , iel    , iok    , izone
integer          inc    , iccocg , iwarnp , imligp , nswrgp
integer          mode   , icla   , ipcla  , ivar0
integer          ivart
integer          idverl
integer          iflux(nozrdm)
double precision epsrgp, climgp, extrap
double precision aa, bb, ckmin, unspi, xlimit, flunmn
double precision flux(nozrdm)
double precision vv, sf, xlc, xkmin, pp

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbrs, rovsdt
double precision, allocatable, dimension(:) :: dcp
double precision, allocatable, dimension(:) :: ckmel
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:,:) :: tempk
double precision, allocatable, dimension(:) :: coefap, coefbp
double precision, allocatable, dimension(:) :: cofafp, cofbfp
double precision, allocatable, dimension(:) :: flurds, flurdb

double precision, dimension(:), pointer :: tparo, bqinci
double precision, dimension(:), pointer :: bxlam, bepa, beps, bfnet

integer    ipadom
data       ipadom /0/
save       ipadom

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

! Allocate temporary arrays for the radiative equations resolution
allocate(viscf(nfac), viscb(ndimfb))
allocate(smbrs(ncelet), rovsdt(ncelet))

! Allocate specific arrays for the radiative transfert module
allocate(tempk(ncelet,nrphas))
allocate(coefap(ndimfb), coefbp(ndimfb))
allocate(cofafp(ndimfb), cofbfp(ndimfb))
allocate(flurds(nfac), flurdb(ndimfb))

! Allocate work arrays
allocate(ckmel(ncelet))

! Mapa field arrays
call field_get_val_s(itparo, tparo)
call field_get_val_s(iqinci, bqinci)
call field_get_val_s(ixlam, bxlam)
call field_get_val_s(iepa, bepa)
call field_get_val_s(ieps, beps)
call field_get_val_s(ifnet, bfnet)

!===============================================================================
! 1. INITIALISATIONS GENERALES
!===============================================================================

!---> Number of passes
ipadom = ipadom + 1
if (ipadom.gt.1 .and. mod(ntcabs,nfreqr).ne.0) return

write(nfecra,1000)

!---> Constants initialization
unspi = 1.d0/pi

!---> Index of thermal variable
ivart = isca(iscalt)

!=============================================================================
! 3.1 Absorption coefficient of environment semitransparent
!=============================================================================

!--> Initialization to a non-admissible value for testing after usray3
do iel = 1, ncel
  propce(iel,ipproc(icak(1))) = -grand
enddo

!--> Absorption coefficient for different modules

! Warning: for the approximation P-1, the absorption coefficient is required
!          for boundary conditions

if (ippmod(iphpar).ge.2) then

  call ppcabs(rtp, propce)
  !==========

  !---> ckmel stores temporarly the absorbption coefficient
  !     of gaz-particle mixing

  if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0  ) then

    do iel = 1, ncel
      ckmel(iel) = propce(iel,ipproc(icak(1)))
    enddo

    if (ippmod(iccoal).ge.0 ) then
      do icla = 1,nclacp
        ipcla = 1+icla
        do iel = 1, ncel
          ckmel(iel) = ckmel(iel)                                   &
                     + ( propce(iel,ipproc(ix2(icla)))              &
                       * propce(iel,ipproc(icak(ipcla))) )
        enddo
      enddo
    else
      do icla = 1,nclafu
        ipcla = 1+icla
        do iel = 1, ncel
          ckmel(iel) = ckmel(iel)                                   &
                     + ( rtpa(iel,isca(iyfol(icla)))                &
                       * propce(iel,ipproc(icak(ipcla))) )
        enddo
      enddo
    endif

    do iel = 1, ncel
      propce(iel,ipproc(icak(1))) = ckmel(iel)
    enddo
  endif

else


  !---> Reading of User datas

  !   - Interface Code_Saturne
  !     ======================

  if (iihmpr.eq.1) then

    call uiray3(propce(1,ipproc(icak(1))), ncel, imodak)

    if (iirayo.eq.2 .and. ippmod(iphpar).le.1 .and. ipadom.le.3) then
      sf = 0.d0
      vv = 0.d0

      ! Compute the caracteristic length of the computational domain
      do ifac = 1, nfabor
        sf = sf + sqrt(surfbo(1,ifac)**2 +                      &
                       surfbo(2,ifac)**2 +                      &
                       surfbo(3,ifac)**2 )
      enddo
      if (irangp.ge.0) then
        call parsom(sf)
      endif

      do iel = 1, ncel
        vv = vv + volume(iel)
      enddo
      if (irangp.ge.0) then
        call parsom(vv)
      endif

      xlc = 3.6d0 * vv / sf

      !  Clipping on ck
      xkmin = 1.d0 / xlc

      iok = 0
      do iel = 1, ncel
        if (propce(iel,ipproc(icak(1))).lt.xkmin) then
          iok = iok +1
        endif
      enddo

      ! Warning if the optical thickness is too big
      pp = xnp1mx/100.0d0
      if (dble(iok).gt.pp*dble(ncel)) then
        write(nfecra,6000) xkmin, dble(iok)/dble(ncel)*100.d0,  &
                           xnp1mx
      endif
    endif

  endif

  call usray3 &
  !==========
( nvar   , nscal  , iappel ,                                     &
  itypfb ,                                                       &
  izfrad ,                                                       &
  dt     ,                                                       &
  propce(1,ipproc(icak(1))))

endif

!--> General checking

!---> P-1: check that ck is strictly greater than 0
if (iirayo.eq.2) then

  ckmin = propce(1,ipproc(icak(1)))
  do iel = 1, ncel
    ckmin = min(ckmin,propce(iel,ipproc(icak(1))))
  enddo
  if (ckmin.lt.0.d0) then
    write(nfecra,2020)
    call csexit (1)
  endif

!---> Dom:  check that ck is greater than 0
else if (iirayo.eq.1) then

  ckmin = propce(1,ipproc(icak(1)))
  do iel = 1, ncel
    ckmin = min(ckmin,propce(iel,ipproc(icak(1))))
  enddo
  if (ckmin.lt.0.d0) then
    write(nfecra,2010) ckmin
    call csexit (1)
  endif

endif

!---> Check of a transparent case
idverl = idiver

aa = zero
do iel = 1, ncel
  aa = aa + propce(iel,ipproc(icak(1)))
enddo
if (irangp.ge.0) then
  call parmax(aa)
endif
if (aa.le.epzero) then
  write(nfecra,1100)
  idverl = -1
endif

!=============================================================================
! 4. Temperature storing (in Kelvin) in tempk(iel, irphas)
!=============================================================================

if (idverl.ge.0) then

  !---> Temperature transport
  if (itherm.eq.1) then

    ! iscalt is in Celsius
    if (itpscl.eq.2) then
      do iel = 1, ncel
        tempk(iel,1) = rtpa(iel,ivart) + tkelvi
      enddo
    else
      do iel = 1, ncel
        tempk(iel,1) = rtpa(iel,ivart)
      enddo
    endif

  !---> Enthalpy transport (flurdb is a temporary array)
  else if (itherm.eq.2) then

    mode = 1

    if (ippmod(iphpar).le.1) then

      call usray4 &
      !==========
      ( nvar   , nscal  ,                                            &
        mode   ,                                                     &
        itypfb ,                                                     &
        dt     ,                                                     &
        tparo  , flurdb , tempk(1,1)  )

    else

      call ppray4 &
      !==========
    ( mode   ,                                                       &
      itypfb ,                                                       &
      rtp    , rtpa   , propce ,                                     &
      tparo  , flurdb , tempk(1,1)  )

    endif

    if (ippmod(iccoal).ge.0) then

      ! Particules' temperature
      do icla = 1, nclacp
        ipcla = 1+icla
        do iel = 1, ncel
          tempk(iel,ipcla) = propce(iel,ipproc(itemp2(icla)))
        enddo
      enddo

    ! Fuel
    else if (ippmod(icfuel).ge.0) then

      do icla = 1, nclafu
        ipcla = 1+icla
        do iel = 1, ncel
          tempk(iel,ipcla) = propce(iel,ipproc(itemp2(icla)))
        enddo
      enddo

    endif

  else
    write(nfecra,3500) itherm
    call csexit (1)
  endif

  !---> on se sert de propce(iel,ipproc(itsri(1))) comme un auxiliaire pour
  !       stocker stephn*ck*tempk**4 ici
  !       plus bas on justifiera le nom.

  if (ippmod(icod3p).eq.-1 .and. ippmod(icoebu).eq.-1) then

    ! Rayonnement standard, flamme CP ou fuel
    do iel = 1, ncel
      propce(iel,ipproc(itsri(1))) = stephn  *             &
       propce(iel,ipproc(icak(1)))*(tempk(iel,1)**4)
    enddo

  else

    ! Flamme de diffusion ou flamme de premelange
    do iel = 1, ncel
      propce(iel,ipproc(itsri(1))) = stephn  *             &
       propce(iel,ipproc(icak(1)))*propce(iel,ipproc(it4m))
    enddo

  endif

  ! Coal
  if (ippmod(iccoal).ge.0) then
    do icla = 1, nclacp
      ipcla = 1+icla
      do iel = 1, ncel
        propce(iel,ipproc(itsri(ipcla))) =  stephn  *           &
          propce(iel,ipproc(icak(ipcla)))*(tempk(iel,ipcla)**4)
      enddo
    enddo

  ! Fuel
  else if (ippmod(icfuel).ge.0) then
    do icla = 1, nclafu
      ipcla = 1+icla
      do iel = 1, ncel
        propce(iel,ipproc(itsri(ipcla))) =  stephn  *           &
          propce(iel,ipproc(icak(ipcla)))*(tempk(iel,ipcla)**4)
      enddo
    enddo
  endif

else
  do iel = 1, ncel
    propce(iel,ipproc(itsri(1))) = zero
  enddo
endif

!===============================================================================
! 5.1 Radiative P-1 model
!===============================================================================

if (iirayo.eq.2) then

  !--> Terme source explicite de l'equation sur Theta4

  do iel = 1, ncel
    smbrs(iel) = 3.d0 * propce(iel,ipproc(icak(1))) *      &
       ( tempk(iel,1) ** 4) * volume(iel)
  enddo

  ! Tenir compte de l'absorption des particules

  ! Coal
  if (ippmod(iccoal).ge.0) then
    do icla = 1, nclacp
      ipcla = 1+icla
      do iel = 1,ncel
        smbrs(iel) = smbrs(iel)                               &
                   + (3.d0*propce(iel,ipproc(ix2(icla)))      &
                     * propce(iel,ipproc(icak(ipcla)))        &
                     * (tempk(iel,ipcla)**4) * volume(iel) )
      enddo
    enddo

  ! Fuel
  else if (ippmod(icfuel).ge.0) then
    do icla = 1, nclafu
      ipcla = 1+icla
      do iel = 1,ncel
        smbrs(iel) = smbrs(iel)                               &
                   + (3.d0*rtpa(iel,isca(iyfol(icla)))        &
                     * propce(iel,ipproc(icak(ipcla)))        &
                     * (tempk(iel,ipcla)**4) * volume(iel) )
      enddo
    enddo
  endif

  !--> Terme source implicite de l'equation sur Theta4
  do iel = 1, ncel
    rovsdt(iel) =  3.d0*propce(iel,ipproc(icak(1)))*volume(iel)
  enddo

  ! Tenir compte de l'absorption des particules

  ! Coal
  if (ippmod(iccoal).ge.0) then
    do icla = 1, nclacp
      ipcla = 1+icla
      do iel = 1,ncel
        rovsdt(iel) = rovsdt(iel)                                      &
                    + (3.d0*propce(iel,ipproc(ix2(icla)))              &
                      * propce(iel,ipproc(icak(ipcla))) * volume(iel) )
      enddo
    enddo

  ! Fuel
  else if (ippmod(icfuel).ge.0) then
    do icla = 1, nclafu
      ipcla = 1+icla
      do iel = 1,ncel
        rovsdt(iel) = rovsdt(iel)                                      &
                    + (3.d0*rtpa(iel,isca(iyfol(icla)))                &
                      * propce(iel,ipproc(icak(ipcla))) * volume(iel) )
      enddo
    enddo

  endif

  !--> Inverse du coefficient de diffusion de l'equation sur Theta4
  !       A priori ckmel contient deja la bonne info, mais pour plus de
  !       securite  on le re-remplit

  do iel = 1, ncel
    ckmel(iel) = propce(iel,ipproc(icak(1)))
  enddo

  ! Tenir compte de l'absorption des particules

  ! Coal
  if (ippmod(iccoal).ge.0) then
    do icla = 1, nclacp
      ipcla = 1+icla
      do iel = 1,ncel
        ckmel(iel) = ckmel(iel)                                  &
                   + ( propce(iel,ipproc(ix2(icla)))             &
                     * propce(iel,ipproc(icak(ipcla))) )
      enddo
    enddo

  ! Fuel
  else if (ippmod(icfuel).ge.0) then
    do icla = 1, nclafu
      ipcla = 1+icla
      do iel = 1, ncel
        ckmel(iel) = ckmel(iel)                                  &
                   + ( rtpa(iel,isca(iyfol(icla)))               &
                     * propce(iel,ipproc(icak(ipcla))) )
      enddo
    enddo
  endif

  ! Update Boundary condiction coefficients

  call raycll &
  !==========
  ( itypfb ,                                                       &
    izfrad ,                                                       &
    coefap , coefbp ,                                              &
    cofafp , cofbfp ,                                              &
    tparo  , bqinci , beps   ,                                     &
    propce(1,ipproc(icak(1))), ckmel )

  ! Solving

  call raypun &
  !==========
( itypfb ,                                                       &
  coefap , coefbp ,                                              &
  cofafp , cofbfp ,                                              &
  flurds , flurdb ,                                              &
  viscf  , viscb  ,                                              &
  smbrs  , rovsdt ,                                              &
  propce(1,ipproc(iabs(1))),propce(1,ipproc(iemi(1))),           &
  propce(1,ipproc(itsre(1))) , propce(1,ipproc(iqx))  ,          &
  propce(1,ipproc(iqy))   , propce(1,ipproc(iqz))  ,             &
  bqinci, beps , tparo  ,                                        &
  ckmel    )

!===============================================================================
! 5.2 Solving of the radiative transfert equation
!===============================================================================

else if (iirayo.eq.1) then

  !--> Terme source explicite de l'equation sur la luminance
  do iel = 1, ncel
    smbrs(iel) = propce(iel,ipproc(itsri(1)))*volume(iel)*unspi
  enddo

  ! Coal
  if (ippmod(iccoal).ge.0) then
    do icla = 1,nclacp
      ipcla = 1+icla
      do iel = 1,ncel
        smbrs(iel) = smbrs(iel)                                 &
                + propce(iel,ipproc(ix2(icla)))                 &
                 *propce(iel,ipproc(itsri(ipcla)))*volume(iel)  &
                 *unspi
      enddo
    enddo

  ! Fuel
  elseif (ippmod(icfuel).ge.0) then
    do icla = 1,nclafu
      ipcla = 1+icla
      do iel = 1,ncel
        smbrs(iel) = smbrs(iel)                                 &
                    + rtpa(iel,isca(iyfol(icla)))               &
                 *propce(iel,ipproc(itsri(ipcla)))*volume(iel)  &
                 *unspi
      enddo
    enddo

  endif

  !--> Terme source implicite de l'equation sur la luminance
  !      KL + div(LS) = KL0 integre sur le volume de controle
  do iel = 1, ncel
    rovsdt(iel) = propce(iel,ipproc(icak(1))) * volume(iel)
  enddo

  ! Coal
  if (ippmod(iccoal).ge.0) then
    do icla = 1,nclacp
      ipcla = 1+icla
      do iel = 1,ncel
        rovsdt(iel) = rovsdt(iel)                                      &
                    + propce(iel,ipproc(ix2(icla)))                    &
                      * propce(iel,ipproc(icak(ipcla))) * volume(iel)
      enddo
    enddo

  ! Fuel
  elseif (ippmod(icfuel).ge.0) then
    do icla = 1,nclafu
      ipcla = 1+icla
      do iel = 1,ncel
        rovsdt(iel) = rovsdt(iel)                                      &
                    + rtpa(iel,isca(iyfol(icla)))                      &
                      * propce(iel,ipproc(icak(ipcla))) * volume(iel)
      enddo
    enddo

  endif

  ! Update Boundary condiction coefficients

  call raycll &
  !==========
  ( itypfb ,                                                       &
    izfrad ,                                                       &
    coefap , coefbp ,                                              &
    cofafp , cofbfp ,                                              &
    tparo  , bqinci , beps   ,                                     &
    propce(1,ipproc(icak(1))), ckmel )

  ! Solving

  call raysol &
  !==========
 ( coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
   flurds , flurdb ,                                              &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt ,                                              &
   propce(1,ipproc(iabs(1))),propce(1,ipproc(iemi(1))) ,          &
   propce(1,ipproc(itsre(1))),propce(1,ipproc(iqx))         ,     &
   propce(1,ipproc(iqy))    , propce(1,ipproc(iqz))   ,           &
   bqinci , bfnet )

endif

!===============================================================================
! 5.3 Storing of the integral of the luminance for the Lagrangian module
!===============================================================================

!  Si dans le module lagrangien on resout une equation de la temperature
!    sur les particules (iphyla=1 et itpvar=1) ou si les particules
!    sont des grains de charbon (iphyla=2), on a besoin de
!                                     /    ->  ->
!    l'integrale de la luminance SA= /  L( X , S ). DOMEGA
!                                   /4.PI
!  On stocke cette variable quelque soit le choix des options

do iel = 1,ncel
  propce(iel,ipproc(ilumin)) = propce(iel,ipproc(itsre(1)))
enddo

!===============================================================================
! 6. Net radiative flux at walls: compuation and integration
!===============================================================================

!--> Initialization to a non-admissible value for testing after usray5
do ifac = 1,nfabor
  bfnet(ifac) = -grand
enddo

!---> Reading of User datas
call usray5 &
!==========
( nvar   , nscal  ,                                              &
  itypfb ,                                                       &
  izfrad ,                                                       &
  dt     ,                                                       &
  coefap , coefbp ,                                              &
  cofafp , cofbfp ,                                              &
  tparo  , bqinci ,                                              &
  bfnet  , bxlam  , bepa   , beps   ,                            &
  propce(1,ipproc(icak(1)))  )

!---> Check flunet
iok = 0
xlimit = -grand*0.1d0
flunmn = grand

do ifac = 1,nfabor
  if (bfnet(ifac).le.xlimit) then
    iok = iok + 1
    flunmn = min(flunmn,bfnet(ifac))
    write(nfecra,4000)ifac,izfrad(ifac),itypfb(ifac)
  endif
enddo

if (iok.ne.0) then
  write(nfecra,4100) flunmn
  call csexit (1)
endif

!--> Integration du flux net sur les differentes zones de frontieres
!     IFLUX sert en parallele pour reperer les zones existantes

do izone = 1, nozrdm
  flux(izone) = 0.d0
  iflux(izone) = 0
enddo
do ifac = 1,nfabor
  izone = izfrad(ifac)
  flux(izone) = flux(izone) + bfnet(ifac)*surfbn(ifac)
  iflux(izone) = 1
enddo
if(irangp.ge.0) then
  call parrsm(nozarm,flux )
  call parimx(nozarm,iflux)
endif


write(nfecra,5000)
write(nfecra,5010)
do izone = 1, nozarm
  if(iflux(izone).eq.1) then
    write(nfecra,5020) izone,flux(izone)
  endif
enddo
write(nfecra,5000)


!--> Integration de la densite de flux net aux frontieres

aa = zero
do ifac = 1,nfabor
  aa =  aa + bfnet(ifac) * surfbn(ifac)
enddo
if(irangp.ge.0) then
  call parsom(aa)
endif
write(nfecra,5030) aa

!===============================================================================
! 7. Implicit and explicit radiative source terms
!===============================================================================


!===============================================================================
! 7.1 Semi-analitical radiative source termes
!===============================================================================

if (idverl.ge.0) then

  do iel = 1, ncel

    !--> part d'absorption du terme source explicite
    propce(iel,ipproc(iabs(1))) = propce(iel,ipproc(icak(1))) &
                                * propce(iel,ipproc(itsre(1)))

    !--> part d'emission du terme source explicite
    propce(iel,ipproc(iemi(1))) = -4.d0*propce(iel,ipproc(itsri(1)))

  enddo

  ! Combustion CP : On rajoute la contribution des particules
  if (ippmod(iccoal).ge.0) then
    do icla = 1, nclacp
      ipcla = 1+icla
      do iel = 1, ncel
        ! Fluid
        propce(iel,ipproc(iabs(1))) = propce(iel,ipproc(iabs(1)))             &
                                    + propce(iel,ipproc(ix2(icla)))           &
                                      * propce(iel,ipproc(icak(ipcla)))       &
                                      * propce(iel,ipproc(itsre(1)))

        propce(iel,ipproc(iemi(1))) = propce(iel,ipproc(iemi(1)))             &
                                    - 4.0d0*propce(iel,ipproc(ix2(icla)))     &
                                      * propce(iel,ipproc(itsri(ipcla)))
        ! Particles
        propce(iel,ipproc(iabs(ipcla))) = propce(iel,ipproc(icak(ipcla)))     &
                                        * propce(iel,ipproc(itsre(1)))
        propce(iel,ipproc(iemi(ipcla))) = - 4.0d0*propce(iel,ipproc(itsri(ipcla)))
        propce(iel,ipproc(itsre(ipcla))) = propce(iel,ipproc(iabs(ipcla)))    &
                                         + propce(iel,ipproc(iemi(ipcla)))
      enddo
    enddo

  ! Combustion Fuel : On rajoute la contribution des particules
  elseif (ippmod(icfuel).ge.0) then
    do icla = 1, nclafu
      ipcla = 1+icla
      do iel = 1, ncel
        ! Fluid
        propce(iel,ipproc(iabs(1))) = propce(iel,ipproc(iabs(1)))             &
                                    + rtpa(iel,isca(iyfol(icla)))             &
                                      * propce(iel,ipproc(icak(ipcla)))       &
                                      * propce(iel,ipproc(itsre(1)))

        propce(iel,ipproc(iemi(1))) = propce(iel,ipproc(iemi(1)))             &
                                    - 4.0d0*rtpa(iel,isca(iyfol(icla)))       &
                                      * propce(iel,ipproc(itsri(ipcla)))
        ! Particles
        propce(iel,ipproc(iabs(ipcla))) = propce(iel,ipproc(icak(ipcla)))     &
                                        * propce(iel,ipproc(itsre(1)))
        propce(iel,ipproc(iemi(ipcla))) = - 4.0d0*propce(iel,ipproc(itsri(ipcla)))
        propce(iel,ipproc(itsre(ipcla))) = propce(iel,ipproc(iabs(ipcla)))    &
                                         + propce(iel,ipproc(iemi(ipcla)))
      enddo
    enddo

  endif

  !--> Premiere methode pour le calcul du terme source explicite :
  !    il est calcule comme la somme des termes d'absorption et d'emission
  !    (il faudra multiplier ce terme par volume(iel) dans covofi->raysca)
  do iel = 1, ncel
    propce(iel,ipproc(itsre(1))) = propce(iel,ipproc(iabs(1)))                &
                                 + propce(iel,ipproc(iemi(1)))
  enddo

  ! Allocate work arrays
  allocate(dcp(ncelet))

  !--> 1/Cp is stored in dcp
  if (icp.gt.0) then
    do iel = 1,ncel
      dcp(iel) = 1.d0/propce(iel,ipproc(icp))
    enddo
  else
    do iel = 1,ncel
      dcp(iel) = 1.d0/cp0
    enddo
  endif

  !--> Terme source implicite,
  !    (il faudra multiplier ce terme par VOLUME(IEL) dans COVOFI->RAYSCA)
  if (ippmod(icod3p).eq.-1 .and. ippmod(icoebu).eq.-1) then
    ! Rayonnement standard, flamme CP ou fuel
    do iel = 1, ncel
      propce(iel,ipproc(itsri(1))) =                         &
       -16.d0*propce(iel,ipproc(icak(1))) *stephn *          &
         (tempk(iel,1)**3) * dcp(iel)
    enddo

  else

    ! Flamme de diffusion ou flamme de premelange
    do iel = 1, ncel
      propce(iel,ipproc(itsri(1))) =                         &
       -16.d0*stephn*propce(iel,ipproc(icak(1)))*            &
           propce(iel,ipproc(it3m)) * dcp(iel)
    enddo

  endif

  deallocate(dcp)

  ! Combustion CP : On rajoute la contribution des particules
  if (ippmod(iccoal).ge.0) then
    do icla = 1, nclacp
      ipcla = 1+icla
      do iel = 1, ncel
        propce(iel,ipproc(itsri(1))) = propce(iel,ipproc(itsri(1)))          &
                                     - 16.d0*propce(iel,ipproc(icak(ipcla))) &
                                       * propce(iel,ipproc(ix2(icla)))       &
                                       * stephn * (tempk(iel,ipcla)**3)      &
                                       / cp2ch(ichcor(icla))
        propce(iel,ipproc(itsri(ipcla))) = -16.d0                            &
                                         * propce(iel,ipproc(icak(ipcla)))   &
                                         * stephn*(tempk(iel,ipcla)**3)      &
                                         / cp2ch(ichcor(icla))
      enddo
    enddo

  ! Combustion FUEL : On rajoute la contribution des particules
  elseif (ippmod(icfuel).ge.0) then
    do icla = 1, nclafu
      ipcla = 1+icla
      do iel = 1, ncel
        propce(iel,ipproc(itsri(1))) = propce(iel,ipproc(itsri(1)))           &
                                     - 16.d0*propce(iel,ipproc(icak(ipcla)))  &
                                       * rtpa(iel,isca(iyfol(icla))) * stephn &
                                       * (tempk(iel,ipcla)**3) / cp2fol
        propce(iel,ipproc(itsri(ipcla))) = -16.d0                             &
                                         * propce(iel,ipproc(icak(ipcla)))    &
                                         * stephn * (tempk(iel,ipcla)**3)     &
                                         / cp2fol
      enddo
    enddo
  endif

else
  do iel = 1, ncel
    propce(iel,ipproc(iabs(1)))  = zero
    propce(iel,ipproc(iemi(1)))  = zero
    propce(iel,ipproc(itsre(1))) = zero
    propce(iel,ipproc(itsri(1))) = zero
  enddo
endif

!===============================================================================
! 7.2 Explicit conservative radiative source termes
!===============================================================================

! coefap and coefbp are NOW Boundary conditions on the divergence

if (idverl.eq.1 .or. idverl.eq.2) then

  ! Allocate a temporary array for gradient computation
  allocate(grad(ncelet,3))

  do ifac = 1,nfabor
    coefbp(ifac) = zero
  enddo

  !--> Calculation of the divergence

  ! En periodique et parallele, echange avant calcul du gradient
  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec &
    !==========
  ( propce(1,ipproc(iqx)), propce(1,ipproc(iqy)), propce(1,ipproc(iqz)) )
  endif

  ! Donnees pour le calcul de la divergence
  inc     = 1
  iccocg  = 1
  imligp  = -1
  iwarnp  = iimlum
  epsrgp  = 1.d-8
  climgp  = 1.5d0
  extrap  = 0.d0
  nswrgp  = 100

  !---> X direction
  do ifac = 1, nfabor
    coefap(ifac) = bfnet(ifac)*surfbo(1,ifac) / surfbn(ifac)
  enddo

  !  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
  !     n'est pas la vitesse ni Rij)
  !    sera a revoir pour la periodicite de rotation
  ivar0 = 0
  call grdcel &
  !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   propce(1,ipproc(iqx))    , coefap , coefbp ,                   &
   grad   )

  do iel = 1,ncel
    propce(iel,ipproc(itsre(1))) = - grad(iel,1)
  enddo

  !---> Y direction
  do ifac = 1, nfabor
    coefap(ifac) = bfnet(ifac)*surfbo(2,ifac) / surfbn(ifac)
  enddo

  !  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
  !     n'est pas la vitesse ni Rij)
  !    sera a revoir pour la periodicite de rotation
  ivar0 = 0
  call grdcel &
  !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   propce(1,ipproc(iqy))    , coefap , coefbp ,                   &
   grad   )

  do iel = 1, ncel
    propce(iel,ipproc(itsre(1))) = propce(iel,ipproc(itsre(1))) - grad(iel,2)
  enddo

  !---> Z direction
  do ifac = 1, nfabor
    coefap(ifac) = bfnet(ifac)*surfbo(3,ifac) / surfbn(ifac)
  enddo

  !  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
  !     n'est pas la vitesse ni Rij)
  !    sera a revoir pour la periodicite de rotation
  ivar0 = 0
  call grdcel &
  !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   propce(1,ipproc(iqz))    , coefap , coefbp ,                   &
   grad   )

  do iel = 1, ncel
    propce(iel,ipproc(itsre(1))) = propce(iel,ipproc(itsre(1))) - grad(iel,3)
  enddo

  ! Free memory
  deallocate(grad)

! Fin du calcul de la divergence
endif


!===============================================================================
! 7.3 Explicite radiative semi-analytical corrected source term
!===============================================================================


if (idverl.eq.2) then

  !---> comparaison des termes sources semi-analytique et conservatif
  aa = zero
  do iel = 1, ncel
    aa = aa + propce(iel,ipproc(itsre(1))) * volume(iel)
  enddo

  bb = zero
  do iel = 1,ncel
    bb = bb                                                                    &
       + (propce(iel,ipproc(iabs(1)))+propce(iel,ipproc(iemi(1))))*volume(iel)
  enddo

  if(irangp.ge.0) then
    call parsom(aa)
    call parsom(bb)
  endif

  aa = aa/bb

  !---> correction du terme source semi-analytique par le conservatif
  do iel = 1,ncel
    propce(iel,ipproc(itsre(1))) = ( propce(iel,ipproc(iabs(1)))              &
                                   + propce(iel,ipproc(iemi(1))))             &
                                 * aa
  enddo

endif

!===============================================================================
! 7.4 Finalization of explicit source terms
!===============================================================================

if (idverl.ge.0) then

  !--> Integration volumique du terme source explicite
  !    Le resultat de cette integration DOIT etre le meme que l'integration
  !    surfacique de la densite de flux net radiatif faite plus haut
  !    si  IDVERL = 1 ou 2

  aa = zero
  do iel = 1, ncel
    aa = aa + propce(iel,ipproc(itsre(1))) * volume(iel)
  enddo

  if(irangp.ge.0) then
    call parsom(aa)
  endif

  write(nfecra,5040) aa
  write(nfecra,5050)
  write(nfecra,5000)

!--> Correction du terme source explicite dans raysca pour permettre un
!    post-processing correct du terme source explicite
!    lorsque la variable transportee est la temperature
!    (pour les calculs en combustion la variable transportee est toujours
!    l'enthalpie)
else
  write(nfecra,5000)
endif

! Free memory
deallocate(viscf, viscb)
deallocate(smbrs, rovsdt)
deallocate(ckmel)
deallocate(tempk)
deallocate(coefap, coefbp)
deallocate(cofafp, cofbfp)
deallocate(flurds, flurdb)

!--------
! Formats
!--------

 1000 FORMAT (/, 3X,'** INFORMATIONS SUR LE TERME SOURCE RADIATIF',/,   &
           3X,'   -----------------------------------------' )
 1100 FORMAT (/, 3X,'   Calcul effectue en rayonnement transparent'  ,/)

 2010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE RAYONNEMENT EST ACTIVE AVEC LE MODELE DOM.           ',/,&
'@      LA VALEUR MINIMALE DU COEFFICIENT D ABSORPTION A EST  ',/,&
'@      EGALE A ', E14.5                                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE RAYONNEMENT EST ACTIVE AVEC LE MODELE P-1.           ',/,&
'@      LE COEFFICIENT D''ABSORBTION DOIT ETRE STRICTEMENT    ',/,&
'@      SUPERIEUR A ZERO.                                     ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3500 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE RAYONNEMENT EST ACTIVE.                              ',/,&
'@                                                            ',/,&
'@    Le modele ITHERM devrait valoir 1 (temperature) ou      ',/,&
'@      2 (enthalpie). On a ITHERM = ',i10                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4000 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT (FLUNET    NON RENSEIGNE)       ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )
 4100 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    LE FLUNET    N''EST PAS RENSEIGNEE POUR CERTAINES       ',/,&
'@        FACES DE BORD                                       ',/,&
'@                                                            ',/,&
'@        Valeur minimale ',E14.5                              ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usray5.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 5000 format('-----------------------------------------------------',   &
          '--------------')

 5010 format('Zone         Flux net radiatif (Watt) (normale',          &
          ' unitaire sortante)')

 5020 format(i6,13x,e10.4)

 5030 format('Flux net radiatif sur toutes les frontieres  Fnet = ',    &
           E10.4,' Watt')

 5040 format('Integrale volumique du terme source radiatif Srad = ',    &
           E10.4,' Watt')

 5050 format('(Si IDIVER = 1 ou 2 alors on doit avoir Srad = -Fnet)')

 6000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT APPROXIMATION P-1  (RAYDOM)     ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    LA LONGUEUR OPTIQUE DU MILIEU SEMI-TRANSPARENT          ',/,&
'@      DOIT AU MOINS ETRE DE L''ORDRE DE L''UNITE POUR ETRE  ',/,&
'@      DANS LE DOMAINE D''APPLICATION DE L''APPROXIMATION P-1',/,&
'@    CELA NE SEMBLE PAS ETRE LE CAS ICI.                     ',/,&
'@                                                            ',/,&
'@    LE COEFFICIENT D''ABSORPTION MINIMUM POUR ASSURER CETTE ',/,&
'@      LONGUEUR OPTIQUE EST XKMIN = ',E10.4                   ,/,&
'@    CETTE VALEUR N''EST PAS ATTEINTE POUR ', E10.4,'%       ',/,&
'@      DES CELLULES DU MAILLAGE.                             ',/,&
'@    LE POURCENTAGE DE CELLULES DU MAILLAGE POUR LESQUELLES  ',/,&
'@      ON ADMET QUE CETTE CONDITION SOIT VIOLEE EST IMPOSE   ',/,&
'@      PAR DEFAUT OU DANS USINI1 A XNP1MX = ', E10.4,'%      ',/,&
'@                                                            ',/,&
'@    Verifier les valeurs du coefficient d''absorption CK    ',/,&
'@      dans l''interface ou le modifier dans USRAY3.         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

end subroutine

