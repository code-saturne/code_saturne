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

subroutine dttvar &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iwarnp ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce ,                            &
   coefa  , coefb  , ckupdc , smacel )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DU PAS DE TEMPS LOCAL
! AFFICHAGE DES NOMBRES DE COURANT + FOURIER MINIMUM, MAXIMUM
! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)

! Sous programme utilise dans le cas une seule phase (ou
! si seule la phase 1 pilote le pas de temps)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,nvar)    !    !     !  source de masse                               !
!                  !    !     ! pour ivar=ipr, smacel=flux de masse            !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use cplsat
use cstnum
use cstphy
use optcal
use entsor
use parall
use ppppar
use ppthch
use ppincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iwarnp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)

! Local variables

character*8      cnom

integer          ifac, iel, icfmax, icfmin, idiff0, iconv0, isym, flid
integer          modntl
integer          ipcvis, ipcvst
integer          iflmas, iflmab
integer          icou, ifou , icoucf
integer          inc, iccocg
integer          nswrgp, imligp
integer          iivar
integer          nbrval, nclptr
integer          ipccou, ipcfou, ntcam1

double precision epsrgp, climgp, extrap
double precision cfmax,cfmin, w1min, w2min, w3min
double precision unpvdt, rom, dtloc
double precision xyzmax(3), xyzmin(3), vmin(1), vmax(1)
double precision dtsdtm
double precision hint

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:) :: wcf
double precision, allocatable, dimension(:) :: cofbft, coefbt, coefbr
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1, w2, w3, dtsdt0
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, crom
!===============================================================================

!===============================================================================
! 0.  INITIALISATION
!===============================================================================

! Pointers to the mass fluxes
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

! Allocate temporary arrays for the time-step resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(dam(ncelet))
allocate(cofbft(nfabor), coefbt(nfabor))

! Allocate other arrays, depending on user options
if (ippmod(icompf).ge.0) then
  allocate(wcf(ncelet))
endif

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))

ipcvis  = ipproc(iviscl)
ipcvst  = ipproc(ivisct)
call field_get_val_s(icrom, crom)
 call field_get_val_s(ibrom, brom)
ipccou  = ipproc(icour)
ipcfou  = ipproc(ifour)

if (ntlist.gt.0) then
  modntl = mod(ntcabs,ntlist)
elseif (ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
  modntl = 0
else
  modntl = 1
endif

if (                                                             &
   .not. (iconv(iu).ge.1.and.                                    &
           (iwarnp.ge.2.or.modntl.eq.0)) .and.                   &
   .not. (idiff(iu).ge.1.and.                                    &
           (iwarnp.ge.2.or.modntl.eq.0)) .and.                   &
   .not. (ippmod(icompf).ge.0.and.                               &
           (iwarnp.ge.2.or.modntl.eq.0)) .and.                   &
   .not. (idtvar.eq.-1.or.idtvar.eq.1.or.idtvar.eq.2.or.         &
           ((iwarnp.ge.2.or.modntl.eq.0).and.                    &
             (idiff(iu).ge.1.or.iconv(iu).ge.1                   &
                               .or.ippmod(icompf).ge.0)))        &
   ) then

  return

endif

!===============================================================================
! 1.  CALCUL DE LA LIMITATION EN COMPRESSIBLE
!===============================================================================

!     On commence par cela afin de disposer de VISCF VISCB comme
!       tableaux de travail.

  if (ippmod(icompf).ge.0) then

    call cfdttv                                                   &
    !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce ,                            &
   ckupdc , smacel ,                                              &
   wcf    ,                                                       &
   viscf  , viscb  , cofbft )

  endif


!===============================================================================
! 2. Compute the diffusivity at the faces
!===============================================================================

!     On s'en sert dans les divers matrdt suivants

!     "VITESSE" DE DIFFUSION FACETTE

if (idiff(iu).ge. 1) then
  do iel = 1, ncel
    w1(iel) = propce(iel,ipcvis) + idifft(iu)*propce(iel,ipcvst)
  enddo
  call viscfa(imvisf, w1, viscf, viscb)
  !==========
else
  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo
endif

!===============================================================================
! 3.  CONDITION LIMITE POUR MATRDT
!===============================================================================

do ifac = 1, nfabor

  if (bmasfl(ifac).lt.0.d0) then
    iel = ifabor(ifac)
    hint = idiff(iu)*(  propce(iel,ipcvis)                         &
                      + idifft(iu)*propce(iel,ipcvst))/distb(ifac)
    coefbt(ifac) = 0.d0
    cofbft(ifac) = hint
  else
    coefbt(ifac) = 1.d0
    cofbft(ifac) = 0.d0
  endif
enddo

!===============================================================================
! 4.  ALGORITHME INSTATIONNAIRE
!===============================================================================

if (idtvar.ge.0) then

!===============================================================================
! 4.1  PAS DE TEMPS VARIABLE A PARTIR DE COURANT ET FOURIER IMPOSES
!===============================================================================

  ! On calcule le pas de temps thermique max (meme en IDTVAR=0, pour affichage)
  ! DTTMAX = 1/SQRT(MAX(0+,gradRO.g/RO) -> W3

  if (iptlro.eq.1) then

    ! Allocate a temporary array for the gradient calculation
    allocate(grad(ncelet,3))
    allocate(coefbr(nfabor))

    do ifac = 1, nfabor
      coefbr(ifac) = 0.d0
    enddo

    nswrgp = nswrgr(ipr)
    imligp = imligr(ipr)
    iwarnp = iwarni(ipr)
    epsrgp = epsrgr(ipr)
    climgp = climgr(ipr)
    extrap = 0.d0

    iivar = 0
    inc   = 1
    iccocg = 1

    call grdcel                                                   &
    !==========
 ( iivar  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   crom, brom, coefbr ,                   &
   grad )

    do iel = 1, ncel
      w3(iel) = (grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz)&
           /crom(iel)
      w3(iel) = 1.d0/sqrt(max(epzero,w3(iel)))

    enddo

    ! Free memory
    deallocate(grad)
    deallocate(coefbr)

  endif

  if (idtvar.eq.1.or.idtvar.eq.2) then

    icou = 0
    ifou = 0

! 4.1.1 LIMITATION PAR LE COURANT
! =============================

    if (coumax.gt.0.d0.and.iconv(iu).ge.1) then

      ! ICOU = 1 marque l'existence d'une limitation par le COURANT
      icou = 1

      ! CONSTRUCTION DE U/DX           (COURANT       ) =W1

      idiff0 = 0

      ! Matrice a priori non symetrique
      isym = 2

      call matrdt &
      !==========
 (iconv(iu), idiff0, isym, coefbt, cofbft, imasfl, bmasfl, viscf, viscb, dam)

      do iel = 1, ncel
        rom = crom(iel)
        w1    (iel) = dam(iel)/(rom*volume(iel))
      enddo

      ! ---> CALCUL DE W1     = PAS DE TEMPS VARIABLE VERIFIANT
      !       LE NOMBRE DE COURANT         MAXIMUM PRESCRIT PAR L'UTILISATEUR

      do iel = 1, ncel
        w1    (iel) = coumax/max(w1    (iel), epzero)
      enddo

      ! PAS DE TEMPS UNIFORME : ON PREND LE MINIMUM DE LA CONTRAINTE

      if (idtvar.eq.1) then
        w1min = grand
        do iel = 1, ncel
          w1min = min(w1min,w1(iel))
        enddo
        if (irangp.ge.0) then
          call parmin (w1min)
          !==========
        endif
        do iel = 1, ncel
          w1(iel) = w1min
        enddo
      endif

    endif

    ! 4.1.2 LIMITATION PAR LE FOURIER
    ! =============================

    if (foumax.gt.0.d0.and.idiff(iu).ge.1) then

      ! IFOU = 1 marque l'existence d'une limitation par le FOURIER
      ifou = 1

      iconv0 = 0
      !                              2
      ! CONSTRUCTION DE      +2.NU/DX  (       FOURIER) =W2

      ! Matrice a priori symetrique
      isym = 1

      call matrdt &
      !==========
 (iconv0, idiff(iu), isym, coefbt, cofbft, imasfl, bmasfl, viscf, viscb, dam)

      do iel = 1, ncel
        rom = crom(iel)
        w2    (iel) = dam(iel)/(rom*volume(iel))
      enddo

      ! ---> CALCUL DE W2     = PAS DE TEMPS VARIABLE VERIFIANT
      !       LE NOMBRE DE         FOURIER MAXIMUM PRESCRIT PAR L'UTILISATEUR

      do iel = 1, ncel
        w2    (iel) = foumax/max(w2    (iel), epzero)
      enddo

      ! ---> PAS DE TEMPS UNIFORME : ON PREND LE MINIMUM DE LA CONTRAINTE

      if (idtvar.eq.1) then
        w2min = grand
        do iel = 1, ncel
          w2min = min(w2min,w2(iel))
        enddo
        if (irangp.ge.0) then
          call parmin (w2min)
          !==========
        endif
        do iel = 1, ncel
          w2(iel) = w2min
        enddo
      endif

    endif

    ! 4.1.3 LIMITATION POUR L'ALGORITHME COMPRESSIBLE
    ! =============================================
    !     Il est important de conserver WCF intact : on le reutilise
    !     plus bas pour l'affichage

    icoucf = 0
    if (coumax.gt.0.d0.and.ippmod(icompf).ge.0) then

      icoucf = 1

      ! ---> CALCUL DE DAM     = PAS DE TEMPS VARIABLE VERIFIANT
      !       LA CONTRAINTE CFL MAXIMUM PRESCRITE PAR L'UTILISATEUR

      do iel = 1, ncel
        dam(iel) = cflmmx/max( wcf(iel), epzero)
      enddo

      ! ---> PAS DE TEMPS UNIFORME : ON PREND LE MINIMUM DE LA CONTRAINTE

      if (idtvar.eq.1) then
        w3min = grand
        do iel = 1, ncel
          w3min = min(w3min,dam(iel))
        enddo
        if (irangp.ge.0) then
          call parmin (w3min)
          !==========
        endif
        do iel = 1, ncel
          dam(iel) = w3min
        enddo
      endif

    endif

    ! 4.1.4 ON PREND LA PLUS CONTRAIGNANTE DES LIMITATIONS
    ! ==================================================
    !    (le minimum des deux si elles existent et
    !     celle qui existe s'il n'en existe qu'une)

    if (icou.eq.1.and.ifou.eq.1) then
      do iel = 1, ncel
        w1(iel) = min(w1(iel),w2(iel))
      enddo
    elseif (icou.eq.0.and.ifou.eq.1) then
      do iel = 1, ncel
        w1(iel) = w2(iel)
      enddo
    endif

    !     En compressible, on prend obligatoirement
    !     en compte la limitation associee à la masse volumique.

    if (icoucf.eq.1) then
      do iel = 1, ncel
        w1(iel) = min(w1(iel),dam(iel))
      enddo
    endif

    ! 4.1.5 ON CALCULE EFFECTIVEMENT LE PAS DE TEMPS
    ! ============================================

    ! --->  MONTEE           PROGRESSIVE DU PAS DE TEMPS
    !              DESCENTE  IMMEDIATE   DU PAS DE TEMPS

    do iel = 1, ncel
      if (w1    (iel).ge.dt(iel)) then
        unpvdt = 1.d0+varrdt
        dt(iel) = min(unpvdt*dt(iel), w1(iel))
      else
        dt(iel) =                     w1(iel)
      endif
    enddo


    ! 4.1.6 ON LIMITE PAR LE PAS DE TEMPS "THERMIQUE" MAX
    ! =================================================
    !     DTTMAX = W3 = 1/SQRT(MAX(0+,gradRO.g/RO)
    !     on limite le pas de temps a DTTMAX

    if (iptlro.eq.1) then

      ! On clippe le pas de temps a DTTMAX (affiche dans ecrlis)

      nclptr = 0

      vmin(1) = dt(1)
      vmax(1) = dt(1)

      do iel = 1, ncel
        vmin(1) = min(vmin(1),dt(iel))
        vmax(1) = max(vmax(1),dt(iel))
        if (dt(iel).gt.w3(iel)) then
          nclptr = nclptr +1
          dt(iel) = w3(iel)
        endif
      enddo

      call log_iteration_clipping('dt (clip/dtrho)', 1, 0, nclptr, vmin, vmax)

      ! ---> PAS DE TEMPS UNIFORME : on reuniformise le pas de temps

      if (idtvar.eq.1) then

        w3min = grand
        do iel = 1, ncel
          w3min = min(w3min,dt(iel))
        enddo
        if (irangp.ge.0) then
          call parmin (w3min)
          !==========
        endif

        do iel = 1, ncel
          dt(iel) = w3min
        enddo
      endif

    endif

    ! 4.1.7 ON CLIPPE LE PAS DE TEMPS PAR RAPPORT A DTMIN ET DTMAX
    ! ==========================================================

    icfmin = 0
    icfmax = 0

    call field_get_id('dt', flid)

    if (idtvar.eq.1) then

      dtloc = dt(1)
      if (dtloc.gt.dtmax) then
        dtloc = dtmax
        icfmax = ncel
      endif
      if (dtloc.lt.dtmin) then
        dtloc = dtmin
        icfmin = ncel
      endif

      ntcam1 = ntcabs - 1
      call cplsyn (ntmabs, ntcam1, dtloc)
      !==========
      if (ntmabs.lt.ntcabs) then
        call csexit(1)
      endif

      call log_iteration_clipping_field(flid, icfmin, icfmax, dt, dt)

      ttcabs = ttcabs + (dtloc - dt(1))
      if (imobil.eq.1) then
        ttcmob = ttcmob + (dtloc - dt(1))
      endif

      do iel = 1, ncel
        dt(iel) = dtloc
      enddo

    else

      vmin(1) = dt(1)
      vmax(1) = dt(1)

      do iel = 1, ncel

        vmin(1) = min(vmin(1),dt(iel))
        vmax(1) = max(vmax(1),dt(iel))

        if (dt(iel).gt.dtmax) then
          icfmax = icfmax +1
          dt(iel) = dtmax
        endif
        if (dt(iel).lt.dtmin) then
          icfmin = icfmin +1
          dt(iel) = dtmin
        endif

      enddo

      call log_iteration_clipping_field(flid, icfmin, icfmax, vmin, vmax)

    endif


    if (iwarnp.ge.2) then
      if (irangp.ge.0) then
        call parcpt (icfmin)
        !==========
        call parcpt (icfmax)
        !==========
      endif
      write (nfecra,1003) icfmin,dtmin,icfmax,dtmax
    endif

  endif

  ! Rapport DT sur DTmax lie aux effets de densite (affichage dans ecrlis)

  if (iptlro.eq.1) then

    dtsdtm = 0.d0

    allocate(dtsdt0(ncel))

    do iel = 1, ncel
      dtsdt0(iel) = dt(iel)/w3(iel)
    enddo

    call log_iteration_add_array('Dt/Dtrho max', 'criterion',       &
                                 MESH_LOCATION_CELLS, .true.,       &
                                 1, dtsdt0)

    deallocate(dtsdt0)

  endif

!===============================================================================
! 4.2  CALCUL DU NOMBRE DE COURANT POUR AFFICHAGE
!===============================================================================

  if (iconv(iu).ge.1 .and. (iwarnp.ge.2.or.modntl.eq.0)) then

    idiff0 = 0
    cnom   =' COURANT'

    ! CONSTRUCTION DE U/DX           (COURANT       ) =W1

    ! MATRICE A PRIORI NON SYMETRIQUE

    isym = 2

    call matrdt &
    !==========
 (iconv(iu), idiff0, isym, coefbt, cofbft, imasfl, bmasfl, viscf, viscb, dam)

    do iel = 1, ncel
      rom = crom(iel)
      w1    (iel) = dam(iel)/(rom*volume(iel))
    enddo

    ! CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax = -grand
    cfmin =  grand
    icfmax= 1
    icfmin= 1

    do iel = 1, ncel
      propce(iel,ipccou) = w1(iel)*dt(iel)
    enddo

    if (iwarnp.ge.2) then

      do iel = 1, ncel
        if (propce(iel,ipccou).le.cfmin) then
          cfmin  = propce(iel,ipccou)
          icfmin = iel
        endif
        if (propce(iel,ipccou).ge.cfmax) then
          cfmax  = propce(iel,ipccou)
          icfmax = iel
        endif
      enddo

      xyzmin(1) = xyzcen(1,icfmin)
      xyzmin(2) = xyzcen(2,icfmin)
      xyzmin(3) = xyzcen(3,icfmin)
      xyzmax(1) = xyzcen(1,icfmax)
      xyzmax(2) = xyzcen(2,icfmax)
      xyzmax(3) = xyzcen(3,icfmax)

      if (irangp.ge.0) then
        nbrval = 3
        call parmnl (nbrval, cfmin, xyzmin)
        !==========
        call parmxl (nbrval, cfmax, xyzmax)
        !==========
      endif

      write(nfecra,1001) cnom, cfmax, xyzmax(1), xyzmax(2), xyzmax(3)
      write(nfecra,1002) cnom, cfmin, xyzmin(1), xyzmin(2), xyzmin(3)

    endif

  endif

!===============================================================================
! 4.3  CALCUL DU NOMBRE DE FOURIER POUR AFFICHAGE
!===============================================================================

  if (idiff(iu).ge.1 .and. (iwarnp.ge.2.or.modntl.eq.0)) then

    iconv0 = 0
    CNOM   =' FOURIER'

    !                              2
    ! CONSTRUCTION DE      +2.NU/DX  (       FOURIER) =W1

    ! MATRICE A PRIORI SYMETRIQUE

    isym = 1

    call matrdt &
    !==========
 (iconv0, idiff(iu), isym, coefbt, cofbft, imasfl, bmasfl, viscf, viscb, dam)

    do iel = 1, ncel
      rom = crom(iel)
      w1    (iel) = dam(iel)/(rom*volume(iel))
    enddo

    ! CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax  = -grand
    cfmin  =  grand
    icfmax = 0
    icfmin = 0

    do iel = 1, ncel
      propce(iel,ipcfou) = w1(iel)*dt(iel)
    enddo

    if (iwarnp.ge.2) then

      do iel = 1, ncel
        if (propce(iel,ipcfou).le.cfmin) then
          cfmin  = propce(iel,ipcfou)
          icfmin = iel
        endif
        if (propce(iel,ipcfou).ge.cfmax) then
          cfmax  = propce(iel,ipcfou)
          icfmax = iel
        endif
      enddo

      xyzmin(1) = xyzcen(1,icfmin)
      xyzmin(2) = xyzcen(2,icfmin)
      xyzmin(3) = xyzcen(3,icfmin)
      xyzmax(1) = xyzcen(1,icfmax)
      xyzmax(2) = xyzcen(2,icfmax)
      xyzmax(3) = xyzcen(3,icfmax)

      if (irangp.ge.0) then
        nbrval = 3
        call parmnl (nbrval, cfmin, xyzmin)
        !==========
        call parmxl (nbrval, cfmax, xyzmax)
        !==========
      endif

      write(nfecra,1001) cnom,cfmax,xyzmax(1),xyzmax(2),xyzmax(3)
      write(nfecra,1002) cnom,cfmin,xyzmin(1),xyzmin(2),xyzmin(3)

    endif

  endif

!===============================================================================
! 4.4  CALCUL DU NOMBRE DE COURANT/FOURIER POUR AFFICHAGE
!===============================================================================

!     En incompressible uniquement (en compressible, on preferera
!       afficher la contrainte liee a la masse volumique)

  if ((iwarnp.ge.2.or.modntl.eq.0)                         &
      .and. (idiff(iu).ge.1.or.iconv(iu).ge.1)             &
      .and. (ippmod(icompf).lt.0))              then

    cnom   =' COU/FOU'

    !                              2
    ! CONSTRUCTION DE U/DX +2.NU/DX  (COURANT +FOURIER) =W1

    ! MATRICE A PRIORI NON SYMETRIQUE

    isym = 1
    if (iconv(iu).gt.0) isym = 2

    call matrdt &
    !==========
 (iconv(iu), idiff(iu), isym, coefbt, cofbft, imasfl, bmasfl, viscf, viscb, dam)

    do iel = 1, ncel
      rom = crom(iel)
      w1    (iel) = dam(iel)/(rom*volume(iel))
    enddo

    ! CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax  = -grand
    cfmin  =  grand
    icfmax = 0
    icfmin = 0

    do iel = 1, ncel
      w2(iel) = w1(iel)*dt(iel)
    enddo

    call log_iteration_add_array('Courant/Fourier', 'criterion',    &
                                 MESH_LOCATION_CELLS, .true.,       &
                                 1, w2)

    if (iwarnp.ge.2) then

      do iel = 1, ncel
        if (w2(iel).le.cfmin) then
          cfmin  = w2(iel)
          icfmin = iel
        endif
        if (w2(iel).ge.cfmax) then
          cfmax  = w2(iel)
          icfmax = iel
        endif
      enddo

      xyzmin(1) = xyzcen(1,icfmin)
      xyzmin(2) = xyzcen(2,icfmin)
      xyzmin(3) = xyzcen(3,icfmin)
      xyzmax(1) = xyzcen(1,icfmax)
      xyzmax(2) = xyzcen(2,icfmax)
      xyzmax(3) = xyzcen(3,icfmax)

      if (irangp.ge.0) then
        nbrval = 3
        call parmnl (nbrval, cfmin, xyzmin)
        !==========
        call parmxl (nbrval, cfmax, xyzmax)
        !==========
      endif

      write(nfecra,1001) cnom,cfmax,xyzmax(1),xyzmax(2),xyzmax(3)
      write(nfecra,1002) cnom,cfmin,xyzmin(1),xyzmin(2),xyzmin(3)

    endif

  endif

!===============================================================================
! 4.5  CALCUL DE LA CONTRAINTE CFL DE LA MASSE VOL. POUR AFFICHAGE
!===============================================================================

  ! En Compressible uniquement

  if ((iwarnp.ge.2.or.modntl.eq.0) .and. (ippmod(icompf).ge.0)) then

    cnom =' CFL/MAS'

    ! CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax  = -grand
    cfmin  =  grand
    icfmax = 0
    icfmin = 0

    do iel = 1, ncel
      w2(iel) = wcf(iel)*dt(iel)
    enddo

    call log_iteration_add_array('CFL / Mass ', 'criterion',        &
                                 MESH_LOCATION_CELLS, .true.,       &
                                 1, w2)

    if (iwarnp.ge.2) then

      do iel = 1, ncel
        if (w2(iel).le.cfmin) then
          cfmin  = w2(iel)
          icfmin = iel
        endif
        if (w2(iel).ge.cfmax) then
          cfmax  = w2(iel)
          icfmax = iel
        endif
      enddo

      xyzmin(1) = xyzcen(1,icfmin)
      xyzmin(2) = xyzcen(2,icfmin)
      xyzmin(3) = xyzcen(3,icfmin)
      xyzmax(1) = xyzcen(1,icfmax)
      xyzmax(2) = xyzcen(2,icfmax)
      xyzmax(3) = xyzcen(3,icfmax)

      if (irangp.ge.0) then
        nbrval = 3
        call parmnl (nbrval, cfmin, xyzmin)
        !==========
        call parmxl (nbrval, cfmax, xyzmax)
        !==========
      endif

      write(nfecra,1001) cnom,cfmax,xyzmax(1),xyzmax(2),xyzmax(3)
      write(nfecra,1002) cnom,cfmin,xyzmin(1),xyzmin(2),xyzmin(3)

    endif

  endif

!===============================================================================
! 5.   ALGORITHME STATIONNAIRE
!===============================================================================
else

  isym = 1
  if (iconv(iu).gt.0) isym = 2

  call matrdt &
  !==========
 ( iconv(iu), idiff(iu), isym,                                  &
   coefb(1,iclrtp(iu,icoef)), coefb(1,iclrtp(iu,icoeff)),       &
   imasfl, bmasfl,                                              &
   viscf, viscb, dt )

  do iel = 1, ncel
    dt(iel) =  relaxv(iu)*crom(iel)                    &
              *volume(iel)/max(dt(iel),epzero)
  enddo

endif

! Free memory
deallocate(viscf, viscb)
deallocate(dam)
deallocate(cofbft)
if (allocated(wcf)) deallocate(wcf)
deallocate(w1, w2, w3)

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1001 FORMAT (/,A8,' MAX= ',E11.4,                                     &
 ' EN ',E11.4,' ',E11.4,' ',E11.4)
 1002 FORMAT (  A8,' MIN= ',E11.4,                                     &
 ' EN ',E11.4,' ',E11.4,' ',E11.4)
 1003 FORMAT (/,'CLIPPINGS DE DT : ',                                  &
                             I10,' A ',E11.4,', ',I10,' A ',E11.4)

#else

 1001 FORMAT (/,A8,' MAX= ',E11.4,                                     &
 ' IN ',E11.4,' ',E11.4,' ',E11.4)
 1002 FORMAT (  A8,' MIN= ',E11.4,                                     &
 ' IN ',E11.4,' ',E11.4,' ',E11.4)
 1003 FORMAT (/,'DT CLIPPING : ',                                      &
                             I10,' A ',E11.4,', ',I10,' A ',E11.4)

#endif

!----
! FIN
!----

return

end subroutine
