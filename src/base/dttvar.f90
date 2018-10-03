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
! Purpose:
! -------

!> \file dttvar.f90
!> \brief Compute the local time step and add the Courant and Fourier number to
!the log.
!>
!> This function has access to the boundary face type, except for the first time
!> step.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     iwarnp        verbosity
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!>                               (see \ref cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!>                               source terms or mass rate (see
!>                               \ref cs_user_mass_source_terms)
!_______________________________________________________________________________

subroutine dttvar &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iwarnp ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
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
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iwarnp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)

! Local variables

character(len=8) :: cnom
character(len=80) :: fname

integer          ifac, iel, icfmax(1), icfmin(1), idiff0, iconv0, isym, flid
integer          modntl
integer          iflmas, iflmab
integer          icou, ifou , icoucf
integer          inc, iccocg
integer          nswrgp, imligp
integer          f_id
integer          nbrval, nclptr
integer          ntcam1
integer          ii

double precision epsrgp, climgp, extrap
double precision cfmax,cfmin, w1min, w2min, w3min
double precision unpvdt, dtloc
double precision xyzmax(3), xyzmin(3), vmin(1), vmax(1)
double precision dtsdtm
double precision hint
double precision mult
double precision prt

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:) :: wcf
double precision, allocatable, dimension(:) :: cofbft, coefbt, coefbr
double precision, dimension(:,:,:), pointer :: coefbv, cofbfv
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1, w2, w3, dtsdt0
double precision, dimension(:), pointer :: imasfl, bmasfl, imasflt, bmasflt
double precision, dimension(:), pointer :: brom, crom, cromt
double precision, dimension(:), pointer :: viscl, visct, cpro_cour, cpro_four

type(var_cal_opt) :: vcopt_u, vcopt_p, vcopt_t

!===============================================================================

!===============================================================================
! 0. Initialisation
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

call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)
call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)
call field_get_val_s(icour, cpro_cour)
call field_get_val_s(ifour, cpro_four)

if (ntlist.gt.0) then
  modntl = mod(ntcabs,ntlist)
elseif (ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
  modntl = 0
else
  modntl = 1
endif

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

if (.not. (vcopt_u%iconv.ge.1 .and. icour.ge.1) .and.              &
    .not. (vcopt_u%idiff.ge.1 .and. ifour.ge.1) .and.              &
    .not. ((vcopt_u%iconv.ge.1 .or.vcopt_u%idiff.ge.1).and.        &
           (iwarnp.ge.2.or.modntl.eq.0)) .and.                     &
    .not. (ippmod(icompf).ge.0.and.                                &
           (iwarnp.ge.2.or.modntl.eq.0)) .and.                     &
    .not. (idtvar.eq.-1.or.idtvar.eq.1.or.idtvar.eq.2.or.          &
           ((iwarnp.ge.2.or.modntl.eq.0).and.                      &
            (vcopt_u%idiff.ge.1.or.vcopt_u%iconv.ge.1.or.          &
             ippmod(icompf).ge.0)))                                &
   ) then

  return

endif

!===============================================================================
! 1. Compute CFL like condition on the time step for positive density for the
!    compressible module
!===============================================================================

if (ippmod(icompf).ge.0) then

  call cfdttv                                                   &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                          &
   icepdc , icetsm , itypsm ,                                   &
   dt     ,                                                     &
   ckupdc , smacel ,                                            &
   wcf    ,                                                     &
   viscf  , viscb  , cofbft )

endif

!===============================================================================
! 2. Compute the diffusivity at the faces
!===============================================================================

if (vcopt_u%idiff.ge.1) then
  do iel = 1, ncel
    w1(iel) = viscl(iel) + vcopt_u%idifft*visct(iel)
  enddo
  call viscfa(imvisf, w1, viscf, viscb)

else
  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo
endif

!===============================================================================
! 3. Boundary condition for matrdt
!===============================================================================

if (idtvar.ge.0) then

  do ifac = 1, nfabor

    if (bmasfl(ifac).lt.0.d0) then
      iel = ifabor(ifac)
      hint = vcopt_u%idiff*(  viscl(iel)                         &
                        + vcopt_u%idifft*visct(iel))/distb(ifac)
      coefbt(ifac) = 0.d0
      cofbft(ifac) = hint
    else
      coefbt(ifac) = 1.d0
      cofbft(ifac) = 0.d0
    endif
  enddo

else

  ! TODO for steady algorithm, check if using the third of the trace
  !      is appropriate, or if a better solution is available
  !      (algorithm was probably broken since velocity components are
  !      coupled)

  mult = 1.d0 / 3.d0

  call field_get_coefb_v (ivarfl(iu), coefbv)
  call field_get_coefbf_v(ivarfl(iu), cofbfv)

  do ifac = 1, nfabor
    coefbt(ifac) = (coefbv(1,1,ifac)+coefbv(2,2,ifac)+coefbv(3,3,ifac)) * mult
    cofbft(ifac) = (cofbfv(1,1,ifac)+cofbfv(2,2,ifac)+cofbfv(3,3,ifac)) * mult
  enddo

endif

!===============================================================================
! 4. Steady algorithm
!===============================================================================

if (idtvar.ge.0) then

!===============================================================================
! 4.1 Variable time step from imposed Courant and Fourier
!===============================================================================

  ! On calcule le pas de temps thermique max (meme en IDTVAR=0, pour affichage)
  ! dttmax = 1/sqrt(max(0+,gradRO.g/RO) -> w3

  if (iptlro.eq.1) then

    ! Allocate a temporary array for the gradient calculation
    allocate(grad(3,ncelet))
    allocate(coefbr(nfabor))

    do ifac = 1, nfabor
      coefbr(ifac) = 0.d0
    enddo

    nswrgp = vcopt_p%nswrgr
    imligp = vcopt_p%imligr
    iwarnp = vcopt_p%iwarni
    epsrgp = vcopt_p%epsrgr
    climgp = vcopt_p%climgr
    extrap = 0.d0

    f_id = -1
    inc   = 1
    iccocg = 1

    call gradient_s &
 ( f_id   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , epsrgp , climgp , extrap ,                            &
   crom, brom, coefbr ,                                           &
   grad )

    do iel = 1, ncel
      w3(iel) = (grad(1,iel)*gx + grad(2,iel)*gy + grad(3,iel)*gz)&
              / crom(iel)
      w3(iel) = 1.d0/sqrt(max(epzero,w3(iel)))

    enddo

    ! Free memory
    deallocate(grad)
    deallocate(coefbr)

  endif

  if (idtvar.eq.1.or.idtvar.eq.2) then

    icou = 0
    ifou = 0

    ! 4.1.1 Courant limitation
    ! ========================

    if (coumax.gt.0.d0.and.vcopt_u%iconv.ge.1) then

      ! ICOU = 1 marque l'existence d'une limitation par le COURANT
      icou = 1

      ! CONSTRUCTION DE U/DX (COURANT) = W1

      idiff0 = 0

      ! Matrice a priori non symetrique
      isym = 2

      call matrdt &
     ( vcopt_u%iconv, idiff0, isym, coefbt, cofbft, imasfl, bmasfl, viscf, &
       viscb, dam)

      do iel = 1, ncel
        w1(iel) = dam(iel)/(crom(iel)*volume(iel))
      enddo

      ! ---> CALCUL DE W1 = PAS DE TEMPS VARIABLE VERIFIANT
      !       LE NOMBRE DE COURANT MAXIMUM PRESCRIT PAR L'UTILISATEUR

      do iel = 1, ncel
        w1(iel) = coumax/max(w1(iel), epzero)
      enddo

      ! PAS DE TEMPS UNIFORME : ON PREND LE MINIMUM DE LA CONTRAINTE

      if (idtvar.eq.1) then
        w1min = grand
        do iel = 1, ncel
          w1min = min(w1min,w1(iel))
        enddo
        if (irangp.ge.0) then
          call parmin(w1min)
        endif
        do iel = 1, ncel
          w1(iel) = w1min
        enddo
      endif

    endif

    ! 4.1.2 Fourier limitation
    ! ========================

    if (foumax.gt.0.d0.and.vcopt_u%idiff.ge.1) then

      ! IFOU = 1 marque l'existence d'une limitation par le FOURIER
      ifou = 1

      iconv0 = 0
      !                              2
      ! CONSTRUCTION DE      +2.NU/DX  (       FOURIER) =W2
      ! Matrice a priori symetrique
      isym = 1

      call matrdt &
      !==========
 (iconv0, vcopt_u%idiff, isym, coefbt, cofbft, imasfl, bmasfl, &
  viscf, viscb, dam)

      do iel = 1, ncel
        w2(iel) = dam(iel)/(crom(iel)*volume(iel))
      enddo


      ! ---> CALCUL DE W2     = PAS DE TEMPS VARIABLE VERIFIANT
      !       LE NOMBRE DE         FOURIER MAXIMUM PRESCRIT PAR L'UTILISATEUR

      do iel = 1, ncel
        w2(iel) = foumax/max(w2(iel), epzero)
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
    !     en compte la limitation associée à la masse volumique.

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
      if (w1(iel).ge.dt(iel)) then
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

    icfmin(1) = 0
    icfmax(1) = 0

    call field_get_id('dt', flid)

    if (idtvar.eq.1) then

      dtloc = dt(1)
      if (dtloc.gt.dtmax) then
        dtloc = dtmax
        icfmax(1) = ncel
      endif
      if (dtloc.lt.dtmin) then
        dtloc = dtmin
        icfmin(1) = ncel
      endif

      ntcam1 = ntcabs - 1
      call cplsyn (ntmabs, ntcam1, dtloc)
      !==========
      if (ntmabs.lt.ntcabs) then
        call csexit(1)
      endif

      call log_iteration_clipping_field(flid, icfmin(1), icfmax(1), dt, dt,icfmin(1), icfmax(1))

      ttcabs = ttcabs + (dtloc - dt(1))
      if (iturbo.eq.2) then
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
          icfmax(1) = icfmax(1) +1
          dt(iel) = dtmax
        endif
        if (dt(iel).lt.dtmin) then
          icfmin(1) = icfmin(1) +1
          dt(iel) = dtmin
        endif

      enddo

      call log_iteration_clipping_field(flid, icfmin(1), icfmax(1), vmin, vmax,icfmin(1), icfmax(1))

    endif


    if (iwarnp.ge.2) then
      if (irangp.ge.0) then
        call parcpt (icfmin(1))
        !==========
        call parcpt (icfmax(1))
        !==========
      endif
      write (nfecra,1003) icfmin(1),dtmin,icfmax(1),dtmax
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

  if (vcopt_u%iconv.ge.1 .and. icour.ge.1) then

    idiff0 = 0
    cnom   =' COURANT'

    ! CONSTRUCTION DE U/DX           (COURANT       ) =W1

    ! MATRICE A PRIORI NON SYMETRIQUE

    isym = 2

    call matrdt &
    !==========
 (vcopt_u%iconv, idiff0, isym, coefbt, cofbft, imasfl, bmasfl, &
  viscf, viscb, dam)

    do iel = 1, ncel
      w1(iel) = dam(iel)/(crom(iel)*volume(iel))
    enddo

    ! CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax = -grand
    cfmin =  grand
    icfmax(1)= 1
    icfmin(1)= 1

    do iel = 1, ncel
      cpro_cour(iel) = w1(iel)*dt(iel)
    enddo

    if (iwarnp.ge.2) then

      do iel = 1, ncel
        if (cpro_cour(iel).le.cfmin) then
          cfmin  = cpro_cour(iel)
          icfmin(1) = iel
        endif
        if (cpro_cour(iel).ge.cfmax) then
          cfmax  = cpro_cour(iel)
          icfmax(1) = iel
        endif
      enddo

      xyzmin(1) = xyzcen(1,max(icfmin(1), 1))
      xyzmin(2) = xyzcen(2,max(icfmin(1), 1))
      xyzmin(3) = xyzcen(3,max(icfmin(1), 1))
      xyzmax(1) = xyzcen(1,max(icfmax(1), 1))
      xyzmax(2) = xyzcen(2,max(icfmax(1), 1))
      xyzmax(3) = xyzcen(3,max(icfmax(1), 1))

      if (irangp.ge.0) then
        nbrval = 3
        call parmnl (nbrval, cfmin, xyzmin)
        call parmxl (nbrval, cfmax, xyzmax)
      endif

      if (icfmin(1).gt.0) then
        write(nfecra,1001) cnom,cfmax,xyzmax(1),xyzmax(2),xyzmax(3)
      else
        write(nfecra,*) cnom, "Too big to be displayed"
      endif

      if (icfmax(1).gt.0) then
        write(nfecra,1002) cnom,cfmin,xyzmin(1),xyzmin(2),xyzmin(3)
      else
        write(nfecra,*) cnom, "Too big to be displayed"
      endif

    endif

  endif

!===============================================================================
! 4.3  CALCUL DU NOMBRE DE FOURIER POUR AFFICHAGE
!===============================================================================

  if (vcopt_u%idiff.ge.1 .and. ifour.ge.1) then

    iconv0 = 0
    cnom   =' FOURIER'

    !                              2
    ! CONSTRUCTION DE      +2.NU/DX  (       FOURIER) =W1

    ! MATRICE A PRIORI SYMETRIQUE

    isym = 1

    call matrdt &
    !==========
 (iconv0, vcopt_u%idiff, isym, coefbt, cofbft, imasfl, bmasfl, &
  viscf, viscb, dam)

    do iel = 1, ncel
      w1(iel) = dam(iel)/(crom(iel)*volume(iel))
    enddo

    ! CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax  = -grand
    cfmin  =  grand
    icfmax(1) = 0
    icfmin(1) = 0

    do iel = 1, ncel
      cpro_four(iel) = w1(iel)*dt(iel)
    enddo

    if (iwarnp.ge.2) then

      do iel = 1, ncel
        if (cpro_four(iel).le.cfmin) then
          cfmin  = cpro_four(iel)
          icfmin(1) = iel
        endif
        if (cpro_four(iel).ge.cfmax) then
          cfmax  = cpro_four(iel)
          icfmax(1) = iel
        endif
      enddo

      xyzmin(1) = xyzcen(1,icfmin(1))
      xyzmin(2) = xyzcen(2,icfmin(1))
      xyzmin(3) = xyzcen(3,icfmin(1))
      xyzmax(1) = xyzcen(1,icfmax(1))
      xyzmax(2) = xyzcen(2,icfmax(1))
      xyzmax(3) = xyzcen(3,icfmax(1))

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

  if (      (vcopt_u%idiff.ge.1.or.vcopt_u%iconv.ge.1) &
      .and. (ippmod(icompf).lt.0)) then

    cnom   =' COU/FOU'

    !                              2
    ! CONSTRUCTION DE U/DX +2.NU/DX  (COURANT +FOURIER) =W1

    ! MATRICE A PRIORI NON SYMETRIQUE

    isym = 1
    if (vcopt_u%iconv.gt.0) isym = 2

    call matrdt &
    !==========
 (vcopt_u%iconv, vcopt_u%idiff, isym, coefbt, cofbft, imasfl, bmasfl, &
  viscf, viscb, dam)

    do iel = 1, ncel
      w1(iel) = dam(iel)/(crom(iel)*volume(iel))
    enddo

    ! CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax  = -grand
    cfmin  =  grand
    icfmax(1) = 0
    icfmin(1) = 0

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
          icfmin(1) = iel
        endif
        if (w2(iel).ge.cfmax) then
          cfmax  = w2(iel)
          icfmax(1) = iel
        endif
      enddo

      xyzmin(1) = xyzcen(1,max(icfmin(1), 1))
      xyzmin(2) = xyzcen(2,max(icfmin(1), 1))
      xyzmin(3) = xyzcen(3,max(icfmin(1), 1))
      xyzmax(1) = xyzcen(1,max(icfmax(1), 1))
      xyzmax(2) = xyzcen(2,max(icfmax(1), 1))
      xyzmax(3) = xyzcen(3,max(icfmax(1), 1))

      if (irangp.ge.0) then
        nbrval = 3
        call parmnl (nbrval, cfmin, xyzmin)
        call parmxl (nbrval, cfmax, xyzmax)
      endif

      if (icfmin(1).gt.0) then
        write(nfecra,1001) cnom,cfmax,xyzmax(1),xyzmax(2),xyzmax(3)
      else
        write(nfecra,*) cnom, "Too big to be displayed"
      endif

      if (icfmax(1).gt.0) then
        write(nfecra,1002) cnom,cfmin,xyzmin(1),xyzmin(2),xyzmin(3)
      else
        write(nfecra,*) cnom, "Too big to be displayed"
      endif

    endif

  endif

!===============================================================================
! 4.5  CALCUL DE LA CONTRAINTE CFL DE LA MASSE VOL. POUR AFFICHAGE
!===============================================================================

  ! En Compressible uniquement

  if (ippmod(icompf).ge.0) then

    cnom =' CFL/MAS'

    ! CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax  = -grand
    cfmin  =  grand
    icfmax(1) = 0
    icfmin(1) = 0

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
          icfmin(1) = iel
        endif
        if (w2(iel).ge.cfmax) then
          cfmax  = w2(iel)
          icfmax(1) = iel
        endif
      enddo

      xyzmin(1) = xyzcen(1,max(icfmin(1), 1))
      xyzmin(2) = xyzcen(2,max(icfmin(1), 1))
      xyzmin(3) = xyzcen(3,max(icfmin(1), 1))
      xyzmax(1) = xyzcen(1,max(icfmax(1), 1))
      xyzmax(2) = xyzcen(2,max(icfmax(1), 1))
      xyzmax(3) = xyzcen(3,max(icfmax(1), 1))

      if (irangp.ge.0) then
        nbrval = 3
        call parmnl (nbrval, cfmin, xyzmin)
        call parmxl (nbrval, cfmax, xyzmax)
      endif

      if (icfmin(1).gt.0) then
        write(nfecra,1001) cnom,cfmax,xyzmax(1),xyzmax(2),xyzmax(3)
      else
        write(nfecra,*) cnom, "Too big to be displayed"
      endif

      if (icfmax(1).gt.0) then
        write(nfecra,1002) cnom,cfmin,xyzmin(1),xyzmin(2),xyzmin(3)
      else
        write(nfecra,*) cnom, "Too big to be displayed"
      endif

    endif

  endif

!===============================================================================
! 4.6 PRINT COMBINED COURANT/FOURIER NUMBER FOR TRANSPORTED SCALARS
!===============================================================================

  do ii = 1, nscal

    ! Get field parameters
    call field_get_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt_t)

    if ( vcopt_t%idiff.ge.1 .or. vcopt_t%iconv.ge.1 ) then

      ! Get mass fluxes
      call field_get_key_int(ivarfl(isca(ii)), kimasf, iflmas)
      call field_get_key_int(ivarfl(isca(ii)), kbmasf, iflmab)
      call field_get_val_s(iflmas, imasflt)
      call field_get_val_s(iflmab, bmasflt)

      ! Get density
      call field_get_key_int(ivarfl(isca(ii)), kromsl, flid)
      if (flid.gt.-1) then
        call field_get_val_s(flid, cromt)
      else
        call field_get_val_s(icrom, cromt)
      endif

      ! Compute diffusivity
      if (vcopt_t%idiff.ge.1) then

        call field_get_key_double(ivarfl(isca(ii)), ksigmas, prt)
        call field_get_key_int(ivarfl(isca(ii)), kivisl, flid)
        if (flid.gt.-1) then
          call field_get_val_s(flid, viscl)
          do iel = 1, ncel
            w1(iel) = viscl(iel) + vcopt_t%idifft*visct(iel)/prt
          enddo
        else
          do iel = 1, ncel
            w1(iel) = visls0(ii) + vcopt_t%idifft*visct(iel)/prt
          enddo
        endif
        call viscfa(imvisf, w1, viscf, viscb)

      else

        do ifac = 1, nfac
          viscf(ifac) = 0.d0
        enddo
        do ifac = 1, nfabor
          viscb(ifac) = 0.d0
        enddo

      endif

      ! Boundary conditions for matrdt
      do ifac = 1, nfabor
        if (bmasflt(ifac).lt.0.d0) then

          iel = ifabor(ifac)
          if (flid.gt.-1) then
            hint = vcopt_t%idiff*( viscl(iel)                                 &
                                 + vcopt_t%idifft*visct(iel)/prt)/distb(ifac)
          else
            hint = vcopt_t%idiff*( visls0(ii)                                 &
                                 + vcopt_t%idifft*visct(iel)/prt)/distb(ifac)
          endif
          coefbt(ifac) = 0.d0
          cofbft(ifac) = hint

        else

          coefbt(ifac) = 1.d0
          cofbft(ifac) = 0.d0

        endif
      enddo

      ! MATRICE A PRIORI NON SYMETRIQUE
      isym = 1
      if (vcopt_t%iconv.gt.0) isym = 2

      call matrdt &
      !==========
 (vcopt_t%iconv, vcopt_t%idiff, isym, coefbt, cofbft, imasflt, bmasflt, &
  viscf, viscb, dam)

      do iel = 1, ncel
        w1(iel) = dam(iel)/(cromt(iel)*volume(iel))
      enddo

      do iel = 1, ncel
        w2(iel) = w1(iel)*dt(iel)
      enddo

      call field_get_name(ivarfl(isca(ii)), fname)
      call log_iteration_add_array(fname(1:15), 'criterion',    &
                                   MESH_LOCATION_CELLS, .true.,       &
                                   1, w2)

    endif

  enddo

!===============================================================================
! 5.   ALGORITHME STATIONNAIRE
!===============================================================================
else

  isym = 1
  if (vcopt_u%iconv.gt.0) isym = 2

  call matrdt &
  !==========
 ( vcopt_u%iconv, vcopt_u%idiff, isym,                          &
   coefbt, cofbft,                                              &
   imasfl, bmasfl,                                              &
   viscf, viscb, dt )

  do iel = 1, ncel
    dt(iel) =  vcopt_u%relaxv*crom(iel)*volume(iel)/max(dt(iel),epzero)
  enddo

endif

! Free memory
deallocate(viscf, viscb)
deallocate(dam)
deallocate(coefbt, cofbft)
if (allocated(wcf)) deallocate(wcf)
deallocate(w1, w2, w3)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1001 format (/,a8,' MAX= ',e11.4, ' EN ',e11.4,' ',e11.4,' ',e11.4)
 1002 format (  a8,' MIN= ',E11.4, ' EN ',e11.4,' ',e11.4,' ',e11.4)
 1003 format (/,'CLIPPINGS DE DT : ',                                  &
                             i10,' A ',e11.4,', ',i10,' A ',e11.4)

#else

 1001 format (/,a8,' MAX= ',e11.4, ' IN ',e11.4,' ',e11.4,' ',e11.4)
 1002 format (  a8,' MIN= ',e11.4, ' IN ',e11.4,' ',e11.4,' ',e11.4)
 1003 format (/,'DT CLIPPING : ',                                      &
                             i10,' A ',e11.4,', ',i10,' A ',e11.4)

#endif

!----
! End
!----

return

end subroutine
