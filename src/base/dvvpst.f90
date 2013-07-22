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

subroutine dvvpst &
!================

 ( nummai , numtyp ,                                              &
   nvar   , nscal  , nvlsta , nvisbr ,                            &
   ncelps , nfacps , nfbrps ,                                     &
   itypps ,                                                       &
   lstcel , lstfac , lstfbr ,                                     &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statce , stativ , statfb ,                   &
   tracel , trafbr )

!===============================================================================
! Purpose:
! --------

! Standard output of variables on post-processing meshes
!   (called after cs_user_extra_operations)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nummai           ! ec ! <-- ! numero du maillage post                        !
! numtyp           ! ec ! <-- ! numero de type de post-traitement              !
!                  !    !     ! (-1: volume, -2: bord, nummai par defaut)      !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nvlsta           ! e  ! <-- ! nombre de variables stat. lagrangien           !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! ncelps           ! e  ! <-- ! nombre de cellules du maillage post            !
! nfacps           ! e  ! <-- ! nombre de faces interieur post                 !
! nfbrps           ! e  ! <-- ! nombre de faces de bord post                   !
! itypps(3)        ! te ! <-- ! indicateur de presence (0 ou 1) de             !
!                  !    !     ! cellules (1), faces (2), ou faces de           !
!                  !    !     ! de bord (3) dans le maillage post              !
! lstcel(ncelps    ! te ! <-- ! liste des cellules du maillage post            !
! lstfac(nfacps    ! te ! <-- ! liste des faces interieures post               !
! lstfbr(nfbrps    ! te ! <-- ! liste des faces de bord post                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! statce           ! tr ! <-- ! statistiques cellules (lagrangien)             !
!  (ncelet,nvlsta) !    !     !                                                !
! stativ           ! tr ! <-- ! statistiques variance cellules (lagrangien)    !
!  (ncelet,nvlsta) !    !     !                                                !
! statfb           ! tr ! <-- ! statistiques faces bord (lagrangien)           !
!  (nfabor,nvisbr) !    !     !                                                !
! tracel(*)        ! tr ! <-- ! tab reel valeurs cellules post                 !
! trafbr(*)        ! tr ! <-- ! tab reel valeurs faces bord post               !
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
use pointe
use entsor
use cstnum
use cstphy
use optcal
use numvar
use parall
use period
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use radiat
use cplsat
use mesh
use field
use post

!===============================================================================

implicit none

! Arguments

integer          nummai , numtyp
integer          nvar   , nscal  , nvlsta , nvisbr
integer          ncelps , nfacps , nfbrps

integer          itypps(3)
integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)

double precision dt(ncelet), rtpa(ncelet,*), rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision statce(ncelet,nvlsta), statfb(nfabor,nvisbr)
double precision stativ(ncelet,nvlsta)
double precision tracel(ncelps*3)
double precision trafbr(nfbrps*3)

! Local variables

character*80     name80

logical          ilved , ientla, ivarpr
integer          inc   , iccocg, nswrgp, imligp, iwarnp
integer          ifac  , iloc  , ivar , iclvar, iclvaf
integer          ira   , idivdt
integer          ipp   , idimt , ii    , kk   , ll, iel
integer          ivarl , ivar0 , iip
integer          iii, ivarl1 , ivarlm , iflu   , ilpd1  , icla
integer          iscal , ipcvsl, ipcvst, iflmab
integer          fldid, fldprv, keycpl, iflcpl
integer          ipccp , ipcrom, ipcsii, keyvis, iflpst, itplus

double precision xcp   , xvsl  , srfbn
double precision visct , flumab, diipbx, diipby, diipbz
double precision epsrgp, climgp, extrap
double precision pcentr

double precision rbid(1)

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: treco
double precision, allocatable, dimension(:) :: wcell
double precision, dimension(:), pointer :: tplusp, tstarp
double precision, dimension(:), pointer :: valsp, coefap, coefbp
double precision, dimension(:,:), pointer :: valvp, cofavp, cofbvp
double precision, dimension(:,:,:), pointer :: cofbtp

!===============================================================================

! Initialize variables to avoid compiler warnings

ipp = 0

!===============================================================================
! 1.1. Fluid domain
!===============================================================================

if (numtyp .eq. -1) then

  !  1.1.2 Automatic additional variables
  !  ------------------------------------

  ! Wall distance (if LES+VanDriest or Rij+Echo or K-w SST)

  if (ineedy.eq.1 .and. abs(icdpar).eq.1) then

    idimt = 1
    ientla = .true.
    ivarpr = .true.

    call post_write_var(nummai, 'DistWall', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, dispar, rbid, rbid)

  endif

  ! Yplus (if LES+VanDriest)

  if (ineedy.eq.1 .and. abs(icdpar).eq.1) then

    if (itytur.eq.4.and.idries.eq.1) then

      idimt = 1
      ientla = .true.
      ivarpr = .true.

      call post_write_var(nummai, 'Yplus', idimt, ientla, ivarpr,  &
                          ntcabs, ttcabs, yplpar, rbid, rbid)

    endif

  endif

  ! Vitesse et pression absolues en cas de calcul en repère relatif

  if (icorio.eq.1) then

    ipcrom = ipproc(irom)

    idimt = 1
    ientla = .true.
    ivarpr = .false.

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      pcentr =   0.5d0*((omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))**2 &
                      + (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))**2 &
                      + (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))**2)

      tracel(iloc) = rtp(iel,ipr) + propce(iel,ipcrom)*pcentr

    enddo

    call post_write_var(nummai, 'Pressure', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, tracel, rbid, rbid)

    idimt = 3
    ientla = .true.
    ivarpr = .false.

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      tracel(1 + (iloc-1)*idimt) = rtp(iel,iu) &
          + (omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))

      tracel(2 + (iloc-1)*idimt) = rtp(iel,iv) &
          + (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))

      tracel(3 + (iloc-1)*idimt) = rtp(iel,iw) &
          + (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))

    enddo

    call post_write_var(nummai, 'Velocity', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, tracel, rbid, rbid)

  endif

  ! Vitesse et pression relatives en cas de calcul en repère fixe

  if (imobil.eq.1) then

    ipcrom = ipproc(irom)

    idimt = 1
    ientla = .true.
    ivarpr = .false.

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      pcentr =   0.5d0*((omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))**2 &
                      + (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))**2 &
                      + (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))**2)

      tracel(iloc) = rtp(iel,ipr) - propce(iel,ipcrom)*pcentr

    enddo

    call post_write_var(nummai, 'Rel Pressure', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, tracel, rbid, rbid)

    idimt = 3
    ientla = .true.
    ivarpr = .false.

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      tracel(1 + (iloc-1)*idimt) = rtp(iel,iu) &
          - (omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))

      tracel(2 + (iloc-1)*idimt) = rtp(iel,iv) &
          - (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))

      tracel(3 + (iloc-1)*idimt) = rtp(iel,iw) &
          - (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))

    enddo

    call post_write_var(nummai, 'Rel Velocity', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, tracel, rbid, rbid)

  endif


!===============================================================================
! 1.2. Boundary
!===============================================================================

else if (numtyp .eq. -2) then

  !  Projection of variables at boundary with no reconstruction
  !  ----------------------------------------------------------

  call field_get_key_id('post_vis', keyvis)
  call field_get_key_id('coupled', keycpl)

  fldprv = -1

  do ivar = 1, nvar  ! Loop on main cell-based variables

    fldid = ivarfl(ivar)

    if (fldid .eq. fldprv) cycle ! already output for multiple components

    fldprv = fldid

    call field_get_key_int(fldid, keyvis, iflpst)

    if (iand(iflpst, 2) .eq. 0) cycle ! nothing to do for this field

    call field_get_dim (fldid, idimt, ilved)
    call field_get_name(fldid, name80(4:80))
    name80(1:3) = 'bc_'


    !  Compute non-reconstructed values at boundary faces

    if (idimt.ne.1) then
      call field_get_key_int(fldid, keycpl, iflcpl)
    else
      iflcpl = 0
    endif

    if (idimt.eq.1) then  ! Scalar

      call field_get_val_s(fldid, valsp)
      call field_get_coefa_s(fldid, coefap)
      call field_get_coefb_s(fldid, coefbp)

      do iloc = 1, nfbrps

        ifac = lstfbr(iloc)
        iel = ifabor(ifac)

        trafbr(iloc) =   coefap(ifac) + coefbp(ifac)*valsp(iel)

      enddo

    else if (iflcpl.eq.0) then  ! Uncoupled vector or tensor

      call field_get_val_v(fldid, valvp)
      call field_get_coefa_v(fldid, cofavp)
      call field_get_coefb_uv(fldid, cofbvp)

      if (.not.ilved) then

        do kk = 0, idimt-1

          do iloc = 1, nfbrps

            ifac = lstfbr(iloc)
            iel = ifabor(ifac)

            trafbr(kk + (iloc-1)*idimt + 1)                      &
                 =   cofavp(ifac,kk+1)                           &
                   + cofbvp(ifac,kk+1)*valvp(iel,kk+1)

          enddo

        enddo

      else ! if interleaved

        do kk = 0, idimt-1

          do iloc = 1, nfbrps

            ifac = lstfbr(iloc)
            iel = ifabor(ifac)

            trafbr(kk + (iloc-1)*idimt + 1)                      &
                 =   cofavp(kk+1,ifac)                           &
                   + cofbvp(kk+1,ifac)*valvp(kk+1,iel)

          enddo

        enddo

      endif

    else ! Coupled vector or tensor

      call field_get_val_v(fldid, valvp)
      call field_get_coefa_v(fldid, cofavp)
      call field_get_coefb_v(fldid, cofbtp)

      if (.not.ilved) then ! in coupled case coefa/coefb interleaved

        do kk = 0, idimt-1

          do iloc = 1, nfbrps

            ifac = lstfbr(iloc)
            iel = ifabor(ifac)

            trafbr(kk + (iloc-1)*idimt + 1) = cofavp(kk+1,ifac)

            do ll = 1, idimt
              trafbr(kk + (iloc-1)*idimt + 1)                    &
                 =   trafbr(kk + (iloc-1)*idimt + 1)             &
                   + cofbtp(kk+1,ll,ifac)*valvp(iel,ll)
            enddo

          enddo

        enddo

      else ! coupled + interleaved case

        do kk = 0, idimt-1

          do iloc = 1, nfbrps

            ifac = lstfbr(iloc)
            iel = ifabor(ifac)

            trafbr(kk + (iloc-1)*idimt + 1) = cofavp(kk+1,ifac)

            do ll = 1, idimt
              trafbr(kk + (iloc-1)*idimt + 1)                    &
                 =   trafbr(kk + (iloc-1)*idimt + 1)             &
                   + cofbtp(kk+1,ll,ifac)*valvp(ll,iel)
            enddo

          enddo

        enddo

      endif

    endif ! test on field dimension and interleaving

    ientla = .true.  ! interleaved result values
    ivarpr = .false. ! defined on work array

    call post_write_var(nummai, trim(name80), idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, rbid, rbid, trafbr)

  enddo ! End of loop on variables

  ! output y+ at the boundary
  ! -------------------------

  if (ipstdv(ipstyp).ne.0) then

    idimt = 1  ! variable dimension

    ! Compute variable on boundary faces

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      trafbr(1 + (iloc-1)*idimt) = yplbr(ifac)
    enddo

    ! Non interleaved values, defined in work array

    ientla = .true.
    ivarpr = .false.

    call post_write_var(nummai, 'Yplus', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, rbid, rbid, trafbr)

  endif ! end of test on output of y+

  ! Handle efforts at boundary
  ! --------------------------------

  if (iand(ipstdv(ipstfo), 1) .ne. 0) then

    ! Compute variable values on boundary faces

    call post_efforts(nfbrps, lstfbr, trafbr)

    idimt = 3        ! variable dimension
    ientla = .true.  ! interleaved values
    ivarpr = .false. ! defined on work array

    call post_write_var(nummai, 'Efforts', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, rbid, rbid, trafbr)

  endif

  if (iand(ipstdv(ipstfo), 2) .ne. 0) then

    ! Compute variable values on boundary faces

    call post_efforts_tangential(nfbrps, lstfbr, trafbr)

    idimt = 3        ! variable dimension
    ientla = .true.  ! interleaved values
    ivarpr = .false. ! defined on work array

    call post_write_var(nummai, 'Tangential Efforts', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, rbid, rbid, trafbr)

  endif

  if (iand(ipstdv(ipstfo), 4) .ne. 0) then

    ! Calcul des valeurs de la variable sur les faces de bord

    call post_efforts_normal(nfbrps, lstfbr, trafbr)

    idimt = 1        ! variable dimension
    ientla = .true.  ! interleaved values
    ivarpr = .false. ! defined on work array

    call post_write_var(nummai, 'Normal Efforts', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, rbid, rbid, trafbr)

  endif

  ! T+ near the boundary
  ! --------------------

  if (ipstdv(ipsttp).ne.0) then

    call field_get_id_try('tplus', itplus)

    if (itplus.ge.0) then

      call field_get_val_s(itplus, tplusp)

      idimt = 1        ! variable dimension
      ientla = .true.  ! interleaved values
      ivarpr = .true.  ! defined on parent array

      call post_write_var(nummai, 'Tplus', idimt, ientla, ivarpr,  &
                          ntcabs, ttcabs, rbid, rbid, tplusp)

    endif ! end of test on presence ot T+

  endif ! end of test on output of y+

  ! Thermal flux at boundary
  ! ------------------------
  !  If working with enthalpy, compute an enthalpy flux

  if (ipstdv(ipstft).ne.0) then

    if (iscalt.gt.0 .and. nscal.gt.0 .and. iscalt.le.nscal) then

      call post_boundary_thermal_flux(nfbrps, lstfbr, rtp, propce, propfb,      &
                                      trafbr)

      idimt = 1        ! variable dimension
      ientla = .true.  ! interleaved values
      ivarpr = .false. ! defined on work array

      call post_write_var(nummai, 'Input thermal flux', idimt, ientla, ivarpr,  &
                          ntcabs, ttcabs, rbid, rbid, trafbr)

    endif

  endif

  ! Temperature at the boundary
  ! ---------------------------

  if (ipstdv(ipsttb).ne.0) then

    idimt = 1        ! variable dimension
    ientla = .true.  ! interleaved values
    ivarpr = .false. ! defined on work array

    ! Compute variable on boundary faces

    call post_boundary_temperature(nfbrps, lstfbr, rtp, propce, propfb,     &
                                   trafbr)

    call post_write_var(nummai, 'Wall temperature', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, rbid, rbid, trafbr)

  endif ! end of test on output of wall temperature

  ! Nusselt at the boundary
  ! -----------------------

  if (ipstdv(ipstnu).ne.0) then

    idimt = 1        ! variable dimension
    ientla = .true.  ! interleaved values
    ivarpr = .false. ! defined on work array

    ! Compute variable on boundary faces

    call post_boundary_nusselt(nfbrps, lstfbr, rtp, propce, propfb,         &
                               trafbr)

    call post_write_var(nummai, 'Wall law Nusselt', idimt, ientla, ivarpr,  &
                        ntcabs, ttcabs, rbid, rbid, trafbr)

  endif ! end of test on output of Nusselt

endif ! end of test on postprocessing mesh number

!===============================================================================
! 2.1. Lagrangian variables
!===============================================================================

if (numtyp .eq. -1) then

  if (iilagr.gt.0 .and. istala.ge.1) then

    ! All standard statistics have dimension 1, and are defined or computed
    ! on the global mesh cells.

    idimt  = 1
    ientla = .true.
    ivarpr = .true.

    allocate(wcell(ncelet))

    iii = nvlsta-nvlsts

    do icla  = 0, nbclst

      ! -> if ICLA = 0: global statistics
      !    if 0 < ICLA =< NBCLST: per group statistics

      do ivarl = 1, nvlsta

        ivarl1 = icla*nvlsta +ivarl
        ivarlm = ivarl1
        ilpd1  = icla*nvlsta +ilpd
        iflu   = 0

        if (ivarl.le.iii) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlag(ivarl)
          else
            write(name80,'(a8,a4,i3)') nomlag(ivarl),'_grp',icla
          endif
        else if (nvlsts.gt.0) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlag(ilvu(ivarl-iii))
          else
            write(name80,'(a8,a4,i3)')                            &
                  nomlag(ilvu(ivarl-iii)),'_grp',icla
          endif
        endif

        call uslaen                                               &
        !==========
          (nvar, nscal, nvlsta,                                   &
           ivarl, ivarl1, ivarlm, iflu, ilpd1, icla,              &
           dt, rtpa, rtp, propce, propfa, propfb,                 &
           statce, stativ, wcell)

        call post_write_var(nummai, trim(name80), idimt, ientla, ivarpr,  &
                            ntcabs, ttcabs, wcell, rbid, rbid)

      enddo

      do ivarl = 1, nvlsta-1

        ivarl1 = icla*(nvlsta-1)+ivarl
        ivarlm = icla*nvlsta+ivarl
        ilpd1  = icla*nvlsta +ilpd
        iflu   = 1

        if (ivarl.le.iii) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlav(ivarl)
          else
            write(name80,'(a8,a4,i3)') nomlav(ivarl),'_grp',icla
          endif
        else if (nvlsts.gt.0) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlav(ilvu(ivarl-iii))
          else
            write(name80,'(a8,a4,i3)')                            &
                 nomlav(ilvu(ivarl-iii)),'_grp',icla
          endif
        endif

        call uslaen                                               &
        !==========
          (nvar, nscal, nvlsta,                                   &
           ivarl, ivarl1, ivarlm, iflu, ilpd1, icla,              &
           dt, rtpa, rtp, propce, propfa, propfb,                 &
           statce, stativ, wcell)

        call post_write_var(nummai, trim(name80), idimt, ientla, ivarpr,  &
                            ntcabs, ttcabs, wcell, rbid, rbid)
      enddo

    enddo

    deallocate(wcell)

  endif

endif

if (numtyp.eq.-2) then

  if (iilagr.gt.0 .and. iensi3.eq.1) then

    iii = nvisbr-nusbor

    do ivarl = 1,nvisbr

      if (ivarl.le.iii) then
        name80 = nombrd(ivarl)
      else if (nusbor.gt.0) then
        name80 = nombrd(iusb(ivarl-iii))
      endif

      if (imoybr(ivarl).eq.3) then

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          if (statfb(ifac,iencnb).gt.seuilf) then
            trafbr(iloc) = statfb(ifac,ivarl)/statfb(ifac,iencnb)
          else
            trafbr(iloc) = 0.d0
          endif
        enddo

      else if (imoybr(ivarl).eq.2) then

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          if (statfb(ifac,inbr).gt.seuilf) then
            trafbr(iloc) = statfb(ifac,ivarl)/statfb(ifac,inbr)
          else
            trafbr(iloc) = 0.d0
          endif
        enddo

      else if (imoybr(ivarl).eq.1) then

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          if (statfb(ifac,inbr).gt.seuilf) then
            trafbr(iloc) = statfb(ifac,ivarl) / tstatp
          else
            trafbr(iloc) = 0.d0
          endif
        enddo

      else

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          if (statfb(ifac,inbr).gt.seuilf) then
            trafbr(iloc) = statfb(ifac,ivarl)
          else
            trafbr(iloc) = 0.d0
          endif
        enddo

      endif

      idimt  = 1
      ientla = .true.
      ivarpr = .false.

      call post_write_var(nummai, trim(name80), idimt, ientla, ivarpr,  &
                          ntcabs, ttcabs, rbid, rbid, trafbr)

    enddo

    !do iloc = 1, nfbrps
    !  ifac = lstfbr(iloc)
    !  trafbr(iloc) = ia(iifrla+ifac-1) !! TODO: ifrlag (cf caltri)
    !enddo

    !idimt  = 1
    !ientla = .true.
    !ivarpr = .false.

    !call post_write_var(nummai, 'lagrangian_boundary_zones', idimt,          &
    !                    ientla, ivarpr, ntcabs, ttcabs, rbid, rbid, trafbr)

  endif
endif
!     Fin du test sur le numero de maillage post.

!===============================================================================
!     2.2. VARIABLES RADIATIVES AUX FRONTIERES
!===============================================================================

if (numtyp.eq.-2) then

  if (iirayo.gt.0) then

    do ivarl = 1,nbrayf

      if (irayvf(ivarl).eq.1) then

        if (ivarl .eq. itparp)      then
          ipp =  ipprob(itparo)
        else if (ivarl .eq. iqincp) then
          ipp = ipprob(iqinci)
        else if (ivarl .eq. ixlamp)  then
          ipp = ipprob(ixlam)
        else if (ivarl .eq. iepap)   then
          ipp = ipprob(iepa)
        else if (ivarl .eq. iepsp)   then
          ipp = ipprob(ieps)
        else if (ivarl .eq. ifnetp)  then
          ipp = ipprob(ifnet)
        else if (ivarl .eq. ifconp) then
          ipp = ipprob(ifconv)
        else if (ivarl .eq. ihconp) then
          ipp = ipprob(ihconv)
        endif

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          trafbr(iloc) = propfb(ifac,ipp)
        enddo

        idimt  = 1
        ientla = .true.
        ivarpr = .false.

        call post_write_var(nummai, trim(nbrvaf(ivarl)), idimt,              &
                            ientla, ivarpr, ntcabs, ttcabs, rbid, rbid, trafbr)

      endif
    enddo

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      trafbr(iloc) = izfrad(ifac)
    enddo

    idimt  = 1
    ientla = .true.
    ivarpr = .false.

    call post_write_var(nummai, 'radiative_boundary_zones', idimt,           &
                        ientla, ivarpr, ntcabs, ttcabs, rbid, rbid, trafbr)

  endif
endif

!===============================================================================
! 2.3. Electric module variables
!===============================================================================

if (numtyp.eq.-1) then

  if (     ippmod(ieljou).ge.1                                      &
      .or. ippmod(ielarc).ge.1                                      &
      .or. ippmod(ielion).ge.1) then

    allocate(grad(ncelet,3))

    if (.true.) then

      ! Gradient of the real potential

      ivar = isca(ipotr)

      inc = 1
      iccocg = 1
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)
      ivar0 = 0

      call field_get_coefa_s(ivarfl(ivar), coefap)
      call field_get_coefb_s(ivarfl(ivar), coefbp)

      call grdcel                                                 &
      !==========
        (ivar0, imrgra, inc, iccocg, nswrgp, imligp,              &
         iwarnp, nfecra, epsrgp, climgp, extrap,                  &
         rtp(1,ivar), coefap, coefbp,                             &
         grad)

      idimt  = 3
      ientla = .false.
      ivarpr = .true.

      call post_write_var(nummai, 'Pot_Gradient_R', idimt, ientla, ivarpr,   &
                          ntcabs, ttcabs, grad, rbid, rbid)

    endif

    ! For Joule Heating by direct conduction:
    !   gradient of the imaginary component of the potential

    if (.true.                                                               &
        .and. (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4)) then

      ivar = isca(ipoti)

      inc = 1
      iccocg = 1
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)

      ivar0 = 0

      call field_get_coefa_s(ivarfl(ivar), coefap)
      call field_get_coefb_s(ivarfl(ivar), coefbp)

      call grdcel                                                 &
      !==========
        (ivar0, imrgra, inc, iccocg, nswrgp, imligp,              &
         iwarnp, nfecra, epsrgp, climgp, extrap,                  &
         rtp(1,ivar), coefap, coefbp,                             &
         grad)

      idimt  = 3
      ientla = .false.
      ivarpr = .true.

      call post_write_var(nummai, 'Pot_Gradient_Im', idimt, ientla, ivarpr,  &
                          ntcabs, ttcabs, grad, rbid, rbid)

    endif

    ! For Joule heating by direct conduction:
    !   imaginary component of the current density

    if (.true.                                                               &
        .and. (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4)) then

      ivar = isca(ipoti)

      ! As in elflux
      ipcsii = ipproc(ivisls(ipoti))

      inc = 1
      iccocg = 1
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)

      ivar0 = 0

      call field_get_coefa_s(ivarfl(ivar), coefap)
      call field_get_coefb_s(ivarfl(ivar), coefbp)

      call grdcel                                                 &
      !==========
        (ivar0, imrgra, inc, iccocg, nswrgp, imligp,              &
         iwarnp, nfecra, epsrgp, climgp, extrap,                  &
         rtp(1,ivar), coefap, coefbp,                             &
         grad)

      do iloc = 1, ncelps
        iel = lstcel(iloc)
        tracel(1 + (iloc-1)*idimt) = -propce(iel,ipcsii)*grad(iel,1)
        tracel(2 + (iloc-1)*idimt) = -propce(iel,ipcsii)*grad(iel,2)
        tracel(3 + (iloc-1)*idimt) = -propce(iel,ipcsii)*grad(iel,3)
      enddo

      idimt  = 3
      ientla = .true.
      ivarpr = .false.

      call post_write_var(nummai, 'Current_Im', idimt, ientla, ivarpr,  &
                          ntcabs, ttcabs, grad, rbid, rbid)

    endif

    ! For electric arcs: electromagnetic field calculation

    if (.true. .and. ippmod(ielarc).ge.2) then

      ! Ax Component

      ivar = isca(ipotva(1))

      inc = 1
      iccocg = 1
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)

      ivar0 = 0

      call field_get_coefa_s(ivarfl(ivar), coefap)
      call field_get_coefb_s(ivarfl(ivar), coefbp)

      call grdcel                                                 &
      !==========
        (ivar0, imrgra, inc, iccocg, nswrgp, imligp,              &
         iwarnp, nfecra, epsrgp, climgp, extrap,                  &
         rtp(1,ivar), coefap, coefbp,                             &
         grad)

      ! B = rot A ( B = curl A)

      do iloc = 1, ncelps
        iel = lstcel(iloc)
        tracel(1 + (iloc-1)*idimt) =  zero
        tracel(2 + (iloc-1)*idimt) =  grad(iel,3)
        tracel(3 + (iloc-1)*idimt) = -grad(iel,2)
      enddo

      ! Ay component

      ivar = isca(ipotva(2))

      inc = 1
      iccocg = 1
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)

      ivar0 = 0

      call field_get_coefa_s(ivarfl(ivar), coefap)
      call field_get_coefb_s(ivarfl(ivar), coefbp)

      call grdcel                                                 &
      !==========
        (ivar0, imrgra, inc, iccocg, nswrgp, imligp,              &
         iwarnp, nfecra, epsrgp, climgp, extrap,                  &
         rtp(1,ivar), coefap, coefbp,                             &
         grad)

      ! B = rot A (B = curl A)

      do iloc = 1, ncelps
        iel = lstcel(iloc)
        tracel(1 + (iloc-1)*idimt) = tracel(1 + (iloc-1)*idimt) - grad(iel,3)
        tracel(3 + (iloc-1)*idimt) = tracel(3 + (iloc-1)*idimt) + grad(iel,1)
      enddo

      ! Az component

      ivar = isca(ipotva(3))

      inc = 1
      iccocg = 1
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)

      ivar0 = 0

      call field_get_coefa_s(ivarfl(ivar), coefap)
      call field_get_coefb_s(ivarfl(ivar), coefbp)

      call grdcel                                                 &
      !==========
        (ivar0, imrgra, inc, iccocg, nswrgp, imligp,              &
         iwarnp, nfecra, epsrgp, climgp, extrap,                  &
         rtp(1,ivar), coefap, coefbp,                             &
         grad   )

      ! B = rot A (B = curl A)

      do iloc = 1, ncelps
        iel = lstcel(iloc)
        tracel(1 + (iloc-1)*idimt) = tracel(1 + (iloc-1)*idimt) + grad(iel,2)
        tracel(2 + (iloc-1)*idimt) = tracel(2 + (iloc-1)*idimt) - grad(iel,1)
      enddo

      idimt  = 3
      ientla = .true.
      ivarpr = .false.

      call post_write_var(nummai, 'Magnetic_field', idimt, ientla, ivarpr,  &
                          ntcabs, ttcabs, tracel, rbid, rbid)

    endif

    ! Calculation of Module and Argument of the complex potential if IELJOU = 4

    if (.true. .and. ippmod(ieljou).eq.4) then

      ivar = isca(ipotr)

      do iloc = 1, ncelps
        iel = lstcel(iloc)
        tracel(iloc) =                                              &
          sqrt( rtp(iel,isca(ipotr))*rtp(iel,isca(ipotr))           &
               +rtp(iel,isca(ipoti))*rtp(iel,isca(ipoti)) )
      enddo

      idimt  = 1
      ientla = .true.
      ivarpr = .false.

      call post_write_var(nummai, 'Pot_Module', idimt, ientla, ivarpr,  &
                          ntcabs, ttcabs, tracel, rbid, rbid)

      ivar = isca(ipotr)

      do iloc = 1, ncelps

        iel = lstcel(iloc)

        if (rtp(iel,isca(ipotr)) .ne. 0.d0) then
          if (rtp(iel,isca(ipotr)) .ge. 0.d0) then
            tracel(iloc) = atan(rtp(iel,isca(ipoti))/rtp(iel,isca(ipotr)))
          else
            if (rtp(iel,isca(ipoti)) .gt. 0.d0) then
              tracel(iloc) = 4.d0*atan(1.d0)                      &
                             + atan(  rtp(iel,isca(ipoti))        &
                                    / rtp(iel,isca(ipotr)))
            else
              tracel(iloc) = -4.d0*atan(1.d0)                     &
                             + atan(  rtp(iel,isca(ipoti))        &
                                    / rtp(iel,isca(ipotr)))
            endif
          endif
        else
          tracel(iloc) = 2.d0*atan(1.d0)
        endif

        if (tracel(iloc) .lt. 0.d0) then
          tracel(iloc) = tracel(iloc) + 8.d0**atan(1.d0)
        endif

      enddo

      idimt  = 1
      ientla = .true.
      ivarpr = .false.

      call post_write_var(nummai, 'Pot_Arg', idimt, ientla, ivarpr,  &
                          ntcabs, ttcabs, tracel, rbid, rbid)

    endif

    ! Free memory
    deallocate(grad)

  endif

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
