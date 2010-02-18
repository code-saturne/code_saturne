!-------------------------------------------------------------------------------

!VERS


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

subroutine ustsma &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ncesmp , iphas  , iappel ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Mass source term for a phase iphas

! The subroutine ustsma is called at three different stages in the code
!  (iappel = 1, 2 or 3)

! iappel = 1
!    Calculation of the number of cells where a mass source term is
!    imposed: ncesmp
!    Called once at the beginnign of the calculation

! iappel = 2
!    Identification of the cells where a mass source term is imposed:
!    array icesmp(ncesmp)
!    Called once at the beginnign of the calculation

! iappel = 3
!    Calculation of the values of the mass source term
!    Called at each time step



! The equation for mass conservation becomes

!           d(rho)/dt + div(rho u) = gamma

! The equation for a variable f becomes

!           d(f)/dt = ..... + gamma*(f_i - f)

!   discretized as

!           rho*(f^(n+1) - f^(n))/dt = .....
!                                    + gamma*(f_i - f^(n+1))

! f_i is the value of f associated to the injecte mass.
! Two options are available:
!   - the mass flux is injected with the local value of variable f
!                           --> f_i = f^(n+1)
!                   (the equation for f is therefore not modified)
!
!   - the mass flux is injected with a specific value for f
!                           --> f_i is specified by the user


! Variables to be specified by the user
! =====================================

!  ncesmp: number of cells where a mass source term is imposed

!  icetsm(ieltsm): identification of the cells where a mass source
!                  term is imposed.
!                  For each cell where a mass source term is imposed
!                  (ielstm in [1;ncesmp]), icetsm(ieltsm) is the
!                  global index number of the corresponding cell
!                  (icestm(ieltsm) in [1;ncel])

!  smacel(ieltsm,ipr(iphas)): value of the injection mass rate gamma (kg/m3/s)
!                             in the ieltsm cell with mass source term

!  itypsm(ieltsm,ivar): type of treatment for variable ivar in the
!                       ieltsm cell with mass source term.
!                     * itypsm = 0 --> injection of ivar at local value
!                     * itypsm = 1 --> injection of ivar at user
!                                      specified value

!  smacel(ieltsm,ivar): specified value for variable ivar associated
!                       to the injected mass in the ieltsm cell with
!                       a mass source term
!                                  except for ivar=ipr(iphas)

!
! Remarks
! =======
!
! - if itypsm(ieltsm,ivar)=0, smacel(ieltsm,ivar) is not used

! - if smacel(ieltsm,ipr(iphas))<0, mass is removed from the system,
!     therefore Code_Saturna automatically considers f_i=f^(n+1),
!     whatever the values of itypsm or smacel specified by the user

! - if a value ivar is not linked to phase iphas for a mass source
!     term is imposed, no source term will be taen into account.

! - if a scalar doesn't evolve following the standard equation
!     d(rho f)/dt + d(rho U f)/dx = ...
!     (alternate convective field for instance), the source term
!     set by this routine will nto be correct (except in case of
!     injection at the local value of the variable). The proper source
!     term should be added directly in ustssc.


! Identification of cells
! =======================
! The selection of cells where to apply the source terms is based on a getcel
! command. For more info on the syntax of the getcel command, refer to the
! user manual or to the comments on the similar command getfbr in the routine
! usclim.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss terms           !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ncssmp           ! i  ! <-- ! number of cells with mass source terms         !
! iphas            ! i  ! <-- ! index number of the current phase              !
! iappel           ! i  ! <-- ! indicates which at which stage the routine is  !
!                  !    !     !  is called                                     !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! maxelt           ! i  ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see uttsma.f90)                              !
! idevel(nideve)   ! ia ! <-- ! integer work array for temporary developpement !
! ituser(nituse    ! ia ! <-- ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfavor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time steps)                        !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate                     !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary developpement    !
! rtuser(nituse    ! ra ! <-- ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "pointe.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "parall.h"
include "period.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas  , iappel

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6)
double precision smacel(ncesmp,nvar)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ieltsm
integer          ifac, iutile, ii
integer          ilelt, nlelt

double precision vent, vent2
double precision dh, ustar2
double precision xkent, xeent
double precision flucel
double precision vtot  , gamma

!===============================================================================

idebia = idbia0
idebra = idbra0

if(iappel.eq.1.or.iappel.eq.2) then

!===============================================================================
! 1. For each phase, one or two calls

!   First call:
!
!       iappel = 1: ncesmp: calculation of the number of cells with
!                             mass source term


!   Second call (if ncesmp>0):
!       iappel = 2: icetsm: index number of cells with mass source terms

! WARNINGS
! ========
!   Do not use smacel in this section (it is set on the third call, iappel=3)

!   Do not use icetsm in this section on the first call (iappel=1)

!   This section (iappel=1 or 2) is only accessed at the beginning of a
!     calculation. Should the localization of the mass source terms evolve
!     in time, the user must identify at the beginning all cells that can
!     potentially becomea mass source term.

!===============================================================================


!  1.1 To be completed by the user: cell selection
!  -----------------------------------------------

! Example 1: No mass source term (default)
  ieltsm = 0


! Example 2 : Mass source term for phase one in the cells that
!              have a boundary face of color 3 and the cells
!              with a coordinate X between 2.5 and 5.
!
!     In this test in two parts, one mut pay attention not to count
!      the cells twice (a cell with a boundary face of color 3 can
!      also have a coordinate X between 2.5 and 5).
!     One should also pay attention that, on the first call, the
!      array icetsm doesn't exist yet. It mustn't be used outside
!      of tests (iappel.eq.2).

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

  iutile = 0
  if(iutile.eq.1) then

    if(iphas.eq.1) then

      ieltsm = 0

!     Cells with coordinate X between 2.5 and 5.
      CALL GETCEL('X > 2.5 and X < 5.0',NLELT,LSTELT)
      do ilelt = 1, nlelt
        ii = lstelt(ilelt)
        ieltsm = ieltsm + 1
        if (iappel.eq.2) icetsm(ieltsm) = ii
      enddo


!     Cells with a boundary face of color 3

      CALL GETFBR('3',NLELT,LSTELT)
      do ilelt = 1, nlelt
        ifac = lstelt(ilelt)
        ii   = ifabor(ifac)
!       The cells that have already been counted above are not
!        counted again.
        if (.not.(xyzcen(1,ii).lt.500.d0.and.                     &
                  xyzcen(1,ii).gt.250.d0)    )then
          ieltsm = ieltsm + 1
          if (iappel.eq.2) icetsm(ieltsm) = ii
        endif
      enddo

    else
      ieltsm = 0
    endif

 endif


!  1.2 Generic subsection: do not modify
!  -------------------------------------

! --- For iappel = 1,
!      Specification of ncesmp. This block is valid for both examples.

  if (iappel.eq.1) then
    ncesmp = ieltsm
  endif

!-------------------------------------------------------------------------------

elseif(iappel.eq.3) then

!===============================================================================

! 2. For each phase with ncesmp > 0 , third call

!       iappel = 3 : itypsm : type of mass source term
!                    smacel : mass source term


! Remark
! ======
! If itypsm(ieltsm,ivar) is set to 1, smacel(ieltsm,ivar) must be set.

!===============================================================================



!  2.1 To be completed by the user: itypsm and smacel
!  --------------------------------------------------

! Example 1: simulation of an inlet condition by mass source terms
!            and printing of the total mass rate for phase 1.

  if(iphas.eq.1) then

    vent = 0.1d0
    vent2 = vent**2
    dh     = 0.5d0
!
! Calculation of the inlet conditions for k and epsilon with standard
!   laws in a circular pipe.
    ustar2 = 0.d0
    xkent  = epzero
    xeent  = epzero

    call keendb                                                   &
    !==========
      ( vent2, dh, ro0(iphas), viscl0(iphas), cmu, xkappa,        &
        ustar2, xkent, xeent )

    flucel = 0.d0
    do ieltsm = 1, ncesmp
      smacel(ieltsm,ipr(iphas)) = 30000.d0
      itypsm(ieltsm,iv(iphas)) = 1
      smacel(ieltsm,iv(iphas)) = vent
      if (itytur(iphas).eq.2) then
        itypsm(ieltsm,ik(iphas)) = 1
        smacel(ieltsm,ik(iphas)) = xkent
        itypsm(ieltsm,iep(iphas)) = 1
        smacel(ieltsm,iep(iphas)) = xeent
      else if (itytur(iphas).eq.3) then
        itypsm(ieltsm,ir11(iphas)) = 1
        itypsm(ieltsm,ir12(iphas)) = 1
        itypsm(ieltsm,ir13(iphas)) = 1
        itypsm(ieltsm,ir22(iphas)) = 1
        itypsm(ieltsm,ir23(iphas)) = 1
        itypsm(ieltsm,ir33(iphas)) = 1
        smacel(ieltsm,ir11(iphas)) = 2.d0/3.d0*xkent
        smacel(ieltsm,ir12(iphas)) = 0.d0
        smacel(ieltsm,ir13(iphas)) = 0.d0
        smacel(ieltsm,ir22(iphas)) = 2.d0/3.d0*xkent
        smacel(ieltsm,ir23(iphas)) = 0.d0
        smacel(ieltsm,ir33(iphas)) = 2.d0/3.d0*xkent
        itypsm(ieltsm,iep(iphas)) = 1
        smacel(ieltsm,iep(iphas)) = xeent
      else if (iturb(iphas).eq.50) then
        itypsm(ieltsm,ik(iphas)) = 1
        smacel(ieltsm,ik(iphas)) = xkent
        itypsm(ieltsm,iep(iphas)) = 1
        smacel(ieltsm,iep(iphas)) = xeent
        itypsm(ieltsm,iphi(iphas)) = 1
        smacel(ieltsm,iphi(iphas)) = 2.d0/3.d0
! There is no mass source term in the equation for f_bar
      else if (iturb(iphas).eq.60) then
        itypsm(ieltsm,ik(iphas)) = 1
        smacel(ieltsm,ik(iphas)) = xkent
        itypsm(ieltsm,iomg(iphas))= 1
        smacel(ieltsm,iomg(iphas))= xeent/cmu/xkent
      endif
      if(nscal.gt.0) then
        do ii = 1, nscal
          if(iphsca(ii).eq.iphas) then
            itypsm(ieltsm,isca(ii)) = 1
            smacel(ieltsm,isca(ii)) = 1.d0
          endif
        enddo
      endif
      flucel = flucel+                                            &
                volume(icetsm(ieltsm))*smacel(ieltsm,ipr(iphas))
    enddo

    if (irangp.ge.0) then
      call parsom (flucel)
    endif

    if (iwarni(ipr(iphas)).ge.1) then
      write(nfecra,1000) iphas, flucel
    endif

  endif

!-------------------------------------------------------------------------------

! Example 2 : simulation of a suction (by a pump for instance) with a
!             total rate of 80 000 kg/s.
!             The suction rate is supposed to be uniformly distributed
!             on all the cells selected above.

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

  iutile = 0
  if(iutile.eq.1) then

    if(iphas.eq.1) then

! Calculation of the total volume of the area where the mass source
!   term is imposed (the case of parallel computing is taken into
!   account with the call to parsom).
      vtot = 0.d0
      do ieltsm = 1, ncesmp
        vtot = vtot + volume(icetsm(ieltsm))
      enddo
      if (irangp.ge.0) then
        call parsom (vtot)
      endif

! The mass suction rate is gamma = -80000/vtot (in kg/m3/s)
! It is set below, with a test for cases where vtot=0. The total
! mass rate is calculated for verification.

      if (vtot.gt.0.d0) then
        gamma = -80000.d0/vtot
      else
        write(nfecra,9000) iphas, vtot
        call csexit (1)
      endif

      flucel = 0.d0
      do ieltsm = 1, ncesmp
        smacel(ieltsm,ipr(iphas)) = gamma
        flucel = flucel+                                          &
                volume(icetsm(ieltsm))*smacel(ieltsm,ipr(iphas))
      enddo

      if (irangp.ge.0) then
        call parsom (flucel)
      endif

      if (iwarni(ipr(iphas)).ge.1) then
        write(nfecra,2000) iphas, flucel, vtot
      endif

    endif

  endif

!-------------------------------------------------------------------------------

endif

!--------
! Formats
!--------

 1000 format(/,'PHASE ',I3,                                             &
         ' : mass rate generated in the domain: ',E14.5,/)

 2000 format(/,'PHASE ',I3,                                             &
         ' : mass flux rate generated in the domain: ',E14.5,/,         &
'                         distributed on the volume: ',E14.5)

 9000 format(/,'PHASE ',I3,                                             &
         ' : error in ustsma                ',/,                        &
'   the volume of the mass suction area is = ',E14.5,/)

!----
! End
!----

return

end subroutine
