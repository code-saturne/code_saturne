!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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

subroutine uslabo &
!================

 ( lndnod ,                                                       &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   kface  , nbpt   , isuivi ,                                     &
   itypfb , itrifb , ifrlag , itepa  , indep  ,                   &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , parbor , vitpar , vitflu , auxl   )

!===============================================================================
! Purpose:
! --------
!
! User subroutine of the Lagrangian particle-tracking module:
! -----------------------------------------------------------
!
! User subroutine (non-mandatory intervention)
!
! User subroutine managing the particle behavior during an interaction between
! a particle and a boundary and recording the boundary statistics.
!
! The user does not need to modify this subroutine in the case of a standard use.
! If he wishes to treat non-standard user-defined interactions, he needs to intervene
! in sections 8 and 10.
!
! The interaction between a particle and a boundary face is treated with respect
! to the information given by the user (value of iusclb per zone) in the subroutine uslag2.
!
! Given the name stored in iusclb and associated to the boundary face kface, the type
! of particle behavior is defined. For a standard use, the value of iusclb can be either:
!
! * ientrl: for a zone where particles are injected into the domain (particle-inlet zone).
! * isortl: for a zone where particle are getting out of the domain (particle-outlet zone).
! * irebol: condition of elastic rebound.
! * idepo1: definitive deposition of the particles; the particle is removed from the calculation
! * idepo2: definitive deposition of the particles; the particle is kept in the calculation
!           (useful only if iensi2 = 1)
! * idepo3: deposition and resuspension possible depending on the flow conditions.
! * iencrl: fouling only for coal particles (iphyla = 2)
!
! Besides, if one wishes to add another kind of non-standard interaction for a zone of
! boundary faces, one must give (in uslag2) in iusclb(kzone) one of the following names:
!

!       JBORD1, JBORD2, JBORD3, JBORD4, JBORD5
!
! And, in the present routine uslabo, the user has to program the behavior of the particles
! for this boundary zone.
!
! CAUTION: At the beginning of the routine, the variable isuivi is initialized with an
! absurd value and MUST be modified before the end of the routine.
!
! The velocities of the the particle and the flow seen must be modified with respect to
! the interactions through the use of the arrays vitpar and vitflu, and MUST NOT be modified
! directly in the ETTP and ETTPA arrays in this routine.
!
! Rule to modify the isuivi parameter:
! ====================================
!
! 1) Set isuivi to 0 if the particle must not be followed in the mesh after its
!    interaction with a boundary face (ex: ientrl, isortl, idepo1, idepo2)

! 2) Set isuivi to 1 if the particle must be followed in the mesh after its
!    interaction with a boundary face (ex: idepo3)

! Remark: During an interaction, the computations of the velocities of the particle
! ------  and of the flow seen are first-order (even if the calculation is second-order
!         elsewhere in the domain)



!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! ntersl           ! i  ! <-- ! number of source terms of return coupling      !
! nvlsta           ! i  ! <-- ! nb of Lagrangian statistical variables         !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! kface            ! i  ! <-- ! number of the interaction face                 !
! nbpt             ! i  ! <-- ! number of the treated particle                 !
! isuivi           ! i  ! <-- ! flag to follow (or not) the particle           !
! itypfb(nfabor)   ! ia ! <-- ! type of the boundary faces                     !
! itrifb(nfabor)   ! ia ! <-- ! indirection array for the sorting of the faces !
! ifrlag           ! ia ! <-- ! number of the boundary face                    !
!   (nfabor)       !    !     ! for the Lagrangian module                      !
! itepa            ! ra ! <-- ! particle information (integer)                 !
! (nbpmax,nivep    !    !     !                                                !
! indep            ! ia ! --> ! for each cell, number of the departure cell    !
!   (nbpmax)       !    !     !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! transported variables at cell centers          !
! (ncelet,*)       !    !     ! (current and previous time step)               !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! coefa, coefb     ! ra ! <-- ! boundary conditions at the boundary faces      !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the previous time step        !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
! parbor(nfabor    ! ra ! <-- ! cumulation of the boundary statistics          !
!    nvisbr)       !    !     !                                                !
! vitpar           ! ra ! <-- ! part. velocity for the treatment of the        !
!   (nbpmax,3)     !    !     ! particle/wall interactions                     !
! vitflu           ! ra ! <-- ! flow velocity for the treatment of the         !
!   (nbpmax,3)     !    !     ! particle/wall interactions                     !
! auxl(nbpmax,3    ! ra ! --- ! work array                                     !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use numvar
use optcal
use entsor
use cstphy
use parall
use period
use lagpar
use lagran
use ppppar
use ppthch
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          lndnod
integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          kface  , nbpt   , isuivi

integer          itypfb(nfabor) , itrifb(nfabor)
integer          ifrlag(nfabor) , itepa(nbpmax,nivep)
integer          indep(nbpmax)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision parbor(nfabor,nvisbr) , auxl(nbpmax,3)
double precision vitpar(nbpmax,3) , vitflu(nbpmax,3)

! Local variables

integer          depch
integer          ip , nfin , kzone , n1 , icha, iok

double precision aa
double precision xp , yp , zp
double precision xq , yq , zq
double precision xk , yk , zk
double precision xpq , ypq , zpq
double precision xnn , ynn , znn
double precision vnorl(1)  , enc3 , viscp , masse
double precision dpinit , dp03 , mp0 , trap , vnorm , ang
double precision energ , energt
double precision uxn   , vyn    , wzn

!===============================================================================

!===============================================================================
! 0.  Memory management
!===============================================================================


!===============================================================================
! 1. Treatment with respect to the type of boundary
!===============================================================================

iok = 0

!--> number of the treated particle

ip = nbpt

!--> indicator of mass flux calculation

depch = 1

!--> Zone of the boundary face to be treated

kzone = ifrlag(kface)

!--> Normalized normale getting out from the face KFACE

aa = 1.d0 / surfbn(kface)
xnn = surfbo(1,kface) * aa
ynn = surfbo(2,kface) * aa
znn = surfbo(3,kface) * aa

!===============================================================================
! 2. Search of the intersection point between the boundary face and the ray.
! The coordinates are stored in XK YK ZK
!===============================================================================

!
!            1)  Equation of a plan, of which normal vector has coordinates (a,b,c):
!            2)  Equation of a line that contains points P and Q:
!                x = XP + (XQ-XP) * AA
!                y = YP + (YQ-YP) * AA
!                z = ZP + (ZQ-ZP) * AA
!                where AA is a parameter that varies in the real ensemble.

!-->We determine the vector PQ:

xp = ettpa(ip,jxp)
yp = ettpa(ip,jyp)
zp = ettpa(ip,jzp)

xq = ettp(ip,jxp)
yq = ettp(ip,jyp)
zq = ettp(ip,jzp)

xpq = xq - xp
ypq = yq - yp
zpq = zq - zp

!-->if the particle has not moved (if it is deposited on the boundary face),
!   it is not treated anymore

if (xpq.eq.0.d0 .and. ypq.eq.0.d0 .and. zpq.eq.0.d0) return

!--> From the equation of the plan of the face and the parametric equation
!    of the ray, the intersection point is determined

aa = xpq * surfbo(1,kface)                                        &
   + ypq * surfbo(2,kface)                                        &
   + zpq * surfbo(3,kface)

if ( aa.eq.0.d0 ) then
  write (nfecra,9010) ip
  nbperr = nbperr + 1
  dnbper = dnbper + tepa(ip,jrpoi)
  isuivi = 0
  itepa(ip,jisor) = 0
  return
endif

aa =                                                              &
     ( surfbo(1,kface) * cdgfbo(1,kface)                          &
     + surfbo(2,kface) * cdgfbo(2,kface)                          &
     + surfbo(3,kface) * cdgfbo(3,kface)                          &
     - surfbo(1,kface) * xp                                       &
     - surfbo(2,kface) * yp                                       &
     - surfbo(3,kface) * zp )                                     &
     / aa

!--> The aa parameter is injected into the equation of the right of the ray to
! get the intersection point of coordinates (XK YK ZK)


xk = xp + xpq * aa
yk = yp + ypq * aa
zk = zp + zpq * aa

!===============================================================================
! 3. If the particle deposits, the number of deposited particles is updated
!===============================================================================

if (iusclb(kzone).eq.idepo1 .or.                                 &
    iusclb(kzone).eq.idepo2 .or.                                 &
    iusclb(kzone).eq.idepo3      ) then

  nbpdep = nbpdep + 1
  dnbdep = dnbdep + tepa(ip,jrpoi)

endif

!===============================================================================
! 3. Departure of the particle from the calculation domain
!    or deposition on a boundary
!===============================================================================

if ( iusclb(kzone).eq.isortl .or.                                 &
     iusclb(kzone).eq.ientrl .or.                                 &
     iusclb(kzone).eq.idepo1      ) then

  isuivi = 0
  itepa(ip,jisor) = 0

!      update of the flow

  deblag(kzone) = deblag(kzone)-tepa(ip,jrpoi)*ettp(ip,jmp)

!--> The particle gets out, but for the Ensight visualization,
!    it is placed correctly at the intersection point
!
!

  ettp(ip,jxp) = xk
  ettp(ip,jyp) = yk
  ettp(ip,jzp) = zk

!===============================================================================
! 4. Deposition of the particle, which remains in memory
!===============================================================================

else if (iusclb(kzone).eq.idepo2) then

!--> The particle does not get out of the domain, it is not treated any more
!    but can still be visualized. The IP number is not reusable.

  isuivi = 0
  itepa(ip,jisor) = -itepa(ip,jisor)
  ettp(ip,jxp) = xk
  ettp(ip,jyp) = yk
  ettp(ip,jzp) = zk

  do n1 = 1,3
     vitpar(ip,n1) = 0.d0
     vitflu(ip,n1) = 0.d0
  enddo

!===============================================================================
! 5. Deposition of the particle, the resuspension is possible
!===============================================================================

else if (iusclb(kzone).eq.idepo3) then

  isuivi = 0
  itepa(ip,jisor) = ifabor(kface)
  ettp(ip,jxp) = xk
  ettp(ip,jyp) = yk
  ettp(ip,jzp) = zk

  do n1 = 1,3
    vitpar(ip,n1) = 0.d0
    vitflu(ip,n1) = 0.d0
  enddo
  do n1 = jup,jwf
    ettpa(ip,n1) = 0.d0
  enddo

!===============================================================================
! 6. Deposition of the particle with DLVO deposition conditions
!===============================================================================

else if (iusclb(kzone).eq.idepfa) then


! Calculation of the criterion

  uxn = ettp(ip,jup)*xnn
  vyn = ettp(ip,jvp)*ynn
  wzn = ettp(ip,jwp)*znn

  energ = 0.5d0*ettp(ip,jmp)*(uxn+vyn+wzn)**2

  energt   = 3.34d-12*ettp(ip,jdp)

  if ( energ .ge. energt )then

! The particle deposits:

    nbpdep = nbpdep + 1
    dnbdep = dnbdep + tepa(ip,jrpoi)

    isuivi = 0
    itepa(ip,jisor) = -itepa(ip,jisor)

    ettp(ip,jxp) = xk
    ettp(ip,jyp) = yk
    ettp(ip,jzp) = zk

    vitpar(ip,1) = 0.d0
    vitpar(ip,2) = 0.d0
    vitpar(ip,3) = 0.d0

    vitflu(ip,1) = 0.d0
    vitflu(ip,2) = 0.d0
    vitflu(ip,3) = 0.d0

  else

! The particle does not deposit:
! It 'rebounds' on the energy barrier:

  isuivi = 1
  itepa(ip,jisor) = ifabor(kface)

!  The mass flux is not calculated
!
    depch = 0

!-->Modification of the starting point

  ettpa(ip,jxp) = xk
  ettpa(ip,jyp) = yk
  ettpa(ip,jzp) = zk

  if (iensi1.eq.1) then
     nfin = 0
     call enslag                                                 &
          !==========
        ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,             &
          nfin   , ip     ,                                        &
          itepa  ,                                                 &
          ettpa  , tepa   )
  endif

  !-->Modification of the arrival point
  !   (the absolute value is intended to avoid the negative scalar products
  !  that may occur due to computer round-off error

  aa = 2.d0 * abs( (xq-xk)*xnn + (yq-yk)*ynn + (zq-zk)*znn )

  ettp(ip,jxp) = xq - aa*xnn
  ettp(ip,jyp) = yq - aa*ynn
  ettp(ip,jzp) = zq - aa*znn

  !--> Modification of the particle velocity at the arrival point

!-->Modification of the particle velocity at the impaction point
!   (like an elastic rebound)

    aa = abs(( vitpar(ip,1)*xnn                                   &
              +vitpar(ip,2)*ynn                                   &
              +vitpar(ip,3)*znn) )*2.d0

    vitpar(ip,1) = vitpar(ip,1) - aa*xnn
    vitpar(ip,2) = vitpar(ip,2) - aa*ynn
    vitpar(ip,3) = vitpar(ip,3) - aa*znn

  !--> Modification of the velocity of the flow seen at the arrival point

    aa = abs( (vitflu(ip,1)*xnn                                   &
             + vitflu(ip,2)*ynn                                   &
             + vitflu(ip,3)*znn) ) * 2.d0

  vitflu(ip,1) = vitflu(ip,1) - aa*xnn
  vitflu(ip,2) = vitflu(ip,2) - aa*ynn
  vitflu(ip,3) = vitflu(ip,3) - aa*znn


  endif

!===============================================================================
! 7. Elastic rebound of the particle on the boundary
!===============================================================================

else if (iusclb(kzone).eq.irebol) then

  isuivi = 1
  itepa(ip,jisor) = ifabor(kface)

!-->Modification of the starting point

  ettpa(ip,jxp) = xk
  ettpa(ip,jyp) = yk
  ettpa(ip,jzp) = zk

    if (iensi1.eq.1) then
      nfin = 0
      call enslag                                                 &
      !==========
       ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,             &
         nfin   , ip     ,                                        &
         itepa  ,                                                 &
         ettpa  , tepa   )
    endif

!-->Modification of the arrival point
!   (the absolute value is intended to avoid the negative scalar products
!  that may occur due to computer round-off error

    aa = 2.d0 * abs( (xq-xk)*xnn + (yq-yk)*ynn + (zq-zk)*znn )

  ettp(ip,jxp) = xq - aa*xnn
  ettp(ip,jyp) = yq - aa*ynn
  ettp(ip,jzp) = zq - aa*znn

!--> Modification of the particle velocity at the arrival point


  aa = abs( (vitpar(ip,1)*xnn                                     &
          +  vitpar(ip,2)*ynn                                     &
          +  vitpar(ip,3)*znn) ) * 2.d0

  vitpar(ip,1) = vitpar(ip,1) - aa*xnn
  vitpar(ip,2) = vitpar(ip,2) - aa*ynn
  vitpar(ip,3) = vitpar(ip,3) - aa*znn

!--> Modification of the velocity of the flow seen at the arrival point

  aa = abs( (vitflu(ip,1)*xnn                                     &
          +  vitflu(ip,2)*ynn                                     &
          +  vitflu(ip,3)*znn) ) * 2.d0

  vitflu(ip,1) = vitflu(ip,1) - aa*xnn
  vitflu(ip,2) = vitflu(ip,2) - aa*ynn
  vitflu(ip,3) = vitflu(ip,3) - aa*znn

!===============================================================================
! 8. Fouling of coal particles
!===============================================================================

else if (iusclb(kzone).eq.iencrl) then

!--> Fouling of the particle, if its properties make it possible
!    and with respect to a probability
!      ICI if  Tp     > TPENC
!          if  VISCP  > VISCREF

  icha = itepa(ip,jinch)

  if ( ettp(ip,jhp).gt.tprenc(icha) ) then

    enc3 = ( (1.d+7 * enc1(icha))/((ettp(ip,jhp)-150.d0)**2) )    &
           + enc2(icha)
      viscp = exp( log(10.d0)*enc3 ) * 0.1d0

      trap = 1.d0
      if ( viscp.gt.visref(icha) ) then
        n1 = 1
        call zufall(n1,vnorl(1))
        trap = 1.d0-visref(icha) / viscp
      endif

!  If VISCP <= VISREF ===> Probability of fouling equal to 1
!  If VISCP  > VISREF ===> Probability equal to TRAP = 1-VISREF/VISCP
!
!                     ===> Fouling if VNORL is between
!                          TRAP et 1.

      if ( viscp.le.visref(icha) .or.                             &
         (viscp.gt.visref(icha) .and. vnorl(1).ge.trap) ) then

! The computation of the mass of coal particles fouled is carried out
! in a following section

         npencr = npencr + 1
         isuivi = 0
         itepa(ip,jisor)  =  0
         ettp(ip,jxp) = xk
         ettp(ip,jyp) = yk
         ettp(ip,jzp) = zk

      endif
    endif

!--> if there is no fouling, then it is an elastic rebound

    if ( itepa(ip,jisor).ne.0 ) then

      isuivi = 1
      itepa(ip,jisor) = ifabor(kface)

!--> Modification of the departure point

    ettpa(ip,jxp) = xk
    ettpa(ip,jyp) = yk
    ettpa(ip,jzp) = zk

      if (iensi1.eq.1) then
        nfin = 0
        call enslag                                               &
        !==========
       ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,             &
         nfin   , ip     ,                                        &
         itepa  ,                                                 &
         ettpa  , tepa   )
      endif

!--> Modification of the arrival point

      aa = 2.d0 * abs((xq-xk)*xnn + (yq-yk)*ynn + (zq-zk)*znn)

    ettp(ip,jxp) = xq - aa*xnn
    ettp(ip,jyp) = yq - aa*ynn
    ettp(ip,jzp) = zq - aa*znn

    endif

  if (itepa(ip,jisor).gt.0) then

!--> Modification of the particle velocity at the arrival point


    aa = abs( (vitpar(ip,1)*xnn                                   &
            +  vitpar(ip,2)*ynn                                   &
            +  vitpar(ip,3)*znn) ) * 2.d0

    vitpar(ip,1) = vitpar(ip,1) - aa*xnn
    vitpar(ip,2) = vitpar(ip,2) - aa*ynn
    vitpar(ip,3) = vitpar(ip,3) - aa*znn

!--> Modification of the velocity of the flow seen at the arrival point

    aa = abs( (vitflu(ip,1)*xnn                                   &
            +  vitflu(ip,2)*ynn                                   &
            +  vitflu(ip,3)*znn) ) * 2.d0

    vitflu(ip,1) = vitflu(ip,1) - aa*xnn
    vitflu(ip,2) = vitflu(ip,2) - aa*ynn
    vitflu(ip,3) = vitflu(ip,3) - aa*znn

    endif

!===============================================================================
! 9. User-defined interaction number 1 : JBORD1
!===============================================================================

!  The following procedure is also valid for JBORD2, JBORD3, JBORD4 et JBORD5
!  The example is given only for JBORD1

!     We first check if we are in the zone of interest:
!      ELSE IF (IUSCLB(KZONE).EQ.JBORD1) THEN

!     if we need to keep on following the particle
!         ISUIVI = 0 OU 1

!     the mesh element of interest
!         ITEPA(IP,JISOR) =

!     modification of the arrival point
!         ETTP(IP,JXP) =
!         ETTP(IP,JYP) =
!         ETTP(IP,JZP) =

!     modification of the particle velocity at the arrival point
!         VITPAR(IP,1) =
!         VITPAR(IP,2) =
!         VITPAR(IP,3) =

!      modification of the velocity of the flow seen at the arrival point
!         VITFLU(IP,1) =
!         VITFLU(IP,2) =
!         VITFLU(IP,3) =


else if (iusclb(kzone).eq.jbord1                                  &
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!    For a standard use, without intervention of the user,
!    we do not wish to go through this part but we want the
!    test  IUSCLB(KZONE).EQ.JBORD1 to be in the us* example
!    and the following source code to be compiled to check for errors
         .and.(0.eq.1)                                            &
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END
                                ) then

!     ----------------------------------------------------
!     Example 1: The particle has 50% probality to definitely deposit
!                 and 50% to bounce back to the flow
!     ----------------------------------------------------


    n1 = 1
    call zufall(n1,vnorl(1))
    trap = 0.5d0

    if (vnorl(1).ge.trap) then

      isuivi = 0
      itepa(ip,jisor)  =  0
      ettp(ip,jxp) = xk
      ettp(ip,jyp) = yk
      ettp(ip,jzp) = zk

    else

      isuivi = 1
      itepa(ip,jisor) = ifabor(kface)

!-->Modification of the departure point

      ettpa(ip,jxp) = xk
      ettpa(ip,jyp) = yk
      ettpa(ip,jzp) = zk

!-->Modification of the arrival point

      aa = 2.d0 * abs((xq-xk)*xnn + (yq-yk)*ynn + (zq-zk)*znn)

      ettp(ip,jxp) = xq - aa*xnn
      ettp(ip,jyp) = yq - aa*ynn
      ettp(ip,jzp) = zq - aa*znn

    endif

!-->No need to treat the particles with ITEPA(IP,JISOR)=0 because
!   they will be removed from the particle list

    if (itepa(ip,jisor).gt.0) then

!-->Modification of the particle velocity at the arrival point

      aa = abs( (vitpar(ip,1)*xnn                                 &
              +  vitpar(ip,2)*ynn                                 &
              +  vitpar(ip,3)*znn) ) * 2.d0

      vitpar(ip,1) = vitpar(ip,1) - aa*xnn
      vitpar(ip,2) = vitpar(ip,2) - aa*ynn
      vitpar(ip,3) = vitpar(ip,3) - aa*znn

!-->Modification of the velocity of the flow seen at the arrival point

      aa = abs( (vitflu(ip,1)*xnn                                 &
              +  vitflu(ip,2)*ynn                                 &
              +  vitflu(ip,3)*znn) ) * 2.d0

      vitflu(ip,1) = vitflu(ip,1) - aa*xnn
      vitflu(ip,2) = vitflu(ip,2) - aa*ynn
      vitflu(ip,3) = vitflu(ip,3) - aa*znn

    endif


!===============================================================================
! 10. Verification and exit if error
!===============================================================================

else
  write (nfecra,9020) kzone
  iok = iok + 1
endif

if (iok.ne.0) then
  call csexit (1)
  !==========
endif

!===============================================================================
! 11. Recording of the particle/boundary interaction if needed
!===============================================================================

! The recording of wall statistics start as soon as the parameter IENSI3
! is set to 1. However, as long as the absolute number of the Lagrangian iteration
! is inferior to NSTBOR, or if the flow is unsteady  (ISTTIO = 0); the array PARBOR
! is reset to 0 before entering this surboutine.

! NPSTF :  number of iteractions of computation of statistics
!          at the unsteady boundaries

! NPSTFT : total number of statistics at the boundaries since the
!          beginning of the computation, included the unsteady part
!         (to be used only for the listing post-processing).
!

! TSTATP : physical duration of the recording of the statistics
!          of the interactions between the particles and the stationary boundaries,
!          if unsteady then it is equal DTP the last Lagrangian time step.
!

!
!
!
!



!    The following lines are only indications:


!     DO IVAR = 1,NVISBR

!      IF (IMOYBR(IVAR).EQ.2) THEN

!         DO IFAC = 1,NFABOR
!           IF (PARBOR(IFAC,INBR).GT.SEUILF) THEN
!             PARBOR(IFAC,IVAR) = PARBOR(IFAC,IVAR) /PARBOR(IFAC,INBR)
!           ELSE
!             PARBOR(IFAC,IVAR) = 0.D0
!           ENDIF
!         ENDDO

!       ELSE IF (IMOYBR(IVAR).EQ.1) THEN

!         DO IFAC = 1,NFABOR
!           IF (PARBOR(IFAC,INBR).GT.SEUILF) THEN
!             PARBOR(IFAC,IVAR) = PARBOR(IFAC,IVAR) / TSTATP
!           ELSE
!             PARBOR(IFAC,IVAR) = 0.D0
!           ENDIF
!         ENDDO
!       ENDIF
!     ENDDO




if ( iensi3.eq.1 ) then

!--> Example of types of interactions about which we want to
!    record information

  if ( iusclb(kzone).eq.irebol .or.                               &
       iusclb(kzone).eq.idepo1 .or.                               &
       iusclb(kzone).eq.idepo2 .or.                               &
       iusclb(kzone).eq.idepo3 .or.                               &
       iusclb(kzone).eq.idepfa ) then

    if (inbrbd.eq.1) then
      parbor(kface,inbr) = parbor(kface,inbr) + tepa(ip,jrpoi)
    endif

    if (iflmbd.eq.1 .and. depch.eq.1) then
        parbor(kface,iflm) = parbor(kface,iflm)                   &
     + ( tepa(ip,jrpoi) * ettp(ip,jmp) /surfbn(kface) )
    endif

    if (iangbd.eq.1) then
      vnorm = ettp(ip,jup) * ettp(ip,jup)                         &
            + ettp(ip,jvp) * ettp(ip,jvp)                         &
            + ettp(ip,jwp) * ettp(ip,jwp)
      vnorm = sqrt( vnorm )
      ang =  ettp(ip,jup) * surfbo(1,kface)                       &
           + ettp(ip,jvp) * surfbo(2,kface)                       &
           + ettp(ip,jwp) * surfbo(3,kface)                       &
           / surfbn(kface)                                        &
           / vnorm
      ang = acos(ang)
      parbor(kface,iang) = parbor(kface,iang) + ang*tepa(ip,jrpoi)
    endif

    if (ivitbd.eq.1) then
      vnorm = ettp(ip,jup) * ettp(ip,jup)                         &
            + ettp(ip,jvp) * ettp(ip,jvp)                         &
            + ettp(ip,jwp) * ettp(ip,jwp)
      vnorm = sqrt( vnorm )
      parbor(kface,ivit) =parbor(kface,ivit) +vnorm*tepa(ip,jrpoi)
    endif

    if (nusbor.gt.0) then
      do n1 = 1,nusbor
        parbor(kface,iusb(n1)) = 0.d0
      enddo
    endif

!--> Particular case of the mass of fouled coal

  else if ( iusclb(kzone).eq.iencrl .and. isuivi.eq.0 ) then

    parbor(kface,inbr) = parbor(kface,inbr) + tepa(ip,jrpoi)

    if (iencbd.eq.1) then

      icha =  itepa(ip,jinch)
      dpinit = tepa(ip,jrd0p)

      dp03 = dpinit * dpinit * dpinit
      mp0  = pi * dp03 * rho0ch(icha) / 6.d0

      masse = ettp(ip,jmch) + ettp(ip,jmck)                       &
                                + xashch(icha) * mp0

      parbor(kface,ienc) = parbor(kface,ienc)                     &
                         + tepa(ip,jrpoi)*masse

    endif

  endif

endif


!===============================================================================
! Archives. This part is left here as it may be useful..
!           Creation of a local referential associated to a boundary face
!===============================================================================

! The local referential (T1,T2,N) is built so that N is
! the normalized normal of the face, and that T1 and T2 belong to the face

!-->1. I know N and PK, I define  T1 so that
!      T1 = PK *  N (cross product)

!     XPK = XK - ETTPA(IP,JXP)
!     YPK = YK - ETTPA(IP,JYP)
!     ZPK = ZK - ETTPA(IP,JZP)

!     XT1 = YPK*ZNN - ZPK*YNN
!     YT1 = ZPK*XNN - XPK*ZNN
!     ZT1 = XPK*YNN - YPK*XNN

!     AA = SQRT(XT1*XT1 + YT1*YT1 + ZT1*ZT1)
!     XT1 = XT1 / AA
!     YT1 = YT1 / AA
!     ZT1 = ZT1 / AA

!-->2. Now I can define T2 = - T1 * N

!     XT2 = YT1*ZNN - ZT1*YNN
!     YT2 = ZT1*XNN - XT1*ZNN
!     ZT2 = XT1*YNN - YT1*XNN

!     AA = SQRT(XT2*XT2 + YT2*YT2 + ZT2*ZT2)
!     XT2 = -XT2 / AA
!     YT2 = -YT2 / AA
!     ZT2 = -ZT2 / AA

!--------
! Formats
!--------

 9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : STOP IN THE EXECUTION OF THE LAGRANGIAN MODULE   ',/,&
'@    =========   (USLABO)                                    ',/,&
'@                                                            ',/,&
'@  the normal of a boundary face is perpendicular to   ',/,&
'@   a PQ ray : impossible                               ',/,&
'@                                                            ',/,&
'@  The particle ',I10,' is ELIMINATED                          ',/,&
'@                                                            ',/,&
'@  Please contact the development team.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : STOP IN THE EXECUTION OF THE LAGRANGIAN MODULE  ',/,&
'@    =========   (USLABO)                                    ',/,&
'@                                                            ',/,&
'@  The type of boundary condition IUSCLB                   ',/,&
'@    is not defined for the boundary NB = ',I10        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be run.                           ',/,&
'@                                                            ',/,&
'@  Please check USLAG2 and USLABO.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return
end subroutine
