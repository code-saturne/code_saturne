!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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


!===============================================================================


subroutine usladp &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   , statis ,                            &
   taup   , tlag   , piil   ,                                     &
   vagaus , gradpr , gradvf ,                                     &
   romp   ,                                                       &
   dppar  , dnxpar , dnypar , dnzpar )

!===============================================================================
! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module :
!   -------------------------------------

!    User subroutine (non-mandatory intervention)

!    For each particle :
!      - Computation of the wall-normal distance
!      - Computation of the normal to the wall

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! ibord            ! ia ! --> ! if nordre=2, contains the number               !
!   (nbpmax)       !    !     ! of the boundary face of part./wall interaction !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! transported variables at the current           !
! (ncelet,*)       !    !     ! and previous time step                         !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the previous time step        !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
! statis           ! ra ! <-- ! cumul for the averages of the volume stats.    !
!(ncelet,nvlsta    !    !     !                                                !
! taup(nbpmax)     ! ra ! <-- ! particle relaxation time                       !
! tlag(nbpmax)     ! ra ! <-- ! relaxation time for the flow                   !
! piil(nbpmax,3    ! ra ! <-- ! term in the sede integration                   !
! tsup(nbpmax,3    ! ra ! <-- ! prediction 1st substep                         !
!                  !    !     ! for the particle velocity                      !
! tsuf(nbpmax,3    ! ra ! <-- ! prediction 1st substep                         !
!                  !    !     ! for the velocity of the flow seen              !
! bx(nbpmax,3,2    ! ra ! <-- ! characteristics of the turbulence              !
! tsfext(nbpmax    ! ra ! <-- ! infos for the return coupling                  !
! vagaus           ! ra ! <-- ! Gaussian random variables                      !
!(nbpmax,nvgaus    !    !     !                                                !
! gradpr(ncel,3    ! ra ! <-- ! pressure gradient                              !
! gradvf(ncel,3    ! ra ! <-- ! flow-velocity gradient                         !
! romp             ! ra ! --- ! particle densite                               !
! fextla           ! ra ! --> ! exterior force field                           !
!(ncelet,3)        !    !     !                                                !
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
use cstnum
use cstphy
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep)

double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp),ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,*)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3)
double precision vagaus(nbpmax,*)
double precision dppar(nbpart)  , dnxpar(nbpart)
double precision dnypar(nbpart) , dnzpar(nbpart)

double precision gradpr(ncelet,3) , gradvf(ncelet,9)
double precision romp(nbpmax)

! Local variables

integer          ip

! User local variables

double precision xnorm

!===============================================================================

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END
!===============================================================================
! 0.  Memory management
!===============================================================================


!===============================================================================
! 1. Example
!===============================================================================

!   This example is unactivated

if (1.eq.0) then

  do ip = 1,nbpart

! Example : for every paricle, we take as
!      - Wall-normal distance: the radius of the duct
!      - Normal: with respect to the radius, by supposing
!                the duct along the Z-axis. Be careful, the
!                convention of Code_Saturne is respected, with
!                the normal oriented from the fluid towards the outside of the domain


    dppar(ip)  = 0.00631942286d0                                  &
                -sqrt( ettp(ip,jxp)*ettp(ip,jxp)                  &
                      +ettp(ip,jyp)*ettp(ip,jyp) )

    xnorm = sqrt( ettp(ip,jxp)*ettp(ip,jxp)                       &
                 +ettp(ip,jyp)*ettp(ip,jyp) )
    dnxpar(ip) = ettp(ip,jxp)/xnorm
    dnypar(ip) = ettp(ip,jyp)/xnorm
    dnzpar(ip) = 0.d0

  enddo

!==============================================================================
! Control: do not modify
!==============================================================================

  do ip = 1,nbpart

    if ( dppar(ip) .le. dparmn ) then
      dppar(ip) = dparmn-dparmn/100.d0
    endif

  enddo

endif

!--------
! Formats
!--------


!----
! End
!----

end subroutine


!===============================================================================


subroutine uslaed &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  , ibord  ,                                              &
   dt     , rtp    , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct , tsvar  , &
   auxl1  , auxl2  , auxl3  )

!===============================================================================
! Purpose :
! ----------

!   Subroutine of the Lagrangian particle-tracking module :
!   -------------------------------------

!     User subroutine (non-mandatory intervention)

!     Integration of the sde for the user-defined variables.
!     The variables are constant by default.


!                                         d T       T - PIP
!     The sde must be of the form:       ----- = - ---------
!                                         d t         Tca


!     T : IIIIeme user-defined variable, given for the ip particle by
!            T = ETTP(IP,JVLS(IIII))
!            T = ETTPA(IP,JVLS(IIII))

!     Tca : Characteristic time for the sde
!           to be prescribed in the array auxl1

!     PIP : Coefficient of the sde (pseudo right member)
!           to be prescribed in the array auxl2
!
!           If the chosen scheme is first order (nordre=1)
!           then, at the first and only passage pip is expressed
!           as a function of the quantities of the previous time step contained in ettpa
!
!           If the chosen scheme is second order (nordre=2)
!           then, at the first passage (nor=1) pip is expressed as
!           a function of the quantities of the previous time step contained in ettpa,
!           and at the second passage (nor=2) pip is expressed as
!           a function of the quantities of the current time step contained in ettp

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! ibord            ! ia ! <-- ! number of the boundary face of part./wall      !
!   (nbpmax)       !    !     ! interaction                                    !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! ra ! <-- ! transported variables at the current           !
! (ncelet,*)       !    !     ! and previous time step                         !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the previous time step        !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
! taup(nbpmax)     ! ra ! <-- ! particle relaxation time                       !
! tlag(nbpmax)     ! ra ! <-- ! relaxation time for the flow                   !
! tempct           ! ra ! <-- ! characteristic thermal time and                !
!  (nbpmax,2)      !    !     ! implicit source term of return coupling        !
! tsvar            ! ra ! <-- ! prediction 1st substep for the ivar variable,  !
! (nbpmax,nvp1)    !    !     ! used for the correction at the 2nd substep     !
!                  !    !     !                                                !
! auxl1(nbpmax)    ! ra ! ---   work array                                     !
! auxl2(nbpmax)    ! ra ! --- ! work array                                     !
! auxl3(nbpmax)    ! ra ! --- ! work array                                     !
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
use cstphy
use cstnum
use optcal
use entsor
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep)  , ibord(nbpmax)

double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision tsvar(nbpmax,nvp1)
double precision auxl1(nbpmax), auxl2(nbpmax), auxl3(nbpmax)

! Local variables

integer          npt , iel , iiii , ipl

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

! We enter this subroutine only if additional variables have been defined in
! uslag1; we must then necessarily define how they are solved.


if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION: STOP IN THE LAGRANGIAN MODULE                  ',/,&
'@    =========                                               ',/,&
'@     THE USER SUBROUTINE uslaed MUST BE FILLED              ',/,&
'@                                                            ',/,&
'@  The calculation will not be run                           ',/,&
'@                                                            ',/,&
'@  Additional variables have been declared in                ',/,&
'@    uslag1 (NVLS=)                                          ',/,&
'@  The subroutine uslaed must be filled to precise           ',/, &
'@    the stochastic differential equation to be solved       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


! FIXME : TODO write a user example


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  Initializations
!===============================================================================


!===============================================================================
! 2. Characteristic time of the current sde
!===============================================================================

! Loop on the additional variables

do iiii = 1,nvls

!      Number of the treated variable in ettp

  ipl = jvls(iiii)

  do npt = 1,nbpart

    if ( itepa(npt,jisor).gt.0 ) then

      iel = itepa(npt,jisor)

!     Characteristic time tca of the differential equation
!     This example must be adapted to the case

      auxl1(npt) = 1.d0

!     Prediction at the first substep
!     This example must be adapted to the case

      if (nor.eq.1) then
        auxl2(npt) = ettpa(npt,ipl)
      else

!     Correction at the second substep
!     This example must be adapted to the case

        auxl2(npt) = ettp(npt,ipl)
      endif

    endif
  enddo

!===============================================================================
! 3. Integration of the variable ipl
!===============================================================================

  call lagitg                                                     &
  !==========
   ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     ipl    ,                                                     &
     itepa(1,jisor)  , ibord  ,                                   &
     ettp   , ettpa  , auxl1  , auxl2  , tsvar  )

enddo

!===============================================================================

!----
! End
!----

end subroutine


!===============================================================================


subroutine uslafe &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  , ibord  ,                                              &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   ettp   , ettpa  , tepa   , statis , stativ ,                   &
   taup   , tlag   , piil   ,                                     &
   tsuf   , tsup   , bx     , tsfext ,                            &
   vagaus , gradpr , gradvf ,                                     &
   romp   , fextla )

!===============================================================================
! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module:
!   -------------------------------------

!    User subroutine (non-mandatory intervention)

!    Management of an external force field acting on the particles
!    It must be prescribed in every cell and be homogeneous to gravity (m/s^2)
!
!    By default gravity and drag force are the only forces acting on the particles
!    (the gravity components gx gy gz are assigned in usini1)
!

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep)   !    !     !                                                !
! ibord(nbpmax)    ! ia ! --> ! if nordre=2, contains the number of the        !
!                  !    !     ! boundary face of particle/wall interaction     !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! transported variables at cell centers for      !
! (ncelet,*)       !    !     ! the current and the previous time step         !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     !                                                !
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     !                                                !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
!  (nbpmax,nvep)   !    !     !                                                !
! statis           ! ra ! <-- ! cumul for the averages of the volume stats.    !
!  (ncelet,nvlsta) !    !     !                                                !
! stativ           ! ra ! <-- ! cumulation for the variance of the volume      !
!  (ncelet,        !    !     ! statistics                                     !
!   nvlsta-1)      !    !     !                                                !
! taup(nbpmax)     ! ra ! <-- ! particle relaxation time                       !
! tlag(nbpmax)     ! ra ! <-- ! relaxation time for the flow                   !
! piil(nbpmax,3)   ! ra ! <-- ! term in the integration of the sde             !
! tsup(nbpmax,3)   ! ra ! <-- ! prediction 1st substep for                     !
!                  !    !     ! the velocity of the particles                  !
! tsuf(nbpmax,3)   ! ra ! <-- ! prediction 1st substep for                     !
!                  !    !     ! the velocity of the flow seen                  !
! bx(nbpmax,3,2)   ! ra ! <-- ! characteristics of the turbulence              !
! tsfext(nbpmax    ! ra ! <-- ! infos for the return coupling                  !
! vagaus           ! ra ! <-- ! Gaussian random variables                      !
!  (nbpmax,nvgaus) !    !     !                                                !
! gradpr(ncel,3)   ! ra ! <-- ! pressure gradient                              !
! gradvf(ncel,3)   ! ra ! <-- ! gradient of the flow velocity                  !
! romp             ! ra ! --- ! particle density                               !
! fextla(ncelet,3) ! ra ! --> ! user external force field (m/s^2)              !
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
use cstnum
use cstphy
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep) , ibord(nbpmax)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,*),stativ(ncelet,*)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3) , bx(nbpmax,3,2)
double precision tsuf(nbpmax,3) , tsup(nbpmax,3)
double precision tsfext(nbpmax)
double precision vagaus(nbpmax,*)
double precision gradpr(ncelet,3) , gradvf(ncelet,9)
double precision romp(nbpmax)
double precision fextla(nbpmax,3)

! Local variables

integer          ip


!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================
!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0. Memory management
!===============================================================================


!===============================================================================
! 1. Example
!===============================================================================

!   This example is unactivated

if (1.eq.0) then


  do ip = 1,nbpart

    fextla(ip,1) = 0.d0
    fextla(ip,2) = 0.d0
    fextla(ip,3) = 0.d0

  enddo


endif

!==============================================================================

!--------
! Formats
!--------


!----
! End
!----

end subroutine


!===============================================================================


subroutine uslain &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nptnew ,                                                       &
   itypfb , itrifb , itepa  , ifrlag , injfac ,                   &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus )

!===============================================================================
! Purpose:
! --------
!
! User subroutine of the Lagrangian particle-tracking module:
! -----------------------------------------
!
! User subroutine (non-mandatory intervention)

! User subroutine for the boundary conditions for the particles
! (inlet and treatment for the other boundaries)
!
! This routine is called after the initialization of the
! ettp, tepa and itepa arrays for the new particles in order to modify them
! to inject new particle profiles.
!

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
! nptnew           ! i  ! <-- ! total number of new particles for all the      !
!                  !    !     ! injection zones                                !
! itrifb(nfabor)   ! ia ! <-- ! indirection for the sorting of the boundary    !
! itypfb(nfabor)   ! ia ! <-- ! type of the boundary faces                     !
! ifrlag(nfabor)   ! ia ! --> ! type of the Lagrangian boundary faces          !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! injfac(npbmax)   ! ia ! <-- ! number of the injection boundary face          !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! transported variables at the previous timestep !
! (ncelet,*)       !    !     !                                                !
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
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
! vagaus           ! ra ! --> ! Gaussian random variables                      !
!(nbpmax,nvgaus    !    !     !                                                !
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
use optcal
use cstnum
use cstphy
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nptnew

integer          itypfb(nfabor) , itrifb(nfabor)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)
integer          injfac(nbpmax)

double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision vagaus(nbpmax,*)

! Local variables

integer          iclas , izone , ifac
integer          ii , ip , npt , npar1 , npar2, ipnorm

! User-defined local variables
! (the dimension of vgauss is 3, but 2 would be sufficient here)

double precision vgauss(3)

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

!     By default, we do not modify them

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

if (nbpnew.eq.0) return

!===============================================================================
! 1. Memory management
!===============================================================================


!===============================================================================
! 2. Initializations
!===============================================================================



!===============================================================================
! 3. Modification of properties of the new particles (injection profiles,
!    position of the injection point, statistical
!    weights, correction of the diameter if the standard-deviation option
!    is activated.)
!===============================================================================

!    These modifications occur after all the initializations related to
!    the particle injection, but before the treatment of the continuous
!    injection: it is thus possible to impose an injection profile with
!    the continous-injection option.
!

!   reinitialization of the counter of the new particles
npt = nbpart

!   for each boundary zone:
do ii = 1,nfrlag
  izone = ilflag(ii)

!       for each class:
  do iclas = 1, iusncl(izone)

!         if new particles must enter the domain:
    if (mod(ntcabs,iuslag(iclas,izone,ijfre)).eq.0) then

      do ip = npt+1 , npt+iuslag(iclas,izone,ijnbp)

!         number of the original boundary face of injection

      ifac = injfac(ip)
!
!-----------------------------------------------------------
!        EXAMPLE OF MODIFICATION OF THE INJECTION VELOCITY
!        WITH RESPECT TO THE INJECTION POSITION
!-----------------------------------------------------------
!    For instance, the user can call his own subroutine that provides
!    the three components of the instantaneous velocities ettp(ip,jup)
!    ettp(ip,jvp) and  ettp(ip,jwp) with respect to  ettp(ip,jzp)
!    (through interpolation for instance). More simply, the user can provide
!    the three components of the instantaneous velocities, under the form
!    of a mean value (taken arbitrarily here equal to (2,0,0) m/s) added
!    to a fluctuating value (equal here to 0,2 m/s for the 1st and 3rd components)
!
!
        ipnorm = 2
        call normalen(ipnorm,vgauss)
        ettp(ip,jup) = 2.d0 + vgauss(1) * 0.2d0
        ettp(ip,jvp) = 0.d0
        ettp(ip,jwp) = 0.d0 + vgauss(2) * 0.2d0

      enddo

      npt = npt + iuslag(iclas,izone,ijnbp)

    endif

  enddo
enddo

!===============================================================================
! 4. SIMULATION OF THE INSTANTANEOUS TURBULENT FLUID FLOW VELOCITIES SEEN
!    BY THE SOLID PARTICLES ALONG THEIR TRAJECTORIES.
!===============================================================================
!
! Entering this subroutine, the ettp(ip,juf) ettp(ip,jvf) and ettp(ip,jwf) arrays
! are filled with the components of the instantaneous velocity (fluctuation + mean value)
! seen by the particles
!
! When the velocity of the flow is modified just above, most of the time
! the user knows only the mean value. In some flow configurations and some
! injection conditions, it may be necessary to reconstruct the fluctuating part.
! That is why the following routine is called.
!
! Caution: this turbulent component must be reconstructed only on the modified
! velocities of the flow seen.
!
! The reconstruction is unactivated here and must be adapted to the case.
!

if ( 1.eq.0 ) then

  npar1 = nbpart+1
  npar2 = nbpart+nbpnew

  call lagipn                                                     &
  !==========
  ( ncelet , ncel   ,                                             &
    nbpmax , nvp    , nvp1   , nvep   , nivep  ,                  &
    npar1  , npar2  ,                                             &
    itepa  ,                                                      &
    rtpa   ,                                                      &
    ettp   , tepa   , vagaus )

endif

!===============================================================================

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine


!===============================================================================


subroutine uslapr &
!================

 ( idvar  , iepart , izone  , iclass ,                            &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   xxpart , yypart , zzpart ,                                     &
   tvpart , uupart , vvpart , wwpart , ddpart , ttpart  ,         &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   )

!===============================================================================
! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module:
!   -------------------------------------

!   User subroutine for the boundary conditions associated to
!   the particles (inlet and treatment of the other boundaries)
!
!   It allows to impose the values of the velocity, the diameter
!   and the temperature for the treated particle.
!
!   if idvar = 0 ==> the volume fraction is retrieved
!   if idvar = 1 ==> the 3 components of the velocity are retrieved
!   if idvar = 2 ==> the diameter is retrieved
!   if idvar = 3 ==> the temperature is retrieved


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idvar            ! i  ! <-- ! type of the value(s) ta calculate              !
! iepart           ! i  ! <-- ! number of the particle cell                    !
! izone            ! i  ! <-- ! number of the particle zone                    !
! iclass           ! i  ! <-- ! number of the particle class                   !
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
! itrifb(nfabor)   ! ia ! <-- ! indirection for the sorting of the             !
! itypfb(nfabor)   ! ia ! <-- ! type of the boundary faces                     !
! ifrlag(nfabor)   ! ia ! --> ! type of the Lagrangian boundary faces          !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! xxpart           !  r ! <-- ! x-coordinate of the particle                   !
! yypart           !  r ! <-- ! y-coordinate of the particle                   !
! zzpart           !  r ! <-- ! z-coordinate of the particle                   !
! tvpart           !  r ! <-- ! value of the volume fraction                   !
! uupart           !  r ! <-- ! x-component of particle velocity               !
! vvpart           !  r ! <-- ! y-component of particle velocity               !
! wwpart           !  r ! <-- ! z-component of particle velocity               !
! ddpart           !  r ! <-- ! particle diameter                              !
! ttpart           !  r ! <-- ! particle temperature                           !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! transported variables at the previous          !
! (ncelet,*)       !    !     ! time step                                      !
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
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
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
use optcal
use cstnum
use cstphy
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments


integer          idvar  , iepart , izone  , iclass

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itypfb(nfabor) , itrifb(nfabor)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)

double precision xxpart , yypart , zzpart
double precision tvpart , uupart , vvpart , wwpart
double precision ddpart , ttpart

double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)

! Local variables


double precision pis6

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START

!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
  !==========
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,&
'@    =========                                               ',/,&
'@     MODULE LAGRANGIEN :                                    ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uslapr DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. Memory management
!===============================================================================


!===============================================================================
! 2. Initialization
!===============================================================================

pis6 = pi / 6.d0

!===============================================================================
! 3. Profile for the volume fraction
!===============================================================================

if (idvar .eq. 0) then

  tvpart = 0.01

endif

!===============================================================================
! 4. Velocity profile
!===============================================================================

if (idvar .eq. 1) then

  uupart = 1.d0
  vvpart = 0.d0
  wwpart = 0.d0

endif

!===============================================================================
! 5. Diameter profile
!===============================================================================

if (idvar .eq. 2) then

  ddpart = 50.d-6

endif


!===============================================================================
! 6. Temperature profile
!===============================================================================

if (idvar .eq. 3) then

  ttpart = 20.d0

endif

!===============================================================================

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine


!===============================================================================


subroutine uslaru &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  ,                                     &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus , croule , auxl  ,                    &
   distpa , distyp )

!===============================================================================
! Purpose:
! --------
!
! User subroutine of the Lagrangian particle-tracking module:
! -----------------------------------------
!
! User subroutine (non-mandatory intervention)

! Calculation of the function of significance for the Russian roulette


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
! itypfb(nfabor)   ! ia ! <-- ! type of the boundary faces                     !
! itrifb(nfabor)   ! ia ! --> ! indirection for the sorting of the             !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep)   !    !     !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! transported variables at cell centers for      !
! (ncelet,*)       !    !     ! the previous timestep                          !
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
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
! vagaus           ! ra ! <-- ! Gaussian random variables                      !
!(nbpmax,nvgaus    !    !     !                                                !
! croule(ncelet    ! ra ! --> ! function of significance for                   !
!                  !    !     ! the Russian roulette                           !
! auxl(nbpmax,3    ! ra ! --- !                                                !
! distpa(ncelet    ! ra ! <-- ! wall-normal distance arrays                    !
! disty(ncelet)    ! ra ! <-- ! y+ distance                                    !
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
use optcal
use entsor
use cstphy
use parall
use period
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itypfb(nfabor) , itrifb(nfabor)
integer          itepa(nbpmax,nivep)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision vagaus(nbpmax,*) , croule(ncelet)
double precision auxl(nbpmax,3)
double precision distpa(ncelet) , distyp(ncelet)

! Local variables

integer          iel
double precision zref

!===============================================================================


!===============================================================================
! 0.  Memory management
!===============================================================================


!===============================================================================
! 1. Default initialization
!---------------------------

!     Caution : the croule parameter is only initialized in this subroutine.
!               Make sure that it is prescribed for every cell.


!===============================================================================

do iel = 1,ncel
  croule(iel) = 1.d0
enddo

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! -1.  If the user does not intervene, croule = 1 everywhere
!===============================================================================

if(1.eq.1) then
  return
endif

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 2. Calculation of a user-defined function of significance
!===============================================================================

!   CAUTION:   the croule array must be filled with positive
!   ^^^^^^^^^  real numbers enabling to weight the importance
!              of some zones with respect to others.
!
!              (the greater croule, the more important the zone)

!              For instance, we can decide that the zone is as important
!              as it is close to a position z=zref; with an importance equal
!              to 1.e-3 near zref, and with an importance equal to 1.e-6
!              far from zref.


zref = 0

do iel = 1,ncel
  croule(iel) = 1.d0/(max( abs(xyzcen(3,iel)-zref),1.d-3 ))
enddo

do iel = 1,ncel
  croule(iel) = max(croule(iel),1.d-6 )
enddo

!===============================================================================

!----
! End
!----

end subroutine


!===============================================================================


subroutine uslatc &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   numpt  , itepa  ,                                              &
   rep    , uvwr   , romf   , romp   , xnul   ,                   &
   xcp    , xrkl   , tauc   ,                                     &
   dt     , rtp    , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   )

!===============================================================================
! Purpose:
! --------
!
! User subroutine of the Lagrangian particle-tracking module:
! -----------------------------------------
!
! User subroutine (non-mandatory intervention)
!
! Modification of the computation of the thermal relaxation time
! of the particles with respect to the chosen formulation of the
! Nusselt number.

! This subroutine being called in a loop on the particle number,
! be careful not to "load" it to heavily..
!
!

!               m   Cp
!                p    p
!      Tau = ---------------
!         c          2
!               PI d    h
!                   p    e

!     Tau  : Thermal relaxation time (value to be computed)
!        c

!     m    : Particle mass
!      p

!     Cp   : Particle specific heat
!       p

!     d    : Particle diameter
!      p

!     h    : Coefficient of thermal exchange
!      e

!  The coefficient of thermal exchange is calculated from a Nusselt number,
!  itself evaluated by a correlation (Ranz-Marshall by default)
!
!

!            h  d
!             e  p
!     Nu = --------  = 2 + 0.55 Re **(0.5) Prt**(0.33)
!           Lambda                p

!     Lambda : Thermal conductivity of the carrier field

!     Re     : Particle Reynolds number
!       p

!     Prt    : Prandtl number

!

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! numpt            ! i  ! <-- !                                                !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! rep              ! r  ! <-- ! particle Reynolds number                       !
!                  !    !     ! rep = uvwr * ettp(numpt,jdp) / xnul            !
! uvwr             ! r  ! <-- ! relative velocity of the particle              !
!                  !    !     ! uvwr = |flow-seen velocity - part. velocity |  !
! romf             ! r  ! <-- ! fluid density at  particle position            !
!                  !    !     !                                                !
! romp             ! r  ! <-- ! particle density                               !
! xnul             ! r  ! <-- ! kinematic viscosity of the fluid at            !
!                  !    !     ! particle position                              !
! xcp              ! r  ! <-- ! specific heat of the fluid at particle         !
!                  !    !     ! position                                       !
! xrkl             ! r  ! <-- ! diffusion coefficient of the fluid at particle !
!                  !    !     ! position                                       !
! tauc             ! r  ! --> ! thermal relaxation time                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! ra ! <-- ! transported variables at cell centers at       !
! (ncelet,*)       !    !     ! the current time step                          !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the previous time step        !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
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
use cstnum
use cstphy
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          numpt

integer          itepa(nbpmax,nivep)

double precision rep    , uvwr   , romf   , romp   , xnul
double precision xcp    , xrkl   , tauc

double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)

! Local variables

integer          ip

! User-defined local variables

double precision prt, fnus

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0. Memory management
!===============================================================================


!===============================================================================
! 1. Initializations
!===============================================================================

ip = numpt

!===============================================================================
! 2. Standard thermal relaxation time
!===============================================================================

!   This example is unactivated, it gives the standard thermal relaxation time
!   as an indication.


if (1.eq.0) then

  prt  = xnul / xrkl

  fnus = 2.d0 + 0.55d0 * rep**0.5d0 * prt**(1.d0/3.d0)

  tauc = ettp(ip,jdp) *ettp(ip,jdp) * romp * ettp(ip,jcp)         &
           / ( fnus * 6.d0 * romf * xcp * xrkl )

endif

!==============================================================================

!--------
! Formats
!--------


!----
! End
!----

end subroutine


!===============================================================================


subroutine uslatp &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   numpt  , itepa  ,                                              &
   rep    , uvwr   , romf   , romp   , xnul   , taup   ,          &
   dt     , rtp    , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   )

!===============================================================================
! Purpose:
! --------
!
! User subroutine of the Lagrangian particle-tracking module:
! -----------------------------------------
!
! User subroutine (non-mandatory intervention)
!
! Modification of the calculation of the particle relaxation time
! with respect to the chosen formulation for the drag coefficient

! This subroutine being called in a loop on the particle number,
! be careful not to "load" it too heavily..
!
!            rho             4 d
!               p               p
!      Tau = ---- --------------------------------
!         p
!            rho   3 C     | U [X (t),t] - V (t) |
!               f     drag    f  p          p

!     Tau  : Particle relaxation time
!        p

!     rho  : Particle density
!        p

!     rho  : Fluid density
!        f

!     C    : Drag coefficient
!      drag

!     d    : Particle diameter
!      p

!     U [X (t),t] : Instantaneous velocity of the flow seen
!      f  p

!     V (t) : Particle velocity
!      p

!
!

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! numpt            ! i  ! <-- !                                                !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! rep              ! r  ! <-- ! particle Reynolds number                       !
!                  !    !     ! rep = uvwr * ettp(numpt,jdp) / xnul            !
! uvwr             ! r  ! <-- ! particle relative velocity                     !
!                  !    !     ! uvwr= |flow-seen velocity - part. velocity|    !
! romf             ! r  ! <-- ! fluid density at  particle position            !
!                  !    !     !                                                !
! romp             ! r  ! <-- ! particle density                               !
! xnul             ! r  ! <-- ! kinematic viscosity of the fluid at            !
!                  !    !     ! particle position                              !
! taup             ! r  ! --> ! particle relaxation time                       !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! ra ! <-- ! transported variables at cells centers         !
! (ncelet,*)       !    !     ! at the previous time step                      !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the previous time step        !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
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
use cstnum
use cstphy
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          numpt

integer          itepa(nbpmax,nivep)

double precision rep    , uvwr   , romf   , romp   , xnul  , taup

double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)

! Local variables

integer          ip
double precision fdr

! User-defined local variables

double precision cd1 , cd2 , dd2
double precision rec1, rec2, rec3, rec4

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0.  Memory management
!===============================================================================


!===============================================================================
! 1. Initializations
!===============================================================================

ip = numpt

!===============================================================================
! 2. Relaxation time with the standard (Wen-Yu) formulation of the drag coefficient
!===============================================================================

! This example is unactivated, it gives the standard relaxation time
! as an indication:

if (1.eq.0) then

  cd1  = 0.15d0
  cd2  = 0.687d0

  if (rep.le.1000) then
      dd2 = ettp(ip,jdp) * ettp(ip,jdp)
      fdr = 18.d0 * xnul * (1.d0 + cd1 * rep**cd2) / dd2
  else
      fdr = (0.44d0 * 3.d0 / 4.d0) * uvwr / ettp(ip,jdp)
  endif

  taup = romp / romf / fdr

endif

!===============================================================================
! 3. Computation of the relaxation time with the drag coefficient of
!    S.A. Morsi and A.J. Alexander, J. of Fluid Mech., Vol.55, pp 193-208 (1972)
!===============================================================================

rec1 =  0.1d0
rec2 =  1.0d0
rec3 =  10.d0
rec4 = 200.d0

dd2 = ettp(ip,jdp) * ettp(ip,jdp)

if ( rep.le.rec1 ) then
  fdr = 18.d0 * xnul / dd2

else if ( rep.le.rec2 ) then
  fdr = 3.d0/4.d0 * xnul / dd2                                     &
      * (22.73d0 + 0.0903d0/rep + 3.69d0*rep )

else if ( rep.le.rec3 ) then
  fdr = 3.d0/4.d0 * xnul / dd2                                     &
      * (29.1667d0 - 3.8889d0/rep + 1.222d0*rep)

else if ( rep.le.rec4 ) then
    fdr = 18.d0*xnul/dd2 *(1.d0 + 0.15d0*rep**0.687d0)

else
   fdr = (0.44d0 * 3.d0 / 4.d0) * uvwr / ettp(ip,jdp)
endif

taup = romp / romf / fdr


!==============================================================================

!--------
! Formats
!--------


!----
! End
!----

end subroutine
