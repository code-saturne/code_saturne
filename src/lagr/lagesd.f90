!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine lagesd &
!================

 ( ifac   , ip     ,                                              &
   taup   , piil   ,                                              &
   vagaus , gradpr , romp,                                        &
   tempf  , romf   , ustar  , lvisq  ,tvisq   ,  depint )

!===============================================================================

! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module:
!   ------------------------------------------------------


!   Deposition submodel:
!
!    1/ Modification of the coordinate system (global ->local)
!    2/ Call of subroutine lagcli
!    3/ Integration of the stochastic differential equations
!    in the 2 directions different from the normal to the boundary face
!    4/ Modification of the coordinate system (local ->global)
!    5/ Update of the particle position
!
!
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ifac             ! e  ! <-- !                                                !
! ip               ! e  ! <-- !                                                !
! taup(nbpart)     ! tr ! <-- ! temps caracteristique dynamique                !
! piil(nbpart,3)   ! tr ! <-- ! terme dans l'integration des eds up            !
! vagaus           ! tr ! <-- ! variables aleatoires gaussiennes               !
!(nbpart,nvgaus)   !    !     !                                                !
! gradpr(3,ncel)   ! tr ! <-- ! gradient de pression                           !
! romp             ! tr ! <-- ! masse volumique des particules                 !
! tempf            !  r ! <-- ! temperature of the fluid (K)                   !
! romf             !  r ! <-- ! density of the fluid                           !
! ustar            !  r ! <-- ! friction velocity                              !
! lvisq            !  r ! <-- ! wall-unit lengthscale                          !
! tvisq            !  r ! <-- ! wall-unit timescale                            !
! depint           !  r ! <-- ! interface location near-wall/core-flow         !
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
use field
use ppincl
use pointe

!===============================================================================

implicit none

! Arguments

integer          ifac   , ip

double precision taup(nbpart)
double precision piil(nbpart,3)
double precision vagaus(nbpart,*)
double precision gradpr(3,ncelet)
double precision romp(nbpart)
double precision ustar

! Local variables

integer          isens  , iel , id, i0
double precision depg(3),depl(3),vpart(3),vvue(3),tempf, romf
double precision vflui(3),vflui1(3),enertur
double precision ggp(3), gdpr(3), piilp(3), tlp,bxp

double precision aa , bb , cc , dd , ee
double precision aux1 , aux2 ,aux3 , aux4 , aux5 , aux6
double precision aux7 , aux8 , aux9 , aux10 , aux11
double precision ter1f , ter2f , ter3f
double precision ter1p , ter2p , ter3p , ter4p , ter5p
double precision ter1x , ter2x , ter3x , ter4x , ter5x
double precision tci , force
double precision gama2 , omegam , omega2
double precision grga2 , gagam , gaome
double precision p11 , p21 , p22 , p31 , p32 , p33
double precision lvisq, tvisq, depint
double precision c0, cl, visccf
double precision energi , dissip , vit(3)
double precision norm_vit , norm

! Local variables for the resuspension model

double precision drag(3)                 ! Hydrodynamic drag on a deposited particle
double precision tordrg(3), tordrg_norm  ! Hydrodynamic torque on a deposited particle

double precision scalax

double precision iner_tor,  cst_1,  cst_4, adh_tor(3),vpart0(3)
double precision kk, kkk

double precision, dimension(:), pointer :: cromf
double precision, dimension(:,:), pointer :: vela
double precision, dimension(:), pointer :: cvara_k
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: viscl

!===============================================================================

! Map field arrays
call field_get_val_prev_v(ivarfl(iu), vela)

if (itytur.eq.2 .or. iturb.eq.50 .or. iturb.eq.60) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
else if (itytur.eq.3) then
  call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
  call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
  call field_get_val_prev_s(ivarfl(ir33), cvara_r33)
endif

!===============================================================================
! 0.  Memory management and Initialization
!===============================================================================

! Initializations to avoid compiler warning
energi = 0.d0
dissip = 0.d0
norm_vit = 0.d0

iel = ipepa(jisor,ip)

! Friction velocity
ifac = ipepa(jdfac,ip)
ustar = uetbor(ifac)

! Constants for the calculation of bxp and tlp
c0   = 2.1d0
cl   = 1.d0 / (0.5d0 + (3.d0/4.d0)*c0)

! Pointer on the density w.r.t the flow

if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then
  call field_get_val_s(iprpfl(ipproc(irom1)), cromf)
else
  call field_get_val_s(icrom, cromf)
endif

call field_get_val_s(iprpfl(iviscl), viscl)

romf = cromf(iel)
visccf = viscl(iel) / romf

norm=sqrt(vela(1,iel)**2 + vela(2,iel)**2 + vela(3,iel)**2)

! Velocity norm w.r.t y+
if (pepa(jryplu,ip).le.5.d0) then
   norm_vit = pepa(jryplu,ip) * ustar
else if ( pepa(jryplu,ip).gt.5.d0.and.pepa(jryplu,ip).le.30.d0) then
   norm_vit = ( -3.05d0 + 5.d0 * log(pepa(jryplu,ip))) * ustar
else if (pepa(jryplu,ip).gt.30.d0.and.pepa(jryplu,ip).lt.100.d0) then
   norm_vit = (2.5d0 * log(pepa(jryplu,ip)) + 5.5d0) * ustar
endif

if (ustar.gt.0.d0) then
  vit(1) = norm_vit * vela(1,iel) / norm
  vit(2) = norm_vit * vela(2,iel) / norm
  vit(3) = norm_vit * vela(3,iel) / norm
else
  vit(1) = 0.d0
  vit(2) = 0.d0
  vit(3) = 0.d0
endif

! Turbulent kinetic energy and dissipation w.r.t y+
if (pepa(jryplu,ip).le.5.d0) then
   energi = 0.1d0 * (pepa(jryplu,ip)**2) * ustar**2
   dissip = 0.2d0 * ustar**4 / visccf
else if (pepa(jryplu,ip).gt.5.d0.and.pepa(jryplu,ip).le.30.d0) then
   energi = ustar**2 / (0.09d0)**0.5
   dissip = 0.2d0 * ustar**4 / visccf
else if (pepa(jryplu,ip).gt.30.d0.and.pepa(jryplu,ip).le.100.d0) then
   energi = ustar**2 / (0.09d0)**0.5
   dissip = ustar**4 / (0.41d0 * pepa(jryplu,ip) * visccf)
endif

!===============================================================================
! 2. Reference frame change:
!---------------------------
! global reference frame --> local reference frame for the boundary face
!===============================================================================

isens = 1

! 2.1 - particle velocity

call  lagprj                                                        &
!===========
     ( isens          ,                                             &
       eptpa(jup,ip)  , eptpa(jvp,ip)  , eptpa(jwp,ip)  ,           &
       vpart(1)       , vpart(2)       , vpart(3)       ,           &
       dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,           &
       dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,           &
       dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )

vpart0(2) = vpart(2)
vpart0(3) = vpart(3)

! 2.2 - flow-seen velocity

call  lagprj                                                        &
!===========
     ( isens          ,                                             &
       eptpa(juf,ip)  , eptpa(jvf,ip)  , eptpa(jwf,ip)  ,           &
       vvue(1)        , vvue(2)        , vvue(3)        ,           &
       dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,           &
       dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,           &
       dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )

! 2.3 - Gravity vector

call  lagprj                                                        &
!===========
    ( isens  ,                                                      &
      gx     , gy     , gz     ,                                    &
      ggp(1) , ggp(2) , ggp(3) ,                                    &
      dlgeo(ifac, 5)  , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,           &
      dlgeo(ifac, 8)  , dlgeo(ifac, 9) , dlgeo(ifac,10) ,           &
      dlgeo(ifac,11)  , dlgeo(ifac,12) , dlgeo(ifac,13) )

! 2.4 - flow velocity

call  lagprj                                                        &
!===========
    ( isens                 ,                                       &
      vit(1)         , vit(2)         , vit(3)         ,            &
      vflui1(1)      , vflui(2)       , vflui(3)       ,            &
      dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,            &
      dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,            &
      dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )

call  lagprj                                                        &
!===========
    ( isens                 ,                                       &
      vela(1,iel)    , vela(2,iel)    , vela(3,iel)    ,            &
      vflui(1)       , vflui1(2)      , vflui1(3)      ,            &
      dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,            &
      dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,            &
      dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )

! 2.5 - pressure gradient

call  lagprj                                                        &
!===========
   ( isens          ,                                               &
     gradpr(1,iel)  , gradpr(2,iel)  , gradpr(3,iel)  ,             &
     gdpr(1)        , gdpr(2)        , gdpr(3)        ,             &
     dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,             &
     dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,             &
     dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13) )

! 2.6 - "piil" term

call  lagprj                                                        &
     !===========
     ( isens          ,                                              &
     piil(ip,1)     , piil(ip,2)     , piil(ip,3)     ,            &
     piilp(1)       , piilp(2)       , piilp(3)       ,            &
     dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,            &
     dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,            &
     dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13) )

! 2.7 - tlag

if (energi.gt.0.d0) then
  tlp = cl * energi / dissip
  tlp = max(tlp,epzero)
else
  tlp = epzero
endif

! 2.8 - bx

bxp = sqrt(c0*dissip)

!===============================================================================
! 3. Integration of the EDS on the particles
!===============================================================================

!  Retrieve of the turbulent kinetic energy

if (itytur.eq.2 .or. iturb.eq.50 .or. iturb.eq.60) then
  enertur = cvara_k(iel)
else if (itytur.eq.3) then
  enertur = 0.5d0*( cvara_r11(iel)                         &
                  + cvara_r22(iel)                         &
                  + cvara_r33(iel) )
endif

call lagcli                                                       &
!==========
   ( ipepa(jimark,ip),                                            &
     tempf        ,                                               &
     lvisq, tvisq,                                                &
     vpart(1)     , vvue(1)   , depl(1) ,                         &
     eptp(jdp,ip) , romp(ip)  , taup(ip),                         &
     pepa(jryplu,ip),pepa(jrinpf,ip), enertur, ggp(1), vflui(1),  &
     gdpr(1), piilp(1), depint )

if (ipepa(jdepo,ip).gt.0) then
   depl(1) = 0.d0
   vpart(1) = 0.d0
endif


!  Integration in the 2 other directions

if (ipepa(jdepo,ip).eq.0) then

   do id = 2,3

      i0 = id - 1

      tci = piilp(id) * tlp + vflui(id)
      force = ( romf * gdpr(id) / romp(ip) + ggp(id) ) * taup(ip)

      aux1 = exp( -dtp / taup(ip))
      aux2 = exp( -dtp / tlp )
      aux3 = tlp / (tlp-taup(ip))

      aux4 = tlp / (tlp+taup(ip))
      aux5 = tlp * (1.d0-aux2)
      aux6 = bxp * bxp * tlp

      aux7 = tlp - taup(ip)
      aux8 = bxp * bxp * aux3**2

      !---> Terms for the trajectory

      aa = taup(ip) * (1.d0 - aux1)
      bb = (aux5 - aa) * aux3
      cc = dtp - aa - bb

      ter1x = aa * vpart(id)
      ter2x = bb * vvue(id)
      ter3x = cc * tci
      ter4x = (dtp - aa) * force

      !---> Terms for the flow-seen velocity

      ter1f = vvue(id) * aux2
      ter2f = tci * (1.d0-aux2)

      !---> Terms for the particles velocity

      dd = aux3 * (aux2 - aux1)
      ee = 1.d0 - aux1

      ter1p = vpart(id) * aux1
      ter2p = vvue(id) * dd
      ter3p = tci * (ee-dd)
      ter4p = force * ee

      !---> (2.3) Coefficients computation for the stochastic integrals:

      gama2  = 0.5d0 * (1.d0 - aux2*aux2 )
      omegam = 0.5d0 * aux4 * ( aux5 - aux2*aa )                      &
           -0.5d0 * aux2 * bb
      omegam = omegam * sqrt(aux6)

      omega2 = aux7                                                   &
           * ( aux7*dtp - 2.d0 * (tlp*aux5-taup(ip)*aa) )      &
           + 0.5d0 * tlp * tlp*aux5 * (1.d0 + aux2)        &
           + 0.5d0 * taup(ip) * taup(ip) * aa * (1.d0+aux1)        &
           - 2.0d0 * aux4 * tlp * taup(ip) * taup(ip)          &
           * (1.d0 - aux1*aux2)

      omega2 = aux8 * omega2

      if (abs(gama2).gt.epzero) then

         p21 = omegam / sqrt(gama2)
         p22 = omega2 - p21**2
         p22 = sqrt( max(zero,p22) )

      else
         p21 = 0.d0
         p22 = 0.d0
      endif

      ter5x = p21 * vagaus(ip,i0) + p22 * vagaus(ip,3+i0)

      !---> Integral for the flow-seen velocity

      p11 = sqrt( gama2*aux6 )
      ter3f = p11*vagaus(ip,i0)

      !---> Integral for particles velocity

      aux9  = 0.5d0 * tlp * (1.d0 - aux2*aux2)
      aux10 = 0.5d0 * taup(ip) * (1.d0 - aux1*aux1)
      aux11 = taup(ip) * tlp * (1.d0 - aux1*aux2)                 &
           / (taup(ip) + tlp)

      grga2 = (aux9 - 2.d0*aux11 + aux10) * aux8
      gagam = (aux9 - aux11) * (aux8 / aux3)
      gaome = ( (tlp - taup(ip)) * (aux5 - aa)                    &
           -  tlp * aux9 - taup(ip) * aux10                    &
           + (tlp + taup(ip)) * aux11 ) * aux8

      if(p11.gt.epzero) then
         p31 = gagam / p11
      else
         p31 = 0.d0
      endif

      if(p22.gt.epzero) then
         p32 = (gaome-p31*p21) / p22
      else
         p32 = 0.d0
      endif

      p33 = grga2 - p31**2 - p32**2
      p33 = sqrt( max(zero,p33) )

      ter5p =  p31 * vagaus(ip,i0)                                   &
           + p32 * vagaus(ip,3+i0)                                 &
           + p33 * vagaus(ip,6+i0)

      !---> trajectory

      depl(id) = ter1x + ter2x + ter3x + ter4x + ter5x

      !---> flow-seen velocity

      vvue(id) = ter1f + ter2f + ter3f

      !---> particles velocity
      vpart(id) = ter1p + ter2p + ter3p + ter4p + ter5p

   enddo

else

   do id = 2,3

      i0 = id - 1

      tci = piilp(id) * tlp + vflui(id)
      force = ( romf * gdpr(id) / romp(ip) + ggp(id) ) * taup(ip)

      aux1 = exp( -dtp / taup(ip))
      aux2 = exp( -dtp / tlp )
      aux3 = tlp / (tlp-taup(ip))

      aux4 = tlp / (tlp+taup(ip))
      aux5 = tlp * (1.d0-aux2)
      aux6 = bxp * bxp * tlp


      !---> Terms for the flow-seen velocity

      ter1f = vvue(id) * aux2
      ter2f = tci * (1.d0-aux2)

      !---> (2.3) Coefficients computation for the stochastic integrals:

      gama2  = 0.5d0 * (1.d0 - aux2*aux2 )

      !---> Integral for the flow-seen velocity

      p11 = sqrt( gama2*aux6 )
      ter3f = p11*vagaus(ip,i0)

      !---> flow-seen velocity

      vvue(id) = ter1f + ter2f + ter3f

   enddo

endif


if (ireent.eq.1) then

   if (ipepa(jdepo,ip).gt.0) then

      ! Resuspension model

      ! Calculation of the hydrodynamic drag and torque
      ! applied on the deposited particle

      drag(1) = 3.d0 * pi * eptp(jdp,ip) * (vvue(1) - vpart(1)) * visccf * romf * 3.39d0
      tordrg(1) = 0.0d0

      do id = 2,3
         drag(id) = 3.d0 * pi * eptp(jdp,ip) * (vvue(id)-vpart(id)) * visccf * romf * 1.7d0
         tordrg(id) = 1.4d0 * drag(id) * eptp(jdp,ip) * 0.5d0
      enddo


      ! Is there direct wall-normal lift-off of the particle ?

      if ((abs(drag(1)).gt.pepa(jfadh,ip)).and.(drag(1).lt.0.d0)) then

         ! The particle is resuspended

         ipepa(jdepo,ip) = 0

         pepa(jfadh,ip) = 0.d0
         pepa(jmfadh,ip) = 0.d0

         ipepa(jnbasg,ip) = 0
         ipepa(jnbasp,ip) = 0

         pepa(jndisp,ip) = 0.d0

         vpart(1) = min(- 1.d0 / eptp(jmp,ip) * abs(drag(1) -  pepa(jfadh,ip)) * dtp, 1.d-3)
         vpart(2) = 0.d0
         vpart(3) = 0.d0

         ! Update of the number and weight of resuspended particles

         nbpres = nbpres + 1
         dnbres = dnbres + pepa( jrpoi,ip)

         parbor(ipepa(jdfac,ip),ires) = parbor(ipepa(jdfac,ip),ires) + pepa(jrpoi,ip)

         parbor(ipepa(jdfac,ip),iflres) = parbor(ipepa(jdfac,ip),iflres) + pepa(jrpoi,ip)  &
                       + ( pepa(jrpoi,ip) * eptp(jmp,ip) / surfbn(ipepa(jdfac,ip)))

         parbor(ipepa(jdfac,ip),iflm) = parbor(ipepa(jdfac,ip),iflm)                       &
                       - ( pepa(jrpoi,ip) * eptp(jmp,ip) / surfbn(ipepa(jdfac,ip)))


         else  ! No direct normal lift-off

         ! Calculation of the norm of the hydrodynamic torque and drag (tangential)

         tordrg_norm = sqrt(tordrg(2)**2 + tordrg(3)**2)

         adh_tor(2) = - pepa(jmfadh,ip) / tordrg_norm * tordrg(2)
         adh_tor(3) = - pepa(jmfadh,ip) / tordrg_norm * tordrg(3)

         do id = 2,3

            iner_tor = (7.d0/5.d0)*eptp(jmp,ip)*(eptp(jdp,ip) * 0.5d0)**2

            cst_4 = 6 * pi * visccf * romf * 1.7 * 1.4 * (eptp(jdp,ip) * 0.5d0)**2

            cst_1 = cst_4 * (eptp(jdp,ip) * 0.5d0) / iner_tor

            vpart0(id) = vpart(id)

            vpart(id) = (vpart0(id) - vvue(id) - adh_tor(id) / cst_4)* exp( - cst_1 * dtp) &
                 + (vvue(id) + adh_tor(id) / cst_4)

         enddo

         scalax = vpart(2) * vvue(2) + vpart(3) * vvue(3)

         if (scalax.gt.0.d0) then

            ! The calculated particle velocity is aligned
            ! with the flow seen
            ! --> The particle starts or keep on rolling

            ipepa(jdepo,ip) = 2

            vpart(1) = 0.0d0

            do id = 2,3

               if (abs(vpart(id)).gt.abs(vvue(id))) then

                  ! The velocity of the rolling particle cannot
                  ! exceed the surrounding fluid velocity

                  vpart(id) = vvue(id)

               endif

               kk = vpart0(id) - vvue(id) - adh_tor(id) / cst_4
               kkk = vvue(id) + adh_tor(id) / cst_4

              depl(id) = ((cst_1 * kkk * dtp + kk) * exp( cst_1 * dtp) - kk) * &
                             exp(- cst_1 * dtp) / cst_1
            enddo

         else

            ! The particle is not set into motion or stops
            ! the flag is set to  10 and velocity and displacement are null

            ipepa(jdepo,ip) = 10
            do id = 2,3
               depl(id) = 0.d0
               vpart(id) = 0.d0
            enddo

         endif ! if (scalax..)

      endif

   endif

else  ! if ireent.eq.0 --> Motionless deposited particle

   if (ipepa(jdepo,ip).gt.0) then

      do id = 2,3
         vpart(id) = 0.d0
         vvue(id) = 0.d0
         depl(id) = 0.d0
      enddo

   endif

endif



!===============================================================================
! 3. Reference frame change:
!---------------------------
!local reference frame for the boundary face  --> global reference frame
!===============================================================================

isens = 2

! 3.1 - Displacement
!
call lagprj                                                       &
!==========
    ( isens   ,                                                   &
      depg(1) , depg(2) , depg(3)     ,                           &
      depl(1) , depl(2) , depl(3)     ,                           &
      dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,          &
      dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,          &
      dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )

! 3.2 - Particle velocity
!
call lagprj                                                       &
!==========
    ( isens          ,                                            &
      eptp(jup,ip)   , eptp(jvp,ip)   , eptp(jwp,ip)   ,          &
      vpart(1)       , vpart(2)       , vpart(3)       ,          &
      dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,          &
      dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,          &
      dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )



! 3.3 - flow-seen velocity

call lagprj                                                       &
!==========
    ( isens          ,                                            &
      eptp(juf,ip)   , eptp(jvf,ip)   , eptp(jwf,ip)   ,          &
      vvue(1)        , vvue(2)        , vvue(3)        ,          &
      dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,          &
      dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,          &
      dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )

!===============================================================================
! 5. Computation of the new particle position
!===============================================================================

   eptp(jxp,ip) = eptp(jxp,ip) + depg(1)
   eptp(jyp,ip) = eptp(jyp,ip) + depg(2)
   eptp(jzp,ip) = eptp(jzp,ip) + depg(3)


!===============================================================================

return
end subroutine lagesd
