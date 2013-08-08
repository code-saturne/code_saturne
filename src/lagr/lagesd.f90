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

subroutine lagesd &
!================

 ( ifac   , ip     ,                                              &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   dlgeo  ,                                                       &
   dt     , rtpa   , propce ,                                     &
   ettp   , ettpa  , tepa   ,                                     &
   statis , taup   , tlag   , piil   ,                            &
   bx     , vagaus , gradpr , gradvf , romp,                      &
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
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nnod             ! e  ! <-- ! nombre de sommets                              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ifabor           ! e  ! <-- !                                                !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
! dlgeo            ! tr ! --> ! tableau contenant les donnees geometriques     !
!(nfabor,ngeol)    !    !     !                                                !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (pas de temps precedent)           !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! ettp             ! tr ! --> ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! tr ! <-- ! cumul des statistiques volumiques              !
!(ncelet,nvlsta    !    !     !                                                !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! piil(nbpmax,3    ! tr ! <-- ! terme dans l'integration des eds up            !
! bx(nbpmax,3,2    ! tr ! <-- ! caracteristiques de la turbulence              !
! vagaus           ! tr ! <-- ! variables aleatoires gaussiennes               !
!(nbpmax,nvgaus    !    !     !                                                !
! gradpr(ncel,3    ! tr ! <-- ! gradient de pression                           !
! gradvf(ncel,3    ! tr ! <-- ! gradient de la vitesse du fluide               !
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
use ppincl
use pointe

!===============================================================================

implicit none

! Arguments

integer          ifac   , ip
integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          itepa(nbpmax,nivep)

double precision dlgeo(nfabor,ngeol)
double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep) , statis(ncelet,*)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3) , bx(nbpmax,3,2)
double precision vagaus(nbpmax,*)
double precision gradpr(ncelet,3) , gradvf(ncelet,9)
double precision romp(nbpmax)
double precision ustar

! Local variables

integer          isens  , iel , mode, id, i0
double precision depg(3),depl(3),vpart(3),vvue(3),tempf, romf
double precision vflui(3),vflui1(3),vdirn,vdirt,vdirtt,vpartl,vvuel,dxl, enertur
double precision ggp(3), gdpr(3), piilp(3), tlp,bxp

double precision aa , bb , cc , dd , ee
double precision aux1 , aux2 ,aux3 , aux4 , aux5 , aux6
double precision aux7 , aux8 , aux9 , aux10 , aux11
double precision ter1f , ter2f , ter3f
double precision ter1p , ter2p , ter3p , ter4p , ter5p
double precision ter1x , ter2x , ter3x , ter4x , ter5x
double precision tci , force, k1
double precision gama2 , omegam , omega2
double precision grga2 , gagam , gaome
double precision p11 , p21 , p22 , p31 , p32 , p33
double precision grav(3)
double precision lvisq, tvisq, depint
double precision c0, cl, visccf
double precision energi , dissip , vit(3)
double precision norm_vit , norm
integer          iromf

! Local variables for the resuspension model

double precision drag(3)                 ! Hydrodynamic drag on a deposited particle
double precision tordrg(3), tordrg_norm  ! Hydrodynamic torque on a deposited particle

double precision omep(3)  ! Angular velocity of the rolling particle
double precision dome(3)  ! Time increment of the angular velocity of the deposited particle

double precision scalax

double precision iner_tor,  cst_1,  cst_4, adh_tor(3),vpart0(3)
double precision kk, kkk

!===============================================================================

!===============================================================================
! 0.  Memory management and Initialization
!===============================================================================

! Initializations to avoid compiler warning
energi = 0.d0
dissip = 0.d0
norm_vit = 0.d0

iel = itepa(ip,jisor)

! Friction velocity
ifac = itepa(ip,jdfac)
ustar = uetbor(ifac)

! Constants for the calculation of bxp and tlp
c0   = 2.1d0
cl   = 1.d0 / (0.5d0 + (3.d0/4.d0)*c0)

! Pointer on the density w.r.t the flow
if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then
  iromf = ipproc(irom1)
else
  iromf = ipproc(irom)
endif

romf = propce(iel,iromf)
visccf = propce(iel,ipproc(iviscl)) / romf

norm=sqrt(rtpa(iel,iu)**2 + rtpa(iel,iv)**2 + rtpa(iel,iw)**2)

! Velocity norm w.r.t y+
if (tepa(ip,jryplu).le.5.d0) then
   norm_vit = tepa(ip,jryplu) * ustar
else if ( tepa(ip,jryplu).gt.5.d0.and.tepa(ip,jryplu).le.30.d0) then
   norm_vit = ( -3.05d0 + 5.d0 * log(tepa(ip,jryplu))) * ustar
else if (tepa(ip,jryplu).gt.30.d0.and.tepa(ip,jryplu).lt.100.d0) then
   norm_vit = (2.5d0 * log(tepa(ip,jryplu)) + 5.5d0) * ustar
endif

vit(1) = norm_vit * rtpa(iel,iu) / norm
vit(2) = norm_vit * rtpa(iel,iv) / norm
vit(3) = norm_vit * rtpa(iel,iw) / norm

! Turbulent kinetic energy and dissipation w.r.t y+
if (tepa(ip,jryplu).le.5.d0) then
   energi = 0.1d0 * (tepa(ip,jryplu)**2) * ustar**2
   dissip = 0.2d0 * ustar**4 / visccf
else if (tepa(ip,jryplu).gt.5.d0.and.tepa(ip,jryplu).le.30.d0) then
   energi = ustar**2 / (0.09d0)**0.5
   dissip = 0.2d0 * ustar**4 / visccf
else if (tepa(ip,jryplu).gt.30.d0.and.tepa(ip,jryplu).le.100.d0) then
   energi = ustar**2 / (0.09d0)**0.5
   dissip = ustar**4 / (0.41d0 * tepa(ip,jryplu) * visccf)
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
       ettpa(ip,jup)  , ettpa(ip,jvp)  , ettpa(ip,jwp)  ,           &
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
       ettpa(ip,juf)  , ettpa(ip,jvf)  , ettpa(ip,jwf)  ,           &
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
      rtpa(iel,iu)   , rtpa(iel,iv)   , rtpa(iel,iw)   ,            &
      vflui(1)       , vflui1(2)      , vflui1(3)       ,            &
      dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,            &
      dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,            &
      dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )

! 2.5 - pressure gradient

call  lagprj                                                        &
     !===========
     ( isens          ,                                              &
     gradpr(iel,1)  , gradpr(iel,2)  , gradpr(iel,3)  ,            &
     gdpr(1)        , gdpr(2)        , gdpr(3)        ,            &
     dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,            &
     dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,            &
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

tlp = cl * energi / dissip
tlp = max(tlp,epzero)

! 2.8 - bx

bxp = sqrt(c0*dissip)

!===============================================================================
! 3. Integration of the EDS on the particles
!===============================================================================

!  Retrieve of the turbulent kinetic energy

if (itytur.eq.2 .or. iturb.eq.50 .or. iturb.eq.60) then
  enertur = rtpa(iel,ik)
else if (itytur.eq.3) then
  enertur = 0.5d0*( rtpa(iel,ir11)                         &
                  + rtpa(iel,ir22)                         &
                  + rtpa(iel,ir33) )
endif

call lagcli                                                       &
!==========
   ( itepa(ip,jimark),                                            &
     tempf        ,                                               &
     romf, ustar, lvisq, tvisq,                                   &
     vpart(1)     , vvue(1)   , depl(1) ,                         &
     ettp(ip,jdp) , romp(ip)  , taup(ip),                         &
     tepa(ip,jryplu),tepa(ip,jrinpf), enertur, ggp(1), vflui(1),  &
     gdpr(1), piilp(1), depint )

if (itepa(ip,jdepo).gt.0) then
   depl(1) = 0.d0
   vpart(1) = 0.d0
endif


!  Integration in the 2 other directions

if (itepa(ip,jdepo).eq.0) then

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

   if (itepa(ip,jdepo).gt.0) then

      ! Resuspension model

      ! Calculation of the hydrodynamic drag and torque
      ! applied on the deposited particle

      drag(1) = 3.d0 * pi * ettp(ip,jdp) * (vvue(1) - vpart(1)) * visccf * romf * 3.39d0
      tordrg(1) = 0.0d0

      do id = 2,3
         drag(id) = 3.d0 * pi * ettp(ip,jdp) * (vvue(id)-vpart(id)) * visccf * romf * 1.7d0
         tordrg(id) = 1.4d0 * drag(id) * ettp(ip,jdp) * 0.5d0
      enddo


      ! Is there direct wall-normal lift-off of the particle ?

      if ((abs(drag(1)).gt.tepa(ip,jfadh)).and.(drag(1).lt.0.d0)) then

         ! The particle is resuspended

         itepa(ip,jdepo) = 0

         tepa(ip,jfadh) = 0.d0
         tepa(ip,jmfadh) = 0.d0

         itepa(ip,jnbasg) = 0
         itepa(ip,jnbasp) = 0

         tepa(ip,jndisp) = 0.d0

         vpart(1) = min(- 1.d0 / ettp(ip,jmp) * abs(drag(1) -  tepa(ip,jfadh)) * dtp, 1.d-3)
         vpart(2) = 0.d0
         vpart(3) = 0.d0

         ! Update of the number and weight of resuspended particles

         nbpres = nbpres + 1
         dnbres = dnbres + tepa(ip, jrpoi)

      else  ! No direct normal lift-off

         ! Calculation of the norm of the hydrodynamic torque and drag (tangential)

         tordrg_norm = sqrt(tordrg(2)**2 + tordrg(3)**2)

         adh_tor(2) = - tepa(ip,jmfadh) / tordrg_norm * tordrg(2)
         adh_tor(3) = - tepa(ip,jmfadh) / tordrg_norm * tordrg(3)

         do id = 2,3

            iner_tor = (7.d0/5.d0)*ettp(ip,jmp)*(ettp(ip,jdp) * 0.5d0)**2

            cst_4 = 6 * pi * visccf * romf * 1.7 * 1.4 * (ettp(ip,jdp) * 0.5d0)**2

            cst_1 = cst_4 * (ettp(ip,jdp) * 0.5d0) / iner_tor

            vpart0(id) = vpart(id)

            vpart(id) = (vpart0(id) - vvue(id) - adh_tor(id) / cst_4)* exp( - cst_1 * dtp) &
                 + (vvue(id) + adh_tor(id) / cst_4)

         enddo

         scalax = vpart(2) * vvue(2) + vpart(3) * vvue(3)

         if (scalax.gt.0.d0) then

            ! The calculated particle velocity is aligned
            ! with the flow seen
            ! --> The particle starts or keep on rolling

            itepa(ip,jdepo) = 2

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

            itepa(ip,jdepo) = 10
            do id = 2,3
               depl(id) = 0.d0
               vpart(id) = 0.d0
            enddo

         endif ! if (scalax..)

      endif

   endif

else  ! if ireent.eq.0 --> Motionless deposited particle

   if (itepa(ip,jdepo).gt.0) then

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
      ettp(ip,jup)   , ettp(ip,jvp)   , ettp(ip,jwp)   ,          &
      vpart(1)       , vpart(2)       , vpart(3)       ,          &
      dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,          &
      dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,          &
      dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )



! 3.3 - flow-seen velocity

call lagprj                                                       &
!==========
    ( isens          ,                                            &
      ettp(ip,juf)   , ettp(ip,jvf)   , ettp(ip,jwf)   ,          &
      vvue(1)        , vvue(2)        , vvue(3)        ,          &
      dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,          &
      dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,          &
      dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )

!===============================================================================
! 5. Computation of the new particle position
!===============================================================================

   ettp(ip,jxp)=ettp(ip,jxp)+depg(1)
   ettp(ip,jyp)=ettp(ip,jyp)+depg(2)
   ettp(ip,jzp)=ettp(ip,jzp)+depg(3)


!===============================================================================

end subroutine
