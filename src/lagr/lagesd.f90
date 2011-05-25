!-------------------------------------------------------------------------------

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

subroutine lagesd &
!================

 ( idbia0 , idbra0 ,                                              &
   ifac   , ip     ,                                              &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  , ia     ,                                              &
   dlgeo  ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   ,                                     &
   statis , taup   , tlag   , piil   ,                            &
   bx     , vagaus , gradpr , gradvf , romp,                      &
   tempf  , romf   , ustar  , lvisq  ,tvisq   ,  depint ,         &
   ra )


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
! idbia0           ! e  ! <-- ! number of first free position in ia            !
! idbra0           ! e  ! <-- ! number of first free position in ra            !
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
! ia(*)            ! tr ! --- ! macro tableau entier                           !
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
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
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
! ra(*)            ! tr ! --- ! macro tableau reel                             !
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
use ppppar
use ppthch
use ppincl
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ifac   , ip
integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          itepa(nbpmax,nivep)  , ia(*)

double precision dlgeo(nfabor,ngeol)
double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep) , statis(ncelet,*)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3) , bx(nbpmax,3,2)
double precision vagaus(nbpmax,*)
double precision gradpr(ncelet,3) , gradvf(ncelet,9)
double precision romp(nbpmax)
double precision ra(*)

! Local variables

integer          idebia , idebra
integer          isens  , iel , mode, id, i0
double precision depg(3),depl(3),vpart(3),vvue(3),tempf, romf
double precision vflui(3),vdirn,vdirt,vdirtt,vpartl,vvuel,dxl, enertur
double precision ggp(3), gdpr(3), piilp(3), tlp(3),bxp(3)

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
double precision grav(3) , ad1, ad2
double precision ustar, lvisq, tvisq, depint

!===============================================================================

!===============================================================================
! 0.  Memory management and Initialization
!===============================================================================

idebia = idbia0
idebra = idbra0

iel = itepa(ip,jisor)


!===============================================================================
! 2. Reference frame change:
!---------------------------
!global reference frame --> in the local reference frame for the boundary face
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

! 2.2 - vu fluid velocity

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

! 2.4 - fluid velocity

call  lagprj                                                        &
!===========
    ( isens                 ,                                       &
      rtpa(iel,iu)   , rtpa(iel,iv)   , rtpa(iel,iw)   ,            &
      vflui(1)       , vflui(2)       , vflui(3)       ,            &
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

call  lagprj                                                        &
!===========
    ( isens          ,                                              &
      tlag(ip,1)     , tlag(ip,2)     , tlag(ip,3)     ,            &
      tlp(1)         , tlp(2)         , tlp(3)         ,            &
      dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,            &
      dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,            &
      dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13) )

tlp(1)=abs(tlp(1))
tlp(2)=abs(tlp(2))
tlp(3)=abs(tlp(3))

if ((abs(tlag(ip,1)-tlag(ip,2)) .le. 1.d-30) .and.                 &
    (abs(tlag(ip,2)-tlag(ip,3)) .le. 1.d-30)) then
  tlp(1)= tlag(ip,1)
  tlp(2)= tlag(ip,1)
  tlp(3)= tlag(ip,1)
endif

! 2.8 - bx

call  lagprj                                                        &
!===========
    ( isens          ,                                              &
      bx(ip,1,nor)   , bx(ip,2,nor)   , bx(ip,3,nor)   ,            &
      bxp(1)         , bxp(2)         , bxp(3)         ,            &
      dlgeo(ifac, 5) , dlgeo(ifac, 6) , dlgeo(ifac, 7) ,            &
      dlgeo(ifac, 8) , dlgeo(ifac, 9) , dlgeo(ifac,10) ,            &
      dlgeo(ifac,11) , dlgeo(ifac,12) , dlgeo(ifac,13)  )

if ((abs(bx(ip,1,nor)-bx(ip,2,nor)) .le. 1.d-30) .and.             &
    (abs(bx(ip,2,nor)-bx(ip,3,nor)) .le. 1.d-30)) then
  bxp(1)= bx(ip,1,nor)
  bxp(2)= bx(ip,1,nor)
  bxp(3)= bx(ip,1,nor)
endif



!===============================================================================
! 3. Integration of the EDS on the particles
!===============================================================================

!  Recovery of the turbulent kinetic energy

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
     romf, ustar, lvisq, tvisq, ifac,                             &
     vpart(1)     , vvue(1)   , depl(1) ,                         &
     ettp(ip,jdp) , romp(ip)  , taup(ip),                         &
     tepa(ip,jryplu),tepa(ip,jrinpf), enertur, ggp(1), vflui(1),  &
     gdpr(1), piilp(1), depint, ra)

!  Integration in the 2 others directions

do id = 2,3

  i0 = id - 1

  tci = piilp(id) * tlp(id) + vflui(id)
  force = ( romf * gdpr(id) / romp(ip) + ggp(id) ) * taup(ip)

  aux1 = exp( -dtp / taup(ip))
  aux2 = exp( -dtp / tlp(id) )
  aux3 = tlp(id) / (tlp(id)-taup(ip))

  aux4 = tlp(id) / (tlp(id)+taup(ip))
  aux5 = tlp(id) * (1.d0-aux2)
  aux6 = bxp(id) * bxp(id) * tlp(id)

  aux7 = tlp(id) - taup(ip)
  aux8 = bxp(id) * bxp(id) * aux3**2

!---> Terms for the trajectory

  aa = taup(ip) * (1.d0 - aux1)
  bb = (aux5 - aa) * aux3
  cc = dtp - aa - bb

  ter1x = aa * vpart(id)
  ter2x = bb * vvue(id)
  ter3x = cc * tci
  ter4x = (dtp - aa) * force

!---> Terms for the vu fluid

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
          * ( aux7*dtp - 2.d0 * (tlp(id)*aux5-taup(ip)*aa) )      &
          + 0.5d0 * tlp (id)* tlp(id)*aux5 * (1.d0 + aux2)        &
          + 0.5d0 * taup(ip) * taup(ip) * aa * (1.d0+aux1)        &
          - 2.0d0 * aux4 * tlp(id) * taup(ip) * taup(ip)          &
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

!---> Integral on the vu fluid velocity

  p11 = sqrt( gama2*aux6 )
  ter3f = p11*vagaus(ip,i0)

!---> Integral on the particles velocity

  aux9  = 0.5d0 * tlp(id) * (1.d0 - aux2*aux2)
  aux10 = 0.5d0 * taup(ip) * (1.d0 - aux1*aux1)
  aux11 = taup(ip) * tlp(id) * (1.d0 - aux1*aux2)                 &
         / (taup(ip) + tlp(id))

  grga2 = (aux9 - 2.d0*aux11 + aux10) * aux8
  gagam = (aux9 - aux11) * (aux8 / aux3)
  gaome = ( (tlp(id) - taup(ip)) * (aux5 - aa)                    &
          -  tlp(id) * aux9 - taup(ip) * aux10                    &
          + (tlp(id) + taup(ip)) * aux11 ) * aux8

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

!---> vu fluid velocity

  vvue(id) = ter1f + ter2f + ter3f

!---> particles velocity

  vpart(id) = ter1p + ter2p + ter3p + ter4p + ter5p

enddo


!===============================================================================
! 3. Reference frame change:
!---------------------------
!local reference frame for the boundary face  --> in the global reference frame
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

! 3.3 - vu fluid velocity

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
