!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine lagdep &
!================

 ( taup   , tlag   , piil   ,                                     &
   bx     , vagaus , gradpr , romp   ,                            &
   fextla , vislen)

!===============================================================================
! Purpose:
! ----------
!
!   Subroutine of the Lagrangian particle-tracking module:
!   ------------------------------------------------------
!
!
!   Deposition sub-model:
!
!   Main subroutine of the submodel
!
!   1/ Calculation of the normalized wall-normal distance of the boundary-cell particles
!   2/ Sorting of the particles with respect to their normalized wall-normal distance
!         * If y^+ > depint : the standard Langevin model is applied
!         * If y^+ < depint : the deposition submodel is applied
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! taup(nbpart)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpart)     ! tr ! <-- ! temps caracteristique fluide                   !
! piil(nbpart,3)   ! tr ! <-- ! terme dans l'integration des eds up            !
! bx(nbpart,3,2)   ! tr ! <-- ! caracteristiques de la turbulence              !
! vagaus           ! tr ! <-- ! variables aleatoires gaussiennes               !
!(nbpart,nvgaus)   !    !     !                                                !
! gradpr(3,ncel)   ! tr ! <-- ! gradient de pression                           !
! romp             ! tr ! <-- ! masse volumique des particules                 !
! fextla           ! tr ! <-- ! champ de forces exterieur                      !
!(ncelet,3)        !    !     !    utilisateur (m/s2)                          !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
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
use mesh
use field
use pointe

!===============================================================================

implicit none

! Arguments

double precision taup(nbpart) , tlag(nbpart,3)
double precision piil(nbpart,3) , bx(nbpart,3,2)
double precision vagaus(nbpart,*)
double precision gradpr(3,ncelet)
double precision romp(nbpart)
double precision fextla(nbpart,3)

! Local variables

integer          iel , ifac, ip , id , i0 , mode

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
double precision grav(3) , romf , vitf
double precision tempf, lvisq, tvisq
double precision ustar, visccf,depint

double precision vislen(nfabor)

double precision, dimension(:), pointer :: cromf
double precision, dimension(:,:), pointer :: vela
double precision, dimension(:), pointer :: viscl

double precision, dimension(:), pointer :: cvara_scalt
double precision, dimension(:), pointer :: cpro_temp1, cpro_temp

!===============================================================================

! Map field arrays
call field_get_val_prev_v(ivarfl(iu), vela)

!===============================================================================
! 1.  Initialization
!===============================================================================

!  Initialize variables to avoid compiler warnings

vitf = 0.d0

grav(1) = gx
grav(2) = gy
grav(3) = gz

! Interface location between near-wall region
! and core of the flow (normalized units)

depint = 100.d0

!  The density pointer according to the flow location

if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then
  call field_get_val_s(iprpfl(ipproc(irom1)), cromf)
else
  call field_get_val_s(icrom, cromf)
endif

if (iscalt.gt.0) call field_get_val_prev_s(ivarfl(isca(iscalt)), cvara_scalt)

if (itemp1.gt.0) call field_get_val_s(iprpfl(itemp1),cpro_temp1)
if (itemp.gt.0) call field_get_val_s(iprpfl(itemp),cpro_temp)

call field_get_val_s(iprpfl(iviscl), viscl)

!===============================================================================
! 2. loop on the particles
!===============================================================================

  do ip = 1,nbpart

    if (ipepa(jisor,ip).gt.0) then

      iel = ipepa(jisor,ip)
      romf = cromf(iel)
      visccf = viscl(iel) / romf

    ! Fluid temperature computation depending on the type of flow

      if (ippmod(iccoal).ge.0 .or.                             &
          ippmod(icpl3c).ge.0 .or.                             &
          ippmod(icfuel).ge.0) then

         tempf = cpro_temp1(iel)

      else if (ippmod(icod3p).ge.0 .or.                        &
               ippmod(icoebu).ge.0 .or.                        &
               ippmod(ielarc).ge.0 .or.                        &
           ippmod(ieljou).ge.0      ) then

         tempf = cpro_temp(iel)

      else if (iscalt.gt.0) then
        if (itherm.eq.1) then
          if (itpscl.eq.2) then
            tempf = cvara_scalt(iel) + tkelvi
          else if (itpscl.eq.1 ) then
            tempf = cvara_scalt(iel)
          endif
        else if (itherm.eq.2) then
          mode = 1
          call usthht(mode,cvara_scalt(iel),tempf)
          !==========
          tempf = tempf+tkelvi
        endif
      else
        tempf = t0
      endif


!=========================================================================
!   If y^+ is greater than the interface location,
!   the standard model is applied
!=========================================================================

     if (pepa(jryplu,ip).gt.depint) then

         ipepa(jimark,ip) = -1

         do id = 1,3

            i0 = id - 1
            if (id.eq.1) vitf = vela(1,iel)
            if (id.eq.2) vitf = vela(2,iel)
            if (id.eq.3) vitf = vela(3,iel)

            tci = piil(ip,id) * tlag(ip,id) + vitf

            if (iadded_mass.eq.0) then
              force =                                                     &
                ( romf * gradpr(id,iel) / romp(ip)                        &
                + grav(id)+ fextla(ip,id)        ) * taup(ip)

              ! Added-mass term?
            else
              force =                                                     &
                ( romf * gradpr(id,iel) / romp(ip)                        &
                  * (1.d0 + 0.5d0*added_mass_const)                       &
                  / (1.d0 + 0.5d0*added_mass_const*romf/romp(ip))         &
                + grav(id)+ fextla(ip,id)        ) * taup(ip)

            endif

            aux1 = exp( -dtp / taup(ip))
            aux2 = exp( -dtp / tlag(ip,id))
            aux3 = tlag(ip,id) / (tlag(ip,id)-taup(ip))

            aux4 = tlag(ip,id) / (tlag(ip,id)+taup(ip))
            aux5 = tlag(ip,id) * (1.d0-aux2)
            aux6 = bx(ip,id,nor) * bx(ip,id,nor) * tlag(ip,id)

            aux7 = tlag(ip,id) - taup(ip)
            aux8 = bx(ip,id,nor) * bx(ip,id,nor) * aux3**2

            !---> trajectory terms

            aa = taup(ip) * (1.d0 - aux1)
            bb = (aux5 - aa) * aux3
            cc = dtp - aa - bb

            ter1x = aa * eptpa(jup+i0,ip)
            ter2x = bb * eptpa(juf+i0,ip)
            ter3x = cc * tci
            ter4x = (dtp - aa) * force

            !---> flow-seen velocity terms

            ter1f = eptpa(juf+i0,ip) * aux2
            ter2f = tci * (1.d0-aux2)

            !---> termes pour la vitesse des particules

            dd = aux3 * (aux2 - aux1)
            ee = 1.d0 - aux1

            ter1p = eptpa(jup+i0,ip) * aux1
            ter2p = eptpa(juf+i0,ip) * dd
            ter3p = tci * (ee-dd)
            ter4p = force * ee

            !---> (2.3) Coefficients computation for the stochastic integral

            !---> Integral for particles position

            gama2  = 0.5d0 * (1.d0 - aux2*aux2 )
            omegam = 0.5d0 * aux4 * ( aux5 - aux2*aa )                  &
                 -0.5d0 * aux2 * bb
            omegam = omegam * sqrt(aux6)

            omega2 = aux7                                               &
                   * (aux7*dtp - 2.d0 * (tlag(ip,id)*aux5-taup(ip)*aa)) &
                   + 0.5d0 * tlag(ip,id) * tlag(ip,id) * aux5           &
                    * (1.d0 + aux2)                                     &
                   + 0.5d0 * taup(ip) * taup(ip) * aa * (1.d0+aux1)     &
                   - 2.d0 * aux4 * tlag(ip,id) * taup(ip) * taup(ip)    &
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

            ter5x = p21 * vagaus(ip,id) + p22 * vagaus(ip,3+id)

            !---> integral for the flow-seen velocity

            p11 = sqrt( gama2*aux6 )
            ter3f = p11*vagaus(ip,id)

            !---> integral for the particles velocity

            aux9  = 0.5d0 * tlag(ip,id) * (1.d0 - aux2*aux2)
            aux10 = 0.5d0 * taup(ip) * (1.d0 - aux1*aux1)
            aux11 = taup(ip) * tlag(ip,id) * (1.d0 - aux1*aux2)         &
                 / (taup(ip) + tlag(ip,id))

            grga2 = (aux9 - 2.d0*aux11 + aux10) * aux8
            gagam = (aux9 - aux11) * (aux8 / aux3)
            gaome = ( (tlag(ip,id) - taup(ip)) * (aux5 - aa)            &
                 - tlag(ip,id) * aux9 - taup(ip) * aux10              &
                 + (tlag(ip,id) + taup(ip)) * aux11) * aux8

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

            ter5p = p31 * vagaus(ip,id)                                &
                 + p32 * vagaus(ip,3+id)                               &
                 + p33 * vagaus(ip,6+id)

            !
            ! Update of the particle state-vector
            !

            eptp(jxp+i0,ip) = eptpa(jxp+i0,ip)                         &
                 + ter1x + ter2x + ter3x + ter4x + ter5x

            eptp(juf+i0,ip) = ter1f + ter2f + ter3f

            eptp(jup+i0,ip) = ter1p + ter2p + ter3p + ter4p + ter5p

         enddo

!===============================================================================
!   Otherwise, the deposition submodel is applied :
!!===============================================================================

         else

            if (pepa(jryplu,ip).lt.pepa(jrinpf,ip)) then

               if ( ipepa(jimark,ip) .lt. 0 ) then
                  ipepa(jimark,ip) = 10
               else
                  ipepa(jimark,ip) = 0
               endif

            else

               !   if (jrinpf < y^+ < depint)

               if ( ipepa(jimark,ip) .lt. 0 ) then
                  ipepa(jimark,ip) = 20
               else if ( ipepa(jimark,ip) .eq. 0 ) then
                  ipepa(jimark,ip) = 30
               endif

            endif

            ifac = ipepa(jdfac,ip)

            ustar = uetbor(ifac)
            lvisq = vislen(ifac)
            if (ustar.gt.0.d0) then
              tvisq = lvisq / ustar
            else
              tvisq = grand
            endif

            call lagesd                                                    &
            ( ifac   , ip     ,                                            &
              taup   , piil   ,                                            &
              vagaus , gradpr , romp   ,                                   &
              tempf  , romf   , ustar  , lvisq  ,tvisq   , depint )

         endif

      endif

   enddo

!----
! End
!----

end subroutine
