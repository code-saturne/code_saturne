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

subroutine lagres &
!================

 ( parbor, nvisbr)

!===============================================================================

! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module:
!   ------------------------------------------------------


!   Calculation of the particle resuspension
!
!
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! parbor           ! tr ! <-->! infos sur interaction des particules           !
!(nfabor,nvisbr)   !    !     !   aux faces de bord                            !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================


!===============================================================================


!===============================================================================
! Module files
!===============================================================================

use paramx
use cstphy
use cstnum
use lagpar
use lagran
use numvar
use optcal
use ppthch
use entsor
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvisbr

double precision parbor(nfabor,nvisbr)

! Local variables

integer ip, ii, iel, ndiam, test_colli
double precision kinetic_energy
double precision adhesion_energ
double precision norm_velocity, norm_face


double precision omep, domep
double precision v_part_t, v_part_t_dt, v_part_inst
double precision sub_dt

double precision, dimension(:), pointer :: cvar_scalt

! ==========================================================================
! 0.    initialization
! ==========================================================================

if (iscalt.gt.0) call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

! ==========================================================================
! 1.    Resuspension sub model
! ==========================================================================

do ip = 1, nbpart

   test_colli = 0

   iel = ipepa(jisor,ip)

   if (ipepa(jdepo,ip).eq.1) then

      ! The particle has just deposited
      ! The adhesion force is calculated

      call lagadh(ip, cvar_scalt(iel), adhesion_energ)

   elseif (ipepa(jdepo,ip).eq.2) then

      ! The particle is rolling
      ! if the number of great asperities
      ! is null it is marked for a possible collision

      if (ipepa(jnbasg,ip).eq.0) then
         test_colli = 1
      endif

      if (pepa(jndisp,ip).gt.eptp(jdp,ip).and. &
           pepa(jndisp,ip).lt. 2.d0 * eptp(jdp,ip)) then

         ! If the particle has a displacement approximately
         ! equal to a diameter, recalculation of the adhesion force

         pepa(jndisp,ip) = 0.d0

         call lagadh(ip, cvar_scalt(iel), adhesion_energ)

            if ((test_colli.eq.1) .and. (ipepa(jnbasg,ip).gt.0)) then

               kinetic_energy = 0.5d0 * eptp(jmp,ip) * (eptp(jup,ip)**2      &
                                                  +    eptp(jvp,ip)**2       &
                                                  +    eptp(jwp,ip)**2)

               if (kinetic_energy.gt.adhesion_energ) then

                  ! The particle is resuspended
                  ! and its kinetic energy is totally converted
                  ! along the wall-normal distance

                  ipepa(jdepo,ip) = 0

                  pepa(jfadh,ip) = 0.d0
                  pepa(jmfadh,ip) = 0.d0

                  ipepa(jnbasg,ip) = 0
                  ipepa(jnbasp,ip) = 0

                  pepa(jndisp,ip) = 0.d0

                  norm_face = surfbn(ipepa(jdfac,ip))

                  norm_velocity = sqrt(eptp(jup,ip)**2 + eptp(jvp,ip)**2 + eptp(jwp,ip)**2)

                  eptp(jup,ip) = - norm_velocity / norm_face * surfbo(1, ipepa(jdfac,ip))
                  eptp(jvp,ip) = - norm_velocity / norm_face * surfbo(2, ipepa(jdfac,ip))
                  eptp(jwp,ip) = - norm_velocity / norm_face * surfbo(3, ipepa(jdfac,ip))

                  ! Update of the number and weight of resuspended particles

                  nbpres = nbpres + 1
                  dnbres = dnbres + pepa( jrpoi,ip)

                  parbor(ipepa(jdfac,ip),ires) = parbor(ipepa(jdfac,ip),ires) + pepa(jrpoi,ip)

                  parbor(ipepa(jdfac,ip),iflres) = parbor(ipepa(jdfac,ip),iflres)                  &
                              + ( pepa(jrpoi,ip) * eptp(jmp,ip) / norm_face)

                  parbor(ipepa(jdfac,ip),iflm) = parbor(ipepa(jdfac,ip),iflm)                   &
                              - ( pepa(jrpoi,ip) * eptp(jmp,ip) / norm_face)


               endif

            endif

      elseif (pepa(jndisp,ip).ge. 2d0 * eptp(jdp,ip)) then

         ndiam = floor(pepa(jndisp,ip) / eptp(jdp,ip))

         ii = 1

         do while ((ii.le.ndiam).and.(ipepa(jdepo,ip).eq.2))

            call lagadh(ip, cvar_scalt(iel), adhesion_energ)

            ! Reconstruct an estimate of the particle velocity
            ! at the current sub-time-step assuming linear variation
            ! (constant acceleration)

            v_part_t = sqrt( eptpa(jup,ip) ** 2                         &
                           + eptpa(jvp,ip) ** 2                         &
                           + eptpa(jwp,ip) ** 2)

            v_part_t_dt = sqrt( eptp(jup,ip) ** 2                       &
                              + eptp(jvp,ip) ** 2                       &
                              + eptp(jwp,ip) ** 2)

            sub_dt = dtp / ndiam

            v_part_inst =  v_part_t + sub_dt * (v_part_t_dt + v_part_t) / dtp

            ! Reconstruct an estimate of the angular velocity
            ! at the current sub-time-step

            omep = v_part_inst / (eptp(jdp,ip) * 0.5d0)

            ! Variation of the angular velocity due to
            ! the update of the adhesion torque

            domep = pepa(jmfadh,ip)                                &
                 /((7.d0/5.d0)*eptp(jmp,ip)*(eptp(jdp,ip) * 0.5d0)**2)

            if ((domep * sub_dt) .gt. omep) then

               ipepa(jdepo,ip) = 10

               eptp(jup,ip) = 0.d0
               eptp(jvp,ip) = 0.d0
               eptp(jwp,ip) = 0.d0

            endif


            if ((test_colli.eq.1) .and. (ipepa(jnbasg,ip).gt.0)) then

               kinetic_energy = 0.5d0 * eptp(jmp,ip) * (eptp(jup,ip)**2      &
                                                  +    eptp(jvp,ip)**2       &
                                                  +    eptp(jwp,ip)**2)


              if (kinetic_energy.gt.adhesion_energ) then

                  ! The particle is resuspended
                  ! and its kinetic energy is totally converted
                  ! along the wall-normal distance

                  ipepa(jdepo,ip) = 0

                  pepa(jfadh,ip) = 0.d0
                  pepa(jmfadh,ip) = 0.d0

                  ipepa(jnbasg,ip) = 0
                  ipepa(jnbasp,ip) = 0

                  pepa(jndisp,ip) = 0.d0

                  norm_face = surfbn(ipepa(jdfac,ip))

                  norm_velocity = sqrt(eptp(jup,ip)**2 + eptp(jvp,ip)**2 + eptp(jwp,ip)**2)

                  eptp(jup,ip) = - norm_velocity / norm_face * surfbo(1, ipepa(jdfac,ip))
                  eptp(jvp,ip) = - norm_velocity / norm_face * surfbo(2, ipepa(jdfac,ip))
                  eptp(jwp,ip) = - norm_velocity / norm_face * surfbo(3, ipepa(jdfac,ip))

                  ! Update of the number and weight of resuspended particles

                  nbpres = nbpres + 1
                  dnbres = dnbres + pepa( jrpoi,ip)

                  parbor(ipepa(jdfac,ip),ires) = parbor(ipepa(jdfac,ip),ires) + pepa(jrpoi,ip)

                  parbor(ipepa(jdfac,ip),iflres) = parbor(ipepa(jdfac,ip),iflres)                  &
                                     + ( pepa(jrpoi,ip) * eptp(jmp,ip) / norm_face)

                  parbor(ipepa(jdfac,ip),iflm) = parbor(ipepa(jdfac,ip),iflm)                   &
                                   - ( pepa(jrpoi,ip) * eptp(jmp,ip) / norm_face)
               endif

               if (ipepa(jnbasg,ip).eq.0) then
                  test_colli = 1

               endif

            endif ! if test_colli

            ii = ii + 1

         enddo ! do while ..

      endif ! if pepa(jndisp,ip)


   endif  ! if jdepo = ...

enddo

end subroutine lagres
