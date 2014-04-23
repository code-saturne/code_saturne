!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

 ( nbpmax , nvp    , nvep   , nivep  ,                            &
   itepa  ,                                                       &
   ettp   , ettpa  , tepa   , rtp , parbor, nvisbr)

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
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
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
use dimens, only: nvar
use cstphy
use cstnum
use lagpar
use lagran
use ppthch
use entsor
use mesh

!===============================================================================

implicit none

! Arguments

integer          nbpmax , nvp    , nvep  , nivep , nvisbr
integer          itepa(nbpmax,nivep)

double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision rtp (ncelet,nflown:nvar)
double precision parbor(nfabor,nvisbr)

! Local variables

integer ip, ii, ndiam, test_colli
double precision kinetic_energy
double precision  adhesion_energ
double precision norm_velocity, norm_face


double precision omep, domep
double precision v_part_t, v_part_t_dt, v_part_inst
double precision sub_dt

! ==========================================================================
! 0.    initialization
! ==========================================================================


! ==========================================================================
! 1.    Resuspension sub model
! ==========================================================================

do ip = 1, nbpart

   test_colli = 0

   if (itepa(ip,jdepo).eq.1) then

      ! The particle has just deposited
      ! The adhesion force is calculated

      call lagadh                                                     &
           ( ip   ,                                                   &
           nbpmax , nvp    , nvep   , nivep  ,                        &
           itepa  ,                                                   &
           ettp   , tepa   , rtp , adhesion_energ)

   elseif (itepa(ip,jdepo).eq.2) then

      ! The particle is rolling
      ! if the number of great asperities
      ! is null it is marked for a possible collision

      if (itepa(ip,jnbasg).eq.0) then
         test_colli = 1
      endif

      if (tepa(ip,jndisp).gt.ettp(ip,jdp).and. &
           tepa(ip,jndisp).lt. 2.d0 * ettp(ip,jdp)) then

         ! If the particle has a displacement approximately
         ! equal to a diameter, recalculation of the adhesion force

         tepa(ip,jndisp) = 0.d0

         call lagadh                                                     &
              ( ip   ,                                                   &
              nbpmax , nvp    , nvep   , nivep  ,                        &
              itepa  ,                                                   &
              ettp   , tepa   , rtp , adhesion_energ)

            if ((test_colli.eq.1) .and. (itepa(ip,jnbasg).gt.0)) then

               kinetic_energy = 0.5d0 * ettp(ip,jmp) * (ettp(ip,jup)**2      &
                                                  +    ettp(ip,jvp)**2       &
                                                  +    ettp(ip,jwp)**2)

               if (kinetic_energy.gt.adhesion_energ) then

                  ! The particle is resuspended
                  ! and its kinetic energy is totally converted
                  ! along the wall-normal distance

                  itepa(ip,jdepo) = 0

                  tepa(ip,jfadh) = 0.d0
                  tepa(ip,jmfadh) = 0.d0

                  itepa(ip,jnbasg) = 0
                  itepa(ip,jnbasp) = 0

                  tepa(ip,jndisp) = 0.d0

                  norm_face = surfbn(itepa(ip,jdfac))

                  norm_velocity = sqrt(ettp(ip,jup)**2 + ettp(ip,jvp)**2 + ettp(ip,jwp)**2)

                  ettp(ip,jup) = - norm_velocity / norm_face * surfbo(1, itepa(ip,jdfac))
                  ettp(ip,jvp) = - norm_velocity / norm_face * surfbo(2, itepa(ip,jdfac))
                  ettp(ip,jwp) = - norm_velocity / norm_face * surfbo(3, itepa(ip,jdfac))

                  ! Update of the number and weight of resuspended particles

                  nbpres = nbpres + 1
                  dnbres = dnbres + tepa(ip, jrpoi)

                  parbor(itepa(ip,jdfac),ires) = parbor(itepa(ip,jdfac),ires) + tepa(ip,jrpoi)

                  parbor(itepa(ip,jdfac),iflres) = parbor(itepa(ip,jdfac),iflres)                  &
                              + ( tepa(ip,jrpoi) * ettp(ip,jmp) / norm_face)

                  parbor(itepa(ip,jdfac),iflm) = parbor(itepa(ip,jdfac),iflm)                   &
                              - ( tepa(ip,jrpoi) * ettp(ip,jmp) / norm_face)


               endif

            endif

      elseif (tepa(ip,jndisp).ge. 2d0 * ettp(ip,jdp)) then

         ndiam = floor(tepa(ip,jndisp) / ettp(ip,jdp))

         ii = 1

         do while ((ii.le.ndiam).and.(itepa(ip,jdepo).eq.2))

            call lagadh                                                     &
                 ( ip   ,                                                   &
                 nbpmax , nvp    , nvep   , nivep  ,                        &
                 itepa  ,                                                   &
                 ettp   , tepa   , rtp , adhesion_energ)

            ! Reconstruct an estimate of the particle velocity
            ! at the current sub-time-step assuming linear variation
            ! (constant acceleration)

            v_part_t = sqrt( ettpa(ip,jup) ** 2                         &
                           + ettpa(ip,jvp) ** 2                         &
                           + ettpa(ip,jwp) ** 2)

            v_part_t_dt = sqrt( ettp(ip,jup) ** 2                       &
                              + ettp(ip,jvp) ** 2                       &
                              + ettp(ip,jwp) ** 2)

            sub_dt = dtp / ndiam

            v_part_inst =  v_part_t + sub_dt * (v_part_t_dt + v_part_t) / dtp

            ! Reconstruct an estimate of the angular velocity
            ! at the current sub-time-step

            omep = v_part_inst / (ettp(ip,jdp) * 0.5d0)

            ! Variation of the angular velocity due to
            ! the update of the adhesion torque

            domep = tepa(ip,jmfadh)                                &
                 /((7.d0/5.d0)*ettp(ip,jmp)*(ettp(ip,jdp) * 0.5d0)**2)

            if ((domep * sub_dt) .gt. omep) then

               itepa(ip,jdepo) = 10

               ettp(ip,jup) = 0.d0
               ettp(ip,jvp) = 0.d0
               ettp(ip,jwp) = 0.d0

            endif


            if ((test_colli.eq.1) .and. (itepa(ip,jnbasg).gt.0)) then

               kinetic_energy = 0.5d0 * ettp(ip,jmp) * (ettp(ip,jup)**2      &
                                                  +    ettp(ip,jvp)**2       &
                                                  +    ettp(ip,jwp)**2)


              if (kinetic_energy.gt.adhesion_energ) then

                  ! The particle is resuspended
                  ! and its kinetic energy is totally converted
                  ! along the wall-normal distance

                  itepa(ip,jdepo) = 0

                  tepa(ip,jfadh) = 0.d0
                  tepa(ip,jmfadh) = 0.d0

                  itepa(ip,jnbasg) = 0
                  itepa(ip,jnbasp) = 0

                  tepa(ip,jndisp) = 0.d0

                  norm_face = surfbn(itepa(ip,jdfac))

                  norm_velocity = sqrt(ettp(ip,jup)**2 + ettp(ip,jvp)**2 + ettp(ip,jwp)**2)

                  ettp(ip,jup) = - norm_velocity / norm_face * surfbo(1, itepa(ip,jdfac))
                  ettp(ip,jvp) = - norm_velocity / norm_face * surfbo(2, itepa(ip,jdfac))
                  ettp(ip,jwp) = - norm_velocity / norm_face * surfbo(3, itepa(ip,jdfac))

                  ! Update of the number and weight of resuspended particles

                  nbpres = nbpres + 1
                  dnbres = dnbres + tepa(ip, jrpoi)

                  parbor(itepa(ip,jdfac),ires) = parbor(itepa(ip,jdfac),ires) + tepa(ip,jrpoi)

                  parbor(itepa(ip,jdfac),iflres) = parbor(itepa(ip,jdfac),iflres)                  &
                                     + ( tepa(ip,jrpoi) * ettp(ip,jmp) / norm_face)

                  parbor(itepa(ip,jdfac),iflm) = parbor(itepa(ip,jdfac),iflm)                   &
                                   - ( tepa(ip,jrpoi) * ettp(ip,jmp) / norm_face)
               endif

               if (itepa(ip,jnbasg).eq.0) then
                  test_colli = 1

               endif

            endif ! if test_colli

            ii = ii + 1

         enddo ! do while ..

      endif ! if tepa(ip,jndisp)


   endif  ! if jdepo = ...

enddo

end subroutine lagres
