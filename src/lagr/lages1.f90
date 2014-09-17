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

subroutine lages1 &
!================

 ( nbpmax ,                                                       &
   rtpa   , propce ,                                              &
   taup   , tlag   , piil   ,                                     &
   bx     , vagaus , gradpr , romp   ,                            &
   brgaus , terbru , fextla )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    INTEGRATION DES EDS PAR UN SCHEMA D'ORDRE 1 (scheO1)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (pas de temps precedent)           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! piil(nbpmax,3    ! tr ! <-- ! terme dans l'integration des eds up            !
! bx(nbpmax,3,2    ! tr ! <-- ! caracteristiques de la turbulence              !
! vagaus           ! tr ! <-- ! variables aleatoires gaussiennes               !
!(nbpmax,nvgaus    !    !     !                                                !
! gradpr(3,ncel)   ! tr ! <-- ! gradient de pression                           !
! romp             ! tr ! <-- ! masse volumique des particules                 !
! fextla           ! tr ! <-- ! champ de forces exterieur                      !
!(ncelet,3)        !    !     !    utilisateur (m/s2)                          !
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
use dimens, only: nvar
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nbpmax

double precision rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3) , bx(nbpmax,3,2)
double precision vagaus(nbpmax,*)
double precision brgaus(nbpmax,*) , terbru(nbpmax)
double precision gradpr(3,ncelet)
double precision romp(nbpmax)
double precision fextla(nbpmax,3)

! Local variables

integer          iel , ip , id , i0 , mode

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
double precision grav(3) , rom , vitf
double precision tempf
double precision ddbr, tix2, tiu2, tixiu
double precision tbrix1, tbrix2, tbriu

double precision, dimension(:), pointer :: cromf
double precision, dimension(:,:), pointer :: vela

!===============================================================================

! Map field arrays
call field_get_val_prev_v(ivarfl(iu), vela)

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

vitf = 0.d0

grav(1) = gx
grav(2) = gy
grav(3) = gz

! Pointeur sur la masse volumique en fonction de l'ecoulement

if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then
  call field_get_val_s(iprpfl(ipproc(irom1)), cromf)
else
  call field_get_val_s(icrom, cromf)
endif

!===============================================================================
! 2. INTEGRATION DES EDS SUR LES PARTICULES
!===============================================================================

do id = 1,3

  i0 = id - 1

  do ip = 1,nbpart

    if (ipepa(jisor,ip).gt.0) then

      iel = ipepa(jisor,ip)

      rom = cromf(iel)

      if (id.eq.1) vitf = vela(1,iel)
      if (id.eq.2) vitf = vela(2,iel)
      if (id.eq.3) vitf = vela(3,iel)

!---> (2.1) Calcul preliminaires :
!     ----------------------------

!      calcul de II*TL+<u> et [(grad<P>/rhop+g)*tau_p+<Uf>] ?

      tci = piil(ip,id) * tlag(ip,id) + vitf

      force =                                                     &
        ( rom * gradpr(id,iel) / romp(ip)                         &
         +grav(id)+ fextla(ip,id)        ) * taup(ip)

!---> (2.2) Calcul des coefficients/termes deterministes :
!     ----------------------------------------------------

!  sur HP-UX : l'option de compil +T genere des messages du type :
!  === MATH LIBRARY ERROR 14: DEXP(X) UNDERFLOW
!  au passage de l'exponentiel si l'exposant est < a -7.0D2


      aux1 = exp( -dtp / taup(ip))
      aux2 = exp( -dtp / tlag(ip,id))
      aux3 = tlag(ip,id) / (tlag(ip,id)-taup(ip))

      aux4 = tlag(ip,id) / (tlag(ip,id)+taup(ip))
      aux5 = tlag(ip,id) * (1.d0-aux2)
      aux6 = bx(ip,id,nor) * bx(ip,id,nor) * tlag(ip,id)

      aux7 = tlag(ip,id) - taup(ip)
      aux8 = bx(ip,id,nor) * bx(ip,id,nor) * aux3**2

!---> termes pour la trajectoire

      aa = taup(ip) * (1.d0 - aux1)
      bb = (aux5 - aa) * aux3
      cc = dtp - aa - bb

      ter1x = aa * eptpa(jup+i0,ip)
      ter2x = bb * eptpa(juf+i0,ip)
      ter3x = cc * tci
      ter4x = (dtp - aa) * force

!---> termes pour le fluide vu

      ter1f = eptpa(juf+i0,ip) * aux2
      ter2f = tci * (1.d0-aux2)

!---> termes pour la vitesse des particules

      dd = aux3 * (aux2 - aux1)
      ee = 1.d0 - aux1

      ter1p = eptpa(jup+i0,ip) * aux1
      ter2p = eptpa(juf+i0,ip) * dd
      ter3p = tci * (ee-dd)
      ter4p = force * ee

!---> (2.3) Calcul des coefficients pour les integrales stochastiques :


!---> integrale sur la position des particules

      gama2  = 0.5d0 * (1.d0 - aux2*aux2 )
      omegam = 0.5d0 * aux4 * ( aux5 - aux2*aa )                  &
              -0.5d0 * aux2 * bb
      omegam = omegam * sqrt(aux6)

      omega2 = aux7                                               &
             * (aux7*dtp - 2.d0 * (tlag(ip,id)*aux5-taup(ip)*aa)) &
           + 0.5d0 * tlag(ip,id) * tlag(ip,id) * aux5             &
             * (1.d0 + aux2)                                      &
           + 0.5d0 * taup(ip) * taup(ip) * aa * (1.d0+aux1)       &
           - 2.d0 * aux4 * tlag(ip,id) * taup(ip) * taup(ip)      &
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

!---> integrale sur la vitesse du fluide vu

      p11 = sqrt( gama2*aux6 )
      ter3f = p11*vagaus(ip,id)

!---> integrale sur la vitesse des particules

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

      ter5p = p31 * vagaus(ip,id)                                 &
            + p32 * vagaus(ip,3+id)                               &
            + p33 * vagaus(ip,6+id)


!---> (2.3) Calcul des Termes dans le cas du mouvement Brownien :


      if ( lamvbr .eq. 1 ) then

        !  Calcul de la temperature du fluide en fonction du type
        !  d'ecoulement

        if ( ippmod(iccoal).ge.0 .or.                             &
             ippmod(icpl3c).ge.0      ) then

          tempf = propce(iel,ipproc(itemp1))

        else if ( ippmod(icod3p).ge.0 .or.                        &
          ippmod(icoebu).ge.0 .or.                                &
          ippmod(ielarc).ge.0 .or.                                &
          ippmod(ieljou).ge.0      ) then

          tempf = propce(iel,ipproc(itemp))

        else if (itherm.eq.1 .and. itpscl.eq.2) then
          tempf = rtpa(iel,isca(iscalt))

        else if (itherm.eq.1 .and. itpscl.eq.1) then
          tempf = rtpa(iel,isca(iscalt))

        else if (itherm.eq.2) then
          mode = 1
          call usthht(mode,rtpa(iel,isca(iscalt)),tempf)
          !==========
          tempf = tempf+tkelvi
        else
          tempf = t0
        endif

        ddbr  = sqrt( 2.d0*kboltz*tempf                           &
                     /(eptp(jmp,ip)*taup(ip)) )
        tix2  = (taup(ip)*ddbr)**2                                &
                *(dtp - taup(ip)*(1.d0-aux1)                      &
                               *(3.d0-aux1) /2.d0 )

        tiu2  = ddbr*ddbr*taup(ip)                                &
               *(1.d0-exp(-2.d0*dtp/taup(ip)))/2.d0

        tixiu = (ddbr*taup(ip)*(1.d0-aux1))**2/2.d0

        tbrix2 = tix2-(tixiu*tixiu)/tiu2
        if ( tbrix2 .gt.0.d0 ) then
          tbrix2 = sqrt(tbrix2)*brgaus(ip,id)
        else
          tbrix2 = 0.d0
        endif

        if ( tiu2 .gt. 0.d0 ) then
          tbrix1 = tixiu/sqrt(tiu2)*brgaus(ip,3+id)
        else
          tbrix1 = 0.d0
        endif

        if ( tiu2 .gt. 0.d0 ) then
          tbriu = sqrt(tiu2)*brgaus(ip,3+id)
          terbru(ip) = sqrt(tiu2)
        else
          tbriu = 0.d0
          terbru(ip) = 0.d0
        endif

      else
        tbrix1 = 0.d0
        tbrix2 = 0.d0
        tbriu  = 0.d0
      endif

!===============================================================================
! 3. FINALISATION DES ECRITURES
!===============================================================================

!---> trajectoire

      eptp(jxp+i0,ip) = eptpa(jxp+i0,ip)                          &
                       + ter1x + ter2x + ter3x + ter4x + ter5x    &
                       + tbrix1 + tbrix2

!---> vitesse fluide vu

      eptp(juf+i0,ip) = ter1f + ter2f + ter3f

!---> vitesse particules

      eptp(jup+i0,ip) = ter1p + ter2p + ter3p + ter4p + ter5p     &
                       + tbriu

    endif

  enddo

enddo

!----
! FIN
!----

return
end subroutine lages1
