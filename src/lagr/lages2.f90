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

subroutine lages2 &
!================

 ( propce ,                                                       &
   taup   , tlag   , piil   ,                                     &
   bx     , tsfext ,                                              &
   vagaus , gradpr ,                                              &
   romp   , brgaus , terbru , fextla )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    INTEGRATION DES EDS PAR UN SCHEMA D'ORDRE 2 (sche2c)

!     Lorsqu'il y a eu interaction avec une face de bord,
!     les calculs de la vitesse de la particule et de
!     la vitesse du fluide vu sont forcement a l'ordre 1
!     (meme si on est a l'ordre 2 par ailleurs).

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! statis           ! tr ! <-- !  cumul des statistiques volumiques             !
!(ncelet,nvlsta    !    !     !                                                !
! taup(nbpart)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpart)     ! tr ! <-- ! temps caracteristique fluide                   !
! piil(nbpart,3)   ! tr ! <-- ! terme dans l'integration des eds up            !
! bx(nbpart,3,2)   ! tr ! <-- ! caracteristiques de la turbulence              !
! tsfext(nbpart)   ! tr ! --> ! infos pour couplage retour dynamique           !
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

!===============================================================================

implicit none

! Arguments

double precision propce(ncelet,*)
double precision taup(nbpart) , tlag(nbpart,3)
double precision piil(nbpart,3) , bx(nbpart,3,2)
double precision tsfext(nbpart)
double precision vagaus(nbpart,*), brgaus(nbpart,*)
double precision gradpr(3,ncelet)
double precision romp(nbpart)
double precision terbru(nbpart)
double precision fextla(nbpart,3)

! Local variables

integer          iel , ip , id , i0
double precision aux0 , aux1 , aux2 , aux3 , aux4 , aux5
double precision aux6 , aux7 , aux8 , aux9 , aux10 , aux11
double precision aux12 , aux14 , aux15 , aux16
double precision aux17 , aux18 , aux19 , aux20
double precision ter1 , ter2 , ter3 , ter4 , ter5
double precision sige , tapn
double precision gamma2 , omegam , omega2
double precision grgam2 , gagam , gaome
double precision p11 , p21 , p22 , p31 , p32 , p33
double precision grav(3) , rom , vitf
double precision tbriu

double precision, dimension(:), pointer :: cromf
double precision, dimension(:,:), pointer :: vel, vela
double precision, dimension(:,:), allocatable :: auxl

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)
call field_get_val_prev_v(ivarfl(iu), vela)

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

allocate(auxl(nbpart,7))

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

!===============================================================================
! 2.1 CALCUL A CHAQUE SOUS PAS DE TEMPS
!===============================================================================

!---> Calcul de tau_p*A_p et de II*TL+<u> :
!     -------------------------------------

do id = 1,3

  do ip = 1,nbpart

    if (ipepa(jisor,ip).gt.0) then

      iel = ipepa(jisor,ip)

      rom = cromf(iel)

      auxl(ip,id) =                                               &
        ( rom *gradpr(id,iel) /romp(ip)                           &
         +grav(id)+fextla(ip,id)       ) * taup(ip)

      if (nor.eq.1) then
        if (id.eq.1) vitf = vela(1,iel)
        if (id.eq.2) vitf = vela(2,iel)
        if (id.eq.3) vitf = vela(3,iel)
      else
        if (id.eq.1) vitf = vel(1,iel)
        if (id.eq.2) vitf = vel(2,iel)
        if (id.eq.3) vitf = vel(3,iel)
      endif
      auxl(ip,id+3) = piil(ip,id) * tlag(ip,id) + vitf

    endif

  enddo

enddo

!===============================================================================
! 2.2 ETAPE DE PREDICTION :
!===============================================================================

if (nor.eq.1) then


!---> Sauvegarde de tau_p^n
!     ---------------------

  do ip = 1,nbpart

    if (ipepa(jisor,ip).gt.0) then

      pepa(jtaux,ip) = taup(ip)

    endif

  enddo

!---> Sauvegarde couplage
!     -------------------

  if (iilagr.eq.2) then

    do ip = 1,nbpart

      if (ipepa(jisor,ip).gt.0) then

        aux0 = -dtp/taup(ip)
        aux1 = exp(aux0)
        tsfext(ip) = taup(ip) * eptp(jmp,ip)                      &
                   * (-aux1 + (aux1-1.d0) / aux0)

      endif

    enddo

  endif


!---> Chargement des termes a t = t_n :
!     ---------------------------------

  do id = 1,3

    i0 = id - 1

    do ip = 1,nbpart

      if (ipepa(jisor,ip).gt.0) then

        aux0 = -dtp / taup(ip)
        aux1 = -dtp / tlag(ip,id)
        aux2 = exp(aux0)
        aux3 = exp(aux1)
        aux4 = tlag(ip,id) / (tlag(ip,id) - taup(ip))
        aux5 = aux3 - aux2

        pepa(jtsuf(id),ip) =   0.5d0 * eptpa(juf+i0,ip) * aux3               &
                             + auxl(ip,id+3) * ( -aux3 + (aux3-1.d0) /aux1)

        ter1 = 0.5d0 * eptpa(jup+i0,ip) * aux2
        ter2 = 0.5d0 * eptpa(juf+i0,ip) * aux4 * aux5
        ter3 = auxl(ip,id+3)*                                     &
              ( - aux2 + ((tlag(ip,id) + taup(ip)) / dtp)         &
             * (1.d0 - aux2)                                      &
             - (1.d0 + tlag(ip,id) / dtp) * aux4 * aux5)
        ter4 = auxl(ip,id) * (-aux2 + (aux2 - 1.d0) / aux0)

        pepa(jtsup(id),ip) = ter1 + ter2 + ter3 + ter4

      endif

    enddo

  enddo

!---> Schema d'Euler :
!     ----------------

  call lages1                                                     &
  !==========
 ( propce ,                                                       &
   taup   , tlag   , piil   ,                                     &
   bx     , vagaus , gradpr , romp   ,                            &
   brgaus , terbru , fextla )

else

!===============================================================================
! 2.2 ETAPE DE CORRECTION :
!===============================================================================

!---> Calcul de Us :
!     --------------

  do id = 1,3

    i0 = id - 1

    do ip = 1,nbpart

      if (ipepa(jisor,ip).gt.0 .and. ipepa(jord1,ip).eq.0) then

        aux0 = -dtp / taup(ip)
        aux1 = -dtp / tlag(ip,id)
        aux2 = exp(aux0)
        aux3 = exp(aux1)
        aux4 = tlag(ip,id) / (tlag(ip,id) - taup(ip))
        aux5 = aux3 - aux2
        aux6 = aux3 * aux3

        ter1 = 0.5d0 * eptpa(juf+i0,ip) * aux3
        ter2 = auxl(ip,id+3) * (1.d0 - (aux3 - 1.d0) / aux1)

        ter3 = -aux6 + (aux6 - 1.d0) / (2.d0 * aux1)
        ter4 = 1.d0 - (aux6 - 1.d0) / (2.d0 * aux1)
        sige = (ter3 * bx(ip,id,nor-1) + ter4 * bx(ip,id,nor))    &
             * (1.d0 / (1.d0 - aux6) )

        ter5 = 0.5d0 * tlag(ip,id) * (1.d0 - aux6)

        eptp(juf+i0,ip) =   pepa(jtsuf(id),ip) + ter1 + ter2      &
                          + sige * sqrt(ter5) * vagaus(ip,id)

!---> Calcul de Up :
!     --------------

        ter1 = 0.5d0 * eptpa(jup+i0,ip) * aux2
        ter2 = 0.5d0 * eptpa(juf+i0,ip) * aux4 * aux5
        ter3 = auxl(ip,id+3)*                                     &
           ( 1.d0 - ((tlag(ip,id) + taup(ip)) /dtp) *(1.d0-aux2)  &
             + (tlag(ip,id) / dtp) * aux4 * aux5)                 &
             + auxl(ip,id) * (1.d0 - (aux2 - 1.d0) / aux0)

        tapn  = pepa(jtaux,ip)
        aux7  = exp(-dtp / tapn)
        aux8  = 1.d0 - aux3 * aux7
        aux9  = 1.d0 - aux6
        aux10 = 1.d0 - aux7 * aux7

        aux11 = tapn / (tlag(ip,id) + tapn)
        aux12 = tlag(ip,id) / (tlag(ip,id) - tapn)
        aux14 = tlag(ip,id) - tapn

        aux15 = tlag(ip,id) * (1.d0-aux3)
        aux16 = tapn * (1.d0 - aux7)
        aux17 = sige * sige * aux12 * aux12

        aux18 = 0.5d0 * tlag(ip,id) * aux9
        aux19 = 0.5d0 * tapn * aux10
        aux20 = tlag(ip,id) * aux11 * aux8

! ---> calcul de la matrice de correlation

        gamma2 = sige * sige * aux18

        grgam2 = aux17 * (aux18 - 2.d0 * aux20 + aux19)

        gagam = sige * sige * aux12 * (aux18 - aux20)

        omega2 = aux17 * (aux14 *                                 &
               (aux14 * dtp - 2.d0 * tlag(ip,id) * aux15          &
               + 2.d0 * tapn * aux16)                             &
               + tlag(ip,id) * tlag(ip,id) * aux18                &
               + tapn * tapn * aux19                              &
               - 2.d0 * tlag(ip,id) * tapn * aux20)

        omegam = aux14 * (1.d0-aux3) - aux18 + tapn * aux11 * aux8
        omegam = omegam * sige * sige * aux12 * tlag(ip,id)

        gaome  = aux17 * (aux14 * (aux15 - aux16) - tlag(ip,id)   &
            * aux18 - tapn * aux19 + tapn * tlag(ip,id) * aux8 )

! ---> simulation du vecteur Gaussien

        p11 = sqrt( max (zero,gamma2) )

        if (p11.gt.epzero) then
          p21 = omegam / p11
          p22 = omega2 - p21 * p21
          p22 = sqrt( max (zero,p22) )
        else
          p21 = 0.d0
          p22 = 0.d0
        endif

        if (p11.gt.epzero) then
          p31 = gagam / p11
        else
          p31 = 0.d0
        endif

        if (p22.gt.epzero) then
          p32 = (gaome - p31 * p21) / p22
        else
          p32 = 0.d0
        endif

        p33 = grgam2 - p31 * p31 - p32 * p32
        p33 = sqrt( max(zero,p33) )

        ter4 = p31 * vagaus(ip,id)                                &
             + p32 * vagaus(ip,3+id) + p33 * vagaus(ip,6+id)

! ---> Calcul des Termes dans le cas du mouvement Brownien

        if ( lamvbr .eq. 1 ) then
          tbriu = terbru(ip)*brgaus(ip,3+id)
        else
          tbriu = 0.d0
        endif

! ---> finalisation de l'ecriture

        eptp(jup+i0,ip) =   pepa(jtsup(id),ip)                    &
                          + ter1 + ter2 + ter3 + ter4 + tbriu

      endif

    enddo

  enddo

endif

deallocate(auxl)

!----
! FIN
!----

return
end subroutine lages2
