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

subroutine rayigc &
!================

 (zbas,zz,pz,zzp,pzp,xa,xda,q,u,tco2,ro)

!================================================================================
!  Purpose:
!  --------

!    Atmospheric module subroutine.

!
!   calcul de l'absorption par le gaz carbonique et par l'ozone
!   dans l'infra-rouge
!
!
!---------------------------------------------------------------------------------
! Arguments
!__________________.____._____.__________________________________________________.
! !    nom    !type!mode!                   role                                 !
!_!___________!____!____!________________________________________________________!
! !  zbas     ! r  ! d  ! altitude absolue au niveau k1                          !
! !  zz       ! r  ! d  ! altitude (z ou zc suivant relief ou pas)               !
! !  pz       ! r  ! d  ! pression normalisee a la pression au sol               !
! !  zzp      ! r  ! d  ! altitude niveau intermediaire                          !
! !  pzp      ! r  ! d  ! idem pz au niveau zzp                                  !
! !  xa       ! r  ! r  ! absorption par le CO2 + O3                             !
! !  xda      ! r  ! r  ! absorption differentielle par le CO2 + O3              !
! !  q        ! r  ! d  ! quantite effective d'absorbant pour la vapeur          !
! !           !    !              ! d'eau
! !           !    !    ! de la pression                                         !
! !  u        ! r  ! d  ! grandeur necessaire au calcul de la                    !
! !           !    !    ! transmission par la vapeur d'eau                       !
! !  tco2     ! r  ! d  ! fonction de la temperature                             !
! !  ro       ! r  ! d  ! masse volumique                                        !
!_!___________!____!____!________________________________________________________!

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
use optcal
use cstphy
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use atincl

!===============================================================================

implicit none

!... declaration des variables externes

double precision zbas
double precision zz,pz,zzp,pzp,xa,xda,q,u,tco2,ro

!... declaration des variables internes

integer          iphas
double precision x1,x2,x3,x4,x1c,x2c,x3c,x4c,y1,y2
double precision tauv,dtauv,xx,exn,exnp1,conco2
double precision uco2,duco2,ao,dao
double precision uo3,duo3,xo1,xo2,xo3,xo4
double precision yo1,yo2,ao3,dao3

data x1,x2,x3,x4/1.33,-0.4572,0.26,0.286/
data x1c,x2c,x3c,x4c/-0.00982,0.0676,0.421,0.01022/
data y1,y2/0.0581,0.0546/
data xo1,xo2,xo3,xo4/0.209,7.e-5,0.436,-0.00321/
data yo1,yo2/0.0749,0.0212/

!===============================================================================

iphas = 1

!    concentration en gaz carbonique
conco2 = 3.5d-2

!   1-calcul de th2o dans la bande 15mu du co2
!   ------------------------------------------
if(u.le.20.d0) then
  tauv = x1 + x2*(u + x4)**x3
  dtauv = ro*q*x2*x3*(u + x4)**(x3 - 1.d0)
else
  tauv = 0.33d0 - 0.2754d0*(log10(u) - 1.3011d0)
  dtauv = -0.2754d0/log(10.d0)*ro*q/u
endif

!   2-calcul de l epaisseur optique pour le co2
!   -------------------------------------------
xx = 1.d0 - 0.0065d0*(zz - zbas)/288.15d0
exn = 0.75d0
exnp1 = exn+1.d0
uco2 = -conco2*288.15d0/0.0065d0/(5.31d0*exnp1)*(pz**exnp1                &
     - pzp**exnp1)*(tkelvi/tco2)**(exn/2.d0)

if(uco2.lt.0.d0) uco2 = -uco2
duco2 = conco2*pz**exnp1/xx
duco2 = duco2*(tkelvi/tco2)**(exn/2.d0)
if(uco2.le.1.d0) then
  ao = x1c + x2c*(uco2 + x4c)**x3c
  dao = duco2*x2c*x3c*(uco2 + x4c)**(x3c - 1.d0)
else
  ao = y1 + y2*log10(uco2)
  dao = y2/log(10.d0)*duco2/uco2
endif

!   3-calcul de l epaisseur optique pour O3
!   ---------------------------------------
uo3 = abs(rayuoz(zz) - rayuoz(zzp))
duo3 = raydoz(zz)
if(uo3.le.0.01d0) then
  ao3 = xo1*(uo3 + xo2)**xo3 + xo4
  dao3 = duo3*xo1*(uo3 + xo2)**(xo3 - 1.d0)
else
  ao3 = yo1 + yo2*log10(uo3)
  dao3 = yo2*duo3/log(10.d0)/uo3
endif

!   4- calcul de l'absorption totale (ozone et co2)
!   -----------------------------------------------
xa = tauv*ao + ao3
xda = tauv*dao + dtauv*ao + dao3

return

!   5- aditionnal functions
!   ------------------------

contains

! 5.1. computes ozones concentration for the altitude zh
!------------------------------------------------------

function rayuoz(zh)

implicit none
double precision, intent(in) :: zh  ! absolute altitude
double precision ::  rayuoz
double precision ::  a, b, c

a = 0.4d0
b = 20000.d0
c = 5000.d0

rayuoz = a*(1.d0 + exp(-b/c))/(1.d0 + exp((zh - b)/c))

end function rayuoz


! 5.2. computes dO3/dz for the altitude zh
!-------------------------------------------------------
function raydoz(zh)

implicit none
double precision, intent(in) :: zh  ! absolute altitude
double precision ::  raydoz
double precision ::  a, b, c

a = 0.4d0
b = 20000.d0
c = 5000.d0

raydoz = -a/c*(exp((zh - b)/c))*(1.d0 + exp(-b/c))                        &
       / ((1.d0 + exp((zh - b)/c))**2.d0)

end function raydoz

end subroutine rayigc
