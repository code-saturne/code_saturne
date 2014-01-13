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

subroutine rayive &
!================

 (tauv,dtauv,qqqq,xqx,qqqqc,xqc,ro)

!==============================================================================
!  Purpose:
!  --------

!    Atmospheric module subroutine.

! *                                                                    *
! *  calcul de l'absorption par le gaz carbonique et par l'ozone       *
! *  dans l'infra-rouge                                                *
! *                                                                    *
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.__________________________________________________.
! !    nom    !type!mode!                   role                                 !
!_!___________!____!____!________________________________________________________!
! ! tauv      ! r  ! r  ! transmission par la vapeur d'eau et                    !
! !           !    !    ! ses dimeres                                            !
! ! dtauv     ! r  ! r  ! derivee de la transmission par rapport                 !
! !           !    !    ! a l'altitude                                           !
! ! qqqq      ! r  ! d  ! quantite effective d'absorbant pour la vapeur          !
! !           !    !    ! d'eau entre le sol et le niveau z                      !
! ! xqx       ! r  ! d  ! quantite effective d'absorbant pour la vapeur          !
! !           !    !    ! d'eau au niveau z'                                     !
! ! qqqqc     ! r  ! d  ! idem qqqq pour les dimeres                             !
! ! xqc       ! r  ! d  ! idem xqx pour les dimeres                              !
! ! ro        ! r  ! d  ! masse volumique au niveau z                            !
!_!___________!____!____!________________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================

!... declaration des variables externes

double precision tauv,dtauv,qqqq,xqx,qqqqc,xqc,ro

!... declaration des variables internes

double precision u,us,xu,xus,za,dza
double precision a0,a1,a2,a3,a4,b0,b1,b2,b3,b4
double precision as0,as1,as2,bs0,bs1,bs2
double precision abs0,dabs0
double precision n1,dn1,d1,dd1,t1,dt1,n2,dn2,d2,dd2,t2,dt2

data a0,a1,a2,a3,a4/7.76192d-7,1.33836d-3,1.66649d-1,2.17686      &
,2.69020/
data b0,b1,b2,b3,b4/7.79097d-7,1.36832d-3,1.79601d-1,2.70573      &
,5.15119/
data as0,as1,as2/1.5075d-2,-3.6185d-2,1.9245d-2/
data bs0,bs1,bs2/1.5075d-2,1.9547d-1,7.5271d-1/

!===============================================================================

u = qqqq/10.d0
us = qqqqc/10.d0
xu = xqx/10.d0
xus = xqc/10.d0

!  emissivite de la vapeur d'eau

if(u.ge.0.01d0) then
  za = 0.24d0*log10(u+0.01d0) + 0.622d0
  dza = 0.24d0/(u + 0.01d0)/log(10.d0)
else
  za = 0.846d0*(u + 3.59d-5)**0.243d0 - 6.9d-2
  dza = 0.846d0*0.243d0*(u + 3.59d-5)**(0.243d0 - 1.d0)
endif

!  emissivite des dimeres

n1 = a0 + u*(a1 + u*(a2 + u*(a3 + u*a4)))
dn1 = a1 + u*(2.d0*a2 + u*(3.d0*a3 + u*4.d0*a4))
d1 = b0 + u*(b1 + u*(b2 + u*(b3 + u*(b4 + u))))
dd1 = b1 + u*(2.d0*b2 + u*(3.d0*b3 + u*(4.d0*b4 + u*5.d0)))
t1 = n1/d1

dt1 = dn1/d1 - n1*dd1/d1/d1

if(us.le.0.5d0) then
  n2 = as0 + us*(as1 + us*as2)
  dn2 = as1 + us*2.d0*as2
  d2 = bs0 + us*(bs1 + us*(bs2 + us))
  dd2 = bs1 + us*(2.d0*bs2 + 3.d0*us)
  t2 = n2/d2

  dt2 = dn2/d2 - n2*dd2/d2/d2
else
  t2 = 0.d0
  dt2 = 0.d0
endif

!  transmission totale

abs0 = za + 0.4614d0*t1*(1.d0 - t2)
tauv = 1.d0 - abs0
dabs0 = (dza*xu + 0.4614d0*(dt1*xu*(1.d0 - t2) - t1*dt2*xus))*ro
dtauv = -dabs0

return
end subroutine rayive
