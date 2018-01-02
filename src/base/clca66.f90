!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine clca66 &
!================

 ( clsyme , eloglo , alpha )

!===============================================================================
! FONCTION :
! --------
!      CALCUL DE LA MATRICE ALPHA (VOIR CNDRIJ)
! POUR LES CONDITIONS AUX LIMITES EN PAROI ET SYMETRIE EN RIJ-EPS

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! clsyme           !  r ! <-- ! indicateur symetrie = 1., paroi = 0.           !
! eloglo           ! tr ! <-- ! matrice de changement de base                  !
!  (6,6)           !    !     !  local -> global                               !
! alpha            ! tr ! --> ! voir cndrij                                    !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision clsyme
double precision eloglo(3,3),alpha(6,6)

! Local variables

integer ii,jj,kk,pp,jj1,jj2

!===============================================================================

! Initialize variables to avoid compiler warnings

jj1 = 0
jj2 = 0

!===============================================================================
! CALCUL DE ALPHA(I,J) POUR I DANS [1,3] ET J DANS [1,3]  : 9 TERMES
!===============================================================================

 do ii = 1,3
   do jj = 1,3

     alpha(ii,jj) =   eloglo(ii,1)**2*eloglo(jj,1)**2 +           &
                      eloglo(ii,2)**2*eloglo(jj,2)**2 +           &
                      eloglo(ii,3)**2*eloglo(jj,3)**2 +           &
   2.d0*clsyme*eloglo(ii,1)*eloglo(ii,3)*eloglo(jj,1)*eloglo(jj,3)

   enddo
 enddo

!===============================================================================
! CALCUL DE ALPHA(I,J) POUR I DANS [1,3] ET J DANS [4,6]  : 9 TERMES
!===============================================================================

 do ii=1,3
   do jj=1,3

     if (jj.eq.1) then
       kk=1
       pp=2
     else if (jj.eq.2) then
       kk=2
       pp=3
     else if (jj.eq.3) then
       kk=1
       pp=3
     endif

     alpha(ii,jj+3) = 2.d0*(                                      &
       eloglo(ii,1)**2*eloglo(kk,1)*eloglo(pp,1) +                &
       eloglo(ii,2)**2*eloglo(kk,2)*eloglo(pp,2) +                &
       eloglo(ii,3)**2*eloglo(kk,3)*eloglo(pp,3) +                &
       clsyme*eloglo(ii,3)*eloglo(ii,1)*                          &
        (eloglo(kk,1)*eloglo(pp,3)+eloglo(kk,3)*eloglo(pp,1)))

   enddo
 enddo

!===============================================================================
! CALCUL DE ALPHA(I,J) POUR I DANS [4,6] ET J DANS [1,3]  : 9 TERMES
!===============================================================================

 do ii=1,3
   do jj=1,3

     if (ii.eq.1) then
       kk=1
       pp=2
     else if (ii.eq.2) then
       kk=2
       pp=3
     else if (ii.eq.3) then
       kk=1
       pp=3
     endif
     alpha(ii+3,jj) =                                             &
       eloglo(kk,1)*eloglo(pp,1)*eloglo(jj,1)**2 +                &
       eloglo(kk,2)*eloglo(pp,2)*eloglo(jj,2)**2 +                &
       eloglo(kk,3)*eloglo(pp,3)*eloglo(jj,3)**2 +                &
       clsyme*eloglo(jj,3)*eloglo(jj,1)*                          &
        (eloglo(kk,1)*eloglo(pp,3)+eloglo(kk,3)*eloglo(pp,1))

   enddo
 enddo


!===============================================================================
! CALCUL DE ALPHA(I,J) POUR I DANS [4,6] ET J DANS [4,6]  : 9 TERMES
!===============================================================================

 do ii=1,3
   do jj=1,3

     if (ii.eq.1) then
       kk=1
       pp=2
     else if (ii.eq.2) then
       kk=2
       pp=3
     else if (ii.eq.3) then
       kk=1
       pp=3
     endif
     if (jj.eq.1) then
       jj1=1
       jj2=2
     else if (jj.eq.2) then
       jj1=2
       jj2=3
     else if (jj.eq.3) then
       jj1=1
       jj2=3
     endif
     alpha(ii+3,jj+3)=                                            &
     2.d0*(eloglo(kk,1)*eloglo(pp,1)*eloglo(jj1,1)*eloglo(jj2,1)+ &
           eloglo(kk,2)*eloglo(pp,2)*eloglo(jj1,2)*eloglo(jj2,2)+ &
           eloglo(kk,3)*eloglo(pp,3)*eloglo(jj1,3)*eloglo(jj2,3))+&
     clsyme*(eloglo(kk,1)*eloglo(pp,3)+eloglo(kk,3)*eloglo(pp,1))*&
         (eloglo(jj1,3)*eloglo(jj2,1)+eloglo(jj1,1)*eloglo(jj2,3))

   enddo
 enddo

 return
 end
