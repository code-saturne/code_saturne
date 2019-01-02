!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
!> \file solmoy.f90
!> \brief Atmospheric soil module - Initialize ground level parameters from land use
!
!> \brief
!>   calcul des coefficients du modele d'interface sol-atmosphere
!>   on connait :
!>-    - la surface occupee par chaque categorie de sol (pourcentage)
!>       dans chaque maille (tableau indsol)
!>-     - les valeurs des coefficients du modele d'interface sol-atm.
!>       pour chaque categorie de sol (tableaux rugdyn,...,csol)
!>
!>    on calcule pour chaque maille :
!>       la moyenne des coefficients des categories de sol
!>                  ponderee par la surface de chaque categorie
!>    ainsi par exemple si on a 3 types de sol (eau, foret, bati)
!>-     albedo(i,j) = albedo(eau  ) * %_surface_eau  (i,j)
!>                 + albedo(foret) * %_surface_foret(i,j)
!>                 + albedo(bati ) * %_surface_bati (i,j)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[out]   ierreu   code error
!-------------------------------------------------------------------------------

subroutine solmoy ( ierreu )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use atincl
use atsoil
use mesh

implicit none

!===============================================================================

integer          ierreu

! Local variables

integer iirugdy,iirugth,iialbed,iiemiss,iicsol,iiveget
integer iic1w,iic2w,iir1,iir2
integer ifac,n
double precision codinv
double precision rugdij,rugtij,albeij,emisij
double precision vegeij,c1wij,c2wij,csolij
double precision r1ij,r2ij
double precision solmax(10),solmea(10),solmin(10)
character(len=12) ::    solnom(10)

!  ================================================================
!     1) initialisations
!  ================================================================

codinv = -999.d0
do ifac = 1, nfmodsol
  solution_sol(ifac)%constantes%rugdyn = codinv
  solution_sol(ifac)%constantes%rugthe = codinv
  solution_sol(ifac)%constantes%albedo = codinv
  solution_sol(ifac)%constantes%emissi = codinv
  solution_sol(ifac)%constantes%vegeta = codinv
  solution_sol(ifac)%constantes%c1w    = codinv
  solution_sol(ifac)%constantes%c2w    = codinv
  solution_sol(ifac)%constantes%csol   = codinv
  solution_sol(ifac)%constantes%r1     = codinv
  solution_sol(ifac)%constantes%r2     = codinv
enddo

iirugdy = 1
iirugth = 2
iialbed = 3
iiemiss = 4
iicsol  = 5
iiveget = 6
iic1w   = 7
iic2w   = 8
iir1    = 9
iir2    = 10

solnom(iirugdy)='z0 dynamique'
solnom(iirugth)='z0 thermique'
solnom(iialbed)='albedo      '
solnom(iiemiss)='emissivite  '
solnom(iicsol )='csol (x1e+6)'
solnom(iiveget)='vegetation  '
solnom(iic1w  )='c1w         '
solnom(iic2w  )='c2w         '
solnom(iir1   )='r1          '
solnom(iir2   )='r2          '

!  ================================================================
!     2) calcul des coefficients pour chaque maille
!  ================================================================



do ifac = 1, nfmodsol
  rugdij = zero
  rugtij = zero
  albeij = zero
  emisij = zero
  csolij = zero
  vegeij = zero
  c1wij  = zero
  c2wij  = zero
  r1ij   = zero
  r2ij   = zero

  !          zrrel = z(2)*(z(km)-CDGFBO(3,IFAC))/z(km)

  do n = 1, nbrsol
    rugdij = rugdij + tab_sol(n)%rugdyn*float(pourcent_sol(ifac,n))/100.d0
    rugtij = rugtij + tab_sol(n)%rugthe*float(pourcent_sol(ifac,n))/100.d0
    albeij = albeij + tab_sol(n)%albedo*float(pourcent_sol(ifac,n))/100.d0
    emisij = emisij + tab_sol(n)%emissi*float(pourcent_sol(ifac,n))/100.d0
    csolij = csolij + tab_sol(n)%csol  *float(pourcent_sol(ifac,n))/100.d0
    vegeij = vegeij + tab_sol(n)%vegeta*float(pourcent_sol(ifac,n))/100.d0
    c1wij  = c1wij  + tab_sol(n)%c1w   *float(pourcent_sol(ifac,n))/100.d0
    c2wij  = c2wij  + tab_sol(n)%c2w   *float(pourcent_sol(ifac,n))/100.d0
    r1ij   = r1ij   + tab_sol(n)%r1    *float(pourcent_sol(ifac,n))/100.d0
    r2ij   = r2ij   + tab_sol(n)%r2    *float(pourcent_sol(ifac,n))/100.d0
  enddo

  solution_sol(ifac)%constantes%rugdyn = rugdij
  solution_sol(ifac)%constantes%rugthe = rugtij
  solution_sol(ifac)%constantes%albedo = albeij
  solution_sol(ifac)%constantes%emissi = emisij
  solution_sol(ifac)%constantes%csol   = csolij
  solution_sol(ifac)%constantes%vegeta = vegeij
  solution_sol(ifac)%constantes%c1w    = c1wij
  solution_sol(ifac)%constantes%c2w    = c2wij
  solution_sol(ifac)%constantes%r1     = r1ij
  solution_sol(ifac)%constantes%r2     = r2ij

  ! Pour temperatures profondes dans un premier temps on initialise a tprini
  solution_sol(ifac)%constantes%tprof = tprini
enddo

!  ================================================================
!     3) controle
!  ================================================================

ierreu = 0
do ifac = 1, nfmodsol
  if(solution_sol(ifac)%constantes%rugdyn .eq. codinv) ierreu = ierreu + 1
  if(solution_sol(ifac)%constantes%rugthe .eq. codinv) ierreu = ierreu + 1
  if(solution_sol(ifac)%constantes%albedo .eq. codinv) ierreu = ierreu + 1
  if(solution_sol(ifac)%constantes%emissi .eq. codinv) ierreu = ierreu + 1
  if(solution_sol(ifac)%constantes%csol   .eq. codinv) ierreu = ierreu + 1
  if(solution_sol(ifac)%constantes%vegeta .eq. codinv) ierreu = ierreu + 1
  if(solution_sol(ifac)%constantes%c1w    .eq. codinv) ierreu = ierreu + 1
  if(solution_sol(ifac)%constantes%c2w    .eq. codinv) ierreu = ierreu + 1
  if(solution_sol(ifac)%constantes%r1     .eq. codinv) ierreu = ierreu + 1
  if(solution_sol(ifac)%constantes%r2     .eq. codinv) ierreu = ierreu + 1
enddo

! impression eventuelle d'un message d'erreur

if(ierreu.ne.0) then
  write(nfecra,9999)
  write(nfecra,9991)ierreu

  ! ou impression de controle

else
  do n = 1, 10
    solmin(n) = +999999.d0
    solmea(n) = 0.d0
    solmax(n) = -999999.d0
  enddo
  do ifac = 1, nfmodsol, 1
    if (solution_sol(ifac)%constantes%rugdyn .gt. solmax(1)) solmax(1) &
         = solution_sol(ifac)%constantes%rugdyn
    if (solution_sol(ifac)%constantes%rugthe .gt. solmax(2)) solmax(2) &
         = solution_sol(ifac)%constantes%rugthe
    if (solution_sol(ifac)%constantes%albedo .gt. solmax(3)) solmax(3) &
         = solution_sol(ifac)%constantes%albedo
    if (solution_sol(ifac)%constantes%emissi .gt. solmax(4)) solmax(4) &
         = solution_sol(ifac)%constantes%emissi
    if (solution_sol(ifac)%constantes%csol   .gt. solmax(5)) solmax(5) &
         = solution_sol(ifac)%constantes%csol
    if (solution_sol(ifac)%constantes%vegeta .gt. solmax(6)) solmax(6) &
         = solution_sol(ifac)%constantes%vegeta
    if (solution_sol(ifac)%constantes%c1w    .gt. solmax(7)) solmax(7) &
         = solution_sol(ifac)%constantes%c1w
    if (solution_sol(ifac)%constantes%c2w    .gt. solmax(8)) solmax(8) &
         = solution_sol(ifac)%constantes%c2w
    if (solution_sol(ifac)%constantes%r1     .gt. solmax(9)) solmax(9) &
         = solution_sol(ifac)%constantes%r1
    if (solution_sol(ifac)%constantes%r2     .gt. solmax(10))solmax(10)&
         = solution_sol(ifac)%constantes%r2
    if (solution_sol(ifac)%constantes%rugdyn .lt. solmin(1)) solmin(1) &
         = solution_sol(ifac)%constantes%rugdyn
    if (solution_sol(ifac)%constantes%rugthe .lt. solmin(2)) solmin(2) &
         = solution_sol(ifac)%constantes%rugthe
    if (solution_sol(ifac)%constantes%albedo .lt. solmin(3)) solmin(3) &
         = solution_sol(ifac)%constantes%albedo
    if (solution_sol(ifac)%constantes%emissi .lt. solmin(4)) solmin(4) &
         = solution_sol(ifac)%constantes%emissi
    if (solution_sol(ifac)%constantes%csol   .lt. solmin(5)) solmin(5) &
         = solution_sol(ifac)%constantes%csol
    if (solution_sol(ifac)%constantes%vegeta .lt. solmin(6)) solmin(6) &
         = solution_sol(ifac)%constantes%vegeta
    if (solution_sol(ifac)%constantes%c1w    .lt. solmin(7)) solmin(7) &
         = solution_sol(ifac)%constantes%c1w
    if (solution_sol(ifac)%constantes%c2w    .lt. solmin(8)) solmin(8) &
         = solution_sol(ifac)%constantes%c2w
    if (solution_sol(ifac)%constantes%r1     .lt. solmin(9)) solmin(9) &
         = solution_sol(ifac)%constantes%r1
    if (solution_sol(ifac)%constantes%r2     .lt. solmin(10))solmin(10)&
         = solution_sol(ifac)%constantes%r2
    solmea(1) = solmea(1) + solution_sol(ifac)%constantes%rugdyn
    solmea(2) = solmea(2) + solution_sol(ifac)%constantes%rugthe
    solmea(3) = solmea(3) + solution_sol(ifac)%constantes%albedo
    solmea(4) = solmea(4) + solution_sol(ifac)%constantes%emissi
    solmea(5) = solmea(5) + solution_sol(ifac)%constantes%csol
    solmea(6) = solmea(6) + solution_sol(ifac)%constantes%vegeta
    solmea(7) = solmea(7) + solution_sol(ifac)%constantes%c1w
    solmea(8) = solmea(8) + solution_sol(ifac)%constantes%c2w
    solmea(9) = solmea(9) + solution_sol(ifac)%constantes%r1
    solmea(10)= solmea(10)+ solution_sol(ifac)%constantes%r2
  enddo
  do n = 1, 10
    solmea(n)= solmea(n)/float(nfmodsol)
  enddo
  write(nfecra,3001)
  write(nfecra,3002)
  n = iirugdy
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iirugth
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iialbed
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iiemiss
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iicsol
  write(nfecra,3003)solnom(n),solmin(n)*1.d6,solmea(n)*1.d6                     &
       ,solmax(n)*1.d6
  n = iiveget
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iic1w
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iic2w
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iir1
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iir2
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  write(nfecra,3004)
endif

! sortie

return

!--------
! formats
!--------

9999 format(//,5x,'%% erreur solmoy: erreur numero  1')
9991 format( 23x,'initialisation incorrecte des coefficients de ',&
     'l''interface sol atmosphere',                         &
             /,23x,'sur',i9,' valeurs il y en a',i9,' non initialisees')
3001 format(//,8x,' ** ========================================= **',&
             /,8x,' ** interface sol atmosphere                  **',      &
             /,8x,' ** valeurs des constantes calculees          **',      &
             /,8x,' ** ========================================= **',/)
3002 format(   ' *            * minimum* moyenne* maximum*')
3003 format(   ' *',a12,'*',3(f8.4,'*'))
3004 format(/,8x,' ** ========================================= **',//)

end subroutine solmoy
