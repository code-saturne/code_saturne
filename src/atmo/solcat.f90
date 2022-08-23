!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
!> \file solcat.f90
!> \brief Soil - atmosphere parameters computed from a "Land use" file

!> \brief ! *   definition des types de sol et des constantes associees
!>   par defaut, on travaille avec un fichier d'occupation du sol
!>   fourni par l'ign.
!>
!>-   le sol est classe soit en 7 categories :
!>    1) eau
!>    2) foret
!>    3) divers
!>    4) sol mineral nu
!>    5) bati diffus
!>    6) bati mixte
!>    7) bati dense
!>-   soit en 5 categories :
!>    1) eau
!>    2) foret
!>    3) divers
!>    4) sol mineral nu
!>    5) bati
!>
!>   l'utilisateur peut modifier :
!>     - les valeurs des constantes prises par defaut
!>        (par exemple la rugosite de la foret)
!>     - les types de sol a utiliser
!>        (dans le cas de donnees ne provenant pas de l'ign)

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     iappel        first pass to set default values,
!>                              second pass to perform some checks and log
!-------------------------------------------------------------------------------
subroutine solcat ( iappel )

!==============================================================================
! Module files
!============================================================================

use paramx
use entsor
use atsoil

implicit none

! Arguments

integer iappel

! Local variables

integer ierreu
integer eau,foret,divers,minral,diffus,mixte,dense,bati
integer n
double precision codinv
integer inityp
character(len=50) :: raison
character(len=10) :: inicat

! Initialisation

inicat = 'xxxxxxxx'
inityp = -9
codinv = -999.d0

! First pass, default values according to the choice of number of soils
!----------------------------------------------------------------------

if (iappel.eq.1) then
  do n = 1, nbrsol
    tab_sol(n)%rugdyn = codinv
    tab_sol(n)%rugthe = codinv
    tab_sol(n)%albedo = codinv
    tab_sol(n)%emissi = codinv
    tab_sol(n)%csol   = codinv
    tab_sol(n)%vegeta = codinv
    tab_sol(n)%c1w    = codinv
    tab_sol(n)%c2w    = codinv
    tab_sol(n)%r1     = codinv
    tab_sol(n)%r2     = codinv
  enddo

  do n = 1, nbrsol
    tab_sol(n)%nomcat = inicat
  enddo

  ! note: si vous utilisez d'autres categories de sol (comme prairie),
  !       rajoutez les a la liste deja existante au lieu de supprimer
  !       celles que vous n'utilisez pas ou encore pire d'utiliser un nom
  !       de categorie pour les coefficients d'une autre

  eau    = inityp
  foret  = inityp
  divers = inityp
  minral = inityp
  diffus = inityp
  mixte  = inityp
  dense  = inityp
  bati   = inityp

  ! cas d'un fichier base sur donnees ign en 7 categories

  if (nbrsol.eq.7) then

    ! definition des types de sol utilises et de l'ordre de rangement qui
    ! doit etre le meme que celui du fichier d'occupation du sol utilise
    ! dans lecsol

    eau    = 1
    foret  = 2
    divers = 3
    minral = 4
    diffus = 5
    mixte  = 6
    dense  = 7

    ! cas d'un fichier base sur donnees ign en 5 categories

  elseif (nbrsol.eq.5) then
    eau    = 1
    foret  = 2
    divers = 3
    minral = 4
    bati   = 5
  endif

  ! Soil categories names

  if (eau    .ne. inityp) tab_sol(eau)%nomcat    = 'water   '
  if (foret  .ne. inityp) tab_sol(foret)%nomcat  = 'forest  '
  if (divers .ne. inityp) tab_sol(divers)%nomcat = 'diverse '
  if (minral .ne. inityp) tab_sol(minral)%nomcat = 'mineral '
  if (diffus .ne. inityp) tab_sol(diffus)%nomcat = 'bg diffu'
  if (mixte  .ne. inityp) tab_sol(mixte )%nomcat = 'bg mixte'
  if (dense  .ne. inityp) tab_sol(dense )%nomcat = 'bg dense'
  if (bati   .ne. inityp) tab_sol(bati  )%nomcat = 'building'

  ! valeurs standard des parametres

  if(eau    .ne. inityp)tab_sol(eau)%rugdyn    = 0.0005d0
  if(foret  .ne. inityp)tab_sol(foret)%rugdyn  = 0.800d0
  if(divers .ne. inityp)tab_sol(divers)%rugdyn = 0.100d0
  if(minral .ne. inityp)tab_sol(minral)%rugdyn = 0.0012d0
  if(diffus .ne. inityp)tab_sol(diffus)%rugdyn = 0.250d0
  if(mixte  .ne. inityp)tab_sol(mixte )%rugdyn = 0.600d0
  if(dense  .ne. inityp)tab_sol(dense )%rugdyn = 1.000d0
  if(bati   .ne. inityp)tab_sol(bati  )%rugdyn = 0.600d0

  if(eau    .ne. inityp)tab_sol(eau)%rugthe    = tab_sol(eau)%rugdyn
  if(foret  .ne. inityp)tab_sol(foret)%rugthe  = tab_sol(foret)%rugdyn*exp(-2.d0)
  if(divers .ne. inityp)tab_sol(divers)%rugthe = tab_sol(divers)%rugdyn*exp(-2.d0)
  if(minral .ne. inityp)tab_sol(minral)%rugthe = tab_sol(minral)%rugdyn*exp(-2.d0)
  if(diffus .ne. inityp)tab_sol(diffus)%rugthe = tab_sol(diffus)%rugdyn*exp(-2.d0)
  if(mixte  .ne. inityp)tab_sol(mixte )%rugthe = tab_sol(mixte )%rugdyn*exp(-2.d0)
  if(dense  .ne. inityp)tab_sol(dense )%rugthe = tab_sol(dense )%rugdyn*exp(-2.d0)
  if(bati   .ne. inityp)tab_sol(bati  )%rugthe = tab_sol(bati  )%rugdyn*exp(-2.d0)

  if(eau    .ne. inityp)tab_sol(eau)%albedo    = 0.08d0
  if(foret  .ne. inityp)tab_sol(foret)%albedo  = 0.16d0
  if(divers .ne. inityp)tab_sol(divers)%albedo = 0.20d0
  if(minral .ne. inityp)tab_sol(minral)%albedo = 0.25d0
  if(diffus .ne. inityp)tab_sol(diffus)%albedo = 0.18d0
  if(mixte  .ne. inityp)tab_sol(mixte )%albedo = 0.18d0
  if(dense  .ne. inityp)tab_sol(dense )%albedo = 0.18d0
  if(bati   .ne. inityp)tab_sol(bati  )%albedo = 0.18d0

  if(eau    .ne. inityp)tab_sol(eau)%emissi    = 0.980d0
  if(foret  .ne. inityp)tab_sol(foret)%emissi  = 0.950d0
  if(divers .ne. inityp)tab_sol(divers)%emissi = 0.940d0
  if(minral .ne. inityp)tab_sol(minral)%emissi = 0.965d0
  if(diffus .ne. inityp)tab_sol(diffus)%emissi = 0.880d0
  if(mixte  .ne. inityp)tab_sol(mixte )%emissi = 0.880d0
  if(dense  .ne. inityp)tab_sol(dense )%emissi = 0.880d0
  if(bati   .ne. inityp)tab_sol(bati  )%emissi = 0.880d0

  if(eau    .ne. inityp)tab_sol(eau)%vegeta    = 0.00d0
  if(foret  .ne. inityp)tab_sol(foret)%vegeta  = 1.00d0
  if(divers .ne. inityp)tab_sol(divers)%vegeta = 1.00d0
  if(minral .ne. inityp)tab_sol(minral)%vegeta = 0.00d0
  if(diffus .ne. inityp)tab_sol(diffus)%vegeta = 0.50d0
  if(mixte  .ne. inityp)tab_sol(mixte )%vegeta = 0.25d0
  if(dense  .ne. inityp)tab_sol(dense )%vegeta = 0.00d0
  if(bati   .ne. inityp)tab_sol(bati  )%vegeta = 0.25d0

  if(eau    .ne. inityp)tab_sol(eau)%csol    =  7.6d-06
  if(foret  .ne. inityp)tab_sol(foret)%csol  = 11.0d-06
  if(divers .ne. inityp)tab_sol(divers)%csol = 11.0d-06
  if(minral .ne. inityp)tab_sol(minral)%csol =  5.0d-06
  if(dense  .ne. inityp)tab_sol(dense )%csol =  3.9d-06
  if(diffus .ne. inityp)                                                     &
       tab_sol(diffus)%csol = tab_sol(foret)%csol*tab_sol(diffus)%vegeta +   &
       tab_sol(dense)%csol*(1.d0-tab_sol(diffus)%vegeta)
  if(mixte  .ne. inityp)                                                     &
       tab_sol(mixte )%csol =  tab_sol(foret )%csol*tab_sol(mixte )%vegeta + &
       tab_sol(dense )%csol*(1.d0-tab_sol(mixte )%vegeta)
  if(bati  .ne. inityp)                                                      &
       tab_sol(bati  )%csol =  tab_sol(foret )%csol*tab_sol(bati  )%vegeta + &
       3.9d-06*(1.d0-tab_sol(bati  )%vegeta)
  if(eau    .ne. inityp)tab_sol(eau)%c1w    = 100.0d0
  if(foret  .ne. inityp)tab_sol(foret)%c1w  = 18.d0*tab_sol(foret)%vegeta + 2.d0
  if(divers .ne. inityp)tab_sol(divers)%c1w = 18.d0*tab_sol(divers)%vegeta + 2.d0
  if(minral .ne. inityp)tab_sol(minral)%c1w = 18.d0*tab_sol(minral)%vegeta + 2.d0
  if(diffus .ne. inityp)tab_sol(diffus)%c1w = 18.d0*tab_sol(diffus)%vegeta + 2.d0
  if(mixte  .ne. inityp)tab_sol(mixte )%c1w = 18.d0*tab_sol(mixte )%vegeta + 2.d0
  if(dense  .ne. inityp)tab_sol(dense )%c1w = 18.d0*tab_sol(dense )%vegeta + 2.d0
  if(bati   .ne. inityp)tab_sol(bati  )%c1w = 18.d0*tab_sol(bati  )%vegeta + 2.d0

  if(eau    .ne. inityp)tab_sol(eau)%c2w    = 1.00d0
  if(foret  .ne. inityp)tab_sol(foret)%c2w  = 0.20d0
  if(divers .ne. inityp)tab_sol(divers)%c2w = 0.20d0
  if(minral .ne. inityp)tab_sol(minral)%c2w = 0.20d0
  if(diffus .ne. inityp)tab_sol(diffus)%c2w = 0.20d0
  if(mixte  .ne. inityp)tab_sol(mixte )%c2w = 0.20d0
  if(dense  .ne. inityp)tab_sol(dense )%c2w = 0.20d0
  if(bati   .ne. inityp)tab_sol(bati  )%c2w = 0.20d0

  if(eau    .ne. inityp)tab_sol(eau)%r1    = 0.d0
  if(foret  .ne. inityp)tab_sol(foret)%r1  = 0.d0
  if(divers .ne. inityp)tab_sol(divers)%r1 = 0.d0
  if(minral .ne. inityp)tab_sol(minral)%r1 = 0.d0
  if(dense  .ne. inityp)tab_sol(dense )%r1 = 30.d0
  if(diffus .ne. inityp)tab_sol(diffus)%r1 = 10.d0
  if(mixte  .ne. inityp)tab_sol(mixte )%r1 = 15.d0
  if(bati   .ne. inityp)tab_sol(bati  )%r1 = 15.d0

  if(eau    .ne. inityp)tab_sol(eau)%r2    = 0.d0
  if(foret  .ne. inityp)tab_sol(foret)%r2  = 0.d0
  if(divers .ne. inityp)tab_sol(divers)%r2 = 0.d0
  if(minral .ne. inityp)tab_sol(minral)%r2 = 0.d0
  if(dense  .ne. inityp)tab_sol(dense )%r2 = 2.0d0
  if(diffus .ne. inityp)tab_sol(diffus)%r2 = 2.0d0/3.d0
  if(mixte  .ne. inityp)tab_sol(mixte )%r2 = 1.d0
  if(bati   .ne. inityp)tab_sol(bati  )%r2 = 1.0d0

endif

! Second pass: log and checks
!----------------------------

if (iappel.eq.2) then
  ! Log
  write(nfecra,2000)
  write(nfecra,2001)
  do n = 1, nbrsol
    write(nfecra,2002) tab_sol(n)%nomcat, tab_sol(n)%rugdyn,tab_sol(n)%rugthe, &
      tab_sol(n)%albedo, tab_sol(n)%emissi,1.d+06*tab_sol(n)%csol,          &
      tab_sol(n)%vegeta, tab_sol(n)%c1w, tab_sol(n)%c2w, tab_sol(n)%r1,     &
      tab_sol(n)%r2
  enddo
  write(nfecra,2012)

  ! controle

  ierreu = nbrsol
  raison = ' Wrong soil cofficients              '
  do n = 1, nbrsol
    if (tab_sol(n)%rugdyn.ne.codinv .and. tab_sol(n)%rugthe.ne.codinv .and.   &
        tab_sol(n)%albedo.ne.codinv .and. tab_sol(n)%emissi.ne.codinv .and.   &
        tab_sol(n)%c1w   .ne.codinv .and. tab_sol(n)%c2w   .ne.codinv .and.   &
        tab_sol(n)%csol  .ne.codinv .and. tab_sol(n)%r1    .ne.codinv .and.   &
        tab_sol(n)%r2    .ne.codinv) ierreu = ierreu - 1
  enddo

  ! Error message if needed
  if (ierreu.ne.0) then
    write(nfecra,9999) ierreu
    write(nfecra,9990) raison

    call csexit(1)
  endif

endif

!--------
! Formats
!--------

9999 format(//,5x,'%% Error in solcat: number of soil with an error = ',i2)
9990 format( 22x,a50)

2000 format(//,   &
            ' Soil-atmosphere interface model',/,   &
            ' Values of tabulated coefficients',/,   &
            ' --------------------------------',/)
2001 format(' Name      z0 dyn   z0 th    albedo   emissi   ', &
            'Cs(e+6)  vegeta   c1w      c2w      ',          &
            'r1       r2')
2002 format(' ',a8,10(f8.4,' '))
2012 format(' --------------------------------',//)

!===============================================================================
! 8. SORTIE
!===============================================================================

return
end subroutine solcat
