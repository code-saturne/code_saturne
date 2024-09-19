!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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
subroutine solcat ( iappel ) &
 bind(C, name='cs_f_solcat')

!==============================================================================
! Module files
!============================================================================

use paramx
use entsor
use atincl
use atsoil

use, intrinsic :: iso_c_binding

implicit none

interface

  subroutine cs_f_atmo_soil_init_arrays(nb_col, p_csol, p_rugdyn, p_rugthe, &
                                        p_albedo, p_emissi, p_vegeta,       &
                                        p_c1w, p_c2w, p_r1, p_r2)           &
    bind(C, name='cs_f_atmo_soil_init_arrays')
    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(in) :: nb_col
    type(c_ptr), intent(out) :: p_csol, p_rugdyn, p_rugthe
    type(c_ptr), intent(out) :: p_albedo, p_emissi, p_vegeta
    type(c_ptr), intent(out) :: p_c1w, p_c2w, p_r1, p_r2
  end subroutine cs_f_atmo_soil_init_arrays

end interface

procedure() :: csexit

! Arguments

integer(c_int), value :: iappel

! Local variables

integer ierreu
integer eau,foret,divers,minral,diffus,mixte,dense,bati
integer n
integer error, n_elts
integer, dimension(:), pointer :: elt_ids
double precision codinv
integer inityp
character(len=50) :: raison
character(len=10) :: nomcat(8)
type(c_ptr) :: c_csol, c_rugdyn, c_rugthe
type(c_ptr) :: c_albedo, c_emissi, c_vegeta
type(c_ptr) :: c_c1w, c_c2w, c_r1, c_r2

DATA nomcat/'','','','','','','',''/

! Initialisation

inityp = -9
codinv = -999.d0

! Get the number of soil zones and their category
! Then allocate the table values
! Only when not called from atmsol.f90 (called from cs_setup.c)
if (iappel.eq.1) then
  call atmo_get_soil_zone(n_elts, nbrsol, elt_ids)

  ! Allocation of table values
  call cs_f_atmo_soil_init_arrays(nbrsol, c_csol, c_rugdyn, c_rugthe,  &
                                  c_albedo, c_emissi, c_vegeta,        &
                                  c_c1w, c_c2w, c_r1, c_r2)

  call c_f_pointer(c_r1, soil_cat_r1, [nbrsol])
  call c_f_pointer(c_r2, soil_cat_r2, [nbrsol])
  call c_f_pointer(c_c1w, soil_cat_c1w, [nbrsol])
  call c_f_pointer(c_c2w, soil_cat_c2w, [nbrsol])
  call c_f_pointer(c_csol, csol, [nbrsol])
  call c_f_pointer(c_rugdyn, rugdyn, [nbrsol])
  call c_f_pointer(c_rugthe, rugthe, [nbrsol])
  call c_f_pointer(c_albedo, soil_cat_albedo, [nbrsol])
  call c_f_pointer(c_emissi, soil_cat_emissi, [nbrsol])
  call c_f_pointer(c_vegeta, soil_cat_vegeta, [nbrsol])
endif

! First pass, default values according to the choice of number of soils
!----------------------------------------------------------------------
if (iappel.eq.1) then
  do n = 1, nbrsol
    csol(n)   = codinv
    rugdyn(n) = codinv
    rugthe(n) = codinv
    soil_cat_albedo(n) = codinv
    soil_cat_emissi(n) = codinv
    soil_cat_vegeta(n) = codinv
    soil_cat_c1w(n)    = codinv
    soil_cat_c2w(n)    = codinv
    soil_cat_r1(n)     = codinv
    soil_cat_r2(n)     = codinv
  enddo

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

  if (eau    .ne. inityp) nomcat(eau)    = 'water'
  if (foret  .ne. inityp) nomcat(foret)  = 'forest'
  if (divers .ne. inityp) nomcat(divers) = 'diverse'
  if (minral .ne. inityp) nomcat(minral) = 'mineral'
  if (diffus .ne. inityp) nomcat(diffus) = 'bg diffu'
  if (mixte  .ne. inityp) nomcat(mixte)  = 'bg mixte'
  if (dense  .ne. inityp) nomcat(dense)  = 'bg dense'
  if (bati   .ne. inityp) nomcat(bati)   = 'building'

  ! note: si vous utilisez d'autres categories de sol (comme prairie),
  !       rajoutez les a la liste deja existante au lieu de supprimer
  !       celles que vous n'utilisez pas ou encore pire d'utiliser un nom
  !       de categorie pour les coefficients d'une autre

  ! valeurs standard des parametres

  if(eau    .ne. inityp)rugdyn(eau)    = 0.0005d0
  if(foret  .ne. inityp)rugdyn(foret)  = 0.800d0
  if(divers .ne. inityp)rugdyn(divers) = 0.100d0
  if(minral .ne. inityp)rugdyn(minral) = 0.0012d0
  if(diffus .ne. inityp)rugdyn(diffus) = 0.250d0
  if(mixte  .ne. inityp)rugdyn(mixte)  = 0.600d0
  if(dense  .ne. inityp)rugdyn(dense)  = 1.000d0
  if(bati   .ne. inityp)rugdyn(bati)   = 0.600d0

  if(eau    .ne. inityp)rugthe(eau)    = rugdyn(eau)
  if(foret  .ne. inityp)rugthe(foret)  = rugdyn(foret)*exp(-2.d0)
  if(divers .ne. inityp)rugthe(divers) = rugdyn(divers)*exp(-2.d0)
  if(minral .ne. inityp)rugthe(minral) = rugdyn(minral)*exp(-2.d0)
  if(diffus .ne. inityp)rugthe(diffus) = rugdyn(diffus)*exp(-2.d0)
  if(mixte  .ne. inityp)rugthe(mixte)  = rugdyn(mixte)*exp(-2.d0)
  if(dense  .ne. inityp)rugthe(dense)  = rugdyn(dense)*exp(-2.d0)
  if(bati   .ne. inityp)rugthe(bati)   = rugdyn(bati)*exp(-2.d0)

  if(eau    .ne. inityp)soil_cat_albedo(eau)    = 0.08d0
  if(foret  .ne. inityp)soil_cat_albedo(foret)  = 0.16d0
  if(divers .ne. inityp)soil_cat_albedo(divers) = 0.20d0
  if(minral .ne. inityp)soil_cat_albedo(minral) = 0.25d0
  if(diffus .ne. inityp)soil_cat_albedo(diffus) = 0.18d0
  if(mixte  .ne. inityp)soil_cat_albedo(mixte)  = 0.18d0
  if(dense  .ne. inityp)soil_cat_albedo(dense)  = 0.18d0
  if(bati   .ne. inityp)soil_cat_albedo(bati)   = 0.18d0

  if(eau    .ne. inityp)soil_cat_emissi(eau)    = 0.980d0
  if(foret  .ne. inityp)soil_cat_emissi(foret)  = 0.950d0
  if(divers .ne. inityp)soil_cat_emissi(divers) = 0.940d0
  if(minral .ne. inityp)soil_cat_emissi(minral) = 0.965d0
  if(diffus .ne. inityp)soil_cat_emissi(diffus) = 0.880d0
  if(mixte  .ne. inityp)soil_cat_emissi(mixte)  = 0.880d0
  if(dense  .ne. inityp)soil_cat_emissi(dense)  = 0.880d0
  if(bati   .ne. inityp)soil_cat_emissi(bati)   = 0.880d0

  if(eau    .ne. inityp)soil_cat_vegeta(eau)    = 0.00d0
  if(foret  .ne. inityp)soil_cat_vegeta(foret)  = 1.00d0
  if(divers .ne. inityp)soil_cat_vegeta(divers) = 1.00d0
  if(minral .ne. inityp)soil_cat_vegeta(minral) = 0.00d0
  if(diffus .ne. inityp)soil_cat_vegeta(diffus) = 0.50d0
  if(mixte  .ne. inityp)soil_cat_vegeta(mixte)  = 0.25d0
  if(dense  .ne. inityp)soil_cat_vegeta(dense)  = 0.00d0
  if(bati   .ne. inityp)soil_cat_vegeta(bati)   = 0.25d0

  if(eau    .ne. inityp)csol(eau)    =  7.6d-06
  if(foret  .ne. inityp)csol(foret)  = 11.0d-06
  if(divers .ne. inityp)csol(divers) = 11.0d-06
  if(minral .ne. inityp)csol(minral) =  5.0d-06
  if(dense  .ne. inityp)csol(dense)  =  3.9d-06
  if(diffus .ne. inityp)                                                     &
       csol(diffus) = csol(foret)*soil_cat_vegeta(diffus) +   &
       csol(dense)*(1.d0-soil_cat_vegeta(diffus))
  if(mixte  .ne. inityp)                                                     &
       csol(mixte) = csol(foret)*soil_cat_vegeta(mixte) + &
       csol(dense)*(1.d0-soil_cat_vegeta(mixte ))
  if(bati  .ne. inityp)                                                      &
       csol(bati) = csol(foret)*soil_cat_vegeta(bati) + &
       3.9d-06*(1.d0-soil_cat_vegeta(bati))
  if(eau    .ne. inityp)soil_cat_c1w(eau)    = 100.0d0
  if(foret  .ne. inityp)soil_cat_c1w(foret)  = 18.d0*soil_cat_vegeta(foret) + 2.d0
  if(divers .ne. inityp)soil_cat_c1w(divers) = 18.d0*soil_cat_vegeta(divers) + 2.d0
  if(minral .ne. inityp)soil_cat_c1w(minral) = 18.d0*soil_cat_vegeta(minral) + 2.d0
  if(diffus .ne. inityp)soil_cat_c1w(diffus) = 18.d0*soil_cat_vegeta(diffus) + 2.d0
  if(mixte  .ne. inityp)soil_cat_c1w(mixte)  = 18.d0*soil_cat_vegeta(mixte) + 2.d0
  if(dense  .ne. inityp)soil_cat_c1w(dense)  = 18.d0*soil_cat_vegeta(dense) + 2.d0
  if(bati   .ne. inityp)soil_cat_c1w(bati)   = 18.d0*soil_cat_vegeta(bati) + 2.d0

  if(eau    .ne. inityp)soil_cat_c2w(eau)    = 1.00d0
  if(foret  .ne. inityp)soil_cat_c2w(foret)  = 0.20d0
  if(divers .ne. inityp)soil_cat_c2w(divers) = 0.20d0
  if(minral .ne. inityp)soil_cat_c2w(minral) = 0.20d0
  if(diffus .ne. inityp)soil_cat_c2w(diffus) = 0.20d0
  if(mixte  .ne. inityp)soil_cat_c2w(mixte)  = 0.20d0
  if(dense  .ne. inityp)soil_cat_c2w(dense)  = 0.20d0
  if(bati   .ne. inityp)soil_cat_c2w(bati)   = 0.20d0

  if(eau    .ne. inityp)soil_cat_r1(eau)    = 0.d0
  if(foret  .ne. inityp)soil_cat_r1(foret)  = 0.d0
  if(divers .ne. inityp)soil_cat_r1(divers) = 0.d0
  if(minral .ne. inityp)soil_cat_r1(minral) = 0.d0
  if(dense  .ne. inityp)soil_cat_r1(dense ) = 30.d0
  if(diffus .ne. inityp)soil_cat_r1(diffus) = 10.d0
  if(mixte  .ne. inityp)soil_cat_r1(mixte) = 15.d0
  if(bati   .ne. inityp)soil_cat_r1(bati)  = 15.d0

  if(eau    .ne. inityp)soil_cat_r2(eau)    = 0.d0
  if(foret  .ne. inityp)soil_cat_r2(foret)  = 0.d0
  if(divers .ne. inityp)soil_cat_r2(divers) = 0.d0
  if(minral .ne. inityp)soil_cat_r2(minral) = 0.d0
  if(dense  .ne. inityp)soil_cat_r2(dense) = 2.0d0
  if(diffus .ne. inityp)soil_cat_r2(diffus) = 2.0d0/3.d0
  if(mixte  .ne. inityp)soil_cat_r2(mixte) = 1.d0
  if(bati   .ne. inityp)soil_cat_r2(bati) = 1.0d0

endif

! Second pass: log and checks
!----------------------------

if (iappel.eq.2) then
  ! Log
  write(nfecra,2000)
  write(nfecra,2001)
  do n = 1, nbrsol
    write(nfecra,2002) nomcat(n), rugdyn(n),rugthe(n), &
      soil_cat_albedo(n), soil_cat_emissi(n),1.d+06*csol(n),          &
      soil_cat_vegeta(n), soil_cat_c1w(n), soil_cat_c2w(n), soil_cat_r1(n),     &
      soil_cat_r2(n)
  enddo
  write(nfecra,2012)

  ! controle

  ierreu = nbrsol
  raison = ' Wrong soil cofficients              '
  do n = 1, nbrsol
    if (rugdyn(n).ne.codinv .and. rugthe(n).ne.codinv .and.   &
        soil_cat_albedo(n).ne.codinv .and. soil_cat_emissi(n).ne.codinv .and.   &
        soil_cat_c1w(n)   .ne.codinv .and. soil_cat_c2w(n)   .ne.codinv .and.   &
        csol(n)  .ne.codinv .and. soil_cat_r1(n)    .ne.codinv .and.   &
        soil_cat_r2(n)    .ne.codinv) ierreu = ierreu - 1
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
            'Cs(e-6)  vegeta   c1w      c2w      ',          &
            'r1       r2')
2002 format(' ',a8,9(f8.4,' '),f8.4)
2012 format(' --------------------------------',//)

!===============================================================================
! 8. SORTIE
!===============================================================================

return
end subroutine solcat
