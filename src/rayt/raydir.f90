!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine raydir
!================

!===============================================================================
! Purpose:
! -------

!    Return a quadrature Sn or Tn

!-------------------------------------------------------------------------------
! Arguments

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use radiat

!===============================================================================

implicit none

! Arguments

! Local variables

double precision  vec(10), poids(5), vect(3)
double precision  teta(6, 3)
double precision  phi1, teta1
integer           ii, jj, kk, nray, jray, iia, iib
integer           nquad, nbpoint, niv, ntri, int1, tri

double precision, dimension(:,:), allocatable :: x3d, y3d, z3d
double precision, dimension(:,:), allocatable :: posnod

!==================================================================
! Initialization
!==================================================================

if (i_quadrature.eq.1) then    ! quadrature S4 : 24 directions

  ndirs = 3

elseif (i_quadrature.eq.2) then ! quadrature S6 : 48 directions

  ndirs = 6

elseif (i_quadrature.eq.3) then ! quadrature S8 : 80 directions

  ndirs = 10

elseif (i_quadrature.eq.4) then    ! Quadrature T2 : 32 directions

  ndirs = 4

elseif (i_quadrature.eq.5) then   ! Quadrature T4 : 128 directions

  ndirs = 16

elseif (i_quadrature.eq.6) then   ! Quadrature Tn : 8 n^2 directions

  ndirs = ndirec * ndirec

endif

call init_quadrature(ndirs)

sx = 0.d0
sy = 0.d0
sz = 0.d0
angsol = 0.d0

vec   = 0.d0
poids = 0.d0

!==================================================================
! Quadrature Sn : n(n+2) directions
!==================================================================

if (i_quadrature.eq.1) then    ! quadrature S4 : 24 directions

  vec(1) = 0.2958759
  vec(2) = 0.9082483
  poids(1) = 0.5235987

  sx(1:2) = vec(1)
  sx(3) = vec(2)
  sy(1) = vec(1)
  sy(2) = vec(2)
  sy(3) = vec(1)
  sz(1) = vec(2)
  sz(2:3) = vec(1)
  angsol(1:ndirs) = poids(1)

elseif (i_quadrature.eq.2) then ! quadrature S6 : 48 directions

  vec(1) = 0.1838670
  vec(2) = 0.6950514
  vec(3) = 0.9656013
  poids(1) = 0.1609517
  poids(2) = 0.3626469

  sx(1:3) = vec(1)
  sx(4:5) = vec(2)
  sx(6) = vec(3)
  sy(1) = vec(1)
  sy(2) = vec(2)
  sy(3) = vec(3)
  sy(4) = vec(1)
  sy(5) = vec(2)
  sy(6) = vec(1)
  sz(1) = vec(3)
  sz(2) = vec(2)
  sz(3) = vec(1)
  sz(4) = vec(2)
  sz(5) = vec(1)
  sz(6) = vec(1)
  angsol(1) = poids(1)
  angsol(2) = poids(2)
  angsol(3) = poids(1)
  angsol(4) = poids(2)
  angsol(5) = poids(2)
  angsol(6) = poids(1)

elseif (i_quadrature.eq.3) then ! quadrature S8 : 80 directions

  vec(1) = 0.1422555
  vec(2) = 0.5773503
  vec(3) = 0.8040087
  vec(4) = 0.9795543
  poids(1) = 0.0992284
  poids(2) = 0.1712359
  poids(3) = 0.4617179

  sx(1:4) = vec(1)
  sx(5:7) = vec(2)
  sx(8:9) = vec(3)
  sx(10) = vec(4)
  sy(1) = vec(1)
  sy(2) = vec(2)
  sy(3) = vec(3)
  sy(4) = vec(4)
  sy(5) = vec(1)
  sy(6) = vec(2)
  sy(7) = vec(3)
  sy(8) = vec(1)
  sy(9) = vec(2)
  sy(10) = vec(1)
  sz(1) = vec(4)
  sz(2) = vec(3)
  sz(3) = vec(2)
  sz(4) = vec(1)
  sz(5) = vec(3)
  sz(6) = vec(2)
  sz(7) = vec(1)
  sz(8) = vec(2)
  sz(9) = vec(1)
  sz(10) = vec(1)
  angsol(1) = poids(2)
  angsol(2) = poids(1)
  angsol(3) = poids(1)
  angsol(4) = poids(2)
  angsol(5) = poids(1)
  angsol(6) = poids(3)
  angsol(7) = poids(1)
  angsol(8) = poids(1)
  angsol(9) = poids(1)
  angsol(10) = poids(2)

!==================================================================
! Quadrature Tn : 8n^2 directions
!==================================================================


elseif (i_quadrature.eq.4) then    ! Quadrature T2 : 32 directions

  vec(1) = 0.2357022604
  vec(2) = 0.9428090416
  vec(3) = 0.5773502692
  poids(1) = 0.5512855984
  poids(2) = 0.3398369095

  sx(1)   = vec(1)
  sx(2)   = vec(2)
  sx(3)   = vec(3)
  sx(4)   = vec(1)
  sy(1:2) = vec(1)
  sy(3)   = vec(3)
  sy(4)   = vec(2)
  sz(1)   = vec(2)
  sz(2)   = vec(1)
  sz(3)   = vec(3)
  sz(4)   = vec(1)
  angsol(1) = poids(2)
  angsol(2) = poids(2)
  angsol(3) = poids(1)
  angsol(4) = poids(2)

elseif (i_quadrature.eq.5) then   ! Quadrature T4 : 128 directions

  vec(1)  = 0.0990147543
  vec(2)  = 0.4923659639
  vec(3)  = 0.2357022604
  vec(4)  = 0.1230914910
  vec(5)  = 0.8616404369
  vec(6)  = 0.6804138174
  vec(7)  = 0.5773502692
  vec(8)  = 0.2721655270
  vec(9)  = 0.9901475430
  vec(10) = 0.9428090416

  poids(1) =0.0526559083
  poids(2) =0.0995720042
  poids(3) =0.0880369928
  poids(4) =0.1320249278
  poids(5) =0.1552108150

  sx(1)   = vec(1)
  sx(2)   = vec(2)
  sx(3)   = vec(3)
  sx(4)   = vec(4)
  sx(5)   = vec(5)
  sx(6)   = vec(6)
  sx(7)   = vec(7)
  sx(8)   = vec(8)
  sx(9)   = vec(4)
  sx(10)  = vec(9)
  sx(11)  = vec(10)
  sx(12)  = vec(5)
  sx(13)  = vec(6)
  sx(14)  = vec(2)
  sx(15)  = vec(3)
  sx(16)  = vec(1)
  sy(1)   = vec(1)
  sy(2)   = vec(4)
  sy(3)   = vec(3)
  sy(4)   = vec(2)
  sy(5)   = vec(4)
  sy(6)   = vec(8)
  sy(7)   = vec(7)
  sy(8)   = vec(6)
  sy(9)   = vec(5)
  sy(10)  = vec(1)
  sy(11)  = vec(3)
  sy(12)  = vec(2)
  sy(13)  = vec(6)
  sy(14)  = vec(5)
  sy(15)  = vec(10)
  sy(16)  = vec(9)
  sz(1)   = vec(9)
  sz(2)   = vec(5)
  sz(3)   = vec(10)
  sz(4)   = vec(5)
  sz(5)   = vec(2)
  sz(6)   = vec(6)
  sz(7)   = vec(7)
  sz(8)   = vec(6)
  sz(9)   = vec(2)
  sz(10)  = vec(1)
  sz(11)  = vec(3)
  sz(12)  = vec(4)
  sz(13)  = vec(8)
  sz(14)  = vec(4)
  sz(15)  = vec(3)
  sz(16)  = vec(1)

  angsol(1) = poids(1)
  angsol(2) = poids(2)
  angsol(3) = poids(3)
  angsol(4) = poids(2)
  angsol(5) = poids(2)
  angsol(6) = poids(4)
  angsol(7) = poids(5)
  angsol(8) = poids(4)
  angsol(9) = poids(2)
  angsol(10) = poids(1)
  angsol(11) = poids(3)
  angsol(12) = poids(2)
  angsol(13) = poids(4)
  angsol(14) = poids(2)
  angsol(15) = poids(3)
  angsol(16) = poids(1)

elseif (i_quadrature.eq.6) then   ! Quadrature Tn : 8 n^2 directions

  nquad = ndirs

  ! calcul de la position x et la position z du centre de tous les sous triangles sur le grand triangle 2D

  nbpoint = 2 * nquad - 1   ! nombre de point max dans un niveau

  allocate(x3d(nquad, nbpoint))
  allocate(y3d(nquad, nbpoint))
  allocate(z3d(nquad, nbpoint))

  z3d = 0.d0
  y3d = 0.d0
  x3d = 0.d0

  ! position z

  do ii = 0, nquad - 1
    niv = nquad - ii
    do jj = 1, 2 * niv - 1, 2
      z3d(niv, jj) = (1.d0 + 3.d0 * DBLE(ii)) / (3.d0 * DBLE(nquad))
    enddo
    do jj = 2, 2 * niv - 1, 2
      z3d(niv, jj) = (2.d0 + 3.d0 * DBLE(ii)) / (3.d0 * DBLE(nquad))
    enddo
  enddo

 ! position y

 ! la position y pour chaque point d'un niveau augmente de 1/(3*(nbpoint+1)/2) dans une meme colonne de triangle.
 ! on remplit donc les points en y par colonne (diagonale allant du haut vers bas gauche)

  jj = 0
  int1 = 0
  do tri = 1, nbpoint
    do niv = 1 + jj, nquad
      y3d(niv, tri) = DBLE(tri + jj) / (3.d0 * (DBLE(nbpoint) + 1.d0) / 2.d0)
    enddo
    int1 = int1 + 1
    if (int1.eq.2) then
      int1 = 0
      jj = jj + 1
    endif
  enddo

  ! position x

  x3d = 1.0d0 - z3d - y3d

  ! reecriture dans sx,sy,sz

  jj = 0
  do niv=1,nquad
    nbpoint = 2 * niv - 1
    do ii = 1, nbpoint
      jj = jj+1
      sx(jj) = x3d(niv, ii)
      sy(jj) = y3d(niv, ii)
      sz(jj) = z3d(niv, ii)
    enddo
  enddo

  ! normalisation des vecteurs

  do ii = 1, ndirs
    vect(1) = sx(ii)
    vect(2) = sy(ii)
    vect(3) = sz(ii)
    call normve(vect)
    sx(ii) = vect(1)
    sy(ii) = vect(2)
    sz(ii) = vect(3)
  enddo

  deallocate(x3d)
  deallocate(y3d)
  deallocate(z3d)

 ! on connait maintenant toutes les directions, on va maintenant calculer les poids de chacune des directions.

 ! calcul de la position des sous noeuds

 !                 1           ! niv 1
 !                 /\
 !                /  \
 !              1/____\2       ! niv 2
 !              /\    /\
 !             /  \  /  \
 !           1/____\/____\3    ! niv 3
 !                 2

      nbpoint = nquad + 1   ! nombre de point max dans un niveau

      allocate(x3d(nquad+1, nbpoint))
      allocate(y3d(nquad+1, nbpoint))
      allocate(z3d(nquad+1, nbpoint))

  ! position z
  niv = 1
  z3d(niv, 1) = 1.d0

  do ii = 2, nquad + 1
    niv = niv + 1
    do jj = 1, niv   ! il y a niv points dans le niveau niv
      z3d(niv, jj) = z3d(niv - 1, 1) - 1.d0 / DBLE(nquad)
    enddo
  enddo

  ! position y
  ! on remplit donc les points en y par colonne (diagonale allant du haut vers bas gauche)

  niv = 1
  y3d(niv, 1) = 0.d0

  do ii = 2,nquad+1
    niv = niv + 1
    do jj = 1, niv  ! Il y a niv points dans le niveau niv
      y3d(niv, jj) = DBLE(jj - 1) * (1.d0 - z3d(niv, 1)) / DBLE(niv - 1)
    enddo
  enddo

  ! position x

  x3d = 1.d0 - z3d - y3d

  ! calcul de l'aire et du poids de chaque triangle

  ntri = 0

  do niv = 1, nquad  ! nombre de niveau de triangle

    ! triangle n°1 du niveau
    allocate(posnod(3, 3))

    posnod(1, 1) = x3d(niv, 1)
    posnod(2, 1) = y3d(niv, 1)
    posnod(3, 1) = z3d(niv, 1)
    posnod(1, 2) = x3d(niv + 1, 1)
    posnod(2, 2) = y3d(niv + 1, 1)
    posnod(3, 2) = z3d(niv + 1, 1)
    posnod(1, 3) = x3d(niv + 1, 2)
    posnod(2, 3) = y3d(niv + 1, 2)
    posnod(3, 3) = z3d(niv + 1, 2)
    ntri = ntri + 1
    call lhuilier(posnod, angsol(ntri))

    deallocate(posnod)

    do ii = 1, niv - 1   ! nombre de double triangle dans le niveau

      allocate(posnod(3,3))

      posnod(1, 1) = x3d(niv, 1 + (ii - 1))
      posnod(2, 1) = y3d(niv, 1 + (ii - 1))
      posnod(3, 1) = z3d(niv, 1 + (ii - 1))
      posnod(1, 2) = x3d(niv, 2 + (ii - 1))
      posnod(2, 2) = y3d(niv, 2 + (ii - 1))
      posnod(3, 2) = z3d(niv, 2 + (ii - 1))
      posnod(1, 3) = x3d(niv + 1, 2 + (ii - 1))
      posnod(2, 3) = y3d(niv + 1, 2 + (ii - 1))
      posnod(3, 3) = z3d(niv + 1, 2 + (ii - 1))
      ntri = ntri + 1
      call lhuilier(posnod, angsol(ntri))

      deallocate(posnod)

      allocate(posnod(3,3))
      posnod(1,1) = x3d(niv, 2 + (ii - 1))
      posnod(2,1) = y3d(niv, 2 + (ii - 1))
      posnod(3,1) = z3d(niv, 2 + (ii - 1))
      posnod(1,2) = x3d(niv + 1, 2 + (ii - 1))
      posnod(2,2) = y3d(niv + 1, 2 + (ii - 1))
      posnod(3,2) = z3d(niv + 1, 2 + (ii - 1))
      posnod(1,3) = x3d(niv + 1, 3 + (ii - 1))
      posnod(2,3) = y3d(niv + 1, 3 + (ii - 1))
      posnod(3,3) = z3d(niv + 1, 3 + (ii - 1))
      ntri = ntri + 1

      call lhuilier(posnod, angsol(ntri))

      deallocate(posnod)

    enddo

  enddo

  deallocate(x3d)
  deallocate(y3d)
  deallocate(z3d)

endif

return

end subroutine


!===============================================================================

subroutine lhuilier (posnod, angleS)
!================

!===============================================================================
!  FONCTION  :
!  ---------

! compute the solid angle associated to a given direction of the Tn quadrature

!-------------------------------------------------------------------------------
! Arguments

double precision posnod(3,3)
double precision angleS

! Local variable
double precision a,b,c,p

!-------------------------------------------------------------------------------


! on renormalise pour les produits scalaire

call normve(posnod(:,1))
call normve(posnod(:,2))
call normve(posnod(:,3))

!longueur des segments du triangle curviligne (R=1)

a = dacos(dot_product(posnod(:, 1), posnod(:, 2)))
b = dacos(dot_product(posnod(:, 2), posnod(:, 3)))
c = dacos(dot_product(posnod(:, 3), posnod(:, 1)))

p = 0.5d0 * (a + b + c) ! périmètre

angleS = 4 * datan(dsqrt(tan(p / 2) * tan((p - a) / 2) * tan((p - b) / 2)*tan((p - c) / 2)))

return

end subroutine

!===============================================================================

subroutine normve (vecteur)
!================

!===============================================================================
!  FONCTION  :
!  ---------

! vector normalization

!-------------------------------------------------------------------------------
! Argument

double precision vecteur(3)

! Local variable

double precision norm

!-------------------------------------------------------------------------------

norm = dsqrt(vecteur(1)**2+vecteur(2)**2+vecteur(3)**2)

vecteur = vecteur / norm

return

end subroutine
