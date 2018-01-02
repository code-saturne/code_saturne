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

!===============================================================================
! Function:
! ---------

!> \file distyp.f90
!>
!> \brief This subroutine computes the dimensionless distance to the wall
!> solving a transport equation.
!>
!> This function solves the following transport equation on \f$ \varia \f$:
!> \f[
!> \dfrac{\partial \varia}{\partial t} + \divs \left( \varia \vect{V} \right)
!>     - \divs \left( \vect{V} \right) \varia = 0
!> \f]
!> where the vector field \f$ \vect{V} \f$ is defined by:
!> \f[
!>  \vect{V} = \dfrac{ \grad y }{\norm{\grad y} }
!> \f]
!> The boundary conditions on \f$ \varia \f$ read:
!> \f[
!>  \varia = \dfrac{u_\star}{\nu} \textrm{ on walls}
!> \f]
!> \f[
!>  \dfrac{\partial \varia}{\partial n} = 0 \textrm{ elsewhere}
!> \f]
!>
!> Then the dimensionless distance is deduced by:
!> \f[
!>  y^+ = y \varia
!> \f]
!>
!>
!> Remarks:
!> - a steady state is looked for.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     itypfb        boundary face types
!> \param[out]    disty         dimensionless distance \f$ y^+ \f$
!_______________________________________________________________________________

subroutine distyp &
 ( itypfb ,                                                       &
   disty  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use coincl
use parall
use period
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)

double precision disty(ncelet)

! Local variables

integer          idtva0, f_id0, f_id  , iconvp, idiffp
integer          ndircp
integer          iescap, iflmb0, itypfl
integer          ifac  , iel   , init
integer          inc   , iccocg, isym  , isweep, infpar
integer          imucpp, idftnp, iswdyp, icvflb
integer          isou  , jsou

integer          ivoid(1)

double precision xnorme, dtminy, dtmaxy, relaxp, thetap, timey
double precision xusnmx, xusnmn, xnorm0
double precision dismax, dismin, usna

double precision rvoid(1)

double precision, allocatable, dimension(:) :: dvarp, smbdp, rovsdp
double precision, allocatable, dimension(:,:) :: q
double precision, allocatable, dimension(:) :: flumas, flumab
double precision, allocatable, dimension(:) :: rom, romb
double precision, allocatable, dimension(:) :: coefap, coefbp
double precision, allocatable, dimension(:,:) :: coefav
double precision, allocatable, dimension(:,:,:) :: coefbv
double precision, allocatable, dimension(:) :: w1, w2
double precision, allocatable, dimension(:) :: dpvar
double precision, dimension(:), pointer :: w_dist
double precision, dimension(:), pointer :: crom, uetbor
double precision, dimension(:), pointer :: viscl

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays for the distance resolution
allocate(dvarp(ncelet), smbdp(ncelet), rovsdp(ncelet))
allocate(q(3,ncelet))
allocate(flumas(nfac), flumab(nfabor))
allocate(rom(ncelet), romb(nfabor))
allocate(coefap(nfabor), coefbp(nfabor))
allocate(coefav(3,nfabor))
allocate(coefbv(3,3,nfabor))
allocate(dpvar(ncelet))

! Allocate work arrays
allocate(w1(ncelet))
allocate(w2(ncelet))

ipass  = ipass + 1

call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)

uetbor => null()

call field_get_id_try('ustar', f_id)
if (f_id.ge.0) then
  call field_get_val_s(f_id, uetbor)
endif

call field_get_id("wall_distance", f_id)
call field_get_val_s(f_id, w_dist)

!===============================================================================
! 2. At the first time step
!===============================================================================

! Au premier pas de temps, on a en general u* = 0 (ou faux)
!   on ne calcule pas y+

! En effet ca prend du temps, d'autant plus que u* est petit, car il
!   alors calculer y+ jusqu'a une grande distance des parois

if(ntcabs.eq.1) then

  do iel = 1, ncel
    disty(iel) = grand
  enddo

  if(iwarny.ge.1) then
    write(nfecra,7000)
  endif

  return

endif

!===============================================================================
! 3. Compute  V = Grad(DISTPA)/|Grad(DISTPA)|
!===============================================================================

! Compute the gradient of the distance to the wall

inc    = 1
iccocg = 1

call field_get_id("wall_distance", f_id)

! Current gradient: iprev = 0
call field_gradient_scalar(f_id, 0, imrgra, inc, iccocg, q)

! Normalization (warning, the gradient may be sometimes equal to 0)
do iel = 1, ncel
  xnorme = max(sqrt(q(1,iel)**2+q(2,iel)**2+q(3,iel)**2),epzero)
  do isou = 1, 3
    q(isou,iel) = q(isou,iel)/xnorme
  enddo
enddo

! Paralellism and periodicity
if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(q)
endif

!===============================================================================
! 4. Compute the flux of V
!===============================================================================

do ifac = 1, nfabor
  romb(ifac) = 1.d0
enddo
do iel = 1, ncelet
  rom(iel)  = 1.d0
enddo

! Le gradient normal de la distance a la paroi vaut -1 en paroi
!   par definition et obeit a un flux nul ailleurs

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
    xnorme = max(surfbn(ifac),epzero**2)
    do isou = 1, 3
      coefav(isou,ifac) = -surfbo(isou,ifac)/xnorme
      do jsou = 1, 3
        coefbv(isou,jsou,ifac) = 0.d0
      enddo
    enddo
  else
    do isou = 1, 3
      coefav(isou,ifac) = 0.d0
      do jsou = 1, 3
        if (isou == jsou) then
          coefbv(isou,jsou,ifac) = 1.d0
        else
          coefbv(isou,jsou,ifac) = 0.d0
        endif
      enddo
    enddo
  endif
enddo

! Calcul du flux de masse

! On ne le met pas a zero en paroi (justement)
iflmb0 = 0
! On l'initialise a 0
init   = 1
! On prend en compte les Dirichlet
inc    = 1
! Il ne s'agit ni de U ni de R
f_id   = -1

itypfl = 1

call inimav                                                       &
 ( f_id   , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgy , imligy ,          &
   iwarny ,                                                       &
   epsrgy , climgy ,                                              &
   rom    , romb   ,                                              &
   q      ,                                                       &
   coefav , coefbv ,                                              &
   flumas , flumab )

!===============================================================================
! 5. Boundary conditions
!===============================================================================

! Dirichlet en u*/nu aux parois, et flux nul ailleurs

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then
    iel = ifabor(ifac)
    coefap(ifac) = uetbor(ifac)*crom(iel)/viscl(iel)
    coefbp(ifac) = 0.0d0
  else
    coefap(ifac) = 0.0d0
    coefbp(ifac) = 1.0d0
  endif
enddo

!===============================================================================
! 6. Compute the time step
!===============================================================================

! On vise un Courant infini (de l'ordre de 1000).

! On calcule avec MATRDT DA = Sigma a S/d
iconvp = 1
idiffp = 0
!     La matrice est non symetrique
isym   = 2

! Warning: no diffusion here, so no need of other Boundary coefficient

call matrdt &
!==========
 ( iconvp , idiffp , isym   ,                                     &
   coefbp , coefbp , flumas , flumab , flumas , flumab , w2     )

! Le Courant est COUMXY = DT w2 / VOLUME
!     d'ou DTMINY = MIN(COUMXY * VOLUME/w2)
! Au cas ou une cellule serait a w2(IEL)=0,
!   on ne la prend pas en compte

! On prend dans QZ un pas de temps variable,
!   si on ne trouve pas de pas de temps, on prend le minimum
!   (cellules sources)
dtminy = grand
dtmaxy = -grand
do iel = 1, ncel
  w1(iel) = -grand
  if(w2(iel).gt.epzero) then
    w1(iel) = coumxy*volume(iel)/w2(iel)
    dtminy  = min(w1(iel),dtminy)
    dtmaxy  = max(w1(iel),dtmaxy)
  endif
enddo
if(irangp.ge.0) then
  call parmin (dtminy)
  call parmax (dtmaxy)
endif
dtminy = max(dtminy,epzero)

do iel = 1, ncel
  if(w1(iel).le.0.d0) then
    w1(iel) = dtminy
  endif
enddo

if(iwarny.ge.2) then
  write(nfecra,2000)dtminy,dtmaxy
endif

!===============================================================================
! 7. Diagonal part of the matrix
!===============================================================================

do iel = 1, ncel
  rovsdp(iel) = volume(iel)*rom(iel)/w1(iel)
enddo

!===============================================================================
! 8. Time loop
!===============================================================================

! Initializations
!=================

! Iterations
isweep = 0

! Temps
timey = 0.d0

! Inconnue
!   Au cas ou on n'atteint pas tout a fait l'etat stationnaire,
!   il faut que le yplus ne soit pas nul dans la zone ou les
!   conditions aux limites n'ont pas ete convectees. On voudrait
!   plutot que yplus y soit maximum.
!   Si on utilise zero ou une valeur negative comme initialisation,
!   on risque de se retrouver avec des valeurs proches de
!   zero issues de la diffusion due au schema upwind au voisinage
!   du front convecte et donc avec des yplus proches de zero
!   n'importe ou.
!   On va donc utiliser la valeur max de u*/nu.

!   A partir du second pas de temps, on a egalement le yplus du pas
!     de temps precedent

! On calcule le min et le max
xusnmx = -grand
xusnmn =  grand
do ifac = 1, nfabor
  if(itypfb(ifac).eq.iparoi .or.                            &
     itypfb(ifac).eq.iparug) then
    xusnmx = max(xusnmx,coefap(ifac))
    xusnmn = min(xusnmn,coefap(ifac))
  endif
enddo
if(irangp.ge.0) then
  call parmax (xusnmx)
  call parmin (xusnmn)
endif

if(ipass.eq.1) then
  do iel = 1, ncelet
    dvarp(iel) = xusnmx
  enddo
else
  do iel = 1, ncel
    usna = disty(iel)/max(w_dist(iel),epzero)
    usna = max(usna,xusnmn)
    usna = min(usna,xusnmx)
    dvarp(iel) = usna
  enddo
endif

! Norme de reference (moyenne des u*/nu)
!   (on divise par le nombre de faces)
xnorm0 = 0.d0
infpar = 0
do ifac = 1, nfabor
  if(itypfb(ifac).eq.iparoi .or.                            &
     itypfb(ifac).eq.iparug) then
    infpar = infpar+1
    xnorm0 = xnorm0 + coefap(ifac)**2
  endif
enddo
if(irangp.ge.0) then
  call parcpt (infpar)
  call parsom (xnorm0)
endif
xnorm0 = xnorm0/dble(infpar)

! To prevent division by 0
if (xnorm0.le.epzero**2) goto 100


! Loops beginning
!=================

do isweep = 1, ntcmxy

  ! Instant (arbitrairement +DTMINY)
  timey = timey + dtminy

  ! -- Echange pour les cas paralleles
  !     a la premiere iteration, c'est inutile (on a fait l'init sur NCELET)

  if(isweep.gt.1.or.ipass.gt.1) then

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(dvarp)
      !==========
    endif

  endif


  ! Save of the solution for convergence test
  !===========================================

  do iel = 1, ncel
    w1(iel) = dvarp(iel)
  enddo

  ! Right hand side
  !=================
  !   Obligatoirement a tous les pas de temps

  do iel = 1, ncel
    smbdp(iel) = 0.d0
  enddo

  ! Solving
  !=========

  ! La variable n'est pas la vitesse ou une composante de Rij
  f_id0= -1
  ! Le cas est convectif, non diffusif
  iconvp = 1
  idiffp = 0
  ! Il y a des Dirichlet (car il y a des parois)
  ndircp = 1
  ! Pas d'estimateurs, ni de multigrille (100 et 10 sont arbitraires)
  iescap = 0
  imucpp = 0
  idftnp = 1
  iswdyp = 0
  nomva0 = 'yplus_wall'
  ! Ordre 1 en temps (etat stationnaire cherche)
  thetap = 1.d0
  ! Pas de stationnaire ni de relaxation -> a modifier eventuellement
  idtva0 = 0
  relaxp = 1.d0
  ! all boundary convective flux with upwind
  icvflb = 0
  ! Warning: no diffusion so no need of other diffusive Boundary coeeficient

  call codits &
  !==========
 ( idtva0 , f_id0  , iconvp , idiffp , ndircp ,                   &
   imrgra , nswrsy , nswrgy , imligy , ircfly ,                   &
   ischcy , isstpy , iescap , imucpp , idftnp , iswdyp ,          &
   iwarny ,                                                       &
   blency , epsily , epsrsy , epsrgy , climgy , extray ,          &
   relaxp , thetap ,                                              &
   dvarp  , dvarp  ,                                              &
   coefap , coefbp ,                                              &
   coefap , coefbp ,                                              &
   flumas , flumab ,                                              &
   flumas , flumab , flumas , flumab , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdp , smbdp  , dvarp  , dpvar  ,                            &
   rvoid  , rvoid  )


  ! Clipping (indispensable si on initialise par u*/nu du pas de
  !==========                                   temps precedent)

  do iel = 1, ncel
    dvarp(iel) = max(dvarp(iel),xusnmn)
    dvarp(iel) = min(dvarp(iel),xusnmx)
  enddo

  ! Stopping test
  !===============

  ! on utilise QY dans lequel la solution precedente a ete sauvee
  ! on souhaite que la variation de l'inconnue sur chaque cellule
  !   soit inferieure a un pourcentage de la valeur moyenne de
  !   u*/nu calculee sur les faces de bord
  ! on limite le test aux cellules qui pourront avoir un interet pour
  !   VanDriest, c'est a dire qu'on ignore celles a y+ > YPLMXY
  !   comme on ne connait pas encore y+, on se base sur min(u*/nu) :
  !   on ignore les cellules a y min(u*/nu) > YPLMXY

  xnorme = -grand
  do iel = 1, ncel
    if(w_dist(iel)*xusnmn.le.yplmxy) then
      xnorme = max(xnorme,(dvarp(iel)-w1(iel))**2)
    endif
  enddo
  if (irangp.ge.0) then
    call parmax (xnorme)
  endif

  if (iwarny.ge.2) then
    write(nfecra,3000)isweep,xnorme,xnorm0,xnorme/xnorm0
  endif

  if (xnorme.le.epscvy*xnorm0) goto 100

enddo

write(nfecra,8000) xnorme, xnorm0, xnorme/xnorm0, ntcmxy

 100  continue


!===============================================================================
! 9. Finalization and printing
!===============================================================================

do iel = 1, ncel
  disty(iel) = dvarp(iel)*w_dist(iel)
enddo

dismax = -grand
dismin =  grand

do iel = 1, ncel
  dismin = min(disty(iel),dismin)
  dismax = max(disty(iel),dismax)
enddo

if (irangp.ge.0) then
  call parmin(dismin)
  call parmax(dismax)
endif

if (iwarny.ge.1) then
  write(nfecra,1000) dismin, dismax, min(isweep,ntcmxy)
endif

! Free memory
deallocate(dvarp, smbdp, rovsdp)
deallocate(q)
deallocate(flumas, flumab)
deallocate(rom, romb)
deallocate(coefap, coefbp)
deallocate(coefav, coefbv)
deallocate(w1, w2)
deallocate(dpvar)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format( &
'                                                             ',/,&
' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,&
'    ------------------------------------                     ',/,&
'                                                             ',/,&
'  Distance+ min = ',E14.5    ,' Distance+ max = ',E14.5       ,/,&
'                                                             ',/,&
'     (Calcul de la distance realise en ',I10   ,' iterations)',/)
 2000 format( &
'                                                             ',/,&
' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,&
'    ------------------------------------                     ',/,&
'                                                             ',/,&
' Yplus:  Dt min = ',E14.5    ,'        Dt max = ',E14.5       ,/)
 3000 format( &
'                                                             ',/,&
' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,&
'    ------------------------------------                     ',/,&
'                                                             ',/,&
' Yplus:  iteration   residu abs.     reference   residu rel. ',/,&
' Yplus: ',I10    ,E14.5        ,E14.5        ,E14.5           ,/)
 7000 format( &
'                                                             ',/,&
' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,&
'    ------------------------------------                     ',/,&
'                                                             ',/,&
'  Elle n''est pas calculee au premier pas de temps           ',/)
 8000 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : Calcul de la distance a la paroi            ',/,&
'@    =========                                               ',/,&
'@     Le systeme iteratif pour le calcul de la distance a la ',/,&
'@       paroi adimensionnelle est imparfaitement converge.   ',/,&
'@                                                            ',/,&
'@          Residu     Reference            Residu relatif    ',/,&
'@  ',2E14.5,12X,E14.5                                         ,/,&
'@                                                            ',/,&
'@     Augmenter la valeur de NTCMXY dans usipsu peut resoudre',/,&
'@       le probleme.                                         ',/,&
'@     La valeur actuelle de cet entier est ',I10              ,/,&
'@                                                            ',/,&
'@     Le calcul se poursuit.                                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format( &
'                                                             ',/,&
' ** DIMENSIONLESS WALL DISTANCE                              ',/,&
'    ---------------------------                              ',/,&
'                                                             ',/,&
'  Min distance+ = ',E14.5    ,' Max distance+ = ',E14.5       ,/,&
'                                                             ',/,&
'     (Distance calculation done in ',I10   ,' iterations)'    ,/)
 2000 format( &
'                                                             ',/,&
' ** DIMENSIONLESS WALL DISTANCE                              ',/,&
'    ---------------------------                              ',/,&
'                                                             ',/,&
' Yplus:  Min dt = ',E14.5    ,'        Max dt = ',E14.5       ,/)
 3000 format( &
'                                                             ',/,&
' ** DIMENSIONLESS WALL DISTANCE                              ',/,&
'    ---------------------------                              ',/,&
'                                                             ',/,&
' Yplus:  iteration   abs. residu     reference   rel. residu ',/,&
' Yplus: ',I10    ,E14.5        ,E14.5        ,E14.5           ,/)
 7000 format( &
'                                                             ',/,&
' ** DIMENSIONLESS WALL DISTANCE                              ',/,&
'    ---------------------------                              ',/,&
'                                                             ',/,&
'  It is not computed at the first time step                  ',/)
 8000 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Dimensionless distance to the wall computation ',/,&
'@    ========                                                ',/,&
'@     The iterative system for the dimensionless distance to ',/,&
'@       the wall calculation has not properly converged.     ',/,&
'@                                                            ',/,&
'@          Residual   Reference            Relative residual ',/,&
'@  ',2E14.5,12X,E14.5                                         ,/,&
'@                                                            ',/,&
'@     Increase the value of NTCMXY in usipsu may solve       ',/,&
'@       the problem.                                         ',/,&
'@     The current value of this integer is ',I10              ,/,&
'@                                                            ',/,&
'@     The calculation will be run.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! End
!----

return
end subroutine
