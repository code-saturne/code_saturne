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
!>
!> Then, Imposition of an amortization of Van Driest type for the LES.
!>        \f$ \nu_T \f$ is absorbed by \f$ (1-\exp(\dfrac{-y^+}{d^+}))^2 \f$
!>        where \f$ d^+ \f$ is set at 26.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     itypfb        boundary face types
!> \param[in]     visvdr        dynamic viscosity in edge cells after
!>                               driest velocity amortization
!_______________________________________________________________________________

subroutine distyp &
 ( itypfb , visvdr)

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
double precision visvdr(ncelet)

! Local variables

integer          idtva0, f_id0, f_id  , iconvp, idiffp
integer          f_id_yplus
integer          ndircp
integer          iflmb0, itypfl
integer          ifac  , iel   , init
integer          inc   , iccocg, isym  , isweep
integer          imucpp, idftnp
integer          nswrgp, nswrsp
integer          icvflb, iescap, imligp, ircflp, iswdyp, isstpp, ischcp, iwarnp
integer          isou  , jsou
integer          iflmas, iflmab

integer          infpar
save             infpar

integer          ivoid(1)

double precision relaxp, blencp, climgp, epsilp, epsrgp, epsrsp, extrap
double precision thetap
double precision xnorme, dtminy, dtmaxy
double precision xusnmx, xusnmn, xnorm0
double precision dismax, dismin, usna
double precision hint, pimp, qimp

double precision rvoid(1)

double precision, allocatable, dimension(:) :: dvarp, smbdp, rovsdp
double precision, allocatable, dimension(:,:) :: q
double precision, allocatable, dimension(:,:) :: coefav
double precision, allocatable, dimension(:,:,:) :: coefbv
double precision, allocatable, dimension(:) :: w1, w2
double precision, allocatable, dimension(:) :: dpvar
double precision, dimension(:), pointer :: w_dist
double precision, pointer, dimension(:)   :: cvar_var
double precision, pointer, dimension(:)   :: cvara_var
double precision, dimension(:), pointer :: crom, uetbor
double precision, dimension(:), pointer :: viscl
double precision, dimension(:), pointer :: visct
double precision, pointer, dimension(:) :: coefap, coefbp
double precision, pointer, dimension(:) :: cofafp, cofbfp
double precision, dimension(:), pointer :: flumas, flumab
type(var_cal_opt) :: vcopt

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

call field_get_id("wall_yplus", f_id_yplus)
call field_get_key_struct_var_cal_opt(f_id_yplus, vcopt)

call field_get_val_s(f_id_yplus, cvar_var)
call field_get_val_prev_s(f_id_yplus, cvara_var)

call field_get_coefa_s( f_id_yplus, coefap)
call field_get_coefb_s( f_id_yplus, coefbp)
call field_get_coefaf_s(f_id_yplus, cofafp)
call field_get_coefbf_s(f_id_yplus, cofbfp)

call field_get_key_int(f_id_yplus, kimasf, iflmas)
call field_get_key_int(f_id_yplus, kbmasf, iflmab)

! Get pointer to the convective mass flux
call field_get_val_s(iflmas, flumas)
call field_get_val_s(iflmab, flumab)

! Number of wall faces
if (ipass.eq.1) then
  infpar = 0
  do ifac = 1, nfabor
    if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
      infpar = infpar+1
    endif
  enddo
  if (irangp.ge.0) then
    call parcpt(infpar)
  endif
endif

! If no wall, no wall distance
if (infpar.eq.0) then
  do iel = 1, ncelet
    cvar_var(iel) = grand
  enddo

  return
endif

!===============================================================================
! 2. At the first time step
!===============================================================================

! Au premier pas de temps, on a en general u* = 0 (ou faux)
!   on ne calcule pas y+

! En effet ca prend du temps, d'autant plus que u* est petit, car il
!   alors calculer y+ jusqu'a une grande distance des parois

if (ntcabs.eq.1) then

  do iel = 1, ncel
    cvar_var(iel) = grand
  enddo

  if (vcopt%iwarni.ge.1) then
    write(nfecra,7000)
  endif

  return

endif

!===============================================================================
! 3. Compute  V = Grad(y)/|Grad(y)|
!===============================================================================

! Compute the gradient of the distance to the wall

inc    = 1
iccocg = 1

call field_get_id("wall_distance", f_id)

! Current gradient: iprev = 0
call field_gradient_scalar(f_id, 0, imrgra, inc, iccocg, q)

! Normalization (warning, the gradient may be sometimes equal to 0)
do iel = 1, ncel
  xnorme = max(sqrt(q(1,iel)**2+q(2,iel)**2+q(3,iel)**2), epzero)
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

! Le gradient normal de la distance a la paroi vaut -1 en paroi
!   par definition et obeit a un flux nul ailleurs

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
    xnorme = max(surfbn(ifac), epzero**2)
    do isou = 1, 3
      coefav(isou,ifac) = - surfbo(isou,ifac)/xnorme
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

! Compute convective mass flux

! Mass flux non-zero at walls
iflmb0 = 0
! Default initilization at 0
init   = 1
! Take Dirichlet into account
inc    = 1
! q=grad(y) is not the Reynolds stress tensor
f_id   = -1

! Convective velocity NOT multiplied by rho
itypfl = 0

epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
iwarnp = vcopt%iwarni

call inimav                                                       &
 ( f_id   , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   rvoid  , rvoid  ,                                              &
   q      ,                                                       &
   coefav , coefbv ,                                              &
   flumas , flumab )

!===============================================================================
! 5. Boundary conditions
!===============================================================================

! Dirichlet u*/nu at walls, homogeneous Neumann elsewhere

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then
    iel = ifabor(ifac)

    ! Dirichlet Boundary Condition
    !-----------------------------

    hint = 1.d0/distb(ifac)
    pimp = uetbor(ifac)*crom(iel)/viscl(iel)

    call set_dirichlet_scalar &
         !====================
       ( coefap(ifac), cofafp(ifac),             &
         coefbp(ifac), cofbfp(ifac),             &
         pimp        , hint        , rinfin )


  else
    ! Neumann Boundary Conditions
    !----------------------------

    hint = 1.d0/distb(ifac)
    qimp = 0.d0

    call set_neumann_scalar &
         !==================
       ( coefap(ifac), cofafp(ifac),             &
         coefbp(ifac), cofbfp(ifac),             &
         qimp        , hint )

  endif
enddo

!===============================================================================
! 6. Compute the time step
!===============================================================================

! A large Courant number is wanted (of the order of 1000).

! On calcule avec MATRDT DA = Sigma a S/d
iconvp = 1
idiffp = 0
! Non symmetric matrix
isym   = 2

! Warning: no diffusion here, so no need of other Boundary coefficient

call matrdt &
!==========
 ( iconvp , idiffp , isym   ,                                     &
   coefbp , coefbp , flumas , flumab , flumas , flumab , w2     )

! Le Courant est coumxy = DT w2 / VOLUME
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

if (vcopt%iwarni.ge.2) then
  write(nfecra,2000) dtminy, dtmaxy
endif

!===============================================================================
! 7. Diagonal part of the matrix
!===============================================================================

do iel = 1, ncel
  rovsdp(iel) = volume(iel)/w1(iel)
enddo

!===============================================================================
! 8. Time loop
!===============================================================================

! Initializations
!=================

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
    usna = cvar_var(iel)/max(w_dist(iel),epzero)
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

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(dvarp)
  endif



  ! Save of the solution for convergence test
  !===========================================

  do iel = 1, ncel
    w1(iel) = dvarp(iel)
  enddo

  ! Right hand side
  !=================
  do iel = 1, ncel
    smbdp(iel) = 0.d0
  enddo

  ! Solving
  !=========
  iconvp = vcopt%iconv
  idiffp = vcopt%idiff
  idftnp = vcopt%idften
  nswrsp = vcopt%nswrsm
  nswrgp = vcopt%nswrgr
  imligp = vcopt%imligr
  ircflp = vcopt%ircflu
  ischcp = vcopt%ischcv
  isstpp = vcopt%isstpc
  ! No error estimate
  iescap = 0
  imucpp = 0
  iswdyp = vcopt%iswdyn
  iwarnp = vcopt%iwarni
  blencp = vcopt%blencv
  epsilp = vcopt%epsilo
  epsrsp = vcopt%epsrsm
  epsrgp = vcopt%epsrgr
  climgp = vcopt%climgr
  extrap = vcopt%extrag
  relaxp = vcopt%relaxv
  thetap = vcopt%thetav
  ! all boundary convective flux with upwind
  icvflb = 0
  init   = 1

  ! There are som Dirichlet BCs
  ndircp = 1
  ! No steady state algo
  idtva0 = 0

  ! Warning: no diffusion so no need of other diffusive Boundary coefficient

  call codits &
  !==========
 ( idtva0 , isweep , f_id_yplus, iconvp , idiffp , ndircp ,       &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   dvarp  , dvarp  ,                                              &
   coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
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

  if (vcopt%iwarni.ge.2) then
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
  cvar_var(iel) = dvarp(iel)*w_dist(iel)
enddo

dismax = -grand
dismin =  grand

do iel = 1, ncel
  dismin = min(cvar_var(iel),dismin)
  dismax = max(cvar_var(iel),dismax)
enddo

if (irangp.ge.0) then
  call parmin(dismin)
  call parmax(dismax)
endif

if (vcopt%iwarni.ge.1) then
  write(nfecra,1000) dismin, dismax, min(isweep,ntcmxy)
endif

!===============================================================================
! 10. Van Driest amortization
!===============================================================================

call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

do iel = 1, ncel
  visct(iel) = visct(iel)*(1.0d0-exp(-cvar_var(iel)/cdries))**2
enddo

! For the wall cells we add the turbulent viscosity which was absorbed
! in clptur and which has served to calculate the boundary conditions
do iel = 1, ncel
  if (visvdr(iel).gt.-900.d0) visct(iel) = visvdr(iel)
enddo

! Free memory
deallocate(dvarp, smbdp, rovsdp)
deallocate(q)
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
