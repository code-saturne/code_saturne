!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

!> \file sootsc.f90
!>
!> \brief Routine physique particuliere : modele de suies a 2 equations.
!>
!> On precise les termes sources pour les scalaires
!> fraction massique de suies et nombre de précurseurs
!> sur un pas de temps
!>
!> Attention :
!>  le traitement des termes sources est different
!>  de celui de ustssc.f
!>
!> On resout rovsdt*d(var) = smbrs
!>
!> rovsdt et smbrs contiennent deja d'eventuels termes sources
!>  utilisateur. il faut donc les incrementer et pas les
!>  ecraser
!>
!> Pour des questions de stabilite, on ne rajoute dans rovsdt
!>  que des termes positifs. il n'y a pas de contrainte pour
!>  smbrs
!>
!> Dans le cas d'un terme source en cexp + cimp*var on doit
!> ecrire :
!>          smbrs  = smbrs  + cexp + cimp*var
!>          rovsdt = rovsdt + max(-cimp,zero)
!>
!> On fournit ici rovsdt et smbrs (ils contiennent RHO*volume)
!>  - smbrs en kg variable/s :
!>     ex : pour la vitesse            kg m/s2
!>          pour les temperatures      kg degres/s
!>          pour les enthalpies        Joules/s
!>  - rovsdt en kg /s
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     iscal         scalar index
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        type of mass source term for the variables
!> \param[in]     izfppp        boundary zone index
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfa        physical properties at interior face centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[in]     coefa, coefb  boundary conditions
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!> \param[in,out] smbrs         explicit right hand side
!> \param[in,out] rovsdt        implicit terms
!_______________________________________________________________________________

subroutine sootsc &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt )

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
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

character*80     chaine
integer          ivar, ipcrom, iel

double precision epsi
parameter       (epsi = 1.d-6)
double precision d1s3, d2s3, cexp, cimp
double precision zetan, zetas, rho, xfu, xm, temp, nn0
double precision ka, kb, kz, kt, chi, po2, wox
double precision aa, bb, cc, taa, tcc, caa, cbb, ccc, dd

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

ivar = isca(iscal)
chaine = nomvar(ipprtp(ivar))
ipcrom = ipproc(irom)

!===============================================================================
! 2. Writtings
!===============================================================================

if (iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!=======================================================================
! --- Moss et al.:
! zeta_s (isca(ifsm)) soot mass fraction zeta_s = (rho_s/rho).f_v
! zeta_n (isca(inpm)) precursor density  zeta_n = n / (rho.No)
!=======================================================================

if (ivar.eq.isca(ifsm).or.ivar.eq.isca(inpm)) then

  ! To be changed for other combustible !FIXME
  ! Methane CH4 (Syed, Stewart and Moss Symposium 1990)
  caa = 6.54d4 !m^3/kg^2.K^0.5.s
  cbb = 1.3d7 ! m^3.K^-1/2.s^-1
  ccc = 0.1d0 ! m^3.kg^-2/3.K^-1/2.s^-1
  taa = 46.1d3 ! K
  tcc = 12.6d3 ! K

  d1s3 = 1.d0/3.d0
  d2s3 = 2.d0/3.d0

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(propce(1,ipproc(itemp)))
    call synsca(rtp(1,ivar))
  endif

  do iel = 1, ncel

    cexp = 0.d0
    cimp = 0.d0

    nn0 = 6.0223d23
    rho = propce(iel,ipproc(irom)) ! Mixture density (kg/m3)
    temp = propce(iel,ipproc(itemp)) ! Temperature

    xm = 1.d0/ (  propce(iel,ipproc(iym(1)))/wmolg(1)             &
                + propce(iel,ipproc(iym(2)))/wmolg(2)             &
                + propce(iel,ipproc(iym(3)))/wmolg(3) )

    xfu = propce(iel,ipproc(iym(1))) * xm / wmolg(1) ! Fuel molar fraction

    ! --- rate of particule nucleation
    aa = caa * rho**2 * temp**0.5d0 * xfu * exp(-taa/temp)

    ! --- coagulation
    bb = cbb * temp**0.5d0

    ! --- surface growth of soot
    cc = ccc * rho * temp**0.5d0 * xfu * exp(-tcc/temp)

    po2 = propce(iel,ipproc(iym(2)))*xm/wmolg(2)*1.d0/4.76d0

    ! --- oxidation
    ka = 20.d0*exp(-15098.d0/temp)
    kb = 4.46d-3*exp(-7650.d0/temp)
    kt = 1.51d5*exp(-48817.d0/temp)
    kz = 21.3d0*exp(2063.d0/temp)

    chi = kb*po2/(kb*po2+kt)

    wox = 1.2d2*( (ka*po2*chi)/(1.d0+kz*po2) + kb*po2*(1.d0-chi) )

    dd = (36.d0*acos(-1.d0)/rosoot**2.d0)**d1s3

   ! -------------------------------------------------------------

    zetas = rtpa(iel,isca(ifsm)) ! fraction massique de suies (SU)
    zetan = rtpa(iel,isca(inpm)) ! densite de precurseurs (SU)

    if (ivar.eq.isca(ifsm)) then

      ! --- Surface growth : quadratic
      if (zetas.gt.epsi) cimp = volume(iel) *                     &
        (  nn0**d1s3 * rho * cc * zetas**(-d1s3) * zetan**d1s3    &
         - rho * dd *nn0**d1s3 *zetan**d1s3 *zetas**(-d1s3)*wox )
      cexp = volume(iel) * (  144.d0*aa )
    endif

    if (ivar.eq.isca(inpm)) then
      cimp = volume(iel) * ( - rho**2.d0 * bb * zetan )
      cexp = volume(iel) * ( aa )
    endif

    smbrs(iel)  = smbrs(iel)  + cexp + cimp*rtpa(iel,ivar)
    rovsdt(iel) = rovsdt(iel) + max(-cimp,0.d0)

  enddo

endif

!--------
! Formats
!--------

 1000 format(' TERMES SOURCES PHYSIQUE PARTICULIERE POUR LA VARIABLE '  &
       ,a8,/)

!----
! End
!----

return

end subroutine
