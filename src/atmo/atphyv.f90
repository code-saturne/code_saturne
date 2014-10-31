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

subroutine atphyv &
!================

   ( propce )

!===============================================================================
! FONCTION :
! --------

! REMPLISSAGE DES VARIABLES PHYSIQUES : Atmospheric Version


! ATTENTION :
! =========

! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)

!  Il FAUT AVOIR PRECISE ICP = 1
!     ==================
!    dans usipsu si on souhaite imposer une chaleur specifique
!    CP variable (sinon: ecrasement memoire).


! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usipsu :
!             . la masse volumique (initialisee a RO0)
!             . la viscosite       (initialisee a VISCL0)
!      - dans usiniv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a la valeur donnee dans usiniv)

! On peut donner ici les lois de variation aux cellules
!     - de la masse volumique                      ROM    kg/m3
!         (et eventuellememt aux faces de bord     ROMB   kg/m3)
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)


! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)


! Il est conseille de ne garder dans ce sous programme que
!    le strict necessaire.



! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________!____!_____!________________________________________________!

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
use cstnum, only: pi
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use atincl
use field
use spefun
use field_operator

!===============================================================================

implicit none

! Arguments

double precision propce(ncelet,*)

! Local variables

integer          ivart, iel
integer          ipctem, ipcliq

double precision xvart, rhum, rscp, pp, zent
double precision lrhum, lrscp
double precision qsl, deltaq
double precision qliq, qwt, tliq, dum

double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:), pointer :: cvar_vart, cvar_totwt

logical activate

! External function

double precision qsatliq
external qsatliq
! call as: qsatliq(temperature,pressure)

!===============================================================================
! 0. INITIALISATIONS A CONSERVER
!===============================================================================

activate = .false.

! Initialize variables to avoid compiler warnings

ivart = -1

! --- Initialisation memoire

! This routine computes the density and the thermodynamic temperature.
! The computations require the pressure profile which is here taken from
! the meteo file. If no meteo file is used, the user should
! give the laws for RHO and T in usphyv.f90

if (imeteo.eq.0) return

!===============================================================================

!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable thermique
!       (et de ses conditions limites)
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

if (iscalt.gt.0) then
  ivart = isca(iscalt)
else
  write(nfecra,9010) iscalt
  call csexit (1)
endif

! --- Masse volumique

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

call field_get_val_s(ivarfl(ivart), cvar_vart)
if (ippmod(iatmos).ge.2) call field_get_val_s(ivarfl(isca(itotwt)), cvar_totwt)

ipctem = ipproc(itempc)

! From potential temperature, compute:
! - Temperature in Celsius
! - Density
! ----------------------

! Computes the perfect gaz constants according to the physics

rhum = rair
rscp = rair/cp0

lrhum = rair
lrscp = rair/cp0

do iel = 1, ncel

  xvart = cvar_vart(iel) !  The thermal scalar is potential temperature

  if (ippmod(iatmos).ge.2) then  ! humid atmosphere
    lrhum = rair*(1.d0 + (rvsra - 1.d0)*cvar_totwt(iel))
    lrscp = (rair/cp0)*(1.d0 + (rvsra - cpvcpa)*                    &
            cvar_totwt(iel))
  endif

  zent = xyzcen(3,iel)

  if (imeteo.eq.0) then
    call atmstd(zent,pp,dum,dum)
  else
    ! Pressure profile from meteo file:
    call intprf &
    !==========
       ( nbmett, nbmetm,                                            &
         ztmet , tmmet , phmet , zent, ttcabs, pp )
  endif

  ! Temperature in Celsius in cell centers:
  ! ---------------------------------------
  ! law: T = theta * (p/ps) ** (Rair/Cp0)

  propce(iel, ipctem) = xvart*(pp/ps)**lrscp
  propce(iel, ipctem) = propce(iel, ipctem) - tkelvi

  !   Density in cell centers:
  !   ------------------------
  !   law:    RHO       =   P / ( Rair * T(K) )

  crom(iel) = pp/(lrhum*xvart)*(ps/pp)**lrscp

enddo

if (ippmod(iatmos).ge.2) then ! humid atmosphere physics
  ipcliq = ipproc(iliqwt)

  if (moddis.eq.1)then ! all or nothing condensation scheme
    call all_or_nothing()
  elseif (moddis.ge.2)then ! gaussian subgrid condensation scheme
    call gaussian()
  endif
endif ! (ippmod(iatmos).ge.2)

! User re-definition:
call usatph ()

!===============================================================================
! FORMATS
!----

9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    APPEL A csexit DANS LE SOUS PROGRAMME atphyv            ',/,&
'@                                                            ',/,&
'@    La variable dont dependent les proprietes physiques ne  ',/,&
'@      semble pas etre une variable de calcul.               ',/,&
'@    En effet, on cherche a utiliser la temperature alors que',/,&
'@      ISCALT = ',I10                                         ,/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usphyv (et le test lors de la     ',/,&
'@      definition de IVART).                                 ',/,&
'@    Verifier la definition des variables de calcul dans     ',/,&
'@      usipsu. Si un scalaire doit jouer le role de la       ',/,&
'@      temperature, verifier que ISCALT a ete renseigne.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return
contains

! *******************************************************************
! *
! *******************************************************************

subroutine all_or_nothing()

lrhum = rhum

do iel = 1, ncel

  zent = xyzcen(3,iel)

  if (imeteo.eq.0) then
    call atmstd(zent,pp,dum,dum)
  else
    ! Pressure profile from meteo file:
    call intprf &
         !   ===========
       ( nbmett, nbmetm,                                            &
         ztmet , tmmet , phmet , zent, ttcabs, pp )
  endif

  xvart = cvar_vart(iel) ! thermal scalar: liquid potential temperature
  tliq = xvart*(pp/ps)**rscp ! liquid temperature
  qwt  = cvar_totwt(iel) !total water content
  qsl = qsatliq(tliq, pp) ! saturated vapor content
  deltaq = qwt - qsl

  if (activate) then
    write(nfecra,*)"atphyv::all_or_nothing::xvart = ",xvart
    write(nfecra,*)"atphyv::all_or_nothing::tliq = ",tliq
    write(nfecra,*)"atphyv::all_or_nothing::qwt = ",qwt
    write(nfecra,*)"atphyv::all_or_nothing::qsl = ",qsl
    write(nfecra,*)"atphyv::all_or_nothing::qwt,qsl,deltaq = ",qwt,qsl,deltaq
    write(nfecra,*)"atphyv::all_or_nothing::zc = ",xyzcen(3,iel)
    write(nfecra,*)"atphyv::all_or_nothing::pp = ",pp
    write(nfecra,*)"atphyv::all_or_nothing::p0 = ",ps
    write(nfecra,*)"atphyv::all_or_nothing::zent = ",zent
  endif

  if (deltaq.le.0.d0) then ! unsaturated air parcel
    lrhum = rair*(1.d0 + (rvsra - 1.d0)*qwt)
    !Celcius temperature of the air parcel
    propce(iel, ipctem) = tliq - tkelvi
    !density of the air parcel
    crom(iel) = pp/(lrhum*tliq)
    !liquid water content
    propce(iel,ipcliq) = 0.d0
    nebdia(iel) = 0.d0
    nn(iel) = 0.d0
  else ! saturated (ie. with liquid water) air parcel
    qliq = deltaq/ &
         (1.d0 + qsl*clatev**2/(rair*rvsra*cp0*tliq**2))
    lrhum = rair*(1.d0 - qliq + (rvsra - 1.d0)*(qwt - qliq))
    ! liquid water content
    propce(iel,ipcliq) = qliq
    ! Celcius temperature of the air parcel
    propce(iel, ipctem) = tliq + (clatev/cp0)*qliq - tkelvi
    ! density
    crom(iel) = pp/(lrhum*(tliq + (clatev/cp0)*qliq))
    nebdia(iel) = 1.d0
    nn(iel) = 0.d0
  endif

enddo ! iel = 1, ncel
end subroutine all_or_nothing

! *******************************************************************
! *
! *******************************************************************

subroutine gaussian()
! subgrid condensation scheme assuming a gaussian distribution for the
! fluctuations of both qw and thetal.
double precision, dimension(:,:), allocatable :: dtlsd
double precision, dimension(:,:), allocatable :: dqsd

double precision a_const
double precision a_coeff
double precision alpha,al
double precision sig_flu ! standard deviation of qw'-alpha*theta'
double precision var_tl,var_q,cov_tlq
double precision q1,qsup

double precision, dimension(:), pointer :: cvar_k, cvar_ep

! rvap = rair*rvsra

allocate(dtlsd(3,ncelet))
allocate(dqsd(3,ncelet))

! computation of grad(thetal)
call grad_thetal(dtlsd)

! computation of grad(qw)
call grad_qw(dqsd)

call field_get_val_s(ivarfl(ik), cvar_k)
call field_get_val_s(ivarfl(iep), cvar_ep)

! -------------------------------------------------------------
! Gradients are used for estimating standard deviations of the
! subgrid fluctuations.
! -------------------------------------------------------------

lrhum = rhum

a_const = 2.d0*cmu/2.3d0
do iel = 1, ncel

  a_coeff = a_const*cvar_k(iel)**3/cvar_ep(iel)**2 ! 2 cmu/c2 * k**3 / eps**2
  var_tl= a_coeff*(dtlsd(1,iel)**2 + dtlsd(2,iel)**2 + dtlsd(3,iel)**2)
  var_q = a_coeff*( dqsd(1,iel)**2 + dqsd(2,iel)**2 + dqsd(3,iel)**2)
  cov_tlq = a_coeff*(  dtlsd(1,iel)*dqsd(1,iel)   &
                     + dtlsd(2,iel)*dqsd(2,iel)   &
                     + dtlsd(3,iel)*dqsd(3,iel))

  zent = xyzcen(3,iel)

  if (imeteo.eq.0) then
    call atmstd(zent,pp,dum,dum)
  else
    ! Pressure profile from meteo file:
    call intprf(nbmett, nbmetm, ztmet , tmmet , phmet , zent, ttcabs, pp)
  endif

  xvart = cvar_vart(iel) ! thermal scalar: liquid potential temperature
  tliq = xvart*(pp/ps)**rscp ! liquid temperature
  qwt  = cvar_totwt(iel) ! total water content
  qsl = qsatliq(tliq, pp) ! saturated vapor content
  deltaq = qwt - qsl
  alpha = (clatev*qsl/(rvap*tliq**2))*(pp/ps)**rscp
  sig_flu = sqrt(var_q + alpha**2*var_tl - 2.d0*alpha*cov_tlq)

  if (sig_flu.lt.1.d-30) sig_flu = 1.d-30
  q1 = deltaq/sig_flu
  al = 1.d0/(1.d0 + qsl*clatev**2/(rair*rvsra*cp0*tliq**2))
  qsup = qsl/sig_flu

  nebdia(iel) = 0.5d0*(1.d0 + ferf(q1/sqrt(2.d0)))

  qliq = (sig_flu                                                               &
        /(1.d0 + qsl*clatev**2/(rvap*cp0*tliq**2)))                             &
        *(nebdia(iel)*q1 + exp(-q1**2/2.d0)/sqrt(2.d0*pi))
  qliq = max(qliq,1d-15)
  nn(iel) = nebdia(iel) - (nebdia(iel)*q1                                       &
          + exp(-q1**2/2.d0)/sqrt(2.d0*pi))*exp(-q1**2/2.d0)/sqrt(2.d0*pi)

  if(qwt.lt.qliq)then
    ! go back to all or nothing
    if (deltaq.le.0.d0) then ! unsaturated air parcel
      lrhum = rair*(1.d0 + (rvsra-1.d0)*qwt)
      !Celcius temperature of the air parcel
      propce(iel, ipctem) = tliq - tkelvi
      !density of the air parcel
      crom(iel) = pp/(lrhum*tliq)
      !liquid water content
      propce(iel,ipcliq) = 0.d0
      nebdia(iel) = 0.d0
      nn(iel) = 0.d0
    else ! saturated (ie. with liquid water) air parcel
      qliq = deltaq / (1.d0 + qsl*clatev**2/(rair*rvsra*cp0*tliq**2))
      lrhum = rair*(1.d0 - qliq + (rvsra - 1.d0)*(qwt - qliq))
      ! liquid water content
      propce(iel,ipcliq) = qliq
      ! Celcius temperature of the air parcel
      propce(iel,ipctem) = tliq+(clatev/cp0)*qliq - tkelvi
      ! density
      crom(iel) = pp/(lrhum*(tliq + (clatev/cp0)*qliq))
      nebdia(iel) = 1.d0
      nn(iel) = 0.d0
    endif
  else ! coherent subgrid diagnostic
    lrhum = rair*(1.d0 - qliq + (rvsra - 1.d0)*(qwt - qliq))
    ! liquid water content
    propce(iel,ipcliq) = qliq
    !Celcius temperature of the air parcel
    propce(iel, ipctem) = tliq + (clatev/cp0)*qliq - tkelvi
    !density
    crom(iel) = pp/(lrhum*(tliq + (clatev/cp0)*qliq))
  endif ! qwt.lt.qliq

enddo

! when properly finished deallocate dtlsd
deallocate(dtlsd)
deallocate(dqsd)

end subroutine gaussian

! *******************************************************************
! *
! *******************************************************************

subroutine grad_thetal(dtlsd)
double precision dtlsd(3,ncelet)

integer    iccocg
integer    iivar
integer    inc
integer    itpp

! Computation of the gradient of the potential temperature

itpp = isca(iscalt)

! options for gradient calculation

iccocg = 1
inc = 1

iivar = itpp

call field_gradient_scalar(ivarfl(iivar), 1, imrgra, inc,           &
                           iccocg,                                  &
                           dtlsd)

end subroutine grad_thetal

! *******************************************************************
! *
! *******************************************************************

subroutine grad_qw(dqsd)
double precision dqsd(3,ncelet)

integer    iccocg
integer    iivar
integer    inc
integer    iqw

! ----------------------------------------------------------------
! now gradient of total humidity
! ----------------------------------------------------------------

iccocg = 1
inc = 1

iqw = isca(itotwt)

iivar = iqw

call field_gradient_scalar(ivarfl(iivar), 1, imrgra, inc,           &
                           iccocg,                                  &
                           dqsd)

end subroutine grad_qw
end subroutine atphyv

