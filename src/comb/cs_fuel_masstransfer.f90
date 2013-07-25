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

subroutine cs_fuel_masstransfer &
!=============================
 ( ncelet , ncel   ,            &
   rtpa   , propce , volume )

!===============================================================================
! FONCTION :
! --------
! CALCUL DES TERMES DE TRANSFERT DE MASSE ENTRE LA PHASE CONTINUE
! ET LA PHASE DISPERSEE
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant precedent)                !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
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
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_fuel_incl

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel

double precision rtpa(ncelet,*), propce(ncelet,*)
double precision volume(ncelet)

! Local variables

integer          iel    , icla
integer          ipcrom , ipcte1 , ipcte2 , ipcro2 , ipcdia
integer          ipcgev , ipcght , ipcyox
integer          ipcvst,ipcvsl,ipccp,ipchgl

double precision xng,xnuss
double precision pparo2 , xdffli , xdfext , xdftot0 , xdftot1
double precision diacka, xuash
double precision dcoke , surf , lambda
!
double precision peva, pref
double precision rhoeva
double precision ueva
double precision deva
double precision obfl

double precision dhet1, dhet2
double precision deva1, deva2
!
!===============================================================================
! 1. INITIALISATIONS ET CALCULS PRELIMINAIRES
!===============================================================================
! --- Initialisation des termes de transfert de masse

do icla = 1, nclafu
  ipcgev = ipproc(igmeva(icla))
  ipcght = ipproc(igmhtf(icla))
  ipchgl = ipproc(ih1hlf(icla))
  do iel = 1, ncel
    propce(iel,ipcgev) = zero
    propce(iel,ipcght) = zero
    propce(iel,ipchgl) = zero
  enddo
enddo

! --- Pointeur

ipcrom = ipproc(irom)
ipcte1 = ipproc(itemp1)
ipcyox = ipproc(iym1(io2))
ipcvst = ipproc(ivisct)
!
  pref = 1.013d0
!===============================================================================
! 2. TERMES SOURCES POUR l'ENTHALPIE LIQUIDE
!===============================================================================
!
! Contribution aux bilans explicite et implicite
! des echanges par diffusion moleculaire
! 6 Lambda Nu / diam**2 / Rho2 * Rho * (T1-T2)
!
do icla = 1, nclafu

  ipcro2 = ipproc(irom2 (icla))
  ipcdia = ipproc(idiam2(icla))
  ipcte2 = ipproc(itemp2(icla))
  ipcght = ipproc(igmhtf(icla))
  ipchgl = ipproc(ih1hlf(icla))

  xnuss = 2.d0
  do iel = 1, ncel
    if ( ivisls(ihm).gt.0 ) then
      ipcvsl = ipproc(ivisls(ihm))
      if ( icp.gt.0 ) then
        ipccp   = ipproc(icp)
        lambda = propce(iel,ipcvsl) * propce(iel,ipccp)
      else
        lambda = propce(iel,ipcvsl) * cp0
      endif
    else
      if ( icp.gt.0 ) then
        ipccp  = ipproc(icp)
        lambda = visls0(ihm) * propce(iel,ipccp)
      else
        lambda = visls0(ihm) * cp0
      endif
    endif
!
    if ( rtpa(iel,isca(iyfol(icla))) .gt. epsifl  .and.                    &
         propce(iel,ipcte1).gt. propce(iel,ipcte2)        ) then
!
       propce(iel,ipchgl) = 6.d0*lambda*xnuss/propce(iel,ipcdia)**2        &
                           /propce(iel,ipcro2)*rtpa(iel,isca(iyfol(icla)))
!
    endif
!
  enddo
!
enddo
!
!===============================================================================
! 3. TRANSFERTS DE MASSE FIOUL
!===============================================================================
!
do icla = 1, nclafu

  ipcro2 = ipproc(irom2 (icla))
  ipcdia = ipproc(idiam2(icla))
  ipcte2 = ipproc(itemp2(icla))
  ipcgev = ipproc(igmeva(icla))
  ipchgl = ipproc(ih1hlf(icla))
  ipcght = ipproc(igmhtf(icla))
!
  do iel = 1, ncel
!
    propce(iel,ipcgev) = zero
    propce(iel,ipcght) = zero
!
    if (rtpa(iel,isca(iyfol(icla))) .gt. epsifl ) then
!
!===============================================================================
! EVAPORATION
!===============================================================================
! Verification sur la masse du fioul liquide.
! a) Si deva1 < deva2 il ne reste plus de fioul liquide.
!    Du coup pas d'evaporation.
! b) Verification sur la plage des temperatures d'evaporation.
! c) Verification sur la temperature de la phase gaz et de la goutte. Il faut
!     Tgaz > Tgoutte
!
    deva1 =  rtpa(iel,isca(iyfol(icla)))                                   &
             /(rtpa(iel,isca(ing(icla)))*rho0fl)
    deva2 =  (pi*(diniin(icla)**3)/6.d0)+(pi*(dinikf(icla)**3)/6.d0)
!
      if ( propce(iel,ipcte2)    .gt. tevap1               .and.           &
           propce(iel,ipcte1)    .gt. propce(iel,ipcte2)   .and.           &
           deva1.gt.deva2                                          ) then
!
! La flux de masse evapore est determinee en fonction d'un profil
!                    dMeva/dTgoutte supposé.
!
      propce(iel,ipcgev) = propce(iel,ipchgl)                              &
                          /( hrfvap + cp2fol*(tevap2-propce(iel,ipcte2)) )
!
      endif
!
!===============================================================================
! COMBUSTION HETEROGENE
!===============================================================================
! Verification sur la masse du coke.
! a) Si deva1.le.deva2 -> Il ne reste plus de fioul liquide. Du coup la particule
!    est constituee de charbon et d'inertes.
! b) Si dhet1.gt.dhet2 il reste du charbon. Si non la particule est constituee
!    que des inertes.
!
    dhet1= rtpa(iel,isca(iyfol(icla)))                                     &
            /(rtpa(iel,isca(ing(icla)))*rho0fl)
    dhet2= pi*(diniin(icla)**3)/6.d0
!
      if (deva1.le.deva2.and.dhet1.gt.dhet2 ) then
!
! On considere la masse du coke restant comme une particule spherique. Le
!   diametre correspendant est dcoke.
      dcoke = ( ( rtpa(iel,isca(iyfol(icla)))                              &
              /(rtpa(iel,isca(ing(icla)))*rho0fl)                          &
              -pi*(diniin(icla)**3)/6.d0  )                                &
               *6.d0/pi )**(1.d0/3.d0)
!
! Calcul de la pression partielle en oxygene (atm)                                                 ---
!   PO2 = RHO1*RR*T*YO2/MO2
!
      pparo2 = propce(iel,ipproc(irom1))*rr*propce(iel,ipcte1)             &
              *propce(iel,ipcyox)/wmole(io2)
      pparo2 = pparo2 / prefth
!
! Coefficient de cinetique chimique de formation de CO
!   en (kg.m-2.s-1.atm(-n))
      xdffli = ahetfl*exp(-ehetfl*4185.d0                                  &
              /(rr*propce(iel,ipcte1)))
!
! Coefficient de diffusion en  (Kg/m2/s/atm) : XDFEXT
! Coefficient global pour n=0.5 en (kg/m2/s) : XDFTOT0
! Coefficient global pour n=1   en (Kg/m2/s) : XDFTOT1
!
      diacka = dcoke/(dinikf(icla))
        if ( diacka .gt. epsifl ) then
        xdfext = 2.53d-7*((propce(iel,ipcte1))**0.75d0)                    &
                / dcoke*2.d0
        xdftot1 = pparo2 / ( 1.d0/xdffli + 1.d0/xdfext )
        xdftot0 = -(xdffli**2)/(2.d0*xdfext**2)+(pparo2*xdffli**2          &
                +(xdffli**4)/(2.d0*xdfext**2))**0.5d0
        else
        xdftot1 = xdffli*pparo2
        xdftot0 = xdffli*pparo2**0.5d0
        endif
!
! Surface de la particule spherique.
      surf = pi*(dcoke**2)
!
! Nombre des particules dans la cellule.
!
      xng   = rtpa(iel,isca(ing(icla)))
!
        if (iofhet.eq.1) then
        propce(iel,ipcght) = - xdftot1*surf*xng
        else
        propce(iel,ipcght) = - xdftot0*surf*xng
        endif

      endif

    endif
!
  enddo
!
enddo
!
!----
! End
!----

return
end subroutine
