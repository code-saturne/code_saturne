!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine ppcabs &
!================

 ( propce , tempk, kgas , agas , agasbo )

!===============================================================================
! FONCTION :
! ----------


!   SOUS-PROGRAMME PHYSIQUES PARTICULIERES

!  DONNE LA VALEUR DU COEFFICIENT D'ABSORPTION POUR
!    LE MELANGE GAZEUX ET LES PARTICULES POUR LE CP.

!  DANS LE CAS DU MODELE P-1 ON VERIFIE QUE LA LONGUEUR OPTIQUE
!    DU MILIEU EST AU MINIMUM DE L'ORDRE DE L'UNITE


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! tempk            ! ra ! <-- ! Gas phase temperature                          !
! kgas             ! ra ! <-- ! Radiation coeffcients of the i different grey  !
!                  !    !     ! gases                                          !
! agas             ! ra ! <-- ! Weights of the of the i different grey gases   !
!                  !    !     ! in cell                                        !
! agasbo           ! ra ! <-- ! Weights of the of the i different grey gases   !
!                  !    !     ! at boundary faces                              !
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
use entsor
use optcal
use cstphy
use cstnum
use parall
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use radiat
use mesh
use field

!===============================================================================

implicit none

! Arguments

double precision propce(ncelet,*)

double precision tempk(ncelet,nrphas)
double precision kgas(ncelet,nwsgg), agas(ncelet,nwsgg), agasbo(nfabor,nwsgg)

! Local variables

integer          iel, ifac, icla, ipck, icha, iok
double precision xm, dd2, vv, sf, xlc, xkmin, pp

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cvar_fsm

!===============================================================================
! 0 - Initialization
!===============================================================================

if (imodak.eq.1.or.imoadf.ge.1.or.imfsck.eq.1) then
  allocate(w1(ncelet), w2(ncelet), w3(ncelet))
endif

call field_get_val_s(icrom, crom)

!===============================================================================
!  1 - COEFFICIENT D'ABSORPTION DU MELANGE GAZEUX (m-1)
!===============================================================================

if ( ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0 ) then

! ----> Combustion gaz : Flamme de diffusion
!                        Flamme de premelange (Modele EBU)

  if (imodak.eq.1) then

    if (isoot.ge.1) call field_get_val_s(ivarfl(isca(ifsm)), cvar_fsm)

    do iel = 1, ncel
      xm = 1.d0/ (  propce(iel,ipproc(iym(1)))/wmolg(1)                 &
                  + propce(iel,ipproc(iym(2)))/wmolg(2)                 &
                  + propce(iel,ipproc(iym(3)))/wmolg(3) )
      w1(iel) = propce(iel,ipproc(iym(3)))*xm/wmolg(3)*xco2
      w2(iel) = propce(iel,ipproc(iym(3)))*xm/wmolg(3)*xh2o

      w3(iel) = 0.d0

      ! Soot model
      if (isoot.eq.0) w3(iel) = Xsoot * propce(iel,ipproc(iym(1))) &
                       * crom(iel) / rosoot
      if (isoot.ge.1) w3(iel) = cvar_fsm(iel) &
                        * crom(iel) / rosoot
    enddo
    call raydak(ncel,ncelet,                                      &
    !==========
      propce(1,ipproc(icak(1))),w1,w2,w3,propce(1,ipproc(itemp)))

    ! the code seems to be good (BS)
    if (ntcabs.eq.ntpabs+1) then
      write(NFECRA,*) ' a verifier '
      write(NFECRA,*) ' a finir   : raydak '
      write(NFECRA,*) ' Le codage est a terminer par le groupe I81'
      write(NFECRA,*) '                     13-10-03 22:38:03      '
      ! call csexit(1)
    endif

  else
    do iel = 1, ncel
      propce(iel,ipproc(icak(1))) = propce(iel,ipproc(ickabs))
    enddo
  endif

else if ( ippmod(iccoal) .ge. 0 ) then

  ! ---->  Coal

  if (imodak.eq.1) then

    do iel = 1, ncel
      ! CO2 volume concentration
      w1(iel) = propce(iel,ipproc(immel))/wmole(ico2)             &
      *propce(iel,ipproc(iym1(ico2)))
      ! H2O volume concentration
      w2(iel) = propce(iel,ipproc(immel))/wmole(ih2o)             &
      *propce(iel,ipproc(iym1(ih2o)))
      ! Soot volume fraction
      w3(iel) = 0.d0
    enddo

    call raydak(ncel,ncelet,                                      &
    !==========
     propce(1,ipproc(icak(1))),w1,w2,w3,propce(1,ipproc(itemp1)))

  else if (imoadf.ge.1) then

    do iel = 1, ncel
      ! CO2 volume concentration
      w1(iel) = propce(iel,ipproc(immel))/wmole(ico2)             &
      *propce(iel,ipproc(iym1(ico2)))
      ! H2O volume concentration
      w2(iel) = propce(iel,ipproc(immel))/wmole(ih2o)             &
      *propce(iel,ipproc(iym1(ih2o)))
      ! Soot volume fraction
      w3(iel) = 0.d0
    enddo

    if (imoadf.eq.1) then
      call radf08 (w1,w2,w3,tempk(1,1),kgas,agas,agasbo)
    else if (imoadf.eq.2) then
      call radf50 (w1,w2,w3,tempk(1,1),kgas,agas,agasbo)
    endif

  else if (imfsck.eq.1) then

    do iel = 1, ncel
      ! CO2 volume concentration
      w1(iel) = propce(iel,ipproc(immel))/wmole(ico2)             &
      *propce(iel,ipproc(iym1(ico2)))
      ! H2O volume concentration
      w2(iel) = propce(iel,ipproc(immel))/wmole(ih2o)             &
      *propce(iel,ipproc(iym1(ih2o)))
      ! Soot volume fraction
      w3(iel) = 0.d0
    enddo

    call rafsck (w1,w2,w3,tempk(1,1),kgas,agas,agasbo)
  else

    do iel = 1, ncel
      propce(iel,ipproc(icak(1))) = ckabs1
    enddo

  endif

else if ( ippmod(icfuel).ge.0 ) then

! ---->  Fuel

  if (imodak.eq.1) then

    do iel = 1,ncel
! concentration volumique en CO2
      w1(iel) = propce(iel,ipproc(immel))/wmole(ico2)             &
               *propce(iel,ipproc(iym1(ico2)))
! concentration volumique en H20
      w2(iel) = propce(iel,ipproc(immel))/wmole(ih2o)             &
               *propce(iel,ipproc(iym1(ih2o)))
! fraction volumique de suies
      w3(iel) = 0.d0

    enddo

    call raydak(ncel,ncelet,                                      &
    !==========
     propce(1,ipproc(icak(1))),w1,w2,w3,propce(1,ipproc(itemp1)))

  else if (imoadf.ge.1) then

    do iel = 1,ncel
! concentration volumique en CO2
      w1(iel) = propce(iel,ipproc(immel))/wmole(ico2)             &
               *propce(iel,ipproc(iym1(ico2)))
! concentration volumique en H20
      w2(iel) = propce(iel,ipproc(immel))/wmole(ih2o)             &
               *propce(iel,ipproc(iym1(ih2o)))
! fraction volumique de suies
      w3(iel) = 0.d0
    enddo

    if (imoadf.eq.1) then
      call radf08 (w1,w2,w3,tempk(1,1),kgas,agas,agasbo)
    elseif (imoadf.eq.2) then
      call radf50 (w1,w2,w3,tempk(1,1),kgas,agas,agasbo)
    endif

  else

    do iel = 1, ncel
      propce(iel,ipproc(icak(1))) = ckabs1
    enddo

  endif

endif

if (imodak.eq.1.or.imoadf.ge.1) deallocate(w1, w2, w3)

!===============================================================================
!  2 - COEFFICIENT D'ABSORPTION DES PARTICULES PAR CLASSE K2/X2 (m-1)
!===============================================================================

! ---->  Charbon

if ( ippmod(iccoal) .ge. 0 ) then

  do icla = 1, nclacp

    ipck = 1 + icla
    icha = ichcor(icla)

    do iel = 1, ncel

! ---> Calcul du diametre des particules

      dd2 = ( xashch(icha)*diam20(icla)**2 +                      &
           ( 1.d0-xashch(icha))                                   &
             *propce(iel,ipproc(idiam2(icla)))**2 )**0.5d0

! ---> Calcul du coeficient d'absorption des particules K2/X2
!         3./2. ROM/(ROM2*DD2)

      propce(iel,ipproc(icak(ipck))) = 1.5d0*crom(iel)            &
                       / ( propce(iel,ipproc(irom2(icla)))*dd2)

    enddo

  enddo

endif

! ---->  Fuel

if ( ippmod(icfuel).ge.0 ) then

  do icla = 1, nclafu

    ipck = 1 + icla

    do iel = 1, ncel

! ---> Calcul du coeficient d'absorption des particules K2/X2
!         3./2. ROM/(ROM2*DD2)

      propce(iel,ipproc(icak(ipck))) =  1.5d0*crom(iel)           &
                 / ( propce(iel,ipproc(irom2(icla)))              &
                    *propce(iel,ipproc(idiam2(icla))) )

    enddo

  enddo

endif

!===============================================================================
!  3 - COEFFICIENT D'ABSORPTION GAZ (Arc Electrique)
!===============================================================================


if ( ippmod(ielarc).ge.1 ) then

  do iel = 1, ncel

! ---> Directement donne par le fichier dp_elec

      propce(iel,ipproc(icak(1))) = propce(iel,ipproc(idrad))

  enddo

endif

!===============================================================================
!  4 - CLIPPING DU COEFFICIENT D'ABSORPTION DANS LA CAS DE L'APPROX P-1
!===============================================================================


!--> MODELE P-1 : Controle standard des valeurs du coefficient
!                 d'absorption. Ce coefficient doit assurer une
!                 longueur optique au minimum de l'ordre de l'unite.
!
!                 !!!!DO NOT USE WITH ADF model!!!!

  if (iirayo.eq.2.and.imoadf.eq.0) then

     allocate(w3(ncelet))

!         Coefficient d'absorption du melange gaz-particules de charbon

     do iel = 1, ncel
       w3(iel) =  propce(iel,ipproc(icak(1)))
     enddo

     if ( ippmod(iccoal).ge.0 ) then
       do icla = 1,nclacp
         ipck = 1+icla
         do iel = 1,ncel
           w3(iel) = w3(iel)                                      &
                    + ( propce(iel,ipproc(ix2(icla)))             &
                      * propce(iel,ipproc(icak(ipck))) )
         enddo
       enddo
     elseif ( ippmod(icfuel).ge.0 ) then
       do icla = 1,nclafu
         ipck = 1+icla
         do iel = 1,ncel
           w3(iel) = w3(iel)                                      &
                    + ( propce(iel,ipproc(iyfol(icla)))           &
                      * propce(iel,ipproc(icak(ipck))) )
         enddo
       enddo
     endif

!         Calcul de la longueur caractéristique XLC du domaine de calcul

    sf = 0.d0
    vv = 0.d0

    do ifac = 1,nfabor
       sf = sf + sqrt(surfbo(1,ifac)**2 + surfbo(2,ifac)**2 + surfbo(3,ifac)**2)
    enddo
    if (irangp.ge.0) then
      call parsom(sf)
      !==========
    endif

    do iel = 1,ncel
       vv = vv + volume(iel)
    enddo
    if (irangp.ge.0) then
      call parsom(vv)
      !==========
    endif

    xlc = 3.6d0 * vv / sf

!         Clipping de CK a XKMIN

    xkmin = 1.d0 / xlc

    iok = 0
    do iel = 1,ncel
      if (w3(iel).lt.xkmin) then
        iok = iok +1
      endif
    enddo
    if (irangp.ge.0) then
      call parcpt(iok)
    endif

!     Arret en fin de pas de temps si epaisseur optique trop grande
    pp = xnp1mx/100.0d0
    if (dble(iok).gt.pp*dble(ncelgb)) then
       write(nfecra,1000) xkmin, dble(iok)/dble(ncelgb)*100.d0, xnp1mx
       istpp1 = 1
    endif

    ! Free memory
    deallocate(w3)

  endif

! -------
! FORMAT
! -------

 1000 format(                                                    &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Radiative module with P-1 approximation        ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    The optical length of the semi-transparent medium       ',/,&
'@      must be at least approximately 1 so as to be in the   ',/,&
'@      application domain of the P-1 approximation.          ',/,&
'@    This does not seem to be the case here.                 ',/,&
'@                                                            ',/,&
'@    The minimum absorption coefficient to ensure this       ',/,&
'@      optical length is xkmin = ',e11.4                     ,/,&
'@    This value is not reached for ', e11.4,'%               ',/,&
'@      of the meshe''s cells.                                ',/,&
'@    The percentage of mesh cells for which we allow this    ',/,&
'@      condition to be violated is fixed by defaul or in     ',/,&
'@      cs_user_parameters.f90 to xnp1mx = ', e11.4,'%        ',/,&
'@                                                            ',/,&
'@    The computation will not be run.                        ',/,&
'@                                                            ',/,&
'@    Check the values of the absorption coefficent Ck        ',/,&
'@      in ppcabs, usray3, or the thermochemistry file.       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return

end subroutine
