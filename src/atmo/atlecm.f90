!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine atlecm &
!================

 ( imode  , tmprom , ztprom , zdprom ,                            &
   xmet   , ymet   , pmer   ,                                     &
   ttprom , qvprom ,                                              &
   uprom  , vprom  , ekprom , epprom,                             &
   rprom  , tpprom , phprom   )

!===============================================================================
!  Purpose:
!  -------


!             Reads the meteo profile data
!             for the atmospheric module
!
!             IMODE = 0 : READING OF DIMENSIONS ONLY
!             IMODE = 1 : READING OF DATA

!             WARNING : ARGUMENTS ARE DEFINED ONLY ON SECOND CALL
!               OF THE ROUTINE (IMODE=1)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use ppppar

use atincl

!===============================================================================

implicit none

! Arguments

integer           imode

double precision  tmprom(nbmetm)
double precision  ztprom(nbmett) , zdprom(nbmetd)
double precision  xmet(nbmetm)   , ymet(nbmetm)  , pmer(nbmetm)
double precision  ttprom(nbmett,nbmetm) , qvprom(nbmett,nbmetm)
double precision  uprom(nbmetd,nbmetm)  , vprom(nbmetd,nbmetm)
double precision  ekprom(nbmetd,nbmetm) , epprom(nbmetd,nbmetm)
double precision  rprom(nbmett,nbmetm)  , tpprom(nbmett,nbmetm)
double precision  phprom(nbmett,nbmetm)

! Local variables

integer itp, ii, ios, k, iphas
integer syear, squant, shour, smin
integer year, quant,hour,minute

double precision ssec, second
double precision sjday, jday
double precision psol,rap,rscp,tmoy

character*80     ccomnt
character*1      csaute

!===============================================================================

if (imode.eq.0) then
  write(nfecra, *) 'reading of dimensions for meteo profiles'
else
  write(NFECRA, *) 'reading of meteo profiles data'
endif

!===============================================================================
! 0. Initialization
!===============================================================================

CSAUTE = '/'

! --> Opens the meteo file
open (unit=impmet, file=ficmet,                                  &
      status='old', form='formatted', access='sequential',       &
      iostat=ios, err=99)
rewind(unit=impmet, err=99)

itp=0

syear=-999
squant=-999
shour=-999
smin=-999
ssec=-999.

!===============================================================================
! 1. Loop on time
!===============================================================================

 100  continue
itp = itp+1

! ---> Read the comments
! ----------------------

 101  read(impmet, '(a80)', err=999, end=906) ccomnt

if(ccomnt(1:1).eq.csaute) go to 101
backspace(impmet)

!===============================================================================
! 2. Read the time of the profile
!===============================================================================

! --> year, quant-day, hour, minute, second  of the profile (TU)

if (imode.eq.0) then
  read(impmet, *, err=999, end=906)
else
  read(impmet, *, err=999, end=906) year, quant, hour, minute, second

  ! --> the date and time of the first meteo profile are taken as the
  !     starting time of the simulation

  if (syear.lt.0) then
    syear=year
    squant=quant
    shour=hour
    smin=minute
    ssec=second
  endif

  !--> Compute the relative time to the starting time of the simulation


  ! --> Compute the julian day for the starting day of the simulation
  !     (julian day at 12h)
  sjday= squant + ((1461 * (syear + 4800 + (1 - 14) / 12)) / 4 +   &
             (367 * (1 - 2 - 12 * ((1 - 14) / 12))) / 12 -         &
             (3 * ((syear + 4900 + (1 - 14) / 12) / 100)) / 4      &
              + 1 - 32075) - 1

  ! --> Compute the julian day for the date of the current profile
  !     (julian day at 12h)
  jday = quant + ((1461 * (year + 4800 + (1 - 14) / 12)) / 4 +     &
             (367 * (1 - 2 - 12 * ((1 - 14) / 12))) / 12 -         &
             (3 * ((year + 4900 + (1 - 14) / 12) / 100)) / 4       &
              + 1 - 32075) - 1

  tmprom(itp) = (jday - sjday)*86400 + (hour - shour)*3600 + (minute - smin)*60 &
              + (second - ssec)

  ! --> check the chronlogical order of profiles

  if (itp.gt.1) then
    if (tmprom(itp).lt.tmprom(itp-1)) then
      write(nfecra, 8000)
      call csexit (1)
      !==========
    endif
  endif

endif

!===============================================================================
! 3. Read the position of the profile
!===============================================================================

 102  read(impmet, '(a80)', err=999, end=999) ccomnt

if(ccomnt(1:1).eq.csaute) go to 102
backspace(impmet)


if (imode.eq.0) then
  read(impmet, *, err=999, end=999)
else
  read(impmet, *, err=999, end=999) xmet(itp), ymet(itp)
endif

!===============================================================================
! 4. Read the sea-level pressure
!===============================================================================

 103  read(impmet, '(a80)', err=999, end=999) ccomnt

if (ccomnt(1:1).eq.csaute) go to 103
backspace(impmet)


if (imode.eq.0) then
  read(impmet, *, err=999, end=999)
else
  read(impmet, *, err=999, end=999) pmer(itp)
endif

!=================================================================================
! 5. Read the temperature and humidity profiles
!=================================================================================

 104  read(impmet, '(a80)', err=999, end=999) ccomnt

if (ccomnt(1:1).eq.csaute) go to 104
backspace(impmet)


if (imode.eq.0) then
  read(impmet, *, err=999, end=999) nbmett

  if (nbmett.le.1) then
    write(nfecra, 8001)
    call csexit (1)
    !==========
  endif

  do ii=1, nbmett
    read (impmet, *, err=999, end=999)
  enddo

else

  read(impmet, *, err=999, end=999)

  do ii=1, nbmett

    ! Altitude, temperature, humidity
    read (impmet, *, err=999, end=999) ztprom(ii),                       &
                                       ttprom(ii, itp), qvprom(ii, itp)

  enddo

endif

!===============================================================================
! 6. Compute hydro pressure profile  (Laplace integration)
!===============================================================================
!  NOTE : BASEE SUR LA PRESSION PMER
!         ET UNE INTEGRATION DE LAPLACE DE BAS EN HAUT

if (imode.eq.1) then
  iphas = 1
  phprom(1, itp)=pmer(itp)
  psol=p0
  rscp=rair/cp0

  do k=2, nbmett
    tmoy=0.5*(ttprom(k-1, itp)+ttprom(k, itp))+tkelvi
    !         rhmoy=rair*(1.+(rvsra-1.)*
    !    &         (qvprom(k-1, itp)+qvprom(k, itp))/2.*ih2o)
    !         rap=-abs(gz)*(zk-zkm1)/rhmoy/tmoy
    rap=-abs(gz)*(ztprom(k)-ztprom(k-1))/rair/tmoy
    phprom(k, itp)=phprom(k-1, itp)*exp(rap)
  enddo

endif

!===============================================================================
! 7. Compute the pot. temperature profile and the density profile
!===============================================================================

if (imode.eq.1) then
  do k=1, nbmett
    ! rhum=rair*(1.d0+(rvsra-1.d0)*qvprom(k, itp)*ih2o)

    !   if ((iscalt.eq.-1).or.(iphysi.eq.0)) then
    ! Constant density
    !     rprom(k, itp)=pmer/(tprom(k, itp)+tkelvi)/rhum
    !   else
    ! Variable density
    !     rprom(k, itp)=phprom(k, itp)/(ttprom(k, itp)+tkelvi)/rhum
    !   endif
    rprom(k, itp)=phprom(k, itp)/(ttprom(k, itp)+tkelvi)/rair
    !    &            *(1.d0+(rvsra-cpvcpa)*qvprom(k, itp)*ih2o)
    tpprom(k, itp)=(ttprom(k, itp)+tkelvi)*                         &
         ((psol/phprom(k, itp))**rscp)
  enddo

endif

!================================================================================
! 8. Read the velocity profile
!================================================================================

 105  read(impmet, '(a80)', err=999, end=999) ccomnt

if (ccomnt(1:1).eq.csaute) go to 105
backspace(impmet)


if (imode.eq.0) then

  read(impmet, *, err=999, end=999) nbmetd

  if (nbmetd.le.1) then
    write(nfecra, 8002)
    call csexit (1)
    !==========
  endif

  do ii=1, nbmetd
    read (impmet, *, err=999, end=999)
  enddo

else

  read(impmet, *, err=999, end=999)

  do ii=1, nbmetd
    !  Altitude, u, v, k, epsilon
    read (impmet, *, err=999, end=999) zdprom(ii),                    &
                                    uprom(ii, itp),  vprom(ii, itp),  &
                                    ekprom(ii, itp), epprom(ii, itp)
  enddo

endif

!================================================================================
! 9. Printings
!================================================================================

if (imode.eq.1) then
  if (itp.eq.1) then
    write(nfecra, *)
    write(nfecra, *) '==================================================='
    write(nfecra, *) 'printing meteo profiles'
  endif
  write(nfecra, *) 'year, quant-day , hour, minute, second'
  write(nfecra, 7995) year, quant, hour, minute, second
 7995   format(1x, i4, i5, 2i4, f8.2)
  write(nfecra, *) 'tmprom(itp)'
  write(nfecra, 7996) tmprom(itp)
 7996   format(1x, f8.2)
  write(nfecra, *) 'zdprom, uprom, vprom, ekprom, epprom'
  do ii=1, nbmetd
    write(nfecra, 7997)                                            &
      zdprom(ii), uprom(ii, itp), vprom(ii, itp), ekprom(ii, itp), epprom(ii, itp)
 7997   format(1x, 3f8.2, 2e10.3)
  enddo
  write(nfecra, *) 'ztprom, ttprom, tpprom, rprom, phprom, qvprom'
  do ii=1, nbmett
    write(nfecra, 7998)                                             &
         ztprom(ii), ttprom(ii, itp), tpprom(ii, itp),              &
         rprom(ii, itp), phprom(ii, itp), qvprom(ii, itp)
 7998     format(1x, 3f8.2, f8.4, f12.3, e10.3)
  enddo

endif

!================================================================================
! 10. End of the loop on time
!================================================================================

goto 100

 906  continue

if (imode.eq.0) nbmetm= itp-1

close(unit=impmet)

! ---
! End
! ---

return

!============================
! XX. Error outputs
!============================

 99   continue
write ( nfecra, 9998 )
call csexit (1)
!==========

 999  continue
write ( nfecra, 9999 )
call csexit (1)
!==========

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)
 8000 format (                                                    &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (atlecm)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE ATMOSPHERIQUE                   ',/,&
'@                                                            ',/,&
'@              Erreur  dans le fichier meteo:                ',/,&
'@      ordre chronologique des profils non respecte          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8001 format (                                                    &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (atlecm)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE ATMOSPHERIQUE                   ',/,&
'@                                                            ',/,&
'@              Erreur  dans le fichier meteo:                ',/,&
'@  le nombre de mesures de temperatures doit etre            ',/,&
'@  superieur a 2                                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8002 format (                                                    &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (atlecm)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE ATMOSPHERIQUE                   ',/,&
'@                                                            ',/,&
'@              Erreur  dans le fichier meteo:                ',/,&
'@  le nombre de mesures de vitesses doit etre                ',/,&
'@  superieur a 2                                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9998 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (atlecm)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE ATMOSPHERIQUE                   ',/,&
'@                                                            ',/,&
'@  Erreur a l''ouverture du fichier meteo                    ',/,&
'@  Verifier le nom du fichier meteo                          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (atlecm)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE ATMOSPHERIQUE                   ',/,&
'@                                                            ',/,&
'@  Erreur a la lecture du fichier meteo.                     ',/,&
'@    Le fichier a ete ouvert mais est peut etre incomplet    ',/,&
'@    ou son format inadapte.                                 ',/,&
'@                                                            ',/,&
'@    year (integer), quantile (integer), hour (integer),    ',/,&
'@          minute (integer), second (dble prec) of the profile',/,&
'@    location of the meteo profile (x,y) (dble prec)         ',/,&
'@    sea level pressure (double precision)                   ',/,&
'@ temperature profile:                                       ',/,&
'@   number of altitudes (integer)                            ',/,&
'@   alt.,temperature  in celcius,humidity in kg/kg (dble prec)',/,&
'@ wind profile:                                              ',/,&
'@   number of altitudes (integer)                            ',/,&
'@   alt.,u,v,k,eps (double precision)                        ',/,&
'@ NO LINE AT THE END OF THE FILE                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 8000 format (                                                    &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA (atlecm)      ',/,&
'@     =======                                                ',/,&
'@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,&
'@                                                            ',/,&
'@              Error in the meteo profile file:              ',/,&
'@      check that the chronogical order of the profiles      ',/,&
'@      are respected                                         ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8001 format (                                                    &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA (atlecm)      ',/,&
'@     =======                                                ',/,&
'@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,&
'@                                                            ',/,&
'@              Error in the meteo profile file:              ',/,&
'@  the number of temperature measurements must be larger     ',/,&
'@  than 2                                                    ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8002 format (                                                    &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA (atlecm)      ',/,&
'@     =======                                                ',/,&
'@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,&
'@                                                            ',/,&
'@              Error in the meteo profile file:              ',/,&
'@  the number of velocity measurements must be larger        ',/,&
'@  than 2                                                    ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9998 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA (atlecm)      ',/,&
'@     =======                                                ',/,&
'@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,&
'@                                                            ',/,&
'@  Error opening the meteo profile file                      ',/,&
'@  check the name of the meteo file                          ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA (atlecm)      ',/,&
'@     =======                                                ',/,&
'@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,&
'@                                                            ',/,&
'@  Error opening the meteo profile file                      ',/,&
'@    The meteo profile file has been opened but its content  ',/,&
'@    is incomplete or under a wrong format                   ',/,&
'@    check the format of the file (see the user guide):      ',/,&
'@                                                            ',/,&
'@    year (integer), quantile (integer), hour (integer),     ',/,&
'@          minute (integer), second (dble prec) of the profile',/,&
'@    location of the meteo profile (x,y) (dble prec)         ',/,&
'@    sea level pressure (double precision)                   ',/,&
'@ temperature profile:                                       ',/,&
'@   number of altitudes (integer)                            ',/,&
'@   alt.,temperature  in celcius,humidity in kg/kg (dble prec)',/,&
'@ wind profile:                                              ',/,&
'@   number of altitudes (integer)                            ',/,&
'@   alt.,u,v,k,eps (double precision)                        ',/,&
'@ NO LINE AT THE END OF THE FILE                             ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

end subroutine
