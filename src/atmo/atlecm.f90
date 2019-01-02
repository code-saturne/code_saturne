!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
!> \file atlecm.f90
!> \brief Reads the meteo profile data for the atmospheric module
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   imode       IMODE = 0 : READING FOR DIMENSIONS ONLY
!>                          IMODE = 1 : READING ACTUAL METEO DATA
!-------------------------------------------------------------------------------
subroutine atlecm ( imode )

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

! a function used in the routine
! for the diagnostic of liquid water content
double precision qsatliq

! Local variables

integer itp, ii, ios, k
integer year, quant,hour,minute, month, day
integer ih2o

double precision second
double precision sjday, jday
double precision rap,rscp,tmoy, rhmoy
double precision ztop, zzmax, tlkelv, pptop, dum
double precision rhum,q0,q1

character(len=80) :: ccomnt,oneline
character(len=1)  :: csaute

!===============================================================================

if (imode.eq.0) then
  write(nfecra, *) 'Reading dimensions for meteo profiles'
else
  write(nfecra, *) 'Reading meteo profiles data'
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
ih2o = 0

if (imode.eq.1) then
  rscp=rair/cp0

  !--> flag to take into account the humidity
  if (ippmod(iatmos).eq.2) ih2o=1

endif

!===============================================================================
! 1. Loop on time
!===============================================================================

100 continue
itp = itp+1

! ---> Read the comments
! ----------------------

101 read(impmet, '(a80)', err=999, end=906) ccomnt

if(ccomnt(1:1).eq.csaute) go to 101
backspace(impmet)

!===============================================================================
! 2. Read the time of the profile
!===============================================================================

! two formats possible :
! --> year, month, day, hour, minute, second  of the profile (UTC)
! --> year, quant-day, hour, minute, second  of the profile (UTC)
! NB: second is real, all others are integers

if (imode.eq.0) then
  read(impmet, *, err=999, end=906)
else
  second = -9999.d0
  read(impmet, '(a80)', err=999, end=906) oneline
  read(oneline, *, err=907, end=907) year, month, day,  &
                                                  hour, minute, second
! --> catch some read errors
  if (month.gt.12.or.day.gt.31) then
    write(nfecra,8005)
    call csexit (1)
  endif
  call comp_quantile(day, month, year, quant)
  goto 908
907  continue
  read(oneline, *, err=999, end=906) year, quant, hour, minute, second
908 continue
! --> catch some read errors
  if (second.lt.0d0.or.quant.gt.366) then
    write(nfecra,8005)
    call csexit (1)
  endif

  ! --> if the date and time are not completed in usati1.f90,
  !     the date and time of the first meteo profile are taken as the
  !     starting time of the simulation

  if (syear.lt.0) then
    syear = year
    squant = quant
    shour = hour
    smin = minute
    ssec = second
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

  tmmet(itp) = (jday - sjday)*86400.d0 + (hour - shour)*3600.d0     &
       + (minute - smin)*60.d0  + (second - ssec)

  ! --> check the chronological order of profiles

  if (itp.gt.1) then
    if (tmmet(itp).lt.tmmet(itp-1)) then
      write(nfecra, 8000)
      call csexit (1)
      !==========
    endif
  endif

endif

!===============================================================================
! 3. Read the position of the profile
!===============================================================================

102 read(impmet, '(a80)', err=999, end=999) ccomnt

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

103 read(impmet, '(a80)', err=999, end=999) ccomnt

if (ccomnt(1:1).eq.csaute) go to 103
backspace(impmet)

if (imode.eq.0) then
  read(impmet, *, err=999, end=999)
else
  read(impmet, *, err=999, end=999) pmer(itp)
endif

!===============================================================================
! 5. Read the temperature and humidity profiles
!===============================================================================

104 read(impmet, '(a80)', err=999, end=999) ccomnt

if (ccomnt(1:1).eq.csaute) go to 104
backspace(impmet)

if (imode.eq.0) then
  read(impmet, *, err=999, end=999) nbmett

  if (nbmett.le.1) then
    write(nfecra, 8001)
    call csexit (1)
    !==========
  endif

  do ii = 1, nbmett
    read (impmet,*,err=999,end=999) zzmax
  enddo

  !-->  Computes nbmaxt:
  !     if 1d radiative model (iatra1 = 1) altitudes are completed up to 11000m
  !     (i.e. nbmaxt > nbmett, used for dimensions of ttmet and qvmet array)
  !     if no radiative model nbmaxt = nbmett
  nbmaxt = nbmett

  if (iatra1.eq.1) then
    ztop = 11000.d0
    zzmax = (int(zzmax)/1000)*1000.d0
    do while(zzmax.le.(ztop-1000.d0))
      zzmax = zzmax + 1000.d0
      nbmaxt = nbmaxt + 1
    enddo
  endif

else

  read(impmet, *, err=999, end=999)

  do ii = 1, nbmett

    !     Altitude, temperature, humidite
    if (ippmod(iatmos).eq.2) then
      read (impmet,*,err=999,end=999) ztmet(ii),                   &
           ttmet(ii,itp),qvmet(ii,itp),&
           ncmet(ii,itp)
    else
      read (impmet,*,err=999,end=999) ztmet(ii),                   &
           ttmet(ii,itp),qvmet(ii,itp)
    endif

    !--> Check the unity of the specific humidity (kg/kg) when used

    if (qvmet(ii,itp).gt.0.1 .or. qvmet(ii,itp).lt.0.) then
      write(nfecra,8003)
      call csexit (1)
    endif
    if (ippmod(iatmos).eq.2) then
      if (ncmet(ii,itp).lt.0.) then
        write(nfecra,8004)
        call csexit (1)
      endif
    endif

  enddo

  !--> If 1D radiative model (iatra1 = 1), complete the temperature and humidity
  !    profiles up to 11000m

  if (iatra1.eq.1) then
    ztop = 11000.d0
    ii = nbmett
    zzmax = (int(ztmet(ii))/1000)*1000.d0
    ii = nbmett
    do while(zzmax.le.(ztop-1000.d0))
      zzmax = zzmax + 1000.d0
      ii = ii + 1
      ztmet(ii) = zzmax
      call atmstd(ztmet(ii), dum, tlkelv, dum) ! standard temperature profile
      ttmet(ii,itp) = tlkelv - tkelvi
      if (iqv0.eq.0) then
        qvmet(ii,itp) = 0.d0   ! standard atmosphere: q=0
      else                     ! use a decreasing exponential
                               ! profile to complete qvmet
        qvmet(ii,itp) = qvmet(nbmett,itp)                       &
             *exp((ztmet(nbmett)-ztmet(ii))/2.5d3)
      endif
      ncmet(ii,itp) = 0.d0
    enddo

  endif

endif

!===============================================================================
! 6. Compute hydro pressure profile  (Laplace integration)
!===============================================================================
! If ihpm = 0 (default): bottom to top Laplace integration based on pressure at
! sea level (pmer(itp))
! If ihpm = 1 (usati1): top to bottom Laplace integration based on pressure at
! the top of the domain (ztmet(nbmaxt)) for the standard atmosphere

if (imode.eq.1) then

  phmet(1, itp) = pmer(itp)
  rscp = rair/cp0

  if (ihpm.eq.0) then
    phmet(1,itp) = pmer(itp)
    do k = 2, nbmaxt
      tmoy = 0.5d0*(ttmet(k-1,itp) + ttmet(k,itp)) + tkelvi

      if(ippmod(iatmos).eq.2) then ! take liquid water into account
        q0 = min( qvmet(k-1,itp), qsatliq( ttmet(k-1,itp) &
            + tkelvi, phmet(k-1,itp)))
        q1 = min( qvmet(k  ,itp), qsatliq( ttmet(k  ,itp) &
            + tkelvi, phmet(k-1,itp)))
        !in q1=.. phmet(k-1,itp) is not a mistake: we can not use phmet(k,itp)
        !since this is what we want to estimate.
      else
        q0 = qvmet(k-1,itp)
        q1 = qvmet(k  ,itp)
      endif

      rhmoy = rair*(1.d0 + (rvsra-1.d0)*                   &
           (q0 + q1)/2.d0*ih2o)
      rap = -abs(gz)*(ztmet(k)-ztmet(k-1))/rhmoy/tmoy
      phmet(k,itp) = phmet(k-1,itp)*exp(rap)
    enddo
  else
    call atmstd (ztmet(nbmaxt), pptop, dum, dum) ! standard profile temperature
    phmet(nbmaxt,itp) = pptop
    do k = nbmaxt-1, 1, -1
      tmoy = 0.5d0*(ttmet(k+1,itp) + ttmet(k,itp)) + tkelvi

      if(ippmod(iatmos).eq.2) then ! take liquid water into account
        q0 = min( qvmet(k,itp), qsatliq( ttmet(k,itp)+tkelvi, phmet(k+1,itp)))
        q1 = min( qvmet(k+1  ,itp), qsatliq( ttmet(k+1,itp)+tkelvi, phmet(k+1,itp)))
        !in q0=.. phmet(k+1,itp) is not a mistake: we can not use phmet(k,itp)
        !since this is what we want to estimate.
      else
        q0=qvmet(k  ,itp)
        q1=qvmet(k+1,itp)
      endif

      rhmoy = rair*(1.d0+(rvsra - 1.d0)*                   &
              (q0 + q1)/2.d0*ih2o)
      rap = abs(gz)*(ztmet(k+1) - ztmet(k))/rhmoy/tmoy
      phmet(k,itp) = phmet(k+1,itp)*exp(rap)
    enddo
  endif

endif

!==============================================================================
! 7. Compute the pot. temperature profile and the density profile
!==============================================================================

if (imode.eq.1) then
  do k = 1, nbmaxt

    rhum = rair*(1.d0+(rvsra-1.d0)*qvmet(k,itp)*ih2o)

    if (ippmod(iatmos).eq.0) then
      !constant density
      rmet(k,itp) = phmet(1,itp)/(ttmet(k,itp) + tkelvi)/rhum
    else
      !variable density
      rmet(k,itp) = phmet(k,itp)/(ttmet(k,itp) + tkelvi)/rhum
    endif
    rscp = (rair/cp0)*(1.d0 + (rvsra-cpvcpa)*qvmet(k,itp)*ih2o)
    tpmet(k,itp) = (ttmet(k,itp)+tkelvi)*((ps/phmet(k,itp))**rscp)
  enddo

endif

!==============================================================================
! 8. Read the velocity profile
!==============================================================================

105 read(impmet, '(a80)', err=999, end=999) ccomnt

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
    read (impmet, *, err=999, end=999) zdmet(ii),                    &
         umet(ii, itp),  vmet(ii, itp),  &
         ekmet(ii, itp), epmet(ii, itp)
  enddo

endif

!===============================================================================
! 9. Printings
!===============================================================================

if (imode.eq.1) then
  if (itp.eq.1) then
    write(nfecra, *)
    write(nfecra, *) '==================================================='
    write(nfecra, *) 'printing meteo profiles'
  endif
  write(nfecra, *) 'year, quant-day, hour, minute, second:'
  write(nfecra, 7995) year, quant, hour, minute, second
7995 format(1x, i4, i5, 2i4, f10.2)
  write(nfecra, *) 'tmmet(itp):'
  write(nfecra, 7996) tmmet(itp)
7996 format(1x, f10.2)
  write(nfecra, *) 'zdmet, umet, vmet, ekmet, epmet:'
  do ii = 1, nbmetd
    write(nfecra, 7997)                                            &
         zdmet(ii), umet(ii, itp), vmet(ii, itp), ekmet(ii, itp), epmet(ii, itp)
7997 format(1x, 3f8.2, 2e10.3)
  enddo
  if (ippmod(iatmos).eq.0) then
    write(nfecra, *) '==================================================='
    write(nfecra, *) 'WARNING : option  constant density                 '
    write(nfecra, *) 'WARNING : thermal profile will be ignored          '
    write(nfecra, *) '==================================================='
  endif
  if (ippmod(iatmos).le.1) then
    write(nfecra,*) 'ztmet,ttmet,tpmet,rmet,phmet,qvmet:'
    do ii = 1, nbmaxt
      write(nfecra,7998)                                            &
           ztmet(ii), ttmet(ii,itp), tpmet(ii,itp),              &
           rmet(ii,itp), phmet(ii,itp), qvmet(ii,itp)
7998  format(1x, 3f8.2,f8.4,f12.3,e10.3)
    enddo
  else
    write(nfecra,*) 'ztmet,ttmet,tpmet,rmet,phmet,qvmet,qsat,nc:'
    do ii = 1, nbmaxt
      write(nfecra,7999)                                            &
           ztmet(ii), ttmet(ii,itp), tpmet(ii,itp),              &
           rmet(ii,itp), phmet(ii,itp), qvmet(ii,itp),           &
           qsatliq(ttmet(ii,itp)+tkelvi , phmet(ii,itp)),         &
           ncmet(ii,itp)
7999  format(1x, 3f8.2,f8.4,f12.3,e10.3,e10.3,e12.5)
    enddo
  endif

endif

!===============================================================================
! 10. End of the loop on time
!===============================================================================

goto 100

906 continue

if (imode.eq.0) nbmetm= itp-1

close(unit=impmet)

! ---
! End
! ---

return

!============================
! XX. Error outputs
!============================

99 continue
write ( nfecra, 9998 )
call csexit (1)
!==========

999 continue
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
8003 format (                                                    &
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                            ',/,&
     '@ @@ATTENTION:VERIFICIATION A L''ENTREE DES DONNEES (atlecm)',/,&
     '@   =========                                               ',/,&
     '@        PHYSIQUE PARTICULIERE ATMOSPHERIQUE                 ',/,&
     '@                                                            ',/,&
     '@     erreur a la lecture du fichier meteo                   ',/,&
     '@  humidites specifiques irrealistes                         ',/,&
     '@  verifier les unites (kg/kg)                               ',/,&
     '@                                                            ',/,&
     '@  Le calcul ne sera pas execute.                            ',/,&
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                             ',/)
8004 format (                                                    &
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                            ',/,&
     '@ @@ATTENTION:VERIFICATION A L''ENTREE DES DONNEES (atlecm)' ,/,&
     '@   =========                                               ',/,&
     '@        PHYSIQUE PARTICULIERE ATMOSPHERIQUE                 ',/,&
     '@                                                            ',/,&
     '@  erreur a la lecture du fichier meteo                      ',/,&
     '@  nombre de gouttes lus <0                                  ',/,&
     '@  verifier le fichier meteo                                 ',/,&
     '@                                                            ',/,&
     '@  Le calcul ne sera pas execute.                            ',/,&
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                             ',/)
8005 format (                                                    &
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                            ',/,&
     '@ @@ATTENTION:VERIFICATION A L''ENTREE DES DONNEES (atlecm)  ',/,&
     '@   =========                                                ',/,&
     '@        PHYSIQUE PARTICULIERE ATMOSPHERIQUE                 ',/,&
     '@                                                            ',/,&
     '@  erreur a la lecture de la date du fichier meteo           ',/,&
     '@  verifier le format (entiers, reel)                        ',/,&
     '@                                                            ',/,&
     '@  Le calcul ne sera pas execute.                            ',/,&
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                             ',/)
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
8003 format (                                                    &
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                            ',/,&
     '@ @@  WARNING:   CHECKING INPUT DATA (atlecm)                ',/,&
     '@     =======                                                ',/,&
     '@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,&
     '@                                                            ',/,&
     '@              Error in the meteo profile file:              ',/,&
     '@  the values for the specific humidity are not realistic    ',/,&
     '@  check the unity (kg/kg)                                   ',/,&
     '@                                                            ',/,&
     '@  The computation will not be run                           ',/,&
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                             ',/)
8004 format (                                                    &
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                            ',/,&
     '@ @@  WARNING:   CHECKING INPUT DATA (atlecm)                ',/,&
     '@     =======                                                ',/,&
     '@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,&
     '@                                                            ',/,&
     '@              Error in the meteo profile file:              ',/,&
     '@  Number of droplets read  <  0                             ',/,&
     '@  The computation will not be run                           ',/,&
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                             ',/)
8005 format (                                                    &
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                            ',/,&
     '@ @@  WARNING:   CHECKING INPUT DATA (atlecm)                ',/,&
     '@     =======                                                ',/,&
     '@      ATMOSPHERIC SPECIFIC PHYSICS                          ',/,&
     '@                                                            ',/,&
     '@              Error in the date of meteo profile file:      ',/,&
     '@  Check the format (integers,real) for the date             ',/,&
     '@  The computation will not be run                           ',/,&
     '@                                                            ',/,&
     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
     '@                                                             ',/)
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

end subroutine atlecm


!===============================================================================
!> \brief Compute the calendar day number from the date
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]       jour        day
!> \param[in]       mois        month
!> \param[in]       annee       year
!> \param[out]      quant       calendar day number
!-------------------------------------------------------------------------------
subroutine comp_quantile(jour,mois,annee,quant)

implicit none

! Arguments
integer jour,mois,annee
integer quant

! Local variables
integer distrib, booll, retrait

distrib = int(mois * 275 / 9.) - 30
booll = int((mois + 9 ) / 12.) !boolean=0 for jan and feb and =1 for march to dec
retrait = 1+int(mod(annee,4)+2)/3 !retrait=2 days non leap years and =1 day for leap years
quant = distrib-(booll*retrait)+jour

return
end subroutine comp_quantile


