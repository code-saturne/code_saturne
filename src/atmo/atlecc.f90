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

!===============================================================================
!  Purpose:
!  --------

!> \file atlecc.f90
!> \brief Reads the chemistry profile data for the atmospheric chemistry
!
!> \brief Reads the chemistry profile data for the atmospheric chemistry
!>-     imode = 0 : reading of dimensions only, imode = 1 : reading of data
!>-     warning : arguments are defined only on second call
!>             of the routine (imode = 1)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     imode         execution mode
!_______________________________________________________________________________

subroutine atlecc (imode )

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use field
use cstnum
use cstphy
use ppppar
use atincl
use numvar
use atchem

implicit none

!===============================================================================

! Arguments

integer           imode

! Local variables

integer f_id
integer itp, ii, ios, k
integer sjday,minute
double precision second
integer year, month, quant, hour, day, jday
character(len=80) :: ccomnt, label, oneline
character(len=1) :: csaute

! altitudes and concentrations of every nespgi species
double precision  zconctemp(nespgi+1)
! names of species defined in the file in case of a user defined chemical scheme
character(len=80), allocatable, dimension(:) ::      namespg

!===============================================================================

if (imode.eq.0) then
  write(nfecra,*) ''
  write(nfecra,*) 'reading of dimensions for concentration profiles'
else
  write(nfecra,*) ''
  write(nfecra,*) 'reading of concentration profiles data'
 endif

!===============================================================================
! 0. initialization
!===============================================================================

csaute = '/'
! --> Opens the concentration profiles file
open ( unit=impmec, file=ficmec, status='old', iostat=ios, err=99 )
rewind ( unit=impmec,err=99 )

itp = 0

!===============================================================================
! 1. loop on time
!===============================================================================

 100  continue
itp = itp+1

! ---> reading the comments
! ---------------------------

 101  read (impmec,'(A80)',err=999,end=906) ccomnt

if (ccomnt(1:1).eq.csaute) go to 101
backspace(impmec)

!===============================================================================
! 2. reading the time of the profile
!===============================================================================

! two formats possible :
! --> year, month, day, hour, minute, second  of the profile (UTC)
! --> year, quant-day, hour, minute, second  of the profile (UTC)
! NB: second is real, all others are integers

if (imode.eq.0) then
  read(impmec, *, err=999, end=906)
else
  second = -9999.d0
  read(impmec, '(a80)', err=999, end=906) oneline
  read(oneline,*,err=907,end=907) year, month, day, hour, minute, second
  ! --> catch some read errors
  if (month.gt.12.or.day.gt.31) then
    write(nfecra,8005)
    call csexit (1)
  endif
  call comp_quantile(day, month, year, quant)
  goto 908
907  continue
  read(oneline,*,err=999,end=906) year, quant, hour, minute, second
908 continue
! --> catch some read errors
  if (second.lt.0d0.or.quant.gt.366) then
    write(nfecra,8005)
    call csexit (1)
  endif

  ! --> if the date and time are not completed in usati1.f90 nor and if no meteo
  ! --> file is given,
  ! --> the date and time of the first concentration profile are taken as the
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
  sjday = squant + ((1461 * (syear + 4800 + (1 - 14) / 12)) / 4 +  &
             (367 * (1 - 2 - 12 * ((1 - 14) / 12))) / 12 -            &
             (3 * ((syear + 4900 + (1 - 14) / 12) / 100)) / 4         &
              + 1 - 32075) - 1

! --> Compute the julian day for the date of the current profile
!     (julian day at 12h)
  jday = quant + ((1461 * (year + 4800 + (1 - 14) / 12)) / 4 +   &
             (367 * (1 - 2 - 12 * ((1 - 14) / 12))) / 12 -           &
            (3 * ((year + 4900 + (1 - 14) / 12) / 100)) / 4          &
              + 1 - 32075) - 1

  tchem(itp) = (jday - sjday)*86400 + (hour - shour)*3600.d0     &
       + (minute - smin)*60.d0  + (second - ssec)


! --> verifying the chronological order of profiles

  if(itp.gt.1) then
    if(tchem(itp).lt.tchem(itp-1)) then
      write(nfecra,8000)
      call csexit (1)
    endif
  endif
endif

!===============================================================================
! 3. reading the position of the profile
!===============================================================================

 102  read (impmec,'(A80)',err=999,end=999) ccomnt

if(ccomnt(1:1).eq.csaute) go to 102
backspace(impmec)

if (imode.eq.0) then
  read(impmec,*,err=999,end=999)
else
  read(impmec,*,err=999,end=999) xchem(itp), ychem(itp)
endif

!===============================================================================
! 4. reading the number of species initialized by the chemistry file
!===============================================================================

 107  read (impmec,'(a80)',err=999,end=999) ccomnt

if (ccomnt(1:1).eq.csaute) go to 107
backspace(impmec)


if (imode.eq.0) then
  read (impmec,*,err=999,end=999) nespgi
  if (nespgi.gt.size(isca_chem)) then
    write(nfecra,8002) size(isca_chem), nespgi
    call csexit (1)
    !==========
  endif
else
  read (impmec,*,err=999,end=999)
endif

!===============================================================================
! 5. reading the species initialized by the chemistry file
!===============================================================================
if (nespgi.ge.1) then

 108  read(impmec,'(a80)',err=999,end=999) ccomnt

  if(ccomnt(1:1).eq.csaute) go to 108
  backspace(impmec)


  if (imode.eq.0) then
    read(impmec,*,err=999,end=999)
  else
    read(impmec,*,err=999,end=999) idespgi(1:nespgi)
  endif

!===============================================================================
! 6. reading the concentration profiles
!===============================================================================

109 read(impmec,'(a80)',err=999,end=999) ccomnt

  if(ccomnt(1:1).eq.csaute) go to 109
  backspace(impmec)

  if (imode.eq.0) then
    read(impmec,*,err=999,end=999) nbchmz
    if(nbchmz.le.1) then
      write(nfecra,8001)
      call csexit (1)
    endif

    do ii = 1, nbchmz
      read (impmec,*,err=999,end=999)
    enddo

  else

    read(impmec,*,err=999,end=999)

    do ii = 1, nbchmz

      !     Altitudes and concentrations of every nespgi scpecies
      read (impmec,*,err=999,end=999) zconctemp

      zproc(ii) = zconctemp(1)
      do k = 2, nespgi+1
        espnum(ii+(itp-1)*nbchmz+(k-2)*nbchmz*nbchim) = zconctemp(k)
      enddo

    enddo
  endif

endif ! fin test nespgi

!===============================================================================
! 7. logging
!===============================================================================

if (imode.eq.1) then
  if (itp.eq.1) then
    write(nfecra, *)
    write(nfecra, *) '==================================================='
    write(nfecra, *) 'printing concentration profiles'
  endif
  write(nfecra, *) 'year, quant-day, hour, minute, second'
  write(nfecra, 7995) year, quant, hour, minute, second
7995 format(1x, i4, i5, 2i4, f10.2)
  write(nfecra, *) 'tchem(itp)'
  write(nfecra, 7996) tchem(itp)
7996 format(1x, f10.2)
  write(nfecra, '(a)', advance='no') 'zproc, '
  do ii = 1, nespgi
    f_id = ivarfl(isca(isca_chem(ii)))
    call field_get_label(f_id, label)
    if (ii .lt. nespgi) then
      write(nfecra, '(a)', advance='no') trim(label)//', '
    else
      write(nfecra, '(a)') trim(label)
    endif
  enddo
  do ii = 1, nbchmz
    write(nfecra, 7797) zproc(ii),                                         &
    (espnum(ii+(itp-1)*nbchmz+(k-1)*nbchmz*nbchim), k = 1, nespgi)
7797 format(1x,f8.2,1x,10(es10.2))
  enddo
endif

!================================================================================
! 11. end of the loop on time
!================================================================================

goto 100

906 continue

if (imode.eq.0) nbchim = itp-1
close(unit=impmec)

! ---
! End
! ---

return

!============================
! 12. error outputs
!============================

 99   continue
write ( nfecra,9998 )
call csexit (1)
!==========

 999  continue
write ( nfecra,9999 )
call csexit (1)
!==========

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)
 8000 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (atlecc)      ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@    MODULE DE CHIMIE (ICHEMISTRY) DEMANDE                   ',/,&
'@                                                            ',/,&
'@              Erreur  dans le fichier chimie :              ',/,&
'@      ordre chronologique des profils non respecte          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8001 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (atlecc)      ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@    MODULE DE CHIMIE (ICHEMISTRY) DEMANDE                   ',/,&
'@                                                            ',/,&
'@              Erreur  dans le fichier chimie :              ',/,&
'@  le nombre de mesures de concentrations doit etre          ',/,&
'@  superieur a 2                                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (atlecc)      ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@    MODULE DE CHIMIE (ICHEMISTRY) DEMANDE                   ',/,&
'@                                                            ',/,&
'@  Le nombre d''especes a initialiser avec le fichier chimie ',/,&
'@  est superieur au nombre de scalaires chimie               ',/,&
'@                                                            ',/,&
'@   Nombre de scalaires chimie       declares : ',I10         ,/,&
'@   Nombre d''especes chimiques a initialiser : ',I10         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8005 format (                                                    &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ATTENTION:VERIFICATION A L''ENTREE DES DONNEES (atlecc)  ',/,&
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
 9998 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (atlecc)      ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@    MODULE DE CHIMIE (ICHEMISTRY) DEMANDE                   ',/,&
'@                                                            ',/,&
'@  Erreur a l''ouverture du fichier chimie                   ',/,&
'@  Verifier le nom du fichier chimie                         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (atlecc)      ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@    MODULE DE CHIMIE (ICHEMISTRY) DEMANDE                   ',/,&
'@                                                            ',/,&
'@  Erreur a la lecture du fichier chimie.                    ',/,&
'@    Le fichier a ete ouvert mais est peut etre incomplet    ',/,&
'@    ou son format inadapte.                                 ',/,&
'@                                                            ',/,&
'@    year (integer), month (integer), day (integer),         ',/,&
'@    hour (integer), minute (integer), second (dble prec)    ',/,&
'@    of the profile                                          ',/,&
'@    OU                                                      ',/,&
'@    year (integer), quantile (integer), hour (integer),     ',/,&
'@    minute (integer), second (dble prec) of the profile     ',/,&
'@    location of the concentration profile (x,y) (dble prec) ',/,&
'@    number of reactions (integer)                           ',/,&
'@    number of species (integer)                             ',/,&
'@    names of the species                                    ',/,&
'@    molar mass of the species                               ',/,&
'@ concentration profiles:                                    ',/,&
'@   number of altitudes (integer)                            ',/,&
'@   alt, concentration for each species                      ',/,&
'@ NO LINE AT THE END OF THE FILE                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 8000 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA (atlecc)      ',/,&
'@     =======                                                ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@              Error in the chemistry profile file:             ',/,&
'@      check that the chronogical order of the profiles      ',/,&
'@      are respected                                         ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8001 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA (atlecc)      ',/,&
'@     =======                                                ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@              Error in the chemistry profile file:             ',/,&
'@  the number of concentrations measurements must be larger  ',/,&
'@  than 2                                                    ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP WHILE READING INPUT DATA (atlecc)        ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@  The number of species to initialize with the chemistry file  ',/,&
'@  is larger than the number of chemistry model scalars      ',/,&
'@                                                            ',/,&
'@   Number of chemistry model scalars: ',I10                  ,/,&
'@   Number of species to initialize: ',I10                    ,/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8005 format (                                                    &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP WHILE READING INPUT DATA (atlecc)        ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@  Error opening the chemistry profile file                  ',/,&
'@  verify input format (integer, real)                       ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                             ',/)
 9998 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA (atlecc)      ',/,&
'@     =======                                                ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@  Error opening the chemistry profile file                  ',/,&
'@  check the name of the chemistry file                      ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA (atlecc)      ',/,&
'@     =======                                                ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@  Error opening the chemistry profile file                     ',/,&
'@    The chemistry profile file has been opened but its content ',/,&
'@    is incomplete or under a wrong format                   ',/,&
'@    check the format of the file (see the user guide):      ',/,&
'@                                                            ',/,&
'@    year (integer), month (integer), day (integer),         ',/,&
'@    hour (integer), minute (integer), second (dble prec)    ',/,&
'@    of the profile                                          ',/,&
'@    OR                                                      ',/,&
'@    year (integer), quantile (integer), hour (integer),     ',/,&
'@    minute (integer), second (dble prec) of the profile     ',/,&
'@    location of the concentration profile (x,y) (dble prec) ',/,&
'@    number of reactions (integer)                           ',/,&
'@    number of species (integer)                             ',/,&
'@    names of the species                                    ',/,&
'@    molar mass of the species                               ',/,&
'@ concentration profiles:                                    ',/,&
'@   number of altitudes (integer)                            ',/,&
'@   alt, concentration for each species                      ',/,&
'@ NO LINE AT THE END OF THE FILE                             ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

end subroutine atlecc
