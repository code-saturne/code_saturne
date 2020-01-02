!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine strini &
!================

 ( dt     )

!===============================================================================
! FONCTION :
! ----------

! INITILISATION DES DONNEES DES STRUCTURES MOBILES EN ALE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
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
use optcal
use cstphy
use entsor
use pointe
use albase
use alstru
use alaste
use parall
use period
use mesh
use post
use field

!===============================================================================

implicit none

! Arguments

double precision dt(ncelet)

! Local variables

integer          n_fields, f_id, flag
integer          ifac  , istr, icompt, ii
integer          mbstru, mbaste

integer          inod
integer          indast

integer, allocatable, dimension(:) :: itrav
integer, allocatable, dimension(:) :: lstfac, idfloc, idnloc

!===============================================================================

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_ast_coupling_initialize(nalimx, epalim) &
    bind(C, name='cs_ast_coupling_initialize')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nalimx
    real(kind=c_double), value :: epalim
  end subroutine cs_ast_coupling_initialize

end interface

!===============================================================================
! 1. INITIALISATION
!===============================================================================

do istr = 1, nstrmx
  dtstr(istr) = dt(1)
enddo

do istr = 1, nastmx
  do ii = 1, 3
    asddlf(ii,istr) = 1
    asddlc(ii,istr) = 1
  enddo
enddo

!     NBSTRU et NBASTE valent -999 si le calcul n'est pas un calcul
!       suite ou s'il est une suite d'un calcul sans ALE
mbstru = nbstru
mbaste = nbaste

!     En ALE on met IHISTR a 1 par defaut
!       (remis a zero sans structure)
ihistr = 1

!===============================================================================
! 2.  RESERVATION DU TABLEAU IDFSTR
!===============================================================================

do ifac = 1, nfabor
  idfstr(ifac) = 0
enddo

!===============================================================================
! 2.  User definition of idfstr
!===============================================================================

! Internal structures

call uistr1 &
( idfstr, mbstru,          &
  aexxst, bexxst, cfopre,  &
  ihistr,                  &
  xstp, xstreq, xpstr )

call usstr1                                                       &
 ( idfstr ,                                                       &
   aexxst , bexxst , cfopre ,                                     &
   xstp   , xpstr  , xstreq )

! External structures: Code_Saturne / code_aster coupling

call uiaste(idfstr, asddlf)
call usaste(idfstr)

!===============================================================================
! 3.  CALCUL DE NBSTRU ET NBASTE
!===============================================================================

! 3.1 STRUCTURES INTERNES :
! -----------------------

nbstru = 0
do ifac = 1, nfabor
  if (idfstr(ifac).gt.nbstru) nbstru = idfstr(ifac)
enddo

if (irangp.ge.0) call parcmx(nbstru)
                 !==========

if (nbstru.gt.nstrmx) then
  write(nfecra,4000)
  call csexit(1)
endif

!     On compare NBSTRU a la valeur eventuelle anterieure

if (mbstru.gt.-999) then
  if (nbstru.ne.mbstru) then
    write(nfecra,4001)mbstru,nbstru
    call csexit(1)
  endif
endif

! 3.2 STRUCTURES EXTERNES : COUPLAGE CODE_SATURNE / CODE_ASTER
! -----------------------

nbaste = 0
do ifac = 1, nfabor
  if (-idfstr(ifac).gt.nbaste) nbaste = -idfstr(ifac)
enddo

if (irangp.ge.0) call parcmx(nbaste)
                 !==========

if (nbaste.gt.nastmx) then
  write(nfecra,4002)
  call csexit(1)
endif

!     On compare NBASTE a la valeur eventuelle anterieure

if (mbaste.gt.-999) then
  if (nbaste.ne.mbaste) then
    write(nfecra,4003)mbaste,nbaste
    call csexit(1)
  endif
endif


!===============================================================================
! 5.  CALCUL ET ENVOI A CODE_ASTER DES PARAMETRES GEOMETRIQUES
!===============================================================================

if (nbaste.gt.0) then

  ! Allocate a work array
  allocate(itrav(nnod))

  do inod = 1, nnod
     itrav(inod) = 0
  enddo

  nbfast = 0
  nbnast = 0

!       Calcul du nombre de faces et noeuds couples avec code_aster
  do ifac = 1, nfabor
    istr = idfstr(ifac)
    if (istr.lt.0) then
      nbfast = nbfast + 1
      do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
        inod = nodfbr(ii)
        itrav(inod) = istr
      enddo
    endif
  enddo
  do inod = 1, nnod
    if (itrav(inod).lt.0) nbnast = nbnast + 1
  enddo

  ! Allocate temporary arrays
  allocate(lstfac(nbfast))
  allocate(idfloc(nbfast), idnloc(nbnast))

  indast = 0
  do ifac = 1, nfabor
    istr = idfstr(ifac)
    if (istr.lt.0) then
      indast = indast + 1
      lstfac(indast) = ifac
      idfloc(indast) = -istr
    endif
  enddo
  nbfast = indast

  indast = 0
  do inod = 1, nnod
    istr = itrav(inod)
    if (istr.lt.0) then
      indast = indast + 1
      idnloc(indast) = -istr
    endif
  enddo
  nbnast = indast

  ! Exchange code_aster coupling parameters
  call cs_ast_coupling_initialize(nalimx, epalim)

  ! Send geometric information to code_aster
  call astgeo(nbfast, lstfac, idfloc, idnloc, almax)

  ! Free memory
  deallocate(lstfac)
  deallocate(idfloc, idnloc)

endif


!===============================================================================
! 6.  MESSAGES D'INFORMATION SUR LE COUPLAGE
!===============================================================================

!     Valeur par defaut et verifiction de IHISTR
if (nbstru.eq.0) ihistr = 0

call field_get_n_fields(n_fields)

icompt = 0
do f_id = 0, n_fields - 1
  call field_get_key_int(f_id, keyvis, flag)
  if (iand(flag, POST_MONITOR).ne.0) icompt = icompt+1
enddo
if(icompt.eq.0 .and. ihistr.eq.0) then
  nthist = -1
  frhist = -1.d0
endif

if (ihistr.ne.0 .and. ihistr.ne.1) then
  write(nfecra,1000)ihistr
  call csexit(1)
endif

!     Si NBSTRU=0 et NBASTE=0, on desalloue IDFSTR et on passe NALIMX a 1
!       si necessaire
if (nbstru.gt.0) then
  write(nfecra,2010) nbstru,alpnmk,betnmk,gamnmk,ihistr
else
  write(nfecra,2000) nbstru
endif
if (nbaste.gt.0) then
  write(nfecra,2012) nbaste
else
  write(nfecra,2002) nbaste
endif
if (nbstru.eq.0.and.nbaste.eq.0) then
  if (nalimx.gt.1) then
    write(nfecra,2001)
    nalimx = 1
  endif
else if (nbstru.gt.0) then
  if (nalimx.eq.1) then
    write(nfecra,2020) aexxst, bexxst, cfopre
  else
    cfopre = 1.d0
    write(nfecra,2030) nalimx, epalim
  endif
endif
write(nfecra,3000)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR D''ECRITURE DES FICHIERS HISTORIQUES DES  ',/,&
'@      STRUCTURES MOBILES NE PEUT VALOIR QUE 0 OU 1.         ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes dans usstru.               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format( &
    /,'MODE COUPLAGE DE STRUCTURES NON ACTIVE  ',/,&
      '                 NBSTRU = ',I10                         ,/)
 2001 format( &
      '               NALIMX INUTILE ET POSITIONNE A 1        ',/)
 2002 format( &
    /,'MODE COUPLAGE CODE_ASTER NON ACTIVE     ',/,&
      '                 NBASTE = ',I10                         ,/)
 2010 format( &
    /,'MODE COUPLAGE DE STRUCTURES ACTIVE      ',/,&
      '                 AVEC NBSTRU = ',I10   ,' STRUCTURE(S) ',/,&
      '                                                       ',/,&
      '               COEFFICIENTS DE NEWMARK :               ',/,&
      '                 ALPNMK = ',E12.4                       ,/,&
      '                 BETNMK = ',E12.4                       ,/,&
      '                 GAMNMK = ',E12.4                       ,/,&
      '                                                       ',/,&
      '               FICHIERS HISTORIQUES DES STRUCTURES :   ',/,&
      '                 IHISTR = ',I4,' ( 1 : active)         ',/)
 2012 format( &
    /,'MODE COUPLAGE CODE_ASTER ACTIVE         ',/,&
      '                 AVEC NBASTE = ',I10   ,' STRUCTURE(S) ',/)
 2020 format( &
    /,'SCHEMA DE COUPLAGE EXPLICITE ACTIVE     ',/,&
      '                                                       ',/,&
      '               COEFFICIENTS DU SCHEMA :                ',/,&
      '                 AEXXST = ',E12.4                       ,/,&
      '                 BEXXST = ',E12.4                       ,/,&
      '                 CFOPRE = ',E12.4                       ,/)
 2030 format( &
    /,'SCHEMA DE COUPLAGE IMPLICITE ACTIVE     ',/,&
      '                                                       ',/,&
      '               NB DE SOUS-ITERATIONS MAX. : ',I10       ,/,&
      '               SEUIL DE CONVERGENCE       : ',E12.4     ,/)

 3000 format( &
'-------------------------------------------------------------',/)
 4000 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES        ',/,&
'@                INTERNES                                    ',/,&
'@                                                            ',/,&
'@    Le nombre de structures definies est superieur au nombre',/,&
'@    maximum autorise NSTRMX :                               ',/,&
'@      Nombre de structures definies      : ',I10             ,/,&
'@      Nombre de structures autorisees    : ',I10             ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Reduire le nombre de structure                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4001 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES MOBILES',/,&
'@                INTERNES                                    ',/,&
'@                                                            ',/,&
'@    Le nombre de structures internes definies est           ',/,&
'@      different de celui du calcul precedent :              ',/,&
'@      Nombre de structures calcul precedent : ',I10          ,/,&
'@      Nombre de structures calcul actuel    : ',I10          ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire ou la specification',/,&
'@      des structures dans usstru.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4002 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES        ',/,&
'@                              EXTERNES (COUPLAGE CODE_ASTER)',/,&
'@                                                            ',/,&
'@    Le nombre de structures definies est superieur au nombre',/,&
'@    maximum autorise NASTMX :                               ',/,&
'@      Nombre de structures definies      : ',I10             ,/,&
'@      Nombre de structures autorisees    : ',I10             ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Reduire le nombre de structure                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4003 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES MOBILES',/,&
'@                             EXTERNES  (COUPLAGE CODE_ASTER)',/,&
'@                                                            ',/,&
'@    Le nombre de structures externes definies est           ',/,&
'@      different de celui du calcul precedent :              ',/,&
'@      Nombre de structures calcul precedent : ',I10          ,/,&
'@      Nombre de structures calcul actuel    : ',I10          ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire ou la specification',/,&
'@      des structures dans usaste.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE TIME MONITORING FILES INDICATOR FOR THE MOBILE      ',/,&
'@      STRUCTURES CAN ONLY TAKE THE VALUES 0 OR 1.           ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify the parameters given in usstru.                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format( &
    /,'COUPLING MODE FOR STRUCTURES NOT ACTIVATED ',/,&
      '              NBSTRU = ',I10                            ,/)
 2001 format( &
      '            NALIMX USELESS AND SET TO 1                ',/)
 2002 format( &
    /,'CODE_ASTER COUPLING MODE NOT ACTIVATED     ',/,&
      '              NBASTE = ',I10                            ,/)
 2010 format( &
    /,'COUPLING MODE FOR STRUCTURES ACTIVATED     ',/,&
      '              WITH NBSTRU = ',I10   ,' STRUCTURE(S)    ',/,&
      '                                                       ',/,&
      '            NEWMARK COEFFICIENTS:                      ',/,&
      '              ALPNMK = ',E12.4                          ,/,&
      '              BETNMK = ',E12.4                          ,/,&
      '              GAMNMK = ',E12.4                          ,/,&
      '                                                       ',/,&
      '            MONITORING FILES FOR STRUCTURES:           ',/,&
      '                 IHISTR = ',I4,' ( 1 : activated)      ',/)
 2012 format( &
    /,'CODE_ASTER COUPLING MODE ACTIVATED         ',/,&
      '              WITH NBASTE = ',I10   ,' STRUCTURE(S)    ',/)
 2020 format( &
    /,'EXPLICIT SCHEME FOR COUPLING ACTIVATED     ',/,&
      '                                                       ',/,&
      '            SCHEME COEFFICIENTS:                       ',/,&
      '              AEXXST = ',E12.4                          ,/,&
      '              BEXXST = ',E12.4                          ,/,&
      '              CFOPRE = ',E12.4                          ,/)
 2030 format(                                                           &
    /,'IMPLICIT SCHEME FOR COUPING ACTIVATED      ',/,&
      '                                                       ',/,&
      '            NB OF MAX INNER ITERATIONS : ',I10          ,/,&
      '            CONVERGENCE THRESHOLD      : ',E12.4        ,/)

 3000 format( &
'-------------------------------------------------------------',/)
 4000 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE INTERNAL STRUCTURES SPECIFICATION ',/,&
'@                                                            ',/,&
'@    The number of defined structures is greater than the    ',/,&
'@      allowed maximum NSTRMX:                               ',/,&
'@      Number of defined structures: ',I10                    ,/,&
'@      Number of allowed structures: ',I10                    ,/,&
'@                                                            ',/,&
'@    The calculation will not be run.                        ',/,&
'@                                                            ',/,&
'@    Decrease the number of structures                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4001 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE INTERNAL MOBILE STRUCTURES        ',/,&
'@             SPECIFICATION                                  ',/,&
'@                                                            ',/,&
'@    The number of defined structures is different from the  ',/,&
'@      previous calculation:                                 ',/,&
'@      Number of structures previous calculation: ',I10       ,/,&
'@      Number of structures current  calculation: ',I10       ,/,&
'@                                                            ',/,&
'@    The calculation will not be run.                        ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restart file or the structures     ',/,&
'@      specifications in usstru.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4002 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE EXTERNAL MOBILE STRUCTURES        ',/,&
'@             SPECIFICATION (CODE_ASTER COUPLING)            ',/,&
'@                                                            ',/,&
'@    The number of defined structures is greater than the    ',/,&
'@      allowed maximum NASTMX:                               ',/,&
'@      Number of defined structures: ',I10                    ,/,&
'@      Number of allowed structures: ',I10                    ,/,&
'@                                                            ',/,&
'@    The calculation will not be run.                        ',/,&
'@                                                            ',/,&
'@    Decrease the number of structures                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4003 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE EXTERNAL MOBILE STRUCTURES        ',/,&
'@             SPECIFICATION (CODE_ASTER COUPLING)            ',/,&
'@                                                            ',/,&
'@    The number of defined structures is different from the  ',/,&
'@      previous calculation:                                 ',/,&
'@      Number of structures previous calculation: ',I10       ,/,&
'@      Number of structures current  calculation: ',I10       ,/,&
'@                                                            ',/,&
'@    The calculation will not be run.                        ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restart file or the structures     ',/,&
'@      specifications in usstru.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! End
!----

end subroutine
