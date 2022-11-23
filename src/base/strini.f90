!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
integer          ifac  , istr, icompt
integer          mbstru, mbaste

integer          indast, verbosity, visualization

integer, allocatable, dimension(:) :: face_ids

!===============================================================================

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_ast_coupling_initialize(verbosity, visualization, &
                                        nalimx, epalim)           &
    bind(C, name='cs_ast_coupling_initialize')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: verbosity, visualization, nalimx
    real(kind=c_double), value :: epalim
  end subroutine cs_ast_coupling_initialize

  subroutine cs_ast_coupling_geometry(n_faces, face_ids, almax) &
    bind(C, name='cs_ast_coupling_geometry')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: n_faces
    integer(c_int), dimension(*), intent(in) :: face_ids
    real(kind=c_double), value :: almax
  end subroutine cs_ast_coupling_geometry

end interface

!===============================================================================
! 1. INITIALISATION
!===============================================================================

do istr = 1, nstrmx
  dtstr(istr) = dt(1)
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

! External structures: code_saturne / code_aster coupling

call uiaste(idfstr)
call usaste(idfstr)

! TODO set verbosity and visualization levels from GUI and user-defined
! functions (or build base structure earlier and allow settings)

verbosity = 1
visualization = 1

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

  nbfast = 0

  ! Number of faces coupled with code_aster
  do ifac = 1, nfabor
    istr = idfstr(ifac)
    if (istr.lt.0) then
      nbfast = nbfast + 1
    endif
  enddo

  ! Allocate temporary arrays
  allocate(face_ids(nbfast))

  indast = 0
  do ifac = 1, nfabor
    istr = idfstr(ifac)
    if (istr.lt.0) then
      indast = indast + 1
      face_ids(indast) = ifac - 1   ! O-based numbering for C
    endif
  enddo
  nbfast = indast

  ! Exchange code_aster coupling parameters
  call cs_ast_coupling_initialize(verbosity, visualization, nalimx, epalim)

  ! Send geometric information to code_aster
  call cs_ast_coupling_geometry(nbfast, face_ids, almax)

  ! Free memory
  deallocate(face_ids)

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
    /,'COUPLING MODE FOR STRUCTURES NOT ACTIVATED'             ,/,&
      '              NBSTRU = ',I10                            ,/)
 2001 format( &
      '            NALIMX USELESS AND SET TO 1'                ,/)
 2002 format( &
    /,'CODE_ASTER COUPLING MODE NOT ACTIVATED'     ,/,&
      '              NBASTE = ',I10                            ,/)
 2010 format( &
    /,'COUPLING MODE FOR STRUCTURES ACTIVATED'     ,/,&
      '              WITH NBSTRU = ',I10   ,' STRUCTURE(S)'    ,/,&
      ''                                                       ,/,&
      '            NEWMARK COEFFICIENTS:'                      ,/,&
      '              ALPNMK = ',E12.4                          ,/,&
      '              BETNMK = ',E12.4                          ,/,&
      '              GAMNMK = ',E12.4                          ,/,&
      ''                                                       ,/,&
      '            MONITORING FILES FOR STRUCTURES:'           ,/,&
      '                 IHISTR = ',I4,' ( 1 : activated)'      ,/)
 2012 format( &
    /,'CODE_ASTER COUPLING MODE ACTIVATED'         ,/,&
      '              WITH NBASTE = ',I10   ,' STRUCTURE(S)'    ,/)
 2020 format( &
    /,'EXPLICIT SCHEME FOR COUPLING ACTIVATED'     ,/,&
      ''                                                       ,/,&
      '            SCHEME COEFFICIENTS:'                       ,/,&
      '              AEXXST = ',E12.4                          ,/,&
      '              BEXXST = ',E12.4                          ,/,&
      '              CFOPRE = ',E12.4                          ,/)
 2030 format( &
    /,'IMPLICIT SCHEME FOR COUPING ACTIVATED'                  ,/,&
      ''                                                       ,/,&
      '            NB OF MAX INNER ITERATIONS : ',I10          ,/,&
      '            CONVERGENCE THRESHOLD      : ',E12.4        ,/)

 3000 format( &
'-------------------------------------------------------------',/)
 4000 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN THE INTERNAL STRUCTURES SPECIFICATION' ,/,&
'@'                                                            ,/,&
'@    The number of defined structures is greater than the'    ,/,&
'@      allowed maximum NSTRMX:'                               ,/,&
'@      Number of defined structures: ',I10                    ,/,&
'@      Number of allowed structures: ',I10                    ,/,&
'@'                                                            ,/,&
'@    The calculation will not be run.'                        ,/,&
'@'                                                            ,/,&
'@    Decrease the number of structures'                       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 4001 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN THE INTERNAL MOBILE STRUCTURES'        ,/,&
'@             SPECIFICATION'                                  ,/,&
'@'                                                            ,/,&
'@    The number of defined structures is different from the'  ,/,&
'@      previous calculation:'                                 ,/,&
'@      Number of structures previous calculation: ',I10       ,/,&
'@      Number of structures current  calculation: ',I10       ,/,&
'@'                                                            ,/,&
'@    The calculation will not be run.'                        ,/,&
'@'                                                            ,/,&
'@    Verify the auxiliary restart file or the structures'     ,/,&
'@      specifications in usstru.'                             ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 4002 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN THE EXTERNAL MOBILE STRUCTURES'        ,/,&
'@             SPECIFICATION (CODE_ASTER COUPLING)'            ,/,&
'@'                                                            ,/,&
'@    The number of defined structures is greater than the'    ,/,&
'@      allowed maximum NASTMX:'                               ,/,&
'@      Number of defined structures: ',I10                    ,/,&
'@      Number of allowed structures: ',I10                    ,/,&
'@'                                                            ,/,&
'@    The calculation will not be run.'                        ,/,&
'@'                                                            ,/,&
'@    Decrease the number of structures'                       ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 4003 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN THE EXTERNAL MOBILE STRUCTURES'        ,/,&
'@             SPECIFICATION (CODE_ASTER COUPLING)'            ,/,&
'@'                                                            ,/,&
'@    The number of defined structures is different from the'  ,/,&
'@      previous calculation:'                                 ,/,&
'@      Number of structures previous calculation: ',I10       ,/,&
'@      Number of structures current  calculation: ',I10       ,/,&
'@'                                                            ,/,&
'@    The calculation will not be run.'                        ,/,&
'@'                                                            ,/,&
'@    Verify the auxiliary restart file or the structures'     ,/,&
'@      specifications in usstru.'                             ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

!----
! End
!----

end subroutine
