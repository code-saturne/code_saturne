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

subroutine tstvec &
!================

 ( ncelet , ncel   , nfac   , nfabor ,                            &
   ifacel , ifabor ,                                              &
   iworkf , ismbs  , ismbv  ,                                     &
   rworkf , rsmbs  , rsmbv )

!===============================================================================
! FONCTION :
! ---------

! VERIFICATION DE LA VECTORISATION APRES RENUMEROTATION

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac  /nfabor    ! e  ! <-- ! nombre total de faces internes/de brd          !
! ifacel           ! te ! <-- ! no des elts voisins d'une face intern          !
! ifabor           ! te ! <-- ! no de l'elt voisin d'une face de bord          !
! nfabor  )        !    !     !                                                !
! iworkf(*         ! te ! --- ! tab de trav de dim max(nfac,nfabor)            !
! ismbs (ncelet    ! te ! --- ! tab de trav pour assemblage scalaire           !
! ismbv (ncelet    ! te ! --- ! tab de trav pour assemblage vectoriel          !
! ipnfaw           ! te ! --- ! tab de trav pour ipnfac                        !
!   (nfac+1)       !    !     !                                                !
! nodfaw           ! te ! --- ! tab de trav pour nodfac                        !
!   (lndfac)       !    !     !                                                !
! ipnfbw           ! te ! --- ! tab de trav pour ipnfbr                        !
!   (nfabor+1)     !    !     !                                                !
! nodfbw           ! te ! --- ! tab de trav pour nodfbr                        !
!   (lndfbr)       !    !     !                                                !
! rworkf(*         ! tr ! --- ! tab de trav de dim max(nfac,nfabor)            !
! rsmbs (ncelet    ! tr ! --- ! tab de trav pour assemblage scalaire           !
! rsmbv (ncelet    ! tr ! --- ! tab de trav pour assemblage vectoriel          !
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
use entsor
use parall

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel, nfac, nfabor
integer          ifacel(2,nfac),ifabor(nfabor)
integer          iworkf(*), ismbs(ncelet), ismbv(ncelet)
double precision rworkf(*), rsmbs(ncelet), rsmbv(ncelet)

! Local variables

integer          ii, jj, iok, ifac
integer          iel, istop


!===============================================================================

! -----> Test d'assemblage sur des entiers

istop = 0

! --- Faces internes

if(ivecti.eq.1) then

  do ifac = 1, nfac
    iworkf(ifac) = 1
  enddo
  do iel = 1, ncelet
    ismbv(iel) = 0
    ismbs(iel) = 0
  enddo

!CDIR NODEP
  do ifac = 1, nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    ismbv(ii) = ismbv(ii) + iworkf(ifac)
    ismbv(jj) = ismbv(jj) + iworkf(ifac)
  enddo

! VECTORISATION NON FORCEE
  do ifac = 1, nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    ismbs(ii) = ismbs(ii) + iworkf(ifac)
    ismbs(jj) = ismbs(jj) + iworkf(ifac)
  enddo

  iok = 0
  do iel = 1, ncel
    if(ismbs(iel).ne.ismbv(iel)) then
      iok = iok - 1
    endif
  enddo

  if(iok.ne.0) then
    write(nfecra,3101)iok
    istop = 1
  endif

endif

! --- Faces de bord

if(ivectb.eq.1) then

  do ifac = 1, nfabor
    iworkf(ifac) = 1
  enddo
  do iel = 1, ncel
    ismbv(iel) = 0
    ismbs(iel) = 0
  enddo

!CDIR NODEP
  do ifac = 1, nfabor
    ii = ifabor(ifac)
    ismbv(ii) = ismbv(ii) + iworkf(ifac)
  enddo

! VECTORISATION NON FORCEE
  do ifac = 1, nfabor
    ii = ifabor(ifac)
    ismbs(ii) = ismbs(ii) + iworkf(ifac)
  enddo

  iok = 0
  do iel = 1, ncel
    if(ismbs(iel).ne.ismbv(iel)) then
      iok = iok - 1
    endif
  enddo

  if(iok.ne.0) then
    write(nfecra,3201)iok
    istop = 1
  endif

endif


! -----> Test d'assemblage sur des reels

! --- Faces internes

if(ivecti.eq.1) then

  do ifac = 1, nfac
    rworkf(ifac) = 1.d0
  enddo
  do iel = 1, ncelet
    rsmbv(iel) = 0.d0
    rsmbs(iel) = 0.d0
  enddo

!CDIR NODEP
  do ifac = 1, nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    rsmbv(ii) = rsmbv(ii) + rworkf(ifac)
    rsmbv(jj) = rsmbv(jj) + rworkf(ifac)
  enddo

! VECTORISATION NON FORCEE
  do ifac = 1, nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    rsmbs(ii) = rsmbs(ii) + rworkf(ifac)
    rsmbs(jj) = rsmbs(jj) + rworkf(ifac)
  enddo

  iok = 0
  do iel = 1, ncel
    if(abs(rsmbs(iel)-rsmbv(iel)).gt.1.d-20) then
      iok = iok - 1
    endif
  enddo

  if(iok.ne.0) then
    write(nfecra,3102)iok
    istop  = 1
  endif

endif

! --- Faces de bord

if(ivectb.eq.1) then

  do ifac = 1, nfabor
    rworkf(ifac) = 1.d0
  enddo
  do iel = 1, ncel
    rsmbv(iel) = 0.d0
    rsmbs(iel) = 0.d0
  enddo

!CDIR NODEP
  do ifac = 1, nfabor
    ii = ifabor(ifac)
    rsmbv(ii) = rsmbv(ii) + rworkf(ifac)
  enddo

! VECTORISATION NON FORCEE
  do ifac = 1, nfabor
    ii = ifabor(ifac)
    rsmbs(ii) = rsmbs(ii) + rworkf(ifac)
  enddo

  iok = 0
  do iel = 1, ncel
    if(abs(rsmbs(iel)-rsmbv(iel)).gt.1.d-20) then
      iok = iok - 1
    endif
  enddo

  if(iok.ne.0) then
    write(nfecra,3202)iok
    istop  = 1
  endif

endif

if(istop.ne.0) then
  write(nfecra,9000)
  call csexit (1)
endif

!===============================================================================
! 6. FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 3101 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS tstvec                           ',/,&
'@    =========                                               ',/,&
'@    PROBLEME D''ASSEMBLAGE D''ENTIERS AUX FACES INTERNES    ',/,&
'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Verifier le maillage.                                     ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3201 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS tstvec                           ',/,&
'@    =========                                               ',/,&
'@    PROBLEME D''ASSEMBLAGE D''ENTIERS AUX FACES DE BORD     ',/,&
'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Verifier le maillage.                                     ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3102 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS tstvec                           ',/,&
'@    =========                                               ',/,&
'@    PROBLEME D''ASSEMBLAGE DE REELS   AUX FACES INTERNES    ',/,&
'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Verifier le maillage.                                     ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3202 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS tstvec                           ',/,&
'@    =========                                               ',/,&
'@    PROBLEME D''ASSEMBLAGE DE REELS   AUX FACES DE BORD     ',/,&
'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Verifier le maillage.                                     ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9000 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS tstvec                           ',/,&
'@    =========                                               ',/,&
'@     CONFIGURATION IMPREVUE A LA RENUMEROTATION DES FACES   ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Verifier le maillage.                                     ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 3101 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN tstvec                                ',/,&
'@    ========                                                ',/,&
'@    PROBLEM IN ASSEMBLING INTEGERS AT THE INTERIOR FACES    ',/,&
'@    ',I10   ,' OCCURRENCES OF THE PROBLEM                   ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify the mesh.                                          ',/,&
'@  Contact the support.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3201 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN tstvec                                ',/,&
'@    ========                                                ',/,&
'@    PROBLEM IN ASSEMBLING INTEGERS AT THE BOUNDARY FACES    ',/,&
'@    ',I10   ,' OCCURRENCES OF THE PROBLEM                   ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify the mesh.                                          ',/,&
'@  Contact the support.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3102 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN tstvec                                ',/,&
'@    ========                                                ',/,&
'@    PROBLEM IN ASSEMBLING REALS AT THE INTERIOR FACES       ',/,&
'@    ',I10   ,' OCCURRENCES OF THE PROBLEM                   ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify the mesh.                                          ',/,&
'@  Contact the support.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3202 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN tstvec                                ',/,&
'@    ========                                                ',/,&
'@    PROBLEM IN ASSEMBLING REALS AT THE BOUNDARY FACES       ',/,&
'@    ',I10   ,' OCCURRENCES OF THE PROBLEM                   ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify the mesh.                                          ',/,&
'@  Contact the support.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9000 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN tstvec                                ',/,&
'@    ========                                                ',/,&
'@     UNEXPECTED CONFIGURATION DURING THE FACES RENUMBERING  ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify the mesh.                                          ',/,&
'@  Contact the support.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

return
end subroutine
