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

subroutine vert1d &
!================

 ( nfabor , nfpt1d , iappel ,                                     &
   ifpt1d , nppt1d , iclt1d ,                                     &
   rgpt1d , eppt1d ,                                              &
   xlmbt1 , rcpt1d , dtpt1d )

!===============================================================================
! FONCTION :
! ----------

!  VERIFICATION DES DONNEES UTILISATEUR DU MODULE THERMIQUE 1D EN PAROI

! IAPPEL = 1 (un seul appel a l'initialisation) :
!             VERIFICATION DU NOMBRE DE CELLULES OU L'ON IMPOSE UNE PAROI
!             VERIFICATION DE ISUIT1

! IAPPEL = 2 (un seul appel a l'initialisation) :
!             VERIFICATION DU REPERAGE DES CELLULES OU L'ON IMPOSE
!                                                               UNE PAROI
!             VERIFICATION DES DONNEES RELATIVE AU MAILLAGE

! IAPPEL = 3 (appel a chaque pas de temps) :
!             VERIFICATION DES TYPES DE C.L. EN PAROI EXTERIEURE
!             VERIFICATION DES VALEURS DES COEFFICIENTS PHYSIQUE DU CALCUL
!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nfpt1d           ! e  ! <-- ! nombre de faces avec module therm 1d           !
! ifpt1d           ! e  ! <-- ! numero de la face en traitement                !
!                  !    !     ! thermique en paroi                             !
! nppt1d           ! e  ! <-- ! nombre de points de discretisation             !
!                  !    !     ! dans la paroi                                  !
! eppt1d           ! r  ! <-- ! epaisseur de la paroi                          !
! rgpt1d           ! r  ! <-- ! raison du maillage                             !
! iclt1d           ! e  ! <-- ! type de condition limite                       !
! xlmbt1           ! r  ! <-- ! diffusivite thermique                          !
! rcpt1d           ! r  ! <-- ! rocp                                           !
! dtpt1d           ! tr ! <-- ! pas de temps                                   !
! iappel           ! e  ! <-- ! indique les donnes a renvoyer                  !
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
use cstnum
use entsor
use optcal

!===============================================================================

implicit none

! Arguments

integer          nfabor , nfpt1d
integer          iappel

integer          ifpt1d(nfpt1d) , nppt1d(nfpt1d) , iclt1d(nfpt1d)

double precision eppt1d(nfpt1d) , rgpt1d(nfpt1d)
double precision xlmbt1(nfpt1d) , rcpt1d(nfpt1d) , dtpt1d(nfpt1d)

! Variables locales
integer          ii, ifac

!===============================================================================




if(iappel.eq.1) then

  if ((nfpt1d.lt.0).or.(nfpt1d.gt.nfabor)) then
    write(nfecra,1000)nfabor,nfpt1d
    call csexit(1)
  endif
  if (isuit1.ne.0 .and. isuit1.ne.1) then
    write(nfecra,1010)isuit1
    call csexit(1)
  endif

else if (iappel.eq.2) then

  do ii = 1, nfpt1d
    if (ifpt1d(ii).lt.0 .or. ifpt1d(ii).gt.nfabor) then
      write(nfecra,2000)nfabor,ii,ifpt1d(ii)
      call csexit(1)
    endif
  enddo

  do ii = 1, nfpt1d
    ifac = ifpt1d(ii)
    if (nppt1d(ii).le.0) then
      write(nfecra,2010)ii,ifac,nppt1d(ii)
      call csexit(1)
    endif
    if (eppt1d(ii).le.0.d0) then
      WRITE(NFECRA,2020)'EPPT1D','EPPT1D',II,EPPT1D(II),IFAC
      call csexit(1)
    endif
    if (rgpt1d(ii).le.0.d0) then
      WRITE(NFECRA,2020)'RGPT1D','RGPT1D',II,RGPT1D(II),IFAC
      call csexit(1)
    endif
  enddo

else if (iappel.eq.3) then

  do ii = 1, nfpt1d
    ifac = ifpt1d(ii)
    if (iclt1d(ii).ne.1 .and. iclt1d(ii).ne.3) then
      WRITE(NFECRA,3000)'ICLT1D','ICLT1D',II,ICLT1D(II),IFAC
      call csexit(1)
    endif
    if (xlmbt1(ii).le.0.d0) then
      WRITE(NFECRA,2020)'XLMBT1','XLMBT1',II,XLMBT1(II),IFAC
      call csexit(1)
    endif
    if (rcpt1d(ii).le.0.d0) then
      WRITE(NFECRA,2020)'RCPT1D','RCPT1D',II,RCPT1D(II),IFAC
      call csexit(1)
    endif
    if (dtpt1d(ii).le.0.d0) then
      WRITE(NFECRA,2020)'DTPT1D','DTPT1D',II,DTPT1D(II),IFAC
      call csexit(1)
    endif
  enddo

endif

!----
! FORMATS
!----

#if defined(_CS_LANG_FR)

 1000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    MODULE THERMIQUE 1D EN PAROI                            ',/,&
'@                                                            ',/,&
'@    NFPT1D DOIT ETRE POSITIF ET INFERIEUR A NFABOR          ',/,&
'@    ON A ICI                                                ',/,&
'@       NFABOR=',I10                                          ,/,&
'@       NFPT1D=',I10                                          ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uspt1d.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    MODULE THERMIQUE 1D EN PAROI                            ',/,&
'@                                                            ',/,&
'@    ISUIT1 DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uspt1d.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    MODULE THERMIQUE 1D EN PAROI                            ',/,&
'@                                                            ',/,&
'@    LE TABLEAU IFPT1D(II) DOIT RENVOYER A UN NUMERO DE FACE ',/,&
'@                                                     DE BORD',/,&
'@    ON A ICI                                                ',/,&
'@       NFABOR            =',I10                              ,/,&
'@       IFPT1D(',I10,   ')=',I10                              ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uspt1d.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2010   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    MODULE THERMIQUE 1D EN PAROI                            ',/,&
'@                                                            ',/,&
'@    LE TABLEAU NPPT1D(II) DOIT RENVOYER A UN ENTIER POSITIF ',/,&
'@    ON A ICI                                                ',/,&
'@       NPPT1D(',I10,   ')=',I10                              ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uspt1d.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2020   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    MODULE THERMIQUE 1D EN PAROI                            ',/,&
'@                                                            ',/,&
'@    LE TABLEAU ',A6,' DOIT RENVOYER A UN REEL > 0           ',/,&
'@    ON A ICI                                                ',/,&
'@       ',A6,'(',I10,   ')=',E14.5                            ,/,&
'@       (FACE DE BORD NUMERO ',I10   ,')                     ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uspt1d.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    MODULE THERMIQUE 1D EN PAROI                            ',/,&
'@                                                            ',/,&
'@    LE TABLEAU ',A6,' NE PEUT CONTENIR QUE LES VALEURS      ',/,&
'@                                                     1 OU 3',/, &
'@    ON A ICI                                                ',/,&
'@       ',A6,'(',I10,   ')=',I10                              ,/,&
'@       (FACE DE BORD NUMERO ',I10   ,')                     ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uspt1d.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION            ',/,&
'@    ========                                                ',/,&
'@    1D-WALL THERMAL MODULE                                  ',/,&
'@                                                            ',/,&
'@    NFPT1D MUST BE POSITIVE AND LOWER THAN NFABOR           ',/,&
'@    ONE HAS HERE                                            ',/,&
'@       NFABOR=',I10                                          ,/,&
'@       NFPT1D=',I10                                          ,/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify uspt1d.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION            ',/,&
'@    ========                                                ',/,&
'@    1D-WALL THERMAL MODULE                                  ',/,&
'@                                                            ',/,&
'@    ISUIT1 MUST BE AN INTEGER EQUAL TO 0 OR 1               ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify uspt1d.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION            ',/,&
'@    ========                                                ',/,&
'@    1D-WALL THERMAL MODULE                                  ',/,&
'@                                                            ',/,&
'@    THE ARRAY IFPT1D(II) MUST GIVE A BOUNDARY FACE NUMBER   ',/,&
'@    ONE HAS HERE                                            ',/,&
'@       NFABOR            =',I10                              ,/,&
'@       IFPT1D(',I10,   ')=',I10                              ,/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify uspt1d.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2010   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION            ',/,&
'@    ========                                                ',/,&
'@    1D-WALL THERMAL MODULE                                  ',/,&
'@                                                            ',/,&
'@    THE ARRAY NPPT1D(II) MUST GIVE A POSITIVE INTEGER       ',/,&
'@    ONE HAS HERE                                            ',/,&
'@       NPPT1D(',I10,   ')=',I10                              ,/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify uspt1d.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2020   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION            ',/,&
'@    ========                                                ',/,&
'@    1D-WALL THERMAL MODULE                                  ',/,&
'@                                                            ',/,&
'@    THE ARRAY ',A6,' MUST GIVE A POSITIVE REAL              ',/,&
'@    ONE HAS HERE                                            ',/,&
'@       ',A6,'(',I10,   ')=',E14.5                            ,/,&
'@       (BOUNDARY FACE NUMBER ',I10   ,')                    ',/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify uspt1d.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION            ',/,&
'@    ========                                                ',/,&
'@    1D-WALL THERMAL MODULE                                  ',/,&
'@                                                            ',/,&
'@    THE ARRAY ',A6,' CAN ONLY TAKE THE VALUES 1 OR 3        ',/,&
'@    ONE HAS HERE                                            ',/,&
'@       ',A6,'(',I10,   ')=',I10                              ,/,&
'@       (BOUNDARY FACE NUMBER ',I10   ,')                    ',/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify uspt1d.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

return

end subroutine

