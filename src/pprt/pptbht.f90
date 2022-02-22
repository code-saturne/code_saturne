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

subroutine pptbht &
!================

 ( ncoel  ,                                                       &
   nomcel , ehcoel , cpcoel , wmolce )

!===============================================================================
! FONCTION :
! --------

!         PHYSIQUES PARTICULIERES

!           CALCUL DE L'ENTHALPIE ET DU CP
!                    A PARTIR DE LA BANDE DE JANAF


! Arguments
!_______________.____._____.________________________________________________.
!    nom        !type!mode !                   role                         !
!_______________!____!_____!________________________________________________!
! ncoel         ! e  ! <-- ! nombre de const. elem.                         !
! nomcel(ngazem)! a  ! <-- ! nom des constituants elementaires              !
! ehcoel        ! tr !  <- ! enthalpie pour chaque constituant              !
! (ngazem,npot) !    !     !                elementaire                     !
! cpcoel        ! tr !  <- ! cp pour chaque constituant                     !
! (ngazem,npot) !    !     !                elementaire                     !
! wmolce(ngazem)! tr !  <- ! masse molaire de chaque constituant            !
!_______________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstphy
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer          ncoel

character(len=12) :: nomcel(ngazem)

double precision ehcoel(ngazem,npot) , cpcoel(ngazem,npot)
double precision wmolce (ngazem)

! Local variables

character(len=256) :: pathdatadir
character(len=40) :: dummy
character(len=12) :: nomesp

integer          ind , iches , indtp , inicff , injcff, impjnf
integer          ne   , nt  , nc , iok
integer          icoeff(ngazem)

double precision cth  , ctc

double precision tlim(3) , wcoeff(2,7) , coeff(ngazem,2,7)

!===============================================================================


!===============================================================================
! 2. LECTURE DU FICHIER DE DONNEES THERMODYNAMIQUES TABLE DE JANAF
!===============================================================================

! Initialisation


do iches= 1, 12
  nomesp(iches:iches)=' '
enddo

do ne = 1 , ngazem
  icoeff(ne) = 0
  do inicff = 1, 2
    do injcff = 1, 7
      coeff(ne,inicff,injcff) = 0.d0
    enddo
  enddo
enddo

do ne = 1 , ncoel
  do nt = 1, npo
    cpcoel(ne,nt)= 0.d0
    ehcoel(ne,nt)= 0.d0
  enddo
enddo

impjnf = impfpp

call csdatadir(len(pathdatadir), pathdatadir)

open(unit=impjnf, file=trim(pathdatadir)// '/data/thch/JANAF',  &
     status='old' , form='formatted')

read(impjnf,'(a)') dummy

! Lecture des domaines de temperature

read (impjnf,*) (tlim(indtp) , indtp=1,3)

! Boucle de lecture des especes chimiques avec stockage partiel

 5    continue

read (impjnf,'(a12,6x,a6)') nomesp, dummy

if (nomesp(1:3).EQ.'END') GOTO 100

read (impjnf,*) (wcoeff(1,injcff), injcff=1,5)
read (impjnf,*) (wcoeff(1,injcff), injcff=6,7),                   &
                (wcoeff(2,injcff), injcff=1,3)
read (impjnf,*) (wcoeff(2,injcff), injcff=4,7)

! On ne stocke les coefficients que si
!  l'espece consideree fait partie de l'exemple

do ne = 1, ncoel
  if ( nomcel(ne).eq.nomesp ) then
    icoeff(ne) = 1
     do inicff = 1, 2
       do injcff = 1, 7
         coeff(ne,inicff,injcff) = wcoeff(inicff,injcff)
       enddo
     enddo
  endif
enddo

goto 5

 100  continue

! Arret de la lecture si tous les renseignements necessaires ont ete
! enregistres

close(impjnf)

! Test et stop eventuel

iok = 0
do ne = 1, ncoel
  if(icoeff(ne).eq.0) then
    iok = iok + 1
    write(nfecra,1000) nomcel(ne)
  endif
enddo
if(iok.ne.0) then
  write(nfecra,1100) iok
  call csexit (1)
  !==========
endif

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PHYSIQUE PARTICULIERE                       ',/,&
'@    =========                                               ',/,&
'@    L''ESPECE ',A12     ,' EST INCONNUE DANS JANAF          ',/,&
'@                                                            ',/,&
'@                                                            '  )
 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE                                   ',/,&
'@                                                            ',/,&
'@    LE FICHIER PARAMETRIQUE FAIT REFERENCE A ',I10           ,/,&
'@       ESPECE(S) INCONNUE(S) DANS JANAF (data/thch).        ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!======================================================================================
! 3. CALCUL DES ENTHALPIES ET DES CP
!======================================================================================

! Constante des gaz parfaits en J/mol/K

do nt = 1,npo

  ! Determination du jeu de coefficients utilises
  if (th(nt) .gt. tlim(2)) then
    ind = 1
  else
    ind = 2
  endif

  do ne = 1, ncoel
    ehcoel(ne,nt)  = coeff(ne,ind,6) + coeff(ne,ind,1) * th(nt)
    cpcoel(ne,nt)  =                   coeff(ne,ind,1)
    cth = th(nt)
    ctc = 1.d0

    ! Dans la table de Janaf, les COEFF sont adimensionnels (CP/R,H/R)
    do nc = 2, 5
      cth = cth * th(nt)
      ctc = ctc * th(nt)
      ehcoel(ne,nt) = ehcoel(ne,nt)                               &
                    + coeff(ne,ind,nc) * cth / dble(nc)
      cpcoel(ne,nt) = cpcoel(ne,nt) + coeff(ne,ind,nc) * ctc
    enddo

    ! Calcul du CP et du H pour chaque espece
    ehcoel(ne,nt) = ehcoel(ne,nt) * cs_physical_constants_r / wmolce(ne)
    cpcoel(ne,nt) = cpcoel(ne,nt) * cs_physical_constants_r / wmolce(ne)
  enddo

enddo

return

end subroutine
