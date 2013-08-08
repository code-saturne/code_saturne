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

subroutine eltssc &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE POUR LE MODULE ELECTRIQUE
! CALCUL DES TERMES SOURCE ELECTRIQUES POUR L'ENTHALPIE
! CALCUL DES TERMES SOURCES POUR LE POTENTIEL VECTEUR

! ATTENTION : LE TRAITEMENT DES TERMES SOURCES EST DIFFERENT
! ---------   DE CELUI DE USTSSC.F

! ON RESOUT ROVSDT*D(VAR) = SMBRS

! ROVSDT ET SMBRS CONTIENNENT DEJA D'EVENTUELS TERMES SOURCES
!  UTILISATEUR. IL FAUT DONC LES INCREMENTER ET PAS LES
!  ECRASER

! POUR DES QUESTIONS DE STABILITE, ON NE RAJOUTE DANS ROVSDT
!  QUE DES TERMES POSITIFS. IL N'Y A PAS DE CONTRAINTE POUR
!  SMBRS

! DANS LE CAS D'UN TERME SOURCE EN CEXP + CIMP*VAR ON DOIT
! ECRIRE :
!          SMBRS  = SMBRS  + CEXP + CIMP*VAR
!          ROVSDT = ROVSDT + MAX(-CIMP,ZERO)

! ON FOURNIT ICI ROVSDT ET SMBRS (ILS CONTIENNENT RHO*VOLUME)
!    SMBRS en kg variable/s :
!     ex : pour la vitesse            kg m/s2
!          pour les temperatures      kg degres/s
!          pour les enthalpies        Joules/s
!    ROVSDT en kg /s

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! itypfb(nfabor    ! te ! --> ! type des faces de bord                         !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
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
use period
use ppppar
use ppthch
use ppincl
use elincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          itypfb(nfabor)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

character*80     chaine
integer          ivar , iel
integer          ipcdc1, ipcdc2, ipcdc3
integer          ipcefj

double precision valmin,valmax

double precision, allocatable, dimension(:) :: w1

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate a temporary array
allocate(w1(ncelet))


! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
chaine = nomvar(ipprtp(ivar))


!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES ET VARIABLES STANDARD ET
!   COMMUNES A TOUTES LES VERSIONS ELECTRIQUES : EFFET JOULE, E, J
!         CAS IELJOU >= 1 ou IELARC >= 1 ou IELION >= 1
!===============================================================================

!   2.1 Terme source pour l'enthalpie :
!  ----------------------------------

if ( ivar.eq.isca(ihm) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  ipcefj = ipproc(iefjou)

!     Pour eviter de prendre en compte des TS trop violents
!     sur les premiers pas de temps (gradient sur le potentiel)
!     on ne prend en compte les termes sources qu'a partir du
!     troisieme pas de temps. Sinon on laisse ROVSDT et SMBRS
!     inchanges (pas de TS si ce n'est ceux de ustssc)

  if ( ntcabs .ge. 3 ) then
    do iel = 1, ncel
      w1(iel) =  propce(iel,ipcefj)*volume(iel)
    enddo

!        Arc electrique : rajout du TS radiatif si l'utilisateur
!        l'a donne dans le fichier dp_elec (purement explicite)

    if( ippmod(ielarc).ge. 1 ) then
      if ( ixkabe .eq. 2 ) then
        do iel = 1, ncel
          w1(iel) = w1(iel)                                       &
                       -propce(iel,ipproc(idrad))*volume(iel)
        enddo
      endif
    endif

    do iel = 1, ncel
       smbrs(iel) = smbrs(iel) + w1(iel)
!            ROVSDT(IEL) = ROVSDT(IEL) + ZERO
    enddo

    if (iwarni(ivar).ge.2) then
      valmin = w1(1)
      valmax = w1(1)
      do iel = 1, ncel
        valmin = min(valmin,w1(iel))
        valmax = max(valmax,w1(iel))
      enddo
      if (irangp.ge.0) then
        call parmax(valmax)
        call parmin(valmin)
      endif
      write(nfecra,2000) valmin,valmax
    endif
  endif
endif

 2000  format(                                                          &
 ' Termes Sources sur H  min= ',E14.5,', max= ',E14.5)
 2010  format(/,                                                  &
 ' Termes Sources sur H  : pour eviter les debuts de calcul ',/,  &
 '   trop violents, on ne prend en compte les termes sources',/,  &
 '   d''effet Joule qu''a partir du troisieme pas de temps. ',/)

!   2.1  Terme source pour les composantes potentiel vecteur : ARC ELECTRIQUE
!  ---------------------------------------------------------

if( ippmod(ielarc).ge. 2 ) then


! -->   Terme source pour les composantes potentiel vecteur
!============================================================

  if ( ivar .eq. isca(ipotva(1)) ) then

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    ipcdc1 = ipproc(idjr(1))
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) +                                   &
            permvi*propce(iel,ipcdc1)*volume(iel)
!            ROVSDT(IEL) = ROVSDT(IEL) + ZERO
    enddo

  else if ( ivar .eq. isca(ipotva(2)) ) then

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    ipcdc2 = ipproc(idjr(2))
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) +                                   &
            permvi*propce(iel,ipcdc2)*volume(iel)
!            ROVSDT(IEL) = ROVSDT(IEL) + ZERO
    enddo

  else if ( ivar .eq. isca(ipotva(3)) ) then

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    ipcdc3 = ipproc(idjr(3))
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) +                                   &
            permvi*propce(iel,ipcdc3)*volume(iel)
!            ROVSDT(IEL) = ROVSDT(IEL) + ZERO
    enddo
  endif

endif

! Free memory
deallocate(w1)

!--------
! FORMATS
!--------

 1000 format(                                                           &
'  Calcul des termes sources pour la variable : ',A8             )



!----
! FIN
!----

return
end subroutine
