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

subroutine cpltss &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp , ia     ,                            &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt , tslagr ,                                     &
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , w10    , w11    ,          &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!   ON PRECISE LES TERMES SOURCES POUR UN SCALAIRE PP
!   SUR UN PAS DE TEMPS

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
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
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
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
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! viscf(nfac)      ! tr ! --- ! tableau de travail    faces internes           !
! viscb(nfabor     ! tr ! --- ! tableau de travail    faces de bord            !
! xam(nfac,2)      ! tr ! --- ! tableau de travail    faces de bord            !
! w1..11(ncelet    ! tr ! --- ! tableau de travail    cellules                 !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use coincl
use cpincl
use ppincl
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          itypfb(nfabor)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          izfppp(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smbrs(ncelet), rovsdt(ncelet)
double precision tslagr(ncelet,*)
double precision viscf(nfac), viscb(nfabor)
double precision xam(nfac,2)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet), w11(ncelet)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia , idebra
integer          ivar   , iel
integer          iscala , icha

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
chaine = nomvar(ipprtp(ivar))



!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES
!===============================================================================

! --> Terme source pour les matieres volatiles legeres

if ( ivar.ge.isca(if1m(1)) .and. ivar.le.isca(if1m(ncharb)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Contribution du TS interfacial aux bilans explicite et implicite

  icha = ivar-isca(if1m(1))+1
  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + tslagr(iel,itsmv1(icha))
!          ROVSDT(IEL) = ROVSDT(IEL) + ZERO
  enddo

endif

! --> Terme source pour les matieres volatiles lourdes

if ( ivar.ge.isca(if2m(1)) .and. ivar.le.isca(if2m(ncharb)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Contribution du TS interfacial pour le bilan explicite

  icha = ivar-isca(if2m(1))+1
  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  +  tslagr(iel,itsmv2(icha))
!          ROVSDT(IEL) = ROVSDT(IEL) + ZERO
  enddo

endif

! --> Terme source pour le traceur 3 (C de la comb. het.)

if ( ivar.eq.isca(if3m) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Contribution du TS interfacial aux bilans explicite et implicite

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + tslagr(iel,itsco)
!          ROVSDT(IEL) = ROVSDT(IEL) + ZERO
  enddo

endif

! --> Terme source pour la variance du traceur 4 (Air)

if ( ivar.eq.isca(if4p2m) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calcul des termes sources explicite et implicite
!      relatif aux echanges interfaciaux entre phases

! -> appel commente => SMBRS et ROVSDT non modifies
!       NUMTRA = 4
!       CALL CPTSVI
!!==========
!     & ( NCELET , NCEL   , NUMTRA ,
!     &   RTP    , PROPCE , VOLUME ,
!     &   SMBRS  , ROVSDT ,
!     &   W1     , W2     ,
!     &   W3 )


! ---- Calcul des termes sources explicite et implicite
!      relatif aux termes de production et de dissipation

!      Pointeur relatif au scalaire associe
!      (0 si pas de scalaire associe)
  iscala = 0

  call cpltsv                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  , iscala ,                                              &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   smbrs  , rovsdt ,                                              &
   viscb  ,                                                       &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     ,                                     &
   ra     )

endif

!--------
! FORMATS
!--------

 1000 format(' TERMES SOURCES PHYSIQUE PARTICULIERE POUR LA VARIABLE '  &
       ,a8,/)

!----
! FIN
!----

return

end subroutine
