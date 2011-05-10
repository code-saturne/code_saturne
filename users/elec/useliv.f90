!-------------------------------------------------------------------------------

!VERS


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

subroutine useliv &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   maxelt , lstelt ,                                              &
   ia     ,                                                       &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  , &
   ra     )

!===============================================================================
! Purpose :
! --------

! variables initialization for specific electric module :
!    the same as USINIV.F for main variables

! Cette routine est appelee en debut de calcul (suite ou non)
!     avant le debut de la boucle en temps

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul (inconnues transportees),
!     la valeur du pas de temps,

! On a repris a titre d'exemple dans cette routine utilisateur
!     l'initialisation choisie par defaut.
!

! Thermodynamic properties and transport coefficients are in
!     PROPCE (at the center of the cells), PROPFA (for internal faces),
!     PROPFB (prop aux faces de bord)
!     Ainsi,
!      PROPCE(IEL,IPPROC(IROM  )) designe ROM   (IEL)
!      PROPCE(IEL,IPPROC(IVISCL)) designe VISCL (IEL)
!      PROPCE(IEL,IPPROC(ICP   )) designe CP    (IEL)
!      PROPCE(IEL,IPPROC(IVISLS(ISCAL))) designe VISLS (IEL ,ISCAL)

!      PROPFA(IFAC,IPPROF(IFLUMA(IVAR ))) designe FLUMAS(IFAC,IVAR)

!      PROPFB(IFAC,IPPROB(IROM  )) designe ROMB  (IFAC)
!      PROPFB(IFAC,IPPROB(IFLUMA(IVAR ))) designe FLUMAB(IFAC,IVAR)

! Thermodynamic properties and transport coefficients modification
!  (ROM, VISCL, VISCLS, CP) will be done by default in  PPPHYV
!     and not here


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! maxelt           !  e ! <-- ! nb max d'elements (cell,fac,fbr)               !
! lstelt(maxelt) te ! --- ! tableau de travail                             !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa coefb      ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
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
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use elincl
use mesh

!===============================================================================

implicit none

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          maxelt, lstelt(maxelt)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ra(*)


! Local variables

integer          idebia, idebra
integer          iel, mode
integer          iesp , idimve

double precision tinit, hinit, coefe(ngazem)

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================
!
! For Joule heating by direct conduction, it stoppes
!   you have tot be sure that the enthalpy function H(T) is the right one
!
if ( ippmod(ieljou).ge.1 ) then

  write(nfecra,9010)
  call csexit (1)

! For electric arc, we continue because the value are given by defauft
!       (H(T) is given from the data file dp_ELE)
elseif(ippmod(ielarc).ge.1) then

  if(ntcabs.eq.1) then
    write(nfecra,9011)
  endif

  return

endif

 9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION : Stop in the definition of Thermal properties  ',/,&
'@    =========                                               ',/,&
'@                      for Electric module                   ',/,&
'@                                                            ',/,&
'@     The user routine uselph has to be completed            ',/,&
'@                                                            ',/,&
'@     This user routine is used to define thermal properties ',/,&
'@     It is unavoidable.                                     ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9011 format(/,                                                   &
' ELECTRIC ARC MODULE : THERMAL PROPERTIES ARE READ IN A FILE',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0. Control output
!===============================================================================

write(nfecra,9001)

!===============================================================================
! 1.  Local variables initialization
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2. Initialization : examples
!      only at the beginning of the calculation
!===============================================================================

if ( isuite.eq.0 ) then

! --> Enthalpy = H(T0) ou 0

!     For electric arc,
!     for the whole compution domain enthalpy is set to H(T0) of the 1st constituant of the gas
!
!     For Joule jeating by direct conduction,
!     enthalpy is set to zero, and the user will enter his H(T) function tabulation.
!
!   --  HINIT calculations

  if ( ippmod(ielarc).ge.1 ) then
    mode = -1
    tinit = t0
    coefe(1) = 1.d0
    if ( ngazg .gt. 1 ) then
      do iesp = 2, ngazg
        coefe(iesp) = 0.d0
      enddo
    endif
    call elthht(mode,ngazg,coefe,hinit,tinit)
  else
    mode = -1
    tinit = t0
    call usthht(mode,hinit,tinit)
  endif

!    -- Entahlpy value

  do iel = 1, ncel
    rtp(iel,isca(ihm)) = hinit
  enddo


! --> Mass fraction  = 1 ou 0

  if ( ngazg .gt. 1 ) then
    do iel = 1, ncel
      rtp(iel,isca(iycoel(1))) = 1.d0
    enddo
    do iesp = 2, ngazg-1
      do iel = 1, ncel
        rtp(iel,isca(iycoel(iesp))) = 0.d0
      enddo
    enddo
  endif


! --> Electric potentials = 0

!     -- Real Component
  do iel = 1, ncel
    rtp(iel,isca(ipotr)) = 0.d0
  enddo

!     -- Imaginary (for Joule heating by direct conduction)
  if ( ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4 ) then
    do iel = 1, ncel
      rtp(iel,isca(ipoti)) = 0.d0
    enddo
  endif

!     -- Vector potential (3D electric arc 3D)
  if ( ippmod(ielarc).ge.2 ) then
    do idimve = 1, ndimve
      do iel = 1, ncel
        rtp(iel,isca(ipotva(idimve))) = 0.d0
      enddo
    enddo
  endif

endif

!----
! FORMATS
!----

 9001 format(/,                                                   &
'                       ELECTRIC MODULE                       ',/,&
'  useliv : variables initialization by user                   ',/,&
'                                                             '  )

!----
! FIN
!----
return
end subroutine
