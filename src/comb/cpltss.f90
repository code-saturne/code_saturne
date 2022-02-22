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

subroutine cpltss &
!================

 ( iscal  ,                                                       &
   itypfb ,                                                       &
   smbrs  , rovsdt , tslagr )

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
! iscal            ! i  ! <-- ! scalar number                                  !
! itypfb(nfabor    ! te ! --> ! type des faces de bord                         !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
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
use lagran
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal

integer          itypfb(nfabor)

double precision smbrs(ncelet), rovsdt(ncelet)
double precision tslagr(ncelet,*)

! Local variables

character(len=80) :: chaine
integer          ivar   , iel
integer          iscala , icha

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

! --- Nom de la variable associee au scalaire a traiter ISCAL
call field_get_label(ivarfl(ivar), chaine)

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES
!===============================================================================

! --> Terme source pour les matieres volatiles legeres

if ( ivar.ge.isca(if1m(1)) .and. ivar.le.isca(if1m(ncharb)) ) then

  if (vcopt%iwarni.ge.1) then
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

  if (vcopt%iwarni.ge.1) then
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

  if (vcopt%iwarni.ge.1) then
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

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calcul des termes sources explicite et implicite
!      relatif aux echanges interfaciaux entre phases

! ---- Calcul des termes sources explicite et implicite
!      relatif aux termes de production et de dissipation

!      Pointeur relatif au scalaire associe
!      (0 si pas de scalaire associe)
  iscala = 0

  call cpltsv                                                     &
  !==========
 ( iscal  , iscala ,                                              &
   itypfb ,                                                       &
   smbrs  , rovsdt )

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
