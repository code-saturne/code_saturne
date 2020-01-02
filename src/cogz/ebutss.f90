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

subroutine ebutss &
!================

 ( iscal  ,                                                       &
   smbrs  , rovsdt )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME PREMELANGE MODELE EBU
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
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

character(len=80) :: chaine
integer          ivar, iel

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_omg
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: cvara_scal

type(var_cal_opt) :: vcopt

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))


! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
call field_get_label(ivarfl(ivar), chaine)

! --- Numero des grandeurs physiques (voir cs_user_boundary_conditions)
call field_get_val_s(icrom, crom)

if ( ivar.eq.isca(iygfm) ) then
  call field_get_val_prev_s(ivarfl(isca(iscal)), cvara_scal)
endif

if (itytur.eq.2.or.iturb.eq.50) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
elseif (itytur.eq.3) then
  call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
  call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
  call field_get_val_prev_s(ivarfl(ir33), cvara_r33)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
elseif (iturb.eq.60) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iomg), cvara_omg)
endif

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES
!===============================================================================

if ( ivar.eq.isca(iygfm) ) then

! ---> Terme source pour la fraction massique moyenne de gaz frais

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---> Calcul de K et Epsilon en fonction du modele de turbulence

  if (itytur.eq.2) then

    do iel = 1, ncel
      w1(iel) = cvara_k(iel)
      w2(iel) = cvara_ep(iel)
    enddo

  elseif (itytur.eq.3) then

    do iel = 1, ncel
      w1(iel) = 0.5d0 *( cvara_r11(iel)                    &
                        +cvara_r22(iel)                    &
                        +cvara_r33(iel) )
      w2(iel) = cvara_ep(iel)
    enddo

  elseif (iturb.eq.50) then

    do iel = 1, ncel
      w1(iel) = cvara_k(iel)
      w2(iel) = cvara_ep(iel)
    enddo

  elseif (iturb.eq.60) then

    do iel = 1, ncel
      w1(iel) = cvara_k(iel)
      w2(iel) = cmu*cvara_k(iel)*cvara_omg(iel)
    enddo

  endif

  do iel = 1, ncel
    if ( w1(iel).gt.epzero .and.                                  &
         w2(iel).gt.epzero       ) then
      w3(iel) = cebu*w2(iel)/w1(iel)                              &
                   *crom(iel)*volume(iel)                &
                   *(1.d0 - cvara_scal(iel))
      smbrs(iel) = smbrs(iel) - cvara_scal(iel)*w3(iel)
      rovsdt(iel) = rovsdt(iel) + max(w3(iel),zero)
    endif
  enddo

endif

! Free memory
deallocate(w1, w2, w3)

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
