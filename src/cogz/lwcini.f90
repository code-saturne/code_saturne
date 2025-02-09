!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

subroutine lwcini () &
  bind(C, name='cs_f_lwcini')

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL
!    POUR LA PHYSIQUE PARTICULIERE : COMBUSTION GAZ MODELE LWC
!    PENDANT DE USINIV.F

! Cette routine est appelee en debut de calcul (suite ou non)
!     avant le debut de la boucle en temps

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul,
!     les valeurs du pas de temps


! On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
!     ou relues d'un fichier suite
! On ne dispose des variables VISCLS, CP (quand elles sont
!     definies) que si elles ont pu etre relues dans un fichier
!     suite de calcul

! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
!     ET PAS ICI

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppthch
use coincl
use ppincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

character(len=80) :: chaine
integer          iel, igg
integer          f_id, ii
double precision hinit, coefg(ngazgm), hair, tinitk
double precision tentm, fmelm
double precision valmax, valmin

double precision, dimension(:), pointer :: cvar_yfm, cvar_yfp2m
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m
double precision, dimension(:), pointer :: cvar_coyfp
double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), pointer :: cvar_scal

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

interface

  subroutine cs_combustion_boundary_conditions_mean_inlet_ebu_lw  &
    (fmm, tkm)                                                           &
    bind(C, name='cs_combustion_boundary_conditions_mean_inlet_ebu_lw')
    use, intrinsic :: iso_c_binding
    implicit none
    real(kind=c_double), intent(out) :: fmm, tkm
  end subroutine cs_combustion_boundary_conditions_mean_inlet_ebu_lw

end interface

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

ipass = ipass + 1

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

call field_get_val_s(iyfm, cvar_yfm)
call field_get_val_s(iyfp2m, cvar_yfp2m)
call field_get_val_s(ifm, cvar_fm)
call field_get_val_s(ifp2m, cvar_fp2m)
if (ippmod(icolwc).ge.2) call field_get_val_s(icoyfp, cvar_coyfp)
if (ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3 .or.      &
    ippmod(icolwc).eq.5) then
  call field_get_val_s(ihm, cvar_scalt)
endif

!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

if ( isuite.eq.0 ) then

! ---> Initialisation au 1er passage avec de l'air a TINITK
!                                    ======================

  if ( ipass.eq.1 ) then

! ----- Temperature du melange : air a TINITK
    tinitk = t0

! ----- Enthalpie de l'air a TINITK
    if ( ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3 .or.        &
         ippmod(icolwc).eq.5  ) then
      coefg(1) = zero
      coefg(2) = 1.d0
      coefg(3) = zero
      hair = cs_gas_combustion_t_to_h(coefg, tinitk)
    endif

! ----- On en profite pour initialiser FRMEL et TGF
!       CAR on n'a pas encore vu usebuc.F

    frmel = zero
    tgf   = 300.d0

    do iel = 1, ncel

! ----- Fraction massique de fuel et sa variance

      cvar_yfm(iel) = fmax
      cvar_yfp2m(iel) = zero

! ----- Fraction de melange et sa variance

      cvar_fm(iel)   = fmax
      cvar_fp2m(iel) = zero

      if ( ippmod(icolwc).ge. 2) then
        cvar_coyfp(iel)   = zero
      endif

! ----- Enthalpie du melange

      if ( ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3 .or.      &
          ippmod(icolwc).eq.5) then
        cvar_scalt(iel) = hair
      endif

    enddo

! ---> Initialisation au 2eme passage

  else if ( ipass.eq.2 ) then

    call cs_combustion_boundary_conditions_mean_inlet_ebu_lw(fmelm, tentm)

! ----- Enthalpie du melange HINIT
    if ( ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3 .or.        &
        ippmod(icolwc).eq.5  ) then
      coefg(1) = fmelm
      coefg(2) = (1.d0-fmelm)
      coefg(3) = zero
      hinit = cs_gas_combustion_t_to_h(coefg, tentm)
    endif

! ----- En periodique et en parallele,
!       il faut echanger ces initialisations

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(cvar_yfm)
      call synsca(cvar_yfp2m)
      call synsca(cvar_fm)
      call synsca(cvar_fp2m)

      if ( ippmod(icolwc).ge.2) then
        call synsca(cvar_coyfp)
      endif

      if (ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3 .or.       &
          ippmod(icolwc).eq.5 ) then
        call synsca(cvar_scalt)
      endif

    endif

!      Impressions de controle

    write(nfecra,2000)

    do ii  = 1, nscapp
      f_id = ivarfl(isca(iscapp(ii)))
      call field_get_val_s(f_id, cvar_scal)
      valmax = -grand
      valmin =  grand
      do iel = 1, ncel
        valmax = max(valmax,cvar_scal(iel))
        valmin = min(valmin,cvar_scal(iel))
      enddo
      call field_get_label(f_id, chaine)
      if (irangp.ge.0) then
        call parmin(valmin)
        call parmax(valmax)
      endif
      write(nfecra,2010)chaine(1:8),valmin,valmax
    enddo

    write(nfecra,2020)

  endif

endif

!----
! FORMATS
!----


 2000 format(                                                           &
'                                                             ',/,&
' ----------------------------------------------------------- ',/,&
'                                                             ',/,&
'                                                             ',/,&
' ** INITIALISATION DES VARIABLES PROPRES AU GAZ (FL PRE LWC) ',/,&
'    -------------------------------------------------------- ',/,&
'           2eme PASSAGE                                      ',/,&
' ---------------------------------                           ',/,&
'  Variable  Valeur min  Valeur max                           ',/,&
' ---------------------------------                           '  )

 2010 format(                                                           &
 2x,     a8,      e12.4,      e12.4                              )

 2020 format(                                                           &
' ---------------------------------                           ',/)

!----
! FIN
!----

return
end subroutine
