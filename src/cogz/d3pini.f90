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

subroutine d3pini &
!================

 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL
!    POUR LA PHYSIQUE PARTICULIERE :
!        COMBUSTION GAZ - FLAMME DE DIFFUSION CHIMIE 3 POINTS
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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
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
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use field

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

character(len=80) :: chaine
integer           iel, igg, mode
integer           iscal, ivar, ii
double precision coefg(ngazgm), hair, tinitk
double precision valmax, valmin

double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m
double precision, dimension(:), pointer :: cvar_npm, cvar_fsm
double precision, dimension(:), pointer :: cvar_scal

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

ipass = ipass + 1

if (ippmod(icod3p).eq.1) call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
call field_get_val_s(ivarfl(isca(ifm)), cvar_fm)
call field_get_val_s(ivarfl(isca(ifp2m)), cvar_fp2m)
if (isoot.eq.1) then
  call field_get_val_s(ivarfl(isca(inpm)), cvar_npm)
  call field_get_val_s(ivarfl(isca(ifsm)), cvar_fsm)
endif

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

if ( isuite.eq.0 ) then

! ---> Initialisation au 1er passage avec de l'air a TINITK
!                                    ======================

  if ( ipass.eq.1 ) then

! ----- Calcul de l'enthalpie de l'air HAIR a TINITK

    tinitk   = t0
    coefg(1) = zero
    coefg(2) = 1.d0
    coefg(3) = zero
    mode     = -1
    call cothht                                                   &
    !==========
  ( mode   , ngazg , ngazgm  , coefg  ,                           &
    npo    , npot   , th     , ehgazg ,                           &
    hair   , tinitk )


    do iel = 1, ncel

! ----- Moyenne et variance du taux de melange

      cvar_fm(iel)   = zero
      cvar_fp2m(iel) = zero

! ----- Enthalpie

      if ( ippmod(icod3p).eq.1 ) then
        cvar_scalt(iel) = hair
      endif

! ---- Soot
      if (isoot.eq.1) then
        cvar_npm(iel) = 0.d0
        cvar_fsm(iel) = 0.d0
      endif
    enddo

! ---> Initialisation au 2eme passage

  else if ( ipass.eq.2 ) then

    do iel = 1, ncel

! ----- Moyenne et variance du taux de melange

      cvar_fm(iel)   = fs(1)
      cvar_fp2m(iel) = zero

! ----- Enthalpie

      if ( ippmod(icod3p).eq.1 ) then
        cvar_scalt(iel) = hinfue*fs(1)+hinoxy*(1.d0-fs(1))
      endif

! ---- Soot
      if (isoot.eq.1) then
        cvar_npm(iel) = 0.d0
        cvar_fsm(iel) = 0.d0
      endif


    enddo

! ----- On donne la main a l'utilisateur

    call cs_user_initialization &
    !==========================
  ( nvar   , nscal  ,                                            &
    dt     )

! ----- En periodique et en parallele,
!       il faut echanger ces initialisations

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(cvar_fm)
      !==========
      call synsca(cvar_fp2m)
      !==========
      if ( ippmod(icod3p).eq.1 ) then
        call synsca(cvar_scalt)
        !==========
      endif
    endif

      ! ---- Soot
      if (isoot.eq.1) then
        call synsca(cvar_npm)
        call synsca(cvar_fsm)
      endif


!      Impressions de controle

    write(nfecra,2000)

    do ii  = 1, nscapp
      iscal = iscapp(ii)
      ivar  = isca(iscal)
      call field_get_val_s(ivarfl(isca(iscal)), cvar_scal)
      valmax = -grand
      valmin =  grand
      do iel = 1, ncel
        valmax = max(valmax,cvar_scal(iel))
        valmin = min(valmin,cvar_scal(iel))
      enddo
      call field_get_label(ivarfl(ivar), chaine)
      if (irangp.ge.0) then
        call parmin(valmin)
        !==========
        call parmax(valmax)
        !==========
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
' ** INITIALISATION DES VARIABLES PROPRES AU GAZ (FL DIF 3PT) ',/,&
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
! Fin
!----

return
end subroutine
