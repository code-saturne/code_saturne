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

subroutine cplphy &
!================

 ( mbrom  , izfppp )

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

! Calcul de RHO de la phase gazeuse


! ATTENTION :
! =========


! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)


!  Il FAUT AVOIR PRECISE icp = 0
!     ==================
!    si on souhaite imposer une chaleur specifique
!    CP variable (sinon: ecrasement memoire).


! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usipsu :
!             . la masse volumique (initialisee a RO0)
!             . la viscosite       (initialisee a VISCL0)
!      - dans usppiv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a la valeur donnee dans usiniv)

! On peut donner ici les lois de variation aux cellules
!     - de la masse volumique                      ROM    kg/m3
!         (et eventuellememt aux faces de bord     ROMB   kg/m3)
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)


! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! mbrom            ! te ! <-- ! indicateur de remplissage de romb              !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          mbrom
integer          izfppp(nfabor)

! Local variables

integer          ntbcpi, ntbcpr
integer          ntbmci, ntbmcr
integer          ntbwoi, ntbwor
integer          iel, icha
integer          izone, ifac
double precision srrom1
double precision wmolme

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5
double precision, allocatable, dimension(:) :: w7, w8
double precision, dimension(:), pointer :: brom,  crom
double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), pointer :: cvar_f1m, cvar_f2m
double precision, dimension(:), pointer :: cvar_f3m, cvar_f4p2m

integer       ipass
data          ipass /0/
save          ipass

!===============================================================================

!===============================================================================
! 0. ON COMPTE LES PASSAGES
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet))
allocate(w7(ncelet), w8(ncelet))

! --- Initialisation memoire


! --- Initialisation des tableaux de travail

do iel = 1, ncel
  w1(iel) = zero
  w2(iel) = zero
  w3(iel) = zero
  w4(iel) = zero
  w5(iel) = zero
  w7(iel) = zero
  w8(iel) = zero
enddo

!===============================================================================
! 2. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE GAZEUSE
!                    VALEURS CELLULES
!                    ----------------
!    TEMPERATURE
!    MASSE VOLUMIQUE
!    CONCENTRATIONS DES ESPECES GAZEUSES
!===============================================================================

! --- Calcul de l'enthalpie du gaz     dans W8
!            de F1M                    dans W2
!            de F2M                    dans W3
!            de F3M                    dans W4
!            de F4M                    dans W5
!            de F4P2M                  dans W7


! ---- W2 = F1M = SOMME(F1M(ICHA))
!      W3 = F2M = SOMME(F2M(ICHA))
!      W4 = F3M
!      W5 = F4M = 1. - F1M - F2M - F3M
!      W7 = F4P2M

call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

do icha = 1, ncharb
  call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m)
  call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m)
  do iel = 1, ncel
    w2(iel) =  w2(iel) + cvar_f1m(iel)
    w3(iel) =  w3(iel) + cvar_f1m(iel)
  enddo
enddo

call field_get_val_s(ivarfl(isca(if3m)), cvar_f3m)
call field_get_val_s(ivarfl(isca(if4p2m)), cvar_f4p2m)
do iel = 1,ncel
  w4(iel) = cvar_f3m(iel)
  w5(iel) = 1.d0 - w2(iel) - w3(iel) -w4(iel)
  w7(iel) = cvar_f4p2m(iel)
  w8(iel) = cvar_scalt(iel)
enddo

! ------ Macro tableau d'entiers TBCPI : NTBCPI
!        Macro tableau de reels  TBCPR : NTBCPR
!        Macro tableau d'entiers TBMCI : NTBMCI
!        Macro tableau de reels  TBMCR : NTBMCR
!        Macro tableau d'entiers TBWOI : NTBWOI
!        Macro tableau de reels  TBWOR : NTBWOR

ntbcpi = 1
ntbcpr = 9
ntbmci = 0
ntbmcr = 2 + 2*ncharb + 4

!  Ce sont en fait X1M, X2M,
!                  F1M(ICHA) et F2M(ICHA) pour chaque charbon
!                  ACHX1F1, ACHX2F2, ACOF1, ACOF2
ntbwoi = 1
ntbwor = 4

call cplph1                                                       &
 ( ncelet , ncel   ,                                              &
   ntbcpi , ntbcpr , ntbmci , ntbmcr , ntbwoi , ntbwor ,          &
   w2     , w3     , w4     , w5     , w7     ,                   &
!  f1m      f2m      f3m      f4m      f4p2m
   w8     ,                                                       &
!  enth
   w1    )
!                          ----
!                 ATTENTION W1 contient RHO1

!===============================================================================
! 3. Relaxation de la masse volumique de la phase gazeuse
!===============================================================================

! --- Calcul de Rho avec relaxation

call field_get_val_s(icrom, crom)

if (ipass.gt.1.or.(isuite.eq.1.and.initro.eq.1)) then
  srrom1 = srrom
else
  srrom1 = 1.d0
endif


do iel = 1, ncel
! ---- Sous relaxation eventuelle a donner dans ppini1.F
  crom(iel) = srrom1*crom(iel)                  &
                     + (1.d0-srrom1)*w1(iel)
enddo


!===============================================================================
! 4. CALCUL DE RHO DE LA PHASE GAZEUSE

!                             VALEURS FACES
!                             -------------
!===============================================================================

mbrom = 1
call field_get_val_s(ibrom, brom)
call field_get_val_s(icrom, crom)

! ---> Masse volumique au bord pour toutes les faces
!      Les faces d'entree seront recalculees.

do ifac = 1, nfabor
  iel = ifabor(ifac)
  brom(ifac) = crom(iel)
enddo

! ---> Masse volumique au bord pour les faces d'entree UNIQUEMENT
!     Le test sur IZONE sert pour les reprises de calcul

if ( ipass.gt.1 .or. isuite.eq.1 ) then
  do ifac = 1, nfabor

    izone = izfppp(ifac)
    if(izone.gt.0) then
      if ( ientat(izone).eq.1 ) then
        wmolme = (1.d0+xsi) / (wmole(io2)+xsi*wmole(in2))
        brom(ifac) = p0                           &
                             /(wmolme*cs_physical_constants_r*timpat(izone))
      endif
    endif

  enddo
endif

! Free memory
deallocate(w1, w2, w3)
deallocate(w4, w5)
deallocate(w7, w8)

!----
! FIN
!----

return
end subroutine
