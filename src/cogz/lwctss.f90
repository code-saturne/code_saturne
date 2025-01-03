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

subroutine lwctss &
 (f_id, smbrs, rovsdt)

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME PREMELANGE MODELE LWC
!   ON PRECISE LES TERMES SOURCES POUR UN SCALAIRE PP
!   SUR UN PAS DE TEMPS

! ATTENTION : LE TRAITEMENT DES TERMES SOURCES EST DIFFERENT
! ---------   DE CELUI DE cs_user_source_terms.

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
! f_id             ! i  ! <-- ! field id (scalar variable)                     !
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
use coincl
use ppincl
use mesh
use field
use field_operator
use pointe

!===============================================================================

implicit none

! Arguments

integer          f_id

double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

integer          iel, idirac
integer          inc , iprev
integer          krvarfl

double precision sum, epsi
double precision tsgrad, tschim, tsdiss
double precision turb_schmidt, rvarfl

double precision, allocatable, dimension(:,:) :: gradf, grady
double precision, allocatable, dimension(:) :: w10, w11
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: visct
double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_omg
double precision, dimension(:,:), pointer :: cvara_rij
double precision, dimension(:), pointer :: cvara_scal
double precision, dimension(:), pointer :: cvara_yfm, cvara_fm
type(pmapper_double_r1), dimension(:), pointer :: cpro_fmel, cpro_fmal
type(pmapper_double_r1), dimension(:), pointer :: cpro_tscl, cpro_rhol

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

epsi   = 1.0d-10

call field_get_key_id("variance_dissipation", krvarfl)
call field_get_key_double(f_id, krvarfl, rvarfl)

! ---
call field_get_val_s(icrom, crom)
call field_get_val_s(ivisct, visct)

call field_get_val_prev_s(f_id, cvara_scal)
call field_get_val_prev_s(iyfm, cvara_yfm)
call field_get_val_prev_s(ifm, cvara_fm)

if (itytur.eq.2.or.iturb.eq.50) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
elseif (itytur.eq.3) then
  call field_get_val_prev_v(ivarfl(irij), cvara_rij)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
elseif (iturb.eq.60) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iomg), cvara_omg)
endif

allocate(cpro_fmel(ndirac))
allocate(cpro_fmal(ndirac))
allocate(cpro_tscl(ndirac))
allocate(cpro_rhol(ndirac))

do idirac = 1, ndirac
  call field_get_val_s(ifmel(idirac), cpro_fmel(idirac)%p)
  call field_get_val_s(ifmal(idirac), cpro_fmal(idirac)%p)
  call field_get_val_s(irhol(idirac), cpro_rhol(idirac)%p)
  call field_get_val_s(itscl(idirac), cpro_tscl(idirac)%p)
enddo

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES
!===============================================================================

if (f_id.eq.iyfm) then

! ---> Terme source pour la fraction massique moyenne de fuel

  do iel = 1, ncel
      sum = zero
      do idirac = 1, ndirac
        sum  = sum + cpro_rhol(idirac)%p(iel)                   &
           *cpro_tscl(idirac)%p(iel)*volume(iel)
      enddo

! terme implicite

      if (cvara_scal(iel).gt.epsi) then
        rovsdt(iel) = rovsdt(iel) + max(-sum/cvara_scal(iel),zero)
      endif

! terme explicite

       smbrs(iel) =  smbrs(iel) + sum

  enddo

endif

! ---> Terme source pour la variance de la fraction massique moyenne de fuel

if (f_id.eq.iyfp2m) then

  do iel = 1, ncel
    sum = zero
    do idirac = 1, ndirac
      sum  = sum + (cpro_tscl(idirac)%p(iel)*volume(iel)        &
        *(cpro_fmal(idirac)%p(iel) - cvara_yfm(iel))           &
             *cpro_rhol(idirac)%p(iel))
    enddo
    smbrs(iel) = smbrs(iel) + sum
  enddo

endif

! ---> Terme source pour la covariance

if (f_id.eq.icoyfp) then

  ! Allocate a temporary array for gradient computation
  allocate(gradf(3,ncelet), grady(3,ncelet))

  ! Allocate work arrays
  allocate(w10(ncelet), w11(ncelet))

! --- Calcul du gradient de F
!     =======================

  do iel = 1, ncel
    w10(iel) = cvara_fm(iel)
  enddo

  iprev = 1
  inc = 1

  call field_gradient_scalar(ifm, iprev, inc, gradf)

! --- Calcul du gradient de Yfuel
!     ===========================

  do iel = 1, ncel
    w11(iel) = cvara_yfm(iel)
  enddo

  iprev = 1
  inc = 1

  call field_gradient_scalar(iyfm, iprev, inc, grady)

! --- Calcul du terme source
!     ======================

! ---> Calcul de K et Epsilon en fonction du modele de turbulence

! ---- TURBULENCE

  if (itytur.eq.2) then

    do iel = 1, ncel
      w10(iel) = cvara_k(iel)
      w11(iel) = cvara_ep(iel)
    enddo

  elseif (itytur.eq.3) then

    do iel = 1, ncel
      w10(iel) = ( cvara_rij(1,iel)                          &
                  +cvara_rij(2,iel)                          &
                  +cvara_rij(3,iel) ) / 2.d0
      w11(iel) = cvara_ep(iel)
    enddo

  elseif (iturb.eq.50) then

    do iel = 1, ncel
      w10(iel) = cvara_k(iel)
      w11(iel) = cvara_ep(iel)
    enddo

  elseif (iturb.eq.60) then

    do iel = 1, ncel
      w10(iel) = cvara_k(iel)
      w11(iel) = cmu*cvara_k(iel)*cvara_omg(iel)
    enddo

  endif

  call field_get_key_double(f_id, ksigmas, turb_schmidt)

  do iel=1,ncel

!  A confirmer :
!   Le terme de dissipation devrait etre implicite
!   Dans le terme de dissipation, il manque une constante Cf
!   Peut-elle etre consideree egale a 1 ?
!   Verifier le signe du terme de production
!-
! terme implicite


    w11(iel) = w11(iel)/(w10(iel)*rvarfl)*volume(iel)*crom(iel)
    rovsdt(iel) = rovsdt(iel) + max(w11(iel),zero)

! terme de gradient

    tsgrad =  (2.0d0                                              &
         * visct(iel)/(turb_schmidt)                              &
         * (  gradf(1,iel)*grady(1,iel)                           &
            + gradf(2,iel)*grady(2,iel)                           &
            + gradf(3,iel)*grady(3,iel) ))                        &
         * volume(iel)


! terme de dissipation

    tsdiss = -w11(iel) * cvara_scal(iel)

! terme de chimique

    tschim = zero
    do idirac = 1, ndirac
      tschim =   tschim                                         &
               + (cpro_tscl(idirac)%p(iel)                      &
                 *(cpro_fmel(idirac)%p(iel)-cvara_fm(iel))      &
                  *volume(iel))*cpro_rhol(idirac)%p(iel)
    enddo

! --> Somme des termes

    smbrs(iel) = smbrs(iel) + tschim + tsgrad + tsdiss

  enddo

  ! Free memory
  deallocate(gradf, grady)
  deallocate(w10, w11)

endif

deallocate(cpro_fmel, cpro_fmal)
deallocate(cpro_tscl, cpro_rhol)

!----
! FIN
!----

return

end subroutine
