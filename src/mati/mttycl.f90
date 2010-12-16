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

subroutine mttycl &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   itypfb , itrifb , icodcl , isostd ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl , frcxt  ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

!===============================================================================
! FONCTION :
! --------

! TRAITEMENT DES CODES DE CONDITIONS AUX LIMITES DE MATISSE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! itypfb           ! ia ! <-- ! boundary face types                            !
!  (nfabor, nphas) !    !     !                                                !
! itrifb(nfabor    ! te ! --> ! tab d'indirection pour tri des faces           !
!  nphas)          !    !     !                                                !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! isostd           ! te ! --> ! indicateur de sortie standard                  !
!    (nfabor+1)    !    !     !  +numero de la face de reference               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! frcxt(ncelet,    ! tr ! <-- ! force exterieure generant la pression          !
!   3,nphas)       !    !     !  hydrostatique                                 !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! rijipb           ! tr ! --- ! tab de trav pour valeurs en iprime             !
! (nfabor,6   )    !    !     !  des rij au bord                               !
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
use dimens, only: ndimfb
use numvar
use optcal
use cstnum
use cstphy
use entsor
use parall
use pointe
use matiss
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas

integer          icodcl(nfabor,nvar)
integer          itypfb(nfabor,nphas) , itrifb(nfabor,nphas)
integer          isostd(nfabor+1,nphas)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision rcodcl(nfabor,nvar,3)
double precision frcxt(ncelet,3,nphas)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,3)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          iphas , ifac  , ifml  , icoul
double precision dbm   , roe

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! --- Gestion memoire

idebia = idbia0
idebra = idbra0

! --- Une seule phase

iphas = 1

!===============================================================================
! 2.  BOUCLE SUR LES FACES DE BORD
!===============================================================================

do ifac = 1, nfabor

! --- Couleur de la face de bord

  ifml  = ifmfbr(ifac  )
  icoul = iprfml(ifml,1)


!===============================================================================
! 3.  ENTREE : classique en convection forcee, a pression imposee sinon
!===============================================================================

  if(icoul.eq.icmtfi) then


! --- Convection forcee
    if(icofor.eq.1)then

!     Type de base : IENTRE
      itypfb(ifac,iphas) = ientre

!     Temperature exterieure prise a TINIT
      rcodcl(ifac,isca(itaamt),1) = tinit
      rcodcl(ifac,isca(itpcmt),1) = tinit
      rcodcl(ifac,isca(itppmt),1) = tinit

!     Vitesses transverses nulles
      rcodcl(ifac,iu(iphas),1) = 0.d0
      rcodcl(ifac,iv(iphas),1) = 0.d0

!     Vitesse debitante W est calculee a partir du debit massique par
!       W = - debit massique / (masse volumique * surface)
!     . On suppose la vitesse verticale descendante (d'ou le signe).
!     . ROE est la masse volumique exterieure evaluee a partir de
!       la temperature TINIT (degres C) en utilisant la loi d'etat.
!     . La surface de l'entree du maillage est calculee par la formule
!       NPTRAN*PTRRES*EPCHEM/RCONVE. En effet, la surface au bas de la
!       cheminee est NPTRAN*PTRRES*EPCHEM et on divise par le rapport
!       du convergent represente sur le maillage (>=1) pour obtenir
!       la surface d'entree. (en 2D, il faut RCONVE=1, sinon, on impose
!       un debit plus grand que le debit souhaite).
!     . La correction du debit massique par FRDTRA correspond au
!       rapport d'echelle transverse eventuel entre le maillage et la
!       realite    .

      roe = rrfmat*(trfmat+tkelvi)/(tinit+tkelvi)
      dbm = debmas/frdtra
      rcodcl(ifac,iw(iphas),1) =                                  &
             -dbm*rconve/(roe*nptran*ptrres*epchem)


! --- Convection naturelle (ICOFOR=0)
    else

!     Type de base indefini
!       (dans une version ulterieure de Code_Saturne, on disposera de
!        sorties a pression imposee)
      itypfb(ifac,iphas) = iindef

!     Dirichlet sur les temperatures
      icodcl(ifac,isca(itaamt)  ) = 1
      rcodcl(ifac,isca(itaamt),1) = tinit
      icodcl(ifac,isca(itpcmt)  ) = 1
      rcodcl(ifac,isca(itpcmt),1) = tinit
      icodcl(ifac,isca(itppmt)  ) = 1
      rcodcl(ifac,isca(itppmt),1) = tinit

!     Dirichlet nul sur les vitesses transverses (nul par symetrie)
      icodcl(ifac,iu(iphas)  ) = 1
      rcodcl(ifac,iu(iphas),1) = 0.d0
      icodcl(ifac,iv(iphas)  ) = 1
      rcodcl(ifac,iv(iphas),1) = 0.d0

!     Neumann homogene sur la vitesse debitante (on impose la pression)
      icodcl(ifac,iw(iphas)  ) = 3

!     Pression nulle
      icodcl(ifac,ipr(iphas)  ) = 1
      rcodcl(ifac,ipr(iphas),1) = p0(iphas)

    endif


!===============================================================================
! 4.  SORTIE (pression imposee)
!===============================================================================

  elseif(icoul.eq.icmtfo) then

!     Type de base : sortie libre
    itypfb(ifac,iphas) = isolib

!     Temperature imposee en cas de reentree
    rcodcl(ifac,isca(itaamt),1) = tcrit
    rcodcl(ifac,isca(itpcmt),1) = tinit
    rcodcl(ifac,isca(itppmt),1) = tinit

!     Pression imposee (ecart hydrostatique par rapport a l'entree
!                       + decalage eventuel DPVENT)
    icodcl(ifac,ipr(iphas)  ) = 1
    rcodcl(ifac,ipr(iphas),1) = p0(iphas) + dpvent                &
         - dabs(gz) * ro0(iphas) * (hcheva - hchali)

!===============================================================================
! 5.  LE RESTE : glissement et Neumann homogene
!===============================================================================

! --- Sol
  elseif(icoul.eq.icmtfg) then
    itypfb(ifac,iphas)   = isymet

! --- Plafond
  elseif(icoul.eq.icmtfc) then
    itypfb(ifac,iphas)   = isymet

! --- Symetries
  elseif(icoul.eq.icmtfs) then
    itypfb(ifac,iphas)   = isymet

! --- Parois
  elseif(icoul.eq.icmtfw) then
    itypfb(ifac,iphas)   = isymet

  endif

enddo


!----
! FORMATS
!----

!----
! FIN
!----

return

end subroutine
