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

subroutine uscfth &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   iccfth , imodif , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   sorti1 , sorti2 , gamagr , xmasmr ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ---------

! ROUTINE UTILISATEUR POUR ENTREE DES PARAMETRES THERMODYNAMIQUES


!    CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE
!    =============================================



! CALCUL DE LA PROPRIETE THERMODYNAMIQUE VOULUE
! EN FONCTION DES PARAMETRES DONNES ET DE LA THERMO CHOISIE


! LOIS THERMODYNAMIQUES IMPLEMENTEES :
! ====================================

!  1. GAZ PARFAIT (fournir la masse molaire XMASML)

!  2. GAZ PARFAIT A GAMMA VARIABLE (à adapter)

!  3. VAN DER WAALS (pas encore implemente)


! CALCULS IMPLEMENTES :
! =====================

!  VARIABLE     P   rho   T    e    h    s  epsilon-CvT  C.L.
!      CODE     1    2    3    4    5    6      7         9


!  CODE DU CALCUL : *VARIABLE i : i

!                                      di|
!                   *DERIVEE PARTIELLE --| : ijk
!                                      dj|k

!                   *CALCUL DES VARIABLES A PARTIR DE i ET j : ij

!  Et pour un reperage automatique, on utilise les nombres premiers
!           P   rho   T   e   s
!           2     3   5   7   13  (* 10000 pour les calculs aux cellules)



!  -> OPTIONS DE CALCUL                           : ICCFTH = -1

!  -> INITIALISATIONS PAR DEFAUT                  : ICCFTH = 0

!  -> CALCUL DE GAMMA                             : ICCFTH = 1

!  -> VERIFICATION DE rho                         : ICCFTH = -2

!  -> VERIFICATION DE E                           : ICCFTH = -4

!  -> CALCUL DE T ET e EN FONCTION DE P ET rho    : ICCFTH = 12 ou  60000

!  -> CALCUL DE rho ET e EN FONCTION DE P ET T    : ICCFTH = 13 ou 100000

!  -> CALCUL DE rho ET T EN FONCTION DE P ET e    : ICCFTH = 14 ou 140000

!  -> CALCUL DE P ET e EN FONCTION DE rho ET T    : ICCFTH = 23 ou 150000

!  -> CALCUL DE P ET T EN FONCTION DE rho ET e    : ICCFTH = 24 ou 210000

!                2    dP |
!  -> CALCUL DE c  = ----|                        : ICCFTH = 126
!                    drho|s

!                      dP|
!  -> CALCUL DE beta = --|                        : ICCFTH = 162
!                      ds|rho

!                    de|
!  -> CALCUL DE Cv = --|                          : ICCFTH = 432
!                    dT|rho

!  -> CALCUL DE L'ENTROPIE                        : ICCFTH = 6


!  -> CALCUL DE epsilon - Cv.T                    : ICCFTH = 7


!  -> CALCUL DES CONDITIONS AUX LIMITES
!              - SYMETRIE                         : ICCFTH = 90
!              - PAROI                            : ICCFTH = 91
!              - ENTREE                           : ICCFTH = 92
!              - SORTIE                           : ICCFTH = 93
!              - T ET e EN FONCTION DE P ET rho   : ICCFTH = 912 ou  60900
!              - rho ET e EN FONCTION DE P ET T   : ICCFTH = 913 ou 100900
!              - rho ET T EN FONCTION DE P ET e   : ICCFTH = 914 ou 140900
!              - P ET e EN FONCTION DE rho ET T   : ICCFTH = 923 ou 150900
!              - P ET T EN FONCTION DE rho ET e   : ICCFTH = 924 ou 210900


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! iccfth           ! e  ! <-- ! numero du type de calcul demande               !
! imodif           ! e  ! <-- ! si calcul aux cellules : modif de rtp          !
!                  !    !     !   quand imodif > 0 (en particulier             !
!                  !    !     !   automatise l'init (voir uscfxi)              !
!                  !    !     ! si calcul aux faces : numero de face           !
!                  !    !     !   automatise l'init (voir uscfxi)              !
! iphas            ! i  ! <-- ! phase number                                   !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp,rtpa         ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa coefb      ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! sorti1,2(*)      ! tr ! --> ! variables de sortie                            !
!                  !    !     !   (inutilise si iccfth.lt.0)                   !
! gamagr(*)        ! tr ! --> ! constante gamma equivalent du gaz              !
!                  !    !     !   (inutilise si iccfth.lt.0)                   !
!                  !    !     !   (premiere case utilisee en g.p.)             !
! xmasmr(*)        ! tr ! --> ! masse molaire des constituants du gaz          !
!                  !    !     !   (inutilise si iccfth.lt.0)                   !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "pointe.h"
include "entsor.h"
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          iccfth   , imodif , iphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*),propfa(nfac,*),propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision sorti1(*), sorti2(*), gamagr(*), xmasmr(*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          iiph   , ifac0
integer          ierr
integer          iel    , ifac   , ivar
integer          ipriph , irhiph , itkiph , ieniph
integer          iuiph  , iviph  , iwiph
integer          iclp   , iclr   , iclt   , icle
integer          iclu   , iclv   , iclw
integer          iutile
double precision gamagp , xmasml , enint
double precision xmach  , xmachi , xmache , dxmach

integer          npmax
parameter (npmax = 1000)
double precision cstgr(npmax)

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0. INITIALISATION
!     Pas d'intervention utilisateur requise.
!===============================================================================

!     Pointeurs mémoire
idebia = idbia0
idebra = idbra0

!     Indicateur d'erreur (stop si non nul)
ierr   = 0

!     Position des variables
if(iccfth.ge.0.or.iccfth.le.-2) then
  ipriph = ipr(iphas)
  irhiph = isca(irho  (iphas))
  itkiph = isca(itempk(iphas))
  ieniph = isca(ienerg(iphas))
  iuiph = iu(iphas)
  iviph = iv(iphas)
  iwiph = iw(iphas)
  iclp = iclrtp(ipriph,icoef)
  iclr = iclrtp(irhiph,icoef)
  iclt = iclrtp(itkiph,icoef)
  icle = iclrtp(ieniph,icoef)
  iclu = iclrtp(iuiph,icoef)
  iclv = iclrtp(iviph,icoef)
  iclw = iclrtp(iwiph,icoef)
endif


!     Si calcul aux cellules,
!       indique, si > 0, que RTP doit etre modifié
!     Si calcul aux faces,
!       indique le numéro de la face considérée
ifac0 = imodif

!===============================================================================
! 1. CHOIX DE LA THERMODYNAMIQUE POUR CHAQUE PHASE
!     Intervention utilisateur requise.
!===============================================================================

! --> IEOS = 1 : Gaz Parfait avec Gamma constant
! --> IEOS = 2 : Gaz Parfait avec Gamma variable
! --> IEOS = 3 : Equation d'etat de Van Der Waals

do iiph = 1, nphas
  ieos(iiph) = 1
enddo


!            NE PAS OUBLIER DE RENSEIGNER LES PARAMETRES
!                  DANS LA THERMODYNAMIQUE CHOISIE
!            ===========================================


!===============================================================================
! 2.  GAZ PARFAIT
!===============================================================================

if(ieos(iphas).eq.1) then

!===============================================================================
! 2.1. PARAMETRES A RENSEIGNER PAR L'UTILISATEUR
!     Intervention utilisateur requise.
!===============================================================================

!     Pour chaque phase

!     Masse molaire du gaz (en kg/mol)

  if(iccfth.ge.0) then

    if(iphas.eq.1) then

      xmasml = 28.8d-3

    endif

  endif

!===============================================================================
! 2.2. LOIS IMPLEMENTEES PAR DEFAUT
!     Pas d'intervention utilisateur requise.
!===============================================================================

!     CALCUL DE LA CONSTANTE GAMAGP
!     =============================

!     On la suppose supérieure ou égale à 1.
!     On la calcule à chaque appel, même si cela peut paraître cher, dans
!       par cohérence avec le cas gamma variable, pour lequel, on n'a pas
!       prévu de la conserver. Avec un save et un test, on pourrait
!       s'affranchir du calcul si nécessaire.

  if(iccfth.gt.0) then

    gamagp = 1.d0 + rr/(xmasml*cp0(iphas)-rr)

    if(gamagp.lt.1.d0) then
      write(nfecra,1010) gamagp
      call csexit (1)
    endif

!     On renvoie gamma si demande
    if(iccfth.eq.1) then
      gamagr(1) = gamagp
    endif

  endif


!     OPTIONS DE CALCUL  : Cp ET Cv CONSTANTS (hypothèse gaz parfait)
!     =================

!     La valeur de CP0 est à fournir dans usini1. La valeur de CV0
!       est calculée plus bas (à partir de CP0).


  if(iccfth.eq.-1) then

    icp(iphas) = 0
    icv(iphas) = 0


!     INITIALISATIONS PAR DEFAUT (avant uscfxi)
!     ==========================
!     T0 est forcément positif (vérifié dans verini)

  elseif(iccfth.eq.0) then

    cv0(iphas) = cp0(iphas) - rr/xmasml

    if ( isuite .eq. 0 ) then
      do iel = 1, ncel
        rtp(iel,irhiph) = p0(iphas)*xmasml/(rr*t0(iphas))
        rtp(iel,ieniph) = cv0(iphas)*t0(iphas)
      enddo
    endif


!     VERIFICATION DE rho
!     ===================

  elseif(iccfth.eq.-2) then


! --- Clipping si rho est exactement nul + write + stop
!       En effet, si c'est le cas, on va nécessairement se planter
!       dans la thermo (bon, d'accord, pas dans tous les cas, mais
!       pour le moment, je présume que les fluides à traiter seront
!       assez classiques, avec des pbs à rho = 0)
!       Appelé en fin de resolution de la masse volumique (apres
!       clipping classique et avant communication parallele).

    ierr = 0
    do iel = 1, ncel
      if(rtp(iel,irhiph).le.0.d0) then
        rtp(iel,irhiph) = epzero
        ierr = ierr + 1
      endif
    enddo
    if(irangp.ge.0) then
      call parcpt (ierr)
    endif
    if(ierr.gt.0) then
      ntmabs = ntcabs
      write(nfecra,8000)ierr, epzero
    endif


!     VERIFICATION DE e
!     ===================

  elseif(iccfth.eq.-4) then


! --- Clipping si e est exactement nul + write + stop
!       En effet, si c'est le cas, on va nécessairement se planter
!       dans la thermo (bon, d'accord, pas dans tous les cas, mais
!       pour le moment, je présume que les fluides à traiter seront
!       assez classiques, avec des pbs à e = 0)

    ierr = 0
    do iel = 1, ncel
      enint = rtp(iel,ieniph)                                     &
               - 0.5d0*( rtp(iel,iuiph)**2                        &
                       + rtp(iel,iviph)**2                        &
                       + rtp(iel,iwiph)**2 )
      if(enint.le.0.d0) then
        rtp(iel,ieniph) = epzero                                  &
               + 0.5d0*( rtp(iel,iuiph)**2                        &
                       + rtp(iel,iviph)**2                        &
                       + rtp(iel,iwiph)**2 )
        ierr = ierr + 1
      endif
    enddo
    if(irangp.ge.0) then
      call parcpt (ierr)
    endif
    if(ierr.gt.0) then
      ntmabs = ntcabs
      write(nfecra,8100)ierr, epzero
    endif


!     CALCUL DE T ET e EN FONCTION DE P ET rho :
!     ========================================

  elseif(iccfth.eq.12.or.iccfth.eq.60000) then

!     Test des valeurs de rho
    ierr = 0
    do iel = 1, ncel
      if(rtp(iel,irhiph).le.0.d0) then
        write(nfecra,3010)rtp(iel,irhiph),iel
      endif
    enddo
!     Si pb on s'arrete, car rho a été fourni par l'utilisateur
!       (initialisation erronée sans doute)
    if(ierr.eq.1) then
      call csexit (1)
    endif

    do iel = 1, ncel
!     Temperature
      sorti1(iel) = xmasml*rtp(iel,ipriph)/(rr*rtp(iel,irhiph))
!     Energie totale
      sorti2(iel) = cv0(iphas)*sorti1(iel)                        &
           + 0.5d0*( rtp(iel,iuiph)**2 + rtp(iel,iviph)**2        &
                                       + rtp(iel,iwiph)**2 )
    enddo

!     Affectation a RTP
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,itkiph) = sorti1(iel)
        rtp(iel,ieniph) = sorti2(iel)
      enddo
    endif


!     CALCUL DE rho ET e EN FONCTION DE P ET T :
!     ========================================

  elseif(iccfth.eq.13.or.iccfth.eq.100000) then

!     Test des valeurs de T
    ierr = 0
    do iel = 1, ncel
      if(rtp(iel,itkiph).le.0.d0) then
        write(nfecra,2010)rtp(iel,itkiph),iel
      endif
    enddo
!     Si pb on s'arrete, car T a été fourni par l'utilisateur
!       (initialisation erronée par exemple, non en K ?)
    if(ierr.eq.1) then
      call csexit (1)
    endif

    do iel = 1, ncel
!     Masse volumique
      sorti1(iel) = xmasml*rtp(iel,ipriph)/(rr*rtp(iel,itkiph))
!     Energie totale
      sorti2(iel) = cv0(iphas)*rtp(iel,itkiph)                    &
           + 0.5d0*( rtp(iel,iuiph)**2 + rtp(iel,iviph)**2        &
                                       + rtp(iel,iwiph)**2 )
    enddo

!     Affectation a RTP
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,irhiph) = sorti1(iel)
        rtp(iel,ieniph) = sorti2(iel)
      enddo
    endif


!     CALCUL DE rho ET T EN FONCTION DE P ET e :
!     ========================================

  elseif(iccfth.eq.14.or.iccfth.eq.140000) then

    do iel = 1, ncel
!     Energie interne (évite de diviser par T)
      enint = rtp(iel,ieniph)                                     &
               - 0.5d0*( rtp(iel,iuiph)**2                        &
                       + rtp(iel,iviph)**2                        &
                       + rtp(iel,iwiph)**2 )
!     Masse volumique
      sorti1(iel) = rtp(iel,ipriph) / ( (gamagp-1.d0) * enint )
!     Temperature
      sorti2(iel) = xmasml * (gamagp-1.d0) * enint / rr
    enddo

!     Affectation a RTP
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,irhiph) = sorti1(iel)
        rtp(iel,itkiph) = sorti2(iel)
      enddo
    endif


!     CALCUL DE P ET e EN FONCTION DE rho ET T :
!     ========================================

  elseif(iccfth.eq.23.or.iccfth.eq.150000) then

    do iel = 1, ncel
!     Pression
      sorti1(iel) = rtp(iel,irhiph)*rtp(iel,itkiph)*rr/xmasml
!     Energie totale
      sorti2(iel) = cv0(iphas)*rtp(iel,itkiph)                    &
           + 0.5d0*( rtp(iel,iuiph)**2 + rtp(iel,iviph)**2        &
                                       + rtp(iel,iwiph)**2 )
    enddo

!     Affectation a RTP
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,ipriph) = sorti1(iel)
        rtp(iel,ieniph) = sorti2(iel)
      enddo
    endif


!     CALCUL DE P ET T EN FONCTION DE rho ET e :
!     ========================================

  elseif(iccfth.eq.24.or.iccfth.eq.210000) then

    do iel = 1, ncel
!     Energie interne (évite de diviser par T)
      enint = rtp(iel,ieniph)                                     &
               - 0.5d0*( rtp(iel,iuiph)**2                        &
                       + rtp(iel,iviph)**2                        &
                       + rtp(iel,iwiph)**2 )
!     Pression
      sorti1(iel) = (gamagp-1.d0) * rtp(iel,irhiph) * enint
!     Temperature
      sorti2(iel) = xmasml * (gamagp-1.d0) * enint / rr
    enddo

!       Affectation a RTP
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,ipriph) = sorti1(iel)
        rtp(iel,itkiph) = sorti2(iel)
      enddo
    endif

!                2                            2          P
!     CALCUL DE c  EN FONCTION DE P ET rho : c  = gamma*---
!     ====================================              rho

  elseif(iccfth.eq.126) then

!     Test des valeurs de rho (on peut éliminer ce test pour optimiser)
    iutile = 0
    if(iutile.eq.1) then
      ierr = 0
      do iel = 1, ncel
        if(rtp(iel,irhiph).le.0.d0) then
          write(nfecra,4010)rtp(iel,irhiph),iel
        endif
      enddo
      if(ierr.eq.1) then
        call csexit (1)
      endif
    endif

    do iel = 1, ncel
      sorti1(iel) = gamagp * rtp(iel,ipriph) / rtp(iel,irhiph)
    enddo

!                                                    gamma
! CALCUL DE beta EN FONCTION DE P ET rho : beta = rho
! ======================================

  elseif(iccfth.eq.162) then

!     Test des valeurs de rho (on peut éliminer ce test pour optimiser)
    iutile = 0
    if(iutile.eq.1) then
      ierr = 0
      do iel = 1, ncel
        if(rtp(iel,irhiph).lt.0.d0) then
          write(nfecra,4020)rtp(iel,irhiph),iel
        endif
      enddo
      if(ierr.eq.1) then
        call csexit (1)
      endif
    endif

    do iel = 1, ncel
      sorti1(iel) = rtp(iel,irhiph)**gamagp
    enddo


!     CALCUL DE LA CHALEUR MASSIQUE A VOLUME CONSTANT
!     ===============================================

!     elle est constante           !

!                                                           P
!     CALCUL DE L'ENTROPIE EN FONCTION DE P ET rho : s = --------
!     ============================================          gamma
!                                                        rho
  elseif(iccfth.eq.6) then

!     Test des valeurs de rho (on peut éliminer ce test pour optimiser)
    ierr = 0
    do iel = 1, ncel
      if(rtp(iel,irhiph).le.0.d0) then
        write(nfecra,4030)rtp(iel,irhiph),iel
      endif
    enddo
    if(ierr.eq.1) then
      call csexit (1)
    endif

    do iel = 1, ncel
      sorti1(iel) = rtp(iel,ipriph) / (rtp(iel,irhiph)**gamagp)
    enddo


!     CALCUL DE epsilon - Cv.T : epsilon - Cv.T = 0
!     ========================

  elseif(iccfth.eq.7) then

! --- A l'interieur du domaine

    do iel = 1, ncel
      sorti1(iel) = 0.d0
    enddo

! --- Sur les bords

    do ifac = 1, nfabor
      sorti2(ifac) = 0.d0
    enddo


!     CALCUL DES CONDITIONS AUX LIMITES : sur la face IFAC = IFAC0
!     =================================

!     PAROI
!     -----

  elseif(iccfth.eq.91) then

    ifac = ifac0
    iel  = ifabor(ifac)

! --- Calcul du nombre de Mach normal a la frontiere

    xmach =                                                       &
         ( rtp(iel,iuiph)*surfbo(1,ifac)                          &
         + rtp(iel,iviph)*surfbo(2,ifac)                          &
         + rtp(iel,iwiph)*surfbo(3,ifac) ) / ra(isrfbn+ifac-1)    &
         / sqrt( gamagp * rtp(iel,ipriph) / rtp(iel,irhiph) )

! --- Pression

!     On utilise un Neumann en pression. Bien que cela ne permette
!       pas d'utiliser Rusanov, on espere que cela stabilise.
!     On ajoute un test pour eviter de basculer
!       de detente a choc d'un pas de temps à l'autre

!     Detente
    if(xmach.lt.0.d0.and.coefb(ifac,iclp).le.1.d0) then

      if(xmach.gt.2.d0/(1.d0-gamagp)) then
        coefb(ifac,iclp) = (1.d0 + (gamagp-1.d0)/2.d0 * xmach)    &
             ** (2.d0*gamagp/(gamagp-1.d0))
!         Si on detend trop fort : Dirichlet nul en pression
!           (la valeur de COEFB ici est utilisee comme indicateur
!            et sera retraitée plus tard dans cfxtcl)
      else
        coefb(ifac,iclp) = rinfin
      endif

!     Choc
    elseif(xmach.gt.0.d0.and.coefb(ifac,iclp).ge.1.d0) then

      coefb(ifac,iclp) = 1.d0 + gamagp*xmach                      &
            *( (gamagp+1.d0)/4.d0*xmach                           &
                + sqrt(1.d0 + (gamagp+1.d0)**2/16.d0*xmach**2) )

!     Oscillation ou Mach nul
    else
      coefb(ifac,iclp) = 1.d0
    endif


!     SYMETRIE
!     --------

  elseif(iccfth.eq.90) then

    ifac = ifac0
    iel  = ifabor(ifac)

! --- On impose un Neumann Homogene en Pression par defaut
!         Pas d'intervention requise


! --- ENTREE SUBSONIQUE, RHO ET U DONNES (subsonique par hypothèse)
!     ----------------------------------

!     Contrairement à la conception initiale, on impose un Dirichlet
!       explicite sur la pression et non pas un flux (on se base cpdt
!       sur la meme valeur physique).
!     Ceci a l'avantage de permettre l'utilisation de Rusanov pour
!       lisser les conditions d'entree.
!     En outre, ceci permet de ne pas avoir à gérer le remplissage de
!       coefb ici (c'est secondaire, car il faut le faire pour la paroi).
!     Si on a des oscillations en temps, on pourra essayer d'eviter
!       le basculement de détente à choc d'un pas de temps à l'autre.
!     La pertinence de cette démarche reste à démontrer.

  elseif(iccfth.eq.92) then

    ifac = ifac0
    iel  = ifabor(ifac)

! --- Calcul du nombre de Mach entrant normalement a la frontiere

    xmachi =                                                      &
         ( rtp(iel,iuiph)*surfbo(1,ifac)                          &
         + rtp(iel,iviph)*surfbo(2,ifac)                          &
         + rtp(iel,iwiph)*surfbo(3,ifac) ) / ra(isrfbn+ifac-1)    &
         / sqrt( gamagp * rtp(iel,ipriph) / rtp(iel,irhiph) )
    xmache =                                                      &
         ( coefa(ifac,iclu)*surfbo(1,ifac)                        &
         + coefa(ifac,iclv)*surfbo(2,ifac)                        &
         + coefa(ifac,iclw)*surfbo(3,ifac) ) /ra(isrfbn+ifac-1)   &
         / sqrt( gamagp * rtp(iel,ipriph) / rtp(iel,irhiph) )
    dxmach = xmachi - xmache

! --- Pression : Detente
    if(dxmach.le.0.d0) then

      if(dxmach.gt.2.d0/(1.d0-gamagp)) then
        coefa(ifac,iclp) = rtp(iel,ipriph)*                       &
             ( (1.d0 + (gamagp-1.d0)*0.50d0*dxmach)               &
               ** (2.d0*gamagp/(gamagp-1.d0))    )
      elseif(dxmach.le.2.d0/(1.d0-gamagp) ) then
        coefa(ifac,iclp) = 0.d0
      endif

! --- Pression : Choc
    else
      coefa(ifac,iclp) = rtp(iel,ipriph)*                         &
           (  1.d0 + gamagp*dxmach                                &
           *( (gamagp+1.d0)*0.25d0*dxmach                         &
           + sqrt(1.d0 + (gamagp+1.d0)**2/16.d0*dxmach**2) )  )
    endif

! --- Energie totale
    coefa(ifac,iclp) = rtp(iel,ipriph)
    coefa(ifac,icle) =                                            &
         coefa(ifac,iclp)/((gamagp-1.d0)*coefa(ifac,iclr))        &
         + 0.5d0*(coefa(ifac,iclu)**2                             &
                + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)


! --- ENTREE SUBSONIQUE, RHO*U ET RHO*U*H DONNES (subsonique par hypothèse)
!     ------------------------------------------

!     Ceci reste a completer :
  elseif(iccfth.eq.94) then

    ifac = ifac0
    iel  = ifabor(ifac)

    write(nfecra,7000)

    call csexit (1)
    !==========

!     P   : a calculer avec un newton
!     u et rho en fonction de P
!     E en fonction de H
!     (je l'ai ecrit, reste à le coder ...)


!     SORTIE SUBSONIQUE (subsonique par hypothèse)
!     -----------------

  elseif(iccfth.eq.93) then

    ifac = ifac0
    iel  = ifabor(ifac)

! --- Detente :
    if(coefa(ifac,iclp).le.rtp(iel,ipriph)) then

!       Masse volumique
      coefa(ifac,iclr) = rtp(iel,irhiph)                          &
           * (coefa(ifac,iclp)/rtp(iel,ipriph))**(1.d0/gamagp)

!       Vitesse
      coefa(ifac,iclu) = rtp(iel,iuiph)                           &
           + 2.d0/(gamagp-1.d0)                                   &
           * sqrt(gamagp*rtp(iel,ipriph)/rtp(iel,irhiph))         &
           * (1.d0-(coefa(ifac,iclp)/rtp(iel,ipriph)              &
                        )**((gamagp-1.d0)/(2.d0*gamagp)))         &
           * surfbo(1,ifac)/ra(isrfbn+ifac-1)

      coefa(ifac,iclv) = rtp(iel,iviph)                           &
           + 2.d0/(gamagp-1.d0)                                   &
           * sqrt( gamagp*rtp(iel,ipriph)/rtp(iel,irhiph))        &
           * (1.d0-(coefa(ifac,iclp)/rtp(iel,ipriph)              &
                        )**((gamagp-1.d0)/(2.d0*gamagp)))         &
           * surfbo(2,ifac)/ra(isrfbn+ifac-1)

      coefa(ifac,iclw) = rtp(iel,iwiph)                           &
           + 2.d0/(gamagp-1.d0)                                   &
           * sqrt( gamagp*rtp(iel,ipriph)/rtp(iel,irhiph))        &
           * (1.d0-(coefa(ifac,iclp)/rtp(iel,ipriph)              &
                        )**((gamagp-1.d0)/(2.d0/gamagp)))         &
           * surfbo(3,ifac)/ra(isrfbn+ifac-1)

!       Energie totale
      coefa(ifac,icle) =                                          &
           coefa(ifac,iclp)/((gamagp-1.d0)*coefa(ifac,iclr))      &
           + 0.5d0*(coefa(ifac,iclu)**2                           &
                  + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)

! --- Choc :
    else

!       Masse volumique
      coefa(ifac,iclr) = rtp(iel,irhiph)                          &
           * ( (gamagp+1.d0)*coefa(ifac,iclp)                     &
             + (gamagp-1.d0)*rtp(iel,ipriph) )                    &
           / ( (gamagp-1.d0)*coefa(ifac,iclp)                     &
             + (gamagp+1.d0)*rtp(iel,ipriph) )

!       Vitesse
      coefa(ifac,iclu) = rtp(iel,iuiph)                           &
           - (coefa(ifac,iclp)-rtp(iel,ipriph))                   &
           * sqrt(2.d0/                                           &
                  (rtp(iel,irhiph)                                &
                   *((gamagp+1.d0)*coefa(ifac,iclp)               &
                    +(gamagp-1.d0)*rtp(iel,ipriph) )))            &
           * surfbo(1,ifac)/ra(isrfbn+ifac-1)

      coefa(ifac,iclv) = rtp(iel,iviph)                           &
           - (coefa(ifac,iclp)-rtp(iel,ipriph))                   &
           * sqrt(2.d0/                                           &
                  (rtp(iel,irhiph)                                &
                   *((gamagp+1.d0)*coefa(ifac,iclp)               &
                    +(gamagp-1.d0)*rtp(iel,ipriph) )))            &
           * surfbo(2,ifac)/ra(isrfbn+ifac-1)

      coefa(ifac,iclw) = rtp(iel,iwiph)                           &
           - (coefa(ifac,iclp)-rtp(iel,ipriph))                   &
           * sqrt(2.d0/                                           &
                  (rtp(iel,irhiph)                                &
                   *((gamagp+1.d0)*coefa(ifac,iclp)               &
                    +(gamagp-1.d0)*rtp(iel,ipriph) )))            &
           * surfbo(3,ifac)/ra(isrfbn+ifac-1)

!       Energie totale
      coefa(ifac,icle) =                                          &
           coefa(ifac,iclp)/((gamagp-1.d0)*coefa(ifac,iclr))      &
           + 0.5d0*(coefa(ifac,iclu)**2                           &
                  + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)

    endif


!     T ET e EN FONCTION DE P ET rho
!     ------------------------------
!      on a donné des valeurs positives strictement (hypothèse)

  elseif(iccfth.eq.912.or.iccfth.eq.60900) then

    ifac = ifac0
    iel  = ifabor(ifac)

!     Temperature
    coefa(ifac,iclt) =                                            &
         xmasml*coefa(ifac,iclp)/(rr*coefa(ifac,iclr))

!     Energie totale
    coefa(ifac,icle) =                                            &
         cv0(iphas)*coefa(ifac,iclt)                              &
         + 0.5d0*( coefa(ifac,iclu)**2                            &
                 + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2 )


!     rho ET e EN FONCTION DE P ET T
!     ------------------------------

  elseif(iccfth.eq.913.or.iccfth.eq.100900) then

    ifac = ifac0
    iel  = ifabor(ifac)

!     Masse volumique
    coefa(ifac,iclr) =                                            &
         xmasml*coefa(ifac,iclp)/(rr*coefa(ifac,iclt))

!     Energie totale
    coefa(ifac,icle) =                                            &
         cv0(iphas)*coefa(ifac,iclt)                              &
         + 0.5d0*( coefa(ifac,iclu)**2                            &
                 + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2 )


!     rho ET T EN FONCTION DE P ET e
!     ------------------------------

  elseif(iccfth.eq.914.or.iccfth.eq.140900) then

    ifac = ifac0
    iel  = ifabor(ifac)

!     Masse volumique
    coefa(ifac,iclr) = coefa(ifac,iclp)/( (gamagp-1.d0)*          &
         (coefa(ifac,icle)                                        &
         - 0.5d0*( coefa(ifac,iclu)**2                            &
                 + coefa(ifac,iclv)**2                            &
                 + coefa(ifac,iclw)**2 ) ) )

!     Temperature
    coefa(ifac,iclt)=                                             &
         xmasml*coefa(ifac,iclp)/(rr*coefa(ifac,iclr))


!     P ET e EN FONCTION DE rho ET T
!     ------------------------------

  elseif(iccfth.eq.923.or.iccfth.eq.150900) then

    ifac = ifac0
    iel  = ifabor(ifac)

!     Pression
    coefa(ifac,iclp) = coefa(ifac,iclr)*rr/xmasml                 &
                                       *coefa(ifac,iclt)

!     Energie totale
    coefa(ifac,icle) = cv0(iphas) * coefa(ifac,iclt)              &
         + 0.5d0*( coefa(ifac,iclu)**2                            &
                 + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2 )


!     P ET T EN FONCTION DE rho ET e
!     ------------------------------

  elseif(iccfth.eq.924.or.iccfth.eq.210900) then

    ifac = ifac0
    iel  = ifabor(ifac)

!     Pression
    coefa(ifac,iclp) = (gamagp-1.d0)*coefa(ifac,iclr)             &
          *( coefa(ifac,icle)                                     &
            - 0.5d0*( coefa(ifac,iclu)**2                         &
                    + coefa(ifac,iclv)**2                         &
                    + coefa(ifac,iclw)**2 ) )


!     Temperature
    coefa(ifac,iclt)=                                             &
         xmasml*coefa(ifac,iclp)/(rr*coefa(ifac,iclr))


! --- Fin de test sur GAZ PARFAIT
  endif


!===============================================================================
! 2.  GAZ PARFAIT AVEC GAMMA VARIABLE
!===============================================================================

elseif(ieos(iphas).eq.2) then

!===============================================================================

!     PARAMETRES
!     ==========

!-------------------------------------------------------------------------------

!     EXEMPLES  (A recopier et a adapter apres 'PARAMETRES')
  if(0.eq.1) then
!-------------------------------------------------------------------------------

! --> EXEMPLE 1 : GAZ PARFAIT A 3 CONSTITUANTS

!     Phase
  if(iphas.eq.1) then

!     Masses molaires des constituants (en kg/mol)
    cstgr(1)  = 18.d-3
    cstgr(2)  = 32.d-3
    cstgr(3)  = 28.d-3

!     Calcul de la masse molaire du melange au centre des cellules
    if(iccfth.gt.0) then
      do iel = 1, ncel
        xmasmr(iel) = 1.d0 / ( rtp(iel,isca(1))/cstgr(1)          &
                             + rtp(iel,isca(2))/cstgr(2)          &
                             + rtp(iel,isca(3))/cstgr(3) )
      enddo


!     Calcul du GAMMA equivalent au centre des cellules
      do iel = 1, ncel
        gamagr(iel) = propce(iel,ipproc(icp(iphas)))              &
           / ( propce(iel,ipproc(icp(iphas))) - rr/xmasmr(iel) )
      enddo
    endif

  endif

!-------------------------------------------------------------------------------

!     FIN DES EXEMPLES
  endif
!-------------------------------------------------------------------------------

!===============================================================================


! TEST DE VALIDITE DE GAMAGR : GAMAGR >= 1
! ==========================

  ierr = 0

  do iel = 1, ncel
    if(iccfth.gt.0 .and. gamagr(iel).lt.1.d0) then
      ierr = 1
      write(nfecra,1020) iel, gamagr(iel)
    endif
  enddo

  if(ierr.eq.1) then
    call csexit (1)
  endif


! OPTIONS DE CALCUL  : Cp ET Cv VARIABLES
! =================

  if(iccfth.eq.-1) then

    icp(iphas) = 1
    cp0(iphas) = epzero
    icv(iphas) = 1
    cv0(iphas) = epzero


! INITIALISATIONS PAR DEFAUT :
! ==========================

  elseif(iccfth.eq.0) then

    do iel = 1, ncel
      propce(iel,ipproc(icp(iphas))) = cp0(iphas)
      propce(iel,ipproc(icv(iphas))) =                            &
           cp0(iphas) - rr/xmasmr(iel)
      rtp(iel,irhiph) = p0(iphas)*xmasmr(iel)/rr/t0(iphas)
      rtp(iel,ieniph) =                                           &
           propce(iel,ipproc(icv(iphas)))*t0(iphas)
    enddo


! CALCUL DE T ET e EN FONCTION DE P ET rho :
! ========================================

  elseif(iccfth.eq.12) then

    do iel = 1, ncel

!     Temperature
      sorti1(iel) =                                               &
           xmasmr(iel)/rr*rtp(iel,ipriph)/rtp(iel,irhiph)

!     Energie totale
      sorti2(iel) = propce(iel,ipproc(icv(iphas)))*sorti1(iel)    &
    + 0.5d0*( rtp(iel,iuiph)**2                                   &
           + rtp(iel,iviph)**2 + rtp(iel,iwiph)**2 )

    enddo


    if(imodif.gt.0) then
!       Affectation directe a RTP
      do iel = 1, ncel
        rtp(iel,itkiph) = sorti1(iel)
        rtp(iel,ieniph) = sorti2(iel)
      enddo
    endif


! CALCUL DE rho ET e EN FONCTION DE P ET T :
! ========================================

  elseif(iccfth.eq.13) then

    do iel = 1, ncel

!     Masse volumique
      sorti1(iel) =                                               &
           xmasmr(iel)/rr*rtp(iel,ipriph)/rtp(iel,itkiph)

!     Energie totale
      sorti2(iel) =                                               &
           propce(iel,ipproc(icv(iphas)))*rtp(iel,itkiph)         &
    + 0.5d0*( rtp(iel,iuiph)**2                                   &
           + rtp(iel,iviph)**2 + rtp(iel,iwiph)**2 )

    enddo


    if(imodif.gt.0) then
!       Affectation directe a RTP
      do iel = 1, ncel
        rtp(iel,irhiph) = sorti1(iel)
        rtp(iel,ieniph) = sorti2(iel)
      enddo
    endif


! CALCUL DE rho ET T EN FONCTION DE P ET e :
! ========================================

  elseif(iccfth.eq.14) then

    do iel = 1, ncel

!     Masse volumique
      sorti1(iel) =                                               &
           rtp(iel,ipriph)/(gamagr(iel)-1.d0)/( rtp(iel,ieniph)   &
  - 0.5d0*( rtp(iel,iuiph)**2                                     &
           + rtp(iel,iviph)**2 + rtp(iel,iwiph)**2 ) )

!     Temperature
      sorti2(iel) = xmasmr(iel)/rr*rtp(iel,ipriph)/sorti1(iel)

    enddo


    if(imodif.gt.0) then
!       Affectation directe a RTP
      do iel = 1, ncel
        rtp(iel,irhiph) = sorti1(iel)
        rtp(iel,itkiph) = sorti2(iel)
      enddo
    endif


! CALCUL DE P ET e EN FONCTION DE rho ET T :
! ========================================

  elseif(iccfth.eq.23) then

    do iel = 1, ncel

!     Pression
      sorti1(iel) =                                               &
           rtp(iel,irhiph)*rr/xmasmr(iel)*rtp(iel,itkiph)

!     Energie totale
      sorti2(iel) =                                               &
           propce(iel,ipproc(icv(iphas)))*rtp(iel,itkiph)         &
    + 0.5d0*( rtp(iel,iuiph)**2                                   &
           + rtp(iel,iviph)**2 + rtp(iel,iwiph)**2 )

    enddo


    if(imodif.gt.0) then
!       Affectation directe a RTP
      do iel = 1, ncel
        rtp(iel,ipriph) = sorti1(iel)
        rtp(iel,ieniph) = sorti2(iel)
      enddo
    endif


! CALCUL DE P ET T EN FONCTION DE rho ET e :
! ========================================

  elseif(iccfth.eq.24) then

    do iel = 1, ncel

!     Pression
      sorti1(iel) =                                               &
           (gamagr(iel)-1.d0)*rtp(iel,irhiph)*( rtp(iel,ieniph)   &
  - 0.5d0*( rtp(iel,iuiph)**2                                     &
           + rtp(iel,iviph)**2 + rtp(iel,iwiph)**2 ) )

!     Temperature
      sorti2(iel) = xmasmr(iel)/rr*sorti1(iel)/rtp(iel,irhiph)

    enddo

    if(imodif.gt.0) then
!       Affectation directe a RTP
      do iel = 1, ncel
        rtp(iel,ipriph) = sorti1(iel)
        rtp(iel,itkiph) = sorti2(iel)
      enddo
    endif

!            2                            2          P
! CALCUL DE c  EN FONCTION DE P ET rho : c  = gamma*---
! ====================================              rho

  elseif(iccfth.eq.126) then

    do iel = 1, ncel

! --- Test de positivite de P
      if(rtp(iel,ipriph).lt.0.d0) then
        write(nfecra,1110) iel , rtp(iel,ipriph)
        ierr = 1

! --- Test de positivite de rho
      elseif(rtp(iel,irhiph).le.0.d0) then
        write(nfecra,1120) iel , rtp(iel,irhiph)
        ierr = 1

      else

! --- Calcul
        sorti1(iel) =                                             &
             gamagr(iel) * rtp(iel,ipriph) / rtp(iel,irhiph)

      endif

    enddo

    if(ierr.eq.1) call csexit (1)

!                                                    gamma
! CALCUL DE beta EN FONCTION DE P ET rho : beta = rho
! ======================================

  elseif(iccfth.eq.162) then

    do iel = 1, ncel

! --- Test de positivite de rho
      if(rtp(iel,irhiph).lt.0.d0) then
        write(nfecra,1220) iel , rtp(iel,irhiph)
        ierr = 1

      else

! --- Calcul
        sorti1(iel) = rtp(iel,irhiph)**gamagr(iel)

      endif

    enddo

    if(ierr.eq.1) call csexit (1)

!                           R
! CALCUL DE Cv : Cv = Cp - ---
! ============              M

  elseif(iccfth.eq.432) then

    do iel = 1, ncel

      sorti1(iel) = propce(iel,ipproc(icp(iphas)))-rr/xmasmr(iel)

    enddo

    if(ierr.eq.1) call csexit (1)

!                                                       P
! CALCUL DE L'ENTROPIE EN FONCTION DE P ET rho : s = --------
! ============================================          gamma
!                                                    rho
  elseif(iccfth.eq.6) then

    do iel = 1, ncel

! --- Test de positivite de P
      if(rtp(iel,ipriph).lt.0.d0) then
        write(nfecra,1310) iel , rtp(iel,ipriph)
        ierr = 1

! --- Test de positivite de rho
      elseif(rtp(iel,irhiph).le.0.d0) then
        write(nfecra,1320) iel , rtp(iel,irhiph)
        ierr = 1

      else

! --- Calcul
        sorti1(iel) =                                             &
             rtp(iel,ipriph) / (rtp(iel,irhiph)**gamagr(iel))

      endif

    enddo

    if(ierr.eq.1) call csexit (1)


! CALCUL DE epsilon - Cv.T : epsilon - Cv.T = 0
! ========================

  elseif(iccfth.eq.7) then

! --- A l'interieur du domaine

    do iel = 1, ncel

      sorti1(iel) = 0.d0

    enddo

! --- Sur les bords

    do ifac = 1, nfabor

      sorti2(ifac) = 0.d0

    enddo

    if(ierr.eq.1) call csexit (1)


! CALCUL DES CONDITIONS AUX LIMITES :
! =================================

!     PAROI/SYMETRIE
!     --------------
  elseif(iccfth.eq.91) then

    ifac = ifac0
    iel  = ifabor(ifac)

! --- Calcul du nombre de Mach normal a la frontiere

    xmach = ( rtp(iel,iuiph)*surfbo(1,ifac)                       &
           + rtp(iel,iviph)*surfbo(2,ifac)                        &
           + rtp(iel,iwiph)*surfbo(3,ifac) ) / ra(isrfbn+ifac-1)  &
         / sqrt( gamagr(iel)*rtp(iel,ipriph)/rtp(iel,irhiph) )

! --- Pression / Entropie : Detente

    coefa(ifac,iclp) = 0.d0

    if(xmach.le.0.d0 .and. xmach.gt.2.d0/(1.d0-gamagr(iel))) then
      coefb(ifac,iclp) = (1.d0 + (gamagr(iel)-1.d0)/2.d0 * xmach) &
           ** (2.d0*gamagr(iel)/(gamagr(iel)-1.d0))
      coefb(ifac,iclt) = 1.d0

    elseif(xmach.le.2.d0/(1.d0-gamagr(iel)) ) then
      coefb(ifac,iclp) = 0.d0
      coefb(ifac,iclt) = 1.d0

! --- Pression / Entropie : Choc

    else
      coefb(ifac,iclp) = 1.d0 + gamagr(iel)*xmach                 &
            *( (gamagr(iel)+1.d0)/4.d0*xmach                      &
           + sqrt(1.d0 + (gamagr(iel)+1.d0)**2/16.d0*xmach**2) )
      coefb(ifac,iclt) = coefb(ifac,iclp)/(1.d0-coefb(ifac,iclp)) &
          / rtp(iel,ipriph) * ( rtp(iel,irhiph)                   &
              * (rtp(iel,iuiph)**2                                &
                +rtp(iel,iviph)**2+rtp(iel,iwiph)**2)             &
              + rtp(iel,ipriph) *(1.d0-coefb(ifac,iclp)) )
    endif

! --- Energie totale : Epsilon sup

    coefa(ifac,icle) = 0.d0

    if(ierr.eq.1) call csexit (1)

!     ENTREE
!     ------
  elseif(iccfth.eq.92) then

    ifac = ifac0
    iel  = ifabor(ifac)

! --- Calcul du nombre de Mach entrant normalement a la frontiere

    xmachi = ( rtp(iel,iuiph)*surfbo(1,ifac)                      &
         + rtp(iel,iviph)*surfbo(2,ifac)                          &
         + rtp(iel,iwiph)*surfbo(3,ifac) )/ra(isrfbn+ifac-1)      &
         / sqrt(gamagr(iel)*rtp(iel,ipriph)/rtp(iel,irhiph))
    xmache = ( coefa(ifac,iclu)*surfbo(1,ifac)                    &
         + coefa(ifac,iclv)*surfbo(2,ifac)                        &
         + coefa(ifac,iclw)*surfbo(3,ifac) )/ra(isrfbn+ifac-1)    &
         / sqrt(gamagr(iel)*rtp(iel,ipriph)/rtp(iel,irhiph))
    dxmach = xmachi - xmache

! --- Pression : Detente
    if(dxmach.le.0.d0) then

      if(dxmach.gt.2.d0/(1.d0-gamagr(iel))) then
        coefa(ifac,iclp) = rtp(iel,ipriph)*                       &
             ( (1.d0 + (gamagr(iel)-1.d0)*0.50d0*dxmach)          &
               ** (2.d0*gamagr(iel)/(gamagr(iel)-1.d0))  )
      elseif(dxmach.le.2.d0/(1.d0-gamagr(iel)) ) then
        coefa(ifac,iclp) = 0.d0
      endif

! --- Pression : Choc
    else
      coefa(ifac,iclp) = rtp(iel,ipriph)*                         &
           (  1.d0 + gamagr(iel)*dxmach                           &
           *( (gamagr(iel)+1.d0)*0.25d0*dxmach                    &
           + sqrt(1.d0 + (gamagr(iel)+1.d0)**2/16.d0              &
                                           *dxmach**2) )  )
    endif

! --- Energie totale

    coefa(ifac,iclp) = rtp(iel,ipriph)

    coefa(ifac,icle) =                                            &
         coefa(ifac,iclp)/((gamagr(iel)-1.d0)*coefa(ifac,iclr))   &
         + 0.5d0*(coefa(ifac,iclu)**2                             &
                + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)

!     SORTIE
!     ------
  elseif(iccfth.eq.93) then

    ifac = ifac0
    iel  = ifabor(ifac)

! --- Calcul du nombre de Mach sortant normalement a la frontiere

    xmach = ( rtp(iel,iuiph)*surfbo(1,ifac)                       &
           + rtp(iel,iviph)*surfbo(2,ifac)                        &
           + rtp(iel,iwiph)*surfbo(3,ifac) ) / ra(isrfbn+ifac-1)  &
         / sqrt(gamagr(iel)*rtp(iel,ipriph)/rtp(iel,irhiph))

! --- Sortie supersonique : Dirichlet sur toutes les variables
!     -------------------
    if(xmach.ge.1.d0) then
      do ivar = 1, nvar
        coefa(ifac,iclrtp(ivar,icoef)) = rtp(iel,ivar)
      enddo

!           Entropie
      coefa(ifac,iclt) =                                          &
           rtp(iel,ipriph)/rtp(iel,irhiph)**gamagr(iel)

! --- Sortie subsonique
!     -----------------
    elseif(xmach.lt.1.d0 .and. xmach.ge.0.d0) then

! --- Detente :
      if(coefa(ifac,iclp).le.rtp(iel,ipriph)) then

!       Masse volumique
        coefa(ifac,iclr) = rtp(iel,irhiph)                        &
             * (coefa(ifac,iclp)/rtp(iel,ipriph))                 &
                **(1.d0/gamagr(iel))

!       Vitesse
        coefa(ifac,iclu) = rtp(iel,iuiph)                         &
             + 2.d0/(gamagr(iel)-1.d0)                            &
 * sqrt( gamagr(iel) * rtp(iel,ipriph) / rtp(iel,irhiph) )        &
             * ( 1.d0                                             &
 - (coefa(ifac,iclp)/rtp(iel,ipriph))                             &
               **((gamagr(iel)-1.d0)/2.d0/gamagr(iel)) )          &
 * surfbo(1,ifac) / ra(isrfbn+ifac-1)

        coefa(ifac,iclv) = rtp(iel,iviph)                         &
             + 2.d0/(gamagr(iel)-1.d0)                            &
 * sqrt( gamagr(iel) * rtp(iel,ipriph) / rtp(iel,irhiph) )        &
             * ( 1.d0                                             &
 - (coefa(ifac,iclp)/rtp(iel,ipriph))                             &
               **((gamagr(iel)-1.d0)/2.d0/gamagr(iel)) )          &
 * surfbo(2,ifac) / ra(isrfbn+ifac-1)

        coefa(ifac,iclw) = rtp(iel,iwiph)                         &
             + 2.d0/(gamagr(iel)-1.d0)                            &
 * sqrt( gamagr(iel) * rtp(iel,ipriph) / rtp(iel,irhiph) )        &
             * ( 1.d0                                             &
 - (coefa(ifac,iclp)/rtp(iel,ipriph))                             &
               **((gamagr(iel)-1.d0)/2.d0/gamagr(iel)) )          &
 * surfbo(3,ifac) / ra(isrfbn+ifac-1)

!       Energie totale
        coefa(ifac,icle) = coefa(ifac,iclp)                       &
 /( (gamagr(iel)-1.d0)*coefa(ifac,iclr) )                         &
              + 0.5d0*(coefa(ifac,iclu)**2                        &
                     + coefa(ifac,iclv)**2                        &
                     + coefa(ifac,iclw)**2)

!       Entropie
        coefa(ifac,iclt) = coefa(ifac,iclp)                       &
                             /coefa(ifac,iclr)**gamagr(iel)

! --- Choc :
      else

!       Masse volumique
        coefa(ifac,iclr) = rtp(iel,irhiph)                        &
 * ( (gamagr(iel)+1.d0)*coefa(ifac,iclp)                          &
   + (gamagr(iel)-1.d0)*rtp(iel,ipriph) )                         &
 / ( (gamagr(iel)-1.d0)*coefa(ifac,iclp)                          &
   + (gamagr(iel)+1.d0)*rtp(iel,ipriph) )

!       Vitesse
        coefa(ifac,iclu) = rtp(iel,iuiph)                         &
 - (coefa(ifac,iclp)-rtp(iel,ipriph))*sqrt(2.d0/rtp(iel,irhiph)   &
 / ( (gamagr(iel)+1.d0)*coefa(ifac,iclp)                          &
   + (gamagr(iel)-1.d0)*rtp(iel,ipriph) ))                        &
 * surfbo(1,ifac) / ra(isrfbn+ifac-1)

        coefa(ifac,iclv) = rtp(iel,iviph)                         &
 - (coefa(ifac,iclp)-rtp(iel,ipriph))*sqrt(2.d0/rtp(iel,irhiph)   &
 / ( (gamagr(iel)+1.d0)*coefa(ifac,iclp)                          &
   + (gamagr(iel)-1.d0)*rtp(iel,ipriph) ))                        &
 * surfbo(2,ifac) / ra(isrfbn+ifac-1)

        coefa(ifac,iclw) = rtp(iel,iwiph)                         &
 - (coefa(ifac,iclp)-rtp(iel,ipriph))*sqrt(2.d0/rtp(iel,irhiph)   &
 / ( (gamagr(iel)+1.d0)*coefa(ifac,iclp)                          &
   + (gamagr(iel)-1.d0)*rtp(iel,ipriph) ))                        &
 * surfbo(3,ifac) / ra(isrfbn+ifac-1)

!       Energie totale
        coefa(ifac,icle) = coefa(ifac,iclp)                       &
 /( (gamagr(iel)-1.d0)*coefa(ifac,iclr) )                         &
     + 0.5d0*(coefa(ifac,iclu)**2                                 &
            + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)

!       Entropie
        coefa(ifac,iclt) = coefa(ifac,iclp)                       &
                             /coefa(ifac,iclr)**gamagr(iel)

      endif

    else
      WRITE(NFECRA,*) 'ICCFTH = ',ICCFTH,'  MACH = ',XMACH
      ierr = 1
    endif

    if(ierr.eq.1) call csexit (1)


!     T ET e EN FONCTION DE P ET rho
!     ------------------------------

  elseif(iccfth.eq.912.or.iccfth.eq.60900) then

    ifac = ifac0
    iel  = ifabor(ifac)

!     Temperature
    coefa(ifac,iclt) = xmasmr(iel)/rr*coefa(ifac,iclp)            &
                                        /coefa(ifac,iclr)

!     Energie totale
    coefa(ifac,icle) = propce(iel,ipproc(icv(iphas)))             &
               * coefa(ifac,iclt) + 0.5d0*( coefa(ifac,iclu)**2   &
                    + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)


!     rho ET e EN FONCTION DE P ET T
!     ------------------------------

  elseif(iccfth.eq.913.or.iccfth.eq.100900) then

    ifac = ifac0
    iel  = ifabor(ifac)

!     Masse volumique
    coefa(ifac,iclr) = xmasmr(iel)/rr*coefa(ifac,iclp)            &
                                       /coefa(ifac,iclt)

!     Energie totale
    coefa(ifac,icle) = propce(iel,ipproc(icv(iphas)))             &
               * coefa(ifac,iclt) + 0.5d0*( coefa(ifac,iclu)**2   &
                    + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)


!     rho ET T EN FONCTION DE P ET e
!     ------------------------------

  elseif(iccfth.eq.914.or.iccfth.eq.140900) then

    ifac = ifac0
    iel  = ifabor(ifac)

!     Masse volumique
    coefa(ifac,iclr) = coefa(ifac,iclp)/(gamagr(iel)-1.d0)        &
           / (coefa(ifac,icle) - 0.5d0*( coefa(ifac,iclu)**2      &
                 + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2 ))

!     Temperature
    coefa(ifac,iclt)= xmasmr(iel)/rr*coefa(ifac,iclp)             &
                                       /coefa(ifac,iclr)


!     P ET e EN FONCTION DE rho ET T
!     ------------------------------

  elseif(iccfth.eq.923.or.iccfth.eq.150900) then

    ifac = ifac0
    iel  = ifabor(ifac)

!     Pression
    coefa(ifac,iclp) = coefa(ifac,iclr)*rr/xmasmr(iel)            &
                                       *coefa(ifac,iclt)

!     Energie totale
    coefa(ifac,icle) = propce(iel,ipproc(icv(iphas)))             &
               * coefa(ifac,iclt) + 0.5d0*( coefa(ifac,iclu)**2   &
                    + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)


!     P ET T EN FONCTION DE rho ET e
!     ------------------------------

  elseif(iccfth.eq.924.or.iccfth.eq.210900) then

    ifac = ifac0
    iel  = ifabor(ifac)

!     Pression
    coefa(ifac,iclp) = (gamagr(iel)-1.d0)*coefa(ifac,iclr)        &
          *( coefa(ifac,icle) - 0.5d0*( coefa(ifac,iclu)**2       &
                + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2 ) )


!     Temperature
    coefa(ifac,iclt)= xmasmr(iel)/rr*coefa(ifac,iclp)             &
                                       /coefa(ifac,iclr)


! --- Fin de test sur GAZ PARFAIT A GAMMA VARIABLE
  endif

! --- Fin de test sur les thermos
endif


!--------
! FORMATS
!--------

 1010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    GAMMA DOIT ETRE UN REEL SUPERIEUR OU EGAL A 1           ',/,&
'@    IL A POUR VALEUR ',E12.4                                 ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    GAMMA DOIT ETRE UN REEL SUPERIEUR OU EGAL A 1           ',/,&
'@    IL EST NEGATIF OU NUL DANS LA CELLULE ',I10              ,/,&
'@    IL A POUR VALEUR ',E12.4                                 ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,&
'@    =========   CALCUL DE LA MASSE VOLUMIQUE                ',/,&
'@                                                            ',/,&
'@    LA TEMPERATURE DOIT ETRE UN REEL POSITIF STRICTEMENT    ',/,&
'@    ELLE VAUT ',E12.4   ,' DANS LA CELLULE ',I10             ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,&
'@    =========   CALCUL DE LA TEMPERATURE                    ',/,&
'@                                                            ',/,&
'@    LA MASSE VOLUMIQUE DOIT ETRE UN REEL POSITIF STRICTEMENT',/,&
'@    ELLE VAUT ',E12.4   ,' DANS LA CELLULE ',I10             ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 4010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,&
'@    =========   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''C2''',/,&
'@                                                            ',/,&
'@    LA MASSE VOLUMIQUE DOIT ETRE UN REEL POSITIF STRICTEMENT',/,&
'@    ELLE VAUT ',E12.4   ,' DANS LA CELLULE ',I10             ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,&
'@    =========   CALCUL DE LA VARIABLE THERMO        ''BETA''',/,&
'@                                                            ',/,&
'@    LA MASSE VOLUMIQUE DOIT ETRE UN REEL POSITIF            ',/,&
'@    ELLE VAUT ',E12.4   ,' DANS LA CELLULE ',I10             ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,&
'@    =========   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''S'' ',/,&
'@                                                            ',/,&
'@    LA MASSE VOLUMIQUE DOIT ETRE UN REEL POSITIF STRICTEMENT',/,&
'@    ELLE VAUT ',E12.4   ,' DANS LA CELLULE ',I10             ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 7000 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : CONDITION A LA LIMITE NON DISPONIBLE USCFTH ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    La condition a la limite de type debit et debit         ',/,&
'@      enthalpique imposes n''est pas disponible dans la     ',/,&
'@      version courante.                                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Modifier uscfcl.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8000 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MASSE VOLUMIQUE NEGATIVE OU NULLE (USCFTH)  ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    rho est negatif ou nul en ',I10   ,' cellules           ',/,&
'@    On le limite a ',E12.4                                   ,/,&
'@    et on met fin au calcul.                                ',/,&
'@                                                            ',/,&
'@    Si on souhaite que le calcul se poursuive, on peut      ',/,&
'@      forcer un clipping standard par valeur inferieure en  ',/,&
'@      renseignant la variable SCAMIN associee au scalaire   ',/,&
'@      representant la masse volumique : ISCA(IRHO(IPHAS))   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8100 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ENERGIE INTERNE NEGATIVE OU NULLE (USCFTH)  ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    e-0.5*u*u est negatif ou nul en ',I10   ,' cellules     ',/,&
'@    On le limite a ',E12.4                                   ,/,&
'@    et on met fin au calcul.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)



! formats a eliminer apres mise au point de gama variable



 1110 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USTHMO                        ',/,&
'@    =========   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''C2''',/,&
'@                                                            ',/,&
'@    P DOIT ETRE UN REEL POSITIF                             ',/,&
'@    IL EST NEGATIF DANS LA CELLULE ',I10                     ,/,&
'@    IL A POUR VALEUR ',E12.4                                 ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1120 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USTHMO                        ',/,&
'@    =========   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''C2''',/,&
'@                                                            ',/,&
'@    RHO DOIT ETRE UN REEL POSITIF STRICTEMENT               ',/,&
'@    IL EST NEGATIF OU NUL DANS LA CELLULE ',I10              ,/,&
'@    IL A POUR VALEUR ',E12.4                                 ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1220 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USTHMO                        ',/,&
'@    =========   CALCUL DE LA VARIABLE THERMO ''BETA''       ',/,&
'@                                                            ',/,&
'@    RHO DOIT ETRE UN REEL POSITIF                           ',/,&
'@    IL EST NEGATIF DANS LA CELLULE ',I10                     ,/,&
'@    IL A POUR VALEUR ',E12.4                                 ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1310 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USTHMO                        ',/,&
'@    =========   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''S'' ',/,&
'@                                                            ',/,&
'@    P DOIT ETRE UN REEL POSITIF                             ',/,&
'@    IL EST NEGATIF DANS LA CELLULE ',I10                     ,/,&
'@    IL A POUR VALEUR ',E12.4                                 ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1320 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS USTHMO                        ',/,&
'@    =========   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''S'' ',/,&
'@                                                            ',/,&
'@    RHO DOIT ETRE UN REEL POSITIF STRICTEMENT               ',/,&
'@    IL EST NEGATIF OU NUL DANS LA CELLULE ',I10              ,/,&
'@    IL A POUR VALEUR ',E12.4                                 ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!----
! FIN
!----

return

end subroutine
