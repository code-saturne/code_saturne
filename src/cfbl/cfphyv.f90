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

subroutine cfphyv &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ibrom  , izfppp ,          &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : COMPRESSIBLE SANS CHOC

! Calcul des proprietes physiques variables


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
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! nphmx            ! e  ! <-- ! nphsmx                                         !
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
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!   (nphmx   )     !    !     !                                                !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...3(ncelet    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
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
include "entsor.h"
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
integer          nideve , nrdeve , nituse , nrtuse , nphmx

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr), ibrom(nphmx)
integer          izfppp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra, ifinia
integer          iphas , iel
integer          ifac
integer          iirom , iiromb
integer          maxelt, ils

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0


!===============================================================================
! 2. ON DONNE LA MAIN A L'UTILISATEUR
!===============================================================================

maxelt = max(ncelet, nfac, nfabor)
ils    = idebia
ifinia = ils + maxelt
CALL IASIZE('CFPHYV',IFINIA)

iuscfp = 1
call uscfpv                                                       &
!==========
 ( ifinia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser , ra     )

!     Si IUSCFP = 0, l'utilisateur n'a pas inclus le ss pgm uscfpv dans
!       ses sources. C'est une erreur si Cp, Cv ou Lambda est variable.
!     On se contente de faire le test au premier passage.
if(ipass.eq.0) then
  ipass = ipass + 1
  do iphas = 1, nphas
    if((ivisls(itempk(iphas)).gt.0.or.                            &
        icp(iphas).gt.0.or.icv(iphas).gt.0).and.iuscfp.eq.0) then
      write(nfecra,1000)                                          &
           ivisls(itempk(iphas)),icp(iphas),icv(iphas)
      call csexit (1)
      !==========
    endif
  enddo
endif

!===============================================================================
! 3. MISE A JOUR DE LAMBDA/CV
!===============================================================================

! On a vérifié auparavant que CV0 était non nul.
! Si CV variable est nul, c'est une erreur utilisateur. On fait
!     un test à tous les passages (pas optimal), sachant que pour
!     le moment, on est en gaz parfait avec CV constant : si quelqu'un
!     essaye du CV variable, ce serait dommage que cela lui explose à la
!     figure pour de mauvaises raisons.
! Si IVISLS(IENERG(IPHAS)).EQ.0, on a forcement IVISLS(ITEMPK(IPHAS)).EQ.0
!     et ICV(IPHAS).EQ.0, par construction de IVISLS(IENERG(IPHAS)) dans
!     le sous-programme cfvarp

do iphas = 1, nphas

  if(ivisls(ienerg(iphas)).gt.0) then

    if(ivisls(itempk(iphas)).gt.0) then

      do iel = 1, ncel
        propce(iel,ipproc(ivisls(ienerg(iphas)))) =               &
             propce(iel,ipproc(ivisls(itempk(iphas))))
      enddo

    else
      do iel = 1, ncel
        propce(iel,ipproc(ivisls(ienerg(iphas)))) =               &
             visls0(itempk(iphas))
      enddo

    endif

    if(icv(iphas).gt.0) then

      do iel = 1, ncel
        if(propce(iel,ipproc(icv(iphas))).le.0.d0) then
          write(nfecra,2000)iel,propce(iel,ipproc(icv(iphas)))
          call csexit (1)
          !==========
        endif
      enddo

      do iel = 1, ncel
        propce(iel,ipproc(ivisls(ienerg(iphas)))) =               &
             propce(iel,ipproc(ivisls(ienerg(iphas))))            &
             / propce(iel,ipproc(icv(iphas)))
      enddo

    else

      do iel = 1, ncel
        propce(iel,ipproc(ivisls(ienerg(iphas)))) =               &
             propce(iel,ipproc(ivisls(ienerg(iphas))))            &
             / cv0(iphas)
      enddo

    endif

  else

    visls0(ienerg(iphas)) = visls0(itempk(iphas))/cv0(iphas)

  endif


enddo


!===============================================================================
! 3. MISE A JOUR DE ROM et ROMB :
!     On ne s'en sert a priori pas, mais la variable existe
!     On a ici des valeurs issues du pas de temps précédent (y compris
!       pour les conditions aux limites) ou issues de valeurs initiales
!     L'échange pério/parall sera fait dans phyvar.
!===============================================================================

do iphas = 1, nphas

  iirom  = ipproc(irom  (iphas))
  iiromb = ipprob(irom  (iphas))

  do iel = 1, ncel
    propce(iel,iirom)  = rtpa(iel,isca(irho(iphas)))
  enddo

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    propfb(ifac,iiromb) =                                         &
         coefa(ifac,iclrtp(isca(irho(iphas)),icoef))              &
         + coefb(ifac,iclrtp(isca(irho(iphas)),icoef))            &
         * rtpa(iel,isca(irho(iphas)))
  enddo

enddo


!--------
! FORMATS
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION (MODULE COMPRESSIBLE)  ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  Une ou plusieurs des propriétés suivantes a été déclarée  ',/,&
'@    variable (repérée ci-dessous par un indicateur non nul) ',/,&
'@    et une loi doit être fournie dans uscfpv.               ',/,&
'@         propriété                               indicateur ',/,&
'@     - conductivité thermique                    ',I10       ,/,&
'@     - capacité calorifique à pression constante ',I10       ,/,&
'@     - capacité calorifique à volume constant    ',I10       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Renseigner uscfpv ou déclarer les propriétés constantes et',/,&
'@    uniformes (uscfx2 pour la conductivité thermique,       ',/,&
'@    uscfth pour les capacités calorifiques).                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION (MODULE COMPRESSIBLE)  ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  La capacité calorifique à volume constant présente (au    ',/,&
'@    moins) une valeur négative ou nulle :                   ',/,&
'@    cellule ',I10,   '  Cv = ',E18.9                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfpv.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!----
! FIN
!----

return
end subroutine
