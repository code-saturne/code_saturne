!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine cfmsvl &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iscal  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , w10    , w11    , w12    ,          &
   wflmas , wflmab ,                                              &
   coefu  ,                                                       &
   rdevel , rtuser ,                                              &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR LA MASSE VOLUMIQUE SUR UN PAS DE TEMPS
!   (ALGORITHME COMPRESSIBLE EN RHO, U, E)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iscal            ! i  ! <-- ! scalar number                                  !
! itspdv           ! e  ! <-- ! calcul termes sources prod et dissip           !
!                  !    !     !  (0 : non , 1 : oui)                           !
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
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
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
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! w1..12(ncelet    ! tr ! --- ! tableau de travail                             !
! wflmas(nfac)     ! tr ! --- ! tableau de w flux de masse aux faces           !
! wflmab(nfabor    ! tr ! --- ! tableau de w flux de masse aux bords           !
! coefu(nfabo,3    ! tr ! --- ! tableau de travail cl de la qdm                !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "numvar.f90"
include "entsor.f90"
include "optcal.f90"
include "cstphy.f90"
include "cstnum.f90"
include "pointe.f90"
include "period.f90"
include "parall.f90"
include "ppppar.f90"
include "ppthch.f90"
include "ppincl.f90"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse
integer          iscal

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbrs(ncelet)
double precision rovsdt(ncelet)
double precision w1(ncelet) , w2(ncelet) , w3(ncelet)
double precision w4(ncelet) , w5(ncelet) , w6(ncelet)
double precision w7(ncelet) , w8(ncelet) , w9(ncelet)
double precision w10(ncelet), w11(ncelet), w12(ncelet)
double precision wflmas(nfac), wflmab(nfabor)
double precision coefu(nfabor,3)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          ivar  , iphas
integer          ifac  , iel
integer          init  , inc   , iccocg, isqrt , ii, jj, iii
integer          iclvar, iclvaf
integer          iflmas, iflmab
integer          ippvar, ipp   , iphydp
integer          nswrgp, imligp, iwarnp
integer          istatp, iconvp, idiffp, ireslp, ndircp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
double precision epsrgp, climgp, extrap, blencp, epsilp
double precision epsrsp
double precision sclnor

integer          iccfth, imodif
integer          iij
integer          iwfabg, iwfbbg
double precision dijpfx, dijpfy, dijpfz, pond
double precision diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz
double precision diipbx, diipby, diipbz
double precision pip   , pjp   , thetv, relaxp

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

iwfabg = idebra
iwfbbg = iwfabg+nfac
idebra = iwfbbg+nfabor
CALL RASIZE('CFMSVL',IDEBRA)

! --- Numero de phase associee au scalaire traite
iphas  = iphsca(iscal)

! --- Numero de variable de calcul et de post associe au scalaire traite
ivar   = isca(iscal)
ippvar = ipprtp(ivar)

! --- Numero des conditions aux limites
iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

! --- Flux de masse associe a l'energie
iflmas = ipprof(ifluma(isca(ienerg(iphas))))
iflmab = ipprob(ifluma(isca(ienerg(iphas))))

chaine = nomvar(ippvar)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!===============================================================================
! 2. TERMES SOURCES
!===============================================================================

! --> Initialisation

do iel = 1, ncel
  smbrs(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo


!     TERME SOURCE DE MASSE
!     =====================

if (ncesmp.gt.0) then
  do ii = 1, ncesmp
    iel = icetsm(ii)
    smbrs(iel) = smbrs(iel) + smacel(iel,ipr(iphas))*volume(iel)
  enddo
endif


!     TERME INSTATIONNAIRE
!     ====================

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel) + istat(ivar)*(volume(iel)/dt(iel))
enddo

!===============================================================================
! 3. CALCUL DU "FLUX DE MASSE" ET DE LA "VISCOSITE" AUX FACES
!===============================================================================

!     Ici VISCF et VISCB sont deux tableaux de travail.
!     On calcule WFLMAS et WFLMAB, WFABGS , WFBBGS

call cfmsgs                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iscal  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   wflmas , wflmab , ra(iwfabg) , ra(iwfbbg) ,                    &
!        ------   ------   ------   ------
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , w10    , w11    , w12    ,          &
   viscf  , viscb  , coefu  , xam    ,                            &
   rdevel , rtuser ,                                              &
   ra     )


!     Calcul du gradient de rho pour la reconstruction de rho
!       (contribution du terme convectif)

ircflp = ircflu(ivar)

if(ircflp.gt.0) then

  inc    = 1
  iccocg = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iphydp = 0
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtpa(1,ivar)    ,                                              &
   coefa(1,iclrtp(ivar,icoef)) , coefb(1,iclrtp(ivar,icoef)) ,    &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

  else
    do ii = 1, ncelet
      w1(ii) = 0.d0
      w2(ii) = 0.d0
      w3(ii) = 0.d0
    enddo
  endif


!     Au bord, on écrase WFLMAB pour imposer le débit souhaité.
!     Si on explicite le terme de convection (choix retenu par défaut,
!       seul testé), pas de problème.
!     Si on implicite le terme de convection, on n'impose pas
!       nécessairement le bon débit (cela dépend du traitement de
!       la convection au bord dans bilsc2 et de la condition à la
!       limite sur rho : ainsi, si l'on se place sur une sortie pour
!       laquelle la condition n'est pas valeur bord = valeur interne,
!       la convection étant traitée en upwind, c'est la valeur interne
!       de rho qui intervient dans le flux de bord et non pas
!       la valeur de bord comme supposé ci-dessous.

if(iconv(ivar).le.0) then

  if(ircflp.gt.0) then

    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      iij = idijpf-1+3*(ifac-1)
      dijpfx = ra(iij+1)
      dijpfy = ra(iij+2)
      dijpfz = ra(iij+3)

      pond   = ra(ipond-1+ifac)

      diipfx = cdgfac(1,ifac) - (xyzcen(1,ii)+                    &
               (1.d0-pond) * dijpfx)
      diipfy = cdgfac(2,ifac) - (xyzcen(2,ii)+                    &
               (1.d0-pond) * dijpfy)
      diipfz = cdgfac(3,ifac) - (xyzcen(3,ii)+                    &
               (1.d0-pond) * dijpfz)
      djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj)+                    &
                   pond  * dijpfx
      djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj)+                    &
                   pond  * dijpfy
      djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj)+                    &
                   pond  * dijpfz

      pip = rtpa(ii,ivar)                                         &
           + ircflp*(w1(ii)*diipfx+w2(ii)*diipfy+w3(ii)*diipfz)
      pjp = rtpa(jj,ivar)                                         &
           + ircflp*(w1(jj)*djjpfx+w2(jj)*djjpfy+w3(jj)*djjpfz)

      wflmas(ifac) = -0.5d0*                                      &
           ( pip          *(wflmas(ifac)+abs(wflmas(ifac)))       &
           + pjp          *(wflmas(ifac)-abs(wflmas(ifac))))
    enddo

  else
    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      wflmas(ifac) = -0.5d0*                                      &
           ( rtpa(ii,ivar)*(wflmas(ifac)+abs(wflmas(ifac)))       &
           + rtpa(jj,ivar)*(wflmas(ifac)-abs(wflmas(ifac))))
    enddo
  endif

  do ifac = 1, nfabor
    wflmab(ifac) = -propfb(ifac,iflmab)
  enddo

  init = 0
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
              ifacel,ifabor,wflmas,wflmab,smbrs)

  do ifac = 1, nfac
    ra(iwfabg+ifac-1) = - ra(iwfabg+ifac-1)
  enddo
  do ifac = 1, nfabor
    ra(iwfbbg+ifac-1) = - ra(iwfbbg+ifac-1)
  enddo
  init = 0
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
              ifacel,ifabor,ra(iwfabg),ra(iwfbbg),smbrs)

else

  if(ircflp.gt.0) then

    do ifac = 1, nfabor
      ii = ifabor(ifac)

      iii = idiipb-1+3*(ifac-1)
      diipbx = ra(iii+1)
      diipby = ra(iii+2)
      diipbz = ra(iii+3)

      pip = rtpa(ii,ivar)                                         &
           +ircflp*(w1(ii)*diipbx+w2(ii)*diipby+w3(ii)*diipbz)

      wflmab(ifac) = -propfb(ifac,iflmab)/                        &
           ( coefa(ifac,iclrtp(ivar,icoef))                       &
           + coefb(ifac,iclrtp(ivar,icoef))*pip           )
    enddo

  else
    do ifac = 1, nfabor
      ii = ifabor(ifac)
      wflmab(ifac) = -propfb(ifac,iflmab)/                        &
           ( coefa(ifac,iclrtp(ivar,icoef))                       &
           + coefb(ifac,iclrtp(ivar,icoef))*rtpa(ii,ivar) )
    enddo
  endif

endif


!     On calcule VISCF et VISCB

call cfmsvs                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , iscal  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   viscf  , viscb  ,                                              &
!        ------   ------
   w1     , w9     , w10    ,                                     &
   rdevel , rtuser ,                                              &
   ra     )


!     On annule la viscosité au bord afin que la contribution
!       au flux de bord soit exactement WFLMAB (dans lequel on a mis
!       le flux de masse souhaité). Si on n'annule pas, on risque
!       d'obtenir des contributions non nulles de la partie diffusive,
!       sauf si on a imposé une condition de Neumann homogene sur rho.
!       Pour le moment, on prefere prendre des precautions (la
!       modification de VISCB se traduit simplement par la modification
!       de la matrice portant sur les incréments, ou encore par une
!       une condition de Neumann homogène sur les incréments).

do ifac = 1, nfabor
  viscb (ifac) = 0.d0
enddo

!===============================================================================
! 4. RESOLUTION
!===============================================================================

istatp = istat (ivar)
iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
ndircp = ndircl(ivar)
nitmap = nitmax(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
iescap = 0
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
ipp    = ippvar
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)
thetv  = thetav(ivar)

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     wflmas          , wflmab          ,          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbrs  , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! 5. IMPRESSIONS ET CLIPPINGS
!===============================================================================

! --- Clipping aux bornes définies par l'utilisateur ou par defaut
!       (par défaut, pas de borne contraignate)

!     Valeur bidon
iii = 1

call clpsca                                                       &
!==========
 ( ncelet , ncel   , nvar   , nscal  , iscal  ,                   &
   propce , rtp(1,iii)      , rtp    )

! --- Traitement utilisateur pour gestion plus fine des bornes
!       et actions correctives éventuelles.
  iccfth = -2
  imodif = 0
  call uscfth                                                     &
  !==========
 ( idebia , idebra ,                                              &
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
   w7     , w8     , w9     , w10    ,                            &
   rdevel , rtuser , ra     )


! --- Bilan explicite (voir codits : on enleve l'increment)

if (iwarni(ivar).ge.2) then
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                       &
            - istat(ivar)*(volume(iel)/dt(iel))                   &
                *(rtp(iel,ivar)-rtpa(iel,ivar))                   &
                * max(0,min(nswrsm(ivar)-2,1))
  enddo
  isqrt = 1
  call prodsc(ncelet,ncel,isqrt,smbrs,smbrs,sclnor)
  write(nfecra,1200)chaine(1:8) ,sclnor
endif

!===============================================================================
! 6. COMMUNICATION DE RHO
!===============================================================================

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtp(1,ivar))
  !==========
endif

!     On ne remplit pas PROPCE et PROPFB ici, car on veut disposer de
!       rho à l'instant précédent pour résoudre en qdm et energie
!     On modifiera PROPCE et PROPFB après resolution de l'energie.


!===============================================================================
! 7. CALCUL DU FLUX DE MASSE ACOUSTIQUE AUX FACES
!===============================================================================

! Ce flux est stocke en tant que flux de masse associe a l'energie
inc    = 1
iccocg = 1

call cfbsc3                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtp(1,ivar)     , coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     wflmas          , wflmab          ,          &
   viscf  , viscb  ,                                              &
   propfa(1,iflmas), propfb(1,iflmab),                            &
!        ----------------  ----------------
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )


!     Si ICONV = 0, le terme convectif n'est pas traite par CFBSC3
!       il faut donc le rajouter a la main
!     Noter egalement que si ICONV = 0, PROPFB contient zero au sortir
!       de cfbsc3 et qu'on lui ajoute donc le flux de masse WFLMAB
!       (cohérent avec le flux imposé au bord et avec un signe négatif,
!        car WFLMAB etait, ci-dessus, utilise au second membre)
if(iconv(ivar).le.0) then
  do ifac = 1, nfac
    propfa(ifac,iflmas) = propfa(ifac,iflmas) - wflmas(ifac)      &
                                              - ra(iwfabg+ifac-1)
  enddo
  do ifac = 1, nfabor
    propfb(ifac,iflmab) = propfb(ifac,iflmab) - wflmab(ifac)      &
                                              - ra(iwfbbg+ifac-1)
  enddo
endif


!===============================================================================
! 8. ACTUALISATION DE LA PRESSION
!===============================================================================
!                               Pred      n+1  n
! On utilise l'equation d'etat P    =P(rho   ,e )

! --- Calcul de P au centre des cellules et actualisation de RTP
if(igrdpp(iphas).gt.0) then

  iccfth = 24
  imodif = 0
  call uscfth                                                     &
  !==========
 ( idebia , idebra ,                                              &
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
   rtp(1,ipr(iphas))        , w8     , w9     , w10    ,          &
!        -----------------
   rdevel , rtuser , ra     )

!===============================================================================
! 9. COMMUNICATION DE LA PRESSION
!===============================================================================

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ipr(iphas)))
    !==========
  endif

endif

!     Pas de vérification ni de clipping de la pression, car on
!       considere que la masse volumique et l'énergie ont été vérifiees,
!       sont correctes et donc que la pression l'est aussi (elle
!       vient d'être calculée avec un sous-pgm utilisateur).

!--------
! FORMATS
!--------

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A8                        ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)

!----
! FIN
!----

return
end subroutine
