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

subroutine raycli &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , isvhb  , isvtb  ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , itrifb , itypfb ,                                     &
   izfrad , isothm ,                                              &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  , hbord  , tbord  ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   text   , tint   , tempk  ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  1) Calcul des temperatures de paroi
!  2) Mise a jours des conditions aux limites de la variable
!     energetique

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
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! isvhb            ! e  ! <-- ! indicateur de sauvegarde des                   !
!                  !    !     !  coefficients d'echange aux bords              !
! isvtb            ! e  ! <-- ! indicateur de sauvegarde des                   !
!                  !    !     !  temperatures aux bords                        !
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
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! itrifb           ! ia ! <-- ! indirection for boundary faces ordering        !
!  (nfabor, nphas) !    !     !                                                !
! itypfb           ! ia ! --> ! boundary face types                            !
!  (nfabor, nphas) !    !     !                                                !
! izfrad(nfabor    ! te ! <-- ! numero de zone des faces de bord               !
! isothm(nfabor    ! te ! <-- ! type de condition de paroi                     !
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
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! hbord            ! tr ! --> ! coefficients d'echange aux bords               !
! (nfabor)         !    !     !                                                !
! tbord            ! tr ! --> ! temperature aux bords           i              !
! (nfabor)         !    !     !                                                !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! text (nfabor     ! tr ! --> ! temperature de bord externe                    !
! tint (nfabor     ! tr ! --> ! temperature de bord interne                    !
! tempk(ncelet)    ! tr ! --> ! temperature en kelvin                          !
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
include "pointe.h"
include "entsor.h"
include "parall.h"
include "ihmpre.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "radiat.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          isvhb  , isvtb

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          idevel(nideve), ituser(nituse), ia(*)
integer          izfrad(nfabor),isothm(nfabor)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision hbord(nfabor),tbord(nfabor)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)

double precision tempk(ncelet)
double precision text(nfabor), tint(nfabor)

double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, iel, ideb, ivart, iscat
integer          mode, iok, ifvu, ii, izonem, izone
integer          maxelt, idbia1, ils

double precision tmin , tmax   , tx
double precision cpp, xmtk

integer    ipacli
data       ipacli /0/
save       ipacli

!===============================================================================
!===============================================================================
! 0 - GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

maxelt = max(ncelet,nfac,nfabor)

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

!---> NUMERO DE PASSAGE RELATIF

ipacli = ipacli + 1
ideb = 0

!---> VALEURS MIN ET MAX ADMISSIBLES POUR LA TEMPERATURE DE PAROI
!         EN KELVIN

tmin = 0.d0
tmax = grand + tkelvi

!---> COEFF DE RELAX

!      TX est strictement superieur a 0 et inferieur ou egal a 1

!      Pour calculer la temperature de paroi, on calcule un increment
!      de temperature DeltaT entre l'etape courante n et l'etape
!      precedente n-1, puis on calcule :
!           n    n-1                                 n-1
!          T  = T    + DeltaT si le rapport DeltaT/T    =< TX, sinon

!           n    n-1                      n-1             n-1
!          T  = T    * (1 + TX *((DeltaT/T   ) / |DeltaT/T   |))

tx = 0.1d0


!---> INITIALISATIONS PAR DEFAUT BIDON

do ifac = 1,nfabor
  izfrad(ifac) = -1
  isothm(ifac) = -1
  propfb(ifac,ipprob(ixlam)) = -grand
  propfb(ifac,ipprob(iepa))  = -grand
  propfb(ifac,ipprob(ieps))  = -grand
  text  (ifac) = -grand
  tint  (ifac) = -grand
enddo


!===============================================================================
! 2. SI PAS DE FICHIER SUITE ALORS INITIALISATION AU PREMIER PASSAGE
!    DE TPAROI ET QINCID :
!      LECTURE DE L'INITIALISATION DE TPAROI A TINT
!      QINCID EST INITIALISE A STEPHN*TINT**4 (SI ON INITIALISE QINCID
!      A ZERO, ON AURA UN DEFICIT SUR LA CONDITION LIMITE DE LUMINANCE
!      AUX PAROIS AU 1er PAS DE TEMPS EN DOM)
!===============================================================================

if (ipacli.eq.1 .and. isuird.eq.0) then

! Indicateur : si non suite et premier pas de temps.
    ideb = 1

      do iel = 1,ncelet
        propce(iel,ipproc(itsri(1))) = zero
        propce(iel,ipproc(itsre(1))) = zero
      enddo

      do ifac = 1,nfabor
        propfb(ifac,ipprob(ihconv)) = zero
        propfb(ifac,ipprob(ifconv)) = zero
      enddo

!     On utilise TBORD comme auxiliaire pour l'appel a USRAY2
!       pour être sur que TPAROI ne sera pas modifié
!       (puisqu'on a TBORD libre)
!     On utilise FLUNET comme auxiliaire pour l'appel a USRAY2
!       pour être sur que QINCID ne sera pas modifié
!       (puisqu'on a FLUNET libre)

      do ifac = 1,nfabor
        tbord(ifac)      = zero
        propfb(ifac,ipprob(ifnet)) = zero
      enddo

!         - Interface Code_Saturne
!           ======================

      if (iihmpr.eq.1) then

!---> NUMERO DU SCALAIRE ET DE LA VARIABLE THERMIQUE
        ivart = isca(iscalt(irapha))

        call uiray2                                               &
        !==========
       ( itypfb, iparoi, iparug, ivart , izfrad,                  &
         isothm, itpimp, ipgrno, iprefl, ifgrno, ifrefl,          &
         nozppm, nfabor, nvar,                                    &
         propfb(1,ipprob(ieps)), propfb(1,ipprob(iepa)),          &
         tint, text,                                              &
         propfb(1,ipprob(ixlam)), rcodcl)

      endif

      ils    = idebia
      idbia1 = ils + maxelt
      CALL IASIZE('RAYCLI',IDBIA1)

      call usray2                                                 &
      !==========
 ( idbia1 , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , irapha  ,                                    &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   maxelt , ia(ils),                                              &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , izfrad , isothm ,                                     &
   tmin   , tmax   , tx     ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &
   tbord  , propfb(1,ipprob(ifnet))  , propfb(1,ipprob(ihconv))  ,&
   propfb(1,ipprob(ifconv)),                                      &
   propfb(1,ipprob(ixlam)) , propfb(1,ipprob(iepa)) ,             &
   propfb(1,ipprob(ieps))  ,                                      &
   text   , tint           ,                                      &
   ra     )

      write(nfecra,1000)

! Tparoi en Kelvin et QINCID en W/m2
      do ifac = 1,nfabor
        propfb(ifac,ipprob(itparo)) = tint(ifac)
        propfb(ifac,ipprob(iqinci)) = stephn*tint(ifac)**4
        if ( itypfb(ifac,irapha).eq.iparoi .or.                   &
             itypfb(ifac,irapha).eq.iparug ) then
          propfb(ifac,ipprob(itparo)) = tint(ifac)
          propfb(ifac,ipprob(iqinci)) = stephn*tint(ifac)**4
        else
          propfb(ifac,ipprob(itparo)) = 0.d0
          propfb(ifac,ipprob(iqinci)) = 0.d0
        endif
      enddo

! Fin détection premier passage
endif

!===============================================================================
! 3. PHASE
!===============================================================================

!---> NUMERO DU SCALAIRE ET DE LA VARIABLE THERMIQUE
  iscat = iscalt(irapha)
  ivart = isca(iscalt(irapha))

!===============================================================================
! 3.1 DONNEES SUR LES FACES FRONTIERES
!===============================================================================

!     On utilise TBORD comme auxiliaire pour l'appel a USRAY2
!       pour être sur que TPAROI ne sera pas modifié
!       (puisqu'on a TBORD libre)
!     On utilise FLUNET comme auxiliaire pour l'appel a USRAY2
!       pour être sur que QINCID ne sera pas modifié
!       (puisqu'on a FLUNET libre)

  do ifac = 1,nfabor
    tbord (ifac)     = propfb(ifac,ipprob(itparo))
    propfb(ifac,ipprob(ifnet)) = propfb(ifac,ipprob(iqinci))
  enddo

!     - Interface Code_Saturne
!       ======================

  if (iihmpr.eq.1) then

    call uiray2                                                   &
    !==========
  ( itypfb, iparoi, iparug, ivart , izfrad,                       &
    isothm, itpimp, ipgrno, iprefl, ifgrno, ifrefl,               &
    nozppm, nfabor, nvar,                                         &
    propfb(1,ipprob(ieps)), propfb(1,ipprob(iepa)), tint, text,   &
    propfb(1,ipprob(ixlam)), rcodcl)

  endif

  ils    = idebia
  idbia1 = ils + maxelt
  CALL iasize('raycli', idbia1)

  call usray2                                                     &
  !==========
 ( idbia1 , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , irapha  ,                                    &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   maxelt , ia(ils),                                              &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , izfrad , isothm ,                                     &
   tmin   , tmax   , tx     ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &
   tbord  , propfb(1,ipprob(ifnet)) ,  propfb(1,ipprob(ifconv))  ,&
   propfb(1,ipprob(ifconv)) , propfb(1,ipprob(ixlam)),            &
   propfb(1,ipprob(iepa))   , propfb(1,ipprob(ieps)) ,            &
   text   , tint   ,                                              &
   ra     )

!===============================================================================
! 3.2 CONTROLE DES DONNEES UTILISATEUR
!===============================================================================

!--> Arret si le numero de zone est non renseigne ou mal renseigne

  iok = 0

  do ifac = 1, nfabor
    if (izfrad(ifac).le.0.or.izfrad(ifac).gt.nozrdm) then
      iok = iok + 1
      write(nfecra,2000)ifac,nozrdm,izfrad(ifac)
    endif
  enddo

  if(iok.ne.0) then
    call csexit (1)
    !==========
  endif

! --> On construit une liste des numeros des zones frontieres.
!           (liste locale au processeur, en parallele)
!     Stop si depassement.

  nzfrad = 0
  do ifac = 1, nfabor
    ifvu = 0
    do ii = 1, nzfrad
      if (ilzrad(ii).eq.izfrad(ifac)) then
        ifvu = 1
      endif
    enddo
    if(ifvu.eq.0) then
      nzfrad = nzfrad + 1
      if(nzfrad.le.nbzrdm) then
        ilzrad(nzfrad) = izfrad(ifac)
      else
        write(nfecra,2001) nbzrdm
        write(nfecra,2002)(ilzrad(ii),ii=1,nbzrdm)
        call csexit (1)
        !==========
      endif
    endif
  enddo

! ---> Plus grand numero de zone atteint

  izonem = 0
  do ii = 1, nzfrad
    izone = ilzrad(ii)
    izonem = max(izonem,izone)
  enddo
  if(irangp.ge.0) then
    call parcmx(izonem)
    !==========
  endif
  nozarm = izonem




! On verra si ca coute cher ou non.
!   Pour le moment on le fait tout le temps.
!        IF(IWARNI(IVART).GE.-1.OR.IPACLI.LE.3) THEN
  if(1.eq.1) then

    iok = 0

!--> Si en paroi ISOTHM non renseignee : stop
    do ifac = 1, nfabor
      if( (itypfb(ifac,irapha).eq.iparoi  .or.                    &
           itypfb(ifac,irapha).eq.iparug) .and.                   &
           isothm(ifac)  .eq.-1    ) then
        iok = iok + 1
        write(nfecra,2110) irapha,ifac,izfrad(ifac)
      endif
    enddo

!--> Si ISOTHM renseignee en non paroi : stop
    do ifac = 1, nfabor
      if( itypfb(ifac,irapha).ne.iparoi .and.                     &
          itypfb(ifac,irapha).ne.iparug .and.                     &
          isothm(ifac)  .ne.-1         ) then
        iok = iok + 1
        write(nfecra,2111)                                        &
             irapha,ifac,izfrad(ifac),isothm(ifac)
      endif
    enddo

!--> Si valeur physique erronee : stop
    do ifac = 1, nfabor
      if(isothm(ifac).eq.itpimp ) then
        if(propfb(ifac,ipprob(ieps)) .lt.0.d0.or.                 &
            propfb(ifac,ipprob(ieps)).gt.1.d0.or.                 &
           tint(ifac).le.0.d0                      ) then
          iok = iok + 1
          write(nfecra,2120) irapha,ifac,izfrad(ifac),            &
               propfb(ifac,ipprob(ieps)),                         &
                              tint(ifac)
        endif
      elseif(isothm(ifac).eq.ipgrno ) then
        if(propfb(ifac,ipprob(ieps)) .lt.0.d0.or.                 &
            propfb(ifac,ipprob(ieps)).gt.1.d0.or.                 &
           propfb(ifac,ipprob(ixlam)).le.0.d0.or.                 &
           propfb(ifac,ipprob(iepa)) .le.0.d0.or.                 &
           text(ifac).le.0.d0.or.                                 &
           tint(ifac).le.0.d0                      ) then
          iok = iok + 1
          write(nfecra,2130) irapha,ifac,izfrad(ifac),            &
               propfb(ifac,ipprob(ieps)) ,                        &
               propfb(ifac,ipprob(ixlam)),                        &
               propfb(ifac,ipprob(iepa)) ,                        &
               text(ifac),tint(ifac)
        endif
      elseif(isothm(ifac).eq.iprefl ) then
        if(propfb(ifac,ipprob(ixlam)).le.0.d0.or.                 &
           propfb(ifac,ipprob(iepa)) .le.0.d0.or.                 &
           text(ifac).le.0.d0.or.                                 &
           tint(ifac).le.0.d0                      ) then
          iok = iok + 1
          write(nfecra,2140) irapha,ifac,izfrad(ifac),            &
                             propfb(ifac,ipprob(ixlam))    ,      &
                             propfb(ifac,ipprob(iepa))     ,      &
               text(ifac),tint(ifac)
        endif
      elseif(isothm(ifac).eq.ifgrno ) then
        if(propfb(ifac,ipprob(ieps)).lt.0.d0.or.                  &
           propfb(ifac,ipprob(ieps)).gt.1.d0.or.                  &
           tint(ifac).le.0.d0                      ) then
          iok = iok + 1
          write(nfecra,2150) irapha,ifac,izfrad(ifac),            &
               propfb(ifac,ipprob(ieps)),                         &
                              tint(ifac)
        endif
      elseif(isothm(ifac).eq.ifrefl ) then
        if(tint(ifac).le.0.d0                      ) then
          iok = iok + 1
          write(nfecra,2160) irapha,ifac,izfrad(ifac),            &
                              tint(ifac)
        endif
      elseif(isothm(ifac).ne.-1) then
          iok = iok + 1
          write(nfecra,2170) irapha,ifac,izfrad(ifac),            &
                             isothm(ifac)
      endif
    enddo

!--> Si valeur renseignee sans raison : stop
    do ifac = 1, nfabor
     if(isothm(ifac).eq.itpimp ) then
        if(propfb(ifac,ipprob(ixlam)).gt.0.d0.or.                 &
           propfb(ifac,ipprob(iepa))  .gt.0.d0.or.                &
           text(ifac).gt.0.d0                      ) then
          iok = iok + 1
          write(nfecra,2220) irapha,ifac,izfrad(ifac),            &
               propfb(ifac,ipprob(ixlam)),                        &
               propfb(ifac,ipprob(iepa)) ,text(ifac)
        endif
      elseif(isothm(ifac).eq.iprefl ) then
        if(propfb(ifac,ipprob(ieps)).ge.0.d0             ) then
          iok = iok + 1
          write(nfecra,2240) irapha,ifac,izfrad(ifac),            &
               propfb(ifac,ipprob(ieps))
        endif
      elseif(isothm(ifac).eq.ifgrno ) then
        if(propfb(ifac,ipprob(ixlam)).gt.0.d0.or.                 &
           propfb(ifac,ipprob(iepa)) .gt.0.d0.or.                 &
           text(ifac).gt.0.d0                      ) then
          iok = iok + 1
          write(nfecra,2250) irapha,ifac,izfrad(ifac),            &
               propfb(1,ipprob(ixlam)),propfb(1,ipprob(iepa)),    &
               text(ifac)
        endif
      elseif(isothm(ifac).eq.ifrefl ) then
        if(propfb(ifac,ipprob(ieps)) .ge.0.d0.or.                 &
           propfb(ifac,ipprob(ixlam)).gt.0.d0.or.                 &
           propfb(ifac,ipprob(iepa)) .gt.0.d0.or.                 &
           text(ifac).gt.0.d0                      ) then
          iok = iok + 1
          write(nfecra,2260) irapha,ifac,izfrad(ifac),            &
               propfb(ifac,ipprob(ieps)) ,                        &
               propfb(ifac,ipprob(ixlam)),                        &
               propfb(ifac,ipprob(iepa)) ,text(ifac)
        endif
      endif
    enddo

!--> Stop si erreur
    if(iok.ne.0) then
      call csexit (1)
      !==========
    endif

  endif


!===============================================================================
! 3.2 COMPLETION DES DONNEES UTILISATEUR
!===============================================================================

! ICODCL et EPS (quand il est nul)

  do ifac = 1, nfabor
    if(    isothm(ifac).eq.itpimp ) then
      icodcl(ifac,ivart) = 5
    elseif(isothm(ifac).eq.ipgrno ) then
      icodcl(ifac,ivart) = 5
    elseif(isothm(ifac).eq.iprefl ) then
      icodcl(ifac,ivart) = 5
      propfb(ifac,ipprob(ieps)) = 0.d0
    elseif(isothm(ifac).eq.ifgrno ) then
      icodcl(ifac,ivart) = 5
    elseif(isothm(ifac).eq.ifrefl ) then
      icodcl(ifac,ivart) = 3
      propfb(ifac,ipprob(ieps)) = 0.d0
    endif
  enddo


!===============================================================================
! 4. STOCKAGE DE LA TEMPERATURE (en Kelvin) dans TEMPK(IEL)
!===============================================================================

  if (abs(iscsth(iscat)).eq.1) then

!---> ON REMPLIT TEMPK

    if (iscsth(iscat).eq.-1) then
      do iel = 1, ncel
        tempk(iel) = rtpa(iel,ivart) + tkelvi
      enddo
    else
      do iel = 1, ncel
        tempk(iel) = rtpa(iel,ivart)
      enddo
    endif

  else if (iscsth(iscat).eq.2) then

!---> LECTURES DES DONNEES UTILISATEURS (TBORD est un auxiliaire)

    mode = 1

    if (ippmod(iphpar).le.1) then

      call usray4                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , irapha  ,                                    &
   nideve , nrdeve , nituse , nrtuse ,                            &
   mode   ,                                                       &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,irapha) ,&
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   propfb(1,ipprob(itparo)) , tbord  , tempk  ,                   &
!                                   Resultat : T en K
   rdevel , rtuser ,                                              &
   ra     )

    else

      call ppray4                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , irapha  ,                                    &
   nideve , nrdeve , nituse , nrtuse ,                            &

   mode   ,                                                       &

   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,irapha) ,&
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   propfb(1,ipprob(itparo)) , tbord  , tempk  ,                   &
!                                   Resultat : T en K
   rdevel , rtuser ,                                              &
   ra     )

    endif

  endif

!===============================================================================
! 5.  CALCUL DES TEMPERATURES DE PAROIS
!===============================================================================


! DANS TOUS LES CAS HFCONV CONTIENT Lambda * Hturb / distance
!   (HFCONV : W/(m2 K) ; Hturb est sans dimension)
!  (au premier passage, il est nul)

!--> CALCUL DU FLUX CONVECTIF
!      Par flux convectif, on entend bien sur
!        flux convectif parallele a la paroi,
!        on suppose que la paroi est etanche...
!      Le flux est calcule dans condli clptur, sauf au premier
!        passage sans suite de calcul, puisque raycli est appele avant.


if (ideb.eq.1) then

  do ifac = 1,nfabor
    if (isothm(ifac).ne.-1) then
      propfb(ifac,ipprob(ifconv)) =                               &
      propfb(ifac,ipprob(ihconv))*(tempk(ifabor(ifac))-           &
      propfb(ifac,ipprob(itparo)))
    endif
  enddo

endif


!--> Les cas ou il faut calculer TPAROI sont, au premier passage sans suite
!      des cas a temperature imposee TPAROI = TINT

  if (ideb.eq.1) then

    do ifac = 1,nfabor
      if (isothm(ifac).eq.ipgrno .or.                             &
          isothm(ifac).eq.iprefl .or.                             &
          isothm(ifac).eq.ifgrno    ) then
        isothm(ifac) = itpimp
      endif
    enddo

  endif

  if(ideb.eq.0) then

    call raypar                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , irapha  ,                                    &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,irapha) ,&
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &

   icodcl , isothm , izfrad ,                                     &

   idevel , ituser , ia     ,                                     &

   tmin   , tmax   , tx     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  ,                                              &
   rdevel , rtuser ,                                              &

   propfb(1,ipprob(itparo)) , propfb(1,ipprob(iqinci)) ,          &
   text   , tint   ,                                              &
   propfb(1,ipprob(ixlam))  , propfb(1,ipprob(iepa))   ,          &
   propfb(1,ipprob(ieps))   , propfb(1,ipprob(ihconv)) ,          &
   propfb(1,ipprob(ifconv)) , tempk  ,                            &

   ra     )

  endif


!===============================================================================
! 6.  CHANGEMENT DES CONDITIONS LIMITES UTILISATEUR
!===============================================================================

!===============================================================================
! 6.1  LA VARIABLE TRANSPORTEE EST LA TEMPERATURE
!===============================================================================

  if (abs(iscsth(iscat)).eq.1) then

    if(iscsth(iscat).eq.-1) then
      xmtk = -tkelvi
    else
      xmtk = 0.d0
    endif

    do ifac = 1,nfabor

      if (isothm(ifac).eq.itpimp .or.                             &
          isothm(ifac).eq.ipgrno .or.                             &
          isothm(ifac).eq.ifgrno    ) then
        rcodcl(ifac,ivart,1) = propfb(ifac,ipprob(itparo))+xmtk
        rcodcl(ifac,ivart,2) = rinfin
        rcodcl(ifac,ivart,3) = 0.d0

      else if (isothm(ifac).eq.iprefl) then
        rcodcl(ifac,ivart,1) = text(ifac)+xmtk
        rcodcl(ifac,ivart,2) = propfb(ifac,ipprob(ixlam))/        &
                                propfb(ifac,ipprob(iepa))
        rcodcl(ifac,ivart,3) = 0.d0

      else if (isothm(ifac).eq.ifrefl) then
        icodcl(ifac,ivart) = 3
        rcodcl(ifac,ivart,1) = 0.d0
        rcodcl(ifac,ivart,2) = rinfin
      endif

    enddo

!===============================================================================
! 6.2  LA VARIABLE TRANSPORTEE EST L'ENTHALPIE
!===============================================================================

  elseif (iscsth(iscat).eq.2) then

!---> LECTURES DES DONNEES UTILISATEURS
!     ON CONVERTIT TPAROI EN ENTHALPIE DE BORD, STOCKEE DANS FLUNET,
!     QUI EST UTILISE COMME AUXILIAIRE

    mode = 0

    do ifac = 1,nfabor
      if (isothm(ifac).eq.itpimp.or.                              &
          isothm(ifac).eq.ipgrno.or.                              &
          isothm(ifac).eq.ifgrno    ) then
        mode = -1
      endif
    enddo

    if (mode.eq.-1) then

      if (ippmod(iphpar).le.1) then

        call usray4                                               &
        !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , irapha  ,                                    &
   nideve , nrdeve , nituse , nrtuse ,                            &
   mode   ,                                                       &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,irapha) ,&
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   propfb(1,ipprob(itparo)) , propfb(1,ipprob(ifnet))  ,          &
   tempk  ,                                                       &
!                          HPAROI
   rdevel , rtuser ,                                              &
   ra     )

      else

        call ppray4                                               &
        !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , irapha  ,                                    &
   nideve , nrdeve , nituse , nrtuse ,                            &

   mode   ,                                                       &

   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,irapha) ,&
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   propfb(1,ipprob(itparo)) , propfb(1,ipprob(ifnet))  ,          &
   tempk  ,                                                       &
!                          HPAROI
   rdevel , rtuser ,                                              &
   ra     )

      endif

    endif

    mode = 0

    do ifac = 1,nfabor
      if (isothm(ifac).eq.iprefl) then
        mode = -1
      endif
    enddo

    if (mode.eq.-1) then

      if (ippmod(iphpar).le.1) then

        call usray4                                               &
        !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , irapha  ,                                    &
   nideve , nrdeve , nituse , nrtuse ,                            &
   mode   ,                                                       &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,irapha) ,&
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   text   , tbord  , tempk  ,                                     &
!                       HEXT
   rdevel , rtuser ,                                              &
   ra     )

      else

        call ppray4                                               &
        !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , irapha  ,                                    &
   nideve , nrdeve , nituse , nrtuse ,                            &

   mode   ,                                                       &

   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,irapha) ,&
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   text   , tbord  , tempk  ,                                     &
!                       HEXT
   rdevel , rtuser ,                                              &
   ra     )

      endif

    endif

    do ifac = 1,nfabor

      if (isothm(ifac).eq.itpimp.or.                              &
          isothm(ifac).eq.ipgrno.or.                              &
          isothm(ifac).eq.ifgrno    ) then
        rcodcl(ifac,ivart,1) = propfb(ifac,ipprob(ifnet))
        rcodcl(ifac,ivart,2) = rinfin
        rcodcl(ifac,ivart,3) = 0.d0

      else if (isothm(ifac).eq.iprefl) then

        if (icp(irapha).gt.0) then
          iel = ifabor(ifac)
          cpp = propce(iel,ipproc(icp(irapha)))
        else
          cpp = cp0(irapha)
        endif

        rcodcl(ifac,ivart,1) = tbord(ifac)
        rcodcl(ifac,ivart,2) =  propfb(ifac,ipprob(ixlam))     &
                             / (propfb(ifac,ipprob(iepa))*cpp)
        rcodcl(ifac,ivart,3) = 0.d0

      else if (isothm(ifac).eq.ifrefl) then
        icodcl(ifac,ivart) = 3
        rcodcl(ifac,ivart,1) = 0.d0
        rcodcl(ifac,ivart,2) = rinfin
      endif

    enddo

  endif

!--------
! FORMATS
!--------

 1000 FORMAT (/, 3X,'** INFORMATIONS SUR LE MODULE DE RAYONNEMENT ',/,  &
           3X,'   ------------------------------------------',/,  &
           3X,' Initialisation de la temperature de paroi   ',/,  &
           3X,' (TPAROI) avec le profil utilisateur (TINTP) ',/,  &
           3X,' et du flux incident aux parois (QINCID).    ',/)

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES ',/,&
'@                                                            ',/,&
'@  Le numero de zone associee a la face ',I10   ,' doit etre ',/,&
'@    un entier strictement positif et inferieur ou egal a    ',/,&
'@    NOZRDM = ',I10                                           ,/,&
'@  Ce numero (IZFRDP(IFAC)) vaut ici ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites dans usray2.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@  Le nombre maximal de zones frontieres qui peuvent etre    ',/,&
'@    definies par l''utilisateur est NBZRDM = ',I10           ,/,&
'@    Il a ete depasse.                                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites dans usray2.          ',/,&
'@                                                            ',/,&
'@  Les NBZRDM premieres zones frontieres                     ',/,&
'@    portent ici les numeros suivants :                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2002 format(i10)

 2110 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    ISOTHP DOIT ETRE RENSEIGNE SUR TOUTES LES FACES DE PAROI',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Il ne l''a pas ete pour la face ',I10                      ,/,&
'@                    zone         ',I10                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2111 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    ISOTHP A ETE RENSEIGNE SUR UNE FACE NON PAROI           ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Sur la face ',I10   ,', zone  ',I10   ,', ISOTHP a ete    ',/,&
'@    renseigne dans usray2 (ISOTHP = ',I10   ,') alors que   ',/,&
'@    la face n''a pas ete declaree de type IPAROI ou IPARUG  ',/,&
'@    dans usclim.                                            ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2 et usclim.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2120 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP  DOIT ETRE UN REEL INCLUS DANS [0.; 1.]            ',/,&
'@    TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF             ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = ITPIMP       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@  TINTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2130 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP DOIT ETRE UN REEL INCLUS DANS [0.; 1.]             ',/,&
'@    XLAMP, EPAP, TINTP, TEXTP DOIVENT ETRE DES REELS        ',/,&
'@                                      STRICTEMENT POSITIFS  ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IPGRNO       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,&
'@  TEXTP = ',E14.5    ,' TINTP = ',E14.5                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2140 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP, TINTP, TEXTP DOIVENT ETRE DES REELS        ',/,&
'@                                      STRICTEMENT POSITIFS  ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IPREFL       ',/,&
'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,&
'@  TEXTP = ',E14.5    ,' TINTP = ',E14.5                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2150 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP  DOIT ETRE UN REEL INCLUS DANS [0.; 1.]            ',/,&
'@    TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF             ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFGRNO       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@  TINTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2160 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@  TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF               ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFREFL       ',/,&
'@  TINTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2170 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@   VALEUR NON ADMISSIBLE DE ISOTHP                          ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = ',I10         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2220 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP ET TEXTP NE DOIVENT PAS ETRE RENSEIGNES     ',/,&
'@                                     AVEC ISOTHP = ITPIMP   ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = ITPIMP       ',/,&
'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,&
'@  TEXTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2240 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP NE DOIT PAS ETRE RENSEIGNE AVEC ISOTHP = IPREFL    ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IPREFL       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2250 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP, TEXTP NE DOIVENT PAS ETRE RENSEIGNES       ',/,&
'@                                       AVEC ISOTHP = IFGRNO ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFGRNO       ',/,&
'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,&
'@  TEXTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2260 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP, TEXTP NE DOIVENT PAS ETRE RENSEIGNES       ',/,&
'@                                       AVEC ISOTHP = IFREFL ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFREFL       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,&
'@  TEXTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return

end subroutine
