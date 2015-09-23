!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine raycli &
!================

 ( nvar   , nscal  ,                                              &
   isvhb  , isvtb  ,                                              &
   icodcl , itrifb , itypfb ,                                     &
   izfrad ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  , hbord  )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! isvhb            ! e  ! <-- ! indicateur de sauvegarde des                   !
!                  !    !     !  coefficients d'echange aux bords              !
! isvtb            ! e  ! <-- ! indicateur de sauvegarde des                   !
!                  !    !     !  temperatures aux bords                        !
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
! itypfb           ! ia ! --> ! boundary face types                            !
! izfrad(nfabor    ! te ! <-- ! numero de zone des faces de bord               !
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
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use ihmpre
use ppppar
use ppthch
use ppincl
use radiat
use dimens, only: ndimfb
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          isvhb  , isvtb

integer          icodcl(ndimfb,nvarcl)
integer          itrifb(ndimfb), itypfb(ndimfb)
integer          izfrad(ndimfb)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision rcodcl(ndimfb,nvarcl,3)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision hbord(ndimfb)

! Local variables
integer          ifac, iel, ideb, ivart
integer          mode, iok, ifvu, ii, izonem, izone

double precision tmin , tmax   , tx
double precision xmtk

integer, allocatable, dimension(:) :: isothm

double precision, allocatable, dimension(:) :: tempk, thwall
double precision, allocatable, dimension(:) :: text, tint

integer    ipacli
data       ipacli /0/
save       ipacli

!===============================================================================
! 1. Initializations
!===============================================================================

! Allocate temporary arrays
allocate(isothm(nfabor))
allocate(tempk(ncelet), thwall(ndimfb))
allocate(text(nfabor), tint(nfabor))

!---> NUMERO DE PASSAGE RELATIF
ipacli = ipacli + 1
ideb = 0

!---> Min and Max values of the temperature (in Kelvin)
tmin = 0.d0
tmax = grand + tkelvi

!---> COEFF DE RELAX

!     tx est strictement superieur a 0 et inferieur ou egal a 1

!     Pour calculer la temperature de paroi, on calcule un increment
!     de temperature DeltaT entre l'etape courante n et l'etape
!     precedente n-1, puis on calcule :
!          n    n-1                                 n-1
!         T  = T    + DeltaT si le rapport DeltaT/T    =< tx, sinon

!          n    n-1                      n-1             n-1
!         T  = T    * (1 + tx *((DeltaT/T   ) / |DeltaT/T   |))

tx = 0.1d0

!---> Default initialization
do ifac = 1, nfabor
  izfrad(ifac) = -1
  isothm(ifac) = -1
  propfb(ifac,ipprob(ixlam)) = -grand
  propfb(ifac,ipprob(iepa))  = -grand
  propfb(ifac,ipprob(ieps))  = -grand
  text  (ifac) = -grand
  tint  (ifac) = -grand
enddo

! Index of the thermal variable
ivart = isca(iscalt)

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
  !       pour etre sur que TPAROI ne sera pas modifie
  !       (puisqu'on a TBORD libre)
  !     On utilise FLUNET comme auxiliaire pour l'appel a USRAY2
  !       pour etre sur que QINCID ne sera pas modifie
  !       (puisqu'on a FLUNET libre)

  do ifac = 1, nfabor
    thwall(ifac) = zero
    propfb(ifac,ipprob(ifnet)) = zero
  enddo

  ! - Interface Code_Saturne
  !   ======================

  if (iihmpr.eq.1) then

    !---> NUMERO DU SCALAIRE ET DE LA VARIABLE THERMIQUE

    call uiray2 &
    !==========
   ( itypfb, iparoi, iparug, ivart , izfrad,                  &
     isothm, itpimp, ipgrno, iprefl, ifgrno, ifrefl,          &
     nozppm, nfabor, nvar,                                    &
     propfb(1,ipprob(ieps)), propfb(1,ipprob(iepa)),          &
     tint, text,                                              &
     propfb(1,ipprob(ixlam)), rcodcl)

  endif

  call usray2 &
  !==========
( nvar   , nscal  ,                                              &
  itypfb ,                                                       &
  icodcl , izfrad , isothm ,                                     &
  tmin   , tmax   , tx     ,                                     &
  dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
  thwall , propfb(1,ipprob(ifnet))  , propfb(1,ipprob(ihconv))  ,&
  propfb(1,ipprob(ifconv)),                                      &
  propfb(1,ipprob(ixlam)) , propfb(1,ipprob(iepa)) ,             &
  propfb(1,ipprob(ieps))  ,                                      &
  text   , tint   )

  write(nfecra,1000)

  ! Tparoi en Kelvin et QINCID en W/m2
  do ifac = 1, nfabor
    propfb(ifac,ipprob(itparo)) = tint(ifac)
    propfb(ifac,ipprob(iqinci)) = stephn*tint(ifac)**4
    if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then
      propfb(ifac,ipprob(itparo)) = tint(ifac)
      propfb(ifac,ipprob(iqinci)) = stephn*tint(ifac)**4
    else
      propfb(ifac,ipprob(itparo)) = 0.d0
      propfb(ifac,ipprob(iqinci)) = 0.d0
    endif
  enddo

endif

!===============================================================================
! 3.1 DONNEES SUR LES FACES FRONTIERES
!===============================================================================

!     On utilise flunet comme auxiliaire pour l'appel a USRAY2
!       pour etre sur que QINCID ne sera pas modifie
!       (puisqu'on a flunet libre)

do ifac = 1, nfabor
  thwall (ifac) = propfb(ifac,ipprob(itparo))
  propfb(ifac,ipprob(ifnet)) = propfb(ifac,ipprob(iqinci))
enddo

!     - Interface Code_Saturne
!       ======================

if (iihmpr.eq.1) then

  call uiray2 &
  !==========
( itypfb, iparoi, iparug, ivart , izfrad,                       &
  isothm, itpimp, ipgrno, iprefl, ifgrno, ifrefl,               &
  nozppm, nfabor, nvar,                                         &
  propfb(1,ipprob(ieps)), propfb(1,ipprob(iepa)), tint, text,   &
  propfb(1,ipprob(ixlam)), rcodcl)

endif

call usray2 &
!==========
( nvar   , nscal  ,                                              &
  itypfb ,                                                       &
  icodcl , izfrad , isothm ,                                     &
  tmin   , tmax   , tx     ,                                     &
  dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
  thwall , propfb(1,ipprob(ifnet)) ,  propfb(1,ipprob(ihconv)) , &
  propfb(1,ipprob(ifconv)) , propfb(1,ipprob(ixlam)),            &
  propfb(1,ipprob(iepa))   , propfb(1,ipprob(ieps)) ,            &
  text   , tint   )

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
  if (ifvu.eq.0) then
    nzfrad = nzfrad + 1
    if (nzfrad.le.nbzrdm) then
      ilzrad(nzfrad) = izfrad(ifac)
    else
      write(nfecra,2001) nbzrdm
      write(nfecra,2002)(ilzrad(ii),ii=1,nbzrdm)
      call csexit (1)
    endif
  endif
enddo

! ---> Plus grand numero de zone atteint

izonem = 0
do ii = 1, nzfrad
  izone = ilzrad(ii)
  izonem = max(izonem,izone)
enddo
if (irangp.ge.0) then
  call parcmx(izonem)
endif
nozarm = izonem

! On verra si ca coute cher ou non.
!   Pour le moment on le fait tout le temps.
!        IF(IWARNI(IVART).GE.-1.OR.IPACLI.LE.3) THEN
if (1.eq.1) then

  iok = 0

  !--> Si en paroi ISOTHM non renseignee : stop
  do ifac = 1, nfabor
    if ((itypfb(ifac).eq.iparoi  .or.                            &
         itypfb(ifac).eq.iparug) .and.                           &
         isothm(ifac).eq.-1) then
      iok = iok + 1
      write(nfecra,2110) ifac,izfrad(ifac)
    endif
  enddo

  !--> Si ISOTHM renseignee en non paroi : stop
  do ifac = 1, nfabor
    if (itypfb(ifac).ne.iparoi .and.                             &
        itypfb(ifac).ne.iparug .and.                             &
        isothm(ifac)  .ne.-1         ) then
      iok = iok + 1
      write(nfecra,2111) ifac, izfrad(ifac), isothm(ifac)
    endif
  enddo

  !--> Si valeur physique erronee : stop
  do ifac = 1, nfabor
    if (isothm(ifac).eq.itpimp) then
      if (propfb(ifac,ipprob(ieps)).lt.0.d0.or.                  &
          propfb(ifac,ipprob(ieps)).gt.1.d0.or.                  &
          tint(ifac).le.0.d0                      ) then
        iok = iok + 1
        write(nfecra,2120) ifac,izfrad(ifac),                    &
                           propfb(ifac,ipprob(ieps)), tint(ifac)
      endif
    elseif (isothm(ifac).eq.ipgrno) then
      if (propfb(ifac,ipprob(ieps)) .lt.0.d0.or.                 &
          propfb(ifac,ipprob(ieps)) .gt.1.d0.or.                 &
          propfb(ifac,ipprob(ixlam)).le.0.d0.or.                 &
          propfb(ifac,ipprob(iepa)) .le.0.d0.or.                 &
          text(ifac).le.0.d0.or.                                 &
          tint(ifac).le.0.d0                      ) then
        iok = iok + 1
        write(nfecra,2130) ifac,izfrad(ifac),                    &
                           propfb(ifac,ipprob(ieps)) ,           &
                           propfb(ifac,ipprob(ixlam)),           &
                           propfb(ifac,ipprob(iepa)) ,           &
                           text(ifac),tint(ifac)
      endif
    elseif (isothm(ifac).eq.iprefl) then
      if (propfb(ifac,ipprob(ixlam)).le.0.d0.or.                 &
          propfb(ifac,ipprob(iepa)) .le.0.d0.or.                 &
          text(ifac).le.0.d0.or.                                 &
          tint(ifac).le.0.d0                      ) then
        iok = iok + 1
        write(nfecra,2140) ifac,izfrad(ifac),                    &
                           propfb(ifac,ipprob(ixlam)),           &
                           propfb(ifac,ipprob(iepa)),            &
                           text(ifac),tint(ifac)
      endif
    elseif (isothm(ifac).eq.ifgrno) then
      if (propfb(ifac,ipprob(ieps)).lt.0.d0.or.                  &
          propfb(ifac,ipprob(ieps)).gt.1.d0.or.                  &
          tint(ifac).le.0.d0                      ) then
        iok = iok + 1
        write(nfecra,2150) ifac, izfrad(ifac),                   &
             propfb(ifac,ipprob(ieps)), tint(ifac)
      endif
    elseif (isothm(ifac).eq.ifrefl) then
      if (tint(ifac).le.0.d0) then
        iok = iok + 1
        write(nfecra,2160) ifac, izfrad(ifac), tint(ifac)
      endif
    elseif (isothm(ifac).ne.-1) then
        iok = iok + 1
        write(nfecra,2170) ifac, izfrad(ifac), isothm(ifac)
    endif
  enddo

  !--> Si valeur renseignee sans raison : stop
  do ifac = 1, nfabor
   if (isothm(ifac).eq.itpimp) then
      if (propfb(ifac,ipprob(ixlam)).gt.0.d0.or.                &
          propfb(ifac,ipprob(iepa)) .gt.0.d0.or.                &
          text(ifac).gt.0.d0                      ) then
        iok = iok + 1
        write(nfecra,2220) ifac,izfrad(ifac),                   &
                           propfb(ifac,ipprob(ixlam)),          &
                           propfb(ifac,ipprob(iepa)), text(ifac)
      endif
    elseif (isothm(ifac).eq.iprefl) then
      if (propfb(ifac,ipprob(ieps)).ge.0.d0) then
        iok = iok + 1
        write(nfecra,2240) ifac, izfrad(ifac), propfb(ifac,ipprob(ieps))
      endif
    elseif (isothm(ifac).eq.ifgrno) then
      if (propfb(ifac,ipprob(ixlam)).gt.0.d0.or.                &
          propfb(ifac,ipprob(iepa)) .gt.0.d0.or.                &
          text(ifac).gt.0.d0                      ) then
        iok = iok + 1
        write(nfecra,2250) ifac,izfrad(ifac),                   &
                           propfb(1,ipprob(ixlam)),             &
                           propfb(1,ipprob(iepa)), text(ifac)
      endif
    elseif (isothm(ifac).eq.ifrefl) then
      if(propfb(ifac,ipprob(ieps)) .ge.0.d0.or.                 &
         propfb(ifac,ipprob(ixlam)).gt.0.d0.or.                 &
         propfb(ifac,ipprob(iepa)) .gt.0.d0.or.                 &
         text(ifac).gt.0.d0                      ) then
        iok = iok + 1
        write(nfecra,2260) ifac,izfrad(ifac),                   &
                           propfb(ifac,ipprob(ieps)),           &
                           propfb(ifac,ipprob(ixlam)),          &
                           propfb(ifac,ipprob(iepa)), text(ifac)
      endif
    endif
  enddo

  if (iok.ne.0) then
    call csexit (1)
  endif

endif

!===============================================================================
! 3.2 COMPLETION DES DONNEES UTILISATEUR
!===============================================================================

! ICODCL et EPS (quand il est nul)

do ifac = 1, nfabor
  if (isothm(ifac).eq.itpimp) then
    icodcl(ifac,ivart) = 5
  elseif (isothm(ifac).eq.ipgrno) then
    icodcl(ifac,ivart) = 5
  elseif (isothm(ifac).eq.iprefl) then
    icodcl(ifac,ivart) = 5
    propfb(ifac,ipprob(ieps)) = 0.d0
  elseif (isothm(ifac).eq.ifgrno) then
    icodcl(ifac,ivart) = 5
  elseif (isothm(ifac).eq.ifrefl) then
    icodcl(ifac,ivart) = 3
    propfb(ifac,ipprob(ieps)) = 0.d0
  endif
enddo

!===============================================================================
! 4. STOCKAGE DE LA TEMPERATURE (en Kelvin) dans TEMPK(IEL)
!===============================================================================

if (abs(iscsth(iscalt)).eq.1) then

  !---> ON REMPLIT TEMPK

  if (iscsth(iscalt).eq.-1) then
    do iel = 1, ncel
      tempk(iel) = rtpa(iel,ivart) + tkelvi
    enddo
  else
    do iel = 1, ncel
      tempk(iel) = rtpa(iel,ivart)
    enddo
  endif

  elseif (iscsth(iscalt).eq.2) then

    !---> LECTURES DES DONNEES UTILISATEURS (TBORD est un auxiliaire)

    mode = 1

    if (ippmod(iphpar).le.1) then

      call usray4 &
      !==========
 ( nvar   , nscal  ,                                              &
   mode   ,                                                       &
   itypfb ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   propfb(1,ipprob(itparo)) , thwall , tempk  )
      ! Resultat : T en K

    else

      call ppray4 &
      !==========
 ( nvar   , nscal  ,                                              &
   mode   ,                                                       &
   itypfb ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   propfb(1,ipprob(itparo)) , thwall , tempk  )
      ! Resultat : T en K

    endif

  endif

!===============================================================================
! 5. CALCUL DES TEMPERATURES DE PAROIS
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

  do ifac = 1, nfabor
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

  call raypar &
  !==========
 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   icodcl , isothm , izfrad ,                                     &
   tmin   , tmax   , tx     ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  ,                                              &
   propfb(1,ipprob(itparo)) , propfb(1,ipprob(iqinci)) ,          &
   text   , tint   ,                                              &
   propfb(1,ipprob(ixlam))  , propfb(1,ipprob(iepa))   ,          &
   propfb(1,ipprob(ieps))   , propfb(1,ipprob(ihconv)) ,          &
   propfb(1,ipprob(ifconv)) , tempk  )

endif

!===============================================================================
! 6.  CHANGEMENT DES CONDITIONS LIMITES UTILISATEUR
!===============================================================================

!===============================================================================
! 6.1  LA VARIABLE TRANSPORTEE EST LA TEMPERATURE
!===============================================================================

if (abs(iscsth(iscalt)).eq.1) then

  if (iscsth(iscalt).eq.-1) then
    xmtk = -tkelvi
  else
    xmtk = 0.d0
  endif

  do ifac = 1, nfabor

    if (isothm(ifac).eq.itpimp .or.                             &
        isothm(ifac).eq.ipgrno .or.                             &
        isothm(ifac).eq.ifgrno) then
      rcodcl(ifac,ivart,1) = propfb(ifac,ipprob(itparo))+xmtk
      rcodcl(ifac,ivart,2) = rinfin
      rcodcl(ifac,ivart,3) = 0.d0

    elseif (isothm(ifac).eq.iprefl) then
      rcodcl(ifac,ivart,1) = text(ifac)+xmtk
      rcodcl(ifac,ivart,2) = propfb(ifac,ipprob(ixlam))/        &
                              propfb(ifac,ipprob(iepa))
      rcodcl(ifac,ivart,3) = 0.d0

    elseif (isothm(ifac).eq.ifrefl) then
      icodcl(ifac,ivart) = 3
      rcodcl(ifac,ivart,1) = 0.d0
      rcodcl(ifac,ivart,2) = rinfin
    endif

  enddo

!===============================================================================
! 6.2  LA VARIABLE TRANSPORTEE EST L'ENTHALPIE
!===============================================================================

elseif (iscsth(iscalt).eq.2) then

  !---> LECTURES DES DONNEES UTILISATEURS
  !     ON CONVERTIT TPAROI EN ENTHALPIE DE BORD, STOCKEE DANS FLUNET,
  !     QUI EST UTILISE COMME AUXILIAIRE

  mode = 0

  do ifac = 1, nfabor
    if (isothm(ifac).eq.itpimp.or.                              &
        isothm(ifac).eq.ipgrno.or.                              &
        isothm(ifac).eq.ifgrno) then
      mode = -1
    endif
  enddo

  if (mode.eq.-1) then

    if (ippmod(iphpar).le.1) then

      call usray4 &
      !==========
    ( nvar   , nscal  ,                                              &
      mode   ,                                                       &
      itypfb ,                                                       &
      dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
      propfb(1,ipprob(itparo)) , propfb(1,ipprob(ifnet))  ,          &
      tempk  )
      ! HPAROI

    else

      call ppray4 &
      !==========
    ( nvar   , nscal  ,                                              &
      mode   ,                                                       &
      itypfb ,                                                       &
      dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
      coefa  , coefb  ,                                              &
      propfb(1,ipprob(itparo)) , propfb(1,ipprob(ifnet))  ,          &
      tempk  )
      ! HPAROI

    endif

  endif

  mode = 0

  do ifac = 1, nfabor
    if (isothm(ifac).eq.iprefl) then
      mode = -1
    endif
  enddo

  if (mode.eq.-1) then

    if (ippmod(iphpar).le.1) then

      call usray4 &
      !==========
    ( nvar   , nscal  ,                                              &
      mode   ,                                                       &
      itypfb ,                                                       &
      dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
      text   , thwall , tempk  )
      ! HEXT

    else

      call ppray4 &
      !==========
    ( nvar   , nscal  ,                                              &
      mode   ,                                                       &
      itypfb ,                                                       &
      dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
      coefa  , coefb  ,                                              &
      text   , thwall , tempk  )
      ! HEXT

    endif

  endif

  do ifac = 1, nfabor

    if (isothm(ifac).eq.itpimp.or.                              &
        isothm(ifac).eq.ipgrno.or.                              &
        isothm(ifac).eq.ifgrno) then
      rcodcl(ifac,ivart,1) = propfb(ifac,ipprob(ifnet))
      rcodcl(ifac,ivart,2) = rinfin
      rcodcl(ifac,ivart,3) = 0.d0

    elseif (isothm(ifac).eq.iprefl) then

      rcodcl(ifac,ivart,1) = thwall(ifac)
      ! hext
      rcodcl(ifac,ivart,2) =  propfb(ifac,ipprob(ixlam))     &
                           / (propfb(ifac,ipprob(iepa)))
      rcodcl(ifac,ivart,3) = 0.d0

    elseif (isothm(ifac).eq.ifrefl) then
      icodcl(ifac,ivart) = 3
      rcodcl(ifac,ivart,1) = 0.d0
      rcodcl(ifac,ivart,2) = rinfin
    endif

  enddo

endif

! Free memory
deallocate(isothm)
deallocate(tempk,thwall)
deallocate(text, tint)

!--------
! Formats
!--------

 1000 format (/, &
 3X,'** INFORMATIONS SUR LE MODULE DE RAYONNEMENT'             ,/,&
3X,'   ------------------------------------------'             ,/,&
3X,' Initialisation de la temperature de paroi'                ,/,&
3X,' (TPAROI) avec le profil utilisateur (TINTP)'              ,/,&
3X,' et du flux incident aux parois (QINCID).'                 ,/)

 2000 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : RAYONNEMENT'                                 ,/,&
'@    ========='                                               ,/,&
'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES' ,/,&
'@'                                                            ,/,&
'@  Le numero de zone associee a la face ',I10   ,' doit etre' ,/,&
'@    un entier strictement positif et inferieur ou egal a'    ,/,&
'@    NOZRDM = ',I10                                           ,/,&
'@  Ce numero (IZFRDP(IFAC)) vaut ici ',I10                    ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@'                                                            ,/,&
'@  Verifier les conditions aux limites dans usray2.'          ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2001 format(                                                     &
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

 2110 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    ISOTHP DOIT ETRE RENSEIGNE SUR TOUTES LES FACES DE PAROI',/,&
'@                                                            ',/,&
'@  Il ne l''a pas ete pour la face ',I10                      ,/,&
'@                    zone         ',I10                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2111 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    ISOTHP A ETE RENSEIGNE SUR UNE FACE NON PAROI           ',/,&
'@                                                            ',/,&
'@  Sur la face ',I10   ,', zone  ',I10   ,', ISOTHP a ete    ',/,&
'@    renseigne dans usray2 (ISOTHP = ',I10   ,') alors que   ',/,&
'@    la face n''a pas ete declaree de type IPAROI ou IPARUG  ',/,&
'@    dans cs_user_boundary_conditions.                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2 et cs_user_boundary_conditions.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2120 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP  DOIT ETRE UN REEL INCLUS DANS [0.; 1.]            ',/,&
'@    TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF             ',/,&
'@                                                            ',/,&
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

 2130 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP DOIT ETRE UN REEL INCLUS DANS [0.; 1.]             ',/,&
'@    XLAMP, EPAP, TINTP, TEXTP DOIVENT ETRE DES REELS        ',/,&
'@                                      STRICTEMENT POSITIFS  ',/,&
'@                                                            ',/,&
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
 2140 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP, TINTP, TEXTP DOIVENT ETRE DES REELS        ',/,&
'@                                      STRICTEMENT POSITIFS  ',/,&
'@                                                            ',/,&
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
 2150 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP  DOIT ETRE UN REEL INCLUS DANS [0.; 1.]            ',/,&
'@    TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF             ',/,&
'@                                                            ',/,&
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
 2160 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@  TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF               ',/,&
'@                                                            ',/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFREFL       ',/,&
'@  TINTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2170 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@   VALEUR NON ADMISSIBLE DE ISOTHP                          ',/,&
'@                                                            ',/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = ',I10         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2220 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP ET TEXTP NE DOIVENT PAS ETRE RENSEIGNES     ',/,&
'@                                     AVEC ISOTHP = ITPIMP   ',/,&
'@                                                            ',/,&
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
 2240 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP NE DOIT PAS ETRE RENSEIGNE AVEC ISOTHP = IPREFL    ',/,&
'@                                                            ',/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IPREFL       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2250 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP, TEXTP NE DOIVENT PAS ETRE RENSEIGNES       ',/,&
'@                                       AVEC ISOTHP = IFGRNO ',/,&
'@                                                            ',/,&
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
 2260 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP, TEXTP NE DOIVENT PAS ETRE RENSEIGNES       ',/,&
'@                                       AVEC ISOTHP = IFREFL ',/,&
'@                                                            ',/,&
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
! End
!----

return

end subroutine
