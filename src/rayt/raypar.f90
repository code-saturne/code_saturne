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

subroutine raypar &
!================

 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   icodcl , isothp , izfrap ,                                     &
   tmin   , tmax   , tx     ,                                     &
   dt     , rtp    , rtpa   , propce , propfb , rcodcl ,          &
   coefa  , coefb  ,                                              &
   tparop , qincip , textp  , tintp  ,                            &
   xlamp  , epap   , epsp   ,                                     &
   hfconp , flconp , tempkp )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
!   -----------------------------------------

!   CALCUL DES TEMPERATURES DE PAROIS AVEC UN BILAN DE FLUX

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! isothp(nfabor    ! te ! <-- ! liste des frontieres isothermes                !
! izfrap(nfabor    ! te ! <-- ! numero de zone des faces de bord               !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
! tparop(nfabor    ! tr ! --> ! temperature de paroi en kelvin                 !
! qincip(nfabor    ! tr ! <-- ! densite de flux radiatif aux bords             !
! textp(nfabor)    ! tr ! <-- ! temperature de bord externe                    !
!                  !    !     ! en degres celcius                              !
! tintp(nfabor)    ! tr ! <-- ! temperature de bord interne                    !
!                  !    !     ! en degres celcius                              !
! xlamp(nfabor)    ! tr ! <-- ! coefficient de conductivite thermique          !
!                  !    !     ! des facettes de paroi (w/m/k)                  !
! epap(nfabor)     ! tr ! <-- ! epaisseur des facettes de paroi (m)            !
! epsp(nfabor)     ! tr ! <-- ! emissivite des facettes de bord                !
! hfconp(nfabor    ! tr ! <-- ! coefficient d'echange fluide aux               !
!                  !    !     ! faces de bord                                  !
! flconp(nfabor    ! tr ! <-- ! densite de flux convectif aux faces            !
! tempkp(ncelet    ! tr ! <-- ! temperature en kelvin                          !
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
use parall
use optcal
use cstphy
use cstnum
use ppppar
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          itypfb(nfabor)

integer          isothp(nfabor), izfrap(nfabor)
integer          icodcl(nfabor,nvarcl)

double precision tmin , tmax , tx
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvarcl,3)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision tparop(nfabor), qincip(nfabor)
double precision textp(nfabor), tintp(nfabor)
double precision xlamp(nfabor), epap(nfabor), epsp(nfabor)
double precision hfconp(nfabor) , flconp(nfabor)
double precision tempkp(ncelet)

! Local variables


integer          ivart, ifac, iel, izone
integer          ifacmx, ifacmn
integer          n1min,n1max,nrelax,nmoins,nplus
integer          iitpim,iipgrn,iipref,iifgrn,iifref
integer          indtp(nozrdm)
integer          nbrval, indtpm
integer          nb1int  ,nb2int  ,nbrrdp
parameter       (nb1int=5,nb2int=5,nbrrdp=5)
integer          inttm1(nb1int),inttm2(nb2int)

double precision esl,epp,sigt3,sigt4
double precision tpmin,tpmax,rapp,rapmax,abrapp
double precision qcmax,qrmax
double precision qcmin,qrmin
double precision qrayt,qconv,qinci,detep
double precision srfbn,und0,tp4
double precision xtpmax,ytpmax,ztpmax,xtpmin,ytpmin,ztpmin
double precision tzomax(nozrdm),tzomin(nozrdm),tzomoy(nozrdm)
double precision flunet(nozrdm),radios(nozrdm),surft(nozrdm)
double precision rdptmp(nbrrdp)


!===============================================================================

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================


und0 = 1.d0

ivart = isca(iscalt)


tpmax = -grand
tpmin =  grand

rapmax = 0.d0
n1min  = 0
n1max  = 0
nrelax = 0
nmoins = 0
nplus  = 0
iitpim = 0
iipgrn = 0
iipref = 0
iifgrn = 0
iifref = 0

ifacmx = 0
ifacmn = 0

!===============================================================================
! 2. Decompte du nombre de couleurs differentes de facettes de bord
!    On suppose qu'on en a au moins une
!===============================================================================


do izone = 1, nozrdm
  indtp (izone) =      0
  tzomax(izone) = -grand
  tzomin(izone) =  grand
enddo

!===============================================================================
! 3. Calcul des temperatures de paroi
!===============================================================================

!     3.1) parois isothermes
!     ----------------------


do ifac = 1, nfabor
  if (isothp(ifac).eq.itpimp) then
    izone = izfrap(ifac)

!      Calcul
    tparop(ifac)  = tintp(ifac)

!      Max-Min
    qconv  = flconp(ifac)
    qinci  = qincip(ifac)
    sigt4  = stephn*tparop(ifac)**4
    epp    = epsp(ifac)
    qrayt  = epp *( qinci - sigt4 )

    if (tpmax.le.tparop(ifac)) then
      ifacmx = ifac
      tpmax = tparop(ifac)
      qcmax = qconv
      qrmax = qrayt
    endif

    if (tpmin.ge.tparop(ifac)) then
      ifacmn = ifac
      tpmin = tparop(ifac)
      qcmin = qconv
      qrmin = qrayt
    endif

    tzomax(izone) = max(tzomax(izone),tparop(ifac))
    tzomin(izone) = min(tzomin(izone),tparop(ifac))

!      Reperage pour impression
    iitpim = 1
    indtp(izone) = itpimp


  endif
enddo


!     3.2) parois grises ou noires (non reflechissantes)
!     --------------------------------------------------

do ifac = 1, nfabor
  if (isothp(ifac).eq.ipgrno) then
    izone = izfrap(ifac)

!       Calcul
    esl    = epap(ifac) /xlamp(ifac)
    qconv  = flconp(ifac)
    qinci  = qincip(ifac)
    epp    = epsp(ifac)

    sigt3  = stephn*tparop(ifac)**3
    qrayt  = epp *( qinci -sigt3*tparop(ifac))

    detep  = ( esl*(qconv+qrayt) - (tparop(ifac)-textp(ifac)) )   &
           / ( 1.d0+ 4.d0*esl*epp*sigt3 + esl*hfconp(ifac) )

    rapp   = detep /tparop(ifac)
    abrapp = abs(rapp)

!       Relaxation
    if (abrapp.ge.tx) then
      nrelax = nrelax +1
      tparop(ifac) = tparop(ifac) * (1.d0+ tx*sign(und0,rapp))
    else
      tparop(ifac) = tparop(ifac) + detep
    endif

    rapmax = max(rapmax,abrapp)

    if (rapp.le.0.d0) then
      nmoins = nmoins+1
    else
      nplus  = nplus +1
    endif

!       Clipping
    if (tparop(ifac).lt.tmin) then
      n1min = n1min +1
      tparop(ifac) = tmin
    endif
    if (tparop(ifac).gt.tmax) then
      n1max = n1max +1
      tparop(ifac) = tmax
    endif

!       Max-Min
    if (tpmax.le.tparop(ifac)) then
      ifacmx = ifac
      tpmax = tparop(ifac)
      qcmax = qconv
      qrmax = qrayt
    endif

    if (tpmin.ge.tparop(ifac)) then
      ifacmn = ifac
      tpmin = tparop(ifac)
      qcmin = qconv
      qrmin = qrayt
    endif

    tzomax(izone) = max(tzomax(izone),tparop(ifac))
    tzomin(izone) = min(tzomin(izone),tparop(ifac))

!       Reperage pour impression
    iipgrn = 1
    indtp(izone) = ipgrno

  endif
enddo


!     3.3) parois reflechissantes
!     ---------------------------


do ifac = 1, nfabor
  if (isothp(ifac).eq.iprefl) then
    izone = izfrap(ifac)

!       Calcul
    esl    = epap(ifac) /xlamp(ifac)
    qconv  = flconp(ifac)

    detep  = ( esl*qconv - (tparop(ifac)-textp(ifac)) )           &
           / ( 1.d0 + esl*hfconp(ifac) )

    rapp = detep /tparop (ifac)
    abrapp = abs(rapp)

!       Relaxation
    if (abrapp.ge.tx) then
      nrelax = nrelax +1
      tparop(ifac) = tparop(ifac) * (1.d0+ tx*sign(und0,rapp))
    else
      tparop(ifac) = tparop(ifac) + detep
    endif

    rapmax = max(rapmax,abrapp)

    if (rapp.le.0.d0) then
      nmoins = nmoins+1
    else
      nplus  = nplus +1
    endif

!       Clipping
    if (tparop(ifac).lt.tmin) then
      n1min = n1min +1
      tparop(ifac) = tmin
    endif
    if (tparop(ifac).gt.tmax) then
      n1max = n1max +1
      tparop(ifac) = tmax
    endif

!       Max-Min
    if (tpmax.le.tparop(ifac)) then
      ifacmx = ifac
      tpmax = tparop(ifac)
      qcmax = qconv
      qrmax = 0.d0
    endif

    if (tpmin.ge.tparop(ifac)) then
      ifacmn = ifac
      tpmin = tparop(ifac)
      qcmin = qconv
      qrmin = 0.d0
    endif

    tzomax(izone) = max(tzomax(izone),tparop(ifac))
    tzomin(izone) = min(tzomin(izone),tparop(ifac))

!       Reperage pour impression
    iipref = 1
    indtp(izone) = iprefl

  endif
enddo



!     3.4) parois a flux de conduction impose non reflechissante
!          si le flux impose est nul, la paroi est adiabatique
!          (la transmition de chaleur par rayonnement est
!          equilibree avec celle de la convection)


do ifac = 1, nfabor
  if (isothp(ifac).eq.ifgrno) then
    izone = izfrap(ifac)

!       Calcul
    qconv  = flconp(ifac)
    qinci  = qincip(ifac)
    epp    = epsp(ifac)

    sigt3  = stephn*tparop(ifac)**3
    qrayt  = epp *( qinci -sigt3*tparop(ifac))

    detep  = ( qconv + qrayt - rcodcl(ifac,ivart,3) )             &
           / ( 4.d0*epp*sigt3 + hfconp(ifac) )

    rapp   = detep /tparop (ifac)
    abrapp = abs(rapp)

!       Relaxation
    if (abrapp.ge.tx) then
      nrelax = nrelax +1
      tparop(ifac) = tparop(ifac) * (1.d0+ tx*sign(und0,rapp))
    else
      tparop(ifac) = tparop(ifac) + detep
    endif

    rapmax = max(rapmax,abrapp)

    if (rapp.le.0.d0) then
      nmoins = nmoins+1
    else
      nplus  = nplus +1
    endif

!       Clipping
    if (tparop(ifac).lt.tmin) then
      n1min = n1min +1
      tparop(ifac) = tmin
    endif
    if (tparop(ifac).gt.tmax) then
      n1max = n1max +1
      tparop(ifac) = tmax
    endif

!       Max-Min
    if (tpmax.le.tparop(ifac)) then
      ifacmx = ifac
      tpmax = tparop(ifac)
      qcmax = qconv
      qrmax = qrayt
    endif

    if (tpmin.ge.tparop(ifac)) then
      ifacmn = ifac
      tpmin = tparop(ifac)
      qcmin = qconv
      qrmin = qrayt
    endif

    tzomax(izone) = max(tzomax(izone),tparop(ifac))
    tzomin(izone) = min(tzomin(izone),tparop(ifac))

!       Reperage pour impression
    iifgrn = 1
    indtp(izone) = ifgrno

  endif
enddo


!     3.5) parois a flux de conduction imposee et reflechissante
!          equivalent a impose un flux total au fluide



do ifac = 1, nfabor
  if (isothp(ifac).eq.ifrefl) then
    izone = izfrap(ifac)

!       Calcul
    iel = ifabor(ifac)

    tparop(ifac)= hfconp(ifac)*tempkp(iel)-rcodcl(ifac,ivart,3)
    tparop(ifac)= tparop(ifac)/max(hfconp(ifac),epzero)

!       Clipping
    if (tparop(ifac).lt.tmin) then
      n1min = n1min +1
      tparop(ifac) = tmin
    endif
    if (tparop(ifac).gt.tmax) then
      n1max = n1max +1
      tparop(ifac) = tmax
    endif

!       Max-Min
    qconv  = flconp(ifac)
    qrayt  = 0.d0

    if (tpmax.le.tparop(ifac)) then
      ifacmx = ifac
      tpmax = tparop(ifac)
      qcmax = qconv
      qrmax = qrayt
    endif

    if (tpmin.ge.tparop(ifac)) then
      ifacmn = ifac
      tpmin = tparop(ifac)
      qcmin = qconv
      qrmin = qrayt
    endif

    tzomax(izone) = max(tzomax(izone),tparop(ifac))
    tzomin(izone) = min(tzomin(izone),tparop(ifac))

!       Reperage pour impression
    iifref = 1
    indtp(izone) = ifrefl

  endif
enddo


!===============================================================================
! 4. IMPRESSIONS
!===============================================================================



! Test pour savoir s'il y a des zones de paroi
!  (donc si des impressions sont necessaires)

if (irangp.ge.0) then
  call parimx(nozarm,indtp )
  !==========
endif

indtpm = 0
do izone = 1, nozrdm
  if (indtp(izone).ne.0) then
    indtpm = 1
  endif
enddo

! Si il y a des parois

if(indtpm.gt.0) then

!   Si on veut imprimer en niveau 1 au moins

  if (iimpar.ge.1) then

!      Calcul de TZOMOY FLUNET RADIOS SURFT

    do izone = 1, nozrdm
      tzomoy(izone) = zero
      flunet(izone) = zero
      radios(izone) = zero
      surft (izone) = zero
    enddo
    do ifac = 1, nfabor
      srfbn = surfbn(ifac)
      izone = izfrap(ifac)
      if (indtp(izone).ne.0) then
        tp4 = tparop(ifac)**4
        tzomoy(izone)= tzomoy(izone) + tparop(ifac)*srfbn
        flunet(izone)= flunet(izone)                              &
             + epsp(ifac) *(qincip(ifac) -stephn*tp4 )*srfbn
        radios(izone)= radios(izone)                              &
             - ( epsp(ifac)       *stephn*tp4                     &
             +(1.d0-epsp(ifac))*qincip(ifac)      )*srfbn
        surft (izone) = surft(izone) + srfbn
      endif
    enddo

    if(irangp.ge.0) then
      call parrsm(nozarm,tzomoy)
      !==========
      call parrsm(nozarm,flunet)
      !==========
      call parrsm(nozarm,radios)
      !==========
      call parrsm(nozarm,surft )
      !==========
    endif

    do izone = 1, nozrdm
      if (indtp(izone).ne.0) then
        tzomoy(izone) = tzomoy(izone)/surft(izone)
        radios(izone) = radios(izone)/surft(izone)
      endif
    enddo


!      Determination de la TPMAX TPMIN et des grandeurs associees

    if(ifacmx.gt.0) then
      xtpmax = xyzcen(1,ifabor(ifacmx))
      ytpmax = xyzcen(2,ifabor(ifacmx))
      ztpmax = xyzcen(3,ifabor(ifacmx))
    else
      xtpmax = 0.d0
      ytpmax = 0.d0
      ztpmax = 0.d0
    endif
    if(ifacmn.gt.0) then
      xtpmin = xyzcen(1,ifabor(ifacmn))
      ytpmin = xyzcen(2,ifabor(ifacmn))
      ztpmin = xyzcen(3,ifabor(ifacmn))
    else
      xtpmin = 0.d0
      ytpmin = 0.d0
      ztpmin = 0.d0
    endif

    if(irangp.ge.0) then

      rdptmp(1) = xtpmax
      rdptmp(2) = ytpmax
      rdptmp(3) = ztpmax
      rdptmp(4) = qcmax
      rdptmp(5) = qrmax
      nbrval = nbrrdp
      call parmxl(nbrval,tpmax,rdptmp)
      !==========
      xtpmax =  rdptmp(1)
      ytpmax =  rdptmp(2)
      ztpmax =  rdptmp(3)
      qcmax  =  rdptmp(4)
      qrmax  =  rdptmp(5)

      rdptmp(1) = xtpmin
      rdptmp(2) = ytpmin
      rdptmp(3) = ztpmin
      rdptmp(4) = qcmin
      rdptmp(5) = qrmin
      nbrval = nbrrdp
      call parmnl(nbrval,tpmin,rdptmp)
      !==========
      xtpmin =  rdptmp(1)
      ytpmin =  rdptmp(2)
      ztpmin =  rdptmp(3)
      qcmin  =  rdptmp(4)
      qrmin  =  rdptmp(5)

    endif


!      Determination des compteurs et autres

    if(irangp.ge.0) then

      call parmax(rapmax)
      !==========

      inttm2(1) = nmoins
      inttm2(2) = nplus
      inttm2(3) = n1min
      inttm2(4) = n1max
      inttm2(5) = nrelax
      nbrval = nb2int
      call parism(nbrval,inttm2)
      !==========
      nmoins = inttm2(1)
      nplus  = inttm2(2)
      n1min  = inttm2(3)
      n1max  = inttm2(4)
      nrelax = inttm2(5)

      call parrmx(nozarm,tzomax)
      !==========
      call parrmn(nozarm,tzomin)
      !==========

      inttm1(1) = iitpim
      inttm1(2) = iipgrn
      inttm1(3) = iipref
      inttm1(4) = iifgrn
      inttm1(5) = iifref
      nbrval    = nb1int
      call parimx(nbrval,inttm1)
      !==========
      iitpim = inttm1(1)
      iipgrn = inttm1(2)
      iipref = inttm1(3)
      iifgrn = inttm1(4)
      iifref = inttm1(5)

    endif


!      Impressions

    write(nfecra,1000)
    write(nfecra,1010)

    if (nrelax.gt.0) then
      write(nfecra,2010) tx*100.d0 ,nrelax
      write(nfecra,1010)
    endif

    if ( n1min.gt.0 .or. n1max.gt.0 ) then
      write(nfecra,2020)
      write(nfecra,2030) n1min
      write(nfecra,2040) n1max
      write(nfecra,1010)
    endif

    if (rapmax.gt.0 .or. nmoins.gt.0 .or. nplus.gt.0 ) then
      write(nfecra,2050) rapmax * 100.d0
      write(nfecra,2060) nmoins
      write(nfecra,2070) nplus
      write(nfecra,1010)
    endif

    if (iitpim.eq.1) then
      write(nfecra,3020)
      do izone = 1, nozrdm
        if (indtp(izone).eq.itpimp) then
          write(nfecra,3000) izone,                               &
               tzomax(izone)-tkelvi,                              &
               tzomin(izone)-tkelvi,                              &
               tzomoy(izone)-tkelvi,                              &
               flunet(izone)
        endif
      enddo
      write(nfecra,1010)
    endif

    if (iipgrn.eq.1) then
      write(nfecra,3010)
      do izone = 1, nozrdm
        if (indtp(izone).eq.ipgrno) then
          write(nfecra,3000) izone,                               &
               tzomax(izone)-tkelvi,                              &
               tzomin(izone)-tkelvi,                              &
               tzomoy(izone)-tkelvi,                              &
               flunet(izone)
        endif
      enddo
      write(nfecra,1010)
    endif

    if (iipref.eq.1) then
      write(nfecra,3030)
      do izone = 1, nozrdm
        if (indtp(izone).eq.iprefl) then
          write(nfecra,3000) izone,                               &
               tzomax(izone)-tkelvi,                              &
               tzomin(izone)-tkelvi,                              &
               tzomoy(izone)-tkelvi,                              &
               flunet(izone)
        endif
      enddo
      write(nfecra,1010)
    endif

    if (iifgrn.eq.1) then
      write(nfecra,3040)
      do izone = 1, nozrdm
        if (indtp(izone).eq.ifgrno) then
          write(nfecra,3000) izone,                               &
               tzomax(izone)-tkelvi,                              &
               tzomin(izone)-tkelvi,                              &
               tzomoy(izone)-tkelvi,                              &
               flunet(izone)
        endif
      enddo
      write(nfecra,1010)
    endif

    if (iifref.eq.1) then
      write(nfecra,3050)
      do izone = 1, nozrdm
        if (indtp(izone).eq.ifrefl) then
          write(nfecra,3000) izone,                               &
               tzomax(izone)-tkelvi,                              &
               tzomin(izone)-tkelvi,                              &
               tzomoy(izone)-tkelvi,                              &
               flunet(izone)
        endif
      enddo
      write(nfecra,1010)
    endif

  endif


!   Si on veut imprimer EN PLUS en niveau 2
!                       =======

  if (iimpar.ge.2) then

    write(nfecra,5010) tpmax-tkelvi
    write(nfecra,5020) xtpmax, ytpmax, ztpmax
    write(nfecra,5030) qcmax
    write(nfecra,5040) qrmax
    write(nfecra,5050) tpmin-tkelvi
    write(nfecra,5060) xtpmin, ytpmin, ztpmin
    write(nfecra,5070) qcmin
    write(nfecra,5080) qrmin

  endif

endif

! -------
! FORMATS
! -------

 1000 FORMAT (/, 3X,'** INFORMATIONS SUR LA TEMPERATURE DES PAROIS',/,  &
           3X,'   ------------------------------------------')

 1010 format('------------------------------------'                     &
         ,'-----------------------------------')

 2010 format('ATTENTION, temperature de paroi relaxee a ',G7.2,'%'      &
      ,' en (',I8,' points)')

 2020 format('ATTENTION, temperature de paroi CLIPPE AU MIN-MAX :')

 2030 format('NOMBRE DE POINTS CLIPPE AU MINIMUM : ',I8)

 2040 format('NOMBRE DE POINTS CLIPPE AU MAXIMUM : ',I8)

 2050 format('Variation maximale : ',G9.4,'%')

 2060 format('Temperature de paroi diminuant  : ',I8,' faces de bord')

 2070 format('Temperature de paroi augmentant : ',I8,' faces de bord')

 3000 format(i10,8x,e10.4,4x,e10.4,4x,e10.4,4x,e10.4)

 3010 format('Grises ou noires Temp max (C)  '                          &
      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')

 3020 format('Profils imposes  Temp max (C)  '                          &
      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')

 3030 format('Parois a EPS=0   Temp max (C)  '                          &
      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')

 3040 format('Flux imp EPS!=0  Temp max (C)  '                          &
      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')

 3050 format('Flux imp EPS=0   Temp max (C)  '                          &
      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')

 5010 format(/,12X,'TEMPERATURE PAROI MAXIMALE (Degre celsius) = ',     &
         g15.7)
 5020 format(15X,' AU POINT X Y Z =  ',E10.4,2X,E10.4,2X,E10.4)

 5030 format(15X,' FLUX CONVECTIF  = ',G15.7)

 5040 format(15X,' FLUX RADIATIF  = ',G15.7)

 5050 format(/,12X,'TEMPERATURE PAROI MINIMALE (Degre celsius) = ',     &
         g15.7)

 5060 format(15X,' AU POINT X Y Z =  ',E10.4,2X,E10.4,2X,E10.4)

 5070 format(15X,' FLUX CONVECTIF  = ',G15.7)

 5080 format(15X,' FLUX RADIATIF  = ',G15.7,/)


! ---
! FIN
! ---

return

end subroutine

!   ---------------------------------------------------------------------------------
!   Variation maximale : .1484E-13%
!   Temperature de paroi diminuant  :       33
!   Temperature de paroi augmentant :       27
!   ---------------------------------------------------------------------------------
!   Zones de parois  Temp max (C)    Temp min (C)    Temp moy (C)    Flux Net (W)
!           15        0.7941E+03      0.7410E+03      0.7688E+03      0.7688E+03
!   ---------------------------------------------------------------------------------
!   Profils imposes  Temp max (C)    Temp min (C)    Temp moy (C)    Flux Net (W)
!           15        0.7941E+03      0.7410E+03      0.7688E+03      0.7688E+03
!   ---------------------------------------------------------------------------------

