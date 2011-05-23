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

subroutine raydom &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   izfrad ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  Enveloppe principale du module de resolution de l'equation
!  des transferts radiatifs

!  Deux methodes sont disponibles :

!    1) La methode : "Discretes Ordinates Methods" (DOM)
!    2) L'approximation P-1 (recommande uniquement pour le CP)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfrad(nfabor    ! te ! <-- ! numero de zone des faces de bord               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use entsor
use optcal
use cstphy
use cstnum
use pointe
use parall
use period
use lagpar
use lagdim
use lagran
use ppppar
use ppthch
use fuincl
use ppincl
use cpincl
use radiat
use ihmpre
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          itypfb(nfabor)
integer          izfrad(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision ra(*)

! Local variables

integer          idebia , idebra
integer          iappel
integer          ifac   , iel    , iok    , izone
integer          inc    , iccocg , iwarnp , imligp , nswrgp
integer          mode   , icla   , ipcla  , ivar0
integer          iscat  , ivart
integer          iflux(nozrdm)
double precision epsrgp, climgp, extrap
double precision aa, bb, ckmin, unspi, xlimit, cofrmn, flunmn
double precision flux(nozrdm)
double precision vv, sf, xlc, xkmin, pp

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:,:) :: xam
double precision, allocatable, dimension(:) :: drtp, smbrs, rovsdt
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: w10
double precision, allocatable, dimension(:,:) :: tempk
double precision, allocatable, dimension(:) :: cofrua, cofrub
double precision, allocatable, dimension(:) :: flurds, flurdb


integer    ipadom
data       ipadom /0/
save       ipadom

!==============================================================================
!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

! Allocate temporary arrays for the radiative equations resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(dam(ncelet), xam(nfac,2))
allocate(drtp(ncelet), smbrs(ncelet), rovsdt(ncelet))

! Allocate specific arrays for the radiative transfert module
allocate(tempk(ncelet,nrphas))
allocate(cofrua(nfabor), cofrub(nfabor))
allocate(flurds(nfac), flurdb(nfabor))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet), w6(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))
allocate(w10(ncelet))

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INITIALISATIONS GENERALES
!===============================================================================

!---> NUMERO DE PASSAGE RELATIF

ipadom = ipadom + 1
if (ipadom.gt.1 .and. mod(ntcabs,nfreqr).ne.0) return

write(nfecra,1000)

!---> INITIALISATION DES CONSTANTES

unspi = 1.d0/pi

!===============================================================================
! 2. PHASE
!===============================================================================

!---> NUMERO DU SCALAIRE ET DE LA VARIABLE THERMIQUE

  iscat = iscalt
  ivart = isca(iscalt)

!===============================================================================
! 3.1 COEFFICIENT D'ABSORPTION DU MILIEU SEMI-TRANSPARENT
!===============================================================================

!--> INITIALISATION NON ADMISSIBLE POUR TEST APRES USRAY3

  do iel = 1,ncel
    propce(iel,ipproc(icak(1))) = -grand
  enddo

!--> COEFFICIENT D'ABSORPTION POUR LA PHYSIQUE PARTICULIERE

!  ATTENTION  : DANS LE CAS DE L'APPROXIMATION P-1, LE COEFFICIENT
!    D'ABSORPTION EST UTILISE POUR CALCULER LES CONDITIONS AUX
!    LIMITES DE L'EQUATION A RESOUDRE.

  if (ippmod(iphpar).ge.2) then

    call ppcabs                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                             &
   itypfb ,                                                      &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     ,                                     &
   ra    )

!-----> W10 SERT A STOCKER TEMPORAIREMENT
!       LE COEFFICIENT D'ABSORPTION DU MELANGE GAZ-PARTICULE

    if (ippmod(icp3pl).ge.0 .or. ippmod(icfuel).ge.0) then

      do iel = 1,ncel
        w10(iel) = propce(iel,ipproc(icak(1)))
      enddo

      if (ippmod(icp3pl).ge.0 ) then
        do icla = 1,nclacp
          ipcla = 1+icla
          do iel = 1,ncel
            w10(iel) = w10(iel)                                   &
                   + ( propce(iel,ipproc(ix2(icla)))              &
                   *   propce(iel,ipproc(icak(ipcla))) )
          enddo
        enddo
      else if ( ippmod(icfuel) .ge.0 ) then
        do icla = 1,nclafu
          ipcla = 1+icla
          do iel = 1,ncel
            w10(iel) = w10(iel)                                   &
                     + ( rtpa(iel,isca(iyfol(icla)))              &
                       * propce(iel,ipproc(icak(ipcla))) )
          enddo
        enddo
      endif

      do iel = 1,ncel
        propce(iel,ipproc(icak(1))) = w10(iel)
      enddo
    endif

  else


!---> LECTURES DES DONNEES UTILISATEURS

!       - Interface Code_Saturne
!         ======================

    if (iihmpr.eq.1) then

      call uiray3(propce(1,ipproc(icak(1))), ncel, imodak)
      !==========

      if (iirayo.eq.2 .and. ippmod(iphpar).le.1 .and. ipadom.le.3) then
        sf = 0.d0
        vv = 0.d0

!         Calcul de la longueur caracteristique du domaine de calcul
        do ifac = 1,nfabor
          sf = sf + sqrt(surfbo(1,ifac)**2 +                      &
                         surfbo(2,ifac)**2 +                      &
                         surfbo(3,ifac)**2 )
        enddo
        if (irangp.ge.0) then
          call parsom(sf)
          !==========
        endif

        do iel = 1,ncel
          vv = vv + volume(iel)
        enddo
        if (irangp.ge.0) then
          call parsom(vv)
          !==========
        endif

        xlc = 3.6d0 * vv / sf

!             Clipping pour la variable CK

        xkmin = 1.d0 / xlc

        iok = 0
        do iel = 1,ncel
          if (propce(iel,ipproc(icak(1))).lt.xkmin) then
            iok = iok +1
          endif
        enddo

!     Alerte si epaisseur optique trop grande

        pp = xnp1mx/100.0d0
        if (dble(iok).gt.pp*dble(ncel)) then
          write(nfecra,6000) xkmin, dble(iok)/dble(ncel)*100.d0,  &
                             xnp1mx
        endif
      endif

    endif

    call usray3                                                   &
    !==========
 ( nvar   , nscal  , iappel ,                                     &
   itypfb ,                                                      &
   izfrad ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   propce(1,ipproc(icak(1))),                                &
   ra     )

  endif

!--> VERIFICATIONS GENERALES :

!---> P-1 : VERIFICATION QUE CK DU MILIEU EST STRICTEMENT
!       SUPERIEUR A ZERO POUR TOUTES CELLULES

  if (iirayo.eq.2) then

    ckmin = propce(1,ipproc(icak(1)))
    do iel = 1, ncel
      ckmin = min(ckmin,propce(iel,ipproc(icak(1))))
    enddo
    if (ckmin.lt.0.d0) then
      write(nfecra,2020)
      call csexit (1)
      !==========
    endif

  else if (iirayo.eq.1) then

!---> DOM : VERIFICATION QUE CK DU MILIEU EST SUPERIEUR OU EGAL A 0

    ckmin = propce(1,ipproc(icak(1)))
    do iel = 1, ncel
      ckmin = min(ckmin,propce(iel,ipproc(icak(1))))
    enddo
    if (ckmin.lt.0.d0) then
      write(nfecra,2010) ckmin
      call csexit (1)
      !==========
    endif

  endif

!---> VERIFICATION D'UN CAS TRANSPARENT

  aa = zero
  do iel = 1,ncel
    aa = aa + propce(iel,ipproc(icak(1)))
  enddo
  if (irangp.ge.0) then
    call parmax(aa)
    !==========
  endif
  if (aa.le.epzero) then
    write(nfecra,1100)
    idiver = -1
  endif

!===============================================================================
! 3.2 CONDITIONS AUX LIMITES POUR LES EQNS DE LA DOM ET DE L'APPROX P-1

!     REMPLISSAGE DES CONDITIONS AUX LIMITES POUR LA LUMINANCE
!     (TABLEAUX  COFRUA & COFRUB)
!===============================================================================

!-----> INITIALISATIONS NON ADMISSIBLES POUR TEST APRES USRAY5

  do ifac = 1,nfabor
    cofrua(ifac) = -grand
    cofrub(ifac) = -grand
  enddo

  iappel = 1

  call usray5                                                     &
  !==========
 ( nvar   , nscal  , iappel ,                                     &
   itypfb ,                                                      &
   izfrad ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cofrua , cofrub ,                                              &
   propfb(1,ipprob(itparo)) , propfb(1,ipprob(iqinci)) ,          &
   propfb(1,ipprob(ifnet))  , propfb(1,ipprob(ixlam))  ,          &
   propfb(1,ipprob(iepa))   , propfb(1,ipprob(ieps))   ,          &
   propce(1,ipproc(icak(1)))    ,                            &
   ra     )


!---> BLINDAGE POUR UN DIRICHLET SUR LA LUMINANCE
!      (UNIQUEMENT POUR LA METHODE DOM, SANS OBJET POUR L'APPROX P-1)

  if (iirayo.eq.1) then
    do ifac = 1,nfabor
      cofrub(ifac) = zero
    enddo
  endif

!---> VERIFICATION DU REMPLISSAGE DE COFRUA & COFRUB

!     Attention : dans le cas de l'approx P-1 la valeur de COFRUA peut
!     etre grande (de l'ordre de TPAROI**4), d'ou la valeur de COFRMN

  iok = 0
  xlimit = -grand*0.1d0
! GRAND n'est pas assez grand...
  cofrmn = rinfin

  do ifac = 1,nfabor
    if (cofrua(ifac).le.xlimit) then
      iok = iok + 1
      cofrmn = min(cofrmn,cofrua(ifac))
      write(nfecra,3000)ifac,izfrad(ifac),itypfb(ifac)
    endif
  enddo

  if (iok.ne.0) then
    write(nfecra,3100) cofrmn
    call csexit (1)
    !==========
  endif

  cofrmn = rinfin

  if (iirayo.eq.2) then

    do ifac = 1,nfabor
      if (cofrub(ifac).le.xlimit) then
        iok = iok + 1
        cofrmn = min(cofrmn,cofrub(ifac))
        write(nfecra,3000)ifac,izfrad(ifac),itypfb(ifac)
      endif
    enddo

    if (iok.ne.0) then
      write(nfecra,3200) cofrmn
      call csexit (1)
      !==========
    endif

  endif

!===============================================================================
! 4. STOCKAGE DE LA TEMPERATURE (en Kelvin) dans TEMPK(IEL,IPH)
!===============================================================================

  if (idiver.ge.0) then

    if(abs(iscsth(iscat)).eq.1) then

!---> TRANSPORT DE LA TEMPERATURE

      if (iscsth(iscat).eq.-1) then
        do iel = 1, ncel
          tempk(iel,1) = rtpa(iel,ivart) + tkelvi
        enddo
      else
        do iel = 1, ncel
          tempk(iel,1) = rtpa(iel,ivart)
        enddo
      endif

!---> TRANSPORT DE L'ENTHALPIE (FLURDB est un auxiliaire)

    else if (iscsth(iscat).eq.2) then

      mode = 1

      if (ippmod(iphpar).le.1) then

        call usray4                                               &
        !==========
 ( nvar   , nscal  ,                                             &
   mode   ,                                                       &
   itypfb ,                                                      &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   propfb(1,ipprob(itparo)) , flurdb , tempk(1,1)  ,              &
   ra     )

      else

        call ppray4                                               &
        !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                             &
   mode   ,                                                       &
   itypfb ,                                                      &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   propfb(1,ipprob(itparo)) , flurdb , tempk(1,1)  ,              &
   ra     )

      endif

      if ( ippmod(icp3pl).ge.0 ) then

! Temperature des particules

        do icla = 1, nclacp
          ipcla = 1+icla
          do iel = 1,ncel
            tempk(iel,ipcla) = propce(iel,ipproc(itemp2(icla)))
          enddo
        enddo

!         Fuel

      else if ( ippmod(icfuel).ge.0 ) then

        do icla = 1, nclafu
          ipcla = 1+icla
          do iel = 1,ncel
            tempk(iel,ipcla) = propce(iel,ipproc(itemp3(icla)))
          enddo
        enddo

      endif

    else
      write(nfecra,3500)iscat,iscsth(iscat)
      call csexit (1)
      !==========
    endif

!---> ON SE SERT DE PROPCE(IEL,IPPROC(ITSRI(1))) COMME UN AUXILIAIRE POUR
!       STOCKER STEPHN*CK*TEMPK**4 ICI
!       PLUS BAS ON JUSTIFIERA LE NOM.

    if ( ippmod(icod3p).eq.-1 .and.                               &
         ippmod(icoebu).eq.-1       ) then

!           Rayonnement standard, flamme CP ou fuel

      do iel = 1,ncel
        propce(iel,ipproc(itsri(1))) = stephn  *             &
         propce(iel,ipproc(icak(1)))*(tempk(iel,1)**4)
      enddo

    else

!           Flamme de diffusion ou flamme de premelange

      do iel = 1,ncel
        propce(iel,ipproc(itsri(1))) = stephn  *             &
         propce(iel,ipproc(icak(1)))*propce(iel,ipproc(it4m))
      enddo

    endif

!     Charbon

    if ( ippmod(icp3pl).ge.0 ) then
      do icla = 1,nclacp
        ipcla = 1+icla
        do iel = 1,ncel
          propce(iel,ipproc(itsri(ipcla))) =  stephn  *           &
            propce(iel,ipproc(icak(ipcla)))*(tempk(iel,ipcla)**4)
        enddo
      enddo

!       Fuel

    else if ( ippmod(icfuel).ge.0 ) then
      do icla = 1,nclafu
        ipcla = 1+icla
        do iel = 1,ncel
          propce(iel,ipproc(itsri(ipcla))) =  stephn  *           &
            propce(iel,ipproc(icak(ipcla)))*(tempk(iel,ipcla)**4)
        enddo
      enddo
    endif

  else
    do iel = 1,ncel
      propce(iel,ipproc(itsri(1))) = zero
    enddo
! fin de IF (IDIVER.GE.0) THEN
  endif

!===============================================================================
! 5.1 MODELE DE RAYONNEMENT P-1
!===============================================================================

  if (iirayo.eq.2) then

!--> Terme source explicite de l'equation sur Theta4

    do iel = 1, ncel
      smbrs(iel) = 3.d0 * propce(iel,ipproc(icak(1))) *      &
         ( tempk(iel,1) ** 4) * volume(iel)
    enddo

! Tenir compte de l'absorption des particules

!       Charbon

    if ( ippmod(icp3pl).ge.0 ) then
      do icla = 1,nclacp
        ipcla = 1+icla
        do iel = 1,ncel
          smbrs(iel) = smbrs(iel)                                 &
                         + ( 3.d0 * propce(iel,ipproc(ix2(icla))) &
                         *   propce(iel,ipproc(icak(ipcla)))      &
                         *  (tempk(iel,ipcla)**4) * volume(iel) )
        enddo
      enddo

!       FUEL

    else if ( ippmod(icfuel).ge.0 ) then
      do icla = 1,nclafu
        ipcla = 1+icla
        do iel = 1,ncel
          smbrs(iel) = smbrs(iel)                                 &
                      +( 3.d0 * rtpa(iel,isca(iyfol(icla)))       &
                         *   propce(iel,ipproc(icak(ipcla)))      &
                         *  (tempk(iel,ipcla)**4) * volume(iel) )
        enddo
      enddo
    endif

!--> Terme source implicite de l'equation sur Theta4

    do iel = 1, ncel
      rovsdt(iel) =  3.d0 * propce(iel,ipproc(icak(1)))      &
                        * volume(iel)
    enddo

! Tenir compte de l'absorption des particules

!        Charbon

    if ( ippmod(icp3pl).ge.0 ) then
      do icla = 1,nclacp
        ipcla = 1+icla
        do iel = 1,ncel
          rovsdt(iel) = rovsdt(iel)                               &
                          + (3.d0 * propce(iel,ipproc(ix2(icla))) &
                          *  propce(iel,ipproc(icak(ipcla)))      &
                          * volume(iel) )
        enddo
      enddo

!        Fuel

    else if ( ippmod(icfuel).ge.0 ) then
      do icla = 1,nclafu
        ipcla = 1+icla
        do iel = 1,ncel
          rovsdt(iel) = rovsdt(iel)                               &
                       + (3.d0*rtpa(iel,isca(iyfol(icla)))        &
                          *  propce(iel,ipproc(icak(ipcla)))      &
                          * volume(iel) )
        enddo
      enddo

    endif

!--> Inverse du coefficient de diffusion de l'equation sur Theta4
!       A priori W10 contient deja la bonne info, mais pour plus de
!       securite  on le re-remplit

    do iel = 1, ncel
      w10(iel) =  propce(iel,ipproc(icak(1)))
    enddo

! Tenir compte de l'absorption des particules

!        Charbon

    if ( ippmod(icp3pl).ge.0 ) then
      do icla = 1,nclacp
        ipcla = 1+icla
        do iel = 1,ncel
          w10(iel) = w10(iel)                                     &
                    + ( propce(iel,ipproc(ix2(icla)))             &
                      * propce(iel,ipproc(icak(ipcla)))   )
        enddo
      enddo

!        Fuel

    else if ( ippmod(icfuel).ge.0 ) then
      do icla = 1,nclafu
        ipcla = 1+icla
        do iel = 1,ncel
          w10(iel) = w10(iel)                                     &
                    + ( rtpa(iel,isca(iyfol(icla)))               &
                      * propce(iel,ipproc(icak(ipcla)))   )
        enddo
      enddo

    endif

    call raypun                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                             &
   itypfb ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cofrua , cofrub ,                                              &
   flurds , flurdb ,                                              &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &

   propce(1,ipproc(iabs(1))),propce(1,ipproc(iemi(1))), &
   propce(1,ipproc(itsre(1))) , propce(1,ipproc(iqx))  ,     &
   propce(1,ipproc(iqy))   , propce(1,ipproc(iqz))  ,             &
   propfb(1,ipprob(iqinci)), propfb(1,ipprob(ieps)) ,             &
   propfb(1,ipprob(itparo)),                                      &

   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , w10    ,                   &
   ra     )

!===============================================================================
! 5.2 RESOLUTION DE L'EQUATION DES TRANSFERTS RADIATIFS
!===============================================================================

  else if (iirayo.eq.1) then

!--> Terme source explicite de l'equation sur la luminance


    do iel = 1, ncel
      smbrs(iel) = propce(iel,ipproc(itsri(1))) *volume(iel) &
                                                    *unspi
    enddo

!       Charbon

    if ( ippmod(icp3pl).ge.0 ) then
      do icla = 1,nclacp
        ipcla = 1+icla
        do iel = 1,ncel
          smbrs(iel) = smbrs(iel)                                 &
                  + propce(iel,ipproc(ix2(icla)))                 &
                   *propce(iel,ipproc(itsri(ipcla)))*volume(iel)  &
                   *unspi
        enddo
      enddo

!       Fuel

    elseif ( ippmod(icfuel).ge.0 ) then
      do icla = 1,nclafu
        ipcla = 1+icla
        do iel = 1,ncel
          smbrs(iel) = smbrs(iel)                                 &
                      + rtpa(iel,isca(iyfol(icla)))               &
                   *propce(iel,ipproc(itsri(ipcla)))*volume(iel)  &
                   *unspi
        enddo
      enddo

    endif

!--> Terme source implicite de l'equation sur la luminance
!      KL + div(LS) = KL0 integre sur le volume de controle

    do iel = 1, ncel
      rovsdt(iel) = propce(iel,ipproc(icak(1))) * volume(iel)
    enddo

!        Charbon

    if ( ippmod(icp3pl).ge.0 ) then
      do icla = 1,nclacp
        ipcla = 1+icla
        do iel = 1,ncel
          rovsdt(iel) = rovsdt(iel) +                             &
                          propce(iel,ipproc(ix2(icla)))           &
                          * propce(iel,ipproc(icak(ipcla)))       &
                          * volume(iel)
        enddo
      enddo

!        Fuel

    elseif ( ippmod(icfuel).ge.0 ) then
      do icla = 1,nclafu
        ipcla = 1+icla
        do iel = 1,ncel
          rovsdt(iel) = rovsdt(iel)                               &
                       + rtpa(iel,isca(iyfol(icla)))              &
                          * propce(iel,ipproc(icak(ipcla)))       &
                          * volume(iel)
        enddo
      enddo

    endif

    call raysol                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                             &
   itypfb ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cofrua , cofrub ,                                              &
   flurds , flurdb ,                                              &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , w10    ,                   &
   propce(1,ipproc(iabs(1))),propce(1,ipproc(iemi(1))) ,&
   propce(1,ipproc(itsre(1))),propce(1,ipproc(iqx))         ,&
   propce(1,ipproc(iqy))    , propce(1,ipproc(iqz))   ,           &
   propfb(1,ipprob(iqinci) ), propfb(1,ipprob(ifnet)) ,           &
   ra     )


  endif

!===============================================================================
! 5.3 STOCKAGE DE L'INTEGRALE DE LA LUMINANCE POUR LE LANGRANGIEN
!===============================================================================

!  Si dans le module lagrangien on resout une equation de la temperature
!    sur les particules (IPHYLA=1 et ITPVAR=1) ou si les particules
!    sont des grains de charbon (IPHYLA=2), on a besoin de
!                                     /    ->  ->
!    l'integrale de la luminance SA= /  L( X , S ). DOMEGA
!                                   /4.PI
!  On stocke cette variable quelque soit le choix des options

  do iel = 1,ncel
    propce(iel,ipproc(ilumin)) = propce(iel,ipproc(itsre(1)))
  enddo

!===============================================================================
! 6. FLUX NET RADIATIF AUX PAROIS : CALCUL ET INTEGRATION
!===============================================================================


!---> INITIALISATION NON ADMISSIBLE POUR TEST APRES USRAY5

  do ifac = 1,nfabor
    propfb(ifac,ipprob(ifnet)) = -grand
  enddo

!---> LECTURES DES DONNEES UTILISATEURS

  iappel = 2

  call usray5                                                     &
  !==========
 ( nvar   , nscal  , iappel ,                                     &
   itypfb ,                                                      &
   izfrad ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cofrua , cofrub ,                                              &
   propfb(1,ipprob(itparo)) , propfb(1,ipprob(iqinci)) ,          &
   propfb(1,ipprob(ifnet))  , propfb(1,ipprob(ixlam))  ,          &
   propfb(1,ipprob(iepa))   , propfb(1,ipprob(ieps))   ,          &
   propce(1,ipproc(icak(1)))  ,                              &
   ra     )

!---> VERIFICATION DE FLUNET

  iok = 0
  xlimit = -grand*0.1d0
  flunmn = grand

  do ifac = 1,nfabor
    if (propfb(ifac,ipprob(ifnet)).le.xlimit) then
      iok = iok + 1
      flunmn = min(flunmn,propfb(ifac,ipprob(ifnet)))
      write(nfecra,4000)ifac,izfrad(ifac),itypfb(ifac)
    endif
  enddo

  if (iok.ne.0) then
    write(nfecra,4100) flunmn
    call csexit (1)
    !==========
  endif

!--> Integration du flux net sur les differentes zones de frontieres
!     IFLUX sert en parallele pour reperer les zones existantes

  do izone = 1, nozrdm
    flux(izone) = 0.d0
    iflux(izone) = 0
  enddo
  do ifac = 1,nfabor
    izone = izfrad(ifac)
    flux(izone) = flux(izone) + propfb(ifac,ipprob(ifnet))*surfbn(ifac)
    iflux(izone) = 1
  enddo
  if(irangp.ge.0) then
    call parrsm(nozarm,flux )
    call parimx(nozarm,iflux)
  endif


  write(nfecra,5000)
  write(nfecra,5010)
  do izone = 1, nozarm
    if(iflux(izone).eq.1) then
      write(nfecra,5020) izone,flux(izone)
    endif
  enddo
  write(nfecra,5000)


!--> Integration de la densite de flux net aux frontieres

  aa = zero
  do ifac = 1,nfabor
    aa =  aa + propfb(ifac,ipprob(ifnet)) * surfbn(ifac)
  enddo
  if(irangp.ge.0) then
    call parsom(aa)
  endif
  write(nfecra,5030) aa

!===============================================================================
! 7. TERMES SOURCES RADIATIFS IMPLICITE ET EXPLICITE
!===============================================================================


!===============================================================================
! 7.1 TERMES SOURCES RADIATIFS SEMI-ANALYTIQUES
!===============================================================================

  if (idiver.ge.0) then

!--> On stocke dans le tableau de travail W9 le CP
!    Attention : il faut conserver W9 dans la suite de la routine,
!    car son contenu est utilise plus loin

    if (icp.gt.0) then
      do iel = 1,ncel
        w9(iel) = 1.d0/propce(iel,ipproc(icp))
      enddo
    else
      do iel = 1,ncel
        w9(iel) = 1.d0/cp0
      enddo
    endif

    do iel = 1,ncel

!--> part d'absorption du terme source explicite

      propce(iel,ipproc(iabs(1))) =                          &
propce(iel,ipproc(icak(1)))*propce(iel,ipproc(itsre(1)))

!--> part d'emission du terme source explicite

      propce(iel,ipproc(iemi(1))) = -4.d0 *                  &
             propce(iel,ipproc(itsri(1)))

    enddo

! Combustion CP : On rajoute la contribution des particules
    if ( ippmod(icp3pl).ge.0 ) then
      do icla = 1,nclacp
        ipcla = 1+icla
        do iel = 1,ncel
!         Fluide
          propce(iel,ipproc(iabs(1))) =                      &
             propce(iel,ipproc(iabs(1)))                     &
           + propce(iel,ipproc(ix2(icla))) *                      &
 propce(iel,ipproc(icak(ipcla)))*propce(iel,ipproc(itsre(1)))
!
          propce(iel,ipproc(iemi(1))) =                      &
                                propce(iel,ipproc(iemi(1)))  &
                          - 4.0d0*propce(iel,ipproc(ix2(icla)))   &
                          * propce(iel,ipproc(itsri(ipcla)))
!         Particule
          propce(iel,ipproc(iabs(ipcla))) =                       &
                 propce(iel,ipproc(icak(ipcla))) *                &
                 propce(iel,ipproc(itsre(1)))
          propce(iel,ipproc(iemi(ipcla))) = - 4.0d0 *             &
                 propce(iel,ipproc(itsri(ipcla)))
          propce(iel,ipproc(itsre(ipcla))) =                      &
                 propce(iel,ipproc(iabs(ipcla))) +                &
                   propce(iel,ipproc(iemi(ipcla)))
        enddo
      enddo

! Combustion Fuel : On rajoute la contribution des particules

    elseif ( ippmod(icfuel).ge.0 ) then
      do icla = 1,nclafu
        ipcla = 1+icla
        do iel = 1,ncel
!         Fluide
          propce(iel,ipproc(iabs(1))) =                      &
             propce(iel,ipproc(iabs(1)))                     &
           + rtpa(iel,isca(iyfol(icla))) *                        &
 propce(iel,ipproc(icak(ipcla)))*propce(iel,ipproc(itsre(1)))
!
          propce(iel,ipproc(iemi(1))) =                      &
                                propce(iel,ipproc(iemi(1)))  &
                          - 4.0d0*rtpa(iel,isca(iyfol(icla)))     &
                          * propce(iel,ipproc(itsri(ipcla)))
!         Particule
          propce(iel,ipproc(iabs(ipcla))) =                       &
                 propce(iel,ipproc(icak(ipcla))) *                &
                 propce(iel,ipproc(itsre(1)))
          propce(iel,ipproc(iemi(ipcla))) = - 4.0d0 *             &
                 propce(iel,ipproc(itsri(ipcla)))
          propce(iel,ipproc(itsre(ipcla))) =                      &
                 propce(iel,ipproc(iabs(ipcla))) +                &
                   propce(iel,ipproc(iemi(ipcla)))
        enddo
      enddo

    endif

!--> Premiere methode pour le calcul du terme source explicite :
!    il est calcule comme la somme des termes d'absorption et d'emission
!    (il faudra multiplier ce terme par VOLUME(IEL) dans COVOFI->RAYSCA)

    do iel = 1,ncel
      propce(iel,ipproc(itsre(1))) =                         &
 propce(iel,ipproc(iabs(1)))+propce(iel,ipproc(iemi(1)))
    enddo

!--> Terme source implicite,
!    (il faudra multiplier ce terme par VOLUME(IEL) dans COVOFI->RAYSCA)

  if ( ippmod(icod3p).eq.-1 .and.                                 &
       ippmod(icoebu).eq.-1       ) then

!         Rayonnement standard, flamme CP ou fuel

    do iel = 1,ncel
      propce(iel,ipproc(itsri(1))) =                         &
       -16.d0*propce(iel,ipproc(icak(1))) *stephn *          &
         (tempk(iel,1)**3) * w9(iel)
    enddo

  else

!         Flamme de diffusion ou flamme de premelange

    do iel = 1,ncel
      propce(iel,ipproc(itsri(1))) =                         &
       -16.d0*stephn*propce(iel,ipproc(icak(1)))*            &
           propce(iel,ipproc(it3m))  *  w9(iel)
    enddo

  endif


! Combustion CP : On rajoute la contribution des particules

    if ( ippmod(icp3pl).ge.0 ) then
      do icla = 1,nclacp
        ipcla = 1+icla
        do iel = 1,ncel
          propce(iel,ipproc(itsri(1))) =                     &
                    propce(iel,ipproc(itsri(1))) -           &
                    16.d0*propce(iel,ipproc(icak(ipcla))) *       &
   propce(iel,ipproc(ix2(icla))) * stephn * (tempk(iel,ipcla)**3) &
                    / cp2ch(ichcor(icla))
          propce(iel,ipproc(itsri(ipcla))) =                      &
               -16.d0*propce(iel,ipproc(icak(ipcla))) *stephn     &
                  *(tempk(iel,ipcla)**3) / cp2ch(ichcor(icla))
        enddo
      enddo

! Combustion FUEL : On rajoute la contribution des particules

    elseif ( ippmod(icfuel).ge.0 ) then
      do icla = 1,nclafu
        ipcla = 1+icla
        do iel = 1,ncel
          propce(iel,ipproc(itsri(1))) =                     &
                    propce(iel,ipproc(itsri(1))) -           &
                    16.d0*propce(iel,ipproc(icak(ipcla))) *       &
     rtpa(iel,isca(iyfol(icla))) * stephn * (tempk(iel,ipcla)**3) &
                 /cp2fol
          propce(iel,ipproc(itsri(ipcla))) =                      &
               -16.d0*propce(iel,ipproc(icak(ipcla))) *stephn     &
                  *(tempk(iel,ipcla)**3) / cp2fol
        enddo
      enddo
    endif

  else
    do iel = 1,ncel
      propce(iel,ipproc(iabs(1)))  = zero
      propce(iel,ipproc(iemi(1)))  = zero
      propce(iel,ipproc(itsre(1))) = zero
      propce(iel,ipproc(itsri(1))) = zero
    enddo
  endif

!===============================================================================
! 7.2 TERME SOURCE RADIATIF EXPLICITE CONSERVATIF
!===============================================================================

! A partir d'ici COFRUA et COFRUB deviennent les CL pour la divergence

  if (idiver.eq.1 .or. idiver.eq.2) then

    do ifac = 1,nfabor
      cofrub(ifac) = zero
    enddo

!--> calcul de la divergence

    !    En periodique et parallele, echange avant calcul du gradient
    if (irangp.ge.0.or.iperio.eq.1) then
      call synvec &
      !==========
    ( propce(1,ipproc(iqx)), propce(1,ipproc(iqy)), propce(1,ipproc(iqz)) )
    endif


!    Donnees pour le calcul de la divergence

    inc     = 1
    iccocg  = 1
    imligp  = -1
    iwarnp  = iimlum
    epsrgp  = 1.d-8
    climgp  = 1.5d0
    extrap  = 0.d0
    nswrgp  = 100

!------->>direction X

    do ifac = 1,nfabor
      cofrua(ifac) = propfb(ifac,ipprob(ifnet)) *surfbo(1,ifac) / surfbn(ifac)
    enddo

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
!    sera a revoir pour la periodicite de rotation
    ivar0 = 0
    call grdcel                                                   &
    !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   propce(1,ipproc(iqx))    , cofrua , cofrub ,                   &
   w1     , w2     , w3     ,                                     &
   ra     )

    do iel = 1,ncel
      propce(iel,ipproc(itsre(1))) = - w1(iel)
    enddo

!------->>direction Y

    do ifac = 1,nfabor
      cofrua(ifac) = propfb(ifac,ipprob(ifnet)) *surfbo(2,ifac) / surfbn(ifac)
    enddo

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
!    sera a revoir pour la periodicite de rotation
    ivar0 = 0
    call grdcel                                                   &
    !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   propce(1,ipproc(iqy))    , cofrua , cofrub ,                   &
   w1     , w2     , w3     ,                                     &
   ra     )

    do iel = 1,ncel
      propce(iel,ipproc(itsre(1))) = propce(iel,ipproc(itsre(1))) - w2(iel)
    enddo

!------->>direction Z

    do ifac = 1,nfabor
      cofrua(ifac) = propfb(ifac,ipprob(ifnet)) *surfbo(3,ifac) / surfbn(ifac)
    enddo

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
!    sera a revoir pour la periodicite de rotation
    ivar0 = 0
    call grdcel                                                   &
    !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   propce(1,ipproc(iqz))    , cofrua , cofrub ,                   &
   w1     , w2     , w3     ,                                     &
   ra     )

    do iel = 1,ncel
      propce(iel,ipproc(itsre(1))) =                         &
          propce(iel,ipproc(itsre(1))) - w3(iel)
    enddo

! Fin du calcul de la divergence
  endif


!===============================================================================
! 7.3 TERME SOURCE RADIATIF EXPLICITE SEMI-ANALYTIQUE CORRIGE
!===============================================================================


  if (idiver.eq.2) then

!---> comparaison des termes sources semi-analytique et conservatif

    aa = zero
    do iel = 1,ncel
      aa = aa + propce(iel,ipproc(itsre(1))) * volume(iel)
    enddo

    bb = zero
    do iel = 1,ncel
      bb = bb + (propce(iel,ipproc(iabs(1))) +               &
            propce(iel,ipproc(iemi(1)))) * volume(iel)
    enddo

    if(irangp.ge.0) then
      call parsom(aa)
      call parsom(bb)
    endif

    aa = aa/bb

!---> correction du terme source semi-analytique par le conservatif

    do iel = 1,ncel
      propce(iel,ipproc(itsre(1))) =                         &
       (propce(iel,ipproc(iabs(1)))  +                       &
         propce(iel,ipproc(iemi(1)))) * aa
    enddo

  endif

!===============================================================================
! 7.4 FINALISATION DU TERME SOURCE EXPLICITE
!===============================================================================

  if (idiver.ge.0) then

!--> Integration volumique du terme source explicite
!    Le resultat de cette integration DOIT etre le meme que l'integration
!    surfacique de la densite de flux net radiatif faite plus haut
!    si  IDIVER = 1 ou 2

    aa = zero
    do iel = 1,ncel
      aa = aa + propce(iel,ipproc(itsre(1))) * volume(iel)
    enddo
    if(irangp.ge.0) then
      call parsom(aa)
    endif
    write(nfecra,5040) aa
    write(nfecra,5050)
    write(nfecra,5000)

!--> Correction du terme source explicite dans raysca pour permettre un
!    post-processing correct du terme source explicite
!    lorsque la variable transportee est la temperature
!    (pour les calculs en combustion la variable transportee est toujours
!    l'enthalpie)

  else
    write(nfecra,5000)
  endif

! Free memory
deallocate(viscf, viscb)
deallocate(dam, xam)
deallocate(drtp, smbrs, rovsdt)
deallocate(w1, w2, w3)
deallocate(w4, w5, w6)
deallocate(w7, w8, w9)
deallocate(w10)
deallocate(tempk)
deallocate(cofrua, cofrub)
deallocate(flurds, flurdb)

!--------
! FORMATS
!--------

 1000 FORMAT (/, 3X,'** INFORMATIONS SUR LE TERME SOURCE RADIATIF',/,   &
           3X,'   -----------------------------------------' )
 1100 FORMAT (/, 3X,'   Calcul effectue en rayonnement transparent'  ,/)

 2010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE RAYONNEMENT EST ACTIVE AVEC LE MODELE DOM.           ',/,&
'@      LA VALEUR MINIMALE DU COEFFICIENT D ABSORPTION A EST  ',/,&
'@      EGALE A ', E14.5                                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE RAYONNEMENT EST ACTIVE AVEC LE MODELE P-1.           ',/,&
'@      LE COEFFICIENT D''ABSORBTION DOIT ETRE STRICTEMENT    ',/,&
'@      SUPERIEUR A ZERO.                                     ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@                CONDITIONS AUX LIMITES MAL RENSEIGNEES      ',/,&
'@                                                            ',/,&
'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )
 3100 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    LES COEFFICIENTS DE CONDITIONS AUX LIMITES (COFRUA)     ',/,&
'@    NE SONT PAS RENSEIGNES POUR CERTAINES                   ',/,&
'@        FACES DE BORD                                       ',/,&
'@                                                            ',/,&
'@        Valeur minimale COFRUA ',E14.5                       ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usray5.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3200 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    LES COEFFICIENTS DE CONDITIONS AUX LIMITES (COFRUB)     ',/,&
'@    NE SONT PAS RENSEIGNES POUR CERTAINES                   ',/,&
'@        FACES DE BORD                                       ',/,&
'@                                                            ',/,&
'@        Valeur minimale COFRUB ',E14.5                       ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usray5.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3500 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE RAYONNEMENT EST ACTIVE.                              ',/,&
'@                                                            ',/,&
'@    Le scalaire ',I10   ,' devrait etre la temperature ou   ',/,&
'@      l''enthalpie. On a ISCSTH = ',I10                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4000 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT (FLUNET    NON RENSEIGNE)       ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )
 4100 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    LE FLUNET    N''EST PAS RENSEIGNEE POUR CERTAINES       ',/,&
'@        FACES DE BORD                                       ',/,&
'@                                                            ',/,&
'@        Valeur minimale ',E14.5                              ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usray5.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 5000 format('-----------------------------------------------------',   &
          '--------------')

 5010 format('Zone         Flux net radiatif (Watt) (normale',          &
          ' unitaire sortante)')

 5020 format(i6,13x,e10.4)

 5030 format('Flux net radiatif sur toutes les frontieres  Fnet = ',    &
           E10.4,' Watt')

 5040 format('Integrale volumique du terme source radiatif Srad = ',    &
           E10.4,' Watt')

 5050 format('(Si IDIVER = 1 ou 2 alors on doit avoir Srad = -Fnet)')

 6000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT APPROXIMATION P-1  (RAYDOM)     ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    LA LONGUEUR OPTIQUE DU MILIEU SEMI-TRANSPARENT          ',/,&
'@      DOIT AU MOINS ETRE DE L''ORDRE DE L''UNITE POUR ETRE  ',/,&
'@      DANS LE DOMAINE D''APPLICATION DE L''APPROXIMATION P-1',/,&
'@    CELA NE SEMBLE PAS ETRE LE CAS ICI.                     ',/,&
'@                                                            ',/,&
'@    LE COEFFICIENT D''ABSORPTION MINIMUM POUR ASSURER CETTE ',/,&
'@      LONGUEUR OPTIQUE EST XKMIN = ',E10.4                   ,/,&
'@    CETTE VALEUR N''EST PAS ATTEINTE POUR ', E10.4,'%       ',/,&
'@      DES CELLULES DU MAILLAGE.                             ',/,&
'@    LE POURCENTAGE DE CELLULES DU MAILLAGE POUR LESQUELLES  ',/,&
'@      ON ADMET QUE CETTE CONDITION SOIT VIOLEE EST IMPOSE   ',/,&
'@      PAR DEFAUT OU DANS USINI1 A XNP1MX = ', E10.4,'%      ',/,&
'@                                                            ',/,&
'@    Verifier les valeurs du coefficient d''absorption CK    ',/,&
'@      dans l''interface ou le modifier dans USRAY3.         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

end subroutine

