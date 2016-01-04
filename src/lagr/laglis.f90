!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine laglis &
!================

 ( nvlsta , nvisbr ,                                                           &
   statis , parbor )

!===============================================================================
! Purpose:
! --------
!   Subroutine of the Lagrangian particle-tracking module :
!   -------------------------------------------------------
!
!   Writes the information about the particle-tracking calculation
!   in the main listing file.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! statis           ! tr ! <-- ! cumul pour les moyennes des                    !
!  (ncelet,nvlsta) !    !     !   statistiques volumiques                      !
! parbor           ! tr ! <-- ! infos sur interaction des particules           !
!  (nfabor,nvisbr) !    !     !   aux faces de bord                            !
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
use cstnum
use entsor
use parall
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvlsta , nvisbr

double precision statis(ncelet,nvlsta)
double precision parbor(nfabor,nvisbr)

! Local variables

integer          ifac , iel , ivf
integer          ivff , iflu , icla , ii , nb, nbfr, nbrcel
double precision aa , gmax , gmin
character        chcond*16

double precision, allocatable, dimension(:) :: tabvr,tabvrfou

integer nbpartall, nbpoutall, nbperrall, nbpdepall, npencrall
integer nbpresall, nbptotall, npclonall, npcsupall, nbpertall
integer npkillall

double precision dnbparall, dnbperall, dnbpouall
double precision dnbdepall, dnpencall, dnbpnwall, dnbresall
double precision dnpcloall, dnpcsuall, dnpkilall
double precision debloc(2)

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

ipass = ipass + 1

! Initialize variables to avoid compiler warnings

nbrcel = 0

if (nbpart.ne.0) then
  aa = 100.d0 / dble(nbpart)
else
  aa = 0.d0
endif

!===============================================================================
! 2. AFFICHAGE LISTING
!===============================================================================

write (nfecra,1000)

! Parallelism management

nbpartall = nbpart
nbpoutall = nbpout
nbperrall = nbperr
nbpdepall = nbpdep
npencrall = npencr
nbpresall = nbpres
nbptotall = nbptot
nbpertall = nbpert

dnbparall = dnbpar
dnbpouall = dnbpou
dnbperall = dnbper
dnbdepall = dnbdep
dnpencall = dnpenc
dnbpnwall = dnbpnw
dnbresall = dnbres

if (irangp.ge.0) then

   call parcpt(nbpartall)
   call parcpt(nbpoutall)
   call parcpt(nbperrall)
   call parcpt(nbpdepall)
   call parcpt(npencrall)
   call parcpt(nbpresall)
   call parcpt(nbptotall)
   call parcpt(nbpertall)

   call parsom(dnbparall)
   call parsom(dnbpouall)
   call parsom(dnbperall)
   call parsom(dnbdepall)
   call parsom(dnpencall)
   call parsom(dnbpnwall)
   call parsom(dnbresall)

endif

! Number of particles

write(nfecra,1003)
write(nfecra,1010) iplas , iplar
write(nfecra,1003)
write(nfecra,1020)
write(nfecra,1003)
write(nfecra,1031) nbpnew, dnbpnwall
if (iroule.ge.1) then
  npcsupall = npcsup
  dnpcsuall = dnpcsu
  npclonall = npclon
  dnpcloall = dnpclo
  npkillall = npkill
  dnpkilall = dnpkil
  if (irangp.ge.0) then
    call parcpt(npcsupall)
    call parcpt(npclonall)
    call parcpt(npkillall)
    call parsom(dnpcsuall)
    call parsom(dnpcloall)
    call parsom(dnpkilall)
  endif
  write(nfecra,1037) npcsupall, dnpcsuall
  write(nfecra,1032) npclonall, dnpcloall
  write(nfecra,1034) npkillall, dnpkilall
endif
if (iphyla.eq.2 .and. iencra.eq.1) then
  write(nfecra,1038) npencrall, dnpencall
endif
write(nfecra,1033) nbpoutall-nbperr, (dnbpou-dnbper)
write(nfecra,1039) nbpdepall, dnbdepall

if (ireent.gt.0) then
write(nfecra,1040) nbpresall, dnbresall
endif

write(nfecra,1035) nbperrall, dnbperall
write(nfecra,1036) nbpartall, dnbparall
if (nbptotall.gt.0) then
  write(nfecra,1050) (nbpertall*100.d0)/dble(nbptotall)
  write(nfecra,1001)
endif

! Flow rate for each zone

write(nfecra,7000)

nbfr = 0
do ii = 1,nfrlag
  if (ilflag(ii) .gt. nbfr) nbfr = ilflag(ii)
enddo
if (irangp.ge.0) then
  call parcmx(nbfr)
endif

do nb = 1,nbfr
  debloc(1) = 0.d0
  debloc(2) = 0.d0
  do ii = 1, nfrlag
    if (ilflag(ii) .eq. nb) then
      debloc(1) = 1.d0
      debloc(2) = deblag(nb)
    endif
  enddo
  if (irangp.ge.0) then
    call parrsm(2,debloc)
  endif
  if (debloc(1) > 0.5d0) then
    if (iusclb(nb) .eq. ientrl) then
      chcond = 'INLET'
    else if (iusclb(nb) .eq. irebol) then
      chcond = 'REBOUND'
    else if (iusclb(nb) .eq. isortl) then
      chcond = 'OUTLET'
    else if (iusclb(nb) .eq. idepo1 .or. iusclb(nb) .eq. idepo2) then
      chcond = 'DEPOSITION'
    else if (iusclb(nb) .eq. iencrl) then
      chcond = 'FOULING'
    else if (iusclb(nb) .eq. idepfa) then
      chcond = 'DLVO CONDITIONS'
    else if (iusclb(nb) .eq. isymtl) then
      chcond = 'DLVO CONDITIONS'
    else
      chcond = 'USER'
    endif
    write(nfecra,7001) nb,debloc(2)/dtp,chcond
  endif
enddo

write(nfecra,1001)

! Volumic statistics

if (istala.eq.1) then
  write(nfecra,2000)
  write(nfecra,1003)
  write(nfecra,2005) idstnt
  if (iplas.ge.idstnt) then
    if (isttio.eq.0) then
       write(nfecra,2010) npstt
    endif
    if (isttio.eq.1 .and. iplas.lt.nstist) then
      write(nfecra,2020) npstt
      write(nfecra,2030) nstist
    else if (isttio.eq.1 .and. iplas.ge.nstist) then
      write(nfecra,2020) npstt
      write(nfecra,2040) npst
    endif
    write(nfecra,1003)
    if (nvlsta.gt.0) then
      write(nfecra,3010)

      ! Allocate a work array
      allocate(tabvr(ncelet))

      ! Calculation of the averages
      do ivf = 1, nvlsta

        ivff = ivf
        icla = 0
        iflu = 0

        gmin = grand
        gmax = -grand

        call uslaen                                               &
        !==========
        ( nvlsta ,                                                &
          ivff   , ivff   , ivff   , iflu   , ilpd   , icla   ,   &
          tabvr  )

        if ((ivf.ne.ilfv).and.(ivf.ne.ilpd)) then

          nbrcel = 0

          do iel = 1,ncel
            if (statis(iel,ivf).gt.seuil) then
              gmax = max (gmax, tabvr(iel))
              gmin = min (gmin, tabvr(iel))
              nbrcel = nbrcel + 1
            endif
          enddo
          if (irangp.ge.0) then
            call parcpt(nbrcel)
          endif

          if (nbrcel.eq.0) then
            gmax =  0.d0
            gmin =  0.d0
          endif

        else

          do iel = 1,ncel
            gmax = max (gmax, tabvr(iel))
            gmin = min (gmin, tabvr(iel))
          enddo

        endif

        if (irangp.ge.0) then
          call parmin(gmin)
          call parmax(gmax)
        endif

        write(nfecra,3020) nomlag(ivf),  gmin, gmax

      enddo

      ! Free memory
      deallocate(tabvr)
    endif
  endif
  write(nfecra,1001)

endif

! Boundary statistics

if (iensi3.eq.1) then

  write(nfecra,5000)
  write(nfecra,1003)
  if (isttio.eq.1) then
    if (iplas.ge.nstbor) then
      write(nfecra,5020) npstf
    else
      write(nfecra,5010) nstbor
    endif
 endif
 write(nfecra,5030) npstft
 write(nfecra,1003)

  if (nvisbr.gt.0) then
    write(nfecra,6000)
    if (nvisbr.gt.1) then

      ! Allocate a work array
      allocate(tabvr(nfabor))

      do ifac = 1,nfabor
        if (parbor(ifac,inbr).gt.seuilf) then
          tabvr(ifac) = 1.d0 / parbor(ifac,inbr)
        else
          tabvr(ifac) = 0.d0
        endif
      enddo

    endif

    if (iencnbbd.eq.1) then
      allocate(tabvrfou(nfabor))
      do ifac = 1,nfabor
        if (parbor(ifac,iencnb).gt.seuilf) then
          tabvrfou(ifac) = 1.d0 / parbor(ifac,iencnb)
        else
          tabvrfou(ifac) = 0.d0
        endif
      enddo
    endif

    do ivf = 1, nvisbr

      ivff = ivf
      call lagstf                                                 &
      !==========
       ( ncelet , nfabor , nvisbr ,                               &
         ivff   ,                                                 &
         gmin   , gmax   ,                                        &
         parbor , tabvr  , tabvrfou )
      if (irangp.ge.0) then
        call parmin(gmin)
        call parmax(gmax)
      endif
      write(nfecra,6010) nombrd(ivf),  gmin, gmax
    enddo
    ! Free memory
    if (allocated(tabvr)) deallocate(tabvr)
    if (allocated(tabvrfou)) deallocate(tabvrfou)
    write(nfecra,1001)
  endif

endif


 ! Information about two-way coupling

if (iilagr.eq.2) then

  if (isttio.eq.0) then
    write(nfecra,4000)
    write(nfecra,1002)

  else if (isttio.eq.1) then
    write(nfecra,4010)
    write(nfecra,1002)

    if (iplas.lt.nstits) then
      write(nfecra,4020) nstits
    else if (iplas.ge.nstist) then
      write(nfecra,4030) npts
    endif

  endif

  write(nfecra,4050) vmax
  write(nfecra,4060) tmamax
  write(nfecra,4070) ntxerr

  write(nfecra,1001)

endif

!===============================================================================

!--------
! FORMATS
!--------


1000 format(3X,'** INFORMATION ON THE LAGRANGIAN CALCULATION',/3X,         &
          '   ------------------------------------------')

1001 format('-----------------------------------------------------',   &
          '----------')

1002 format('   ---------------------------------------------------',  &
          '-----')

1003 format('  ')

1010 format('Lagrangian iteration n° (absolute/relative) : ',          &
          I10,' /',I10)

1020 format('   For the current iteration, number of particles',/,     &
          '   (with and without statistical weight) :')
1031 format('ln  newly injected                           ',I8,3X,E14.5)
1032 format('ln  new by cloning                           ',I8,3X,E14.5)
1033 format('ln  out, or deposited and eliminated         ',I8,3X,E14.5)
1034 format('ln  eliminated by russian roulette           ',I8,3X,E14.5)
1035 format('ln  lost in the location stage               ',I8,3X,E14.5)
1036 format('ln  total number at the end of the time step ',I8,3X,E14.5)
1037 format('ln  which have undergone cloning             ',I8,3X,E14.5)
1038 format('ln  coal particles fouled                    ',I8,3X,E14.5)
1039 format('ln  deposited                                ',I8,3X,E14.5)
1040 format('ln  resuspended                              ',I8,3X,E14.5)
1050 format('% of lost particles (restart(s) included) :  ',e11.4)

2000 format('   Volume statistics :')

2005 format('Start of calculation from absolute Lagrangian iteration n°:   ', i10)

2010 format('Number of iterations in unsteady statistics: ', i10)

2020 format('Total number of iterations in the statistics:', i10)

2030 format('Start of steady-state statistics from Lagrangian iteration n°: ', i10,' )')

2040 format('Number of iterations in steady-state statistics :', i9)

3010 format('                            Min value    Max value    ')

3020 format('lc  ',A20,2X,E12.5,2X,E12.5,2X,E12.5)

4000 format('   Unsteady two-way coupling source terms:')

4010 format('   Two-way coupling source terms:')

4020 format('Reset of the source terms (Start of steady-state at:): ',I10,')')

4030 format('Number of iterations for the steady-state source terms:',I10)

4050 format('Maximum particle volume fraction : ',E14.5)

4060 format('Maximum particle mass fraction :  ',E14.5)

4070 format('Number of cells with a part. volume fraction greater than 0.8 :',I10)

5000 format('   Boundary statistics :')

5010 format('Start of steady-state statistics from Lagrangian iteration n°: ', I8,')')

5020 format('Number of iterations in steady-state statistics: ',I10)

5030 format('Total number of iterations in the statistics:',I10)

6000 format('                           Min value    Max value    ')

6010 format('lp  ', A20, 2X, E12.5, 2X, E12.5, 2X, E12.5)

7000 format(3X,'Zone     Mass flow rate(kg/s)      Boundary type    ')

7001 format(2x, i3, 10x, e12.5, 9x, a)

!====
! FIN
!====

end subroutine

