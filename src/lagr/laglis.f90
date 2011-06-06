!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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

subroutine laglis &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , statis , stativ , tslagr , parbor ,          &
   ra     )

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     ECITURE DES INFOS DANS LE LISTING

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! tr ! <-- ! cumul pour les moyennes des                    !
!(ncelet,nvlsta    !    !     !   statistiques volumiques                      !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,ntersl    !    !     !   lagrangien sur la phase porteuse             !
! parbor           ! tr ! <-- ! infos sur interaction des particules           !
!(nfabor,nvisbr    !    !     !   aux faces de bord                            !
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
use numvar
use optcal
use cstnum
use entsor
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep)
integer          ia(*)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifinia, ifinra
integer          ifac , iel , ivf , itabvr , nbrcel
integer          ivff , iflu , icla , ii , nb
double precision aa , bb , gmax , gmin , gmoy
character        chcond*16

double precision, allocatable, dimension(:) :: tabvr

!===============================================================================
!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Initialize variables to avoid compiler warnings

if (nbpart.ne.0) then
  aa = 100.d0 / dble(nbpart)
else
  aa = 0.d0
endif

!===============================================================================
! 2. AFFICHAGE LISTING
!===============================================================================

write (nfecra,1000)

! NOMBRE DE PARTICULES

write(nfecra,1001)
write(nfecra,1010) iplas , iplar
write(nfecra,1001)
write(nfecra,1020)
write(nfecra,1002)
write(nfecra,1031) nbpnew, dnbpnw
if (iroule.ge.1) then
  write(nfecra,1037) npcsup, dnpcsu
  write(nfecra,1032) npclon, dnpclo
  write(nfecra,1034) npkill, dnpkil
endif
if (iphyla.eq.2 .and. iencra.eq.1) then
  write(nfecra,1038) npencr, dnpenc
endif
write(nfecra,1033) nbpout-nbperr, (dnbpou-dnbper)
write(nfecra,1039) nbpdep, dnbdep
write(nfecra,1035) nbperr, dnbper
write(nfecra,1036) nbpart, dnbpar
write(nfecra,1001)
if (nbptot.gt.0) then
  write(nfecra,1050) (nbpert*100.d0)/dble(nbptot)
  write(nfecra,1001)
endif

! DEBIT SUR CHAQUE ZONE

write(nfecra,7000)
write(nfecra,1002)
do ii = 1,nfrlag
  nb = ilflag(ii)

  if ( iusclb(nb) .eq. ientrl) then
     CHCOND = 'ENTREE'
  else if ( iusclb(nb) .eq. irebol) then
     CHCOND = 'REBOND'
  else if ( iusclb(nb) .eq. isortl) then
     CHCOND = 'SORTIE'
  else if ( iusclb(nb) .eq. idepo1 .or.                           &
       iusclb(nb) .eq. idepo2 .or.                                &
       iusclb(nb) .eq. idepo3     ) then
    CHCOND = 'DEPOSITION'
  else if ( iusclb(nb) .eq. iencrl) then
    CHCOND = 'ENCRASSEMENT'
  else if ( iusclb(nb) .eq. idepfa) then
    CHCOND = 'FORCES_CHIMIQUES'
  else
    CHCOND = 'UTILISATEUR'
  endif

  write(nfecra,7001) nb,deblag(nb)/dtp,chcond
enddo
write(nfecra,1001)

! STATISTIQUES VOLUMIQUES

if (istala.eq.1) then
  write(nfecra,2000)
  write(nfecra,1002)
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
    write(nfecra,1001)

    if (nvlsta.gt.0) then
      write(nfecra,3010)
      write(nfecra,1002)

      ! Allocate a work array
      allocate(tabvr(ncelet))

!     MOYENNE

      do ivf = 1, nvlsta

        ivff = ivf
        icla = 0
        iflu = 0

        call uslaen                                               &
        !==========
 ( nvar   , nscal  , nvlsta ,                                     &
   ivff   , ivff   , ivff   , iflu   , ilpd   , icla   ,          &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statis , stativ , tabvr  ,                   &
   ra     )

        nbrcel = 0

        gmax = -grand
        gmin =  grand
        gmoy =  0.d0

        do iel = 1,ncel
          bb = tabvr(iel)
          if (statis(iel,ilpd).gt.seuil) then
            nbrcel = nbrcel + 1
            gmax = max (gmax, bb)
            gmin = min (gmin, bb)
            gmoy = gmoy + bb
          endif
        enddo

        if (nbrcel.gt.0) then
          gmoy = gmoy /dble(nbrcel)
        else
          gmax =  0.d0
          gmin =  0.d0
          gmoy =  0.d0
        endif
        write(nfecra,3020) nomlag(ivf),  gmin, gmax, gmoy

      enddo

      ! Free memory
      deallocate(tabvr)

    endif
  endif

  write(nfecra,1001)

endif

! STATISTIQUES PARIETALES

if (iensi3.eq.1) then

  write(nfecra,5000)
  write(nfecra,1002)
  if (isttio.eq.1) then
    if (iplas.ge.nstbor) then
      write(nfecra,5020) npstf
    else
      write(nfecra,5010) nstbor
    endif
  endif
  write(nfecra,5030) npstft
  write(nfecra,1001)

  if (nvisbr.gt.0) then
    write(nfecra,6000)
    write(nfecra,1002)

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

    do ivf = 1, nvisbr

      ivff = ivf
      call lagstf                                                 &
      !==========
       ( ncelet , nfabor , nvisbr ,                               &
         ivff   ,                                                 &
         gmin   , gmax   , gmoy   ,                               &
         parbor , tabvr  )

      write(nfecra,6010) nombrd(ivf),  gmin, gmax, gmoy

    enddo

    ! Free memory
    if (allocated(tabvr)) deallocate(tabvr)

    write(nfecra,1001)

  endif

endif

! INFO SUR LE COUPLAGE RETOUR

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

 1000 format(3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN',/3X,         &
          '   -------------------------------------')

 1001 format('-----------------------------------------------------',   &
       '----------')

 1002 format('   ---------------------------------------------------',  &
       '-----')

 1010 format('Iters Lagrangiennes absolues/relatives : ',               &
         I10,' /',I10)

 1020 format('   Pour cette iteration, nombre de particules',/,   &
       '   (sans et avec leur poids statistique) :')
 1031 format('ln  nouvelles injectees                ',I8,3X,E14.5)
 1032 format('ln  nouvelles par clonage              ',I8,3X,E14.5)
 1033 format('ln  sorties, ou deposees et supprimees ',I8,3X,E14.5)
 1034 format('ln  eliminees par roulette russe       ',I8,3X,E14.5)
 1035 format('ln  perdues par erreur de reperage     ',I8,3X,E14.5)
 1036 format('ln  total restantes en fin de passage  ',I8,3X,E14.5)
 1037 format('ln  qui ont subit le clonage           ',I8,3X,E14.5)
 1038 format('ln  de charbon encrassees              ',I8,3X,E14.5)
 1039 format('ln  deposees                           ',I8,3X,E14.5)

 1050 format('% de particules perdues (suites comprises) : ',E10.4)

 2000 format('   Statistiques volumiques :')

 2005 format('Debut des statistiques a l''iteration absolue :   ',I10)

 2010 format('Nombre d''iterations dans les stats instationnaires : ',  &
        i10)

 2020 format('Nombre d''iterations total dans les statistiques :',I10)

 2030 format('RAZ des stats (debut du calcul stationnaire : ',I10,' )')

 2040 format('Nombre d''iterations dans les stats stationnaires :',I9)

 3010 format('   Stats volumiques  Valeur min    Valeur max    ',       &
       'Valeur moy')

 3020 format('lc  ',A13,2X,E12.5,2X,E12.5,2X,E12.5)

 4000 format('   Termes sources de couplage-retour instationnaires :')

 4010 format('   Termes sources de couplage-retour :')

 4020 format('RAZ des TS (debut du calcul stationnaire : ',I10,')')

 4030 format('Nombre d''iterations dans les TS stationnaires :',I10)

 4050 format('Taux volumiques max de particules : ',E14.5)

 4060 format('Taux massiques max de particules :  ',E14.5)

 4070 format('Nbr de cellules qui ont un taux volumique > 0.8 :',I10)

 5000 format('   Statistiques aux frontieres :')

 5010 format('RAZ des stats aux frontieres (debut stationnaire : ',     &
        I8,')')

 5020 format('Nbr d''iters des stats frontieres stationnaires : ',I10)

 5030 format('Nbr d''iters total dans les stats aux frontieres :',I10)

 6000 format('   Stats frontieres  Valeur min    Valeur max    ',       &
       'Valeur moy')

 6010 format('lp  ', A13, 2X, E12.5, 2X, E12.5, 2X, E12.5)

 7000 format(3X,'ZONE          DEBIT(kg/s)       TYPE CL    ')

 7001 format(2x, i3, 10x, e12.5, 9x, a16)

!====
! FIN
!====

end subroutine

!-----------------------------------------------------------------------
!Iterations Lagrangiennes absolues/relatives : 1234567 /1234567
!-----------------------------------------------------------------------
!   Pour cette iteration, nombre particules :
!   ---------------------------------------------------------
!ln   nouvelles injectees                12345678  0.3067e+04
!ln   qui ont ete clonees                12345678
!ln   nouvelles par clonage              12345678
!ln   sortantes du domaine               12345678
!ln   eliminees par roulette russe       12345678
!ln   perdues par erreur de reperage     12345678
!ln   total restantes en fin de passage  12345678
!-----------------------------------------------------------------------
!Poucentage de particules perdues pour ce calcul : 1234567890
!-----------------------------------------------------------------------
!   Statistiques volumiques :
!   ---------------------------------------------------------
!Debut des stat a l'iteration Lagrangienne absolue : 1234567
!Nombre de passages dans les stat instationnaires : 1234567
!Nombre de passages total dans les stat : 1234567
!RAZ des statistiques stationnaires (debut du cumul : 1234567)
!Nombre de passages dans les stat stationnaires : 1234567
!-----------------------------------------------------------------------
!   Stat volumiques   Valeur min    Valeur max    Valeur moy
!   ---------------------------------------------------------
!lc   1234567890123  -0.47302e+03  -0.30677e+04  -0.30677e+04
!-----------------------------------------------------------------------
!   Termes sources de couplage-retour
!   ---------------------------------------------------------
!Taux volumiques max de particules :
!Taux massiques max de particules :
!Nbr de cellules qui ont un taux volumique > 0.8 :
!-----------------------------------------------------------------------
!   Statistiques aux frontieres :
!   ---------------------------------------------------------
!RAZ des stat aux frontieres (debut stationnaire :   1234567890)
!Nbr d''iter des stat aux frontieres stationnaires : 1234567890
!Nbr d''iter total dans les stat aux frontieres :    1234567890
!-----------------------------------------------------------------------
!   Stat volumiques   Valeur min    Valeur max    Valeur moy
!   ---------------------------------------------------------
!lp   nombreImpacts       -0.47302e+03  -0.30677e+04
!-----------------------------------------------------------------------

