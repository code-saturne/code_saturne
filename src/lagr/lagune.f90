!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine lagune &
!================

 ( lndnod ,                                                       &
   nvar   , nscal  ,                                              &
   dt     , propce )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!   Sous-programme principal du module de modelisation Lagrangienne
!   des ecoulements diphasiques a inclusions dispersees.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use entsor
use cstphy
use cstnum
use parall
use period
use pointe
use lagpar
use lagdim, only: ntersl, nvlsta, nvisbr
use lagran
use mesh
use field
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

integer          lndnod
integer          nvar   , nscal

double precision dt(ncelet)
double precision propce(ncelet,*)

! Local variables

integer          ip     , npt    , iok
integer          iel    , ivar, ifld
integer          npar1  , npar2
integer          modntl
integer          iprev

integer          ifac

double precision visccf, romf
double precision ustarmoy, surftot, surfb

double precision, allocatable, dimension(:) :: taup
double precision, allocatable, dimension(:,:) :: tlag, piil
double precision, allocatable, dimension(:,:,:) :: bx
double precision, allocatable, dimension(:,:) :: tempct
double precision, allocatable, dimension(:) :: tsfext
double precision, allocatable, dimension(:) :: cpgd1, cpgd2, cpght
double precision, allocatable, dimension(:) :: terbru
double precision, allocatable, dimension(:,:) :: gradpr
double precision, allocatable, dimension(:,:,:) :: gradvf
double precision, allocatable, dimension(:) :: croule
double precision, allocatable, dimension(:) :: w1, w2, w3

double precision, allocatable, save, dimension(:) :: vislen

double precision, allocatable, dimension(:) :: energt

double precision, allocatable, dimension(:) :: tempp

double precision, dimension(:), pointer :: cromf
double precision, dimension(:), pointer :: viscl

double precision, dimension(:), pointer :: cvar_scalt

integer ii
integer nbpartall, nbpper

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
! 0.  GESTION MEMOIRE ET COMPTEUR DE PASSAGE
!===============================================================================

! Allocate temporary arrays
allocate(gradpr(3,ncelet))
allocate(w1(ncelet), w2(ncelet), w3(ncelet))

! Allocate other arrays depending on user options
if (modcpl.gt.0) then
  allocate(gradvf(3,3,ncelet))
endif
if (iroule.eq.1) then
  allocate(croule(ncelet))
endif
if (idlvo.eq.1) then
  allocate(energt(nfabor))
  if (iclogst.eq.1 .or.  irough .eq.1) then
    allocate(tempp(nfabor))
  endif
endif

ipass = ipass + 1

if ((idepst.eq.1).and.(ipass.eq.1)) then
  allocate(vislen(nfabor))
  do ifac = 1, nfabor
    vislen(ifac) = grand
  enddo
endif

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

iplar = iplar + 1
iplas = iplas + 1

!-->Sur Champ fige Lagrangien :
!   values at previous time step = values at current time step
!   Rem : cette boucle pourrait etre faite au 1er passage
!         mais la presence de cs_user_extra_operations incite a la prudence...

if (iilagr.eq.3) then
  ifld = -1
  do ivar = 1,nvar
    if (ivarfl(ivar) .ne. ifld) then
      ifld = ivarfl(ivar)
      call field_current_to_previous(ifld)
    endif
  enddo
endif

!-->au premier passage relatif :

if (iplar.eq.1) then

!      Connectivite cellules -> faces + Alloc. structures en C

  call lagbeg                                                     &
  !==========
 ( nlayer , iphyla , idepst , irough , ireent , iclogst,          &
   nvls   , nbclst , icocel , itycel ,                            &
   jisor  , jisora , jirka  , jord1  ,                            &
   jrval  , jrpoi  , jrtsp  , jdp    , jmp    ,                   &
   jxp    , jyp    , jzp    , jup    , jvp    , jwp    ,          &
   juf    , jvf    , jwf    , jtaux  , jryplu ,                   &
   jrinpf , jdfac  , jimark ,                                     &
   jtp    , jhp    , jtf    , jmwat  , jmch   , jmck   ,          &
   jcp    , jrdck  , jrd0p  , jinch  , jrhock ,                   &
   jreps  , jdepo  , jnbasg , jnbasp , jfadh  , jmfadh ,          &
   jndisp , jclst  , jvls   )

  call lagr_update_pointers

  ! first initializations

  nbpart = 0;
  dnbpar = 0.d0

! --> if the deposition model is activated

  if (idepst.ge.1) then

     ustarmoy = 0.d0
     surftot = 0.d0

    ! boundary faces data

     call laggeo
     !==========

     if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then
       call field_get_val_s(iprpfl(ipproc(irom1)), cromf)
     else
       call field_get_val_s(icrom, cromf)
     endif

     call field_get_val_s(iprpfl(iviscl), viscl)

     ! Average friction velocity calculation
     do ifac = 1, nfabor

        if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then

           iel = ifabor(ifac)

           surfb = sqrt( surfbo(1,ifac)*surfbo(1,ifac)              &
                      +  surfbo(2,ifac)*surfbo(2,ifac)              &
                      +  surfbo(3,ifac)*surfbo(3,ifac) )

           ! the density pointer according to the flow location

           romf = cromf(iel)
           visccf = viscl(iel) / romf

           if ( uetbor(ifac).gt.1.d-15) then

              ustarmoy = (surftot * ustarmoy +  surfb * uetbor(ifac))   &
                       / (surftot + surfb)
              surftot = surftot +  surfb
              vislen(ifac) = visccf / uetbor(ifac)

           endif

        endif

     enddo

!  Average friction velocity display

     write(nfecra,4100) ustarmoy
!
  endif

endif

! Initialization

nbpnew = 0
npcsup = 0
npclon = 0
npkill = 0
npencr = 0
nbpout = 0
nbperr = 0
nbpdep = 0
nbpres = 0

dnbpnw = 0.d0
dnpcsu = 0.d0
dnpclo = 0.d0
dnpkil = 0.d0
dnpenc = 0.d0
dnbpou = 0.d0
dnbper = 0.d0
dnbdep = 0.d0
dnbres = 0.d0

!===============================================================================
! 1.bis  Initialization for the clogging model
!===============================================================================

if ( iclogst.eq.1 ) then

  if (iscalt.gt.0) call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

  do ifac = 1,nfabor
    iel = ifabor(ifac)

    if (iscalt.gt.0) then

      if (itherm.eq.1 .and. itpscl.eq.2) then
        tempp(ifac) = cvar_scalt(iel) + tkelvi
      else if (itherm.eq.1 .and. itpscl.eq.2) then
        tempp(ifac) = cvar_scalt(iel)
      else if (itherm.eq.2) then
        call usthht(1,cvar_scalt(iel),tempp(ifac))
      endif

    else
      tempp(ifac) = t0
    endif

  enddo

  call cloginit                                                   &
  !===========
  ( cstfar, epsvid, epseau, fion, jamlim, mporos, tempp,          &
    valen , phi1  , phi2  , cstham, dcutof, lambwl, kboltz )

endif

!===============================================================================
! 1.ter  Initialization for the roughness surface model
!===============================================================================

if ( irough .eq. 1 ) then

  if (iscalt.gt.0) call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

  do ifac = 1,nfabor
    iel = ifabor(ifac)

    if (iscalt.gt.0) then

      if (itherm.eq.1 .and. itpscl.eq.2) then
        tempp(ifac) = cvar_scalt(iel) + tkelvi
      else if (itherm.eq.1 .and. itpscl.eq.2) then
        tempp(ifac) = cvar_scalt(iel)
      else if (itherm.eq.2) then
        call usthht(1,cvar_scalt(iel),tempp(ifac))
      endif

    else
      tempp(ifac) = t0
    endif

  enddo

  call roughness_init                                &
  !===========
  ( cstfar, epsvid, epseau, fion, tempp,  valen    , &
    phi1  , phi2  , cstham, dcutof, lambwl, kboltz , &
    espasg , denasp , rayasp , rayasg)

endif

!===============================================================================
! 2.  MISE A JOUR DES NOUVELLES PARTICULES ENTREES DANS LE DOMAINE
!===============================================================================

! At the first time step we initialize particles to
! values at current time step and not at previous time step, because
! values at previous time step = initialization

if ( ntcabs.eq.1 ) then
  ! Use fields at current time step
  iprev = 0
else
  ! Use fields at previous time step
  iprev = 1
endif

call lagent                                                      &
!==========
( lndnod ,                                                       &
  nvar   , nscal  ,                                              &
  ntersl , nvlsta , nvisbr , iprev  ,                            &
  itycel , icocel , dlgeo  ,                                     &
  itypfb , itrifb , ifrlag ,                                     &
  dt     )

!===============================================================================
! 2.1 CALCUL DE LA FONCTION D'IMPORTANCE POUR LA ROULETTE RUSSE
!===============================================================================

if (iroule.ge.1) then

  call uslaru                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb ,                                              &
   dt     ,                                                       &
   croule ,                                                       &
   dispar , yplpar )

  iok = 0
  do iel = 1,ncel
    if (croule(iel).le.0.d0) iok = iok + 1
  enddo
  if (iok.ne.0) then
    write(nfecra,9001)
    call csexit (1)
    !==========
  endif

endif

!===============================================================================
! 3.  GESTION DU TEMPS QUI PASSE...
!===============================================================================

!-->Gestion du pas de temps Lagrangien

dtp = dt(1)

!-->Incrementation du TEMPS COURANT LAGRANGIEN

ttclag = ttclag + dtp

!-->Test pour savoir si le domaine contient des particules

nbpartall = nbpart

if ( irangp .ge. 0 ) then
  call parcpt(nbpartall)
endif

if (nbpartall.eq.0) goto 20

! Record particle's starting cell and rank, and reset order 1 switch

do ip = 1,nbpart
  ipepa(jisora,ip) = ipepa(jisor,ip)
  ipepa(jirka,ip) = irangp
  ipepa(jord1,ip) = 0
enddo

!===============================================================================
! 4.  GRADIENT DE PRESSION ET DE LA VITESSE FLUIDE
!===============================================================================

! At the first time step we initialize particles to
! values at current time step and not at previous time step, because
! values at previous time step = initialization (null gradients)

if (ntcabs.eq.1) then

  call laggra(0, gradpr, gradvf)
  !==========

else

  call laggra(1, gradpr, gradvf)
  !==========

endif

!===============================================================================
! 5. PROGRESSION DES PARTICULES
!===============================================================================

 10   continue

nor = mod(nor,nordre)
nor = nor + 1

! Allocations

allocate(taup(nbpart))
allocate(tlag(nbpart,3))
allocate(piil(nbpart,3))
allocate(bx(nbpart,3,2))
if (iilagr.eq.2) then
  allocate(tsfext(nbpart))
endif
if (iilagr.eq.2 .and. iphyla.eq.2 .and. ltsthe.eq.1) then
  allocate(cpgd1(nbpart))
  allocate(cpgd2(nbpart))
  allocate(cpght(nbpart))
endif
if ((iphyla.eq.1 .and. itpvar.eq.1) .or. iphyla.eq.2) then
  allocate(tempct(nbpart,2))
endif
if (lamvbr.eq.1) then
  allocate(terbru(nbpart))
endif

!---> Recopie des resultats de l'etape precedente :

if (nor.eq.1) then

  do ip = 1,nbpart
    call lagr_current_to_previous(ip)
  enddo

endif

!-----> COMPUTATION OF THE FLUID'S PRESSURE AND VELOCITY GRADIENT
!       AT N+1 (with values at current time step)

if (nor.eq.2 .and. iilagr.ne.3) then

  call laggra(0, gradpr, gradvf)
  !==========

endif

!-----> CALCUL DES CARACTERISTIQUES DES PARTICULES

if (nor.eq.1) then
  ! Use fields at previous time step
  iprev = 1
else
  ! Use fields at current time step
  iprev = 0
endif

call lagcar                                                     &
!==========
 ( nvar   , nscal  ,                                            &
   nvlsta , iprev  ,                                            &
   dt     ,                                                     &
   taup   , tlag   ,                                            &
   piil   , bx     , tempct , statis ,                          &
   gradpr , gradvf , w1     , w2     )


!---> INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES
!     POSITION, VITESSE FLUIDE, VITESSE PARTICULE

call lagesp                                                       &
!==========
   ( nvar   , nscal  ,                                            &
     ntersl , nvlsta , nvisbr ,                                   &
     dt     , propce ,                                            &
     statis , stativ , taup   , tlag   , piil   ,                 &
     bx     , tsfext ,                                            &
     gradpr , gradvf , terbru , vislen)

!---> INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES
!     LIEES AUX PHYSIQUES PARTICULIERES PARTICULAIRES

if (iphyla.eq.1 .or. iphyla.eq.2) then

  if ( nor.eq.1 ) then
    ! Use fields at previous time step
    iprev = 1
  else
    ! Use fields at current time step
    iprev = 0
  endif

  call lagphy                                                   &
  !==========
  ( ntersl , nvlsta , nvisbr ,                                  &
    iprev  , dt     , propce ,                                  &
    taup   , tlag   , tempct ,                                  &
    cpgd1  , cpgd2  , cpght  )

endif

!===============================================================================
! 6.  Couplage Retour - Calcul des termes sources
!===============================================================================

if (iilagr.eq.2 .and. nor.eq.nordre) then

  call lagcou                                                     &
  !==========
   ( ntersl ,                                                     &
     propce ,                                                     &
     taup   , tempct , tsfext ,                                   &
     cpgd1  , cpgd2  , cpght  ,                                   &
     w1     , w2   )

endif

!===============================================================================
! 7.  Calcul de la barrière d'énergie dans le cas DLVO
!===============================================================================

if (idlvo.eq.1) then

   call lagbar (energt)
   !==========

endif

! Deallocate arrays whose size is based on nbpart (which may change next)

deallocate(tlag)
deallocate(taup)
deallocate(piil)
deallocate(bx)
if (iilagr.eq.2) then
  deallocate(tsfext)
endif
if (iilagr.eq.2 .and. iphyla.eq.2 .and. ltsthe.eq.1) then
  deallocate(cpgd1)
  deallocate(cpgd2)
  deallocate(cpght)
endif
if ((iphyla.eq.1 .and. itpvar.eq.1) .or. iphyla.eq.2) then
  deallocate(tempct)
endif
if (lamvbr.eq.1) then
  deallocate(terbru)
endif

!===============================================================================
! 8.  Reperage des particules - Traitement des conditions aux limites
!     pour la position des particules
!===============================================================================

if (nor.eq.1) then

  !--> Si on est en instationnaire, RAZ des statistiques aux frontieres

  if (iensi3.eq.1) then

    if (isttio.eq.0 .or. (isttio.eq.1 .and. iplas.le.nstbor)) then
      tstatp = 0.d0
      npstf = 0
      do ii = 1,nvisbr
        do ifac = 1,nfabor
          parbor(ifac,ii) = 0.d0
        enddo
      enddo
    endif

    tstatp = tstatp + dtp
    npstf  = npstf  + 1
    npstft = npstft + 1

  endif

  call getbdy                                                     &
  !==========
 ( nflagm , nfrlag , injcon , ilflag , iusncl ,                   &
   iusclb , deblag , ifrlag )


  call dplprt                                                     &
  !==========
 ( nordre   , parbor   , iensi3   ,                               &
   inbr     , inbrbd   , iflm     , iflmbd   , iang     ,         &
   iangbd   , ivit     , ivitbd   , iencnb   , iencma   ,         &
   iencdi   , iencck   , iencnbbd , iencmabd , iencdibd ,         &
   iencckbd , inclg    , iscovc   ,                               &
   nusbor   , iusb     , vislen   , dlgeo    , energt   ,         &
   tprenc   , visref   , enc1     , enc2     , tkelvi)

  call lagr_update_pointers

  if (ierr.eq.1) then
    ntmabs = ntcabs
    write (nfecra,1000) ntmabs
    goto 20
  endif

endif


!===============================================================================
! 10.  TEMPS DE SEJOUR
!===============================================================================

if (nor.eq.nordre) then

  do npt = 1,nbpart
    if (ipepa(jisor,npt).ne.0) then
      pepa(jrtsp,npt) = pepa(jrtsp,npt) + dtp
    endif
  enddo

endif

!===============================================================================
! 11.  CALCUL DE L'ADHESION SI MODELE DE REENTRAINEMENT
!===============================================================================

if (ireent.gt.0) then

  call lagres                                                     &
  !==========
 ( parbor, nvisbr)

endif

!===============================================================================
! 11.  CALCUL STATISTIQUES
!===============================================================================

if (nor.eq.nordre .and. istala.eq.1 .and. iplas.ge.idstnt) then

  call lagsta                                                     &
  !==========
 ( nvlsta ,                                                       &
   statis , stativ )

endif

!===============================================================================
! 12.  Equation de Poisson
!===============================================================================

if (nor.eq.nordre .and. ilapoi.eq.1) then
  call lagpoi
  !==========
endif

!===============================================================================
! 13.  Methode de reduction de variances : Clonage/Fusion des particules
!===============================================================================

if ( nor.eq.nordre .and. iroule.ge.1 ) then

  call lagrus(ncelet, ncel, croule)
  !==========

  if (npclon.gt.0) then

    npar1 = nbpart - npclon + 1
    npar2 = nbpart

    ! Use fields at current time step
    iprev = 0

    call lagipn(npar1, npar2, iprev)
    !==========

  endif

endif

!===============================================================================
! 14. UN AUTRE TOUR ?
!===============================================================================

if (nordre.eq.2 .and. nor.eq.1) goto 10

!===============================================================================
! 15. BRANCHEMENT UTILISATEUR POUR MODIF DES VARIABLES EVENTUELLES
!     EN FIN D'ITERATION LAGRANGIENNE
!===============================================================================

call uslast                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   ntersl , nvlsta , nvisbr ,                                     &
   dt     )

!===============================================================================
! 16. Visualisations
!===============================================================================

 20   continue

!===============================================================================
! 17. NOMBRE DE PARTICULES PERDUES (SUITES COMPRISES)
!===============================================================================

nbpper = nbperr
if (irangp .ge. 0) then
  call parcpt(nbpper)
endif
nbpert = nbpert + nbpper

!===============================================================================
! 18. ECRITURE SUR FICHIERS DES INFORMATIONS SUR LE NOMBRE DE PARTICULES
!        - nombre de particules dans le domaine
!        - nombre de particules entrantes
!        - nombre de particules sorties
!        - ...

!===============================================================================

if (ipass.eq.1) then
  modntl = 0
elseif(ntlal.gt.0) then
  modntl = mod(ntcabs,ntlal)
elseif(ntlal.eq.-1.and.ntcabs.eq.ntmabs) then
  modntl = 0
else
  modntl = 1
endif

if (modntl.eq.0) then
  call lagaff
  !==========
endif

! Free memory

deallocate(gradpr)
deallocate(w1, w2, w3)
if (modcpl.gt.0) then
  deallocate(gradvf)
endif
if (iroule.eq.1) then
  deallocate(croule)
endif

if ((idepst.eq.1).and.(ntcabs.eq.ntmabs)) then
   deallocate(vislen)
endif
if (idlvo.eq.1) then
   deallocate(energt)
  if (iclogst.eq.1 .or. irough .eq. 1 ) then
     deallocate(tempp)
  endif
endif

!===============================================================================

!--------
! FORMATS
!--------

 9001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LA TECHNIQUE DE CLONAGE/FUSION DES PARTICULES           ',/,&
'@      EST ENCLENCHEE AVEC UNE FONCTION D''IMPORTANCE        ',/,&
'@      COMPORTANT DES VALEURS NEGATIVES OU NULLES            ',/,&
'@      (LAGUNE).                                             ',/,&
'@                                                            ',/,&
'@    LES ELEMENTS DU TABLEAU CROULE DOIVENT STRICTEMENT      ',/,&
'@      POSITIFS.                                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les valeurs de CROULE dans la subroutine USLARU. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 4100 format(                                                     &
'                                                               '/,&
'   ** LAGRANGIAN MODULE:  '                                     /,&
'   ** deposition submodel  '                                   ,/,&
'      ---------------------------------------------  '         ,/,&
'                                                               '/,&
'                                                               '/,&
'   ** Mean friction velocity  (ustar) =  ',F7.3                ,/,&
'---------------------------------------------------------------  ',/)

!----
! Formats
!----

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'=============================================================',/,&
' Erreur dans le module lagrangien: tentative de terminaison',  /,&
'   ntmabs remis a ', i10,                                      /,&
'=============================================================',/,&
                                                                /)
#else

 1000 format(/,                                                   &
'=============================================================',/,&
' Lagrangian module error: trying to finish cleanly',           /,&
'   ntmabs reset to ', i10,                                   /,&
'=============================================================',/,&
                                                                /)
#endif

!----
! End
!----

end subroutine
