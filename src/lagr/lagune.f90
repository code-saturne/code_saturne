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

subroutine lagune &
!================

 ( lndnod ,                                                       &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   icocel , itycel , ifrlag , itepa  ,                            &
   dlgeo  ,                                                       &
   dt     , rtpa   , rtp    , propce ,                            &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , statis , stativ , tslagr , parbor )

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
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! icocel           ! te ! --> ! connectivite cellules -> faces                 !
!   (lndnod)       !    !     !    face de bord si numero negatif              !
! itycel           ! te ! --> ! connectivite cellules -> faces                 !
!   (ncelet+1)     !    !     !    pointeur du tableau icocel                  !
! ifrlag           ! te ! --> ! numero de zone de la face de bord              !
!   (nfabor)       !    !     !  pour le module lagrangien                     !
! itepa            ! te ! --> ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! dlgeo            ! tr ! --> ! tableau contenant les donnees geometriques     !
! (nfabor,ngeol)   !    !     ! pour le sous-modele de depot                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant et prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp             ! tr ! --> ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettp             ! tr ! --> ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! --> ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! tr ! --> ! moyennes statistiques                          !
!(ncelet,nvlsta    !    !     !                                                !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! tslagr           ! tr ! --> ! terme de couplage retour du                    !
!(ncelet,ntersl    !    !     !   lagrangien sur la phase porteuse             !
! parbor           ! tr ! --> ! infos sur interaction des particules           !
!(nfabor,nvisbr    !    !     !   aux faces de bord                            !
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
use entsor
use cstphy
use cstnum
use parall
use period
use pointe
use lagpar
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
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          icocel(lndnod) , itycel(ncelet+1)
integer          ifrlag(nfabor) , itepa(nbpmax,nivep)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp), ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)
double precision dlgeo(nfabor,ngeol)

! Local variables

integer          ip     , npt    , iok
integer          iel    , ivf
integer          npar1  , npar2
integer          modntl

integer          ifac

double precision visccf, romf
double precision ustarmoy, surftot, surfb

integer, allocatable, dimension(:) :: indep, ibord

double precision, allocatable, dimension(:) :: taup
double precision, allocatable, dimension(:,:) :: tlag, piil
double precision, allocatable, dimension(:,:) :: vagaus
double precision, allocatable, dimension(:,:,:) :: bx
double precision, allocatable, dimension(:,:) :: tsuf, tsup
double precision, allocatable, dimension(:,:) :: tsvar
double precision, allocatable, dimension(:,:) :: tempct
double precision, allocatable, dimension(:) :: tsfext
double precision, allocatable, dimension(:) :: cpgd1, cpgd2, cpght
double precision, allocatable, dimension(:,:) :: brgaus
double precision, allocatable, dimension(:) :: terbru
double precision, allocatable, dimension(:,:) :: gradpr, gradvf
double precision, allocatable, dimension(:) :: croule
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:,:) :: auxl, auxl2

double precision, allocatable, dimension(:,:) :: tslag

double precision, allocatable, save, dimension(:) :: vislen

double precision, allocatable, dimension(:):: energt

double precision, allocatable, dimension(:):: tempp

double precision, dimension(:), pointer :: cromf

integer ii
integer nbpartall, nbpper

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 0.  GESTION MEMOIRE ET COMPTEUR DE PASSAGE
!===============================================================================

! Allocate temporary arrays
allocate(indep(nbpmax), ibord(nbpmax))
allocate(auxl(nbpmax,3))
allocate(taup(nbpmax))
allocate(tlag(nbpmax,3))
allocate(piil(nbpmax,3))
allocate(vagaus(nbpmax,nvgaus))
allocate(tsuf(nbpmax,3))
allocate(tsup(nbpmax,3))
allocate(bx(nbpmax,3,2))
allocate(tsvar(nbpmax,nvp1))
allocate(gradpr(ncelet,3))
allocate(w1(ncelet), w2(ncelet), w3(ncelet))

! Allocate other arrays depending on user options
if ((iphyla.eq.1 .and. itpvar.eq.1) .or. iphyla.eq.2) then
  allocate(tempct(nbpmax,2))
endif
if (iilagr.eq.2) then
  allocate(tsfext(nbpmax))
endif
if (iilagr.eq.2 .and. iphyla.eq.2 .and. ltsthe.eq.1) then
  allocate(cpgd1(nbpmax))
  allocate(cpgd2(nbpmax))
  allocate(cpght(nbpmax))
endif
if (modcpl.gt.0) then
  allocate(gradvf(ncelet,9))
endif
if (iroule.eq.1) then
  allocate(croule(ncelet))
endif
if (lamvbr.eq.1) then
  allocate(brgaus(nbpmax,nbrgau))
  allocate(terbru(nbpmax))
endif
if (nordre.eq.2) then
  allocate(auxl2(nbpmax,7))
endif
if (idlvo.eq.1) then

  allocate(energt(nfabor))
  if (iclogst.eq.1) then
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


!-->Sur Champ fige Lagrangien : RTPA = RTP
!   Rem : cette boucle pourrait etre faite au 1er passage
!         mais la presence de cs_user_extra_operations incite a la prudence...

if (iilagr.eq.3) then
  do ivf = 1,nvar
    do iel = 1,ncel
      rtpa(iel,ivf) = rtp(iel,ivf)
    enddo
  enddo
endif

!-->au premier passage relatif :

if (iplar.eq.1) then

!      Connectivite cellules -> faces + Alloc. structures en C

  call lagbeg                                                     &
  !==========
 ( nbpmax , nlayer , iphyla , idepst , ireent , iclogst  ,        &
   nvls   , nbclst , icocel , itycel ,                            &
   jisor  , jrval  , jrpoi  , jrtsp  , jdp    , jmp    ,          &
   jxp    , jyp    , jzp    , jup    , jvp    , jwp    ,          &
   juf    , jvf    , jwf    , jtaux  , jryplu ,                   &
   jrinpf , jdfac  , jimark ,                                     &
   jtp    , jhp    , jtf    , jmwat  , jmch   , jmck   ,          &
   jcp    , jrdck  , jrd0p  , jinch  , jrhock ,                   &
   jreps  , jdepo  , jnbasg , jnbasp , jfadh  , jmfadh ,          &
   jndisp)

! --> if the deposition model is activated

  if (idepst.ge.1) then

     ustarmoy = 0.d0
     surftot = 0.d0

    ! boundary faces data

     call laggeo(dlgeo)
     !==========

     if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then
       call field_get_val_s(iprpfl(ipproc(irom1)), cromf)
     else
       call field_get_val_s(icrom, cromf)
     endif

     ! Average friction velocity calculation
     do ifac = 1, nfabor

        if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then

           iel = ifabor(ifac)

           surfb = sqrt( surfbo(1,ifac)*surfbo(1,ifac)              &
                      +  surfbo(2,ifac)*surfbo(2,ifac)              &
                      +  surfbo(3,ifac)*surfbo(3,ifac) )

           ! the density pointer according to the flow location

           romf = cromf(iel)
           visccf = propce(iel,ipproc(iviscl)) / romf

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

!===============================================================================
! 1.bis  Initialization for the clogging model
!===============================================================================

if ( iclogst.eq.1 ) then

   do ifac = 1,nfabor
      iel = ifabor(ifac)

      if (iscalt.gt.0) then

         if (iscsth(iscalt).eq.-1) then
            tempp(ifac) = rtp(iel,isca(iscalt)) + tkelvi
         else if (iscsth(iscalt).eq. 1) then
            tempp(ifac) = rtp(iel,isca(iscalt))
         else if (iscsth(iscalt).eq.2) then
            call usthht(1,rtp(iel,isca(iscalt)),tempp(ifac))
         endif

      else
         tempp(ifac) = t0
      endif

   enddo

   call cloginit                                                   &
   !===========
   ( cstfar, epsvid, epseau, fion, jamlim, tempp )

endif


!===============================================================================
! 2.  MISE A JOUR DES NOUVELLES PARTICULES ENTREES DANS LE DOMAINE
!===============================================================================

! Au premier pas de temps on initalise les particules avec RTP et
! non RTPA car RTPA = initialisation

if ( ntcabs.eq.1 ) then

  call lagent                                                     &
  !==========
 ( lndnod ,                                                       &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itycel , icocel , dlgeo  ,                                     &
   itypfb , itrifb , ifrlag , itepa  ,                            &
   dt     , rtp    , propce ,                                     &
   ettp   , tepa   , vagaus )

else

  call lagent                                                     &
  !==========
 ( lndnod ,                                                       &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itycel , icocel , dlgeo  ,                                     &
   itypfb , itrifb , ifrlag , itepa  ,                            &
   dt     , rtpa   , propce ,                                     &
   ettp   , tepa   , vagaus )
endif

!===============================================================================
! 2.1 CALCUL DE LA FONCTION D'IMPORTANCE POUR LA ROULETTE RUSSE
!===============================================================================

if (iroule.ge.1) then

  call uslaru                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa ,                                      &
   dt     , rtpa   , propce ,                                     &
   ettp   , tepa   , vagaus , croule , auxl ,                     &
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

!
if (nbpartall.eq.0) goto 20

!-->On enregistre l'element de depart de la particule

do ip = 1,nbpart
  indep(ip) = itepa(ip,jisor)
enddo

!===============================================================================
! 4.  GRADIENT DE PRESSION ET DE LA VITESSE FLUIDE
!===============================================================================

! Au premier pas de temps on calcul les gradient avec RTP et
! non RTPA car RTPA = initialisation (gradients nuls)

if ( ntcabs.eq.1 ) then

  call laggra(rtp, coefa, coefb, gradpr, gradvf)
  !==========

else

  call laggra(rtpa, coefa, coefb, gradpr, gradvf)
  !==========

endif

!===============================================================================
! 4.  Initialisation des variables aleatoires gaussiennes
!===============================================================================

!---> CALCUL DES TIRAGES ALEATOIRES
!     remarque : NORMALEN est dans le fichier ZUFALL.F
!     ^^^^^^^^

if (idistu.eq.1) then
  do ivf = 1,nvgaus
    call normalen(nbpart, vagaus(1,ivf))
  enddo
else
  do ivf = 1,nvgaus
    do ip = 1,nbpmax
      vagaus(ip,ivf) = 0.d0
    enddo
  enddo
endif

!---> CALCUL DES TIRAGES ALEATOIRES POUR LE MVT BROWNIEN

if ( lamvbr .eq. 1 ) then

  do ivf = 1,nbrgau
    call normalen(nbpart, brgaus(1,ivf))
  enddo

endif

!===============================================================================
! 5. PROGRESSION DES PARTICULES
!===============================================================================

 10   continue

nor = mod(nor,nordre)
nor = nor + 1

!---> Recopie des resultats de l'etape precedente :

if (nor.eq.1) then

  do ivf = 1,nvp
    do ip = 1,nbpart
      ettpa(ip,ivf) = ettp(ip,ivf)
    enddo
  enddo

endif

!-----> CALCUL GRADIENT DE PRESSION ET DE LA VITESSE FLUIDE
!       EN N+1 (avec RTP)

if (nor.eq.2 .and. iilagr.ne.3) then

  call laggra(rtp, coefa, coefb, gradpr, gradvf)
  !==========

endif

!-----> CALCUL DES CARACTERISTIQUES DES PARTICULES

if (nor.eq.1) then

!      sous pas de temps n (avec RTPA)

  call lagcar                                                     &
  !==========
   ( nvar   , nscal  ,                                            &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     nvlsta ,                                                     &
     itepa  ,                                                     &
     dt     , rtpa   , propce ,                                   &
     ettp   , ettpa  , tepa   , taup   , tlag   ,                 &
     piil   , bx     , tempct , statis ,                          &
     gradpr , gradvf , w1     , w2     , auxl(1,1) )

else

!     sous pas de temps n+1 (avec RTP)

  call lagcar                                                     &
  !==========
   ( nvar   , nscal  ,                                            &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     nvlsta ,                                                     &
     itepa  ,                                                     &
     dt     , rtp    , propce ,                                   &
     ettp   , ettpa  , tepa   , taup   , tlag   ,                 &
     piil   , bx     , tempct , statis ,                          &
     gradpr , gradvf , w1     , w2     , auxl(1,1) )

endif


!---> INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES
!     POSITION, VITESSE FLUIDE, VITESSE PARTICULE


call lagesp                                                       &
!==========
   ( nvar   , nscal  ,                                            &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     ntersl , nvlsta , nvisbr ,                                   &
     itepa  , ibord  ,                                            &
     dlgeo  ,                                                     &
     dt     , rtpa   , rtp    , propce ,                          &
     ettp   , ettpa  , tepa   ,                                   &
     statis , stativ , taup   , tlag   , piil   ,                 &
     tsuf   , tsup   , bx     , tsfext ,                          &
     vagaus , gradpr , gradvf , brgaus , terbru ,                 &
     auxl(1,1) , auxl2        , vislen)

!---> INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES
!     LIEES AUX PHYSIQUES PARTICULIERES PARTICULAIRES

if ( iphyla.eq.1 .or. iphyla.eq.2 ) then

  if ( nor.eq.1 ) then
    call lagphy                                                   &
    !==========
    ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                &
      ntersl , nvlsta , nvisbr ,                                  &
      itepa  , ibord  ,                                           &
      dt     , rtpa   , propce ,                                  &
      ettp   , ettpa  , tepa   , taup   , tlag   , tempct ,       &
      tsvar  , auxl   , cpgd1  , cpgd2  , cpght  )
  else
    call lagphy                                                   &
    !==========
    ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                &
      ntersl , nvlsta , nvisbr ,                                  &
      itepa  , ibord  ,                                           &
      dt     , rtp    , propce ,                                  &
      ettp   , ettpa  , tepa   , taup   , tlag   , tempct ,       &
      tsvar  , auxl   , cpgd1  , cpgd2  , cpght  )
  endif

endif

!===============================================================================
! 6.  Couplage Retour - Calcul des termes sources
!===============================================================================

if (iilagr.eq.2 .and. nor.eq.nordre) then

  ! Allocate a temporary array
  allocate(tslag(nbpmax,ntersl))

  call lagcou                                                     &
  !==========
   ( nbpmax , nvp    , nvep   , nivep  ,                          &
     ntersl ,                                                     &
     itepa  , indep  , ibord  ,                                   &
     rtp    , propce ,                                            &
     ettp   , ettpa  , tepa   , taup   ,                          &
     tempct , tsfext , tslagr ,                                   &
     cpgd1  , cpgd2  , cpght  ,                                   &
     tslag  , w1     , w2   ,                                     &
     auxl(1,1) , auxl(1,2)   , auxl(1,3) )

     ! Free memory
     deallocate(tslag)

endif

!===============================================================================
! 7.  Calcul de la barrière d'énergie dans le cas DLVO
!===============================================================================

if (idlvo.eq.1) then

   call lagbar                                                    &
   !==========
  ( rtp , energt )

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
   iusclb , iusmoy , deblag , ifrlag )


  call prtget                                                     &
  !==========
 ( nbpmax , nbpart ,                                              &
   ettp   , ettpa  , itepa  , tepa   ,                            &
   ibord  , indep  )

  call dplprt                                                     &
  !==========
 ( nbpart   , dnbpar   , nordre   , parbor   , iensi3   ,         &
   inbr     , inbrbd   , iflm     , iflmbd   , iang     ,         &
   iangbd   , ivit     , ivitbd   , iencnb   , iencma   ,         &
   iencdi   , iencck   , iencnbbd , iencmabd , iencdibd ,         &
   iencckbd , inclg    , iscovc   ,                               &
   nusbor   , iusb     , vislen   , dlgeo,                        &
   rtp      , iu       , iv       , iw       , energt   ,         &
   tprenc   , visref   , enc1     , enc2     , tkelvi)

  call prtput                                                     &
  !==========
 ( nbpmax , nbpart , dnbpar , nbpout , dnbpou , nbperr , dnbper,  &
   nbpdep , dnbdep , npencr , dnpenc ,                            &
   ettp   , ettpa  , itepa  , tepa   ,                            &
   ibord  )


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
    if ( itepa(npt,jisor).ne.0 ) then
      tepa(npt,jrtsp) = tepa(npt,jrtsp) + dtp
    endif
  enddo

endif

!===============================================================================
! 11.  CALCUL DE L'ADHESION SI MODELE DE REENTRAINEMENT
!===============================================================================

if (ireent.gt.0) then

  call lagres                                                     &
  !==========
 ( nbpmax , nvp    , nvep   , nivep  ,                            &
   itepa  ,                                                       &
   ettp   , ettpa  , tepa )

endif

!===============================================================================
! 11.  CALCUL STATISTIQUES
!===============================================================================

if (nor.eq.nordre .and. istala.eq.1 .and. iplas.ge.idstnt) then

  call lagsta                                                     &
  !==========
 ( nbpmax , nvp    , nvep   , nivep  ,                            &
   nvlsta ,                                                       &
   itepa  ,                                                       &
   ettp   , tepa   , statis , stativ )

endif

!===============================================================================
! 12.  Equation de Poisson
!===============================================================================

if (nor.eq.nordre .and. ilapoi.eq.1) then

  call lagpoi                                                     &
  !==========
 ( nbpmax , nvp    , nivep  ,                                     &
   nvlsta ,                                                       &
   itepa  ,                                                       &
   ettp   , statis )

endif

!===============================================================================
! 13.  Methode de reduction de variances : Clonage/Fusion des particules
!===============================================================================

if ( nor.eq.nordre .and. iroule.ge.1 ) then

  call lagrus                                                     &
  !==========
   ( ncelet , ncel   ,                                            &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     itepa  , indep  ,                                            &
     ettp   , ettpa  , tepa   , croule )

  if (npclon.gt.0) then

    npar1 = nbpart - npclon + 1
    npar2 = nbpart

    call lagipn                                                   &
    !==========
    ( nbpmax , nvp    , nvep   , nivep  ,                         &
      npar1  , npar2  ,                                           &
      itepa  ,                                                    &
      rtp    ,                                                    &
      ettp   , tepa   , vagaus ,                                  &
      icocel , lndnod , itycel , dlgeo  , propce , ifrlag )

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
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   dt     , rtpa   , rtp    , propce ,                            &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct ,          &
   statis , stativ )

!===============================================================================
! 16. Visualisations
!===============================================================================

 20   continue

call ucdprt                                                        &
!==========
 ( nbpmax , nbpart , dnbpar , nbpout , dnbpou , nbperr ,           &
   dnbper , nbpdep , dnbdep , npencr , dnpenc ,                    &
   ettp   , ettpa  , itepa  , tepa   , ibord  , indep  )

!===============================================================================
! 17. NOMBRE DE PARITICULES PERDUES (SUITES COMPRISES)
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

deallocate(indep, ibord)
deallocate(auxl)
deallocate(taup)
deallocate(tlag)
deallocate(piil)
deallocate(vagaus)
deallocate(tsuf)
deallocate(tsup)
deallocate(bx)
deallocate(tsvar)
deallocate(gradpr)
deallocate(w1, w2, w3)
if ((iphyla.eq.1 .and. itpvar.eq.1) .or. iphyla.eq.2) then
  deallocate(tempct)
endif
if (iilagr.eq.2) then
  deallocate(tsfext)
endif
if (iilagr.eq.2 .and. iphyla.eq.2 .and. ltsthe.eq.1) then
  deallocate(cpgd1)
  deallocate(cpgd2)
  deallocate(cpght)
endif
if (modcpl.gt.0) then
  deallocate(gradvf)
endif
if (iroule.eq.1) then
  deallocate(croule)
endif
if (lamvbr.eq.1) then
  deallocate(brgaus)
  deallocate(terbru)
endif
if (nordre.eq.2) then
  deallocate(auxl2)
endif

if ((idepst.eq.1).and.(ntcabs.eq.ntmabs)) then
   deallocate(vislen)
endif
if (idlvo.eq.1) then
   deallocate(energt)
  if (iclogst.eq.1) then
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
