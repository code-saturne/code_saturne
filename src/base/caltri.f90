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

subroutine caltri &
!================

 ( iverif ,                                                       &
   ia     ,                                                       &
   ra     )

!===============================================================================
! Purpose:
! --------

! Main solver subroutine (read, time loop, write)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iverif           ! i  ! <-- ! elementary tests flag                          !
! ia(*)            ! ia ! --- ! main integer work array                        !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________.____._____.________________________________________________.

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use pointe
use optcal
use cstphy
use entsor
use albase
use parall
use period
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use lagpar
use lagdim
use lagran
use vorinc
use ihmpre
use radiat
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments


integer          iverif
integer          ia(*)

double precision ra(*)

! Local variables

integer          ipropc , ipropf , ipropb
integer          icoefa , icoefb
integer          irtp   , irtpa
integer          idt
integer          iisstd , ifrcx

integer          idebia , idebra
integer          ifinia , ifinra , idbia1 , idbra1, idbia2
integer          iditia, iditra
integer          ifnia1 , ifnra1 , ifnia2 , ifnia3, ifnra2
integer          iiii

integer          modhis, iappel, modntl, iisuit, iwarn0
integer          ntsdef, ntcam1
integer          iphas , ivar

integer          iicoce , iityce
integer          iiitep , iitepa , istatc , istatf
integer          iettp  , iettpa , iauxl  , itslag , istatv
integer          itaup  , iitlag , ipiil  , iindep , iibord
integer          ivagau , itsuf  , itsup  , ibx    , iauxl2
integer          ibrgau , itebru
integer          igradp , igradv , icroul
integer          itepct , itsfex , itsvar
integer          icpgd1 , icpgd2 , icpght
integer          ilagia , ilagra , iiwork
integer          iw1    , iw2    , iw3
integer          inod   , idim
integer          itrale , indact , indwri
integer          maxelt , ils

double precision titer1, titer2
double precision tecrf1, tecrf2

double precision, save :: ttchis

character        ficsui*32


!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Initialize first free position in arrays "ia" and "ra"
idebia = 1
idebra = 1

! Initialize random number generator
! (not always necessary, but not at all costly)
call zufalli(0)
!===========

!---> Stop test set to 1 if P-1 radiative module "sees" too many cells
!     with an optical thickness greater than 1 (see ppcabs).
istpp1 = 0

!---> Maximum number of elements for selector.
maxelt = max(ncelet, nfac, nfabor)

!--> Probes output tracking
ttchis = -1.d0

!===============================================================================
! 2. Geometry
!===============================================================================

! (memclg directly initializes pointe.f90)
call memclg                                                       &
!==========
 ( idebia , idebra ,                                              &
   ifinia , ifinra )

!---> Memory will be maintained until the end.
idebia = ifinia
idebra = ifinra

call cregeo (idebia, idebra, ia, ra)
!==========

!---> Memory will be maintained until the end.
idebia = ifinia
idebra = ifinra

!===============================================================================
! 3. End of modules initialization
!===============================================================================

call initi2(idebia, idebra,  ia, ra)
!==========

if (iilagr.gt.0) then

  !--> Compute "lndnod" (lagran.f90)

  ! Integer work array of size ncelet
  iiwork = idebia
  ifinia = iiwork + ncelet
  call iasize ('lagini', ifinia)

  ifinra = idebra

  call lagini                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ncelet , ncel   , nfac   , nfabor ,                            &
   lndnod ,                                                       &
   ifacel , ifabor ,                                              &
   ia(iiwork) ,                                                   &
   ia     , ra     )

endif

!===============================================================================
! 4. Other arrays
!===============================================================================

!---> Memory management

call memtri                                                       &
!==========
 ( idebia , idebra , iverif ,                                     &
   nvar   , nscal  , nphas  ,                                     &
   ncofab , nproce , nprofa , nprofb ,                            &
   iisstd , ifrcx  ,                                              &
   idt    , irtp   , irtpa  , ipropc , ipropf , ipropb ,          &
   icoefa , icoefb ,                                              &
   ifinia , ifinra )

! Memory reservation for additional arrays required by specific physics.

idbia1 = ifinia
idbra1 = ifinra

call memppt                                                       &
!==========
 ( idbia1 , idbra1 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ifinia , ifinra )

!===============================================================================
! 4.1 Memory reservation for semi-transparent radiation module
!===============================================================================

if (iirayo.gt.0) then

  idbia1 = ifinia
  idbra1 = ifinra

  call memra1                                                     &
  !==========
 ( idbia1 , idbra1 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ifinia , ifinra )

endif

!===============================================================================
! 4.2 Memory reservation for Lagrangian module
!===============================================================================

!     Si on ne fait pas de Lagrangien, on initialise
!       quand meme les "pointeurs".

idbia1 = ifinia
idbra1 = ifinra

call memla1                                                       &
!==========
  ( idbia1 , idbra1 ,                                             &
    lndnod ,                                                      &
    nbpmax , nvp    , nvp1   , nvep   , nivep  ,                  &
    ntersl , nvlsta , nvisbr ,                                    &
    iiitep , iicoce , iityce ,                                    &
    iettp  , iettpa , iitepa , istatc , istatv, itslag , istatf , &
    ifinia , ifinra )

!===============================================================================
! 4.3 TESTS ELEMENTAIRES : appel a testel.f90
!===============================================================================

if (iverif.eq.1) then

  write(nfecra, 1000)

  call testel                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   nphas  , nvar   ,                                              &
   ia     ,                                                       &
   ra(irtp) , ra(icoefa) , ra(icoefb) ,                           &
   ra     )

  goto 200

endif

!===============================================================================
! 5. INITIALISATIONS PAR DEFAUT
!===============================================================================

call iniva0                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  , nphas  , ncofab ,                            &
   ia     ,                                                       &
   ra(idt)    , ra(irtp)   , ra(ipropc) , ra(ipropf) , ra(ipropb),&
   ra(icoefa) , ra(icoefb) , ra(ifrcx ) ,                         &
   ra     )

!===============================================================================
! 6. CALCUL SUITE EVENTUEL
!===============================================================================

if (isuite.eq.1) then

  call lecamo                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nnod   ,          &
   nvar   , nscal  , nphas  ,                                     &
   ia     ,                                                       &
   ra(idt)    , ra(irtp)   , ra(ipropc) , ra(ipropf) , ra(ipropb),&
   ra(icoefa) , ra(icoefb) , ra(ifrcx ) ,                         &
   ra     )

  ! Using ALE, geometric parameters must be recalculated
  if (iale.eq.1) then

    do inod = 1, nnod
      do idim = 1, ndim
        xyznod(idim,inod) =   ra(ixyzn0+(inod-1)*ndim+idim-1)     &
                            + ra(idepal+(idim-1)*nnod+inod-1)
      enddo
    enddo

    call algrma
    !==========
    call calgeo                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ia     ,                                                       &
   volmin , volmax , voltot ,                                     &
   ra     )

  endif

endif

!===============================================================================
! 7. INITIALISATIONS (Utilisateur et complementaires)
!    RTP DT ROM ROMB VISCL VISCT VISCLS
!    (TPUCOU en PERIODICITE)
!===============================================================================

call inivar                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  , nphas  , ncofab ,                            &
   ia     ,                                                       &
   ra(idt)    , ra(irtp)   , ra(ipropc) , ra(ipropf) , ra(ipropb),&
   ra(icoefa) , ra(icoefb) , ra(ifrcx ) ,                         &
   ra     )

!===============================================================================
! 8.1 MODULE DE RAYONNEMENT : CALCUL SUITE EVENTUEL
!===============================================================================

if (iirayo.gt.0 .and. isuird.eq.1) then

  call raylec                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor ,                   &
   ia     ,                                                       &
   ra(ipropc) , ra(ipropb) ,                                      &
   ra     )

endif

!===============================================================================
! 8.2 INITIALISATIONS DES PARTICULES POUR LE LAGRANGIEN
!===============================================================================

if (iilagr.gt.0) then

  call laglec                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor ,                   &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   ia(iiitep) , ia ,                                              &
   ra(irtpa)  , ra(ipropc) ,                                      &
   ra(iettp)  , ra(iitepa) , ra(istatc) , ra(istatv) ,            &
   ra(istatf) , ra(itslag) , ra         )

endif

!===============================================================================
! 8.3 INITIALISATIONS POUR LE MODULE THERMIQUE 1D EN PAROI
!===============================================================================
! On suppose que toutes les phases voient la meme temperature de paroi
! USPT1D a un fonctionnement similaire a USKPDC et USTSMA, mais comme
! on ecrit des infos dans un fichier suite, on a besoin d'une partie de
! la memoire meme apres la boucle en temps -> IFPT1D et TPPT1D
!                                            (IFNIA1 et IFNRA1)

idbia1 = ifinia
idbra1 = ifinra

ils    = idbia1
ifnia2 = ils + maxelt
call iasize('caltri',ifnia2)

iphas = 1

!     Premier appel : definition de NFPT1D et ISUIT1
iappel = 1
call uspt1d                                                       &
!==========
 ( ifnia2 , idbra1 ,                                              &
   nvar   , nscal  , nphas  , nfpt1d , iphas  , iappel ,          &
   maxelt , ia(ils),                                              &
   ia(idbia1) , ia(idbia1) , ia(idbia1) ,                         &
   ia     ,                                                       &
   ra(idbra1) , ra(idbra1) , ra(idbra1) ,                         &
   ra(idbra1) , ra(idbra1) , ra(idbra1) ,                         &
   ra(idbra1) , ra(idbra1) , ra(idbra1) ,                         &
   ra(idt)    , ra(irtpa)  ,                                      &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra     )

iappel = 1
call vert1d                                                       &
!==========
 (idbia1     , idbra1     ,                                       &
  nfabor     , nfpt1d     , iappel     ,                          &
  ia(idbia1) , ia(idbia1) , ia(idbia1) , ia     ,                 &
  ra(idbra1) , ra(idbra1) ,                                       &
  ra(idbra1) , ra(idbra1) , ra(idbra1) , ra     )

call memt1d                                                       &
!==========
 ( idbia1 , idbra1 , nfabor , ifnia1 , ifnra1 ,ifnia2 , ifnra2 ,  &
   ifinia , ifinra , ia     , ra     )

! On appelle uspt1d lorqu'il y a sur un processeur au moins des faces de
!     bord avec module thermique 1D.

if (nfpt1t.gt.0) then
! Deuxieme appel : remplissage des tableaux de definition de la geometrie
!            et de l'initialisation (IFPT1D,NPPT1D,EPPT1D,RGPT1D,TPPT1D)
  ils    = ifinia
  ifnia3 = ils + maxelt
  call iasize('caltri',ifnia3)

  iappel = 2
  call  uspt1d                                                    &
  !===========
 ( ifnia3 , ifinra ,                                              &
   nvar   , nscal  , nphas  , nfpt1d , iphas  , iappel ,          &
   maxelt , ia(ils),                                              &
   ia(iifpt1) , ia(inppt1) , ia(iiclt1) ,                         &
   ia     ,                                                       &
   ra(itppt1) , ra(irgpt1) , ra(ieppt1) ,                         &
   ra(itept1) , ra(ihept1) , ra(ifept1) ,                         &
   ra(ixlmt1) , ra(ircpt1) , ra(idtpt1) ,                         &
   ra(idt)    , ra(irtpa)  ,                                      &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra     )

  iappel = 2
  call vert1d                                                     &
  !==========
 (ifinia     , ifinra     ,                                       &
  nfabor     , nfpt1d     , iappel     ,                          &
  ia(iifpt1) , ia(inppt1) , ia(iiclt1) , ia     ,                 &
  ra(irgpt1) , ra(ieppt1) ,                                       &
  ra(ixlmt1) , ra(ircpt1) , ra(idtpt1) , ra     )

!     Calcul du max des NPPT1D (pour les fichiers suite)
  nmxt1d = 0
  do iiii = 1, nfpt1d
    nmxt1d = max(ia(inppt1+iiii-1),nmxt1d)
  enddo
  if (irangp.ge.0) call parcmx(nmxt1d)
                   !==========

  if (isuit1.eq.1) then

    ficsui = '1dwall_module'
    call lect1d                                                   &
    !==========
 ( ficsui     , len(ficsui), nfpt1d     , nfpt1t    ,             &
   nmxt1d     , nfabor     , ia(inppt1) , ia(iifpt1) , ra(ieppt1),&
   ra(irgpt1) , ra(itppt1))

  else
!     Creation du maillage, initialisation de la temperature.

    call mait1d                                                   &
    !==========
 ( nfpt1d, ia(inppt1), ra(ieppt1), ra(irgpt1),ra(itppt1))

  endif

endif
!     Les infos sur l'epaisseur de la paroi, le nombre de points de
!     discretisation et la raison geometrique ont ete transmises a
!     la structure C. Elles sont maintenant inutiles dans le Fortran.
!     -> on libere la memoire.
ifinia = ifnia2
ifinra = ifnra2

!===============================================================================
! 9. TABLEAUX POUR BLC EN TEMPS MAIS A OUBLIER ENSUITE
!===============================================================================

!  En fin de bloc en temps on doit retrouver IFNIA1 et IFNRA1
iditia = ifnia1
iditra = ifnra1

idbia1 = ifinia
idbra1 = ifinra


do iphas = 1, nphas

  iappel = 1

  if (iihmpr.eq.1) then
    call uikpdc &
    !==========
  ( iappel, iphas, ncelet, ncepdc,     &
    ia(idbia1), ra(idbra1) , ra(irtpa) )
  endif

  ils    = idbia1
  idbia2 = ils + maxelt
  call iasize('caltri',idbia2)

  call  uskpdc &
  !===========
( idbia2 , idbra1 ,                                              &
  nvar   , nscal  , nphas  ,                                     &
  ncepdc(iphas) , iphas  , iappel ,                              &
  maxelt , ia(ils),                                              &
  ia(idbia1),                                                    &
  ia     ,                                                       &
  ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
  ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
  ra(icoefa) , ra(icoefb) , ra(idbra1) ,                         &
  ra     )

enddo

call mempdc                                                       &
!==========
 ( idbia1, idbra1, ncelet, ncel,  nphas, ndim, ifinia, ifinra)


! On appelle uskpdc lorqu'il y a sur un processeur au moins des cellules
!     avec terme source de masse.
!     On ne fait que remplir le tableau d'indirection des cellules
!     On appelle cependant uskpdc avec tous les processeurs, au cas ou
!     l'utilisateur aurait mis en oeuvre des operations globales.

do iphas = 1, nphas

  if(ncpdct(iphas).gt.0) then

    iappel = 2

    if (iihmpr.eq.1) then
      call uikpdc &
      !==========
    ( iappel, iphas, ncelet, ncepdc,                  &
      ia(iicepd(iphas)), ra(ickupd(iphas)), ra(irtpa) )
    endif

    ils    = ifinia
    ifnia2 = ils + maxelt
    call iasize('caltri',ifnia2)

    call  uskpdc                                                &
    !===========
  ( ifnia2 , ifinra ,                                              &
    nvar   , nscal  , nphas  ,                                     &
    ncepdc(iphas) , iphas  , iappel ,                              &
    maxelt , ia(ils),                                              &
    ia(iicepd(iphas)),                                             &
    ia     ,                                                       &
    ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
    ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
    ra(icoefa) , ra(icoefb) , ra(ickupd(iphas)) ,                  &
    ra     )

  endif

enddo

idbia1 = ifinia
idbra1 = ifinra

ils    = idbia1
idbia2 = ils + maxelt
call iasize('caltri',idbia2)

do iphas = 1, nphas

  iappel = 1
  call  ustsma                                                    &
  !===========
 ( idbia2 , idbra1 ,                                              &
   nvar   , nscal  , nphas  , ncepdc(iphas)   ,                   &
   ncetsm(iphas) ,   iphas  , iappel ,                            &
   maxelt , ia(ils),                                              &
   ia(iicepd(iphas)) ,                                            &
   ia(idbia1) , ia(idbia1),                                       &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  ,                                      &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) , ra(ickupd(iphas))       , ra(idbra1),&
   ra     )

enddo

call memtsm                                                       &
!==========
     ( idbia1 , idbra1 ,                                          &
       ncelet , ncel   , nvar   , nphas  ,                        &
       ifinia , ifinra )

! On appelle ustsma lorqu'il y a sur un processeur au moins des cellules
!     avec terme source de masse.
!     On ne fait que remplir le tableau d'indirection des cellules
!     On appelle cependant ustsma avec tous les processeurs, au cas ou
!     l'utilisateur aurait mis en oeuvre des operations globales.

do iphas = 1, nphas

  if(nctsmt(iphas).gt.0) then

    ils    = ifinia
    ifnia2 = ils + maxelt
    call iasize('caltri',ifnia2)

    iappel = 2
    call  ustsma                                                  &
    !===========
 ( ifnia2 , ifinra ,                                              &
   nvar   , nscal  , nphas  , ncepdc(iphas)   ,                   &
   ncetsm(iphas) ,   iphas  , iappel ,                            &
   maxelt , ia(ils),                                              &
   ia(iicepd(iphas)) ,                                            &
   ia(iicesm(iphas)) , ia(iitpsm(iphas)),                         &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  ,                                      &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) , ra(ickupd(iphas)), ra(ismace(iphas)),&
   ra     )

  endif

enddo


! -- Methode des vortex pour la L.E.S.
!    (dans verini on s'est deja assure que ITYTUR=4 si IVRTEX=1)

if (ivrtex.eq.1) then

  idbia1 = ifinia
  idbra1 = ifinra

  iphas  = 1
  iappel = 1

!  On met une valeur factice a certains parametres non utilise en IAPPEL=1

  call memvor(idbia1, idbra1, iappel, nfabor, ifinia, ifinra)
  !==========

  call vorin0(nfabor, ia(iirepv))
  !==========

  ils    = ifinia
  ifnia2 = ils + maxelt
  call iasize('caltri',ifnia2)

  call usvort                                                     &
  !==========
 ( ifnia2 , ifinra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iphas  , iappel ,                                              &
   maxelt , ia(ils),                                              &
   ia(iirepv)      ,                                              &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  ,                                      &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra     )

  call vorver ( nfabor , ia(iirepv)  , iappel )
  !==========

  idbia1 = ifinia
  idbra1 = ifinra

! Attention, vorpre reserve de la memoire qu'il faut garder ensuite
!           (-> on utilise IFINIA/IFINRA ensuite)

  call vorpre                                                     &
  !==========
 ( idbia1 , idbra1 , ifinia , ifinra ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ia(iirepv),                                                    &
   ia     ,                                                       &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra     )

endif

! -- Fin de zone Methode des vortex pour la L.E.S.

! -- Structures mobiles en ALE

if (iale.eq.1) then

  idbia1 = ifinia
  idbra1 = ifinra

! Attention, strini reserve de la memoire qu'il faut garder ensuite
!           (-> on utilise IFINIA/IFINRA ensuite)
  call strini                                                     &
  !==========
 ( idbia1 , idbra1 , ifinia , ifinra ,                            &
   ia     ,                                                       &
   ra(idt),                                                       &
   ra     )

endif

! -- Fin de zone Structures mobiles en ALE

! -- Couplage Code_Saturne/Code_Saturne

call cscini                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ia     ,                                                       &
   ra     )


!===============================================================================
! 10. DEBUT DE LA BOUCLE EN TEMPS
!===============================================================================

write(nfecra,2000)

ntcabs = ntpabs
ttcabs = ttpabs

iwarn0 = 1
do ivar = 1, nvar
  iwarn0 = max(iwarn0,iwarni(ivar))
enddo

if(iwarn0.gt.0) then
  write(nfecra,3000)
endif

if(inpdt0.eq.1) then
  ntmabs = ntcabs
endif

!     Nb d'iter ALE (nb relatif a l'execution en cours)
!     Si ITALIN=1, on fait une iteration d'initialisation
!     (si ITALIN=-999, c'est qu'on a fait une suite de calculs
!      sans relire lecamx -> comme si on faisait une suite
!      d'un calcul sans ALE)
if (italin.eq.-999) italin = 1
itrale = 1
if (italin.eq.1) then
  itrale = 0
  write(nfecra,3002) ttcabs
endif

!     En cas de couplage avec SYRTHES, on lit des maintenant l'entete
!     du premier message, au cas ou il s'agisse d'un message de
!     terminaison.

if (itrale.gt.0) then

  ntcam1 = ntcabs - 1

  call tstsyr (ntmabs, ntcam1)
  !==========

  if (ntmabs .eq. ntcam1) then
    call csexit (0)
    !==========
  endif

endif

 100  continue

if (inpdt0.eq.0 .and. itrale.gt.0) then
  ntcabs = ntcabs + 1
  if(idtvar.eq.0.or.idtvar.eq.1) then
    ttcabs = ttcabs + ra(idt)
  else
    ttcabs = ttcabs + dtref
  endif
  if(iwarn0.gt.0) then
    write(nfecra,3001) ttcabs,ntcabs
  endif
endif


!===============================================================================
! 11. AVANCEE EN TEMPS
!===============================================================================


!     On teste la presence de ficstp pour modifier NTMABS le cas echeant
call modpar(ntcabs,ntmabs)
!==========

call dmtmps(titer1)
!==========

!     Synchronisation Syrthes 3, si ITRALE>0
if (itrale.gt.0) then
  call itdsyr(ntcabs,ntmabs)
  !==========
endif

call tridim                                                       &
!==========
 ( ifinia , ifinra , itrale ,                                     &
   nvar   , nscal  , nphas  ,                                     &
   ia(iisstd),                                                    &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(itslag) , ra(icoefa) , ra(icoefb) ,                         &
   ra(ifrcx)  ,                                                   &
   ra     )

!===============================================================================
! 12. CALCUL DES MOYENNES TEMPORELLES (ACCUMULATION)
!===============================================================================


if(inpdt0.eq.0 .and. itrale.gt.0) then
call calmom                                                       &
!==========
     ( ifinia , ifinra , ncel   , ncelet ,                        &
       ia     ,                                                   &
       ra(irtp  ) , ra(idt   ) , ra(ipropc) ,                     &
       ra     )
endif


!===============================================================================
! 13. APPEL DU MODULE LAGRANGIEN
!===============================================================================

if (iilagr.gt.0 .and. inpdt0.eq.0 .and. itrale.gt.0) then

  call memla2                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   iindep , iibord , iettpa , iauxl  , iauxl2 ,                   &
   itaup  , iitlag , ipiil  ,                                     &
   ivagau , itsuf  , itsup  , ibx    ,                            &
   igradp , igradv , icroul ,                                     &
   itepct , itsfex , itsvar ,                                     &
   icpgd1 , icpgd2 , icpght ,                                     &
   ibrgau , itebru ,                                              &
   iw1    , iw2    , iw3    ,                                     &
   ilagia , ilagra  )

  call lagune                                                     &
  !==========
 ( ilagia , ilagra ,                                              &
   lndnod ,                                                       &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   ia(iicoce) , ia(iityce) , ia(iifrla) , ia(iiitep) ,            &
   ia(iindep) , ia(iibord) ,                                      &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(iettp)  , ra(iettpa) , ra(iitepa) , ra(istatc) , ra(istatv),&
   ra(itslag),                                                    &
   ra(istatf) , ra(itaup)  , ra(iitlag) , ra(ipiil)  , ra(ibx  ) ,&
   ra(ivagau) , ra(itsuf ) , ra(itsup ) , ra(itsvar) ,            &
   ra(itepct) , ra(itsfex) ,                                      &
   ra(icpgd1) , ra(icpgd2) , ra(icpght) ,                         &
   ra(igradp) , ra(igradv) , ra(icroul) ,                         &
   ra(ibrgau) , ra(itebru) ,                                      &
   ra(iw1   ) , ra(iw2   ) , ra(iw3   ) , ra(iauxl)  , ra(iauxl2),&
   ra     )

!--> Ici on libere la memoire reserve par MEMLA2
!      (i.e. on oublie ILAGIA et ILAGRA)

endif

!===============================================================================
! 14. BRANCHEMENT UTILISATEUR POUR MODIF DES VARIABLES EVENTUELLES
!===============================================================================

if (itrale.gt.0) then

!       Sortie postprocessing de profils 1D

  if (iihmpr.eq.1) then
    call uiprof                                                   &
    !==========
  ( ncelet , ncel,                                                &
    ntmabs, ntcabs, ttcabs,                                       &
    xyzcen, ra(irtp), ra(ipropc) )
  endif

  ils    = ifinia
  ifnia2 = ils + maxelt
  call iasize('caltri',ifnia2)

  call usproj                                                     &
  !==========
 ( ifnia2 , ifinra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvep   , nivep  , ntersl , nvlsta , nvisbr , &
   maxelt , ia(ils),                                              &
   ia(iiitep),                                                    &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(iettp)  , ra(iettpa) , ra(iitepa) , ra(istatc) , ra(istatv),&
   ra(itslag) , ra(istatf) ,                                      &
   ra     )

endif

!===============================================================================
! 15. MISE A JOUR DU MAILLAGE (ALE)
!===============================================================================

if (iale.eq.1 .and. inpdt0.eq.0) then

  if (itrale.eq.0 .or. itrale.gt.nalinf) then

    call alemaj                                                   &
    !==========
 ( ifinia , ifinra , itrale ,                                     &
   nvar   , nscal  , nphas  ,                                     &
   ia(iimpal)      ,                                              &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(idepal) , ra(ixyzn0) ,                                      &
   ra     )

  endif

endif

!===============================================================================
! 16. TEST D'ARRET PAR MANQUE DE TEMPS
!===============================================================================

call armtps(ntcabs,ntmabs)
!==========

!===============================================================================
! 17. TEST D'ARRET ISSU DU MODULE DE RAYONNEMENT P-1
!===============================================================================
if (istpp1.eq.1) then
  ntmabs = ntcabs
endif

!===============================================================================
! 18. TEST D'ARRET PAR DEMANDE DE SYRTHES
!===============================================================================

!     En cas de couplage, on lit des maintenant l'entete du premier
!     message du pas de temps suivant, au cas ou il s'agisse d'un
!     message de terminaison (pas de test sur ITRALE ici, car
!     il serait sur ITRALE + 1, toujours > 0).

call tstsyr (ntmabs, ntcabs)
!==========

!===============================================================================
! 19. SORTIE EVENTUELLE DU FICHIER SUITE
!===============================================================================

iisuit = 0
if(ntcabs.lt.ntmabs) then
  if(ntsuit.eq.0) then
    ntsdef = max((ntmabs-ntpabs)/4,10)
    if(ntsdef.gt.0) then
      if(mod(ntcabs-ntpabs,ntsdef).eq.0) then
        iisuit = 1
      endif
    endif
  elseif(ntsuit.gt.0) then
    if(mod(ntcabs,ntsuit).eq.0) then
      iisuit = 1
    endif
  endif
  if (itrale.eq.0) iisuit = 0
  if (iwarn0.gt.0 .and. iisuit.eq.1) write(nfecra,3020)ntcabs,ttcabs
else if(ntcabs.eq.ntmabs .and. ntsuit.gt.-2) then
  iisuit = 1
  if(iwarn0.gt.0) write(nfecra,3021)ntcabs,ttcabs
endif

if (iisuit.eq.1) then

  call dmtmps(tecrf1)
  !==========

  call ecrava                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nnod   ,          &
   nvar   , nscal  , nphas  ,                                     &
   ia     ,                                                       &
   xyzcen     , surfac     , surfbo     , cdgfac     , cdgfbo    ,&
   ra(idt)    , ra(irtp)   , ra(ipropc) , ra(ipropf) , ra(ipropb),&
   ra(icoefa) , ra(icoefb) , ra(ifrcx)  ,                         &
   ra     )

  if (nfpt1t.gt.0) then
    ficsui = '1dwall_module'
    call ecrt1d                                                   &
    !==========
 ( ficsui   , len(ficsui), nfpt1d   ,  nmxt1d  ,                  &
   nfabor   , ra(itppt1) , ia(iifpt1))
  endif

  if (ippmod(iaeros).ge.0) then
     ficsui = 'cooling_towers'
     call ecrctw ( ficsui , len(ficsui) )
     !==========
  endif

  if (iilagr.gt.0) then

    call lagout                                                   &
    !==========
 ( ifinia , ifinra ,                                              &
   lndnod ,                                                       &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   ia(iicoce) , ia(iityce) , ia(iiitep) ,                         &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(iettp)  , ra(iitepa) , ra(istatf) ,                         &
   ra(istatc) , ra(istatv) , ra(itslag) ,                         &
   ra     )

  endif

  if (iirayo.gt.0) then
    call rayout                                                   &
    !==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra     )
  endif

  call dmtmps(tecrf2)
  !==========

  if(iwarn0.gt.0) write(nfecra,3022) tecrf2-tecrf1

endif ! iisuit = 1

!===============================================================================
! 20. TEST POUR SAVOIR SI ON SORT UN FICHIER POST OU NON
!===============================================================================

call usnpst                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) , ra(istatc) ,                         &
   ra     )

!===============================================================================
! 21. SORTIE DES FICHIERS POST STANDARDS
!===============================================================================

!     Si ITRALE=0 on desactive tous les writers (car la geometrie n'a pas ete
!       ecrite)
if (itrale.eq.0) then
  indwri = 0
  indact = 0
  call pstact(indwri, indact)
  !==========
endif

call pstvar                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ntcabs ,                                                       &
   nvar   , nscal  , nphas  , nvlsta , nvisbr ,                   &
   ia     ,                                                       &
   ttcabs ,                                                       &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(istatc) , ra(istatv) , ra(istatf) ,                         &
   ra     )

!===============================================================================
! 22. HISTORIQUES
!===============================================================================

if ((nthist.gt.0 .or.frhist.gt.0.d0) .and. itrale.gt.0) then

  modhis = -1
  if (frhist.gt.0.d0) then
    if ((ttcabs - ttchis) .gt. frhist*(1-1.0d-6)) then
      modhis = 1
    endif
  else if (nthist.gt.0 .and. mod(ntcabs, nthist).eq.0) then
    modhis = 1
  endif

  if (modhis.eq.1) then

    ttchis = ttcabs

    call ecrhis(ndim, ncelet, ncel, modhis, xyzcen, ra)
    !==========

    if (iilagr.gt.0) then
      call laghis(ndim, ncelet, ncel, modhis, nvlsta,             &
      !==========
                  xyzcen, volume, ra(istatc), ra(istatv))
    endif

    if (ihistr.eq.1) then
      call strhis(modhis)
      !==========
    endif

  endif

endif

call ushist                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra     )


!===============================================================================
! 23. ECRITURE LISTING TOUTES LES NTLIST ITERATIONS
!===============================================================================

if(ntlist.gt.0) then
  modntl = mod(ntcabs,ntlist)
elseif(ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
  modntl = 0
else
  modntl = 1
endif
if(modntl.eq.0) then
  call ecrlis                                                     &
  !==========
  ( ifinia , ifinra ,                                             &
    nvar   , nphas  , ndim   , ncelet , ncel   ,                  &
    irtp   ,                                                      &
    ia     ,                                                      &
    ra(irtp  ) , ra(irtpa ) , ra(idt ) , volume , xyzcen,         &
    ra     )

  if (iilagr.gt.0) then

    call laglis                                                   &
    !==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   ia(iiitep),                                                    &
   ia     ,                                                       &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(iettp)  , ra(iitepa) , ra(istatc) ,  ra(istatv) ,           &
   ra(itslag) , ra(istatf),                                       &
   ra     )

  endif

endif

call dmtmps(titer2)
!==========

if(iwarn0.gt.0) then
  if (itrale.gt.0) then
    write(nfecra,3010)ntcabs,titer2-titer1
  else
    write(nfecra,3012)titer2-titer1
  endif
endif


!===============================================================================
! 24. FIN DE LA BOUCLE EN TEMPS
!===============================================================================

itrale = itrale + 1

if(ntcabs.lt.ntmabs) goto 100


! LIBERATION DES TABLEAUX INTERMEDIAIRES (PDC+TSM)

ifinia = iditia
ifinra = iditra


!===============================================================================
! 25. FINALISATION HISTORIQUES
!===============================================================================

if(iwarn0.gt.0) then
  write(nfecra,4000)
endif

call dmtmps(tecrf1)
!==========

! Ici on sauve les historiques (si on en a stocke)

modhis = 2
call ecrhis(ndim, ncelet, ncel, modhis, xyzcen, ra)
!==========

if (iilagr.gt.0) then
  call laghis(ndim, ncelet, ncel, modhis , nvlsta,                &
  !==========
              xyzcen, volume, ra(istatc), ra(istatv))
endif
if (ihistr.eq.1) then
  call strhis(modhis)
  !==========
endif

!     LE CAS ECHEANT, ON LIBERE LES STRUCTURES C DU MODULE THERMIQUE 1D
!     ET/OU ON FERME LE LISTING LAGRANGIEN

if (nfpt1d.gt.0) then
  call lbrt1d
  !==========
endif

if (iale.gt.0) then
  call lbrale
  !==========
endif

if (iilagr.gt.0) close(implal)

call dmtmps(tecrf2)
!==========

if(iwarn0.gt.0) then
  write(nfecra,4010)tecrf2-tecrf1
endif

 200  continue

!===============================================================================
! 26. MEMOIRE UTILISEE
!===============================================================================

!     Liberation des structures liees a la lecture du fichier xml

if (iihmpr.eq.1) then

  call memui2
  !==========

  call memui1(ncharb)
  !==========
endif


write(nfecra,7000)

!----
! FORMATS
!----

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'===============================================================',&
                                                              /,/,&
'              FONCTIONNEMENT EN MODE VERIFICATION            ',/,&
'              ===================================            ',/,&
'                                                             ',/,&
' =========================================================== ',/,&
                                                                /,&
                                                                /)

 2000 format(/,/,                                                 &
'===============================================================',&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                       CORPS DU CALCUL                       ',/,&
'                       ===============                       ',/,&
                                                                /,&
                                                                /,&
'===============================================================',&
                                                              /,/,&
                                                                /)
 3000 format(/,                                                   &
'===============================================================',&
 /)
 3001 format(/,' INSTANT ',E18.9,        '   ITERATION NUMERO ',I15,/,  &
' ============================================================= ',&
 /,/)
 3002 format(/,' INSTANT ',E18.9,        '   INITIALISATION ALE ',/,    &
' ============================================================= ',&
 /,/)
 3010 format(/,' TEMPS CPU POUR L''ITERATION ',I15,' :    ',E14.5,/,/,  &
'===============================================================',&
 /)
 3012 format(/,' TEMPS CPU POUR L''INITIALISATION ALE :    ',E14.5,/,/, &
'===============================================================',&
 /)
 3020 format(/,/,                                                 &
 ' Sortie intermediaire de fichiers suite',/,                     &
 '   Sauvegarde a l''iteration ', I10, ', Temps physique ',E14.5,/,/)
 3021 format(/,/,                                                 &
 ' Sortie finale de fichiers suite',/,                     &
 '   Sauvegarde a l''iteration ', I10, ', Temps physique ',E14.5,/,/)
 3022 format(/,/,                                                 &
 ' Temps CPU pour les fichiers suite : ',E14.5,/,/)

 4000 format(/,/,                                                 &
'===============================================================',&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                   ETAPES FINALES DU CALCUL                  ',/,&
'                   ========================                  ',/,&
                                                                /,&
                                                                /,&
' =========================================================== ',/,&
                                                                /,&
                                                                /)
 4010 format(                                                   /,&
 3X,'** TEMPS CPU POUR LES SORTIES FINALES : ',E14.5           ,/,&
 3X,'   ----------------------------------                    ',/)
 7000 format(/,/,                                                 &
' =========================================================== ',/,&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                 FIN DE L''EXECUTION DU CALCUL               ',/,&
'                 ============================                ',/,&
                                                                /,&
                                                                /,&
'===============================================================')

#else

 1000 format(/,                                                   &
'===============================================================',&
                                                              /,/,&
'              RUNNING IN VERIFICATION MODE                   ',/,&
'              ============================                   ',/,&
'                                                             ',/,&
' =========================================================== ',/,&
                                                                /,&
                                                                /)

 2000 format(/,/,                                                 &
'===============================================================',&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                       MAIN CALCULATION                      ',/,&
'                       ================                      ',/,&
                                                                /,&
                                                                /,&
'===============================================================',&
                                                              /,/,&
                                                                /)
 3000 format(/,                                                   &
'===============================================================',&
 /)
 3001 format(/,' INSTANT ',E18.9,        '   TIME STEP NUMBER ',I15,/,  &
' ============================================================= ',&
 /,/)
 3002 format(/,' INSTANT ',E18.9,        '   ALE INITIALIZATION ',/,    &
' ============================================================= ',&
 /,/)
 3010 format(/,' CPU TIME FOR THE TIME STEP  ',I15,':     ',E14.5,/,/,  &
'===============================================================',&
 /)
 3012 format(/,' CPU TIME FOR ALE INITIALIZATION:          ',E14.5,/,/, &
'===============================================================',&
 /)
 3020 format(/,/,                                                 &
 ' Write intermediate restart files',/,                           &
 '   checkpoint at iteration ',    I10,  ', Physical time ',E14.5,/,/)
 3021 format(/,/,                                                 &
 ' Write final restart files',/,                           &
 '   checkpoint at iteration ',    I10,  ', Physical time ',E14.5,/,/)
 3022 format(/,/,                                                 &
 ' CPU time for restart files: ',E14.5,/,/)

 4000 format(/,/,                                                 &
'===============================================================',&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                 FINAL STAGE OF THE CALCULATION              ',/,&
'                 ==============================              ',/,&
                                                                /,&
                                                                /,&
' =========================================================== ',/,&
                                                                /,&
                                                                /)
 4010 format(                                                   /,&
 3X,'** CPU TIME FOR FINAL WRITING: ',E14.5                    ,/,&
 3X,'   ---------------------------                           ',/)
 7000 format(/,/,                                                 &
' =========================================================== ',/,&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                      END OF CALCULATION                     ',/,&
'                      ==================                     ',/,&
                                                                /,&
                                                                /,&
'===============================================================')

#endif

!===============================================================================
! 26. End
!===============================================================================

return
end subroutine
