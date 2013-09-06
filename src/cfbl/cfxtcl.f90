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

subroutine cfxtcl &
!================

 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   rcodcl )

!===============================================================================
! FONCTION :
! --------

!    CONDITIONS AUX LIMITES AUTOMATIQUES

!           COMPRESSIBLE SANS CHOC


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
! itrifb           ! ia ! <-- ! indirection for boundary faces ordering        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

! Arguments

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvarcl)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvarcl,3)

! Local variables

integer          ivar  , ifac  , iel
integer          ii    , iii   , imodif, iccfth
integer          icalep, icalgm
integer          iflmab, ipcrom, ipbrom
integer          ien   , itk
integer          iclp
integer          iclu  , iclv  , iclw
integer          nvarcf

integer          nvcfmx
parameter       (nvcfmx=6)
integer          ivarcf(nvcfmx)

double precision hint  , gammag

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w7
double precision, allocatable, dimension(:) :: wbfb
double precision, allocatable, dimension(:,:) :: bval

double precision, dimension(:), pointer :: bmasfl
double precision, dimension(:), pointer :: coefbp

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet), w6(ncelet))

allocate(w7(nfabor), wbfb(nfabor))
allocate(bval(nfabor,nvar))

ien = isca(ienerg)
itk = isca(itempk)
iclp   = iclrtp(ipr,icoef)
iclu   = iclrtp(iu ,icoef)
iclv   = iclrtp(iv ,icoef)
iclw   = iclrtp(iw ,icoef)

call field_get_key_int(ivarfl(ien), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

ipbrom = ipprob(irom)
ipcrom = ipproc(irom)

!     Liste des variables compressible :
ivarcf(1) = ipr
ivarcf(2) = iu
ivarcf(3) = iv
ivarcf(4) = iw
ivarcf(5) = ien
ivarcf(6) = itk
nvarcf    = 6

call field_get_coefb_s(ivarfl(ipr), coefbp)
do ifac = 1, nfabor
  wbfb(ifac) = coefbp(ifac)
enddo

!     Calcul de epsilon_sup = e - CvT
!       On en a besoin si on a des parois a temperature imposee.
!       Il est calculé aux cellules W5 et aux faces de bord COEFU.
!       On n'en a besoin ici qu'aux cellules de bord : s'il est
!         nécessaire de gagner de la mémoire, on pourra modifier
!         cfther.

icalep = 0
do ifac = 1, nfabor
  if(icodcl(ifac,itk).eq.5) then
    icalep = 1
  endif
enddo
if(icalep.ne.0) then
  iccfth = 7
  imodif = 0
  call cfther                                                    &
  !==========
( nvar   ,                                                       &
  iccfth , imodif ,                                              &
  dt     , rtp    , rtpa   , propce , propfb ,                   &
  w5     , w7     , w3     , w4     , rvoid  , rvoid  )
endif


!     Calcul de gamma (constant ou variable ; pour le moment : cst)
!       On en a besoin pour les entrees sorties avec rusanov

icalgm = 0
do ifac = 1, nfabor
  if ( ( itypfb(ifac).eq.iesicf ) .or.                    &
       ( itypfb(ifac).eq.isopcf ) .or.                    &
       ( itypfb(ifac).eq.ierucf ) .or.                    &
       ( itypfb(ifac).eq.ieqhcf ) ) then
    icalgm = 1
  endif
enddo
if(icalgm.ne.0) then
  iccfth = 1
  imodif = 0
  call cfther                                                    &
  !==========
( nvar   ,                                                       &
  iccfth , imodif ,                                              &
  dt     , rtp    , rtpa   , propce , propfb ,                   &
  w1     , w2     , w6     , w4     , rvoid  , rvoid  )

  if(ieos.eq.1) then
    gammag = w6(1)
  else
!     Gamma doit etre passe a cfrusb ; s'il est variable
!       il est dans le tableau W6 et il faut ajouter
!           GAMMAG = W6(IFABOR(IFAC)) selon IEOS
!       dans la boucle sur les faces.
!     En attendant que IEOS different de 1 soit code, on stoppe
    write(nfecra,7000)
    call csexit (1)
  endif

endif



!     Boucle sur les faces

do ifac = 1, nfabor
  iel = ifabor(ifac)

!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES DE PAROI
!===============================================================================

  if ( itypfb(ifac).eq.iparoi) then

!     Les RCODCL ont ete initialises a -RINFIN pour permettre de
!       verifier ceux que l'utilisateur a modifies. On les remet a zero
!       si l'utilisateur ne les a pas modifies.
!       En paroi, on traite toutes les variables.
    do ivar = 1, nvar
      if(rcodcl(ifac,ivar,1).le.-rinfin*0.5d0) then
        rcodcl(ifac,ivar,1) = 0.d0
      endif
    enddo

!     Le flux de masse est nul

    bmasfl(ifac) = 0.d0

!     Pression :

!       Si la gravite est predominante : pression hydrostatique
!         (approximatif et surtout explicite en rho)

    if(icfgrp.eq.1) then

      icodcl(ifac,ipr) = 3
      hint = dt(iel)/distb(ifac)
      rcodcl(ifac,ipr,3) = -hint                             &
           * ( gx*(cdgfbo(1,ifac)-xyzcen(1,iel))                    &
           + gy*(cdgfbo(2,ifac)-xyzcen(2,iel))                    &
           + gz*(cdgfbo(3,ifac)-xyzcen(3,iel)) )                  &
           * propce(iel,ipcrom)

    else

!       En général : proportionnelle a la valeur interne
!         (Pbord = COEFB*Pi)
!       Si on détend trop : Dirichlet homogene

      iccfth = 91

      call cfther                                                 &
      !==========
 ( nvar   ,                                                       &
   iccfth , ifac   ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   w1     , w2     , w3     , w4     , bval   , wbfb )

!       En outre, il faut appliquer une pre-correction pour compenser
!        le traitement fait dans condli... Si on pouvait remplir COEFA
!        et COEFB directement, on gagnerait en simplicite, mais cela
!        demanderait un test sur IPPMOD dans condli : à voir)

!FIXME with the new cofaf
      icodcl(ifac,ipr) = 1
      if(wbfb(ifac).lt.rinfin*0.5d0.and.                  &
         wbfb(ifac).gt.0.d0  ) then
        hint = dt(iel)/distb(ifac)
        rcodcl(ifac,ipr,1) = 0.d0
        rcodcl(ifac,ipr,2) = hint*(1.d0/wbfb(ifac)-1.d0)
      else
        rcodcl(ifac,ipr,1) = 0.d0
      endif

    endif


!       La vitesse et la turbulence sont traitées de manière standard,
!         dans condli.

!       Pour la thermique, on doit effectuer ici un prétraitement,
!         la variable résolue étant l'energie
!         (energie interne+epsilon sup+energie cinétique). En particulier
!         lorsque la paroi est à température imposée, on prépare le
!         travail de clptur. Hormis l'énergie résolue, toutes les
!         variables rho et s prendront arbitrairement une condition de
!         flux nul (leurs conditions aux limites ne servent qu'à la
!         reconstruction des gradients et il parait délicat d'imposer
!         autre chose qu'un flux nul sans risque de créer des valeurs
!         aberrantes au voisinage de la couche limite)

!       Par défaut : adiabatique
    if(  icodcl(ifac,itk).eq.0.and.                          &
         icodcl(ifac,ien).eq.0) then
      icodcl(ifac,itk) = 3
      rcodcl(ifac,itk,3) = 0.d0
    endif

!       Temperature imposee
    if(icodcl(ifac,itk).eq.5) then

!           On impose la valeur de l'energie qui conduit au bon flux.
!             On notera cependant qu'il s'agit de la condition à la
!               limite pour le flux diffusif. Pour la reconstruction
!               des gradients, il faudra utiliser autre chose.
!               Par exemple un flux nul ou encore toute autre
!               condition respectant un profil : on se calquera sur
!               ce qui sera fait pour la température si c'est possible,
!               sachant que l'energie contient l'energie cinetique,
!               ce qui rend le choix du profil délicat.

      icodcl(ifac,ien) = 5
      if(icv.eq.0) then
        rcodcl(ifac,ien,1) = cv0*rcodcl(ifac,itk,1)
      else
        rcodcl(ifac,ien,1) = propce(iel,ipproc(icv))  &
             *rcodcl(ifac,itk,1)
      endif
      rcodcl(ifac,ien,1) = rcodcl(ifac,ien,1)             &
           + 0.5d0*(rtp(iel,iu)**2+rtp(iel,iv)**2+rtp(iel,iw)**2)          &
           + w5(iel)
!                   ^epsilon sup (cf USCFTH)

!           Les flux en grad epsilon sup et énergie cinétique doivent
!             être nuls puisque tout est pris par le terme de
!             diffusion d'energie.
      ifbet(ifac) = 1

!           Flux nul pour la reconstruction éventuelle de température
      icodcl(ifac,itk) = 3
      rcodcl(ifac,itk,3) = 0.d0

!       Flux impose
    elseif(icodcl(ifac,itk).eq.3) then

!           On impose le flux sur l'energie
      icodcl(ifac,ien) = 3
      rcodcl(ifac,ien,3) = rcodcl(ifac,itk,3)

!           Les flux en grad epsilon sup et énergie cinétique doivent
!             être nuls puisque tout est pris par le terme de
!             diffusion d'energie.
      ifbet(ifac) = 1

!           Flux nul pour la reconstruction éventuelle de température
      icodcl(ifac,itk) = 3
      rcodcl(ifac,itk,3) = 0.d0

    endif


!     Scalaires : flux nul (par defaut dans typecl pour iparoi)


!===============================================================================
! 3.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES DE SYMETRIE
!===============================================================================

  elseif ( itypfb(ifac).eq.isymet ) then

!     Les RCODCL ont ete initialises a -RINFIN pour permettre de
!       verifier ceux que l'utilisateur a modifies. On les remet a zero
!       si l'utilisateur ne les a pas modifies.
!       En symetrie, on traite toutes les variables.
    do ivar = 1, nvar
      if(rcodcl(ifac,ivar,1).le.-rinfin*0.5d0) then
        rcodcl(ifac,ivar,1) = 0.d0
      endif
    enddo

!     Le flux de masse est nul

    bmasfl(ifac) = 0.d0

!     Condition de Pression

    iccfth = 90

    call cfther                                                   &
    !==========
 ( nvar   ,                                                       &
   iccfth , ifac   ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   w1     , w2     , w3     , w4     , bval   ,rvoid )


!     Pression :
!       En général : proportionnelle a la valeur interne
!         (Pbord = COEFB*Pi)
!       Si on détend trop : Dirichlet homogene

!       En outre, il faut appliquer une pre-correction pour compenser le
!        traitement fait dans condli... Si on pouvait remplir COEFA
!        et COEFB directement, on gagnerait en simplicite, mais cela
!        demanderait un test sur IPPMOD dans condli : à voir)

    icodcl(ifac,ipr) = 3
    rcodcl(ifac,ipr,1) = 0.d0
    rcodcl(ifac,ipr,2) = rinfin
    rcodcl(ifac,ipr,3) = 0.d0

!       Toutes les autres variables prennent un flux nul (sauf la vitesse
!         normale, qui est nulle) : par defaut dans typecl pour isymet.

!===============================================================================
! 4.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE/SORTIE
!       ETAPE DE THERMO
!===============================================================================


!===============================================================================
!     4.1 Entree/sortie imposée (par exemple : entree supersonique)
!===============================================================================

  elseif ( itypfb(ifac).eq.iesicf ) then

!     On a
!       - la vitesse,
!       - 2 variables parmi P, rho, T, E (mais pas (T,E)),
!       - la turbulence
!       - les scalaires

!     On recherche la variable a initialiser
!       (si on a donne une valeur nulle, c'est pas adapte : on supposera
!        qu'on n'a pas initialise et on sort en erreur)
    iccfth = 10000
    if(rcodcl(ifac,ipr,1).gt.0.d0) iccfth = 2*iccfth
    if(rcodcl(ifac,itk,1).gt.0.d0) iccfth = 5*iccfth
    if(rcodcl(ifac,ien,1).gt.0.d0) iccfth = 7*iccfth
    if((iccfth.le.70000.and.iccfth.ne.60000).or.                &
         (iccfth.eq.350000)) then
      write(nfecra,1000)iccfth
      call csexit (1)
    endif
    iccfth = iccfth + 900

!     Les RCODCL ont ete initialises a -RINFIN pour permettre de
!       verifier ceux que l'utilisateur a modifies. On les remet a zero
!       si l'utilisateur ne les a pas modifies.
!       On traite d'abord les variables autres que la turbulence et les
!       scalaires passifs : celles-ci sont traitees plus bas.
    do iii = 1, nvarcf
      ivar = ivarcf(iii)
      if(rcodcl(ifac,ivar,1).le.-rinfin*0.5d0) then
        rcodcl(ifac,ivar,1) = 0.d0
      endif
    enddo

!     On calcule les variables manquantes parmi P,rho,T,E
!     COEFA sert de tableau de transfert dans USCFTH

    do ivar = 1, nvar
      bval(ifac,ivar) = rcodcl(ifac,ivar,1)
    enddo

    call cfther                                                   &
    !==========
 ( nvar   ,                                                       &
   iccfth , ifac   ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   w1     , w2     , w3     , w4     , bval   , rvoid )


!     Rusanov, flux de masse et type de conditions aux limites :
!       voir plus bas


!===============================================================================
!     4.2 Sortie supersonique
!===============================================================================

  elseif ( itypfb(ifac).eq.isspcf ) then

!     On impose un Dirichlet égal à la valeur interne pour rho u E
!       (on impose des Dirichlet déduit pour les autres variables).
!       Il est inutile de passer dans Rusanov.
!     Il serait nécessaire de reconstruire ces valeurs en utilisant
!       leur gradient dans la cellule de bord : dans un premier temps,
!       on utilise des valeurs non reconstruites (non consistant mais
!       potentiellement plus stable).
!     On pourrait imposer des flux nuls (a tester), ce qui éviterait
!       la nécessité de reconstruire les valeurs.

!     Les RCODCL ont ete initialises a -RINFIN pour permettre de
!       verifier ceux que l'utilisateur a modifies. On les remet a zero
!       si l'utilisateur ne les a pas modifies.
!       On traite d'abord les variables autres que la turbulence et les
!       scalaires passifs : celles-ci sont traitees plus bas.
    do iii = 1, nvarcf
      ivar = ivarcf(iii)
      if(rcodcl(ifac,ivar,1).le.-rinfin*0.5d0) then
        rcodcl(ifac,ivar,1) = 0.d0
      endif
    enddo

!     Valeurs de rho u E
    propfb(ifac,ipbrom) = propce(iel,irom)
    rcodcl(ifac,iu ,1) = rtp(iel,iu)
    rcodcl(ifac,iv ,1) = rtp(iel,iv)
    rcodcl(ifac,iw ,1) = rtp(iel,iw)
    rcodcl(ifac,ien,1) = rtp(iel,ien)

!     Valeurs de P et s déduites
    iccfth = 924

    do ivar = 1, nvar
      bval(ifac,ivar) = rcodcl(ifac,ivar,1)
    enddo

    call cfther                                                   &
    !==========
 ( nvar   ,                                                       &
   iccfth , ifac   ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   w1     , w2     , w3     , w4     , bval   , rvoid )

!               flux de masse et type de conditions aux limites :
!       voir plus bas


!===============================================================================
!     4.3 Sortie a pression imposee
!===============================================================================

  elseif ( itypfb(ifac).eq.isopcf ) then

!       Sortie subsonique a priori (si c'est supersonique dans le
!         domaine, ce n'est pas pour autant que c'est supersonique
!         à la sortie, selon la pression que l'on a imposée)

!     On utilise un scenario dans lequel on a une 1-détente et un
!       2-contact entrant dans le domaine. On détermine les conditions
!       sur l'interface selon la thermo et on passe dans Rusanov
!       ensuite pour lisser.

!     Si P n'est pas donné, erreur ; on sort aussi en erreur si P
!       négatif, même si c'est possible, dans la plupart des cas ce
!       sera une erreur
    if(rcodcl(ifac,ipr,1).lt.-rinfin*0.5d0) then
      write(nfecra,1100)
      call csexit (1)
    endif

!     Les RCODCL ont ete initialises a -RINFIN pour permettre de
!       verifier ceux que l'utilisateur a modifies. On les remet a zero
!       si l'utilisateur ne les a pas modifies.
!       On traite d'abord les variables autres que la turbulence et les
!       scalaires passifs : celles-ci sont traitees plus bas.
    do iii = 1, nvarcf
      ivar = ivarcf(iii)
      if(rcodcl(ifac,ivar,1).le.-rinfin*0.5d0) then
        rcodcl(ifac,ivar,1) = 0.d0
      endif
    enddo

!     Valeurs de rho, u, E, s
    iccfth = 93

    do ivar = 1, nvar
      bval(ifac,ivar) = rcodcl(ifac,ivar,1)
    enddo

    call cfther                                                   &
    !==========
 ( nvar   ,                                                       &
   iccfth , ifac   ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   w1     , w2     , w3     , w4     , bval   , rvoid )

!     Rusanov, flux de masse et type de conditions aux limites :
!       voir plus bas


!===============================================================================
!     4.4 Entree à rho et U imposes
!===============================================================================

  elseif ( itypfb(ifac).eq.ierucf ) then

!       Entree subsonique a priori (si c'est supersonique dans le
!         domaine, ce n'est pas pour autant que c'est supersonique
!         à l'entree, selon les valeurs que l'on a imposées)

!     On utilise un scenario détente ou choc.
!       On détermine les conditions sur l'interface
!       selon la thermo et on passe dans Rusanov ensuite pour lisser.

!     Si rho et u ne sont pas donnés, erreur

!   FIXME: Implement Q,H boundary condition
    if(rcodcl(ifac,iu ,1).lt.-rinfin*0.5d0.or.               &
         rcodcl(ifac,iv ,1).lt.-rinfin*0.5d0.or.               &
         rcodcl(ifac,iw ,1).lt.-rinfin*0.5d0) then
      write(nfecra,1200)
      call csexit (1)
    endif

!     Les RCODCL ont ete initialises a -RINFIN pour permettre de
!       verifier ceux que l'utilisateur a modifies. On les remet a zero
!       si l'utilisateur ne les a pas modifies.
!       On traite d'abord les variables autres que la turbulence et les
!       scalaires passifs : celles-ci sont traitees plus bas.
    do iii = 1, nvarcf
      ivar = ivarcf(iii)
      if(rcodcl(ifac,ivar,1).le.-rinfin*0.5d0) then
        rcodcl(ifac,ivar,1) = 0.d0
      endif
    enddo

!     Valeurs de P, E, s
    iccfth = 92

    do ivar = 1, nvar
      bval(ifac,ivar) = rcodcl(ifac,ivar,1)
    enddo

    call cfther                                                   &
    !==========
 ( nvar   ,                                                       &
   iccfth , ifac   ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   w1     , w2     , w3     , w4     , bval   , rvoid )

!     Rusanov, flux de masse et type de conditions aux limites :
!       voir plus bas

!===============================================================================
!     4.5 Entree à P et H imposees
!===============================================================================

  elseif ( itypfb(ifac).eq.iephcf ) then

!       Entree subsonique a priori (si c'est supersonique dans le
!         domaine, ce n'est pas pour autant que c'est supersonique
!         à l'entree, selon les valeurs que l'on a imposées)

!     On utilise un scenario détente ou choc.
!       On détermine les conditions sur l'interface
!       selon la thermo.

!     Si P et H ne sont pas donnés, erreur

! rcodcl(ifac,isca(ienerg),1) holds the boundary total enthalpy values prescribed by the user
    if(rcodcl(ifac,ipr ,1).lt.-rinfin*0.5d0.or.               &
         rcodcl(ifac,isca(ienerg) ,1).lt.-rinfin*0.5d0) then
      write(nfecra,1200)
      call csexit (1)
    endif

!     Les RCODCL ont ete initialises a -RINFIN pour permettre de
!       verifier ceux que l'utilisateur a modifies. On les remet a zero
!       si l'utilisateur ne les a pas modifies.
!       On traite d'abord les variables autres que la turbulence et les
!       scalaires passifs : celles-ci sont traitees plus bas.
    do iii = 1, nvarcf
      ivar = ivarcf(iii)
      if(rcodcl(ifac,ivar,1).le.-rinfin*0.5d0) then
        rcodcl(ifac,ivar,1) = 0.d0
      endif
    enddo

!     Valeurs de rho, U, E
    iccfth = 95

    do ivar = 1, nvar
      bval(ifac,ivar) = rcodcl(ifac,ivar,1)
    enddo

    call cfther                                                   &
    !==========
 ( nvar   ,                                                       &
   iccfth , ifac   ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   w1     , w2     , w3     , w4     , bval   , rvoid )

!     flux de masse et type de conditions aux limites :
!     voir plus bas


!===============================================================================
!     4.6 Entree à rho*U et rho*U*H imposes
!===============================================================================

  elseif ( itypfb(ifac).eq.ieqhcf ) then

!       Entree subsonique a priori (si c'est supersonique dans le
!         domaine, ce n'est pas pour autant que c'est supersonique
!         à l'entree, selon les valeurs que l'on a imposées)

!     On utilise un scenario dans lequel on a un 2-contact et une
!       3-détente entrant dans le domaine. On détermine les conditions
!       sur l'interface selon la thermo et on passe dans Rusanov
!       ensuite pour lisser.

!     Si rho et u ne sont pas donnés, erreur
! TODO implement the Q,H boundary condition
    if(rcodcl(ifac,irunh,1).lt.-rinfin*0.5d0) then
      write(nfecra,1300)
      call csexit (1)
    endif

!     Les RCODCL ont ete initialises a -RINFIN pour permettre de
!       verifier ceux que l'utilisateur a modifies. On les remet a zero
!       si l'utilisateur ne les a pas modifies.
!       On traite d'abord les variables autres que la turbulence et les
!       scalaires passifs : celles-ci sont traitees plus bas.
    do iii = 1, nvarcf
      ivar = ivarcf(iii)
      if(rcodcl(ifac,ivar,1).le.-rinfin*0.5d0) then
          rcodcl(ifac,ivar,1) = 0.d0
      endif
    enddo

!     A coder

!     IRUNH = ISCA(IENER)
!     (aliases pour simplifier uscfcl)

    write(nfecra,1301)
    call csexit (1)

!===============================================================================
! 5. CONDITION NON PREVUE
!===============================================================================
!     Stop
  else

    write(nfecra,1400)
    call csexit (1)

! --- Fin de test sur les types de faces
  endif


!===============================================================================
! 6. FIN DU TRAITEMENT DES ENTREE/SORTIES
!     CALCUL DU FLUX DE MASSE,
!     CALCUL DES FLUX DE BORD AVEC RUSANOV (SI BESOIN)
!     TYPE DE C    .L. (DIRICHLET NEUMANN)
!===============================================================================

  if ( ( itypfb(ifac).eq.iesicf ) .or.                    &
       ( itypfb(ifac).eq.isspcf ) .or.                    &
       ( itypfb(ifac).eq.iephcf ) .or.                    &
       ( itypfb(ifac).eq.isopcf ) .or.                    &
       ( itypfb(ifac).eq.ierucf ) .or.                    &
       ( itypfb(ifac).eq.ieqhcf ) ) then

!===============================================================================
!     6.1 Flux de bord Rusanov ou simplement flux de masse
!         Attention a bien avoir calcule gamma pour Rusanov
!===============================================================================

!     Sortie supersonique :
    if ( itypfb(ifac).eq.isspcf ) then

!     Seul le flux de masse est calcule (on n'appelle pas Rusanov)
!       (toutes les variables sont connues)

      bmasfl(ifac) = propfb(ifac,ipbrom) *                                     &
                     ( bval(ifac,iu)*surfbo(1,ifac)                            &
                     + bval(ifac,iv)*surfbo(2,ifac)                            &
                     + bval(ifac,iw)*surfbo(3,ifac) )

!     Entree subsonique

    else if ( itypfb(ifac).eq.ierucf ) then

!     Seul le flux de masse est calcule (on n'appelle pas Rusanov)

      bmasfl(ifac) = propfb(ifac,ipbrom) *                                     &
                     ( bval(ifac,iu)*surfbo(1,ifac)                            &
                     + bval(ifac,iv)*surfbo(2,ifac)                            &
                     + bval(ifac,iw)*surfbo(3,ifac) )



!     Autres entrees/sorties :
    else

!     On calcule des flux par Rusanov (PROPFB)
!       (en particulier, le flux de masse est complete)
!       pour la condition d'entree supersonique seulement

      if ( itypfb(ifac).eq.iesicf ) then

        call cfrusb                                                   &
        !==========
     ( nvar   , nscal  ,                                              &
       ifac   ,                                                       &
       gammag ,                                                       &
       dt     , rtp    , rtpa   , propce , propfb , bval ,            &
       w3     , w4     )

!    Pour les autres types (sortie subsonique, entree QH, entree PH),
!    On calcule des flux analytiques

      else

        call cffana                                                   &
        !==========
      ( nvar   ,  ifac   , propfb , bval )

      endif

    endif

!===============================================================================
!     6.2 Recuperation de COEFA
!===============================================================================

!     On rétablit COEFA dans RCODCL
    do ivar = 1, nvar
      rcodcl(ifac,ivar,1) = bval(ifac,ivar)
    enddo

!===============================================================================
!     6.3 Types de C.L.
!===============================================================================

!     P               : Dirichlet sauf IESICF : Neumann (choix arbitraire)
!     rho, U, E, T    : Dirichlet
!     k, R, eps, scal : Dirichlet/Neumann selon flux de masse

!     Pour P, le Neumann est censé etre moins genant pour les
!       reconstructions de gradient si la valeur de P fournie par
!       l'utilisateur est tres differente de la valeur interne.
!       Le choix est cependant arbitraire.

!     On suppose que par defaut,
!            RCODCL(IFAC,X,1) = utilisateur ou calcule ci-dessus
!            RCODCL(IFAC,X,2) = RINFIN
!            RCODCL(IFAC,X,3) = 0.D0
!       et si ICODCL(IFAC,X) = 3, seul RCODCL(IFAC,X,3) est utilisé


!-------------------------------------------------------------------------------
!     Pression : Dirichlet ou Neumann homogene
!-------------------------------------------------------------------------------

      icodcl(ifac,ipr)   = 13

!-------------------------------------------------------------------------------
!     rho U E T : Dirichlet
!-------------------------------------------------------------------------------

!     Vitesse
    icodcl(ifac,iu)    = 1
    icodcl(ifac,iv)    = 1
    icodcl(ifac,iw)    = 1
!     Energie totale
    icodcl(ifac,ien)   = 1
!     Temperature
    icodcl(ifac,itk)   = 1

!-------------------------------------------------------------------------------
!     turbulence et scalaires passifs : Dirichlet/Neumann selon flux
!-------------------------------------------------------------------------------

!       Dirichlet ou Neumann homogène
!       On choisit un Dirichlet si le flux de masse est entrant et
!       que l'utilisateur a donné une valeur dans RCODCL

    if (bmasfl(ifac).ge.0.d0) then
      if(itytur.eq.2) then
        icodcl(ifac,ik ) = 3
        icodcl(ifac,iep) = 3
      elseif(itytur.eq.3) then
        icodcl(ifac,ir11) = 3
        icodcl(ifac,ir22) = 3
        icodcl(ifac,ir33) = 3
        icodcl(ifac,ir12) = 3
        icodcl(ifac,ir13) = 3
        icodcl(ifac,ir23) = 3
        icodcl(ifac,iep ) = 3
      elseif(iturb.eq.50) then
        icodcl(ifac,ik  ) = 3
        icodcl(ifac,iep ) = 3
        icodcl(ifac,iphi) = 3
        icodcl(ifac,ifb ) = 3
      elseif(iturb.eq.60) then
        icodcl(ifac,ik  ) = 3
        icodcl(ifac,iomg) = 3
      elseif(iturb.eq.70) then
        icodcl(ifac,inusa) = 3
      endif
      if(nscaus.gt.0) then
        do ii = 1, nscaus
          icodcl(ifac,isca(ii)) = 3
        enddo
      endif
    else
      if(itytur.eq.2) then
        if(rcodcl(ifac,ik ,1).gt.0.d0.and.               &
             rcodcl(ifac,iep,1).gt.0.d0) then
          icodcl(ifac,ik ) = 1
          icodcl(ifac,iep) = 1
        else
          icodcl(ifac,ik ) = 3
          icodcl(ifac,iep) = 3
        endif
      elseif(itytur.eq.3) then
        if(rcodcl(ifac,ir11,1).gt.0.d0.and.              &
             rcodcl(ifac,ir22,1).gt.0.d0.and.              &
             rcodcl(ifac,ir33,1).gt.0.d0.and.              &
             rcodcl(ifac,ir12,1).gt.-rinfin*0.5d0.and.     &
             rcodcl(ifac,ir13,1).gt.-rinfin*0.5d0.and.     &
             rcodcl(ifac,ir23,1).gt.-rinfin*0.5d0.and.     &
             rcodcl(ifac,iep ,1).gt.0.d0) then
          icodcl(ifac,ir11) = 1
          icodcl(ifac,ir22) = 1
          icodcl(ifac,ir33) = 1
          icodcl(ifac,ir12) = 1
          icodcl(ifac,ir13) = 1
          icodcl(ifac,ir23) = 1
          icodcl(ifac,iep ) = 1
        else
          icodcl(ifac,ir11) = 3
          icodcl(ifac,ir22) = 3
          icodcl(ifac,ir33) = 3
          icodcl(ifac,ir12) = 3
          icodcl(ifac,ir13) = 3
          icodcl(ifac,ir23) = 3
          icodcl(ifac,iep ) = 3
        endif
      elseif(iturb.eq.50) then
        if(rcodcl(ifac,ik  ,1).gt.0.d0.and.              &
             rcodcl(ifac,iep ,1).gt.0.d0.and.              &
             rcodcl(ifac,iphi,1).gt.0.d0.and.              &
             rcodcl(ifac,ifb ,1).gt.-rinfin*0.5d0 ) then
          icodcl(ifac,ik  ) = 1
          icodcl(ifac,iep ) = 1
          icodcl(ifac,iphi) = 1
          icodcl(ifac,ifb ) = 1
        else
          icodcl(ifac,ik  ) = 3
          icodcl(ifac,iep ) = 3
          icodcl(ifac,iphi) = 3
          icodcl(ifac,ifb ) = 3
        endif
      elseif(iturb.eq.60) then
         if(rcodcl(ifac,ik  ,1).gt.0.d0.and.               &
              rcodcl(ifac,iomg,1).gt.0.d0 ) then
           icodcl(ifac,ik  ) = 1
           icodcl(ifac,iomg) = 1
         else
           icodcl(ifac,ik  ) = 3
           icodcl(ifac,iomg) = 3
         endif
       elseif(iturb.eq.70) then
         if(rcodcl(ifac,inusa,1).gt.0.d0) then
           icodcl(ifac,inusa) = 1
         else
           icodcl(ifac,inusa) = 3
         endif
       endif
       if(nscaus.gt.0) then
         do ii = 1, nscaus
           if(rcodcl(ifac,isca(ii),1).gt.-rinfin*0.5d0) then
             icodcl(ifac,isca(ii)) = 1
           else
             icodcl(ifac,isca(ii)) = 3
           endif
         enddo
       endif
     endif


!     Les RCODCL ont ete initialises a -RINFIN pour permettre de
!       verifier ceux que l'utilisateur a modifies. On les remet a zero
!       si l'utilisateur ne les a pas modifies.
!       On traite la turbulence et les scalaires passifs (pour
!       simplifier la boucle, on traite toutes les variables : les
!       variables du compressible sont donc vues deux fois, mais ce
!       n'est pas grave).
     do ivar = 1, nvar
       if(rcodcl(ifac,ivar,1).le.-rinfin*0.5d0) then
         rcodcl(ifac,ivar,1) = 0.d0
       endif
     enddo


! --- Fin de test sur les faces d'entree sortie
   endif

! --- Fin de boucle sur les faces de bord
 enddo

! Free memory
deallocate(w1, w2, w3)
deallocate(w4, w5, w6)
deallocate(w7)
deallocate(bval)

!----
! FORMATS
!----

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,&
'@    =========                                               ',/,&
'@    Deux variables independantes et deux seulement parmi    ',/,&
'@    P, rho, T et E doivent etre imposees aux bords de type  ',/,&
'@    IESICF dans uscfcl (ICCFTH = ',I10,').                  ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites dans uscfcl.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,&
'@    =========                                               ',/,&
'@    La pression n''a pas ete fournie en sortie a pression   ',/,&
'@    imposée.                                                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites dans uscfcl.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,&
'@    =========                                               ',/,&
'@    La masse volumique ou la vitesse n''a pas été fournie   ',/,&
'@    en entree a masse volumique et vitesse imposee.         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites dans uscfcl.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1300 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,&
'@    =========                                               ',/,&
'@    Le debit massique ou le debit enthalpique n''a pas été  ',/,&
'@    fourni en entree a debit massique et enthalpique imposé.',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites dans uscfcl.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1301 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,&
'@    =========                                               ',/,&
'@    Entree à debit massique et debit enthalpique non prevue ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Contacter l''equipe de developpement pour uscfcl.         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1400 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,&
'@    =========                                               ',/,&
'@    Une condition a la limite ne fait pas partie des        ',/,&
'@      conditions aux limites predefinies en compressible.   ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites dans uscfcl.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,&
'@    =========                                               ',/,&
'@    cfxtcl doit etre modifie pour prendre en compte une loi ',/,&
'@      d''etat a gamma variable. Seul est pris en compte le  ',/,&
'@      cas IEOS = 1                                          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IEOS dans cfther.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
!----
! FIN
!----

return
end subroutine
