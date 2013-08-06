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

subroutine schtmp &
!================

 ( nvar   , nscal  , iappel ,                                     &
   isostd ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  )

!===============================================================================
! FONCTION :
! ----------

! TRAITEMENT DU FLUX DE MASSE, DE LA VISCOSITE, DE LA MASSE
! VOLUMIQUE, DE LA CHALEUR SPECIFIQUE ET DU TABLEAU TSNSA
! DANS LE CAS D'UN THETA SCHEMA

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iappel           ! e  ! <-- ! numero de l'appel (avant ou apres              !
!                  !    !     ! phyvar                                         !
! isostd           ! te ! <-- ! indicateur de sortie standard                  !
!    (nfabor+1)    !    !     !  +numero de la face de reference               !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
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
use pointe
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal  , iappel

integer          isostd(nfabor+1)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

! Local variables

integer          iel    , ifac   , iscal
integer          ipcrom , ipcroa
integer          ipbrom , ipbroa
integer          iflmas , iflmab , iflmba, iflmsa
integer          ipcvis , ipcvst
integer          ipcvsa , ipcvta , ipcvsl
integer          iicp   , iicpa
double precision flux   , theta  , aa, bb, viscos, xmasvo, varcp
double precision, dimension(:), pointer :: imasfl, bmasfl

!===============================================================================


!===============================================================================
! 0. INITIALISATION
!===============================================================================


!===============================================================================
! 1.  AU TOUT DEBUT DE LA BOUCLE EN TEMPS
!===============================================================================

if (iappel.eq.1) then

! --- Application du schema en temps sur le flux de masse
!     Soit F le flux de masse
!       - si istmpf = 2 (schema non standard, theta>0 en fait = 0.5)
!         imasfl           contient F_(n-2+theta) et on y met F(n-1+theta)
!         propfa(1,iflmsa) contient F_(n-1+theta) et
!                                    on y met une extrapolation en n+theta
!       - si istmpf = 1 (schema non standard, theta=0)
!         imasfl           contient deja F_n et
!         propfa(1,iflmsa) n'est pas utilise : on ne fait rien
!       - sinon : istmpf = 0 (schema standard, theta= -999)
!         imasfl           et propfa(1,iflmsa) contiennent tous les deux F(n)
!                                            : on ne fait rien


!     Ordre 2 en temps pour le flux (theta = 0.5) a entrer dans navsto
!       Flux convectif = 2F(n-1+theta)-F(n-2+theta)
!       Flux conserve  =  F(n-1+theta)
!     Au premier pas de temps, l'ancien a ete initialise dans inivar (a 0)
!       en suite de calcul, les deux ont ete relus.

  if (istmpf.eq.2) then
    call field_get_key_int(ivarfl(iu), kimasf, iflmas)
    call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
    call field_get_val_s(iflmas, imasfl)
    call field_get_val_s(iflmab, bmasfl)
    iflmsa = ipprof(ifluaa(iu))
    iflmba = ipprob(ifluaa(iu))
    do ifac = 1 , nfac
      flux                =      imasfl(ifac)
      imasfl(ifac) = 2.d0*imasfl(ifac)            &
           - propfa(ifac,iflmsa)
      propfa(ifac,iflmsa) = flux
    enddo
    do ifac = 1 , nfabor
      flux                =      bmasfl(ifac)
      bmasfl(ifac) = 2.d0*bmasfl(ifac)            &
           - propfb(ifac,iflmba)
      propfb(ifac,iflmba) = flux
    enddo
  endif

!     Les valeurs courantes ecrasent les valeurs anterieures
!       en cas d'extrapolation en temps (theta > 0)
!       Pour RHO, on le fait en double si ICALHY = 1 (et sur NCELET)
!     Au debut du calcul les flux nouveau et ancien ont ete initialises inivar
  if(iroext.gt.0) then
    ipcrom = ipproc(irom  )
    ipcroa = ipproc(iroma )
    do iel = 1, ncelet
      propce(iel,ipcroa) = propce(iel,ipcrom)
    enddo
    ipbrom = ipprob(irom  )
    ipbroa = ipprob(iroma )
    do ifac = 1, nfabor
      propfb(ifac,ipbroa) = propfb(ifac,ipbrom)
    enddo
  endif
  if(iviext.gt.0) then
    ipcvis = ipproc(iviscl)
    ipcvst = ipproc(ivisct)
    ipcvsa = ipproc(ivisla)
    ipcvta = ipproc(ivista)
    do iel = 1, ncel
      propce(iel,ipcvsa) = propce(iel,ipcvis)
      propce(iel,ipcvta) = propce(iel,ipcvst)
    enddo
  endif
  if(icpext.gt.0) then
    if(icp   .gt.0) then
      iicp   = ipproc(icp   )
      iicpa  = ipproc(icpa  )
      do iel = 1, ncel
        propce(iel,iicpa ) = propce(iel,iicp  )
      enddo
    endif
  endif

!     Remarque : si on faisant cette operation pour tous les
!       scalaires, on la ferait plusieurs fois pour les scalaires
!       ayant une variance
  if (nscal.ge.1) then
    do iscal = 1, nscal
      if(ivisls(iscal).gt.0.and.iscavr(iscal).le.0) then
        if(ivsext(iscal).gt.0) then
          ipcvsl = ipproc(ivisls(iscal))
          ipcvsa = ipproc(ivissa(iscal))
          do iel = 1, ncel
            propce(iel,ipcvsa) = propce(iel,ipcvsl)
          enddo
        endif
      endif
    enddo
  endif

  return


!===============================================================================
! 2.  JUSTE APRES PHYVAR (ET DONC AVANT NAVSTO)
!===============================================================================


elseif(iappel.eq.2) then

! 2.1 MISE A JOUR DES VALEURS ANCIENNES
! =====================================

!     On passe ici dans le cas ou l'on suspecte que la valeur portee
!       par les variables "anciennes" n'est pas satisfaisante pour
!       extrapoler.
!     On passe ici au premier pas de temps et lorsque le fichier suite
!       ne comportait pas grandeur requise.

  if(initro.ne.1) then
    initro = 1
    if(iroext.gt.0) then
      ipcrom = ipproc(irom  )
      ipcroa = ipproc(iroma )
      do iel = 1, ncelet
        propce(iel,ipcroa) = propce(iel,ipcrom)
      enddo
      ipbrom = ipprob(irom  )
      ipbroa = ipprob(iroma )
      do ifac = 1, nfabor
        propfb(ifac,ipbroa) = propfb(ifac,ipbrom)
      enddo
    endif
  endif
  if(initvi.ne.1) then
    initvi = 1
    if(iviext.gt.0) then
      ipcvis = ipproc(iviscl)
      ipcvst = ipproc(ivisct)
      ipcvsa = ipproc(ivisla)
      ipcvta = ipproc(ivista)
      do iel = 1, ncel
        propce(iel,ipcvsa) = propce(iel,ipcvis)
        propce(iel,ipcvta) = propce(iel,ipcvst)
      enddo
    endif
  endif
  if(initcp.ne.1) then
    initcp = 1
    if(icpext.gt.0) then
      if(icp   .gt.0) then
        iicp   = ipproc(icp   )
        iicpa  = ipproc(icpa  )
        do iel = 1, ncel
          propce(iel,iicpa ) = propce(iel,iicp  )
        enddo
      endif
    endif
  endif

!     Remarque : si on faisant cette operation pour tous les
!       scalaires, on la ferait plusieurs fois pour les scalaires
!       ayant une variance
  if (nscal.ge.1) then
    do iscal = 1, nscal
      if(initvs(iscal).ne.1) then
        initvs(iscal) = 1
        if(ivisls(iscal).gt.0.and.iscavr(iscal).le.0) then
          if(ivsext(iscal).gt.0) then
            ipcvsl = ipproc(ivisls(iscal))
            ipcvsa = ipproc(ivissa(iscal))
            do iel = 1, ncel
              propce(iel,ipcvsa) = propce(iel,ipcvsl)
            enddo
          endif
        endif
      endif
    enddo
  endif


! 2.2 EXTRAPOLATION DES NOUVELLES VALEURS
! =======================================

! --- Extrapolation de la viscosite et de la masse volumique dans le cas d'un
!     theta schema
!     A partir de Fn-1 et Fn on calcule Fn+theta
!     On conserve les nouvelles valeurs dans l'ancien tableau pour
!     retablir en fin de pas de temps

!     Le calcul pour Rho est fait sur NCELET afin d'economiser un echange.

  if(iroext.gt.0) then
    ipcrom = ipproc(irom  )
    ipcroa = ipproc(iroma )
    theta  = thetro
    do iel = 1, ncelet
      xmasvo = propce(iel,ipcrom)
      propce(iel,ipcrom) = (1.d0+theta) * propce(iel,ipcrom)    &
           -       theta  * propce(iel,ipcroa)
      propce(iel,ipcroa) = xmasvo
    enddo
    ipbrom = ipprob(irom  )
    ipbroa = ipprob(iroma )
    do ifac = 1, nfabor
      xmasvo = propfb(ifac,ipbrom)
      propfb(ifac,ipbrom) = (1.d0+theta) * propfb(ifac,ipbrom)  &
           -       theta  * propfb(ifac,ipbroa)
      propfb(ifac,ipbroa) = xmasvo
    enddo
  endif
  if(iviext.gt.0) then
    ipcvis = ipproc(iviscl)
    ipcvst = ipproc(ivisct)
    ipcvsa = ipproc(ivisla)
    ipcvta = ipproc(ivista)
    theta  = thetvi
    do iel = 1, ncel
      viscos = propce(iel,ipcvis)
      propce(iel,ipcvis) = (1.d0+theta) * propce(iel,ipcvis)    &
           -       theta  * propce(iel,ipcvsa)
      propce(iel,ipcvsa) = viscos
      viscos = propce(iel,ipcvst)
      propce(iel,ipcvst) = (1.d0+theta) * propce(iel,ipcvst)    &
           -       theta  * propce(iel,ipcvta)
      propce(iel,ipcvta) = viscos
    enddo
  endif
  if(icpext.gt.0) then
    if(icp.gt.0) then
      iicp   = ipproc(icp   )
      iicpa  = ipproc(icpa  )
      theta  = thetcp
      do iel = 1, ncel
        varcp  = propce(iel,iicp  )
        propce(iel,iicp  ) = (1.d0+theta) * propce(iel,iicp  )  &
             -       theta  * propce(iel,iicpa )
        propce(iel,iicpa ) = varcp
      enddo
    endif
  endif

!     Remarque : si on faisant cette operation pour tous les
!       scalaires, on la ferait plusieurs fois pour les scalaires
!       ayant une variance ET CE SERAIT FAUX
  if (nscal.ge.1) then
    do iscal = 1, nscal
      if(ivisls(iscal).gt.0.and.iscavr(iscal).le.0) then
        if(ivsext(iscal).gt.0) then
          theta  = thetvs(iscal)
          ipcvsl = ipproc(ivisls(iscal))
          ipcvsa = ipproc(ivissa(iscal))
          if(ipcvsl.gt.0) then
            do iel = 1, ncel
              viscos = propce(iel,ipcvsl)
              propce(iel,ipcvsl) = (1.d0+theta)*propce(iel,ipcvsl)&
                                 -       theta *propce(iel,ipcvsa)
              propce(iel,ipcvsa) = viscos
            enddo
          endif
        endif
      endif
    enddo
  endif

  return



!===============================================================================
! 3.  JUSTE APRES NAVSTO, DANS LES BOUCLES U/P ET ALE
!===============================================================================

elseif(iappel.eq.3) then

!     On traite ici le flux de masse uniquement
!        On suppose qu'il n'y en a qu'un seul.

!     si istmpf = 1 : standard : on ne fait rien
!     si istmpf = 2 : ordre 2 (thetfl > 0 : = 0.5)
!       on calcule f(n+theta) par interpolation a partir
!       de F_(n-1+theta) et f(n+1) et on le met dans imasfl
!     si istmpf = 0 : explicite (thetfl = 0) : on remet F(n) dans
!       imasfl sauf a la derniere iteration (un traitement
!       complementaire sera fait en iappel=4)

!     Dans le cas ou on itere sur navsto, on passe ici
!     - a toutes les iterations si ISTMPF.NE.0
!     - a toutes les iterations sauf la derniere si ISTMPF.EQ.0

  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_s(iflmas, imasfl)
  call field_get_val_s(iflmab, bmasfl)

  if(istmpf.eq.2) then
    iflmsa = ipprof(ifluaa(iu))
    iflmba = ipprob(ifluaa(iu))
    theta  = thetfl
    aa = 1.d0/(2.d0-theta)
    bb = (1.d0-theta)/(2.d0-theta)
    do ifac = 1 , nfac
      imasfl(ifac) = aa * imasfl(ifac)            &
           + bb * propfa(ifac,iflmsa)
    enddo
    do ifac = 1 , nfabor
      bmasfl(ifac) = aa * bmasfl(ifac)            &
           + bb * propfb(ifac,iflmba)
    enddo
  elseif(istmpf.eq.0) then
    iflmsa = ipprof(ifluaa(iu))
    iflmba = ipprob(ifluaa(iu))
    do ifac = 1 , nfac
      imasfl(ifac) = propfa(ifac,iflmsa)
    enddo
    do ifac = 1 , nfabor
      bmasfl(ifac) = propfb(ifac,iflmba)
    enddo
  endif

  return

!===============================================================================
! 3.  JUSTE APRES NAVSTO (ET STRDEP), HORS LES BOUCLES U/P ET ALE
!===============================================================================

elseif(iappel.eq.4) then

!     On traite ici le flux de masse uniquement
!        On suppose qu'il n'y en a qu'un seul.

!     Si istmpf = 1 : standard : on ne fait rien
!     Si istmpf = 2 : ordre 2 (thetfl > 0 : = 0.5)
!       on calcule F(n+theta) par interpolation a partir
!       de F_(n-1+theta) et F(n+1) et on le met dans imasfl)
!     Si istmpf = 0 : explicite (thetfl = 0)
!       On sauvegarde F_(n+1) dans propfa(1,iflmsa),mais on continue
!       les calculs avec F_(n) mis dans imasfl

!     On retablira au dernier appel de schtmp pour istmpf = 0

!     Dans le cas ou on itere sur navsto, on passe ici
!       - a toutes les iterations            si istmpf.ne.0
!       - uniquement a la derniere iteration si istmpf.eq.0
!         (ce faisant, a partir de la deuxieme sous-iteration,
!          le calcul sera fait avec F(n+1) et plus F(n), mais on
!          suppose que l'utilisateur a choisi de faire des sous-iter
!          aussi pour impliciter le flux de masse)

  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_s(iflmas, imasfl)
  call field_get_val_s(iflmab, bmasfl)

  if(istmpf.eq.2) then
    iflmsa = ipprof(ifluaa(iu))
    iflmba = ipprob(ifluaa(iu))
    theta  = thetfl
    aa = 1.d0/(2.d0-theta)
    bb = (1.d0-theta)/(2.d0-theta)
    do ifac = 1 , nfac
      imasfl(ifac) = aa * imasfl(ifac)            &
           + bb * propfa(ifac,iflmsa)
    enddo
    do ifac = 1 , nfabor
      bmasfl(ifac) = aa * bmasfl(ifac)            &
           + bb * propfb(ifac,iflmba)
    enddo
  elseif(istmpf.eq.0) then
    iflmsa = ipprof(ifluaa(iu))
    iflmba = ipprob(ifluaa(iu))
    do ifac = 1 , nfac
      flux = imasfl(ifac)
      imasfl(ifac) = propfa(ifac,iflmsa)
      propfa(ifac,iflmsa) = flux
    enddo
    do ifac = 1 , nfabor
      flux = bmasfl(ifac)
      bmasfl(ifac) = propfb(ifac,iflmba)
      propfb(ifac,iflmba) = flux
    enddo
  endif

  return

!===============================================================================
! 3.  JUSTE APRES SCALAI
!===============================================================================

elseif(iappel.eq.5) then

! 3.1 RETABLISSEMENT POUR LE FLUX DE MASSE
! ========================================

!     On corrige les manipulations sur le flux de masse faites dans
!       l'appel precedent afin d'etre pret pour le pas de temps suivant.

!     Si istmpf = 1 : standard : on ne fait rien
!     Si istmpf = 2 : ordre 2 (thetfl > 0 : = 0.5) : on ne fait rien
!     Si istmpf = 0 : explicite (thetfl = 0)
!       on remet F_(n+1) (stocke dans propfa(1,iflmsa)) dans imasfl
!       de sorte que les deux flux contiennent la meme chose

  if(istmpf.eq.0) then
    call field_get_key_int(ivarfl(iu), kimasf, iflmas)
    call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
    call field_get_val_s(iflmas, imasfl)
    call field_get_val_s(iflmab, bmasfl)
    iflmsa = ipprof(ifluaa(iu))
    iflmba = ipprob(ifluaa(iu))
    do ifac = 1 , nfac
      imasfl(ifac) = propfa(ifac,iflmsa)
    enddo
    do ifac = 1 , nfabor
      bmasfl(ifac) = propfb(ifac,iflmba)
    enddo
  endif

! 3.1 RETABLISSEMENT POUR LES PROPRIETES PHYSIQUES
! ================================================

!     Le calcul pour Rho est fait sur NCELET afin d'economiser un echange.

  if(iroext.gt.0) then
    ipcrom = ipproc(irom  )
    ipcroa = ipproc(iroma )
    do iel = 1, ncelet
      propce(iel,ipcrom) = propce(iel,ipcroa)
    enddo
    ipbrom = ipprob(irom  )
    ipbroa = ipprob(iroma )
    do ifac = 1, nfabor
      propfb(ifac,ipbrom) = propfb(ifac,ipbroa)
    enddo
  endif
  if(iviext.gt.0) then
    ipcvis = ipproc(iviscl)
    ipcvst = ipproc(ivisct)
    ipcvsa = ipproc(ivisla)
    ipcvta = ipproc(ivista)
    do iel = 1, ncel
      propce(iel,ipcvis) = propce(iel,ipcvsa)
      propce(iel,ipcvst) = propce(iel,ipcvta)
    enddo
  endif
  if(icpext.gt.0) then
    if(icp.gt.0) then
      iicp   = ipproc(icp   )
      iicpa  = ipproc(icpa  )
      do iel = 1, ncel
        propce(iel,iicp  ) = propce(iel,iicpa )
      enddo
    endif
  endif

!     Remarque : si on faisant cette operation pour tous les
!       scalaires, on la ferait plusieurs fois pour les scalaires
!       ayant une variance
  if (nscal.ge.1) then
    do iscal = 1, nscal
      if(ivisls(iscal).gt.0.and.iscavr(iscal).le.0) then
        if(ivsext(iscal).gt.0) then
          ipcvsl = ipproc(ivisls(iscal))
          ipcvsa = ipproc(ivissa(iscal))
          do iel = 1, ncel
            propce(iel,ipcvsl) = propce(iel,ipcvsa)
          enddo
        endif
      endif
    enddo
  endif

  return

endif

!===============================================================================

!--------
! FORMATS
!--------

!----
! FIN
!----

end subroutine
