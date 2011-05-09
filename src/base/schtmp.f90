!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

subroutine schtmp &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , iappel ,                            &
   isostd ,                                                       &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! iappel           ! e  ! <-- ! numero de l'appel (avant ou apres              !
!                  !    !     ! phyvar                                         !
! isostd           ! te ! <-- ! indicateur de sortie standard                  !
!    (nfabor+1)    !    !     !  +numero de la face de reference               !
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

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas  , iappel

integer          isostd(nfabor+1,nphas)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ra(*)

! Local variables

integer          idebia , idebra
integer          iel    , ifac   , iscal
integer          iphas  , iuiph
integer          ipcrom , ipcroa
integer          ipbrom , ipbroa
integer          iflmas , iflmab , iflmba, iflmsa
integer          ipcvis , ipcvst
integer          ipcvsa , ipcvta , ipcvsl
integer          iicp   , iicpa
double precision flux   , theta  , aa, bb, viscos, xmasvo, varcp

!===============================================================================


!===============================================================================
! 0. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1.  AU TOUT DEBUT DE LA BOUCLE EN TEMPS
!===============================================================================

if(iappel.eq.1) then

! --- Application du schema en temps sur le flux de masse
!     Soit F le flux de masse
!       - si ISTMPF = 2 (schema non standard, theta>0 en fait = 0.5)
!         PROPFA(1,IFLMAS) contient F_(n-2+theta) et on y met F(n-1+theta)
!         PROPFA(1,IFLMSA) contient F_(n-1+theta) et
!                                    on y met une extrapolation en n+theta
!       - si ISTMPF = 1 (schema non standard, theta=0)
!         PROPFA(1,IFLMAS) contient deja F_n et
!         PROPFA(1,IFLMSA) n'est pas utilise : on ne fait rien
!       - sinon : ISTMPF = 0 (schema standard, theta= -999)
!         PROPFA(1,IFLMAS) et PROPFA(1,IFLMSA) contiennent tous les deux F(n)
!                                            : on ne fait rien


!     Ordre 2 en temps pour le flux (theta = 0.5) a entrer dans navsto
!       Flux convectif = 2F(n-1+theta)-F(n-2+theta)
!       Flux conserve  =  F(n-1+theta)
!     Au premier pas de temps, l'ancien a ete initialise dans inivar (a 0)
!       en suite de calcul, les deux ont ete relus.

  do iphas = 1, nphas
    if(istmpf.eq.2) then
      iuiph  = iu
      iflmas = ipprof(ifluma(iuiph))
      iflmab = ipprob(ifluma(iuiph))
      iflmsa = ipprof(ifluaa(iuiph))
      iflmba = ipprob(ifluaa(iuiph))
      do ifac = 1 , nfac
        flux                =      propfa(ifac,iflmas)
        propfa(ifac,iflmas) = 2.d0*propfa(ifac,iflmas)            &
                                 - propfa(ifac,iflmsa)
        propfa(ifac,iflmsa) = flux
      enddo
      do ifac = 1 , nfabor
        flux                =      propfb(ifac,iflmab)
        propfb(ifac,iflmab) = 2.d0*propfb(ifac,iflmab)            &
                                 - propfb(ifac,iflmba)
        propfb(ifac,iflmba) = flux
      enddo
    endif
  enddo

!     Les valeurs courantes ecrasent les valeurs anterieures
!       en cas d'extrapolation en temps (theta > 0)
!       Pour RHO, on le fait en double si ICALHY = 1 (et sur NCELET)
!     Au debut du calcul les flux nouveau et ancien ont ete initialises inivar
  do iphas = 1, nphas
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
  enddo
  do iphas = 1, nphas
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
  enddo
  do iphas = 1, nphas
    if(icpext.gt.0) then
      if(icp   .gt.0) then
        iicp   = ipproc(icp   )
        iicpa  = ipproc(icpa  )
        do iel = 1, ncel
          propce(iel,iicpa ) = propce(iel,iicp  )
        enddo
      endif
    endif
  enddo

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

  do iphas = 1, nphas
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
  enddo
  do iphas = 1, nphas
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
  enddo
  do iphas = 1, nphas
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
  enddo

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

  do iphas = 1, nphas
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
  enddo

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

!     Si ISTMPF = 1 : standard : on ne fait rien
!     Si ISTMPF = 2 : ordre 2 (THETFL > 0 : = 0.5)
!       on calcule F(n+theta) par interpolation a partir
!       de F_(n-1+theta) et F(n+1) et on le met dans PROPFA(1,IFLMAS)
!     Si ISTMPF = 0 : explicite (THETFL = 0) : on remet F(n) dans
!       PROPFA(1,IFLMAS) sauf a la derniere iteration (un traitement
!       complementaire sera fait en IAPPEL=4)

!     Dans le cas ou on itere sur navsto, on passe ici
!     - a toutes les iterations si ISTMPF.NE.0
!     - a toutes les iterations sauf la derniere si ISTMPF.EQ.0

  do iphas = 1, nphas
    iuiph  = iu
    iflmas = ipprof(ifluma(iuiph))
    iflmab = ipprob(ifluma(iuiph))
    if(istmpf.eq.2) then
      iflmsa = ipprof(ifluaa(iuiph))
      iflmba = ipprob(ifluaa(iuiph))
      theta  = thetfl
      aa = 1.d0/(2.d0-theta)
      bb = (1.d0-theta)/(2.d0-theta)
      do ifac = 1 , nfac
        propfa(ifac,iflmas) = aa * propfa(ifac,iflmas)            &
                            + bb * propfa(ifac,iflmsa)
      enddo
      do ifac = 1 , nfabor
        propfb(ifac,iflmab) = aa * propfb(ifac,iflmab)            &
                            + bb * propfb(ifac,iflmba)
      enddo
    elseif(istmpf.eq.0) then
      iflmsa = ipprof(ifluaa(iuiph))
      iflmba = ipprob(ifluaa(iuiph))
      do ifac = 1 , nfac
        propfa(ifac,iflmas) = propfa(ifac,iflmsa)
      enddo
      do ifac = 1 , nfabor
        propfb(ifac,iflmab) = propfb(ifac,iflmba)
      enddo
    endif
  enddo

  return

!===============================================================================
! 3.  JUSTE APRES NAVSTO (ET STRDEP), HORS LES BOUCLES U/P ET ALE
!===============================================================================

elseif(iappel.eq.4) then

!     On traite ici le flux de masse uniquement
!        On suppose qu'il n'y en a qu'un seul.

!     Si ISTMPF = 1 : standard : on ne fait rien
!     Si ISTMPF = 2 : ordre 2 (THETFL > 0 : = 0.5)
!       on calcule F(n+theta) par interpolation a partir
!       de F_(n-1+theta) et F(n+1) et on le met dans PROPFA(1,IFLMAS)
!     Si ISTMPF = 0 : explicite (THETFL = 0)
!       On sauvegarde F_(n+1) dans PROPFA(1,IFLMSA),mais on continue
!       les calculs avec F_(n) mis dans PROPFA(1,IFLMAS)

!     On retablira au dernier appel de schtmp pour ISTMPF = 0

!     Dans le cas ou on itere sur navsto, on passe ici
!       - a toutes les iterations            si ISTMPF.NE.0
!       - uniquement a la derniere iteration si ISTMPF.EQ.0
!         (ce faisant, a partir de la deuxieme sous-iteration,
!          le calcul sera fait avec F(n+1) et plus F(n), mais on
!          suppose que l'utilisateur a choisi de faire des sous-iter
!          aussi pour impliciter le flux de masse)

  do iphas = 1, nphas
    iuiph  = iu
    iflmas = ipprof(ifluma(iuiph))
    iflmab = ipprob(ifluma(iuiph))
    if(istmpf.eq.2) then
      iflmsa = ipprof(ifluaa(iuiph))
      iflmba = ipprob(ifluaa(iuiph))
      theta  = thetfl
      aa = 1.d0/(2.d0-theta)
      bb = (1.d0-theta)/(2.d0-theta)
      do ifac = 1 , nfac
        propfa(ifac,iflmas) = aa * propfa(ifac,iflmas)            &
                            + bb * propfa(ifac,iflmsa)
      enddo
      do ifac = 1 , nfabor
        propfb(ifac,iflmab) = aa * propfb(ifac,iflmab)            &
                            + bb * propfb(ifac,iflmba)
      enddo
    elseif(istmpf.eq.0) then
      iflmsa = ipprof(ifluaa(iuiph))
      iflmba = ipprob(ifluaa(iuiph))
      do ifac = 1 , nfac
        flux = propfa(ifac,iflmas)
        propfa(ifac,iflmas) = propfa(ifac,iflmsa)
        propfa(ifac,iflmsa) = flux
      enddo
      do ifac = 1 , nfabor
        flux = propfb(ifac,iflmab)
        propfb(ifac,iflmab) = propfb(ifac,iflmba)
        propfb(ifac,iflmba) = flux
      enddo
    endif
  enddo

  return

!===============================================================================
! 3.  JUSTE APRES SCALAI
!===============================================================================

elseif(iappel.eq.5) then

! 3.1 RETABLISSEMENT POUR LE FLUX DE MASSE
! ========================================

!     On corrige les manipulations sur le flux de masse faites dans
!       l'appel precedent afin d'etre pret pour le pas de temps suivant.

!     Si ISTMPF = 1 : standard : on ne fait rien
!     Si ISTMPF = 2 : ordre 2 (THETFL > 0 : = 0.5) : on ne fait rien
!     Si ISTMPF = 0 : explicite (THETFL = 0)
!       on remet F_(n+1) (stocke dans PROPFA(1,IFLMSA)) dans PROPFA(1,IFLMAS)
!       de sorte que les deux flux contiennent la meme chose

  do iphas = 1, nphas
    if(istmpf.eq.0) then
      iuiph  = iu
      iflmas = ipprof(ifluma(iuiph))
      iflmab = ipprob(ifluma(iuiph))
      iflmsa = ipprof(ifluaa(iuiph))
      iflmba = ipprob(ifluaa(iuiph))
      do ifac = 1 , nfac
        propfa(ifac,iflmas) = propfa(ifac,iflmsa)
      enddo
      do ifac = 1 , nfabor
        propfb(ifac,iflmab) = propfb(ifac,iflmba)
      enddo
    endif
  enddo

! 3.1 RETABLISSEMENT POUR LES PROPRIETES PHYSIQUES
! ================================================

!     Le calcul pour Rho est fait sur NCELET afin d'economiser un echange.

  do iphas = 1, nphas
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
  enddo

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
