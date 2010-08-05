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

subroutine mtphyv &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ibrom  ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   ,                                     &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE MATISSE : REMPLISSAGE DES VARIABLES PHYSIQUES
!             COPIE DE LA ROUTINE UTILISATEUR USPHYV



! ATTENTION : (on conserve pour memoire les mises en garde de usphyv)
! =========


! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)


!  Il FAUT AVOIR PRECISE ICP(IPHAS) = 1
!     ==================
!    dans usini1 si on souhaite imposer une chaleur specifique
!    CP variable pour la phase IPHAS (sinon: ecrasement memoire).


!  Il FAUT AVOIR PRECISE IVISLS(Numero de scalaire) = 1
!     ==================
!     dans usini1 si on souhaite une diffusivite VISCLS variable
!     pour le scalaire considere (sinon: ecrasement memoire).




! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usini1 :
!             . la masse volumique (initialisee a RO0(IPHAS))
!             . la viscosite       (initialisee a VISCL0(IPHAS))
!      - dans usiniv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a la valeur donnee dans usiniv)

! On peut donner ici les lois de variation aux cellules
!     - de la masse volumique                      ROM    kg/m3
!         (et eventuellememt aux faces de bord     ROMB   kg/m3)
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)


! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! nphmx            ! e  ! <-- ! nphsmx                                         !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!   (nphmx   )     !    !     !                                                !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "cstnum.f90"
include "paramx.f90"
include "pointe.f90"
include "numvar.f90"
include "optcal.f90"
include "cstphy.f90"
include "entsor.f90"
include "parall.f90"
include "period.f90"
include "matiss.f90"


!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse , nphmx

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr), ibrom(nphmx)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          iel   , icoul , ifml
integer          iphas , ivarta, ivartc
integer          ipcrom, ipcvis, ipcvsl, ipcvsc, ipcvsp
double precision xrtp  , un0   , visct1
double precision trfmtk, ustmat, drfmat
double precision hray  , zeta  , f12   , cxy


!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

! --- Initialisation memoire
idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. REPERAGE DES VARIABLES
!===============================================================================

! --- Une seule phase
iphas  = 1

! --- Reperage des scalaires (température air et de peau des colis)
!     ISCALT est complete en retour de l'IHM

ivarta  = isca(iscalt(iphas))
ivartc  = isca(itpcmt)

! --- Rang de la masse volumique de la phase courante IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCROM

ipcrom = ipproc(irom(iphas))

! --- Rang de la viscosite dynamique de la phase courante IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCVIS

ipcvis = ipproc(iviscl(iphas))

! --- Rang des diffusivites des scalaires
!     dans PROPCE, prop. physiques au centre des elements       : IPCVS*
!     - IPCVSL : scalaire ITAAMT = T air ambiant (scalaire temperature)
!     - IPCVSC : scalaire ITPCMT = T PeauColis
!     - IPCVSP : scalaire ITPPMT = T PeauParoi

ipcvsl = ipproc(ivisls(iscalt(iphas)))
ipcvsc = ipproc(ivisls(itpcmt))
ipcvsp = ipproc(ivisls(itppmt))


!===============================================================================
! 2. GRANDEURS GEOMETRIQUES
!===============================================================================

! --- CXY et F12 (facteur de forme rayonnement)
cxy  = (ptrres + plgres) * 0.5d0
zeta = cxy/dmcont
f12  = 2.d0/pi*( (zeta**2 - 1.d0)**0.5d0                          &
     - zeta + asin(1.d0/zeta) )


!===============================================================================
! 3. MASSE VOLUMIQUE
!===============================================================================

! --- Masse volumique au centre des cellules

!     Utilisation de la loi des gaz parfait
!     (attention, il faut des temperatures en Kelvin)
trfmtk = trfmat+tkelvi
do iel = 1, ncel
  xrtp = rtp(iel,ivarta)
  propce(iel,ipcrom) = trfmtk*rrfmat /(xrtp + tkelvi)
enddo

!===============================================================================
! 4. VISCOSITE DYNAMIQUE ET DIFFUSIVITE DES SCALAIRES
!===============================================================================

! --- Remarque : par securite, on s'assure que la viscosite et la
!       diffusivite ne sont jamais nulles en ajoutant si besoins la
!       viscosite moleculaire de l'air.


do iel = 1, ncel


! --- Vitesse de reference locale : max de la vitesse instantanee et
!       locale et de la vitesse de reference calculee dans mtini1

  un0 = sqrt(rtpa(iel,iu(iphas))**2 + rtpa(iel,iv(iphas))**2      &
       + rtpa(iel,iw(iphas))**2)
  un0 = max(vitref,un0)


! --- Viscosite totale de reference : VISCL0 calculee dans mtini1
!       MAIS en convection naturelle, VITREF est modifie dans mttsns
!         il faut donc alors recalculer VISCL0 ici.
!       Par securite on le recalcule aussi en convection forcee.

!     La viscosite dynamique turbulente est modelisee sous la forme
!       mu_t = rho * u*_ref * D_ref, formule dans laquelle :
!     . u*_ref est la vitesse de frottement deduite de la vitesse de
!              reference VITREF calculee precedemment et d'une intensite
!              turbulente imposee RTURB0
!     . D_ref  est une distance de reference basee sur l'encombrement
!              du milieu

!     La valeur de u*_ref est calculee a partir de l'energie cinetique
!       turbulente k par u*_ref = Cmu**(1/4) k**(1/2). La valeur de k
!       se deduit de l'intensite turbulente I, supposee connue, par
!       k = 3/2 (I V_ref)**2 avec V_ref la vitesse de reference
!       VITREF calculee precedemment.

!     La valeur de D_ref est calculee par D_ref = 0.2(Pt - d) ou Pt est
!       le pas transverse du reseau et d le diametre des conteneurs.

!     Ainsi, mu_t = rho * u*_ref * L_ref
!            . u*_ref = Cmu**(1/4) * (3/2)**(1/2) * I * V_ref
!            . L_ref  = 0.2 * (Pt - d)

!     On ajoute la viscosite moleculaire XMUMAT a mu_t pour obtenir
!       la viscosite dynamique totale VISCL0

  ustmat = cmu**0.25d0 * sqrt(1.5d0) * (rturb0/100.d0) * vitref
  drfmat = 0.2d0 * (ptrres-dmcont)

  viscl0(iphas) = ro0(iphas) * ustmat * drfmat  + xmumat


! --- Viscosite totale VISCT1 calculee comme ci-dessus mais
!     . avec la vitesse de reference UN0 calculee ci-dessus
!     . avec la masse volumique instantanee et locale

!     La viscosite dynamique turbulente est modelisee sous la forme
!       mu_t = rho * u*_ref * D_ref, formule dans laquelle :
!     . u*_ref est la vitesse de frottement deduite de la vitesse de
!              reference UN0    calculee ci-dessus et d'une intensite
!              turbulente imposee RTURB0
!     . D_ref  est un distance de reference basee sur l'encombrement
!              du milieu

!     La valeur de u*_ref est calculee à partir de l'energie cinetique
!       turbulente k par u*_ref = Cmu**(1/4) k**(1/2). La valeur de k
!       se deduit de l'intensite turbulente I, supposee connue, par
!       k = 3/2 (I V_ref)**2 avec V_ref la vitesse de reference
!       UN0    calculee precedemment.

!     La valeur de D_ref est calculee par D_ref = 0.2(Pt - d) ou Pt est
!       le pas transverse du reseau et d le diametre des conteneurs.

!     Ainsi, mu_t = rho * u*_ref * L_ref
!            . u*_ref = Cmu**(1/4) * (3/2)**(1/2) * I * V_ref
!            . L_ref  = 0.2 * (Pt - d)

!     On ajoute la viscosite moleculaire XMUMAT a mu_t pour obtenir
!       la viscosite dynamique totale VISCL0

  ustmat = cmu**0.25d0 * sqrt(1.5d0) * (rturb0/100.d0) * un0
  drfmat = 0.2d0 * (ptrres-dmcont)

  visct1 = propce(iel,ipcrom) * ustmat * drfmat  + xmumat


! --- Reperage des elements selon les zones du maillage

!     Couleur de l'element courant
!     (hypotheses : une couleur et une seule, pas de groupe,
!                   et donc pas de boucle sur les familles)
  ifml  = ifmcel(iel   )
  icoul = iprfml(ifml,1)


! --- Zone de stockage
  if(icoul.eq.icmtst) then

!     Viscosite totale
    propce(iel,ipcvis) = max(viscl0(iphas),visct1)

!     Conductivite totale pour l'air (Nb de Prandlt total = 1.)
    propce(iel,ipcvsl) = propce(iel,ipcvis)

!     Diffusivite equivalente pour la temperature de peau des colis
!     . Au dessus de HRESO ou entreposage en alveole : zero
!     . En dessous de HRESO sans alveole : rayonnement
!       phi lineique = HRAY*RF*(CXY-DMCONT)*gradT
!       avec RF=DMCONT/CXY rapport surfacique
!     On ajoute la viscosite moleculaire
!       (pour eviter zero dans les conditions aux limites)
    if(xyzcen(3,iel).ge.hreso.or.ialveo.eq.1) then
      propce(iel,ipcvsc) = 0.d0
    else
      hray = stephn*f12*emicon*(rtp(iel,ivartc)+tkelvi)**3
      propce(iel,ipcvsc) = hray*(cxy-dmcont)*dmcont/cxy
    endif
    propce(iel,ipcvsc) = propce(iel,ipcvsc) + xmumat

!     Diffusivite equivalente pour la temperature de peau des murs
!       Pas de diffusion pour ce scalaire : traitement par terme source
!       (viscosite moleculaire pour eviter zero)
    propce(iel,ipcvsp) = xmumat


! --- Obstacles (registres amont)

  elseif(icoul.eq.icmtri)then
    propce(iel,ipcvis) = max(viscl0(iphas),visct1)
    propce(iel,ipcvsl) = xmumat
    propce(iel,ipcvsc) = xmumat
    propce(iel,ipcvsp) = xmumat

! --- Obstacles (registres aval)

  elseif(icoul.eq.icmtro) then
    propce(iel,ipcvis) = max(viscl0(iphas),visct1)
    propce(iel,ipcvsl) = xmumat
    propce(iel,ipcvsc) = xmumat
    propce(iel,ipcvsp) = xmumat

! --- Defaut (les autres zones)

  elseif(icoul.eq.icmtdf) then
    propce(iel,ipcvis) = max(viscl0(iphas),visct1)
    propce(iel,ipcvsl) = propce(iel,ipcvis)
    propce(iel,ipcvsc) = xmumat
    propce(iel,ipcvsp) = xmumat

! --- Cheminee d'alimentation (au dessus du convergent eventuel)
!     . on dope la diffusion pour homogeneiser T et U
!       et diminuer les risques de refoulement

  elseif(icoul.eq.icmtci) then
    propce(iel,ipcvis) = max(viscl0(iphas),20.d0*visct1)
    propce(iel,ipcvsl) = propce(iel,ipcvis)
    propce(iel,ipcvsc) = xmumat
    propce(iel,ipcvsp) = xmumat

! --- Cheminee d'evacuation (au dessus du convergent eventuel)
!     . on dope la diffusion pour homogeneiser T et U
!       et diminuer les risques de refoulement

  elseif(icoul.eq.icmtco) then
    propce(iel,ipcvis) = max(viscl0(iphas),20.d0*visct1)
    propce(iel,ipcvsl) = propce(iel,ipcvis)
    propce(iel,ipcvsc) = xmumat
    propce(iel,ipcvsp) = xmumat

! --- Sinon ...

  else
    propce(iel,ipcvis) = max(viscl0(iphas),visct1)
    propce(iel,ipcvsl) = propce(iel,ipcvis)
    propce(iel,ipcvsc) = xmumat
    propce(iel,ipcvsp) = xmumat
  endif

enddo


!===============================================================================
! FORMATS
!===============================================================================

!----
! FIN
!----

return
end subroutine
