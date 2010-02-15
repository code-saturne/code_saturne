!-------------------------------------------------------------------------------

!VERS


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

subroutine uslast &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , itepa  ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct ,          &
   statis , stativ ,                                              &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!       SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!       -----------------------------------

!    SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)

!    MODIFICATIONS UTILSATEUR SUR LES VARIABLES EN FIN D'ITERATION
!    LAGRANGIENNES ET CALCUL DES STATISTIQUES UTILISATEUR
!    SUPPLEMENTAIRES SUR LES PARTICULES

!   POUR LES STATISTIQUES UTILISATEUR SUPPLEMENTAIRES,
!   ON RAPPELLE QUE :

!   ISTTIO = 0 : calcul instationnaire pour le lagrangien
!          = 1 : calcul stationnaire   pour le lagrangien

!   ISTALA : calcul statistiques       si  >= 1 sinon pas de stat

!   ISUIST : suite calcul statistiques si  >= 1 sinon pas de stat

!   IDSTNT : Numero du pas de temps pour debut statistque

!   NSTIST : iteration Lagrangienne du debut calcul stationnaire

!   NPST   : nombre d'iterations de calcul de stat stationnaires

!   NPSTT  : nombre d'iterations total des stats depuis le debut
!            du calcul, partie instationnaire comprise

!   TSTAT  : Temps physique d'enregistrement des stats volumiques
!            stationnaires
!            (en instationnaire TSTAT=DTP le pas de temps Lagrangien)

!-------------------------------------------------------------------------------
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
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfac+1)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (lndfac)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!  (nfabor+1)      !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (lndfbr  )     !    !     !  (optionnel)                                   !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
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
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
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
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa(nbpmax,     ! tr ! <-- ! caracteristiques des particules                !
!       nvep)      !    !     !  aux particules (poids, ...)                   !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpmax,2)      !    !     !                                                !
! statis           ! tr ! <-- ! cumul pour les moyennes des                    !
!(ncelet,nvlsta    !    !     !   statistiques volumiques                      !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! w1..w3(ncelet    ! tr ! --- ! tableaux de travail                            !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "numvar.h"
include "cstnum.h"
include "optcal.h"
include "pointe.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"
include "cstphy.h"
include "ppppar.h"
include "ppthch.h"
include "cpincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse
integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1) , nodfbr(lndfbr)
integer          itepa(nbpmax,nivep)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia , idebra
integer          ifinia, ifinra
integer          npt ,  iel , iphas

integer          ivf , ivff , itabvr , iflu , icla

! Local variables UTILISATEUR

integer          nxlist
parameter       (nxlist=100)

integer          iplan
integer          ii, ind, il
integer          inoeud, irang0, indic
integer          ist(6)

double precision zz(4), zzz(8), xlist(nxlist,8), xyzpt(3)

character        name(8)*4

double precision debm(4)
save             debm

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START

if(istala.eq.1 .and. iplas.ge.idstnt .and. nvlsts.gt.0) then

!     Si l'on passe ici, il faut que l'utilisateur complete
!       l'exemple ci-dessous et l'adapte...

  if(1.eq.1) then
    write(nfecra,9000)nvlsts
    call csexit (1)
  endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS LE MODULE LAGRANGIEN             ',/,&
'@    =========                                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uslast DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Des variables statistiques supplementaires ont ete        ',/,&
'@    demandees dans uslag1 (NVLSTS=',   I10,')               ',/,&
'@  Le sous-programme uslast doit etre complete pour preciser ',/,&
'@    le  calcul de leur cumul.                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

else

!     On entre toujours dans ce sous programme en lagrangien,
!       si on ne souhaite rien y faire, on sort immediatement.

  return

endif

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INITIALISATION
!===============================================================================

iphas = ilphas

!===============================================================================
! 2 - CALCUL DES STATISTIQUES PARTICULAIRES UTILISATEURS
!===============================================================================

!   D'une facon generale, dans cette routine on realise les cumuls
!   de la quantite dont on souhaite faire les statistiques.
!   La moyenne et la variance sont calculees dans la routine
!   USLAEN.F. Ce calcul est le plus souvent obtenu par division
!   des cumuls soit par le temps du cumul stationnaire contenu dans
!   la variable TSTAT, soit par le nombre de particules en poids
!   statistiques. Cette division est appliquee pour chaque ecriture
!   dans le listing et pour les sorties post-processing.

!   Cet exemple est desactive et doit etre adapte au cas traite

if (1.eq.0) then

 if(istala.eq.1 .and. iplas.ge.idstnt .and. nvlsts.gt.0) then

  do npt = 1,nbpart

    if( itepa(npt,jisor).gt.0 ) then

      iel = itepa(npt,jisor)

! -------------------------------------------------
! EXEMPLE 1 : Cumul pour la concentration massique
! -------------------------------------------------

      statis(iel,ilvu(1)) = statis(iel,ilvu(1))                   &
        + tepa(npt,jrpoi) *ettp(npt,jmp)

      stativ(iel,ilvu(1)) = stativ(iel,ilvu(1))                   &
        + tepa(npt,jrpoi) *ettp(npt,jmp) *ettp(npt,jmp)

    endif

  enddo

 endif

endif

!===============================================================================
! 3 - CALCUL UTILISATEUR DU DEBIT MASSIQUE DE PARTICULES SUR 4 PLANS
!===============================================================================

!   Cet exemple est desactive et doit etre adapte au cas traite

if (1.eq.0) then

  zz(1) = 0.1d0
  zz(2) = 0.15d0
  zz(3) = 0.20d0
  zz(4) = 0.25d0

!   Si on est en instationnaire, ou si le debut des stat stationnaires
!   n'est pas encore atteint, toutes les statistiques sont remises a
!   zero a chaque pas de temps avant d'entrer dans ce sous-programme.

  if(isttio.eq.0 .or. npstt.le.nstist) then
    do iplan = 1,4
      debm(iplan) = 0.d0
    enddo
  endif

  do iplan = 1,4

    do npt = 1,nbpart

      if(itepa(npt,jisor).gt.0) then

        iel = itepa(npt,jisor)

        if( ettp(npt,jxp).gt.zz(iplan) .and.                      &
            ettpa(npt,jxp).le.zz(iplan)      ) then
          debm(iplan) = debm(iplan) +tepa(npt,jrpoi)*ettp(npt,jmp)
        endif

      endif

    enddo
  enddo

  do iplan = 1,4
    write(nfecra,1001)iplan,debm(iplan)/tstat
  enddo

 1001   format(' Debit massique particulaire en Z(',I10,') : ',E14.5)

endif


!===============================================================================
! 4 - EXTRACTION DE STATISTIQUES VOLUMIQUES EN FIN DE CALCUL
!===============================================================================

!   Cet exemple est desactive et doit etre adapte au cas traite

if (1.eq.0) then

  if(ntcabs.eq.ntmabs) then

    zzz(1) = 0.005d0
    zzz(2) = 0.025d0
    zzz(3) = 0.050d0
    zzz(4) = 0.075d0
    zzz(5) = 0.100d0
    zzz(6) = 0.150d0
    zzz(7) = 0.200d0
    zzz(8) = 0.250d0

    NAME(1) = 'XB01'
    NAME(2) = 'XB05'
    NAME(3) = 'XB10'
    NAME(4) = 'XB15'
    NAME(5) = 'XB20'
    NAME(6) = 'XB30'
    NAME(7) = 'XB40'
    NAME(8) = 'XB50'

    ist(1) = ilvx
    ist(2) = ilvz
    ist(3) = ilfv
    ist(4) = ilpd

    npts = nxlist

    ifinia = idebia
    itabvr = idebra
    ifinra = itabvr + ncelet
    CALL RASIZE('USLAST',IFINRA)
    !==========

    do iplan = 1,8

!     Pour le fichier ci-dessous :
!       l'utilisateur verifiera qu'il n'a pas laisse ouverte l'unite
!       IMPUSR(1), dans un autre sous-programme utilisateur
      OPEN(FILE=NAME(IPLAN),UNIT=IMPUSR(1),FORM='formatted')

      xyzpt(1) = zzz(iplan)

      do ivf = 1,4

        ivff = ist(ivf)
        icla = 0
        iflu = 0

        call uslaen                                               &
        !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivff   , ivff   , ivff   , iflu   , ilpd   , icla   ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statis , stativ , ra(itabvr) ,               &
   rdevel , rtuser , ra     )

        ind = 0
        do ii = 1, npts

          xyzpt(2) = 0.d0
          xyzpt(3) = float(ii-1)/float(npts-1)*150.d-3

          call findpt                                             &
          !==========
          (ncelet, ncel, xyzcen,                                  &
           xyzpt(1), xyzpt(2), xyzpt(3), inoeud, irang0)

          indic = ituser(inoeud)
          ituser(inoeud) = 1
          if (indic.eq.1) then
            ind = ind +1
            xlist(ind,1) = xyzcen(1,inoeud)
            xlist(ind,2) = xyzcen(3,inoeud) * (1.d3 / 5.d0)
            xlist(ind,ivf+2) = ra(itabvr+inoeud-1)
          endif
        enddo
      enddo

      do il = 1, ind
        WRITE (IMPUSR(1),'(8E13.5)') (XLIST(IL,II), II=1,6)
      enddo

      close(impusr(1))

    enddo

  endif

endif



!===============================================================================

!====
! FIN
!====

return

end subroutine
