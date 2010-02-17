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

subroutine usray5 &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  , iappel ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr , izfrdp ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cofrua , cofrub ,                                              &
   w1     , w2     , w3     , w4     , w5     ,  w6     ,         &
   tparoi , qincid , flunet , xlam   , epa    , eps     ,  ck   , &
   rdevel , rtuser ,                                              &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
!   -----------------------------------------

!   Ce sous-programme est appele 2 fois pour chaque phase IPHAS
!   pour laquelle il faut faire un calcul de rayonnement
!   semi-transparent.



!  1. PREMIER APPEL (IAPPEL = 1)
!  =============================



!    1.1 Conditions aux limites pour la luminance
!    --------------------------------------------

!       Il faut completer COFRUA qui fournit la luminance au bord
!         selon le type de frontiere (condition de Dirichlet).
!       La luminance est la puissance surfacique par unite
!         d'angle solide


!       Par exemple, on a


! 1/ Paroi grise : rayonnement isotrope.
!                                    4
!                      eps.sig.tparoi         (1-eps).qincid
!        cofrua   =    --------------    +    --------------
!                            pi                     pi
!  luminance de bord   emission propre         flux reflechi.

!     (eps=1 : paroi noire ; eps=0 : paroi reflechissante )

! Pour une paroi , la luminance (i.e l'energie rayonnee par la
! paroi ) comprend son emission propre  et l'energie reflechie.
!  (CF LE SOUS-PROGRAMME usray2)


! 2/ Milieu libre : luminance rentrante nulle

!        cofrua   =   0.D0

!    (si l'utilisateur a plus d'informations, il peut ameliorer
!     la situation)



!    L'exemple fourni ci-apres represente les conditions "par defaut"
!      et suffit generalement (si l'utilisateur a utilise
!      les types standard de conditions aux limites dans usclim)



!    1.2 Conditions aux limites pour le modele P-1
!    ---------------------------------------------




!  2. DEUXIEME APPEL (IAPPEL = 2)
!  =============================

!      La densite de flux net radiatif doit etre calculee
!        de maniere coherente avec les conditions aux limites
!        de la luminance. La densite de flux net radiatif est
!DONF         le bilan entre le rayonnement qu'une face de bord emet
!        (et non pas le rayonnement qu'elle reflechit) et celui
!        qu'elle absorbe.

!      L'exemple fourni est coherent avec l'exemple fourni pour
!        conditions aux limites sur la luminance au premier appel
!        et suffi donc en general.


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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
! iphas            ! i  ! <-- ! phase number                                   !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! itypfb           ! ia ! <-- ! boundary face types                            !
!  (nfabor, nphas) !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! izfrdp(nfabor    ! te ! <-- ! numero de zone pour les faces de bord          !
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
! cofrua,cofrub    ! tr ! --> ! conditions aux limites aux                     !
!(nfabor)          !    !     !    faces de bord pour la luminances            !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
! tparoi(nfabor    ! tr ! <-- ! temperature de paroi en kelvin                 !
! qincid(nfabor    ! tr ! <-- ! densite de flux radiatif aux bords             !
! flunet(nfabor    ! tr ! --> ! densite de flux net radiatif                   !
! ck (ncelet)      ! tr ! --> ! coefficient d'absorption du milieu             !
!                  !    !     ! (nul si transparent)                           !
! xlam(nfabor)     ! tr ! <-- ! coefficient de conductivite thermique          !
!                  !    !     ! des facettes de paroi (w/m/k)                  !
! epa (nfabor)     ! tr ! <-- ! epaisseur des facettes de paroi (m)            !
! eps (nfabor)     ! tr ! <-- ! emissivite des facettes de bord                !
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
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "pointe.h"
include "parall.h"
include "period.h"
include "ppppar.h"
include "ppthch.h"
include "cpincl.h"
include "ppincl.h"
include "radiat.h"
include "ihmpre.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , iphas  , iappel
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml) , itypfb(nfabor)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr),izfrdp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision cofrua(nfabor), cofrub(nfabor)

double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)

double precision tparoi(nfabor), qincid(nfabor)
double precision xlam(nfabor), epa(nfabor)
double precision eps(nfabor), flunet(nfabor)
double precision ck(ncelet)

double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! Local variables

integer          idebia , idebra , ifac, iok
double precision unspi, xit, distbf

!===============================================================================

!===============================================================================
! 0 - GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

! Indicateur d'arret (pour savoir si des faces ont ete oubliees)
iok = 0

unspi = 1.d0/pi

!===============================================================================
!  1. PREMIER APPEL
!  ================
!===============================================================================

if (iappel.eq.1) then


!===============================================================================
!  1.1 - CONDITIONS AUX LIMITES :
!        MODELE DOM : COFRUA CONTIENT LA LUMINANCE
!        MODELE P-1 : COFRUA ET COFRUB SONT A REMPLIR
!      LES EXEMPLES DONNES ICI SONT LES INITIALISATIONS FAITES
!      PAR DEFAUT ET SONT SUFFISANTES EN GENERAL
!===============================================================================




!      A - MODELE DOM
!      ^^^^^^^^^^^^^^




  if (iirayo.eq.1) then

    do ifac = 1,nfabor

!      1.1.1 - SYMETRIE :
!              ----------
!          REFLEXION TOTALE DU RAYONNEMENT ( EPS=0 )
!          -----------------------------------------

      if (itypfb(ifac).eq.isymet) then

        cofrua(ifac) = qincid(ifac) * unspi


!      1.1.2 - PAROIS 'FLUIDES' : LUMINANCES RENTRANTES "NULLES"
!              (ATTENTION LOGIQUE DIFFERENTE DU MODELE P-1)
!          -------------------------------------------------

      else if (itypfb(ifac).eq.ientre                             &
          .or. itypfb(ifac).eq.isolib) then

        cofrua(ifac) = epzero


!      1.1.3. - PAROIS 'SOLIDES' DE TEMPERATURE TPAROI ET D'EMISSIVITE EPS
!               ----------------------------------------------------------

      else if (itypfb(ifac).eq.iparoi                             &
          .or. itypfb(ifac).eq.iparug) then

        cofrua(ifac)  = eps(ifac)*stephn*(tparoi(ifac)**4)*unspi  &
                          + (1.d0-eps(ifac))* qincid(ifac)*unspi

      else

!      1.1.4 - SI DES FACES N'ONT PAS ETE TRAITEES, IL FAUT S'ARRETER
!              ------------------------------------------------------

!           ==============================================

!             CONSERVER IMPERATIVEMENT LE TEST D'ARRET

!           ==============================================

        write (nfecra,1000) ifac,izfrdp(ifac),itypfb(ifac)
        iok = iok + 1
      endif

    enddo




!   B - MODELE P-1
!   ^^^^^^^^^^^^^^





  else if (iirayo.eq.2) then

    do ifac = 1,nfabor

!      1.1.1 - SYMETRIE ET PAROI REFLECHISSANTE (EPS = 0) :
!              CONDITION DE FLUX NUL
!              ------------------------------------------

      if (itypfb(ifac).eq.isymet     .or.                         &
         ((itypfb(ifac).eq.iparoi.or.                             &
           itypfb(ifac).eq.iparug).and.eps(ifac).eq.0d0)) then

        cofrua(ifac) = 0.d0
        cofrub(ifac) = 1.d0


!      1.1.2 - PAROIS 'FLUIDES' : CONDITION DE FLUX NUL
!              (ATTENTION LOGIQUE DIFFERENTE DU MODELE DOM)
!              --------------------------------------------

      else if (itypfb(ifac).eq.ientre                             &
          .or. itypfb(ifac).eq.isolib) then

        cofrua(ifac) = 0.d0
        cofrub(ifac) = 1.d0


!      1.1.3 - PAROIS 'SOLIDES' DE TEMPERATURE TPAROI ET D'EMISSIVITE EPS
!              (EPS NON NUL)
!              ----------------------------------------------------------

      else if (itypfb(ifac).eq.iparoi .or.                        &
               itypfb(ifac).eq.iparug ) then

        distbf = ra(idistb-1+ifac)

        xit = 1.5d0 *distbf *ck(ifabor(ifac))                     &
            * (2.d0 /(2.d0-eps(ifac)) -1.d0)

        cofrub(ifac) = 1.d0 / (1.d0 + xit)
        cofrua(ifac) = xit * tparoi(ifac)**4 * cofrub(ifac)

      else

!      1.1.4 - SI DES FACES N'ONT PAS ETE TRAITEES, IL FAUT S'ARRETER
!              ------------------------------------------------------

!           ==============================================

!             CONSERVER IMPERATIVEMENT LE TEST D'ARRET

!           ==============================================

        write (nfecra,1000) ifac,izfrdp(ifac),itypfb(ifac)
      iok = iok + 1
    endif

  enddo

  endif

  if (iok.ne.0) then
    write (nfecra,1100) iphas
    call csexit (1)
    !==========
  endif

!===============================================================================
!  2 - DEUXIEME APPEL
!  ===================
!===============================================================================

else if (iappel.eq.2) then

!===============================================================================
!  2.1 - DENSITE DE FLUNET RADIATIF AUX DIFFERENTES FRONTIERES
!      L'EXEMPLE DONNE ICI EST L'INITIALISATION FAITE PAR DEFAUT
!===============================================================================

!    DANS LE CAS OU LES CONDITIONS A LA LIMITES CI-DESSUS
!      AURAIENT ETE MODIFIEES, IL EST NECESSAIRE DE MODIFIER
!      LA MANIERE DONT EST CALCULE LA DENSITE DE FLUX NET RADIATIF,
!      DE MANIERE COHERENTE.
!    LA REGLE EST LA SUIVANTE :
!      LA DENSITE DE FLUX NET EST UN BILAN ENTRE CE QU'UNE FACE
!      DE BORD EMET COMME RAYONNEMENT (ET NON CE QU'ELLE REFLECHIT)
!      ET CE QU'ELLE ABSORBE  (ORIENTATION DE LA NORMALE SORTANTE)
!      AINSI, SI UNE PAROI CHAUFFE LE FLUIDE, FLUNET < 0




  do ifac = 1,nfabor

    if (itypfb(ifac).eq.iparoi .or.                               &
        itypfb(ifac).eq.iparug) then

!      2.1.1 - PAROIS 'SOLIDES' DE TEMPERATURE TPAROI ET D'EMISSIVITE EPS
!              ----------------------------------------------------------

      flunet(ifac) =                                              &
      eps(ifac) *(qincid(ifac) - stephn*tparoi(ifac)**4)


!      2.1.2 - SYMETRIE :
!              ----------
!          REFLEXION TOTALE DU RAYONNEMENT ( FLUNET = 0 )
!          ----------------------------------------------
    else if (itypfb(ifac).eq.isymet) then

      flunet(ifac)= zero


!      2.1.3 - PAROIS 'FLUIDES'
!              ----------------

    else if (itypfb(ifac).eq.ientre                               &
        .or. itypfb(ifac).eq.isolib) then

      if (iirayo.eq.1) then

      flunet(ifac)= qincid(ifac) -pi*cofrua(ifac)

      else if (iirayo.eq.2) then

        flunet(ifac)= 0.d0

      endif


!      2.1.4 - SI DES FACES N'ONT PAS ETE TRAITEES, IL FAUT S'ARRETER
!              ------------------------------------------------------
    else

!           ==============================================

!             CONSERVER IMPERATIVEMENT LE TEST D'ARRET

!           ==============================================

      write (nfecra,2000) ifac,izfrdp(ifac),itypfb(ifac)
      iok = iok + 1

    endif

  enddo


  if (iok.ne.0) then
    write (nfecra,2100) iphas
    call csexit (1)
    !==========
  endif


endif

! -------
! FORMATS
! -------

 1000 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@                CONDITIONS AUX LIMITES NON RENSEIGNEES      ',/,&
'@                                                            ',/,&
'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )

 1100 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    LES CONDITIONS AUX LIMITES NE SONT PAS RENSEIGNEES POUR ',/,&
'@     CERTAINES FACES DE BORD (Phase ',I10   ,')             ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usray3.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT (FLUNET    NON RENSEIGNE)       ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )

 2100 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    LE FLUNET    N''EST PAS RENSEIGNEE POUR CERTAINES       ',/,&
'@        FACES DE BORD (Phase ',I10   ,')                    ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usray3.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)



end subroutine
