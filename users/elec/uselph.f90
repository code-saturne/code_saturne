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

subroutine uselph &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   maxelt , lstelt ,                                              &
   ibrom  , izfppp ,                                              &
   idevel , ituser , ia     ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!   REMPLISSAGE DES VARIABLES PHYSIQUES POUR LE MODULE ELECTRIQUE

!     ----> Effet Joule
!     ----> Arc Electrique
!     ----> Conduction Ionique

!      1) Masse Volumique
!      2) Viscosite moleculaire
!      3) Chaleur massique Cp
!      4) Lambda/Cp moleculaire
!      4) Diffusivite moleculaire



! ATTENTION :
! =========


! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)

! Pour le module electrique, toutes les proprietes physiques sont
!   supposees variables et contenues dans le tableau PROPCE
!   (meme si elles sont physiquement constantes)


! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usini1 :
!             . la masse volumique (initialisee a RO0(IPHAS))
!             . la viscosite       (initialisee a VISCL0(IPHAS))
!      - dans usiniv/useliv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a l
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


! Il est conseille de ne garder dans ce sous programme que
!    le strict necessaire.


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


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
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! nphmx            ! e  ! <-- ! nphsmx                                         !
! maxelt           ! i  ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!   (nphmx   )     !    !     !                                                !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
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
use cstphy
use entsor
use ppppar
use ppthch
use ppincl
use elincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse , nphmx

integer          maxelt, lstelt(maxelt)
integer          ibrom(nphmx)
integer          izfppp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          iel   , iphas
integer          ipcrom, ipcvis, ipccp , ipcvsl, ipcsig
integer          mode

double precision tp
double precision xkr   , xbr
double precision rom0  , temp0 , dilar , aa    , bb    , cc
double precision srrom1, rhonp1

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 0 - INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0

ipass = ipass + 1

iphas = 1

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================


!     En Joule, on s'arrete : il faut que l'utilisateur
!       donne les proprietes physiques
if ( ippmod(ieljou).ge.1 ) then

  write(nfecra,9010)
  call csexit (1)

!     En Arc on continue car on a un fichier de donnees
!       Un message indique que l'utilisateur n'a rien fourni
elseif(ippmod(ielarc).ge.1) then

  if(ipass.eq.1) then
    write(nfecra,9011)
  endif

  return

endif

 9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES PROP. PHYSIQUES   ',/,&
'@    =========                                               ',/,&
'@                      MODULE ELECTRIQUE                     ',/,&
'@                                                            ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uselph DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@     Ce sous-programme utilisateur permet de definir les    ',/,&
'@       proprietes physiques. Il est indispensable.          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9011 format(/,                                                   &
' Module arc electrique: pas d''intervention utilisateur pour ',/,&
'                          le calcul des proprietes physiques.',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!     Message au premier passage pour indiquer que l'utilisateur a
!       rapatrie le sous-programme.
if(ipass.eq.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 1 - EFFET JOULE
!===============================================================================

if ( ippmod(ieljou).ge.1 ) then


!     Attention, dans les modules electriques, la chaleur massique, la
!       conductivite thermique et la conductivite electriques sont
!       toujours dans le tableau PROPCE
!       qu'elles soient physiquement variables ou non.

!       On n'utilisera donc PAS les variables
!          =====================
!                                CP0(IPHAS), VISLS0(ISCALT(IPHAS))
!                                VISLS0(IPOTR) et VISLS0(IPOTI)

!       Informatiquement, ceci se traduit par le fait que
!                                ICP(IPHAS)>0, IVISLS(ISCALT(IPHAS))>0,
!                                IVISLS(IPOTR)>0 et IVISLS(IPOTI)>0





!       Calcul de la temperature a partir de l'enthalpie
!       ------------------------------------------------

!       Ceci depend largement des choix utilisateur en
!         matiere de loi H-T (T en Kelvin)

!       On demande de fournir cette loi dans le sous programme usthht
!          (USERS/base/usthht.F)
!           usthht fournit en particulier un exemple d'interpolation
!            a partir d'une tabulation utilisateur
!           usthht en mode T->H sera utilise pour l'initialisation
!            de l'enthalpie dans useliv.

!       MODE = 1 : H=RTP(IEL,ISCA(IHM)) -> T=PROPCE(IEL,IPPROC(ITEMP))
  mode = 1

  do iel = 1, ncel
    call usthht (mode,                                            &
         rtp(iel,isca(ihm)),propce(iel,ipproc(itemp)))
  enddo


!       Masse volumique au centre des cellules
!       --------------------------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la masse volumique ici
!       en renseignant PROPCE(IEL,IPCROM)
!       (meme si elle est uniforme ou constante).


!     Masse Vol : RO = ROM0 / (1+DILAR*(T-T0)
!         (Choudhary) semblable a Plard (HE-25/94/017)

!          avec sous-relaxation (sauf au premier pas de temps)

  temp0  = 300.d0
  rom0   = 2500.d0
  dilar  = 7.5d-5
  if(ntcabs.gt.1) then
    srrom1 = srrom
  else
    srrom1 = 0.d0
  endif

  ipcrom = ipproc(irom(iphas))
  do iel = 1, ncel
    rhonp1 = rom0 /                                               &
            (1.d0+ dilar * (propce(iel,ipproc(itemp))-temp0) )
    propce(iel,ipcrom) =                                          &
         srrom1*propce(iel,ipcrom)+(1.d0-srrom1)*rhonp1
  enddo


!       Viscosite moleculaire dynamique en kg/(m s)
!        ------------------------------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la viscosite ici
!       en renseignant PROPCE(IEL,IPCVIS)
!       (meme si elle est uniforme ou constante).


!     Viscosite : MU = EXP((AA/T-BB)-CC)
!          (Choudhary)
!      Plard (HE-25/94/017) ; limite a 1173K par C Delalondre

  ipcvis = ipproc(iviscl(iphas))
  aa     = 10425.d0
  bb     =   500.d0
  cc     =-6.0917d0

  do iel = 1, ncel
    if ( propce(iel,ipproc(itemp)) .gt. 1173.d0 ) then
      tp = propce(iel,ipproc(itemp))
    else
      tp= 1173.d0
    endif
    propce(iel,ipcvis) = exp( (aa/(tp-bb))+cc )
  enddo


!       Chaleur specifique J/(kg degres)
!       --------------------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la chaleur massique ici
!       en renseignant PROPCE(IEL,IPCPP)
!       (meme si elle est uniforme ou constante).


!        CP = 1381 (Choudhary)
!          coherent avec Plard (HE-25/94/017)

  ipccp  = ipproc(icp(iphas))
  do iel = 1, ncel
    propce(iel,ipccp) = 1381.d0
  enddo


!       Lambda/Cp en kg/(m s)
!       ---------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la conductivite ici
!       en renseignant PROPCE(IEL,IPCVSL)
!       (meme si elle est uniforme ou constante).


!         Lambda
!          On suppose Cp renseigne au prealable.

!          Plard (HE-25/94/017)

  ipcvsl = ipproc(ivisls(iscalt(iphas)))

  do iel = 1, ncel
    xbr = 85.25d0                                                 &
         -5.93d-2*(propce(iel,ipproc(itemp))-tkelvi)              &
         +2.39d-5*(propce(iel,ipproc(itemp))-tkelvi)**2
    xkr = 16.d0*stephn*(1.4d0)**2*(propce(iel,ipproc(itemp)))**3  &
         /(3.d0*xbr)

    propce(iel,ipcvsl) = 1.73d0 + xkr
  enddo

! --- On utilise CP calcule  dans PROPCE ci dessus
  do iel = 1, ncel
    propce(iel,ipcvsl) = propce(iel,ipcvsl)/propce(iel,ipccp)
  enddo


!       Conductivite electrique en S/m
!       ==============================

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la conductivite ici
!       en renseignant PROPCE(IEL,IPCSIG)
!       (meme si elle est uniforme ou constante).


!         SIGMA  (Plard HE-25/94/017)

  ipcsig = ipproc(ivisls(ipotr))
  do iel = 1, ncel
    propce(iel,ipcsig) =                                          &
         exp(7.605d0-7200.d0/propce(iel,ipproc(itemp)))
  enddo

!     La conductivite electrique pour le potentiel imaginaire est
!       toujours implicitement prise egale a la conductivite
!       utilisee pour le potentiel reel.
!       IL NE FAUT PAS la renseigner.

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
!     Conductivite electrique imaginaire :
!     La conductivite reelle et imaginaire sont dans le meme tableau.
!       Ce choix est fait en dur dans varpos.
!       Les pointeurs pour les deux existent quand meme.
!     Sinon, on pourrait faire ceci :
  if(1.eq.0) then
    if ( ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4 ) then
      do iel = 1, ncel
        propce(iel,ipproc(ivisls(ipoti))) =                       &
             propce(iel,ipproc(ivisls(ipotr)))
      enddo
    endif
  endif
!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!       Diffusivite variable a l'exclusion de l'enthalpie et du potentiel
!       -----------------------------------------------------------------
!     Pour le moment, il n'y a pas d'autres scalaires et
!                                                  on ne fait donc rien

endif

!===============================================================================
! 2 - ARC ELECTRIQUE
!===============================================================================

!     Les proprietes physiques sont a priori fournies par fichier
!       de donnees. IL n'y a donc rien a faire ici.

!      IF ( IPPMOD(IELARC).GE.1 ) THEN
!      ENDIF


!===============================================================================
! 3 - CONDUCTION IONIQUE
!===============================================================================

!     CETTE OPTION N'EST PAS ACTIVABLE

!--------
! FORMATS
!--------

 1000 format(/,                                                   &
' Module electrique: intervention utilisateur pour        ',/,    &
'                      le calcul des proprietes physiques.',/)

!----
! FIN
!----

return
end subroutine
