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

subroutine scalai &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   tslagr , coefa  , coefb  ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR LES SCALAIRES SUR UN PAS DE TEMPS

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
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
use entsor
use optcal
use cstphy
use cstnum
use pointe
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use elincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision tslagr(ncelet,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ra(*)

! Local variables

integer          idebia, idebra, ifinia
integer          iscal, ivar, iel
integer          ii, iisc, itspdv, icalc, iappel
integer          ispecf

double precision, allocatable, dimension(:) :: dtr
double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbrs, rovsdt

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Allocate temporary arrays for the species resolution
allocate(dtr(ncelet))
allocate(viscf(nfac), viscb(nfabor))
allocate(smbrs(ncelet), rovsdt(ncelet))

idebia = idbia0
idebra = idbra0

ipass = ipass + 1

!===============================================================================
! 2. TRAITEMENT DES SCALAIRES A PHYSIQUE PARTICULIERE
!    Cette section sera dediee aux traitement particuliers.
!    On donne un exemple qui n'est que la recopie (a ISCAPP pres) de
!      la section 3 sur les scalaires utilisateurs.
!===============================================================================

if (ippmod(iphpar).ge.1) then

! ---> Initialisation des RTP  a partir des CL

!   On initialise RTP (et pas RTPA) car RTPA ne sert pas
!      dans codits.

  call ppinv2                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  , &
   ra     )

!     On pourra eviter les bugs en initialisant aussi RTPA (sur NCELET)
!     (d'ailleurs, on s'en sert ensuite dans le cpflux etc)
  if(ipass.eq.1.and.isuite.eq.0) then
    if(nscapp.gt.0) then
      do ii = 1, nscapp
        iscal = iscapp(ii)
        ivar  = isca(iscal)
        do iel = 1, ncelet
          rtpa (iel,ivar) = rtp (iel,ivar)
        enddo
      enddo
    endif
  endif

! ---> Calculs TS relatifs a la physique du charbon
!      GMDEV1, GMDEV2, GMHET, GMDCH

  if ( ippmod(icp3pl).ne.-1 ) then

    call cpflux                                                   &
    !==========
   ( idebia , idebra , ncelet , ncel   ,                          &
     rtpa   , propce , volume ,                                   &
     ra     )

  endif

! ---> Calculs TS relatifs a la physique du fuel
!      GAMEVA, GAMHTF

  if ( ippmod(icfuel).ne.-1 ) then

    call fuflux                                                   &
    !==========
   ( idebia , idebra , ncelet , ncel   ,                          &
     rtpa   , propce , volume ,                                   &
     ra     )

  endif

!    ATTENTION : POUR LE CLIPPING AVEC ICLP = 1, IL FAUT

!                QUE LES SCALAIRES SOIENT RESOLUS AVANT
!                LEUR VARIANCE ASSOCIEE






! ---> Boucle sur les scalaires physiques particulieres.
!      On peut imaginer a la place des resolutions couplees.
!      Ici, on ne donne qu'un exemple.

  do ii = 1, nscapp

    iscal = iscapp(ii)
    ivar  = isca(iscal)

! ---> Pas de temps (avec facteur multiplicatif eventuel)

    if(cdtvar(ivar).ne.1.d0) then
      do iel = 1, ncel
        dtr(iel) = dt(iel)*cdtvar(ivar)
      enddo
    else
      do iel = 1, ncel
        dtr(iel) = dt(iel)
      enddo
    endif

!     Schema compressible sans choc :
! ---> Traitement special pour la masse volumique,
!                     la temperature et l'energie
!     L'indicateur ISPECF sera non nul si l'on ne doit pas resoudre
!       le scalaire plus bas avec covofi.

    ispecf = 0

    if ( ippmod(icompf).ge.0 ) then

      if ( iscal.eq.irho .or. iscal.eq.itempk ) then
        ispecf = 1
      elseif ( iscal.eq.ienerg ) then
        ispecf = 2
      endif

! ---> Masse volumique : deja resolue
! ---> Temperature     : n'est pas une variable transportee
! ---> Enthalpie       : a resoudre

      if(ispecf.eq.2) then

        call cfener &
        !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ncepdc , ncetsm ,                                              &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt ,                                              &
   ra     )

      endif

    endif

!     Pour le compressible, on ne resout pas celles deja resolues ou
!       qui ne doivent pas l'etre
!     Pour les autres physiques, on resout les scalaires dans l'ordre
!       (pas de tests a faire)

    if(ispecf.eq.0) then

! ---> Variances et scalaires

!   iscavr dira scalaire/variance
!   itspdv dira si calcul des termes de prod et dissip supplementaires
!         ou non (1 : oui, 0 : non)

!         si iscavr = 0
!           scalaire
!           itspdv = 0
!         sinon
!           variance
!           si iscavr > 0 et iscavr < nscal+1
!             itspdv = 1
!           sinon
!             pour le moment, on s'arrete
!             a terme, les combustionnistes pourront donner leur propre
!               grandeur scalaire associee a la variance et
!               eventuellement reconstruite en dehors de covofi.
!               il faudra alors peut etre declarer des tableaux
!               de travail suppl
!           fin si
!         fin si

      iisc = iscal
      if(iscavr(iisc).eq.0) then
        itspdv = 0
      elseif(iscavr(iisc).gt.0.and.iscavr(iisc).le.nscal) then
        itspdv = 1
      else
        write(nfecra,9000)iisc,iisc,nscal,iscavr(iisc)
        call csexit (1)
      endif


! ---> Appel a covofi pour la resolution

      call covofi                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ncepdc , ncetsm ,                                              &
   iisc   , itspdv ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dtr    , rtp    , rtpa   , propce , propfa , propfb , tslagr , &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt ,                                              &
   ra     )


! ---> Versions Electriques
!             Effet Joule
!             Arc Electrique
!             Conduction ionique

!     On calcule ici j, E, j.E reels et imagimaires

      if ( ippmod(ieljou).ge.1 .or.                               &
           ippmod(ielarc).ge.1 .or.                               &
           ippmod(ielion).ge.1       ) then


!     On utilise le  fait que les scalaires sont dans l'ordre
!       H, PotR, [PotI], [A] pour faire le calcul de j, E et j.E
!       apres  la determination de PotR [et PotI].

        icalc = 0
!            On peut y aller apres PotR si on est en arc
!                                         ou en Joule sans PotI

        if(ippmod(ielarc).ge.1.or.ippmod(ieljou).eq.1             &
             .or.ippmod(ieljou).eq.3) then
          if(iscal.eq.ipotr) then
            icalc = 1
          endif
        endif
!     On y va apres PotI si on est en Joule avec PotI
        if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
          if(iscal.eq.ipoti) then
            icalc = 1
          endif
        endif

        if(icalc.eq.1) then

!     Calcul de j, E et j.E
          iappel = 1

          call elflux                                             &
          !==========
 ( idebia , idebra , iappel ,                                     &
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , viscf  , viscb  ,                            &
   ra     )


!     Recalage des variables electriques j, j.E (et Pot, E)

          if ( ielcor .eq.1  .and. ntcabs .gt. 1 ) then

            call uselrc                                           &
            !==========
  ( nvar   , nscal  ,                                             &
    ia     ,                                                      &
    dt     , rtpa   , rtp    , propce , propfa , propfb ,         &
    coefa  , coefb  ,                                             &
    ra     )

          endif

        endif

      endif

    endif


! ---> Fin de la Boucle sur les scalaires physiques particulieres.
  enddo
endif


!     On calcule ici A, B, jxB

if ( ippmod(ielarc).ge.1       ) then

!     On utilise le  fait que les scalaires sont dans l'ordre
!       H, PotR, [PotI], [A] pour faire le calcul de A, B, jxB
!       apres la determination et le recalage de j
  iappel = 2

  call elflux                                                     &
  !==========
 ( idebia , idebra , iappel ,                                     &
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , viscf  , viscb  ,                            &
   ra     )

endif

!===============================================================================
! 3. TRAITEMENT DES SCALAIRES UTILISATEURS STANDARD
!     On voit que meme s'ils sont numerotes en premier, on peut les
!       traiter en dernier si on veut.
!     On peut imaginer aussi de by-passer cette phase si le modele a
!       physique particuliere le demande.
!===============================================================================

if(nscaus.gt.0) then

! ---> Boucle sur les scalaires utilisateur.

  do ii = 1, nscaus

    iscal = ii
    ivar  = isca(iscal)

! ---> Pas de temps (avec facteur multiplicatif eventuel)

    if(cdtvar(ivar).ne.1.d0) then
      do iel = 1, ncel
        dtr(iel) = dt(iel)*cdtvar(ivar)
      enddo
    else
      do iel = 1, ncel
        dtr(iel) = dt(iel)
      enddo
    endif


! ---> Variances et scalaires

!   iscavr dira scalaire/variance
!   itspdv dira si calcul des termes de prod et dissip supplementaires
!         ou non (1 : oui, 0 : non)

!         si iscavr = 0
!           scalaire
!           itspdv = 0
!         sinon
!           variance
!           si iscavr > 0 et iscavr < nscal+1
!             itspdv = 1
!           sinon
!             pour le moment, on s'arrete
!             a terme, les combustionnistes pourront donner leur propre
!               grandeur scalaire associee a la variance et
!               eventuellement reconstruite en dehors de covofi.
!               il faudra alors peut etre declarer des tableaux
!               de travail suppl
!           fin si
!         fin si

    iisc = iscal
    if(iscavr(iisc).eq.0) then
      itspdv = 0
    elseif(iscavr(iisc).gt.0.and.iscavr(iisc).le.nscal) then
      itspdv = 1
    else
      write(nfecra,9000)iisc,iisc,nscal,iscavr(iisc)
      call csexit (1)
    endif


! ---> Appel a covofi pour la resolution

    call covofi                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ncepdc , ncetsm ,                                              &
   iisc   , itspdv ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dtr    , rtp    , rtpa   , propce , propfa , propfb , tslagr , &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt ,                                              &
   ra     )






! ---> Fin de la Boucle sur les scalaires utilisateurs.
  enddo

endif

! Free memory
deallocate(dtr)
deallocate(viscf, viscb)
deallocate(smbrs, rovsdt)

!===============================================================================
! 4.  FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA RESOLUTION DES SCALAIRES         ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE NUMERO ',I10                                    ,/,&
'@    ISCAVR(',I10   ,') DOIT ETRE UN ENTIER                  ',/,&
'@      POSITIF OU NUL ET                                     ',/,&
'@      INFERIEUR OU EGAL A NSCAL = ',I10                      ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Si ISCAVR(I) est nul, le scalaire I n est pas une variance',/,&
'@  Si ISCAVR(I) est positif, le scalaire I est une variance :',/,&
'@    il s agit de la variance des fluctuations du scalaire J ',/,&
'@    dont le numero est ISCAVR(I)                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT WHILE SOLVING SCALARS EQUATIONS          ',/,&
'@    ========                                                ',/,&
'@    SCALAR NUMBER ',I10                                      ,/,&
'@    ISCAVR(',I10   ,') MUST BE A POSITIVE OR NULL INTEGER   ',/,&
'@      AND LOWER OR EQUAL THAN NSCAL = ', I10                 ,/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculaton will not be run.                           ',/,&
'@                                                            ',/,&
'@  If ISCAVR(I) is null, the scalar I is not a variance      ',/,&
'@  If ISCAVR(I) is positive, the scalar I is a variance:     ',/,&
'@    it is the variance of the fluctuations of the scalar J  ',/,&
'@    whose number is ISCAVR(I)                               ',/,&
'@                                                            ',/,&
'@  Verify usini1.                                            ',/,&
'@  Contact support.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! FIN
!----
return

end subroutine
