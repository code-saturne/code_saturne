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

subroutine cfphyv &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : COMPRESSIBLE SANS CHOC

! Calcul des proprietes physiques variables


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!        !    !     !                                                !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...3(ncelet    ! tr ! --- ! tableau de travail                             !
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
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          ibrom
integer          izfppp(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra, ifinia
integer          iel
integer          ifac
integer          iirom , iiromb
integer          maxelt, ils

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0


!===============================================================================
! 2. ON DONNE LA MAIN A L'UTILISATEUR
!===============================================================================

maxelt = max(ncelet, nfac, nfabor)
ils    = idebia
ifinia = ils + maxelt
call iasize('cfphyv',ifinia)

iuscfp = 1
call uscfpv                                                       &
!==========
 ( ifinia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   maxelt , ia(ils),                                              &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     ,                                     &
   ra     )

!     Si IUSCFP = 0, l'utilisateur n'a pas inclus le ss pgm uscfpv dans
!       ses sources. C'est une erreur si Cp, Cv ou Lambda est variable.
!     On se contente de faire le test au premier passage.
if(ipass.eq.0) then
  ipass = ipass + 1
  if((ivisls(itempk).gt.0.or.                            &
       icp.gt.0.or.icv.gt.0).and.iuscfp.eq.0) then
    write(nfecra,1000)                                          &
         ivisls(itempk),icp,icv
    call csexit (1)
    !==========
  endif
endif

!===============================================================================
! 3. MISE A JOUR DE LAMBDA/CV
!===============================================================================

! On a vérifié auparavant que CV0 était non nul.
! Si CV variable est nul, c'est une erreur utilisateur. On fait
!     un test à tous les passages (pas optimal), sachant que pour
!     le moment, on est en gaz parfait avec CV constant : si quelqu'un
!     essaye du CV variable, ce serait dommage que cela lui explose à la
!     figure pour de mauvaises raisons.
! Si IVISLS(IENERG).EQ.0, on a forcement IVISLS(ITEMPK).EQ.0
!     et ICV.EQ.0, par construction de IVISLS(IENERG) dans
!     le sous-programme cfvarp

if(ivisls(ienerg).gt.0) then

  if(ivisls(itempk).gt.0) then

    do iel = 1, ncel
      propce(iel,ipproc(ivisls(ienerg))) =               &
           propce(iel,ipproc(ivisls(itempk)))
    enddo

  else
    do iel = 1, ncel
      propce(iel,ipproc(ivisls(ienerg))) =               &
           visls0(itempk)
    enddo

  endif

  if(icv.gt.0) then

    do iel = 1, ncel
      if(propce(iel,ipproc(icv)).le.0.d0) then
        write(nfecra,2000)iel,propce(iel,ipproc(icv))
        call csexit (1)
        !==========
      endif
    enddo

    do iel = 1, ncel
      propce(iel,ipproc(ivisls(ienerg))) =               &
           propce(iel,ipproc(ivisls(ienerg)))            &
           / propce(iel,ipproc(icv))
    enddo

  else

    do iel = 1, ncel
      propce(iel,ipproc(ivisls(ienerg))) =               &
           propce(iel,ipproc(ivisls(ienerg)))            &
           / cv0
    enddo

  endif

else

  visls0(ienerg) = visls0(itempk)/cv0

endif

!===============================================================================
! 3. MISE A JOUR DE ROM et ROMB :
!     On ne s'en sert a priori pas, mais la variable existe
!     On a ici des valeurs issues du pas de temps précédent (y compris
!       pour les conditions aux limites) ou issues de valeurs initiales
!     L'échange pério/parall sera fait dans phyvar.
!===============================================================================

iirom  = ipproc(irom  )
iiromb = ipprob(irom  )

do iel = 1, ncel
  propce(iel,iirom)  = rtpa(iel,isca(irho))
enddo

do ifac = 1, nfabor
  iel = ifabor(ifac)
  propfb(ifac,iiromb) =                                         &
       coefa(ifac,iclrtp(isca(irho),icoef))              &
       + coefb(ifac,iclrtp(isca(irho),icoef))            &
       * rtpa(iel,isca(irho))
enddo

!--------
! FORMATS
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION (MODULE COMPRESSIBLE)  ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  Une ou plusieurs des propriétés suivantes a été déclarée  ',/,&
'@    variable (repérée ci-dessous par un indicateur non nul) ',/,&
'@    et une loi doit être fournie dans uscfpv.               ',/,&
'@         propriété                               indicateur ',/,&
'@     - conductivité thermique                    ',I10       ,/,&
'@     - capacité calorifique à pression constante ',I10       ,/,&
'@     - capacité calorifique à volume constant    ',I10       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Renseigner uscfpv ou déclarer les propriétés constantes et',/,&
'@    uniformes (uscfx2 pour la conductivité thermique,       ',/,&
'@    uscfth pour les capacités calorifiques).                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION (MODULE COMPRESSIBLE)  ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  La capacité calorifique à volume constant présente (au    ',/,&
'@    moins) une valeur négative ou nulle :                   ',/,&
'@    cellule ',I10,   '  Cv = ',E18.9                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfpv.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!----
! FIN
!----

return
end subroutine
