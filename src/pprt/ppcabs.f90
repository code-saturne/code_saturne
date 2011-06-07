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

subroutine ppcabs &
!================

 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     )

!===============================================================================
! FONCTION :
! ----------


!   SOUS-PROGRAMME PHYSIQUES PARTICULIERES

!  DONNE LA VALEUR DU COEFFICIENT D'ABSORPTION POUR
!    LE MELANGE GAZEUX ET LES PARTICULES POUR LE CP.

!  DANS LE CAS DU MODELE P-1 ON VERIFIE QUE LA LONGUEUR OPTIQUE
!    DU MILIEU EST AU MINIMUM DE L'ORDRE DE L'UNITE


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...3(ncelet    ! tr ! --- ! tableau de travail                             !
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
use parall
use ppppar
use ppthch
use coincl
use cpincl
use fuincl
use ppincl
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          itypfb(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet), w2(ncelet), w3(ncelet)



! Local variables

integer          iel, ifac, icla, ipck, icha, iok
double precision xm, d2, vv, sf, xlc, xkmin, pp

!===============================================================================

!===============================================================================
! 0 - GESTION MEMOIRE
!===============================================================================


!===============================================================================
!  1 - COEFFICIENT D'ABSORPTION DU MELANGE GAZEUX (m-1)
!===============================================================================

if ( ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0 ) then

! ----> Combustion gaz : Flamme de diffusion
!                        Flamme de premelange (Modele EBU)

  if (imodak.eq.1) then

    do iel = 1, ncel
      xm = 1.d0/ (  propce(iel,ipproc(iym(1)))/wmolg(1)                 &
                  + propce(iel,ipproc(iym(2)))/wmolg(2)                 &
                  + propce(iel,ipproc(iym(3)))/wmolg(3) )
      w1(iel) = propce(iel,ipproc(iym(3)))*xm/wmolg(3)*xco2
      w2(iel) = propce(iel,ipproc(iym(3)))*xm/wmolg(3)*xh2o
      w3(iel) = 0.d0
    enddo
    call raydak(ncel,ncelet,                                      &
    !==========
      propce(1,ipproc(icak(1))),w1,w2,w3,propce(1,ipproc(itemp)))

    write(NFECRA,*) ' a verifier '
    write(NFECRA,*) ' a finir   : raydak '
    write(NFECRA,*) ' Le codage est a terminer par le groupe I81'
    write(NFECRA,*) '                     13-10-03 22:38:03      '
    call csexit(1)

  else
    do iel = 1, ncel
      propce(iel,ipproc(icak(1))) = propce(iel,ipproc(ickabs))
    enddo
  endif

else if ( ippmod(icp3pl).ge.0 ) then

! ---->  Charbon

  if (imodak.eq.1) then

    do iel = 1,ncel
! concentration volumique en CO2
      w1(iel) = propce(iel,ipproc(immel))/wmole(ico2)             &
               *propce(iel,ipproc(iym1(ico2)))
! concentration volumique en H20
      w2(iel) = propce(iel,ipproc(immel))/wmole(ih2o)             &
               *propce(iel,ipproc(iym1(ih2o)))
! fraction volumique de suies
      w3(iel) = 0.d0

    enddo

    call raydak(ncel,ncelet,                                      &
    !==========
     propce(1,ipproc(icak(1))),w1,w2,w3,propce(1,ipproc(itemp1)))

  else
    do iel = 1, ncel
      propce(iel,ipproc(icak(1))) = ckabs1
    enddo
  endif

else if ( ippmod(icfuel).ge.0 ) then

! ---->  Fuel

  if (imodak.eq.1) then

    do iel = 1,ncel
! concentration volumique en CO2
      w1(iel) = propce(iel,ipproc(immel))/wmole(ico2)             &
               *propce(iel,ipproc(iym1(ico2)))
! concentration volumique en H20
      w2(iel) = propce(iel,ipproc(immel))/wmole(ih2o)             &
               *propce(iel,ipproc(iym1(ih2o)))
! fraction volumique de suies
      w3(iel) = 0.d0

    enddo

    call raydak(ncel,ncelet,                                      &
    !==========
     propce(1,ipproc(icak(1))),w1,w2,w3,propce(1,ipproc(itemp1)))

  else
    do iel = 1, ncel
      propce(iel,ipproc(icak(1))) = ckabs1
    enddo
  endif

endif


!===============================================================================
!  2 - COEFFICIENT D'ABSORPTION DES PARTICULES PAR CLASSE K2/X2 (m-1)
!===============================================================================

! ---->  Charbon

if ( ippmod(icp3pl).ge.0 ) then

  do icla = 1, nclacp

    ipck = 1 + icla
    icha = ichcor(icla)

    do iel = 1, ncel

! ---> Calcul du diametre des particules

      d2 = ( xashch(icha)*diam20(icla)**2 +                       &
           ( 1.d0-xashch(icha))                                   &
             *propce(iel,ipproc(idiam2(icla)))**2 )**0.5d0

! ---> Calcul du coeficient d'absorption des particules K2/X2
!         3./2. ROM/(ROM2*D2)

      propce(iel,ipproc(icak(ipck))) =                            &
                         1.5d0*propce(iel,ipproc(irom))    &
                       / ( propce(iel,ipproc(irom2(icla)))*d2)

    enddo

  enddo

endif

! ---->  Fuel

if ( ippmod(icfuel).ge.0 ) then

  do icla = 1, nclafu

    ipck = 1 + icla

    do iel = 1, ncel

! ---> Calcul du coeficient d'absorption des particules K2/X2
!         3./2. ROM/(ROM2*D2)

      propce(iel,ipproc(icak(ipck))) =                            &
                         1.5d0*propce(iel,ipproc(irom))    &
                 / ( propce(iel,ipproc(irom3(icla)))              &
                    *propce(iel,ipproc(idiam3(icla))) )

    enddo

  enddo

endif

!===============================================================================
!  3 - COEFFICIENT D'ABSORPTION GAZ (Arc Electrique)
!===============================================================================


if ( ippmod(ielarc).ge.1 ) then

  do iel = 1, ncel

! ---> Directement donne par le fichier dp_elec

      propce(iel,ipproc(icak(1))) = propce(iel,ipproc(idrad))

  enddo

endif

!===============================================================================
!  4 - CLIPPING DU COEFFICIENT D'ABSORPTION DANS LA CAS DE L'APPROX P-1
!===============================================================================


!--> MODELE P-1 : Controle standard des valeurs du coefficient
!                 d'absorption. Ce coefficient doit assurer une
!                 longueur optique au minimum de l'ordre de l'unite.


  if (iirayo.eq.2) then

!         Coefficient d'absorption du melange gaz-particules de charbon

     do iel = 1, ncel
       w3(iel) =  propce(iel,ipproc(icak(1)))
     enddo

     if ( ippmod(icp3pl).ge.0 ) then
       do icla = 1,nclacp
         ipck = 1+icla
         do iel = 1,ncel
           w3(iel) = w3(iel)                                      &
                    + ( propce(iel,ipproc(ix2(icla)))             &
                      * propce(iel,ipproc(icak(ipck))) )
         enddo
       enddo
     elseif ( ippmod(icfuel).ge.0 ) then
       do icla = 1,nclafu
         ipck = 1+icla
         do iel = 1,ncel
           w3(iel) = w3(iel)                                      &
                    + ( propce(iel,ipproc(iyfol(icla)))           &
                      * propce(iel,ipproc(icak(ipck))) )
         enddo
       enddo
     endif

!         Calcul de la longueur caractéristique XLC du domaine de calcul

    sf = 0.d0
    vv = 0.d0

    do ifac = 1,nfabor
       sf = sf + sqrt(surfbo(1,ifac)**2 + surfbo(2,ifac)**2 + surfbo(3,ifac)**2)
    enddo
    if (irangp.ge.0) then
      call parsom(sf)
      !==========
    endif

    do iel = 1,ncel
       vv = vv + volume(iel)
    enddo
    if (irangp.ge.0) then
      call parsom(vv)
      !==========
    endif

    xlc = 3.6d0 * vv / sf

!         Clipping de CK a XKMIN

    xkmin = 1.d0 / xlc

    iok = 0
    do iel = 1,ncel
      if (w3(iel).lt.xkmin) then
        iok = iok +1
      endif
    enddo
    if (irangp.ge.0) then
      call parcpt(iok)
    endif

!     Arret en fin de pas de temps si epaisseur optique trop grande
    pp = xnp1mx/100.0d0
    if (dble(iok).gt.pp*dble(ncelgb)) then
       write(nfecra,1000) xkmin, dble(iok)/dble(ncelgb)*100.d0, xnp1mx
       istpp1 = 1
    endif

  endif

! -------
! FORMAT
! -------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT APPROXIMATION P-1 (PPCABS)      ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    LA LONGUEUR OPTIQUE DU MILIEU SEMI-TRANSPARENT          ',/,&
'@      DOIT AU MOINS ETRE DE L''ORDRE DE L''UNITE POUR ETRE  ',/,&
'@      DANS LE DOMAINE D''APPLICATION DE L''APPROXIMATION P-1',/,&
'@    CELA NE SEMBLE PAS ETRE LE CAS ICI.                     ',/,&
'@                                                            ',/,&
'@    LE COEFFICIENT D''ABSORPTION MINIMUM POUR ASSURER CETTE ',/,&
'@      LONGUEUR OPTIQUE EST XKMIN = ',E10.4                   ,/,&
'@    CETTE VALEUR N''EST PAS ATTEINTE POUR ', E10.4,'%       ',/,&
'@      DES CELLULES DU MAILLAGE.                             ',/,&
'@    LE POURCENTAGE DE CELLULES DU MAILLAGE POUR LESQUELLES  ',/,&
'@      ON ADMET QUE CETTE CONDITION SOIT VIOLEE EST IMPOSE   ',/,&
'@      PAR DEFAUT OU DANS USINI1 A XNP1MX = ', E10.4,'%      ',/,&
'@                                                            ',/,&
'@    Le calcul est interrompu.                               ',/,&
'@                                                            ',/,&
'@    Verifier les valeurs du coefficient d''absorption CK    ',/,&
'@      dans PPCABS, USRAY3 ou Fichier thermochimie.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return

end subroutine
