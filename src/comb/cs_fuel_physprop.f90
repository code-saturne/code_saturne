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

subroutine cs_fuel_physprop &
!==========================

 ( mbrom  , izfppp ,                                              &
   rtp    , rtpa   , propce )

!===============================================================================
! FONCTION :
! --------
! ROUTINE PHYSIQUE PARTICULIERE : COMBUSTION CHARBON PULVERISE
! Calcul de RHO du melange
!
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! mbrom            ! te ! <-- ! indicateur de remplissage de romb              !
!   (nphmx   )     !    !     !                                                !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field
!===============================================================================

implicit none

! Arguments

integer          mbrom
integer          izfppp(nfabor)

double precision rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)

! Local variables

integer          iel, icla, ipcro2
integer          izone, ifac
integer          iromf , ioxy , nbclip1,nbclip2

double precision x1sro1, x2sro2, srrom1, uns1pw
double precision x2tot, wmolme, unsro1
double precision ff3min,ff3max,valmin,valmax

integer          ipass
data             ipass /0/
save             ipass

integer          iok1,iok2,iok3
double precision , dimension ( : )     , allocatable :: f1m,f2m,f3m,f4m,f5m
double precision , dimension ( : )     , allocatable :: f6m,f7m,f8m,f9m
double precision , dimension ( : )     , allocatable :: enth1 , x2 ,fvp2m
double precision , dimension ( : )     , allocatable :: xoxyd,enthox
double precision, dimension(:), pointer ::  brom, crom
!===============================================================================
!
!===============================================================================
! 0. ON COMPTE LES PASSAGES
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

!===============================================================================
! Deallocation dynamic arrays
!----
allocate(f1m(1:ncelet),f2m(1:ncelet),f3m(1:ncelet),              STAT=iok1)
allocate(f4m(1:ncelet),f5m(1:ncelet),                            STAT=iok1)
allocate(f6m(1:ncelet),f7m(1:ncelet),f8m(1:ncelet),f9m(1:ncelet),STAT=iok2)
allocate(enth1(1:ncel),x2(1:ncel)   ,fvp2m(1:ncel),              STAT=iok3)
!----
if ( iok1 > 0 .or. iok2 > 0 .or. iok3 > 0) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_fuel_physprop            '
  call csexit(1)
endif
if ( ieqnox .eq. 1 ) then
!----
  allocate(xoxyd(1:ncelet),enthox(1:ncelet),STAT=iok1)
!----
  if ( iok1 > 0 ) then
    write(nfecra,*) ' Memory allocation error inside:         '
    write(nfecra,*) '   cs_fuel_physprop for xoxyd and enthox '
  endif
endif
!===============================================================================

!     Pointeur sur masse volumique du gaz aux cellules
iromf = ipproc(irom1)
!
!===============================================================================
! 2. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE DISPERSEE
!                    VALEURS CELLULES
!                    ----------------
!    FRACTION MASSIQUE DE SOLIDE
!    DIAMETRE
!    MASSE VOLUMIQUE
!===============================================================================
!
call cs_fuel_physprop2 ( ncelet , ncel , rtp , propce )
!=====================

!===============================================================================
! 3. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE GAZEUSE
!                    VALEURS CELLULES
!                    ----------------
!    TEMPERATURE
!    MASSE VOLUMIQUE
!    CONCENTRATIONS DES ESPECES GAZEUSES
!===============================================================================

! --- Calcul de l'enthalpie du gaz     enth1
!        de F1M
!        de F2M
!        de F3M                    dans W3=1-F1M-F2M-F4M-F5M-F6M-F7M-F8M-F9M
!        de F4M
!        de F5M
!        de F6M
!        de F7M
!        de F8M
!        de F9M
!        de FVP2M
!
! Initialisation des fm et de x2 a 0
! f1m est toujours egal zero. Dans lecontexte de la combustion fioul il y a qu'un
! seul combustible gazeuse.
f1m( : ) = 0.d0
! f2m = ifvap
f2m( : ) = 0.d0
! f3m, f4m, f5m = Oxydants
f3m( : ) = 0.d0
f4m( : ) = 0.d0
f5m( : ) = 0.d0
! Vapeur d'eau
f6m( : ) = 0.d0
! Combustion heterogene
f7m( : ) = 0.d0
! f8m, f9m est toujours egal zero.
f8m( : ) = 0.d0
f9m( : ) = 0.d0
! fraction massique de la phase liquide.
x2 ( : ) = 0.d0
!
! - nclacp = nclafu
! - X2 est calcule directement en fonction de iyfol
do icla = 1, nclafu
  do iel = 1, ncel
    x2(iel) =  x2(iel) + rtp(iel,isca(iyfol(icla)))
  enddo
enddo

do iel = 1, ncel
  f2m(iel) =  f2m(iel) + rtp(iel,isca(ifvap))
enddo

if ( ieqnox .eq. 1 ) then
  do iel = 1, ncel
    xoxyd(iel)= (1.d0-x2(iel))-f1m(iel)-f2m(iel)
  enddo
endif

ff3min = 1.d+20
ff3max =-1.d+20
nbclip1= 0
nbclip2= 0
valmin = 1.d+20
valmax =-1.d+20
!
do iel = 1, ncel
  uns1pw = 1.d0/(1.d0-x2(iel))

! -Pour l'instant, la variable <noxyd> n'existe pas dans la version fioul car
!  elle n'est pas declarer dans fulecd.f90. Il faut qu'on l'ajoute dans ce
!  subroutine (mot clef: oxycombustion).
  if ( noxyd .ge. 2 ) then
    f4m(iel) = rtp(iel,isca(if4m))
    if ( noxyd .eq. 3 ) then
      f5m(iel) = rtp(iel,isca(if5m))
    endif
  endif

! - Scalaire relatif a la combustion heterogene.
  f7m(iel) =  rtp(iel,isca(if7m))

! La variance transportee.
  fvp2m(iel) = rtp(iel,isca(ifvp2m))

! Unites: kg scalaires / kg gas
  f1m(iel)  = f1m(iel)    *uns1pw
  f2m(iel)  = f2m(iel)    *uns1pw
  f4m(iel)  = f4m(iel)    *uns1pw
  f5m(iel)  = f5m(iel)    *uns1pw
  f6m(iel)  = f6m(iel)    *uns1pw
  f7m(iel)  = f7m(iel)    *uns1pw
  f8m(iel)  = f8m(iel)    *uns1pw
  f9m(iel)  = f9m(iel)    *uns1pw

  fvp2m(iel)= fvp2m(iel)*uns1pw

  f3m(iel) = 1.d0                                        &
           -( f1m(iel)+f2m(iel)+f4m(iel)+f5m(iel)        &
             +f6m(iel)+f7m(iel)+f8m(iel)+f9m(iel) )

  ff3max = max(ff3max,f3m(iel))
  ff3min = min(ff3min,f3m(iel))
!
  if ( ieqnox .eq. 1 ) then
    enthox(iel) = rtp(iel,isca(ihox))/xoxyd(iel)
  endif
!
enddo
!
if ( irangp .ge. 0 ) then
  call parmin(ff3min)
  call parmax(ff3max)
  call parcpt(nbclip1)
  call parcpt(nbclip2)
  call parmin(valmin)
  call parmax(valmax)
endif
write(nfecra,*) ' Values of f3 min and max: ',FF3MIN,FF3MAX
if ( nbclip1 .gt. 0 ) then
  write(nfecra,*) ' Clipping phase gas variance in min:',nbclip1,valmin
endif
if ( nbclip2 .gt. 0 ) then
  write(nfecra,*) ' Clipping phase gas variance in max:',nbclip1,valmin
endif

! ---- Enthalpie du gaz H1
enth1( : ) =0.D0
do icla = 1, nclafu
  do iel = 1, ncel
    enth1(iel) =  enth1(iel) + rtp(iel,isca(ih2(icla)))
  enddo
enddo
do iel = 1, ncel
  enth1(iel) = (rtp(iel,isca(ihm))-enth1(iel))/ ( 1.d0-x2(iel) )
enddo

call cs_fuel_physprop1                 &
!=================================
 ( ncelet , ncel   ,                                      &
   f1m    , f2m    , f3m    , f4m    , f5m    ,           &
   f6m    , f7m    , f8m    , f9m    , fvp2m  ,           &
   enth1  , enthox ,                                      &
   rtp    , propce , propce(1,iromf)   )

!===============================================================================
! 4. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE DISPERSEE
!                    VALEURS CELLULES
!                    ----------------
!    TEMPERATURE
!===============================================================================

! --- Transport d'H2

call  cs_fuel_thfieldconv2 ( ncelet , ncel , rtp , propce )
!=========================

!===============================================================================
! 5. CALCUL DES PROPRIETES PHYSIQUES DU MELANGE
!                    VALEURS CELLULES
!                    ----------------
!    MASSE VOLUMIQUE
!===============================================================================
! --- Calcul de Rho du melange : 1/Rho = X1/Rho1 + Somme(X2/Rho2)
!     On sous relaxe quand on a un rho n a disposition, ie
!       a partir du deuxieme passage ou
!       a partir du premier passage si on est en suite de calcul et
!         qu'on a relu la masse volumique dans le fichier suite.

call field_get_val_s(icrom, crom)

if (ipass.gt.1.or.(isuite.eq.1.and.initro.eq.1)) then
  srrom1 = srrom
else
  srrom1 = 0.d0
endif

do iel = 1, ncel
  x2sro2 = zero

  do icla = 1, nclafu
    ipcro2 = ipproc(irom2(icla))
    x2sro2 = x2sro2 + rtp(iel,isca(iyfol(icla))) / propce(iel,ipcro2)
  enddo
  x1sro1 = (1.d0-x2(iel)) / propce(iel,iromf)
! ---- Sous relaxation eventuelle a donner dans ppini1.F
  crom(iel) = srrom1*crom(iel)                  &
                     + (1.d0-srrom1)/(x1sro1+x2sro2)
enddo


!===============================================================================
! 6. CALCUL DE RHO DU MELANGE
!                      VALEURS FACETTES
!                      ----------------
!===============================================================================

mbrom = 1
call field_get_val_s(ibrom, brom)
call field_get_val_s(icrom, crom)
! ---> Masse volumique au bord pour toutes les facettes
!      Les facettes d'entree seront recalculees.
!
do ifac = 1, nfabor
  iel = ifabor(ifac)
  brom(ifac) = crom(iel)
enddo

! ---> Masse volumique au bord pour les facettes d'entree UNIQUEMENT
!     Le test sur IZONE sert pour les reprises de calcul

if ( ipass.gt.1 .or. isuite.eq.1 ) then
  do ifac = 1, nfabor

    izone = izfppp(ifac)
    if(izone.gt.0) then
      if ( ientat(izone).eq.1 .or. ientfl(izone).eq.1 ) then
        x2sro2 = zero
        x2tot  = zero
        do icla = 1, nclafu
          x2sro2 = x2sro2 + x20(izone,icla)/rho0fl
          x2tot  = x2tot  + x20(izone,icla)
        enddo

        ioxy = inmoxy(izone)
        wmolme =( oxyo2(ioxy)+oxyn2(ioxy)                         &
                 +oxyh2o(ioxy)+oxyco2(ioxy))                      &
               /( wmole(io2) *oxyo2(ioxy)                         &
                 +wmole(in2) *oxyn2(ioxy)                         &
                 +wmole(ih2o)*oxyh2o(ioxy)                        &
                 +wmole(ico2)*oxyco2(ioxy) )

        unsro1 = (wmolme*rr*timpat(izone)) / p0
        x1sro1 = (1.d0-x2tot) * unsro1
        brom(ifac) = 1.d0 / (x1sro1+x2sro2)
      endif
    endif

  enddo
endif
!--------
! Formats
!--------

!===============================================================================
! Deallocation dynamic arrays
!----
deallocate(f1m,f2m,f3m,f4m,f5m,STAT=iok1)
deallocate(f6m,f7m,f8m,f9m,    STAT=iok2)
deallocate(enth1,x2,fvp2m,     STAT=iok3)
!----
if ( iok1 > 0 .or. iok2 > 0 .or. iok3 > 0 ) then
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '     cs_fuel_physprop              '
  call csexit(1)
endif
if ( ieqnox .eq. 1 ) then
!----
  deallocate(xoxyd,enthox)
!----
  IF ( iok1 > 0 ) THEN
    write(nfecra,*) ' Memory deallocation error inside:      '
    write(nfecra,*) '   cs_fuel_physprop for xoxyd and enthox'
    call csexit(1)
  endif
endif
!===============================================================================

!----
! End
!----
return
end subroutine
