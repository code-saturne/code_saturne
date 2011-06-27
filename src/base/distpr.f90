!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine distpr &
!================

 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   distpa )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DE LA DISTANCE A LA PAROI PAR INVERSION DE LA SOLUTION 3D
!   DE L'EQUATION DE DIFFUSION D'UN SCALAIRE

!  On resout
!    div[grad(T)] = -1
!      avec :
!      T(bord)   = 1 en paroi
!      grad(T).n = 0 ailleurs

!  La distance a la paroi vaut alors :

!   d ~ -|grad(T)|+[grad(T).grad(T)+2.T]^(1/2)


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! distpa(ncelet    ! tr ! --> ! tab des distances a la paroi                   !
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
use ppppar
use coincl
use parall
use period
use mltgrd
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          itypfb(nfabor)

double precision distpa(ncelet)

! Local variables


integer          ndircp, iconvp, idiffp, isym
integer          ipol  , ireslp, ipp
integer          niterf, icycle, ncymxp, nitmfp
integer          iinvpe
integer          isqrt , iel   , ifac
integer          inc   , iccocg, ivar
integer          isweep, nittot, idtva0
integer          ibsize

double precision relaxp, thetap, rnorm, residu, rnoini
double precision dismax, dismin

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: coefad, coefbd
double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:,:) :: xam
double precision, allocatable, dimension(:) :: rtpdp, smbdp, rovsdp
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w7, w8, w9

!===============================================================================



!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Allocate temporary arrays for the species resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(coefad(nfac), coefbd(nfabor))
allocate(dam(ncelet), xam(nfac,2))
allocate(rtpdp(ncelet), smbdp(ncelet), rovsdp(ncelet))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))

! Initialize variables to avoid compiler warnings

rnoini = 0.d0

!     Memoire

!     Nombre d'iteration totale pour l'inversion
nittot = 0

!     La distance a la paroi est initialisee a 0 pour la reconstruction

do iel = 1, ncelet
  distpa(iel) = 0.d0
  rtpdp(iel)  = 0.d0
enddo

!===============================================================================
! 2.CONDITIONS LIMITES
!===============================================================================

!     Conditions aux limites pour le scalaire resolu T
!       Dirichlet a 0 en paroi
!       Neumann nul ailleurs
!     On test aussi la pressence d'un Dirichlet

ndircp = 0

do ifac = 1, nfabor
  if(itypfb(ifac).eq.iparoi .or.                            &
     itypfb(ifac).eq.iparug) then
    coefad(ifac) = 0.0d0
    coefbd(ifac) = 0.0d0
    ndircp = 1
  else
    coefad(ifac) = 0.0d0
    coefbd(ifac) = 1.0d0
  endif
enddo

!===============================================================================
! 3. PREPARATION DU SYSTEME A RESOUDRE
!===============================================================================

! -- Diagonale

do iel = 1, ncel
  rovsdp(iel) = 0.d0
enddo

! -- Diffusion aux faces

do iel = 1, ncel
  w1(iel) = 1.d0
enddo

call viscfa                                                       &
!==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

iconvp = 0
idiffp = 1
isym   = 1
thetap = 1.d0

call matrix                                                       &
!==========
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp ,                                     &
   isym   , nfecra ,                                              &
   thetap ,                                                       &
   ifacel , ifabor ,                                              &
   coefbd , rovsdp ,                                              &
   viscf  , viscb  , viscf  , viscb  ,                            &
   dam    , xam    )

! -- Second membre

do iel = 1, ncel
  smbdp(iel)  = volume(iel)
enddo
!===============================================================================
! 4. RESOLUTION DU SYSTEME
!===============================================================================
ipp = 1
NOMVAR(IPP) = 'DisParoi'
ipol   = 0
ireslp = 0
! Pas de multigrille (NCYMXP,NITMFP sont arbitraies)
ncymxp = 100
nitmfp = 10
! Periodicite
iinvpe = 0
if(iperio.eq.1) iinvpe = 1
isqrt = 1
ibsize = 1

! -- Boucle de reconstruction

! Si NSWRSY = 1, on doit faire 2 inversions
do isweep = 0, nswrsy

  call prodsc(ncelet,ncel,isqrt,smbdp,smbdp,rnorm)
  if (iwarny.ge.2) then
     write(nfecra,5000) nomvar(ipp),isweep,rnorm
  endif
  if (isweep.le.1) rnoini = rnorm
! Test de convergence
  if(rnorm.le.10.d0*epsily*rnoini) goto 100

  do iel = 1, ncelet
    rtpdp(iel) = 0.d0
  enddo

  call invers                                                     &
  !==========
 ( nomvar(ipp)     , isym   , ibsize ,                            &
   ipol   , ireslp , nitmay , imgrpy ,                            &
   ncymxp , nitmfp ,                                              &
   iwarny , nfecra , niterf , icycle , iinvpe ,                   &
   epsily , rnorm  , residu ,                                     &
   dam    , xam    , smbdp  , rtpdp  )

  nittot = nittot + niterf
  do iel = 1, ncel
    distpa(iel) = distpa(iel) + rtpdp(iel)
  enddo

!    - Echange pour le parallelisme

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtpdp)
    !==========
  endif

  if(isweep.lt.nswrsy) then
    inc    = 0
    iccocg = 1
    ivar = 0
    idtva0 = 0
    relaxp = 1.d0

    call bilsc2                                                   &
    !==========
 ( nvar   , nscal  ,                                              &
   idtva0 , ivar   , iconvp , idiffp , nswrgy , imligy , ircfly , &
   ischcy , isstpy , inc    , imrgra , iccocg ,                   &
   ipp    , iwarny ,                                              &
   blency , epsrgy , climgy , extray , relaxp , thetap ,          &
   rtpdp  , rtpdp  , coefad , coefbd , coefad , coefbd ,          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   smbdp  )

  endif
enddo

 100  continue

! On travail ensuite sur RTPDP pour calculer DISTPA
do iel=1,ncel
  rtpdp(iel)  = distpa(iel)
enddo

!===============================================================================
! 5. CALCUL DE LA DISTANCE A LA PAROI
!===============================================================================

! Allocate a temporary array for the gradient calculation
allocate(grad(ncelet,3))

!    - Echange pour le parallelisme et pour la periodicite

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtpdp)
  !==========
endif


!    - Calcul du gradient

inc    = 1
iccocg = 1
ivar   = 0

call grdcel                                                       &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgy , imligy ,          &
   iwarny , nfecra , epsrgy , climgy , extray ,                   &
   rtpdp  , coefad , coefbd ,                                     &
   grad   )

do iel = 1, ncel
  w1(iel) = grad(iel,1)**2.d0+grad(iel,2)**2.d0+grad(iel,3)**2.d0
  if(w1(iel)+2.d0*rtpdp(iel).gt.0.d0) then
    distpa(iel) = - sqrt(w1(iel)) + sqrt(w1(iel)+2.d0*rtpdp(iel))
  else
    write(nfecra,8000)iel, xyzcen(1,iel),xyzcen(2,iel),xyzcen(3,iel)
  endif
enddo

! Free memory
deallocate(grad)

!===============================================================================
! 6. CALCUL DES BORNES ET IMPRESSIONS
!===============================================================================
  dismax = -grand
  dismin =  grand

  do iel = 1, ncel
      dismin = min(distpa(iel),dismin)
      dismax = max(distpa(iel),dismax)
  enddo

  if (irangp.ge.0) then
    call parmin(dismin)
    call parmax(dismax)
  endif

!     Impressions

  write(nfecra,1000)dismin, dismax, nittot


! Free memory
deallocate(viscf, viscb)
deallocate(coefad, coefbd)
deallocate(dam, xam)
deallocate(rtpdp, smbdp, rovsdp)
deallocate(w1, w2, w3)
deallocate(w7, w8, w9)

!===============================================================================
! 7. FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'                                                             ',/,&
' ** DISTANCE A LA PAROI                                      ',/,&
'    -------------------                                      ',/,&
'                                                             ',/,&
'   Distance min = ',E14.5    ,'  Distance max = ',E14.5       ,/,&
'                                                             ',/,&
'     (Calcul de la distance realise en ',I10   ,' iterations)',/)

 5000 format(1X,A8,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6)

 8000   format(                                                         &
'@                                                            ',/,&
'@ @@ ATTENTION : Calcul de la distance a la paroi            ',/,&
'@    =========                                               ',/,&
'@  La variable associee ne converge pas a la cellule ',I10    ,/,&
'@       Coord X      Coord Y      Coord Z                    ',/,&
'@ ',3E13.5                                                    ,/)

#else

 1000 format(                                                           &
'                                                             ',/,&
' ** WALL DISTANCE                                            ',/,&
'    -------------                                            ',/,&
'                                                             ',/,&
'  Min distance = ',E14.5    ,' Max distance = ',E14.5         ,/,&
'                                                             ',/,&
'     (Distance calculation done in ',I10   ,' iterations)'    ,/)

 5000 format(1X,A8,' : SWEEP = ',I5,' RIGHT HAND SIDE NORM = ',E14.6)

 8000   format(                                                         &
'@                                                            ',/,&
'@ @@ WARNING: Wall distance calculation                      ',/,&
'@    ========                                                ',/,&
'@  The associated variable does not converge in cell ',I10    ,/,&
'@       Coord X      Coord Y      Coord Z                    ',/,&
'@ ',3E13.5                                                    ,/)

#endif

return
end subroutine
