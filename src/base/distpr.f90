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

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   itypfb ,                                                       &
   ia     ,                                                       &
   distpa ,                                                       &
   viscf  , viscb  , dam    , xam    , smbdp  , rovsdp ,          &
   rtpdp  , coefad , coefbd ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , w7     , &
   w8     , w9     ,                                              &
   ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! itypfb           ! ia ! <-- ! boundary face types                            !
! ia(*)            ! ia ! --- ! main integer work array                        !
! distpa(ncelet    ! tr ! --> ! tab des distances a la paroi                   !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! smbrp(ncelet)    ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdp(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! rtpdp(ncelet)    ! tr ! --- ! var de travail du sclaire diffuse              !
! drtp(ncelet)     ! tr ! --- ! tableau de travail pour increment              !
! coefad,coefbd    ! tr ! --- ! conditions aux limites aux                     !
!  (nfabor)        !    !     ! faces de bord du scalaire diffuse              !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
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
use ppppar
use coincl
use parall
use period
use mltgrd
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas

integer          itypfb(nfabor)
integer          ia(*)

double precision distpa(ncelet), viscf (nfac)  , viscb (nfabor)
double precision dam   (ncelet), xam   (nfac,2)
double precision smbdp (ncelet), rovsdp(ncelet)
double precision rtpdp (ncelet)
double precision coefad(nfabor), coefbd(nfabor)
double precision w1    (ncelet), w2    (ncelet), w3    (ncelet)
double precision w4    (ncelet), w5    (ncelet), w6    (ncelet)
double precision w7    (ncelet), w8    (ncelet), w9    (ncelet)
double precision ra(*)

! Local variables


integer          idebia, idebra
integer          ndircp, iconvp, idiffp, isym
integer          ipol  , ireslp, ipp
integer          niterf, icycle, ncymxp, nitmfp
integer          iinvpe
integer          isqrt , iel   , ifac
integer          inc   , iccocg, iphydp, ivar
integer          isweep, nittot, idtva0

double precision relaxp, thetap, rnorm, residu, rnoini
double precision dismax, dismin
!===============================================================================



!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

rnoini = 0.d0

!     Memoire
idebia = idbia0
idebra = idbra0

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
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  ,                                              &
   ra     )

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
 ( nomvar(ipp)     , isym   , ipol   , ireslp , nitmay , imgrpy , &
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
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   idtva0 , ivar   , iconvp , idiffp , nswrgy , imligy , ircfly , &
   ischcy , isstpy , inc    , imrgra , iccocg ,                   &
   ipp    , iwarny ,                                              &
   blency , epsrgy , climgy , extray , relaxp , thetap ,          &
   ia     ,                                                       &
   rtpdp  , rtpdp  , coefad , coefbd , coefad , coefbd ,          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   smbdp  ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

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

!    - Echange pour le parallelisme et pour la periodicite

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtpdp)
  !==========
endif


!    - Calcul du gradient

inc    = 1
iccocg = 1
iphydp = 0
ivar   = 0

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   ivar   , imrgra , inc    , iccocg , nswrgy , imligy , iphydp , &
   iwarny , nfecra , epsrgy , climgy , extray ,                   &
   ia     ,                                                       &
   w1     , w2     , w3     ,                                     &
   rtpdp  , coefad , coefbd ,                                     &
   w4     , w5     , w6     ,                                     &
   w7     , w8     , w9     ,                                     &
   ra     )

do iel = 1, ncel
  w1(iel) = w4(iel)**2.d0+w5(iel)**2.d0+w6(iel)**2.d0
  if(w1(iel)+2.d0*rtpdp(iel).gt.0.d0) then
    distpa(iel) = - sqrt(w1(iel))                                 &
                  + sqrt(w1(iel)+2.d0*rtpdp(iel))
  else
    write(nfecra,8000)iel, xyzcen(1,iel)                          &
                     ,xyzcen(2,iel),xyzcen(3,iel)
  endif
enddo

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
