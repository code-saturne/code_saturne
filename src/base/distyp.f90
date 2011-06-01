!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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

subroutine distyp &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   ia     ,                                                       &
   distpa , propce , disty  ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DE LA DISTANCE A LA PAROI PAR UTILISATION D'UNE EQUATION
! DE TRANSPORT

! On resout dT/dt + grad (T U) = Gamma Ti

!    avec U = grad(distance a la paroi)/||grad(distance a la paroi)||
!         Gamma = div U
!         Ti = T (puits ou injection a la valeur locale)
!    et aux parois
!         T = u*/nu
!    et ailleurs
!         grad(T).n = 0

!    on cherche un etat stationnaire

!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! ia(*)            ! ia ! --- ! main integer work array                        !
! distpa(ncelet    ! tr ! <-- ! tab des distances a la paroi                   !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! disty(ncelet)    ! tr ! --> ! distance y+                                    !
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
use pointe, only: uetbor
use ppppar
use coincl
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          itypfb(nfabor)
integer          ia(*)

double precision distpa(ncelet),propce(ncelet,*)
double precision disty(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          idtva0, ivar  , iconvp, idiffp
integer          ndircp, ireslp
integer          iescap, iflmb0, imaspe
integer          ncymxp, nitmfp, ipp
integer          ifac  , iel   , ipcvis, init  , ipcrom
integer          inc   , iccocg, isym  , ntcont, infpar

double precision xnorme, dtminy, dtmaxy, relaxp, thetap, timey
double precision xusnmx, xusnmn, xnorm0
double precision dismax, dismin, usna

double precision rvoid(1)

double precision, allocatable, dimension(:) :: rtpdp, smbdp, rovsdp
double precision, allocatable, dimension(:) :: qx, qy, qz
double precision, allocatable, dimension(:) :: flumas, flumab
double precision, allocatable, dimension(:) :: rom, romb
double precision, allocatable, dimension(:) :: coefax, coefay, coefaz
double precision, allocatable, dimension(:) :: coefbx, coefby, coefbz
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w2

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Allocate temporary arrays for the distance resolution
allocate(rtpdp(ncelet), smbdp(ncelet), rovsdp(ncelet))
allocate(qx(ncelet), qy(ncelet), qz(ncelet))
allocate(flumas(nfac), flumab(nfabor))
allocate(rom(nfac), romb(nfabor))
allocate(coefax(nfabor), coefay(nfabor), coefaz(nfabor))
allocate(coefbx(nfabor), coefby(nfabor), coefbz(nfabor))

! Allocate work arrays
allocate(w2(ncelet))

idebia = idbia0
idebra = idbra0

ipass  = ipass + 1

ipcrom = ipproc(irom  )
ipcvis = ipproc(iviscl)

!===============================================================================
! 2. AU PREMIER PAS DE TEMPS
!===============================================================================

!     Au premier pas de temps, on a en general u* = 0 (ou faux)
!       on ne calcule pas y+

!     En effet ca prend du temps, d'autant plus que u* est petit, car il
!       alors calculer y+ jusqu'a une grande distance des parois

if(ntcabs.eq.1) then

  do iel = 1, ncel
    disty(iel) = grand
  enddo

  if(iwarny.ge.1) then
    write(nfecra,7000)
  endif

  return

endif

!===============================================================================
! 3. CALCUL DE Q=Grad(DISTPA)/|Grad(DISTPA)|
!===============================================================================

!     La distance a la paroi vaut 0 en paroi
!       par definition et obeit a un flux nul ailleurs

do ifac = 1, nfabor
  if(itypfb(ifac).eq.iparoi .or.                            &
     itypfb(ifac).eq.iparug ) then
    coefax(ifac) = 0.0d0
    coefbx(ifac) = 0.0d0
  else
    coefax(ifac) = 0.0d0
    coefbx(ifac) = 1.0d0
  endif
enddo

! Allocate a temporary array for the gradient calculation
allocate(grad(ncelet,3))

!     Calcul du gradient

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(distpa)
  !==========
endif

inc    = 1
iccocg = 1
ivar   = 0

call grdcel                                                       &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgy , imligy ,          &
   iwarny , nfecra , epsrgy , climgy , extray ,                   &
   ia     ,                                                       &
   distpa , coefax , coefbx ,                                     &
   grad   ,                                                       &
   ra     )


!     Normalisation (attention, le gradient peut etre nul, parfois)

do iel = 1, ncel
  xnorme = max(sqrt(grad(iel,1)**2+grad(iel,2)**2+grad(iel,3)**2),epzero)
  qx(iel) = grad(iel,1)/xnorme
  qy(iel) = grad(iel,2)/xnorme
  qz(iel) = grad(iel,3)/xnorme
enddo

! Free memory
deallocate(grad)

!===============================================================================
! 4. CALCUL DU FLUX DE Q ET DE GAMMA = div(Q)
!===============================================================================


do ifac = 1, nfabor
  romb(ifac) = 1.d0
enddo
do iel = 1, ncelet
  rom(iel)  = 1.d0
enddo


!     Le gradient normal de la distance a la paroi vaut -1 en paroi
!       par definition et obeit a un flux nul ailleurs

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
    xnorme = max(surfbn(ifac),epzero**2)
    coefax(ifac) = -surfbo(1,ifac)/xnorme
    coefbx(ifac) = 0.d0
    coefay(ifac) = -surfbo(2,ifac)/xnorme
    coefby(ifac) = 0.d0
    coefaz(ifac) = -surfbo(3,ifac)/xnorme
    coefbz(ifac) = 0.d0
  else
    coefax(ifac) = 0.d0
    coefbx(ifac) = 1.d0
    coefay(ifac) = 0.d0
    coefby(ifac) = 1.d0
    coefaz(ifac) = 0.d0
    coefbz(ifac) = 1.d0
  endif
enddo


!     Parallelisme et periodicite en preparation du calcul du flux

if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(qx, qy, qz)
  !==========
endif


!     Calcul du flux de masse

!     On ne le met pas a zero en paroi (justement)
iflmb0 = 0
!     On l'initialise a 0
init   = 1
!     On prend en compte les Dirichlet
inc    = 1
!     On recalcule les gradients complets
iccocg = 1
!     Calcul de flux std (pas de divrij)
imaspe = 1
!     Il ne s'agit ni de U ni de R
ivar = 0

call inimas                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ivar   , ivar   , ivar   , imaspe ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgy , imligy , &
   iwarny , nfecra ,                                              &
   epsrgy , climgy , extray ,                                     &
   ia     ,                                                       &
   rom    , romb   ,                                              &
   qx     , qy     , qz     ,                                     &
   coefax , coefay , coefaz , coefbx , coefby , coefbz ,          &
   flumas , flumab ,                                              &
   ra     )



!     A partir d'ici, QX, QY et QZ sont des tableaux de travail.



!     Calcul de la divergence du flux de masse

init = 1

call divmas                                                       &
!==========
 ( ncelet , ncel   , nfac   , nfabor , init   , nfecra ,          &
   ifacel , ifabor , flumas , flumab , qx)

!===============================================================================
! 5.CONDITIONS LIMITES
!===============================================================================

!     Dirichlet en u*/nu aux parois, et flux nul ailleurs

do ifac = 1, nfabor
  if(itypfb(ifac).eq.iparoi .or.                            &
     itypfb(ifac).eq.iparug) then
    iel = ifabor(ifac)
    coefax(ifac) = uetbor(ifac)*propce(iel,ipcrom)/propce(iel,ipcvis)
    coefbx(ifac) = 0.0d0
  else
    coefax(ifac) = 0.0d0
    coefbx(ifac) = 1.0d0
  endif
enddo

!===============================================================================
! 6.CALCUL DU PAS DE TEMPS
!===============================================================================

!     On vise un Courant infini (de l'ordre de 1000).

!     On calcule avec MATRDT DA = Sigma a S/d
iconvp = 1
idiffp = 0
!     La matrice est non symetrique
isym   = 2

call matrdt                                                       &
!==========
 ( idebia , idebra ,                                              &
   iconvp , idiffp , isym   ,                                     &
   ia     ,                                                       &
   coefbx , flumas , flumab , flumas , flumab , w2     ,          &
   ra     )

!     Le Courant est COUMXY = DT W2 / VOLUME
!         d'ou DTMINY = MIN(COUMXY * VOLUME/W2)
!     Au cas ou une cellule serait a W2(IEL)=0,
!       on ne la prend pas en compte

!     On prend dans QZ un pas de temps variable,
!       si on ne trouve pas de pas de temps, on prend le minimum
!       (cellules sources)
dtminy = grand
dtmaxy = -grand
do iel = 1, ncel
  qz(iel) = -grand
  if(w2(iel).gt.epzero) then
    qz(iel) = coumxy*volume(iel)/w2(iel)
    dtminy  = min(qz(iel),dtminy)
    dtmaxy  = max(qz(iel),dtmaxy)
  endif
enddo
if(irangp.ge.0) then
  call parmin (dtminy)
  call parmax (dtmaxy)
endif
dtminy = max(dtminy,epzero)

do iel = 1, ncel
  if(qz(iel).le.0.d0) then
    qz(iel) = dtminy
  endif
enddo

if(iwarny.ge.2) then
  write(nfecra,2000)dtminy,dtmaxy
endif

!===============================================================================
! 7. DIAGONALE DE LA MATRICE
!===============================================================================

do iel = 1, ncel
  rovsdp(iel) = volume(iel)*rom(iel)/qz(iel)-qx(iel)
enddo

!===============================================================================
! 8. BOUCLE EN TEMPS
!===============================================================================

!     Initialisations
!    =================

!     Iterations
ntcont = 0

!     Temps
timey = 0.d0



!     Inconnue
!       Au cas ou on n'atteint pas tout a fait l'etat stationnaire,
!       il faut que le yplus ne soit pas nul dans la zone ou les
!       conditions aux limites n'ont pas ete convectees. On voudrait
!       plutot que yplus y soit maximum.
!       Si on utilise zero ou une valeur negative comme initialisation,
!       on risque de se retrouver avec des valeurs proches de
!       zero issues de la diffusion due au schema upwind au voisinage
!       du front convecte et donc avec des yplus proches de zero
!       n'importe ou.
!       On va donc utiliser la valeur max de u*/nu.

!       A partir du second pas de temps, on a egalement le yplus du pas
!         de temps precedent

!     On calcule le min et le max
xusnmx = -grand
xusnmn =  grand
do ifac = 1, nfabor
  if(itypfb(ifac).eq.iparoi .or.                            &
     itypfb(ifac).eq.iparug) then
    xusnmx = max(xusnmx,coefax(ifac))
    xusnmn = min(xusnmn,coefax(ifac))
  endif
enddo
if(irangp.ge.0) then
  call parmax (xusnmx)
  call parmin (xusnmn)
endif

if(ipass.eq.1) then
  do iel = 1, ncelet
    rtpdp(iel) = xusnmx
  enddo
else
  do iel = 1, ncel
    usna = disty(iel)/max(distpa(iel),epzero)
    usna = max(usna,xusnmn)
    usna = min(usna,xusnmx)
    rtpdp(iel) = usna
  enddo
endif

!     Norme de reference (moyenne des u*/nu)
!       (on divise par le nombre de faces)
xnorm0 = 0.d0
infpar = 0
do ifac = 1, nfabor
  if(itypfb(ifac).eq.iparoi .or.                            &
     itypfb(ifac).eq.iparug) then
    infpar = infpar+1
    xnorm0 = xnorm0 + coefax(ifac)**2
  endif
enddo
if(irangp.ge.0) then
  call parcpt (infpar)
  call parsom (xnorm0)
endif
xnorm0 = xnorm0/dble(infpar)

!     Pour ne pas diviser par zero
if(xnorm0.le.epzero**2) goto 100


!     Debut des iterations
!    ======================

do ntcont = 1, ntcmxy

!     Instant (arbitrairement +DTMINY)
  timey = timey + dtminy

! -- Echange pour les cas paralleles
!       a la premiere iteration, c'est inutile (on a fait l'init sur NCELET)

  if(ntcont.gt.1.or.ipass.gt.1) then

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtpdp)
      !==========
    endif

  endif


!     Sauvegarde de la solution pour test de convergence
!    ====================================================

  do iel = 1, ncel
    qy(iel) = rtpdp(iel)
  enddo

!     Second membre
!    ===============
!       Obligatoirement a tous les pas de temps

  do iel = 1, ncel
    smbdp(iel) = qx(iel)*rtpdp(iel)
  enddo

!     Resolution
!    ============

!     La variable n'est pas la vitesse ou une composante de Rij
  ivar = 0
!     Le cas est convectif, non diffusif
  iconvp = 1
  idiffp = 0
!     Il y a des Dirichlet (car il y a des parois)
  ndircp = 1
!     On resout par la methode automatique
  ireslp = -1
!     Pas d'estimateurs, ni de multigrille (100 et 10 sont arbitraires)
  iescap = 0
  ncymxp = 100
  nitmfp = 10
!     La case 1 est une poubelle
  ipp    = 1
  NOMVAR(IPP) = 'YplusPar'
!     Ordre 1 en temps (etat stationnaire cherche)
  thetap = 1.d0
!     Pas de stationnaire ni de relaxation -> a modifier eventuellement
  idtva0 = 0
  relaxp = 1.d0

  call codits                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   idtva0 , ivar   , iconvp , idiffp , ireslp , ndircp , nitmay , &
   imrgra , nswrsy , nswrgy , imligy , ircfly ,                   &
   ischcy , isstpy , iescap ,                                     &
   imgrpy , ncymxp , nitmfp , ipp    , iwarny ,                   &
   blency , epsily , epsrsy , epsrgy , climgy , extray ,          &
   relaxp , thetap ,                                              &
   ia     ,                                                       &
   rtpdp  , rtpdp  ,                                              &
   coefax , coefbx , coefax , coefbx , flumas , flumab ,          &
   flumas , flumab , flumas , flumab ,                            &
   rovsdp , smbdp  , rtpdp  ,                                     &
   rvoid  ,                                                       &
   ra     )


!     Clipping (indispensable si on initialise par u*/nu du pas de
!    ==========                                   temps precedent)

  do iel = 1, ncel
    rtpdp(iel) = max(rtpdp(iel),xusnmn)
    rtpdp(iel) = min(rtpdp(iel),xusnmx)
  enddo

!     Test d'arret
!    ==============

!     on utilise QY dans lequel la solution precedente a ete sauvee
!     on souhaite que la variation de l'inconnue sur chaque cellule
!       soit inferieure a un pourcentage de la valeur moyenne de
!       u*/nu calculee sur les faces de bord
!     on limite le test aux cellules qui pourront avoir un interet pour
!       VanDriest, c'est a dire qu'on ignore celles a y+ > YPLMXY
!       comme on ne connait pas encore y+, on se base sur min(u*/nu) :
!       on ignore les cellules a y min(u*/nu) > YPLMXY

!        XNORME = 0.D0
!        DO IEL = 1, NCEL
!          XNORME = XNORME + (RTPDP(IEL)-QY(IEL))**2
!        ENDDO
!        IF(IRANGP.GE.0) THEN
!          CALL PARSOM (XNORME)
!        ENDIF
!        XNORME = XNORME/DBLE(NCELGB)

  xnorme = -grand
  do iel = 1, ncel
    if(distpa(iel)*xusnmn.le.yplmxy) then
      xnorme = max(xnorme,(rtpdp(iel)-qy(iel))**2)
    endif
  enddo
  if(irangp.ge.0) then
    call parmax (xnorme)
  endif

  if(iwarny.ge.2) then
    write(nfecra,3000)ntcont,xnorme,xnorm0,xnorme/xnorm0
  endif

  if(xnorme.le.epscvy*xnorm0) goto 100

enddo

write(nfecra,8000)xnorme,xnorm0,xnorme/xnorm0,ntcmxy

 100  continue

!===============================================================================
! 9. CALCUL DE YPLUS ET IMPRESSIONS
!===============================================================================


do iel = 1, ncel
  disty(iel) = rtpdp(iel)*distpa(iel)
enddo

dismax = -grand
dismin =  grand

do iel = 1, ncel
  dismin = min(disty(iel),dismin)
  dismax = max(disty(iel),dismax)
enddo

if (irangp.ge.0) then
  call parmin(dismin)
  call parmax(dismax)
endif

if(iwarny.ge.1) then
  write(nfecra,1000)dismin, dismax, min(ntcont,ntcmxy)
endif


! Free memory
deallocate(rtpdp, smbdp, rovsdp)
deallocate(qx, qy, qz)
deallocate(flumas, flumab)
deallocate(rom, romb)
deallocate(coefax, coefay, coefaz)
deallocate(coefbx, coefby, coefbz)
deallocate(w2)


#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'                                                             ',/,&
' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,&
'    ------------------------------------                     ',/,&
'                                                             ',/,&
'  Distance+ min = ',E14.5    ,' Distance+ max = ',E14.5       ,/,&
'                                                             ',/,&
'     (Calcul de la distance realise en ',I10   ,' iterations)',/)
 2000 format(                                                           &
'                                                             ',/,&
' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,&
'    ------------------------------------                     ',/,&
'                                                             ',/,&
' Yplus:  Dt min = ',E14.5    ,'        Dt max = ',E14.5       ,/)
 3000 format(                                                           &
'                                                             ',/,&
' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,&
'    ------------------------------------                     ',/,&
'                                                             ',/,&
' Yplus:  iteration   residu abs.     reference   residu rel. ',/,&
' Yplus: ',I10    ,E14.5        ,E14.5        ,E14.5           ,/)
 7000 format(                                                           &
'                                                             ',/,&
' ** DISTANCE A LA PAROI ADIMENSIONNELLE                      ',/,&
'    ------------------------------------                     ',/,&
'                                                             ',/,&
'  Elle n''est pas calculee au premier pas de temps           ',/)
 8000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : Calcul de la distance a la paroi            ',/,&
'@    =========                                               ',/,&
'@     Le systeme iteratif pour le calcul de la distance a la ',/,&
'@       paroi adimensionnelle est imparfaitement converge.   ',/,&
'@                                                            ',/,&
'@          Residu     Reference            Residu relatif    ',/,&
'@  ',2E14.5,12X,E14.5                                         ,/,&
'@                                                            ',/,&
'@     Augmenter la valeur de NTCMXY dans usini1 peut resoudre',/,&
'@       le probleme.                                         ',/,&
'@     La valeur actuelle de cet entier est ',I10              ,/,&
'@                                                            ',/,&
'@     Le calcul se poursuit.                                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                           &
'                                                             ',/,&
' ** DIMENSIONLESS WALL DISTANCE                              ',/,&
'    ---------------------------                              ',/,&
'                                                             ',/,&
'  Min distance+ = ',E14.5    ,' Max distance+ = ',E14.5       ,/,&
'                                                             ',/,&
'     (Distance calculation done in ',I10   ,' iterations)'    ,/)
 2000 format(                                                           &
'                                                             ',/,&
' ** DIMENSIONLESS WALL DISTANCE                              ',/,&
'    ---------------------------                              ',/,&
'                                                             ',/,&
' Yplus:  Min dt = ',E14.5    ,'        Max dt = ',E14.5       ,/)
 3000 format(                                                           &
'                                                             ',/,&
' ** DIMENSIONLESS WALL DISTANCE                              ',/,&
'    ---------------------------                              ',/,&
'                                                             ',/,&
' Yplus:  iteration   abs. residu     reference   rel. residu ',/,&
' Yplus: ',I10    ,E14.5        ,E14.5        ,E14.5           ,/)
 7000 format(                                                           &
'                                                             ',/,&
' ** DIMENSIONLESS WALL DISTANCE                              ',/,&
'    ---------------------------                              ',/,&
'                                                             ',/,&
'  It is not computed at the first time step                  ',/)
 8000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Dimensionless distance to the wall computation ',/,&
'@    ========                                                ',/,&
'@     The iterative system for the dimensionless distance to ',/,&
'@       the wall calculation has not properly converged.     ',/,&
'@                                                            ',/,&
'@          Residual   Reference            Relative residual ',/,&
'@  ',2E14.5,12X,E14.5                                         ,/,&
'@                                                            ',/,&
'@     Increase the value of NTCMXY in usini1 may solve       ',/,&
'@       the problem.                                         ',/,&
'@     The current value of this integer is ',I10              ,/,&
'@                                                            ',/,&
'@     The calculation will be run.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif



return
end subroutine
