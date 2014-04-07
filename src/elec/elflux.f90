!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine elflux &
!================

 ( iappel , propce )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE POUR LE MODULE ELECTRIQUE

!              CALCULS DES VARIABLES PROPCE
!        ELLES SONT UTILISEES POUR LE CALCUL DES TS

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iappel           ! e  ! <-- ! numero d'appel                                 !
!                  !    !     ! 1 : j, e, j.e                                  !
!                  !    !     ! 2 : b,    jxb                                  !
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
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use elincl
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

integer          iappel

double precision propce(ncelet,*)

! Local variables

integer          iel
integer          ipcefj, ipcsig, ipcsii
integer          ipcla1, ipcla2, ipcla3
integer          ipcdc1, ipcdc2, ipcdc3
integer          ipcdi1, ipcdi2, ipcdi3
integer          iprev , inc   , iccocg
integer          ivar  , modntl

double precision vrmin, vrmax, var

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:,:,:) :: gradv

!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(grad(3 ,ncelet))
allocate(gradv(3, 3 ,ncelet))

! Initialize variables to avoid compiler warnings

ipcsii = 0
ipcdi1 = 0
ipcdi2 = 0
ipcdi3 = 0
ipcla1 = 0
ipcla2 = 0
ipcla3 = 0

! --- Numero des grandeurs physiques
ipcsig = ipproc(ivisls(ipotr))

if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
  ipcsii = ipproc(ivisls(ipoti))
endif

ipcefj = ipproc(iefjou)
ipcdc1 = ipproc(idjr(1))
ipcdc2 = ipproc(idjr(2))
ipcdc3 = ipproc(idjr(3))

if (ippmod(ieljou).eq.4) then
  ipcdi1 = ipproc(idji(1))
  ipcdi2 = ipproc(idji(2))
  ipcdi3 = ipproc(idji(3))
endif

! --- Necessite d'une impression (impression si MODNTL=0)
if (ntlist.gt.0) then
  modntl = mod(ntcabs,ntlist)
elseif (ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
  modntl = 0
else
  modntl = 1
endif

!===============================================================================
! 1. PREMIER APPEL : J, E => J.E
!===============================================================================

if (iappel.eq.1) then

!===============================================================================
! 1.1 PRISE EN COMPTE DES TERMES SOURCES ET VARIABLES STANDARD ET
!   COMMUNES A TOUTES LES VERSIONS ELECTRIQUES : EFFET JOULE, E, J
!         CAS IELJOU >= 1 ou IELARC >= 1 ou IELION >= 1
!===============================================================================

!    Pour toutes les versions electriques :
!      on doit calculer le terme joule (partie reelle) = j . E
!      produit de la densite de courant par le champ electrique


!   2.1 Calcul du grad (potR)
!  ---------------------------

  ivar = isca(ipotr)

  iprev = 0
  inc = 1
  iccocg = 1

  call field_gradient_scalar(ivarfl(ivar), iprev, imrgra, inc,     &
                             iccocg,                               &
                             grad)

!   2.2 Calcul du champ electrique E = - grad (potR)
!  -------------------------------------------------

!   2.3 Calcul de la densite de courant j = sig E  :
!  -------------------------------------------------
!                                   PROPCE(IEL,IPCSIG)


  if (ippmod(ieljou).ge.1 .or. ippmod(ielarc).ge.1) then
    do iel = 1, ncel
      propce(iel,ipcdc1)= - propce(iel,ipcsig) * grad(1,iel)
      propce(iel,ipcdc2)= - propce(iel,ipcsig) * grad(2,iel)
      propce(iel,ipcdc3)= - propce(iel,ipcsig) * grad(3,iel)
    enddo
  endif

!   2.4 Calcul de l'Effet Joule j . E : PROPCE( .,IPPROC(IEFJOU) )
!  -----------------------------------------------------------------
!                                                           sig E.E
  do iel = 1, ncel
    propce(iel,ipcefj)=                                           &
         propce(iel,ipcsig)*(grad(1,iel)**2+grad(2,iel)**2+grad(3,iel)**2)
  enddo


!   2.5 On imprime les extrema de E et j
!  -------------------------------------

  if (modntl.eq.0) then

    write(nfecra,1000)

    ! Grad PotR = -E
    var    = grad(1,1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var    = grad(1,iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    write(nfecra,1010)'Gr_PotRX',vrmin,vrmax

    var    = grad(2,1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var    = grad(2,iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    write(nfecra,1010)'Gr_PotRY',vrmin,vrmax

    var    = grad(3,1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var    = grad(3,iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    write(nfecra,1010)'Gr_PotRZ',vrmin,vrmax

    var    = -propce(1,ipcsig) * grad(1,1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var = -propce(iel,ipcsig) * grad(1,iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    write(nfecra,1010)'Cour_ReX',vrmin,vrmax

    var    = -propce(1,ipcsig) * grad(2,1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var = -propce(iel,ipcsig) * grad(2,iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    write(nfecra,1010)'Cour_ReY',vrmin,vrmax

    var    = -propce(1,ipcsig) * grad(3,1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var = -propce(iel,ipcsig) * grad(3,iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    write(nfecra,1010)'Cour_ReZ',vrmin,vrmax

  endif



!===============================================================================
! 1.2 PRISE EN COMPTE ET AJOUT DES TERMES SOURCES ET VARIABLES
!    RELATIVES A LA PRISE EN COMPTE DU POTENTIEL COMPLEXE
!                           CAS IELJOU = 2
!===============================================================================

  if (ippmod(ieljou).ge.2 .or. ippmod(ieljou).eq.4) then


!   3.1 Calcul du grad (potI) :
!  ----------------------------

    ivar = isca(ipoti)

    iprev = 0
    inc = 1
    iccocg = 1

    call field_gradient_scalar(ivarfl(ivar), iprev, imrgra, inc,     &
                               iccocg,                               &
                               grad)

!   3.2 Calcul du champ electrique Ei = - grad (potI) :
!  -------------------------------------------------

!   3.3 Partie imaginaire de la densite de courant :
!  ------------------------------------------------
!                                   PROPCE(IEL,IPCSIG)

  if ( ippmod(ieljou).eq.4 ) then
    do iel = 1, ncel
      propce(iel,ipcdi1)= -propce(iel,ipcsig)*grad(1,iel)
      propce(iel,ipcdi2)= -propce(iel,ipcsig)*grad(2,iel)
      propce(iel,ipcdi3)= -propce(iel,ipcsig)*grad(3,iel)
    enddo
  endif


!   3.4 Effet Joule total : PROPCE( .,IPPROC(IEFJOU) )
!  ----------------------------------------------------

    do iel = 1, ncel

!             ajout de la partie imaginaire et ...
      propce(iel,ipcefj) = propce(iel,ipcefj)                     &
         + propce(iel,ipcsii)*(grad(1,iel)**2+grad(2,iel)**2+grad(3,iel)**2)
!             .    ..division par 2
      propce(iel,ipcefj) = 0.5d0*propce(iel,ipcefj)

    enddo



!   3.5 On imprime les extrema de E et j
!  -------------------------------------

    if (modntl.eq.0) then

!     Grad PotI = -Ei
      var    = grad(1,1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = grad(1,iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      write(nfecra,1010)'Gr_PotIX',vrmin,vrmax

      var    = grad(2,1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = grad(2,iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      write(nfecra,1010)'Gr_PotIY',vrmin,vrmax

      var    = grad(3,1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = grad(3,iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      write(nfecra,1010)'Gr_PotIZ',vrmin,vrmax

!     j=sigma E
      var    = -propce(1,ipcsii) * grad(1,1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var = -propce(iel,ipcsii) * grad(1,iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      write(nfecra,1010)'Cour_ImX',vrmin,vrmax

      var    = -propce(1,ipcsii) * grad(2,1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var = -propce(iel,ipcsii) * grad(2,iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      write(nfecra,1010)'Cour_ImY',vrmin,vrmax

      var    = -propce(1,ipcsii) * grad(3,1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var = -propce(iel,ipcsii) * grad(3,iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      write(nfecra,1010)'Cour_ImZ',vrmin,vrmax

      write(nfecra,1001)

    endif

  endif

endif
!     Fin du test IAPPEL = 1

!===============================================================================
! 2. DEUXIEME APPEL : A, B, JXB
!===============================================================================

if (iappel.eq.2) then

!===============================================================================
! 2.1 PRISE EN COMPTE ET AJOUT DES TERMES SOURCES ET VARIABLES RELATIVES
!    A LA PRISE EN COMPTE DE L'ARC ELECTRIQUE
!                           CAS IELARC = 1 OU 2
!===============================================================================


  if (ippmod(ielarc).ge.1) then
    ipcla1 = ipproc(ilapla(1))
    ipcla2 = ipproc(ilapla(2))
    ipcla3 = ipproc(ilapla(3))
  endif

!   4.1 ARC ELECTRIQUE DIT 3D : IELARC = 2
!  ----------------------------------------
!         ON PASSE PAR LA RESOLUTION D'EQUATIONS DE POISSON SUR LES
!        -----------------------------------------------------------
!           COMPOSANTES DU POTENTIEL VECTEUR
!          ----------------------------------

  if (ippmod(ielarc).ge.2) then

! --> Calcul des composantes du champ magnetique B = (W1, W2, W3)
!===================================================


!    Sur A

    ivar = isca(ipotva)

    iprev = 0
    inc = 1
    iccocg = 1

    call field_gradient_vector(ivarfl(ivar), iprev, imrgra, inc,     &
                               gradv)

    ! B = rot A

    do iel = 1, ncel
      w1(iel)=  gradv(3, 2, iel) - gradv(2, 3, iel)
      w2(iel)=  gradv(1, 3, iel) - gradv(3, 1, iel)
      w3(iel)=  gradv(2, 1, iel) - gradv(1, 2, iel)
    enddo

!   4.2 CAS D'ARC AXISYMETRIQUE : IELARC = 1
!  ------------------------------------------
!         ON PEUT UTILISER LE THEOREME D'AMPERE
!        ---------------------------------------

  else if (ippmod(ielarc).eq.1) then

!      Calcul du Champ magnetique calcule :
!===========================================

!       B = / j d S

!       Cette version n'a pas ete developpe pour l'instant :
!         - elle ne fonctionne que dans le cas d'un arc axisymetrique
!         - elle demande donc un reperage des plans perpendiculaires
!           a l'arc
!         - elle necessite la connaissance des coordonnees de chaque
!           point du maillage

    write(nfecra,2000)
    call csexit (1)

  endif

!   4.3 ARC ELECTRIQUE : CALCUL DES FORCES DE LAPLACE j X B
!  ---------------------------------------------------------
!         POUR TOUS LES CAS D'ARCS ELECTRIQUES
!        --------------------------------------

  if (ippmod(ielarc) .ge. 1) then
    do iel = 1, ncel
      propce(iel,ipcla1)= propce(iel,ipcdc2) * w3(iel)            &
                         -propce(iel,ipcdc3) * w2(iel)
      propce(iel,ipcla2)= propce(iel,ipcdc3) * w1(iel)            &
                         -propce(iel,ipcdc1) * w3(iel)
      propce(iel,ipcla3)= propce(iel,ipcdc1) * w2(iel)            &
                         -propce(iel,ipcdc2) * w1(iel)
    enddo
  endif


!   4.4 Impression de B en arc electrique
!  --------------------------------------

  if (ippmod(ielarc) .ge. 2) then

    if (modntl.eq.0) then

      write(nfecra,1000)

!     B=rot A
      var    = w1(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = w1(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      write(nfecra,1010)'Ch_MagX ',vrmin,vrmax

      var    = w2(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = w2(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      write(nfecra,1010)'Ch_MagY ',vrmin,vrmax

      var    = w3(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = w3(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      write(nfecra,1010)'Ch_MagZ ',vrmin,vrmax

      write(nfecra,1001)

    endif

  endif

endif
!     Fin du test IAPPEL = 2

! Free memory
deallocate(w1, w2, w3)
deallocate(grad)
deallocate(gradv)

!--------
! Formats
!--------

 1000 format(/,                                                   &
'-----------------------------------------                    ',/,&
'   Variable         Minimum       Maximum                    ',/,&
'-----------------------------------------                    '  )
 1010 format(                                                           &
'v  ',A8,'    ',E12.5,'  ',E12.5                                 )
 1001 format(                                                           &
'-----------------------------------------                    '  )

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DANS ELTSSC            ',/,&
'@    =========                                               ',/,&
'@     LA VERSION ARC ELECTRIQUE AVEC THEOREME D''AMPERE      ',/,&
'@       (IELARC = 1)    N''EST PAS DISPONIBLE                ',/,&
'@                                                            ',/,&
'@            VEUILLEZ UTILISER LA VERSION DITE 3D            ',/,&
'@                      IELARC = 2                            ',/,&
'@   AVEC EQUATION DE TRANSPORT SUR LE POTENTIEL VECTEUR A    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return
end subroutine
