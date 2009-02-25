!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine autmgr &
!================


    ( igr    , isym   , iagmax , nagmax ,                         &
      ncelf  , ncelfe , nfacf  , iwarnp ,                         &
      ifacef ,                                                    &
      daf    , xaf    , surfaf , volumf , xyzfin ,                &
      irscel ,                                                    &
      indic  , inombr , irsfac , indicf , w1     , w2 )

!===============================================================================
! FONCTION :
! ----------

!  MULTIGRILLE ALGEBRIQUE :
!  CONSTRUCTION D'UN NIVEAU DE MAILLAGE GROSSIER A PARTIR
!  DU NIVEAU SUPERIEUR SUIVANT CRITERE AUTOMATIQUE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nfecra           ! e  ! <-- ! numero du fichier d'impressions                !
! isym             ! e  ! <-- ! indicateur = 1 matrice sym                     !
!                  !    !     !            = 2 matrice non sym                 !
! igr              ! e  ! <-- ! niveau du maillage grossier                    !
! ncelf            ! e  ! <-- ! nombre d'elements maillage fin                 !
! ncelfe           ! e  ! <-- ! nombre d'elements etendus fin                  !
! nfacf            ! e  ! <-- ! nombre de faces internes maill. fin            !
! iwarnp           ! e  ! <-- ! niveau d'impression                            !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! ifacef           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfacf)       !    !     !  du maillage fin                               !
! daf(ncelfe)      ! tr ! <-- ! diagonale matrice maillage fin                 !
! xaf              ! tr ! <-- ! extradiagonale matrice maillage fin            !
! (nfacf, isym)    !    !     !                                                !
! surfaf           ! tr ! <-- ! surfaces faces internes maillage fin           !
! (3, nfacf)       !    !     !                                                !
! volumf           ! tr ! <-- ! volumes des cellules du maillage fin           !
! (ncelfe)         !    !     !                                                !
! xyzfin           ! tr ! <-- ! centres des cellules du maillage fin           !
! (3, ncelfe)      !    !     !                                                !
! irscel           ! te ! --> ! cellule fine -> cellule grossiere              !
!  (ncelfe)        !    !     !                                                !
! indic(ncelfe)    ! te ! --- ! tableau de travail                             !
! inombr           ! te ! --- ! tableau de travail                             !
!  (ncelfe)        !    !     !                                                !
! irsfac           ! te ! --- ! face fine -> face grossiere                    !
!  (nfacf)         !    !     !  (tableau de travail)                          !
! indicf(nfacf)    ! te ! --- ! indicateur de regroupement des faces           !
! icelfa           ! te ! --- ! connectivite cellules->faces mailla-           !
!  (2*nfacf)       !    !     ! ge fin                                         !
! icelce           ! te ! --- ! connectivite cellules->cellules                !
!  (2*nfacf)       !    !     ! voisines du maillage fin                       !
! rw(ncelf)        ! tr ! --- ! tableau de travail                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "entsor.h"
include "optcal.h"
include "parall.h"

!===============================================================================


!===============================================================================

! Arguments

integer          igr, isym, iagmax, nagmax
integer          ncelf, ncelfe, nfacf
integer          iwarnp

integer          ifacef(2, nfacf)
integer          irscel(ncelfe)
integer          indic(ncelfe), inombr(ncelfe)
integer          indicf(nfacf), irsfac(nfacf)

double precision daf(ncelfe), xaf(nfacf,isym)
double precision surfaf(3,nfacf), volumf(ncelfe)
double precision xyzfin(3,ncelfe)
double precision w1(ncelfe), w2(ncelfe)


! VARIABLES LOCALES

double precision critr, epslon

integer          ncelg, icel, ifac , ifac1, icelg
integer          nfacn,nfacnr,npass,npasmx
integer          inditt, noaglo, ngros, incvoi
integer          ihist(10)
integer          i, j, imin, imax


!===============================================================================

!     PARAMETRES PAR DEFAUT

epslon = +1.d-6
ngros  = 8
npasmx = 10
incvoi = 1

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

do icel = 1, ncelfe
  indic(icel) = -1
  irscel(icel) = 0
  inombr(icel) = 1
enddo
do ifac = 1, nfacf
  indicf(ifac) = ifac
  irsfac(ifac) = ifac
enddo

!     CALCUL DU CARDINAL (NOMBRE DE VOISINS DE CHAQUE CELLULE -1)
do ifac = 1, nfacf
  i = ifacef(1,ifac)
  j = ifacef(2,ifac)
  indic(i) = indic(i) + 1
  indic(j) = indic(j) + 1
enddo

ncelg  = 0
nfacnr = nfacf
npass  = 0
noaglo = ncelf

 100  continue

!     LES PASSES

npass = npass+1
nfacn = nfacnr
iagmax= iagmax +1
iagmax= min(iagmax, nagmax)

do ifac=1,nfacn
  irsfac(ifac) = indicf(ifac)
  indicf(ifac) = 0
enddo
if (nfacn .lt. nfacf) then
  do ifac = nfacn+1, nfacf
    indicf(ifac) = 0
    irsfac(ifac) = 0
  enddo
endif

if (iwarnp .gt. 3) then
  WRITE(NFECRA,*) '    autmgr.F : passe ', NPASS,                 &
                  'NFACNR = ', NFACNR, ' NOAGLO = ',NOAGLO
endif

!     INCREMENTATION DU NOMBRE DE VOISINS

do icel = 1, ncelf
  indic(icel) = indic(icel) + incvoi
enddo

!     INITIALISATION DE FACES NON ELIMINEES

nfacnr = 0

!     BOUCLE SUR LES FACES NON ELIMINEES

do ifac1 = 1, nfacn

  ifac = irsfac(ifac1)
  i = ifacef(1,ifac)
  j = ifacef(2,ifac)

!       ON NE CONSIDERE PAS LES FACES EN FRONTIERE PARALLELE OU
!       PERIODIQUE, POUR NE PAS AGGLOMERER LA GRILLE A TRAVERS
!       CE TYPE DE FRONTIERES (CE QUI DEMANDERAIT UNE CONSTRUCTION
!       PLUS COMPLEXE PUIS CHANGERAIT LE SCEMA DE COMMUNICATION).

  if (i.le.ncelf .and. j.le.ncelf) then

    inditt = 0
    critr  = (daf(i)/indic(i))*(daf(j)/indic(j))                  &
             /( xaf(ifac,1)*xaf(ifac,isym))

    if (       critr.lt.(1.d0-epslon)                             &
         .and. irscel(i)*irscel(j).le.0) then

      if (irscel(i).gt.0 .and. irscel(j).le.0) then
        if(inombr(irscel(i)) .le. iagmax) then
          irscel(j) = irscel(i)
          inombr(irscel(i)) = inombr(irscel(i)) +1
          inditt = inditt +1
        endif
      else if (irscel(i).le.0 .and. irscel(j).gt.0) then
        if (inombr(irscel(j)).le.iagmax) then
          irscel(i) = irscel(j)
          inombr(irscel(j)) = inombr(irscel(j)) + 1
          inditt = inditt +1
        endif
      else if (irscel(i).le.0 .and. irscel(j).le.0) then
        ncelg = ncelg+1
        irscel(i) = ncelg
        irscel(j) = ncelg
        inombr(ncelg) = inombr(ncelg) +1
        inditt = inditt +1
      endif

    endif

    if (inditt.ne.0 .and.inditt.ne.1) then
      WRITE(NFECRA,*)' Bug dans autmgr.F, arret '
      call csexit(1)
    endif

    if (inditt.eq.0 .and. irscel(i)*irscel(j).le.0) then
      nfacnr = nfacnr +1
      indicf(nfacnr) = ifac
    endif

  endif

enddo

!     CONTROLE DU NOMBRE DE CELLULES GROSSIERES CREE

noaglo = 0
do icel=1,ncelf
  if (irscel(icel).le.0) noaglo = noaglo+1
enddo

!     PASSES SUIVANTES SI AGGLOMERATION INSUFFISANTE

if (noaglo.gt.0) then
  if ((ncelg+noaglo)*ngros .ge. ncelf) then
    if (npass.lt.npasmx .and. nfacnr.gt.0) then
      goto 100
    endif
  endif
endif

!     ON TERMINE L'ASSEMBLAGE

do icel = 1, ncelf
  if (irscel(icel).le.0) then
    ncelg = ncelg+1
    irscel(icel) = ncelg
  endif
enddo

!     CONTROLE DIVERS ET VARIES

imax = 0
imin = 2*ncelf
do icelg =1,ncelg
  imax = max(imax, inombr(icelg))
  imin = min(imin, inombr(icelg))
enddo

if (irangp .ge. 0) then
  call parcmn(imin)
  call parcmx(imax)
endif

if (iwarnp.gt.3) then

  WRITE(NFECRA,*) '    autmgr.F : INOMBR MIN = ', IMIN,           &
                  ' MAX = ', IMAX, ' CIBLE = ', NAGMAX
  WRITE(NFECRA,*) '      histogramme '
  noaglo=imax-imin+1
  if (noaglo.gt.0) then
    if (noaglo.gt.10) then
      WRITE(NFECRA,*) ' IHIST MAL DIMENSIONNE DANS autmgr.F'
      call csexit(1)
    endif
    do i = 1, noaglo
      ihist(i) = 0
    enddo
    do icelg = 1, ncelg
      do i = 1, noaglo
        if (inombr(icelg).eq.(imin+i-1))then
          ihist(i)=ihist(i)+1
        endif
      enddo
    enddo
    if (irangp .ge. 0) then
      call parism(noaglo, ihist)
    endif
    do i = 1, noaglo
      epslon = 100.d0*ihist(i)/ncelg
      WRITE(NFECRA,*) '        regroupement ',IMIN+I-1,'  =  ',   &
                      EPSLON,' % '
    enddo
  endif

endif

do icel = 1, ncelf
  indic(icel) = 0
enddo
do icel = 1, ncelf
  icelg = irscel(icel)
  indic(icelg) = indic(icelg) +1
enddo

i=0
j=2*ncelf
noaglo = 0
do icelg = 1, ncelg
  i = max(i, indic(icelg))
  j = min(j, indic(icelg))
  noaglo = noaglo + indic(icelg)
enddo

if (irangp .ge. 0) then
  call parcmn(j)
  call parcmx(j)
endif

if (iwarnp.gt.3) then
  WRITE(NFECRA,*) '    autmgr.F : agglomeration MIN = ', J,       &
                  ' MAX= ',I
endif

if (noaglo .ne. ncelf) then
  WRITE(NFECRA,*) 'BUG DANS autmgr.f, Contacter l''assistance.'
  call csexit(1)
endif

!==============================================================================

!--------
! FORMATS
!--------


!----
! FIN
!----

return
end

