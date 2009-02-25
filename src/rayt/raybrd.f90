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

subroutine raybrd &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
   nvar   , nscal  , nphas  , iph    ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   itypfb , itrifb , izfrad ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   surfbn ,                                                       &
   tparoi , qincid , xlam   , epa    , eps    ,                   &
   flunet , flconv , hfconv , w1     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
!   -----------------------------------------

! Sortie au format CASE de grandeurs
! sur les faces de frontieres du domaine de calcul.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nnod             ! e  ! <-- ! nombre de sommets                              !
! lndnod           ! e  ! <-- ! longueur du tableau icocel (optionnel          !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr (optionnel          !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! iph              ! e  ! <-- ! numero de la phase courante associee           !
!                  !    !     ! au rayonnement                                 !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
!                  !    !     ! le module lagrangien                           !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! nfml  ,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas)          !    !     !                                                !
! itrifb(nfabor    ! te ! --> ! tab d'indirection pour tri des faces           !
!  nphas)          !    !     !                                                !
! izfrad(nfabor    ! te ! <-- ! numero de zone des faces de bord               !
! isothm(nfabor    ! te ! <-- ! type de condition de paroi                     !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds (optionnel)              !
! (ndim,nnod)      !    !     !                                                !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! surfbn(nfabor    ! tr ! <-- ! surface des faces de bord                      !
! tparoi(nfabor    ! tr ! <-- ! temperature de paroi en kelvin                 !
!   ,nphast)       !    !     !                                                !
! qincid(nfabor    ! tr ! <-- ! densite de flux radiatif aux bords             !
!   ,nphast)       !    !     !                                                !
! xlam(nfabor      ! tr ! <-- ! coefficient de conductivite thermique          !
!   ,nphast)       !    !     ! des facettes de paroi (w/m/k)                  !
! epa (nfabor      ! tr ! <-- ! epaisseur des facettes de paroi (m)            !
!   ,nphast)       !    !     !                                                !
! eps (nfabor      ! tr ! <-- ! emissivite des facettes de bord                !
!   ,nphast)       !    !     !                                                !
! flunet(nfabor    ! tr ! <-- ! densite de flux net radiatif aux               !
!   ,nphast)       !    !     ! faces de bord                                  !
! flconv(nfabor    ! tr ! <-- ! densite de flux convectif aux faces            !
!   ,nphast)       !    !     ! de bord                                        !
! hfconv(nfabor    ! tr ! <-- ! coefficient d'echange fluide aux               !
!   ,nphast)       !    !     ! faces de bord                                  !
! w1(nfabor)       ! tr ! <-- ! tableau de travail                             !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
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
include "numvar.h"
include "optcal.h"
include "entsor.h"
include "cstphy.h"
include "pointe.h"
include "period.h"
include "parall.h"
include "radiat.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndnod , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , iph
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          itypfb(nfabor,nphas) , itrifb(nfabor,nphas)
integer          izfrad(nfabor,nphas)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision surfbn(nfabor)
double precision tparoi(nfabor,nphast), qincid(nfabor,nphast)
double precision xlam(nfabor,nphast), epa(nfabor,nphast)
double precision eps(nfabor,nphast), flunet(nfabor,nphast)
double precision flconv(nfabor,nphast), hfconv(nfabor,nphast)
double precision w1(nfabor)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia , idebra

integer          ifac , nbp , ivar , iphas
integer          ii , jj , in , kpt
integer          ipart , ntria3 , nquad4 , nsided

integer          iz , iok , izone , nbzone
integer          izones(nozrdm) , impray

integer          ii1 , ii2 , lpos , n1 , n2
character        fich*80 , name*80 , entete*7

!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1.  A t-on quelquechose a voir sur le maillage frontiere ?
!===============================================================================

iok = 0
iphas = irapha(iph)
do iz = 1,nbrayf
  if (irayvf(iz,iphas).eq.1) iok = iok + 1
enddo
if (iok.eq.0) return

ENTETE = ' '
WRITE(ENTETE,'(A5,I1,A1)')'bord_',IPHAS,'.'

impray = 20

!===============================================================================
! 2.  ECRITURE DU MAILLAGE DE PEAU AU FORMAT ENSIGHT GOLD
!===============================================================================

!--> On stocke le nombre de zones et leur numero associe

nbzone = 0
do iz = 1,nozrdm
  izones(iz) = -1
enddo

do ifac = 1,nfabor

  izone = izfrad(ifac,iph)

  if (izone.gt.0) then

    iok = 0
    do iz = 1,nbzone+1
      if (izone.ne.izones(iz) .and. iok.eq.0) then
        iok = 1
      else if (izones(iz).eq.izone) then
        iok = 2
      endif
    enddo

    if (iok.eq.1) then
      nbzone = nbzone +1
      if (nbzone.gt.nozrdm) then
        write(nfecra,9010)
        return
      endif
      izones(nbzone) = izone
    endif

  endif

enddo

if (nbzone.eq.0) then
  write(nfecra,9020)
  return
endif

!===============================================================================
! 3.  ECRITURE DU MAILLAGE DE PEAU AU FORMAT ENSIGHT GOLD
!===============================================================================

!--> Ouverture

OPEN (IMPRAY,FILE= ENTETE // 'geom',                              &
      STATUS='UNKNOWN',FORM='FORMATTED',                          &
      ACCESS='SEQUENTIAL')

!--> Entete

WRITE(IMPRAY,'(A)') 'geometrie de la frontiere'
WRITE(IMPRAY,'(A)') 'au format ensight gold'
WRITE(IMPRAY,'(A)') 'node id given'
WRITE(IMPRAY,'(A)') 'element id given'

ipart = 0

do iz = 1,nbzone

  izone = izones(iz)

  nquad4 = 0
  ntria3 = 0
  nsided = 0

  do ifac = 1,nfabor
    if (izfrad(ifac,iph).eq.izone) then
      nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
      if (nbp.eq.4) then
        nquad4 = nquad4 + 1
      else if (nbp.eq.3) then
        ntria3 = ntria3 + 1
      else if (nbp.ge.5) then
        nsided = nsided + 1
      endif
    endif
  enddo


!-> Part des triangles

  if (ntria3.gt.0) then

    ipart = ipart+1
    WRITE(IMPRAY,'(A)')    'part'
    WRITE(IMPRAY,'(I10)')  IPART
    WRITE(IMPRAY,'(A,I4)') 'faces triangles zone ',IZ
    WRITE(IMPRAY,'(A)')    'coordinates'
    WRITE(IMPRAY,'(I10)')  NTRIA3*3

    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        if (nbp.eq.3) then
          do in = ipnfbr(ifac),ipnfbr(ifac+1)-1
            WRITE(IMPRAY,'(I10)')  NODFBR(IN)
          enddo
        endif
      endif
    enddo

    do in   = 1,3
      do ifac = 1,nfabor
        if (izfrad(ifac,iph).eq.izone) then
          nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
          if (nbp.eq.3) then
            do ii = ipnfbr(ifac),ipnfbr(ifac+1)-1
              WRITE(IMPRAY,'(E12.5)') XYZNOD(IN,NODFBR(II))
            enddo
          endif
        endif
      enddo
    enddo

    WRITE(IMPRAY,'(A)')  'tria3'
    WRITE(IMPRAY,'(I10)') NTRIA3
    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        IF (NBP.EQ.3) WRITE(IMPRAY,'(I10)') IFAC
      endif
    enddo

    kpt = 0
    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        if (nbp.eq.3) then
          WRITE(IMPRAY,'(3I10)') KPT+1,KPT+2,KPT+3
          kpt = kpt+3
        endif
      endif
    enddo

! Fin ecriture part des triangles
  endif


!-> Part des quadrangles

  if (nquad4.gt.0) then

    ipart = ipart+1
    WRITE(IMPRAY,'(A)')    'part'
    WRITE(IMPRAY,'(I10)')  IPART
    WRITE(IMPRAY,'(A,I4)') 'faces quadrangles zone ',IZ
    WRITE(IMPRAY,'(A)')    'coordinates'
    WRITE(IMPRAY,'(I10)')  NQUAD4*4

    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        if (nbp.eq.4) then
          do in = ipnfbr(ifac),ipnfbr(ifac+1)-1
            WRITE(IMPRAY,'(I10)') NODFBR(IN)
          enddo
        endif
      endif
    enddo

    do in   = 1,3
      do ifac = 1,nfabor
        if (izfrad(ifac,iph).eq.izone) then
          nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
          if (nbp.eq.4) then
            do jj = ipnfbr(ifac),ipnfbr(ifac+1)-1
              WRITE(IMPRAY,'(E12.5)') XYZNOD(IN,NODFBR(JJ))
            enddo
          endif
        endif
      enddo
    enddo

    WRITE(IMPRAY,'(A)')  'quad4'
    WRITE(IMPRAY,'(I10)') NQUAD4
    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        IF (NBP.EQ.4) WRITE(IMPRAY,'(I10)') IFAC
      endif
    enddo

    kpt = 0
    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        if (nbp.eq.4) then
          WRITE(IMPRAY,'(4I10)') KPT+1,KPT+2,KPT+3,KPT+4
          kpt = kpt+4
        endif
      endif
    enddo

! Fin ecriture part des quadrangles
  endif

!--> Geometrie : nsided

  if (nsided.gt.0) then

    ipart = ipart+1
    WRITE(IMPRAY,'(A)')    'part'
    WRITE(IMPRAY,'(I10)')  IPART
    WRITE(IMPRAY,'(A,I4)') 'faces polygones zone ',IZ
    WRITE(IMPRAY,'(A)')    'coordinates'

    kpt = 0
    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        if (nbp.ge.5) kpt = kpt + nbp
      endif
    enddo
    WRITE(IMPRAY,'(I10)')  KPT

    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        if (nbp.ge.5) then
          do in = ipnfbr(ifac),ipnfbr(ifac+1)-1
            WRITE(IMPRAY,'(I10)') NODFBR(IN)
          enddo
        endif
      endif
    enddo

    do in   = 1,3
      do ifac = 1,nfabor
        if (izfrad(ifac,iph).eq.izone) then
          nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
          if (nbp.ge.5) then
            do jj = ipnfbr(ifac),ipnfbr(ifac+1)-1
              WRITE(IMPRAY,'(E12.5)') XYZNOD(IN,NODFBR(JJ))
            enddo
          endif
        endif
      enddo
    enddo

    WRITE(IMPRAY,'(A)')  'nsided'
    WRITE(IMPRAY,'(I10)') NSIDED
    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        IF (NBP.GE.5) WRITE(IMPRAY,'(I10)') IFAC
      endif
    enddo

    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        IF (NBP.GE.5) WRITE(IMPRAY,'(I10)') NBP
      endif
    enddo

    kpt = 0
    do ifac = 1,nfabor
      if (izfrad(ifac,iph).eq.izone) then
        nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
        if (nbp.ge.5) then
          WRITE(IMPRAY,'(100I10)') (KPT+II,II = 1,NBP)
          kpt = kpt+nbp
        endif
      endif
    enddo

! Fin ecriture Nsided
  endif

! Fin boucle sur les zones
enddo

close(impray)

!===============================================================================
! 4.  ECRITURE FICHIERS VARIABLES AU FORMAT ENSIGHT GOLD
!===============================================================================

iphas = irapha(iph)

do ivar = 1,nbrayf

  if (irayvf(ivar,iphas).eq.1) then

    if      (ivar.eq.itparp) then
      do ifac = 1,nfabor
        w1(ifac) = tparoi(ifac,iph)
      enddo
    else if (ivar.eq.iqincp) then
      do ifac = 1,nfabor
        w1(ifac) = qincid(ifac,iph)
      enddo
    else if (ivar.eq.ixlamp) then
      do ifac = 1,nfabor
        w1(ifac) = xlam(ifac,iph)
      enddo
    else if (ivar.eq.iepap) then
      do ifac = 1,nfabor
        w1(ifac) = epa(ifac,iph)
      enddo
    else if (ivar.eq.iepsp) then
      do ifac = 1,nfabor
        w1(ifac) = eps(ifac,iph)
      enddo
    else if (ivar.eq.ifnetp) then
      do ifac = 1,nfabor
        w1(ifac) = flunet(ifac,iph)
      enddo
    else if (ivar.eq.ifconp) then
      do ifac = 1,nfabor
        w1(ifac) = flconv(ifac,iph)
      enddo
    else if (ivar.eq.ihconp) then
      do ifac = 1,nfabor
        w1(ifac) = hfconv(ifac,iph)
      enddo
    endif

    FICH = ' '
    fich = entete
    call verlon (fich,ii1,ii2,lpos)
    call verlon (nbrvaf(ivar,iphas),n1,n2,lpos)
    fich(ii2+1:ii2+lpos) = nbrvaf(ivar,iphas)(n1:n2)
    call verlon (fich,ii1,ii2,lpos)

!--> Ouverture

    open (impray,file=fich(ii1:ii2),                              &
          STATUS='UNKNOWN',FORM='FORMATTED',                      &
          ACCESS='SEQUENTIAL')

!--> Entete

    WRITE(IMPRAY,'(A)')  NBRVAF(IVAR,IPHAS)(N1:N2)

    ipart = 0

    do iz = 1,nbzone

      izone = izones(iz)

      nquad4 = 0
      ntria3 = 0
      nsided = 0

      do ifac = 1,nfabor
        if (izfrad(ifac,iph).eq.izone) then
          nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
          if (nbp.eq.4) then
            nquad4 = nquad4 + 1
          else if (nbp.eq.3) then
            ntria3 = ntria3 + 1
          else if (nbp.ge.5) then
            nsided = nsided + 1
          endif
        endif
      enddo


!--> Variable scalaire


      if (ntria3.gt.0) then
        ipart = ipart +1
        WRITE(IMPRAY,'(A)')    'part'
        WRITE(IMPRAY,'(I10)')  IPART
        WRITE(IMPRAY,'(A)')  'tria3'
        do ifac = 1,nfabor
          if (izfrad(ifac,iph).eq.izone) then
            nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
            IF(NBP.EQ.3) WRITE(IMPRAY,'(E12.5)') W1(IFAC)
          endif
        enddo
      endif

      if (nquad4.gt.0) then
        ipart = ipart +1
        WRITE(IMPRAY,'(A)')    'part'
        WRITE(IMPRAY,'(I10)')  IPART
        WRITE(IMPRAY,'(A)')  'quad4'
        do ifac = 1,nfabor
          if (izfrad(ifac,iph).eq.izone) then
            nbp = ipnfbr(ifac+1) - ipnfbr(ifac)
            IF(NBP.EQ.4) WRITE(IMPRAY,'(E12.5)') W1(IFAC)
          endif
        enddo
      endif

      if (nsided.gt.0) ipart = ipart +1

! Fin boucle sur les zones
    enddo

    close(impray)

  endif

enddo

!===============================================================================
! 5.  ECRITURE FICHIER .CASE AU FORMAT ENSIGHT GOLD
!===============================================================================

OPEN (IMPRAY,FILE = ENTETE // 'CASE',                             &
      STATUS='UNKNOWN',FORM='FORMATTED',                          &
      ACCESS='SEQUENTIAL')

WRITE(IMPRAY,'(A)') 'FORMAT'
WRITE(IMPRAY,'(A)') 'type:     ensight gold'
WRITE(IMPRAY,'(A)') 'GEOMETRY'
WRITE(IMPRAY,'(A)') 'model:    '//ENTETE//'geom'
WRITE(IMPRAY,'(A)') 'VARIABLE'

iphas = irapha(iph)

do ivar = 1,nbrayf

  if (irayvf(ivar,iphas).eq.1) then

    FICH = ' '
    FICH = 'scalar per element:'
    call verlon (fich,ii1,ii2,lpos)
    call verlon (nbrvaf(ivar,iphas),n1,n2,lpos)
    fich(ii2+11:ii2+11+lpos) = nbrvaf(ivar,iphas)(n1:n2)

    NAME = ' '
    name = entete
    call verlon (fich,ii1,ii2,lpos)
    call verlon (name,n1,n2,lpos)
    fich(ii2+2:ii2+1+lpos) = name(n1:n2)

    call verlon (fich,ii1,ii2,lpos)
    call verlon (nbrvaf(ivar,iphas),n1,n2,lpos)
    fich(ii2+1:ii2+lpos) = nbrvaf(ivar,iphas)(n1:n2)

    call verlon (fich,ii1,ii2,lpos)

    WRITE(IMPRAY,'(A)') FICH(II1:II2)

   endif

 enddo

 close(impray)

!--------
! FORMATS
!--------

 9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DE L''ECRITURE DES SORTIES ENSIGHT    ',/,&
'@    =========   SUR LE MAILLAGE FRONTIERE DANS LE MODULE    ',/,&
'@                DE TRANSFERT THERMIQUE RADIATIF             ',/,&
'@                                                            ',/,&
'@    SOUS DIMENSIONNEMENT DE NOZRDM                          ',/,&
'@                                                            ',/,&
'@    Il faut augmenter NOZRDM dans l''include radiat.h       ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit, mais les sorties ensight         ',/,&
'@    ne sont pas effectuees                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DE L''ECRITURE DES SORTIES ENSIGHT    ',/,&
'@    =========   SUR LE MAILLAGE FRONTIERE DANS LE MODULE    ',/,&
'@                DE TRANSFERT THERMIQUE RADIATIF             ',/,&
'@                                                            ',/,&
'@    PAS DE ZONES DE FRONTIERES DETECTEES                    ',/,&
'@                                                            ',/,&
'@    Il faut verifier le contenu du tableau IZFRAD           ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit, mais les sorties ensight         ',/,&
'@    ne sont pas effectuees                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end
