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

subroutine typecl &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   itypfb , itrifb , icodcl , isostd ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl , frcxt  ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   rdevel , rtuser , ra     )

!===============================================================================
! Function :
! --------

! Handle boundary condition type code (itypfb)

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
! lndfac           ! e  ! <-- ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr (optionnel          !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
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
! icodcl           ! te ! <-- ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! isostd           ! te ! --> ! indicateur de sortie standard                  !
!    (nfabor+1)    !    !     !  +numero de la face de reference               !
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
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! frcxt(ncelet,    ! tr ! <-- ! force exterieure generant la pression          !
!   3,nphas)       !    !     !  hydrostatique                                 !
! w1,2,3,4,5,6     ! tr ! --- ! tableaux de travail                            !
!  (ncelet         !    !     !  (calcul du gradient de pression)              !
! rijipb           ! tr ! --- ! tab de trav pour valeurs en iprime             !
! (nfabor,6   )    !    !     !  des rij au bord                               !
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
! Fortran common blocks
!===============================================================================

include "dimfbr.h"
include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstnum.h"
include "cstphy.h"
include "entsor.h"
include "pointe.h"
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itypfb(nfabor,nphas) , itrifb(nfabor,nphas)
integer          isostd(nfabor+1,nphas)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision rcodcl(nfabor,nvar,3)
double precision frcxt(ncelet,3,nphas)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,3)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

character        chaine*80
integer          idebia, idebra
integer          ifac, ivar, iel
integer          iok, inc, iccocg, ideb, ifin, inb, isum, iwrnp
integer          ifrslb(nphsmx), itbslb(nphsmx)
integer          ityp, ii, jj, iphas, iwaru, iflmab
integer          nswrgp, imligp, iwarnp
integer          ipriph, iuiph, iviph, iwiph
integer          ir11ip, ir22ip, ir33ip, ir12ip, ir13ip,ir23ip
integer          ikiph , iepiph, iphiph, ifbiph, iomgip
integer          iprnew, kphas, iii
integer          irangd, iclipr, iiptot
integer          itrois, ifadir, nfadir
double precision pref, epsrgp, climgp, extrap, coefup
double precision diipbx, diipby, diipbz
double precision flumbf, flumty(ntypmx)
double precision ro0iph, p0iph, pr0iph, xxp0, xyp0, xzp0
double precision xyzref(3)

integer          ipass
data             ipass /0/
save             ipass


!===============================================================================

!===============================================================================
! 1.  Initialization
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2.  Check consistency of types given in usclim
!===============================================================================

iok = 0

do iphas = 1, nphas
  do ifac = 1, nfabor
    ityp = itypfb(ifac,iphas)
    if(ityp.le.0.or.ityp.gt.ntypmx) then
      itypfb(ifac,iphas) = 0
      iok = iok + 1
    endif
  enddo
enddo

if (irangp.ge.0) call parcmx(iok)
if(iok.ne.0) then
  call bcderr(nphas, itypfb)
endif

!===============================================================================
! 3.  Sort boundary faces
!===============================================================================


! Count faces of each type (temporarily in ifinty)

do iphas = 1, nphas
  do ii = 1, ntypmx
   ifinty(ii,iphas) = 0
 enddo
enddo

do iphas = 1, nphas
  do ifac = 1, nfabor
    ityp = itypfb(ifac,iphas)
    ifinty(ityp,iphas) = ifinty(ityp,iphas) + 1
  enddo
enddo


! Set start of each group of faces in itrifb (sorted by type): idebty

do iphas = 1, nphas
  do ii = 1, ntypmx
    idebty(ii,iphas) = 1
  enddo
enddo

do iphas = 1, nphas
  do ii = 1, ntypmx-1
    do jj = ii+1, ntypmx
      idebty(jj,iphas) = idebty(jj,iphas) + ifinty(ii,iphas)
    enddo
  enddo
enddo

! Sort faces in itrifb and use the opportunity to correctly set ifinty

do iphas = 1, nphas
  do ii = 1, ntypmx
    ifinty(ii,iphas) = idebty(ii,iphas)-1
  enddo
enddo

do iphas = 1, nphas
  do ifac = 1, nfabor
    ityp = itypfb(ifac,iphas)
    ifin = ifinty(ityp,iphas)+1
    itrifb(ifin,iphas) = ifac
    ifinty(ityp,iphas) = ifin
  enddo
enddo

! Basic check

iok = 0
do iphas = 1, nphas
  do ii = 1, ntypmx-1
    if(ifinty(ii,iphas).ge.idebty(ii+1,iphas)) then
      if (iok.eq.0) iok = ii
    endif
  enddo
  if (irangp.ge.0) call parcmx(iok)
  if (iok.gt.0) then
    ii = iok
    write(nfecra,2010) iphas
    if(ifinty(ii,iphas).ge.idebty(ii+1,iphas)) then
      write(nfecra,2020) (ifinty(jj,iphas),jj=1,ntypmx)
      write(nfecra,2030) (idebty(jj,iphas),jj=1,ntypmx)
      write(nfecra,2040) (itypfb(jj,iphas),jj=1,nfabor)
      write(nfecra,2098) ii,ifinty(ii,iphas),ii+1,idebty(ii+1,iphas)
    else
      write(nfecra,2099) ii,ii+1
    endif
    call csexit (1)
  endif
enddo

iok = 0
do iphas = 1, nphas
  isum = 0
  do ii = 1, ntypmx
    isum = isum + ifinty(ii,iphas) - idebty(ii,iphas) + 1
  enddo
  if (irangp.ge.0) call parcpt (isum)
  if(isum.ne.nfbrgb) then
    write(nfecra,3099) iphas, isum, nfbrgb
    iok = iok + 1
  endif
enddo
if (iok.ne.0) then
  call csexit (1)
  !==========
endif


! ---> On ecrit les types de faces avec la borne inf et sup et le nb
!       pour chaque type de face trouve (tjrs pour les types par defaut)

if(ipass.eq.0.or.iwarni(iu(1)).ge.2) then

  ipass = 1

  write(nfecra,6010)

  do iphas = 1, nphas

    write(nfecra,6011)iphas

    if ( ippmod(icompf).lt.0 ) then

      ii = ientre
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Entree           ', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)
      ii = iparoi
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Paroi lisse      ', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)
      ii = iparug
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Paroi rugueuse   ', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)
      ii = isymet
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Symetrie         ', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

       ii = isolib
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Sortie libre     ', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      ii = iindef
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Indefini         ', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      do ii = 1, ntypmx
        if( ii.ne.ientre .and.                                    &
             ii.ne.iparoi .and.                                   &
             ii.ne.iparug .and.                                   &
             ii.ne.isymet .and.                                   &
             ii.ne.isolib .and.                                   &
             ii.ne.iindef ) then
          inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
          if (irangp.ge.0) call parcpt (inb)
          if(inb.gt.0) then
            WRITE(NFECRA,6020) 'Type utilisateur ', ii,           &
                 inb,                                             &
                 ifinty(ii,iphas)-idebty(ii,iphas)+1,             &
                 idebty(ii,iphas),ifinty(ii,iphas)
          endif
        endif
      enddo

    else

      ii = ieqhcf
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Entree sub. enth.', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      ii = ierucf
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Entree subsonique', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      ii = iesicf
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Entree/Sortie imp', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      ii = isopcf
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Sortie subsonique', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      ii = isspcf
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Sortie supersoniq', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      ii = iparoi
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Paroi lisse      ', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      ii = iparug
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Paroi rugueuse   ', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      ii = isymet
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Symetrie         ', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      ii = iindef
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) call parcpt (inb)
      WRITE(NFECRA,6020) 'Indefini         ', ii,                 &
           inb,                                                   &
           ifinty(ii,iphas)-idebty(ii,iphas)+1,                   &
           idebty(ii,iphas),ifinty(ii,iphas)

      do ii = 1, ntypmx
        if (ii.ne.iesicf .and. &
            ii.ne.isspcf .and. &
            ii.ne.ieqhcf .and. &
            ii.ne.ierucf .and. &
            ii.ne.isopcf .and. &
            ii.ne.iparoi .and. &
            ii.ne.iparug .and. &
            ii.ne.isymet .and. &
            ii.ne.iindef ) then
          inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
          if (irangp.ge.0) call parcpt (inb)
          if(inb.gt.0) then
            WRITE(NFECRA,6020) 'Type utilisateur ',II,            &
                 inb,                                             &
                 ifinty(ii,iphas)-idebty(ii,iphas)+1,             &
                 idebty(ii,iphas),ifinty(ii,iphas)
          endif
        endif
      enddo

    endif

    write(nfecra,6030)

  enddo

endif

!================================================================================
! 4.  rcodcl(.,    .,1) has been initialized as rinfin so as to check what
!     the user has modified. Those not modified are reset to zero here.
!     isolib and ientre are handled later.
!================================================================================

do iphas = 1, nphas
  do ivar=1, nvar
     do ifac = 1, nfabor
        if((itypfb(ifac,iphas) .ne. isolib) .and. &
           (itypfb(ifac,iphas) .ne. ientre) .and. &
           (rcodcl(ifac,ivar,1) .gt. rinfin*0.5d0)) then
           rcodcl(ifac,ivar,1) = 0.d0
       endif
      enddo
   enddo
enddo


!===============================================================================
! 5.  Compute pressure at boundary (in coefu(*,1))
!     (if we need it, that is if there are outlet boudary faces).

!     The loop on phases starts here and ends at the end of the next block.
!===============================================================================

! --- Boucle sur les phases : debut
do iphas = 1, nphas

  ro0iph = ro0  (iphas)
  p0iph  = p0   (iphas)
  pr0iph = pred0(iphas)
  xxp0   = xyzp0(1,iphas)
  xyp0   = xyzp0(2,iphas)
  xzp0   = xyzp0(3,iphas)
  ipriph = ipr (iphas)
  iuiph  = iu  (iphas)
  iviph  = iv  (iphas)
  iwiph  = iw  (iphas)
  if(itytur(iphas).eq.2) then
    ikiph  = ik(iphas)
    iepiph = iep(iphas)
  elseif(itytur(iphas).eq.3) then
    ir11ip = ir11(iphas)
    ir22ip = ir22(iphas)
    ir33ip = ir33(iphas)
    ir12ip = ir12(iphas)
    ir13ip = ir13(iphas)
    ir23ip = ir23(iphas)
    iepiph = iep(iphas)
  elseif(iturb(iphas).eq.50) then
    ikiph  = ik(iphas)
    iepiph = iep(iphas)
    iphiph = iphi(iphas)
    ifbiph = ifb(iphas)
  elseif(iturb(iphas).eq.60) then
    ikiph  = ik (iphas)
    iomgip = iomg(iphas)
  endif

! Check if the pressure (unique) has not been handled already

  iprnew = 1
  if(iphas.gt.1) then
    do kphas = 1, iphas-1
      if(ipr(iphas).eq.ipr(kphas)) then
        iprnew = 0
      endif
    enddo
  endif

! ifrslb = first free standard outlet face (ICODCL not modified)
! itbslb = max of ifrslb on all ranks, standard outlet face presence inidcator

  ifrslb(iphas) = 0
  do ii = ifinty(isolib,iphas), idebty(isolib,iphas), -1
    ifac = itrifb(ii,iphas)
    if (icodcl(ifac,ipriph).eq.0)                                 &
         ifrslb(iphas) = ifac
  enddo
  itbslb(iphas) = ifrslb(iphas)
  if (irangp.ge.0) then
     call parcmx (itbslb(iphas))
  endif

  if ((itbslb(iphas).gt.0) .and. (iprnew.eq.1)) then

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ipriph)
    imligp = imligr(ipriph)
    iwarnp = iwarni(ipriph)
    epsrgp = epsrgr(ipriph)
    climgp = climgr(ipriph)
    extrap = extrag(ipriph)
    iclipr = iclrtp(ipriph,icoef)

    call grdcel                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ipriph , imrgra , inc    , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   frcxt(1,1,iphas), frcxt(1,2,iphas), frcxt(1,3,iphas),          &
   rtpa(1,ipriph)  , coefa(1,iclipr) , coefb(1,iclipr) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )


!  Put in coefu the value at I' or F (depending on iphydr) of the
!  total pressure, computed from P*

    if (iphydr.eq.0) then
      do ifac = 1, nfabor
        ii = ifabor(ifac)
        iii = idiipb-1+3*(ifac-1)
        diipbx = ra(iii+1)
        diipby = ra(iii+2)
        diipbz = ra(iii+3)
        coefu(ifac,1) = rtpa(ii,ipriph)                           &
             + diipbx*w1(ii)+ diipby*w2(ii) + diipbz*w3(ii)       &
             + ro0iph*( gx*(cdgfbo(1,ifac)-xxp0)                  &
                      + gy*(cdgfbo(2,ifac)-xyp0)                  &
                      + gz*(cdgfbo(3,ifac)-xzp0))                 &
             + p0iph - pr0iph
      enddo
    else
      do ifac = 1, nfabor
        ii = ifabor(ifac)
        coefu(ifac,1) = rtpa(ii,ipriph)                           &
             + (cdgfbo(1,ifac)-xyzcen(1,ii))*w1(ii)               &
             + (cdgfbo(2,ifac)-xyzcen(2,ii))*w2(ii)               &
             + (cdgfbo(3,ifac)-xyzcen(3,ii))*w3(ii)               &
             + ro0iph*(  gx*(cdgfbo(1,ifac)-xxp0)                 &
                       + gy*(cdgfbo(2,ifac)-xyp0)                 &
                       + gz*(cdgfbo(3,ifac)-xzp0))                &
             + p0iph - pr0iph
      enddo
    endif

  endif


!===============================================================================
! 6.  Convert to rcodcl and icodcl
!     (if this has not already been set by the user)

!     First, process variables for which a specific treatement is done
!     (pressure, velocity, ...)
!===============================================================================

! 6.1 ENTREE
! ===========

! ---> La pression a un traitement Neumann, le reste Dirichlet
!                                           sera traite plus tard.

  ideb = idebty(ientre,iphas)
  ifin = ifinty(ientre,iphas)

  do ivar = 1, nvar
    if (ivar.eq.ipriph) then
      if(iprnew.eq.1) then
        do ii = ideb, ifin
          ifac = itrifb(ii,iphas)
          if(icodcl(ifac,ivar).eq.0) then
            icodcl(ifac,ivar)   = 3
            rcodcl(ifac,ivar,1) = 0.d0
            rcodcl(ifac,ivar,2) = rinfin
            rcodcl(ifac,ivar,3) = 0.d0
          endif
        enddo
      endif
    endif
  enddo


! 6.2 SORTIE (entrée-sortie libre) (ISOLIB)
! ===================

! ---> La pression a un traitement Dirichlet, les vitesses 9
!        (le reste Neumann, ou Dirichlet si donnée utilisateur,
!        sera traite plus tard)

  if (iphydr.eq.1) then

!     En cas de prise en compte de la pression hydrostatique,
!     on remplit le tableau ISOSTD
!     0 -> pas une face de sortie standard (i.e. pas sortie ou sortie avec CL
!                                                de pression modifiee)
!     1 -> face de sortie libre avec CL de pression automatique.
!     le numero de la face de reference est stocke dans ISOSTD(NFABOR+1,IPHAS)
!     qui est d'abord initialise a -1 (i.e. pas de face de sortie std)
    isostd(nfabor+1,iphas) = -1
    do ifac = 1,nfabor
      isostd(ifac,iphas) = 0
      if ((itypfb(ifac,iphas).eq.isolib).and.                     &
           (icodcl(ifac,ipriph).eq.0)) then
        isostd(ifac,iphas) = 1
      endif
    enddo
  endif

! ---> Pression de recalage (unique, meme s'il y a plusieurs sorties)
!     En cas de prise en compte de la pression hydrostatique, on cherche
!     la face de reference

!   Determination de la pression P I' unique en parallele
!     S'il y a des faces de sortie libre, on cherche le premier
!     proc sur lequel il y en a.

!     On recupere aussi les coordonnees du point de reference pour les
!     sorties, pour calculer PREF ensuite

  if    (itbslb(iphas).gt.0) then

    if (irangp.ge.0) then

!     Indicateur de sorties locales et pointeur sur la premiere face
!       S'il y a des sorties libres standards quelque part
!         et s'il n'y en a pas en local, on affecte a IRANGD la valeur -1
      if(ifrslb(iphas).le.0) then
        irangd = -1
!         et s'il y en a en local, on affecte a IRANGD la valeur IRANGP
      else
        irangd = irangp
      endif

!     Valeur de P I'
!         on prend le numero du dernier proc sur lequel il y en a
      call parcmx(irangd)
!         si c'est le proc courant, on affecte P I' a COEFUP
      if (irangp.eq.irangd) then
        coefup = coefu(ifrslb(iphas),1)
        if (iphydr.eq.1) isostd(nfabor+1,iphas) = ifrslb(iphas)
!         sinon on affecte 0 a COEFUP
      else
        coefup = 0.d0
      endif
!         la somme sur les procs de COEFUP donne donc P I', disponible
!            pour tous les procs (c'est plus simple qu'un bcast pour
!            lequel il faudrait au prealable que tous les proc sachent
!            quel proc envoie).
      call parsom(coefup)

!     Reference de pression pour les sorties
      if (irangp.eq.irangd) then
        xyzref(1) = cdgfbo(1,ifrslb(iphas))
        xyzref(2) = cdgfbo(2,ifrslb(iphas))
        xyzref(3) = cdgfbo(3,ifrslb(iphas))
      else
        xyzref(1) = 0.d0
        xyzref(2) = 0.d0
        xyzref(3) = 0.d0
      endif
      itrois = 3
      call parrsm(itrois,xyzref)

!   Determination de la pression P I' unique en sequentiel
!     on repere la premiere face de sortie libre standard
!                       et on affecte PI' a COEFUP

    else

      coefup = coefu(ifrslb(iphas),1)
      if (iphydr.eq.1) isostd(nfabor+1,iphas) = ifrslb(iphas)
      xyzref(1) = cdgfbo(1,ifrslb(iphas))
      xyzref(2) = cdgfbo(2,ifrslb(iphas))
      xyzref(3) = cdgfbo(3,ifrslb(iphas))
    endif

!     Si l'utilisateur n'a rien specifie, on met IXYZP0 a 2 pour mettre
!     a jour le point de reference
    if (ixyzp0(iphas).eq.-1) ixyzp0(iphas) = 2

  elseif (ixyzp0(iphas).lt.0) then
!     S'il n'y a pas de faces de sortie, on cherche des Dirichlets
!     eventuels specifies par l'utilisateur pour y localiser
!     le point de reference.
    ifadir = -1
    irangd = -1
    do ifac = 1, nfabor
      if (icodcl(ifac,ipriph).eq.1) then
        ifadir = ifac
        irangd = irangp
      endif
    enddo
    nfadir = ifadir
    if (irangp.ge.0) call parcmx(nfadir)
    if (nfadir.gt.0) then
!     on met IXYZP0 a 2 pour mettre a jour le point de reference
      ixyzp0(iphas) = 2
!     en parallele on prend le numero du dernier proc sur lequel il y en a
      if (irangp.ge.0) then
        call parcmx(irangd)
        if (irangp.eq.irangd) then
          xyzref(1) = cdgfbo(1,ifadir)
          xyzref(2) = cdgfbo(2,ifadir)
          xyzref(3) = cdgfbo(3,ifadir)
        else
          xyzref(1) = 0.d0
          xyzref(2) = 0.d0
          xyzref(3) = 0.d0
        endif
        itrois = 3
        call parrsm(itrois,xyzref)
      else
        xyzref(1) = cdgfbo(1,ifadir)
        xyzref(2) = cdgfbo(2,ifadir)
        xyzref(3) = cdgfbo(3,ifadir)
      endif
    endif

  endif


!   Si le point de reference n'a pas ete specifie par l'utilisateur
!   on le change et on decale alors COEFU s'il y a des sorties.
!   La pression totale dans PROPCE est aussi decalee (c'est a priori
!   inutile sauf si l'utilisateur l'utilise dans ustsns par exemple)

  if (ixyzp0(iphas).eq.2) then
    ixyzp0(iphas) = 1
    iiptot = ipproc(iprtot(iphas))
    xxp0 = xyzref(1) - xyzp0(1,iphas)
    xyp0 = xyzref(2) - xyzp0(2,iphas)
    xzp0 = xyzref(3) - xyzp0(3,iphas)
    xyzp0(1,iphas) = xyzref(1)
    xyzp0(2,iphas) = xyzref(2)
    xyzp0(3,iphas) = xyzref(3)
    do iel = 1, ncelet
      propce(iel,iiptot) = propce(iel,iiptot)                     &
           - ro0iph*( gx*xxp0 + gy*xyp0 + gz*xzp0 )
    enddo
    if (itbslb(iphas).gt.0) then
      write(nfecra,8000)iphas,xxp0,xyp0,xzp0
      do ifac = 1, nfabor
        coefu(ifac,1) = coefu(ifac,1)                             &
             - ro0iph*( gx*xxp0 + gy*xyp0 + gz*xzp0 )
      enddo
      coefup = coefup - ro0iph*( gx*xxp0 + gy*xyp0 + gz*xzp0 )
    else
      write(nfecra,8001)iphas,xxp0,xyp0,xzp0
    endif
  elseif (ixyzp0(iphas).eq.-1) then
!     Il n'y a pas de sorties ni de Dirichlet et l'utilisateur n'a
!     rien specifie -> on met IXYZP0 a 0 pour ne plus y toucher, tout
!     en differenciant du cas =1 qui necessitera une ecriture en suite
    ixyzp0(iphas) = 0
  endif

!     La pression totale doit etre recalee en Xref a la valeur
!     Po + rho_0*g.(Xref-X0)
  if (itbslb(iphas).gt.0) then
    xxp0 = xyzp0(1,iphas)
    xyp0 = xyzp0(2,iphas)
    xzp0 = xyzp0(3,iphas)
    pref = p0(iphas)                                              &
         + ro0iph*( gx*(xyzref(1)-xxp0)                           &
                  + gy*(xyzref(2)-xyp0)                           &
                  + gz*(xyzref(3)-xzp0) )                         &
         - coefup
  endif


! ---> Entree/Sortie libre

  ideb = idebty(isolib,iphas)
  ifin = ifinty(isolib,iphas)

  do ivar = 1, nvar
    if (ivar.eq.ipriph) then
      if(iprnew.eq.1) then
        do ii = ideb, ifin
          ifac = itrifb(ii,iphas)
          if(icodcl(ifac,ivar).eq.0) then
            icodcl(ifac,ivar)   = 1
            rcodcl(ifac,ivar,1) = coefu(ifac,1) + pref
            rcodcl(ifac,ivar,2) = rinfin
            rcodcl(ifac,ivar,3) = 0.d0
          endif
        enddo
      endif
    elseif(ivar.eq.iuiph.or.ivar.eq.iviph.or.ivar.eq.iwiph) then
      do ii = ideb, ifin
        ifac = itrifb(ii,iphas)
        if(icodcl(ifac,ivar).eq.0) then
          icodcl(ifac,ivar)   = 9
          rcodcl(ifac,ivar,1) = 0.d0
          rcodcl(ifac,ivar,2) = rinfin
          rcodcl(ifac,ivar,3) = 0.d0
        endif
      enddo
    endif
  enddo


! 6.3 SYMETRIE
! =============

! ---> Les vecteurs et tenseurs ont un traitement particulier
!        le reste Neumann sera traite plus tard

  ideb = idebty(isymet,iphas)
  ifin = ifinty(isymet,iphas)

  do ivar = 1, nvar
    if ( ivar.eq.iuiph.or.ivar.eq.iviph.or.ivar.eq.iwiph.or.      &
       ( itytur(iphas).eq.3.and.                                  &
          (ivar.eq.ir11ip.or.ivar.eq.ir22ip.or.ivar.eq.ir33ip.or. &
           ivar.eq.ir12ip.or.ivar.eq.ir13ip.or.ivar.eq.ir23ip)    &
                                                          ) ) then
      do ii = ideb, ifin
        ifac = itrifb(ii,iphas)
        if(icodcl(ifac,ivar).eq.0) then
          icodcl(ifac,ivar)   = 4
!                RCODCL(IFAC,IVAR,1) = Modifie eventuellement par l'ALE
          rcodcl(ifac,ivar,2) = rinfin
          rcodcl(ifac,ivar,3) = 0.d0
        endif
      enddo
    elseif(ivar.eq.ipriph) then
      if(iprnew.eq.1) then
        do ii = ideb, ifin
          ifac = itrifb(ii,iphas)
          if(icodcl(ifac,ivar).eq.0) then
            icodcl(ifac,ivar)   = 3
            rcodcl(ifac,ivar,1) = 0.d0
            rcodcl(ifac,ivar,2) = rinfin
            rcodcl(ifac,ivar,3) = 0.d0
          endif
        enddo
      endif
    endif
  enddo

! 6.4 PAROI LISSE
! ===============

! ---> La vitesse et les grandeurs turbulentes ont le code 5
!        le reste Neumann sera traite plus tard

  ideb = idebty(iparoi,iphas)
  ifin = ifinty(iparoi,iphas)

  do ivar = 1, nvar
    if ( ivar.eq.iuiph.or.ivar.eq.iviph.or.ivar.eq.iwiph) then
      do ii = ideb, ifin
        ifac = itrifb(ii,iphas)
        if(icodcl(ifac,ivar).eq.0) then
          icodcl(ifac,ivar)   = 5
!               RCODCL(IFAC,IVAR,1) = Utilisateur
          rcodcl(ifac,ivar,2) = rinfin
          rcodcl(ifac,ivar,3) = 0.d0
        endif
      enddo
    elseif (                                                      &
       ( itytur(iphas).eq.2.and.                                  &
          (ivar.eq.ikiph  .or.ivar.eq.iepiph) ).or.               &
       ( itytur(iphas).eq.3.and.                                  &
          (ivar.eq.ir11ip.or.ivar.eq.ir22ip.or.ivar.eq.ir33ip.or. &
           ivar.eq.ir12ip.or.ivar.eq.ir13ip.or.ivar.eq.ir23ip.or. &
           ivar.eq.iepiph)                    ).or.               &
       ( iturb(iphas).eq.50.and.                                  &
          (ivar.eq.ikiph.or.ivar.eq.iepiph.or.ivar.eq.iphiph.or.  &
           ivar.eq.ifbiph)                    ).or.               &
       ( iturb(iphas).eq.60.and.                                  &
          (ivar.eq.ikiph.or.ivar.eq.iomgip)   ) ) then
      do ii = ideb, ifin
        ifac = itrifb(ii,iphas)
        if(icodcl(ifac,ivar).eq.0) then
          icodcl(ifac,ivar)   = 5
          rcodcl(ifac,ivar,1) = 0.d0
          rcodcl(ifac,ivar,2) = rinfin
          rcodcl(ifac,ivar,3) = 0.d0
        endif
      enddo
    elseif(ivar.eq.ipriph) then
      if(iprnew.eq.1) then
        do ii = ideb, ifin
          ifac = itrifb(ii,iphas)
          if(icodcl(ifac,ivar).eq.0) then
            icodcl(ifac,ivar)   = 3
            rcodcl(ifac,ivar,1) = 0.d0
            rcodcl(ifac,ivar,2) = rinfin
            rcodcl(ifac,ivar,3) = 0.d0
          endif
        enddo
      endif
    endif
  enddo

! 6.5 PAROI RUGUEUSE
! ==================

! ---> La vitesse et les grandeurs turbulentes ont le code 6
!      la rugosite est stockee dans rcodcl(..,..,3)
!      le reste Neumann sera traite plus tard (idem paroi lisse)

  ideb = idebty(iparug,iphas)
  ifin = ifinty(iparug,iphas)

  do ivar = 1, nvar
    if ( ivar.eq.iuiph.or.ivar.eq.iviph.or.ivar.eq.iwiph) then
      do ii = ideb, ifin
        ifac = itrifb(ii,iphas)
        if(icodcl(ifac,ivar).eq.0) then
          icodcl(ifac,ivar)   = 6
!               RCODCL(IFAC,IVAR,1) = Utilisateur
          rcodcl(ifac,ivar,2) = rinfin
!               RCODCL(IFAC,IVAR,3) = Utilisateur
        endif
      enddo
    elseif (                                                      &
       ( itytur(iphas).eq.2.and.                                  &
          (ivar.eq.ikiph  .or.ivar.eq.iepiph) ).or.               &
       ( itytur(iphas).eq.3.and.                                  &
          (ivar.eq.ir11ip.or.ivar.eq.ir22ip.or.ivar.eq.ir33ip.or. &
           ivar.eq.ir12ip.or.ivar.eq.ir13ip.or.ivar.eq.ir23ip.or. &
           ivar.eq.iepiph)                    ).or.               &
       ( iturb(iphas).eq.50.and.                                  &
          (ivar.eq.ikiph.or.ivar.eq.iepiph.or.ivar.eq.iphiph.or.  &
           ivar.eq.ifbiph)                    ).or.               &
       ( iturb(iphas).eq.60.and.                                  &
          (ivar.eq.ikiph.or.ivar.eq.iomgip)   ) ) then
      do ii = ideb, ifin
        ifac = itrifb(ii,iphas)
        if(icodcl(ifac,ivar).eq.0) then
          icodcl(ifac,ivar)   = 6
          rcodcl(ifac,ivar,1) = 0.d0
          rcodcl(ifac,ivar,2) = rinfin
          rcodcl(ifac,ivar,3) = 0.d0
        endif
      enddo
    elseif(ivar.eq.ipriph) then
      if(iprnew.eq.1) then
        do ii = ideb, ifin
          ifac = itrifb(ii,iphas)
          if(icodcl(ifac,ivar).eq.0) then
            icodcl(ifac,ivar)   = 3
            rcodcl(ifac,ivar,1) = 0.d0
            rcodcl(ifac,ivar,2) = rinfin
            rcodcl(ifac,ivar,3) = 0.d0
          endif
        enddo
      endif
    endif
  enddo


enddo
! --- Boucle sur les phases : fin


!===============================================================================
! 6.bis  CONVERSION EN RCODCL ICODCL
!   (SI CE DERNIER N'A PAS DEJA ETE RENSEIGNE PAR L'UTILISATEUR)

!     MAINTENANT LES VARIABLES POUR LESQUELLES IL N'EXISTE PAS DE
!    TRAITEMENT PARTICULIER (HORS PRESSION, VITESSE ...)
!===============================================================================

! --- Boucle sur les phases : debut
do iphas = 1, nphas


! 6.1 ENTREE bis
! ===========

! ---> La pression a un traitement Neumann (deja traitee plus haut),
!      La vitesse  Dirichlet. Les scalaires ont un traitement
!     Dirichlet si l'utilisateur fournit une valeur, sinon on utilise
!     Neumann homogene si le flux de masse est sortant (erreur sinon).

  ideb = idebty(ientre,iphas)
  ifin = ifinty(ientre,iphas)

  iok = 0
  do ivar = 1, nvar
    do ii = ideb, ifin
      ifac = itrifb(ii,iphas)
      if(icodcl(ifac,ivar).eq.0) then

        if (ivar.eq.iuiph.or.ivar.eq.iviph.or.ivar.eq.iwiph)      &
             then
          if (rcodcl(ifac,ivar,1).gt.rinfin*0.5d0) then
            itypfb(ifac,iphas) = - abs(itypfb(ifac,iphas))
            if (iok.eq.0) then
              iok = 1
              chaine = nomvar(ipprtp(ivar))
              write(nfecra,6040)ifac,chaine(1:8)
            endif
          else
            icodcl(ifac,ivar) = 1
!                 RCODCL(IFAC,IVAR,1) = Utilisateur
            rcodcl(ifac,ivar,2) = rinfin
            rcodcl(ifac,ivar,3) = 0.d0
          endif

        elseif (rcodcl(ifac,ivar,1).gt.rinfin*0.5d0) then

          flumbf = propfb(ifac,ipprob(ifluma(iuiph)))
          if( flumbf.ge.-epzero) then
            icodcl(ifac,ivar)   = 3
            rcodcl(ifac,ivar,1) = 0.d0
            rcodcl(ifac,ivar,2) = rinfin
            rcodcl(ifac,ivar,3) = 0.d0
          else
            itypfb(ifac,iphas) = - abs(itypfb(ifac,iphas))
            if (iok.lt.2) then
              iok = iok + 2
              chaine = nomvar(ipprtp(ivar))
              write(nfecra,6050)ifac,chaine(1:8)
            endif
          endif
        else
          icodcl(ifac,ivar) = 1
!               RCODCL(IFAC,IVAR,1) = Utilisateur
          rcodcl(ifac,ivar,2) = rinfin
          rcodcl(ifac,ivar,3) = 0.d0
        endif

      endif
    enddo
  enddo

  if (irangp.ge.0) call parcmx(iok)
  if (iok.gt.0) then
    if (iok.eq.1 .or. iok.eq.3) write(nfecra,6060)
    if (iok.eq.2 .or. iok.eq.3) write(nfecra,6070)
    call bcderr(nphas, itypfb)
  endif



! 6.2 SORTIE (entree sortie libre)
! ===================

! ---> La pression a un traitement Dirichlet, les vitesses 9 ont ete
!        traites plus haut.
!      Le reste Dirichlet si l'utilisateur fournit une donnee
!        (flux de masse entrant ou sortant).
!      S'il n'y a pas de donnee utilisateur, on utilise un Neumann homogene
!        (flux entrant et sortant)


! ---> Sortie ISOLIB

  ideb = idebty(isolib,iphas)
  ifin = ifinty(isolib,iphas)

  do ivar = 1, nvar
    do ii = ideb, ifin
      ifac = itrifb(ii,iphas)
      if(icodcl(ifac,ivar).eq.0) then

         if (rcodcl(ifac,ivar,1).gt.rinfin*0.5d0) then
              icodcl(ifac,ivar) = 3
              rcodcl(ifac,ivar,1) = 0.d0
              rcodcl(ifac,ivar,2) = rinfin
              rcodcl(ifac,ivar,3) = 0.d0
          else
              icodcl(ifac,ivar) = 1
!                   RCODCL(IFAC,IVAR,1) = Utilisateur
              rcodcl(ifac,ivar,2) = rinfin
              rcodcl(ifac,ivar,3) = 0.d0
         endif
      endif
    enddo
  enddo



! 6.3 SYMETRIE bis
! =============

! ---> Les vecteurs et tenseurs ont un traitement particulier
!        traite plus haut
!        le reste Neumann

  ideb = idebty(isymet,iphas)
  ifin = ifinty(isymet,iphas)

  do ivar = 1, nvar
    do ii = ideb, ifin
      ifac = itrifb(ii,iphas)
      if(icodcl(ifac,ivar).eq.0) then
        icodcl(ifac,ivar)   = 3
        rcodcl(ifac,ivar,1) = 0.d0
        rcodcl(ifac,ivar,2) = rinfin
        rcodcl(ifac,ivar,3) = 0.d0
      endif
    enddo
  enddo

! 6.4 PAROI LISSE bis
! ===============

! ---> La vitesse et les grandeurs turbulentes ont le code 5
!        traite plus haut
!        le reste Neumann

  ideb = idebty(iparoi,iphas)
  ifin = ifinty(iparoi,iphas)

  do ivar = 1, nvar
    do ii = ideb, ifin
      ifac = itrifb(ii,iphas)
      if(icodcl(ifac,ivar).eq.0) then
        icodcl(ifac,ivar)   = 3
        rcodcl(ifac,ivar,1) = 0.d0
        rcodcl(ifac,ivar,2) = rinfin
        rcodcl(ifac,ivar,3) = 0.d0
      endif
    enddo
  enddo

! 6.5 PAROI RUGUEUSE bis
! ==================

! ---> La vitesse et les grandeurs turbulentes ont le code 6
!        traite plus haut
!        le reste Neumann

  ideb = idebty(iparug,iphas)
  ifin = ifinty(iparug,iphas)

  do ivar = 1, nvar
    do ii = ideb, ifin
      ifac = itrifb(ii,iphas)
      if(icodcl(ifac,ivar).eq.0) then
        icodcl(ifac,ivar)   = 3
        rcodcl(ifac,ivar,1) = 0.d0
        rcodcl(ifac,ivar,2) = rinfin
        rcodcl(ifac,ivar,3) = 0.d0
      endif
    enddo
  enddo

enddo
! --- Boucle sur les phases : fin
!===============================================================================
! 7.  RENFORCEMENT DIAGONALE DE LA MATRICE SI AUCUN POINTS DIRICHLET
!===============================================================================
! On renforce si ISTAT=0 et si l'option est activee (IDIRCL=1)
! Si une de ces conditions est fausse, on force NDIRCL a valoir
! au moins 1 pour ne pas declaer la diagonale.

do ivar = 1, nvar
  ndircl(ivar) = 0
  if ( istat(ivar).gt.0 .or. idircl(ivar).eq.0 ) ndircl(ivar) = 1
enddo

do ivar = 1, nvar
  do ifac = 1, nfabor
    if( icodcl(ifac,ivar).eq.1 .or. icodcl(ifac,ivar).eq.5 ) then
      ndircl(ivar) = ndircl(ivar) +1
    endif
  enddo
  if (irangp.ge.0) call parcpt (ndircl(ivar))
enddo

!===============================================================================
! 8.  ON CALCULE LE FLUX DE MASSE AUX DIFFERENTS TYPES DE FACES
!       ET ON IMPRIME.

!     Ca serait utile de faire l'impression dans ECRLIS, mais attention,
!       on imprime le flux de masse du pas de temps precedent
!       or dans ECRLIS, on imprime a la fin du pas de temps
!       d'ou une petite incoherence possible.
!     D'autre part, ca serait utile de sortir d'autres grandeurs
!       (flux de chaleur par exemple, bilan de scalaires ...)

!===============================================================================

iwaru = -1
do iphas = 1, nphas
  iuiph  = iu(iphas)
  iwaru = max(iwarni(iuiph),iwaru)
enddo
if (irangp.ge.0) call parcmx(iwaru)

if(iwaru.ge.1 .or. mod(ntcabs,ntlist).eq.0                        &
       .or.(ntcabs.le.ntpabs+2).or.(ntcabs.ge.ntmabs-1)) then
  write(nfecra,7010)
endif

do iphas = 1, nphas

  iuiph  = iu(iphas)
  iflmab = ipprob(ifluma(iuiph))

  iwrnp = iwarni(iu(iphas))
  if (irangp.ge.0) call parcmx (iwrnp)
                   !==========

!     On ecrit le flux de masse si IWARNI>0, a la periodicite NTLIST
!     et au deux premiers et deux derniers pas de temps.
  if(iwrnp.ge.1 .or. mod(ntcabs,ntlist).eq.0                      &
       .or.(ntcabs.le.ntpabs+2).or.(ntcabs.ge.ntmabs-1)) then

    do ii = 1, ntypmx
      flumty(ii) = 0.d0
    enddo

    do ii = 1, ntypmx
      ideb = idebty(ii,iphas)
      ifin = ifinty(ii,iphas)
      do jj = ideb, ifin
        ifac = itrifb(jj,iphas)
        flumty(ii) = flumty(ii) + propfb(ifac,iflmab)
      enddo
    enddo


    write(nfecra,7011) iphas

    if (ippmod(icompf).lt.0 ) then

      ii = ientre
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Entree           ',II,INB,FLUMTY(II)
      ii = iparoi
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Paroi lisse      ',II,INB,FLUMTY(II)
      ii = iparug
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Paroi rugueuse   ',II,INB,FLUMTY(II)
      ii = isymet
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Symetrie         ',II,INB,FLUMTY(II)

      ii = isolib
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Sortie libre     ',II,INB,FLUMTY(II)


      ii = iindef
     inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Indefini         ',II,INB,FLUMTY(II)

      do ii = 1, ntypmx
        if( ii.ne.ientre .and.                                    &
            ii.ne.iparoi .and.                                    &
            ii.ne.iparug .and.                                    &
            ii.ne.isymet .and.                                    &
            ii.ne.isolib .and.                                    &
            ii.ne.iindef ) then
          inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
          if (irangp.ge.0) then
            call parcpt (inb)
            call parsom (flumty(ii))
          endif
          if(inb.gt.0) then
            write(nfecra,7020)                                    &
                 'Type utilisateur ',II,INB,FLUMTY(II)
          endif
        endif
      enddo

    else

      ii = ieqhcf
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Entree sub. enth.',II,INB,FLUMTY(II)

      ii = ierucf
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Entree subsonique',II,INB,FLUMTY(II)

      ii = iesicf
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Entree/Sortie imp',II,INB,FLUMTY(II)

      ii = isopcf
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Sortie subsonique',II,INB,FLUMTY(II)

      ii = isspcf
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Sortie supersoniq',II,INB,FLUMTY(II)

      ii = iparoi
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Paroi            ',II,INB,FLUMTY(II)

      ii = isymet
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Symetrie         ',II,INB,FLUMTY(II)

      ii = iindef
      inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
      if (irangp.ge.0) then
        call parcpt (inb)
        call parsom (flumty(ii))
      endif
      WRITE(NFECRA,7020) 'Indefini         ',II,INB,FLUMTY(II)

      do ii = 1, ntypmx
        if( ii.ne.iesicf .and.                                    &
             ii.ne.isspcf .and.                                   &
             ii.ne.ieqhcf .and.                                   &
             ii.ne.ierucf .and.                                   &
             ii.ne.isopcf .and.                                   &
             ii.ne.iparoi .and.                                   &
             ii.ne.isymet .and.                                   &
             ii.ne.iindef ) then
          inb = ifinty(ii,iphas)-idebty(ii,iphas)+1
          if (irangp.ge.0) then
            call parcpt (inb)
            call parsom (flumty(ii))
          endif
          if(inb.gt.0) then
            write(nfecra,7020)                                    &
                 'Type utilisateur ',II,INB,FLUMTY(II)
          endif
        endif
      enddo

    endif

    write(nfecra,7030)

  endif

enddo


!===============================================================================
! FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1099 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@    CONDITIONS AUX LIMITES INCORRECTES OU INCOMPLETES       ',/,&
'@                                                            ',/,&
'@    La valeur du type ITYPFB des conditions aux limites doit',/,&
'@      etre superieure ou egale a          1                 ',/,&
'@        et inferieure ou egale a ',I10                       ,/,&
'@    Une ou plusieurs erreurs sont listees ci-dessus.        ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage du sous-programme de definition des  ',/,&
'@    aux limites.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2010 format(/,' PHASE  : ',I10)
 2020 format(/,'   IFINTY : ',I10)
 2030 format(/,'   IDEBTY : ',I10)
 2040 format(/,'   ITYPFB : ',I10)
 2098 format(/,                                                   &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@    PROBLEME DE TRI DES FACES DE BORD                       ',/,&
'@                                                            ',/,&
'@    IFINTY(',I10   ,') = ',I10                               ,/,&
'@      est superieur a                                       ',/,&
'@    IDEBTY(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Contacter l''assistance.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2099 format(/,                                                   &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@    PROBLEME DE TRI DES FACES DE BORD SUR UN RANG DISTANT   ',/,&
'@                                                            ',/,&
'@    IFINTY(',I10   ,')                                      ',/,&
'@      est superieur a                                       ',/,&
'@    IDEBTY(',I10   ,')                                      ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Contacter l''assistance.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3099 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@    PROBLEME DE TRI DES FACES DE BORD                       ',/,&
'@                                                            ',/,&
'@    Phase          ',I10                                     ,/,&
'@      nombre de faces classees par type = ',I10              ,/,&
'@      nombre de faces de bord  NFABOR   = ',I10              ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Contacter l''assistance.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


 6010 format ( /,/,                                               &
 '   ** INFORMATIONS SUR LE TYPE DE FACES DE BORD',/,             &
 '      -----------------------------------------',/)
 6011 format (                                                    &
'   Phase : ',I4,                                               /,&
'---------------------------------------------------------------',&
'----------',                                                     &
                                                                /,&
'Type de bord           Code    Nb faces (tot/loc)  ',            &
' Index min   Index max',                                         &
                                                                /,&
'---------------------------------------------------------------',&
'----------')
 6020 format (                                                    &
 a17,i10,i10,3i12)
 6030 format(                                                     &
'---------------------------------------------------------------',&
'----------'/)

 6040 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@    CONDITIONS AUX LIMITES INCORRECTES OU INCOMPLETES       ',/,&
'@                                                            ',/,&
'@    Face de bord ',I10                                       ,/,&
'@    Face declaree en entree (ou sortie) a vitesse imposee,  ',/,&
'@      mais la valeur de la vitesse n''a pas ete fournie pour',/,&
'@      la composante ',A8                                     ,/,&
'@    (ce message ne s''affiche que pour la premiere face ou  ',/,&
'@      le probleme est detecte, mais d''autres faces peuvent ',/,&
'@      etre concernees)                                      ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier les conditions aux limites dans l''Interface   ',/,&
'@    ou dans le sous-programme utilisateur correspondant.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 6050 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@    CONDITIONS AUX LIMITES INCORRECTES OU INCOMPLETES       ',/,&
'@                                                            ',/,&
'@    Face de bord ',I10                                       ,/,&
'@    Face declaree en entree (ou sortie) a vitesse imposee   ',/,&
'@      avec un flux rentrant.                                ',/,&
'@    La valeur de la variable ',A8  ,' n''a pas ete specifiee',/,&
'@      (condition de Dirichlet).                             ',/,&
'@    (ce message ne s''affiche que pour la premiere face ou  ',/,&
'@      le probleme est detecte, mais d''autres faces peuvent ',/,&
'@      etre concernees)                                      ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute                           ',/,&
'@                                                            ',/,&
'@    Verifier les conditions aux limites dans l''Interface   ',/,&
'@    ou dans le sous-programme utilisateur correspondant.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 6060 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@    CONDITIONS AUX LIMITES INCORRECTES OU INCOMPLETES       ',/,&
'@                                                            ',/,&
'@    Au moins une face de bord declaree en entree            ',/,&
'@      (ou sortie) a vitesse imposee pour laquelle la valeur ',/,&
'@      de la vitesse n''a pas ete fournie pour toutes les    ',/,&
'@      composantes.                                          ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier les conditions aux limites dans l''Interface   ',/,&
'@    ou dans le sous-programme utilisateur correspondant.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 6070 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,&
'@    =========                                               ',/,&
'@    CONDITIONS AUX LIMITES INCORRECTES OU INCOMPLETES       ',/,&
'@                                                            ',/,&
'@    Au moins une face de bord declaree en entree            ',/,&
'@      (ou sortie) a vitesse imposee avec un flux rentrant   ',/,&
'@      pour laquelle la valeur d''une variable n''a pas ete  ',/,&
'@      specifiee (condition de Dirichlet).                   ',/,&
'@    Le calcul ne sera pas execute                           ',/,&
'@                                                            ',/,&
'@    Verifier les conditions aux limites dans l''Interface   ',/,&
'@    ou dans le sous-programme utilisateur correspondant.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


 7010 format ( /,/,                                               &
 '   ** INFORMATIONS SUR LE FLUX DE MASSE AU BORD',/,             &
 '      -----------------------------------------',/)
 7011 format (                                                    &
'   Phase : ',I4,                                               /,&
'---------------------------------------------------------------',&
                                                                /,&
'Type de bord           Code    Nb faces           Flux de masse',&
                                                                /,&
'---------------------------------------------------------------')
 7020 format (                                                    &
 a17,i10,i12,6x,e18.9)
 7030 format(                                                     &
'---------------------------------------------------------------',&
                                                                /)

 8000 format(/,                                                   &
'PHASE ',I4,' :                                               ',/,&
'Faces de bord d''entree/sortie libre detectees               ',/,&
'Mise a jour du point de reference pour la pression totale    ',/,&
' XYZP0 = ',E14.5,E14.5,E14.5                  ,/)
 8001 format(/,                                                   &
'PHASE ',I4,' :                                               ',/,&
'Faces de bord a Dirichlet de pression impose detectees       ',/,&
'Mise a jour du point de reference pour la pression totale    ',/,&
' XYZP0 = ',E14.5,E14.5,E14.5                  ,/)

!-------------------------------------------------------------------------------

#else

 1099 format(                                                     &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,&
'@    ========'                                                ,/,&
'@    INCORRECT OR INCOMPLETE BOUNDARY CONDITIONS'             ,/,&
'@'                                                            ,/,&
'@    Value of type ITYPFB for boundary conditions must be'    ,/,&
'@          greater or equal to          1'                    ,/,&
'@      and less than or equal to ', I10                       ,/,&
'@    One or more errors are listed above.'                    ,/,&
'@'                                                            ,/,&
'@    The calculation will not be run.'                        ,/,&
'@'                                                            ,/,&
'@    Verify the boundary condition definitions in the'        ,/,&
'@    appropriate user subroutine.'                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

 2010 format(/,' PHASE  : ',I10)
 2020 format(/,'   IFINTY : ',I10)
 2030 format(/,'   IDEBTY : ',I10)
 2040 format(/,'   ITYPFB : ',I10)
 2098 format(/,                                                   &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,&
'@    ========'                                                ,/,&
'@    PROBLEM WITH ORDERING OF BOUNDARY FACES'                 ,/,&
'@'                                                            ,/,&
'@    IFINTY(',I10   ,') = ',I10                               ,/,&
'@      is greater than'                                       ,/,&
'@    IDEBTY(',I10   ,') = ',I10                               ,/,&
'@'                                                            ,/,&
'@    The calculation will not be run.'                        ,/,&
'@'                                                            ,/,&
'@    Contact support.'                                        ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2099 format(/,                                                   &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,&
'@    ========'                                                ,/,&
'@    PROBLEM WITH ORDERING OF BOUNDARY FACES'                 ,/,&
'@    ON A DISTANT RANK.'                                      ,/,&
'@'                                                            ,/,&
'@    IFINTY(',I10   ,')                                      ',/,&
'@      is greater than'                                       ,/,&
'@    IDEBTY(',I10   ,')                                      ',/,&
'@'                                                            ,/,&
'@    The calculation will not be run.'                        ,/,&
'@'                                                            ,/,&
'@    Contact support.'                                        ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3099 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,&
'@    ========'                                                ,/,&
'@    PROBLEM WITH ORDERING OF BOUNDARY FACES'                 ,/,&
'@'                                                            ,/,&
'@    Phase          ',I10                                     ,/,&
'@      number of faces classified by type = ',I10             ,/,&
'@      number of boundary faces (NFABOR)  = ',I10             ,/,&
'@'                                                            ,/,&
'@    The calculation will not be run.'                        ,/,&
'@'                                                            ,/,&
'@    Contact support.'                                        ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


 6010 format ( /,/,                                               &
 '   ** INFORMATION ON BOUNDARY FACES TYPE',/,                    &
 '      ----------------------------------',/)
 6011 format (                                                    &
'   Phase : ',I4,                                               /,&
'---------------------------------------------------------------',&
'----------',                                                     &
                                                                /,&
'Boundary type          Code    Nb faces (tot/loc)  ',            &
' Min index   Max index',                                         &
                                                                /,&
'---------------------------------------------------------------',&
'----------')
 6020 format (                                                    &
 a17,i10,i10,3i12)
 6030 format(                                                     &
'---------------------------------------------------------------',&
'----------'/)

 6040 format(                                                           &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,&
'@    ========'                                                ,/,&
'@    INCORRECT OR INCOMPLETE BOUNDARY CONDITIONS'             ,/,&
'@'                                                            ,/,&
'@    Boundary face ',I10                                      ,/,&
'@    Face declared as inlet (or outlet) with prescribed'      ,/,&
'@      velocity, but velocity value has not been given for'   ,/,&
'@      component ',A8                                         ,/,&
'@    (this message is printed only for the first face on'     ,/,&
'@      which this problem is detected, but other faces may'   ,/,&
'@      be concerned)'                                         ,/,&
'@'                                                            ,/,&
'@    The calculation will not be run.                        ',/,&
'@'                                                            ,/,&
'@    Verify the boundary condition definitions in the GUI'    ,/,&
'@    or in the appropriate user subroutine.'                  ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 6050 format(                                                           &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,&
'@    ========'                                                ,/,&
'@    INCORRECT OR INCOMPLETE BOUNDARY CONDITIONS'             ,/,&
'@'                                                            ,/,&
'@    Boundary face ',I10                                      ,/,&
'@    Face declared as inlet (or outlet) with prescribed'      ,/,&
'@      velocity with entering flow.'                          ,/,&
'@    Value of variable ',A8  ,' has not been specified'       ,/,&
'@      (Dirichlet condition).'                                ,/,&
'@    (this message is printed only for the first face on'     ,/,&
'@      which this problem is detected, but other faces may'   ,/,&
'@      be concerned)'                                         ,/,&
'@'                                                            ,/,&
'@    The calculation will not be run.                        ',/,&
'@'                                                            ,/,&
'@    Verify the boundary condition definitions in the GUI'    ,/,&
'@    or in the appropriate user subroutine.'                  ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 6060 format(                                                           &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,&
'@    ========'                                                ,/,&
'@    INCORRECT OR INCOMPLETE BOUNDARY CONDITIONS'             ,/,&
'@'                                                            ,/,&
'@    At least one boundary face declared as inlet (or'        ,/,&
'@      outlet) with prescribed velocity for which the'        ,/,&
'@      velocity value has not been assigned for all'          ,/,&
'@      components.'                                           ,/,&
'@    The calculation will not be run.                        ',/,&
'@'                                                            ,/,&
'@    Verify the boundary condition definitions in the GUI'    ,/,&
'@    or in the appropriate user subroutine.'                  ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 6070 format(                                                     &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,&
'@    ========'                                                ,/,&
'@    INCORRECT OR INCOMPLETE BOUNDARY CONDITIONS'             ,/,&
'@'                                                            ,/,&
'@    At least one boundary face declared as inlet (or'        ,/,&
'@      outlet) with prescribed velocity with an entering'     ,/,&
'@      flow for which the value of a variable has not been'   ,/,&
'@      specified (Dirichlet condition).'                      ,/,&
'@    The calculation will not be run.                        ',/,&
'@'                                                            ,/,&
'@    Verify the boundary condition definitions in the GUI'    ,/,&
'@    or in the appropriate user subroutine.'                  ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


 7010 format ( /,/,                                               &
 '   ** BOUNDARY MASS FLOW INFORMATION',/,                        &
 '      ------------------------------',/)
 7011 format (                                                    &
'   Phase : ',I4,                                               /,&
'---------------------------------------------------------------',&
                                                                /,&
'Boundary type          Code    Nb faces           Mass flow'   , &
                                                                /,&
'---------------------------------------------------------------')
 7020 format (                                                    &
 a17,i10,i12,6x,e18.9)
 7030 format(                                                     &
'---------------------------------------------------------------',&
                                                                /)

 8000 format(/,                                                   &
'PHASE ',I4,' :'                                               ,/,&
'Boundary faces with free inlet/outlet detected'               ,/,&
'Update of reference point for total pressure'                 ,/,&
' XYZP0 = ',E14.5,E14.5,E14.5                  ,/)
 8001 format(/,                                                   &
'PHASE ',I4,' :'                                               ,/,&
'Boundary faces with pressure Dirichlet condition detected'    ,/,&
'Update of reference point for total pressure'                 ,/,&
' XYZP0 = ',E14.5,E14.5,E14.5                  ,/)

#endif


return
end
