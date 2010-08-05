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

subroutine inimas &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ivar1  , ivar2  , ivar3  , imaspe , iphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgu , imligu , &
   iwarnu , nfecra ,                                              &
   epsrgu , climgu , extrau ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , isympa ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rom    , romb   ,                                              &
   ux     , uy     , uz     ,                                     &
   coefax , coefay , coefaz , coefbx , coefby , coefbz ,          &
   flumas , flumab ,                                              &
   dpdx   , dpdy   , dpdz   , dpdxa  , dpdya  , dpdza  ,          &
   qdmx   , qdmy   , qdmz   , coefqa ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------
!                                                              -->
! INCREMENTATION DU FLUX DE MASSE A PARTIR DU CHAMP VECTORIEL ROM.U
!  .     .          -->   -->
!  m   = m   +(rom* U ) . n
!   ij    ij         ij    ij


! Pour la reconstruction, grad(rho u) est calcule avec des
! conditions aux limites approchees :
!    COEFA(rho u) = ROMB * COEFA(u)
!    COEFB(rho u) = COEFB (u)

! et pour le flux de masse au bord on ecrit
!  FLUMAB = [ROMB*COEFA(u) + ROMB*COEFB(u)*Ui
!                + COEFB(u)*II'.grad(rho u) ].Sij
! ce qui utilise de petites approximations sur les
! non-orthogonalites (cf. notice)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ivar1            ! e  ! <-- ! variable de la direction 1                     !
! ivar2            ! e  ! <-- ! variable de la direction 2                     !
! ivar3            ! e  ! <-- ! variable de la direction 3                     !
! imaspe           ! e  ! <-- ! suivant l'appel de inimas                      !
!                  !    !     ! = 1 si appel de navsto resolp                  !
!                  !    !     ! = 2 si appel de divrij                         !
! iphas            ! i  ! <-- ! phase number                                   !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iflmb0           ! e  ! <-- ! =1 : flux de masse annule sym-paroi            !
! init             ! e  ! <-- ! > 0 : initialisation du flux de masse          !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! imrgra           ! e  ! <-- ! indicateur = 0 gradrc 97                       !
!                  ! e  ! <-- !            = 1 gradmc 99                       !
! iccocg           ! e  ! <-- ! indicateur = 1 pour recalcul de cocg           !
!                  !    !     !              0 sinon                           !
! nswrgu           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligu           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! iwarnu           ! e  ! <-- ! niveau d'impression                            !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! epsrgu           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgu           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrau           ! r  ! <-- ! coef extrap gradient                           !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! isympa           ! te ! <-- ! zero pour annuler le flux de masse             !
! (nfabor     )    !    !     !(symetries et parois avec cl couplees)          !
!                  !    !     ! un sinon                                       !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! rom(ncelet       ! tr ! <-- ! masse volumique aux cellules                   !
! romb(nfabor)     ! tr ! <-- ! masse volumique aux bords                      !
! ux,y,z(ncelet    ! tr ! <-- ! vitesse                                        !
! coefax, b        ! tr ! <-- ! tableaux des cond lim pour ux, uy, uz          !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! dpd. (ncelet     ! tr ! --- ! tableau de travail pour le grad de p           !
! dpd.a(ncelet     ! tr ! --- ! tableau de travail pour le grad de p           !
!                  !    !     !  avec decentrement amont                       !
! qdm.(ncelet)     ! tr ! --- ! tableau de travail pour la qdm                 !
! coefqa(nfab,3    ! tr ! --- ! tableau de travail cl de la qdm                !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "dimfbr.f90"
include "paramx.f90"
include "pointe.f90"
include "period.f90"
include "parall.f90"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          ivar1  , ivar2  , ivar3  , imaspe , iphas
integer          iflmb0 , init   , inc    , imrgra , iccocg
integer          nswrgu , imligu
integer          iwarnu , nfecra
double precision epsrgu , climgu , extrau

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr), isympa(nfabor)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rom(ncelet), romb(nfabor)
double precision ux(ncelet), uy(ncelet), uz(ncelet)
double precision coefax(nfabor), coefay(nfabor), coefaz(nfabor)
double precision coefbx(nfabor), coefby(nfabor), coefbz(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision dpdx (ncelet),dpdy (ncelet),dpdz (ncelet)
double precision dpdxa(ncelet),dpdya(ncelet),dpdza(ncelet)
double precision qdmx(ncelet),qdmy(ncelet),qdmz(ncelet)
double precision coefqa(ndimfb,3)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, ii, jj, iel, iof, iii
integer          iappel, iphydp
double precision pfac,pip,uxfac,uyfac,uzfac
double precision dofx,dofy,dofz,pond
double precision diipbx, diipby, diipbz

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! ---> CALCUL DE LA QTE DE MOUVEMENT


if( init.eq.1 ) then
  do ifac = 1, nfac
    flumas(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    flumab(ifac) = 0.d0
  enddo

elseif(init.ne.0) then
  write(nfecra,1000) init
  call csexit (1)
endif

do iel = 1, ncel
  qdmx(iel) = rom(iel)*ux(iel)
  qdmy(iel) = rom(iel)*uy(iel)
  qdmz(iel) = rom(iel)*uz(iel)
enddo

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(qdmx, qdmy, qdmz)
  !==========
endif

do ifac =1, nfabor
  coefqa(ifac,1) = romb(ifac)*coefax(ifac)
  coefqa(ifac,2) = romb(ifac)*coefay(ifac)
  coefqa(ifac,3) = romb(ifac)*coefaz(ifac)
enddo

!===============================================================================
! 2.  CALCUL DU FLUX DE MASSE SANS TECHNIQUE DE RECONSTRUCTION
!===============================================================================

if( nswrgu.le.1 ) then

!     FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pond = ra(ipond-1+ifac)

    flumas(ifac) =  flumas(ifac)                                  &
     +(pond*qdmx(ii)+(1.d0-pond)*qdmx(jj) )*surfac(1,ifac)        &
     +(pond*qdmy(ii)+(1.d0-pond)*qdmy(jj) )*surfac(2,ifac)        &
     +(pond*qdmz(ii)+(1.d0-pond)*qdmz(jj) )*surfac(3,ifac)

  enddo


!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    uxfac = inc*coefqa(ifac,1) +coefbx(ifac)*romb(ifac)*ux(ii)
    uyfac = inc*coefqa(ifac,2) +coefby(ifac)*romb(ifac)*uy(ii)
    uzfac = inc*coefqa(ifac,3) +coefbz(ifac)*romb(ifac)*uz(ii)

    flumab(ifac) = flumab(ifac)                                   &
     +( uxfac*surfbo(1,ifac)                                      &
                    +uyfac*surfbo(2,ifac) +uzfac*surfbo(3,ifac) )

  enddo

endif


!===============================================================================
! 4.  CALCUL DU FLUX DE MASSE AVEC TECHNIQUE DE RECONSTRUCTION
!        SI LE MAILLAGE EST NON ORTHOGONAL
!===============================================================================

if( nswrgu.gt.1 ) then




!     TRAITEMENT DE LA PERIODICITE SPEFICIQUE A INIMAS AU DEBUT
!     =========================================================

  if(iperot.gt.0) then
    iappel = 1

    call permas                                                   &
    !==========
 ( imaspe , iphas  , iappel ,                                     &
   rom    ,                                                       &
   ra(idudxy) , ra(idrdxy) , ra(iwdudx) , ra(iwdrdx) )

  endif

!     FLUX DE MASSE SUIVANT X
!     =======================

! ---> CALCUL DU GRADIENT

  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar1  , imrgra , inc    , iccocg , nswrgu , imligu , iphydp , &
   iwarnu , nfecra , epsrgu , climgu , extrau ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dpdxa  , dpdxa  , dpdxa  ,                                     &
   qdmx   , coefqa(1,1) , coefbx ,                                &
   dpdx   , dpdy   , dpdz   ,                                     &
!        ------   ------   ------
   dpdxa  , dpdya  , dpdza  ,                                     &
   rdevel , rtuser , ra     )


! ---> FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pond = ra(ipond-1+ifac)

    iof = idofij-1+3*(ifac-1)
!---> DOF = OF
    dofx = ra(iof+1)
    dofy = ra(iof+2)
    dofz = ra(iof+3)

    flumas(ifac) = flumas(ifac)                                   &
         +( pond*qdmx(ii) +(1.d0-pond)*qdmx(jj)                   &
           +0.5d0*( dpdx(ii) +dpdx(jj) )*dofx                     &
           +0.5d0*( dpdy(ii) +dpdy(jj) )*dofy                     &
           +0.5d0*( dpdz(ii) +dpdz(jj) )*dofz    )*surfac(1,ifac)

  enddo

! ---> FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    iii = idiipb-1+3*(ifac-1)
    diipbx = ra(iii+1)
    diipby = ra(iii+2)
    diipbz = ra(iii+3)

    pip = romb(ifac) * ux(ii)                                     &
          +dpdx(ii)*diipbx                                        &
          +dpdy(ii)*diipby +dpdz(ii)*diipbz
    pfac = inc*coefqa(ifac,1) +coefbx(ifac)*pip

    flumab(ifac) = flumab(ifac) +pfac*surfbo(1,ifac)

  enddo


!     FLUX DE MASSE SUIVANT Y
!     =======================

! ---> CALCUL DU GRADIENT

  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar2  , imrgra , inc    , iccocg , nswrgu , imligu , iphydp , &
   iwarnu , nfecra , epsrgu , climgu , extrau ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dpdxa  , dpdxa  , dpdxa  ,                                     &
   qdmy   , coefqa(1,2) , coefby ,                                &
   dpdx   , dpdy   , dpdz   ,                                     &
!        ------   ------   ------
   dpdxa  , dpdya  , dpdza  ,                                     &
   rdevel , rtuser , ra     )


! ---> FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pond = ra(ipond-1+ifac)

    iof = idofij-1+3*(ifac-1)

!---> DOF = OF

    dofx = ra(iof+1)
    dofy = ra(iof+2)
    dofz = ra(iof+3)

    flumas(ifac) = flumas(ifac)                                   &
         +( pond*qdmy(ii) +(1.d0-pond)*qdmy(jj)                   &
           +0.5d0*( dpdx(ii) +dpdx(jj) )*dofx                     &
           +0.5d0*( dpdy(ii) +dpdy(jj) )*dofy                     &
           +0.5d0*( dpdz(ii) +dpdz(jj) )*dofz    )*surfac(2,ifac)

  enddo

! ---> FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    iii = idiipb-1+3*(ifac-1)
    diipbx = ra(iii+1)
    diipby = ra(iii+2)
    diipbz = ra(iii+3)

    pip = romb(ifac) * uy(ii)                                     &
        +dpdx(ii)*diipbx                                          &
        +dpdy(ii)*diipby +dpdz(ii)*diipbz
    pfac = inc*coefqa(ifac,2) +coefby(ifac)*pip

    flumab(ifac) = flumab(ifac) +pfac*surfbo(2,ifac)

  enddo

!     FLUX DE MASSE SUIVANT Z
!     =======================

! ---> CALCUL DU GRADIENT

  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar3  , imrgra , inc    , iccocg , nswrgu , imligu , iphydp , &
   iwarnu , nfecra , epsrgu , climgu , extrau ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dpdxa  , dpdxa  , dpdxa  ,                                     &
   qdmz     , coefqa(1,3) , coefbz ,                              &
   dpdx   , dpdy   , dpdz   ,                                     &
!        ------   ------   ------
   dpdxa  , dpdya  , dpdza  ,                                     &
   rdevel , rtuser , ra     )

!     FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pond = ra(ipond-1+ifac)

    iof = idofij-1+3*(ifac-1)

!---> DOF = OF

    dofx = ra(iof+1)
    dofy = ra(iof+2)
    dofz = ra(iof+3)

    flumas(ifac) = flumas(ifac)                                   &
         +( pond*qdmz(ii) +(1.d0-pond)*qdmz(jj)                   &
           +0.5d0*( dpdx(ii) +dpdx(jj) )*dofx                     &
           +0.5d0*( dpdy(ii) +dpdy(jj) )*dofy                     &
           +0.5d0*( dpdz(ii) +dpdz(jj) )*dofz    )*surfac(3,ifac)

  enddo

! ---> FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    iii = idiipb-1+3*(ifac-1)
    diipbx = ra(iii+1)
    diipby = ra(iii+2)
    diipbz = ra(iii+3)

    pip = romb(ifac) * uz(ii)                                     &
        +dpdx(ii)*diipbx                                          &
        +dpdy(ii)*diipby +dpdz(ii)*diipbz
    pfac = inc*coefqa(ifac,3) +coefbz(ifac)*pip

    flumab(ifac) = flumab(ifac) +pfac*surfbo(3,ifac)

  enddo



!     TRAITEMENT DE LA PERIODICITE SPEFICIQUE A INIMAS A LA FIN
!     =========================================================

  if(iperot.gt.0) then
    iappel = 2

    call permas                                                   &
    !==========
 ( imaspe , iphas  , iappel ,                                     &
   rom    ,                                                       &
   ra(idudxy) , ra(idrdxy) , ra(iwdudx) , ra(iwdrdx) )

  endif




endif


!===============================================================================
! 6.  POUR S'ASSURER DE LA NULLITE DU FLUX DE MASSE AUX LIMITES
!       SYMETRIES PAROIS COUPLEES
!===============================================================================

if(iflmb0.eq.1) then
! FORCAGE DE FLUMAB a 0 pour la vitesse'
  do ifac = 1, nfabor
    if(isympa(ifac).eq.0) then
      flumab(ifac) = 0.d0
    endif
  enddo
endif

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format('INIMAS APPELE AVEC INIT =',I10)

#else

 1000 format('INIMAS CALLED WITH INIT =',I10)

#endif

!----
! FIN
!----

return

end subroutine
