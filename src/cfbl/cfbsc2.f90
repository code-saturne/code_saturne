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

subroutine cfbsc2 &
!================

 ( nvar   , nscal  ,                                              &
   ivar   , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   pvar   , coefap , coefbp , cofafp , cofbfp ,                   &
   flumas , flumab , viscf  , viscb  ,                            &
   smbrp  )

!===============================================================================
! FONCTION :
! ---------

! CALCUL DU BILAN EXPLICITE DE LA VARIABLE PVAR (VITESSE,SCALAIRES)

!                   -- .                  ----->        -->
! SMBRP = SMBRP -  \   m   PVAR +Visc   ( grad PVAR )  . n
!                  /__  ij    ij     ij             ij    ij
!                   j

! ATTENTION : SMBRP DEJA INITIALISE AVANT L'APPEL A BILSCA
!            IL CONTIENT LES TERMES SOURCES EXPLICITES, ETC....

! BLENCP = 1 : PAS D'UPWIND EN DEHORS DU TEST DE PENTE
! BLENCP = 0 : UPWIND
! ISCHCP = 1 : CENTRE
! ISCHCP = 0 : SECOND ORDER

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ivar             ! e  ! <-- ! numero de la variable                          !
! iconvp           ! e  ! <-- ! indicateur = 1 convection, 0 sinon             !
! idiffp           ! e  ! <-- ! indicateur = 1 diffusion , 0 sinon             !
! nswrgp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligp           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! ircflp           ! e  ! <-- ! indicateur = 1 rec flux, 0 sinon               !
! ischcp           ! e  ! <-- ! indicateur = 1 centre , 0 2nd order            !
! isstpp           ! e  ! <-- ! indicateur = 1 sans test de pente              !
!                  !    !     !            = 0 avec test de pente              !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! imrgra           ! e  ! <-- ! indicateur = 0 gradrc 97                       !
!                  ! e  ! <-- !            = 1 gradmc 99                       !
! iccocg           ! e  ! <-- ! indicateur = 1 pour recalcul de cocg           !
!                  !    !     !              0 sinon                           !
! ipp              ! e  ! <-- ! numero de variable pour post                   !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! blencp           ! r  ! <-- ! 1 - proportion d'upwind                        !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! pvar (ncelet     ! tr ! <-- ! variable resolue (instant precedent)           !
! coefap, b        ! tr ! <-- ! tableaux des cond lim pour p                   !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! cofafp, b        ! tr ! <-- ! tableaux des cond lim pour le flux de          !
!   (nfabor)       !    !     !  diffusion de p                                !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! viscf (nfac)     ! tr ! <-- ! visc*surface/dist aux faces internes           !
!                  !    !     !  pour second membre                            !
! viscb (nfabor    ! tr ! <-- ! visc*surface/dist aux faces de bord            !
!                  !    !     !  pour second membre                            !
! smbrp(ncelet     ! tr ! <-- ! bilan au second membre                         !
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
use pointe
use entsor
use parall
use period
use cfpoin, only: ifbrus
use mesh
use numvar, only: ivarfl

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ivar   , iconvp , idiffp , nswrgp , imligp
integer          ircflp , ischcp , isstpp
integer          inc    , imrgra , iccocg
integer          iwarnp , ipp
double precision blencp , epsrgp , climgp, extrap


double precision pvar (ncelet), coefap(nfabor), coefbp(nfabor)
double precision                cofafp(nfabor), cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscf (nfac), viscb (nfabor)
double precision smbrp(ncelet)

! Local variables

character*80     chaine
character*8      cnom
integer          ifac,ii,jj,infac,iel,iupwin, iij, iii, iok
integer          idimtr, irpvar

double precision pfac,pfacd,pip,pjp,flui,fluj,flux
double precision difx,dify,difz,djfx,djfy,djfz,pif,pjf
double precision testi,testj,testij
double precision dpxf,dpyf,dpzf, pi, pj
double precision dcc, ddi, ddj, tesqck
double precision dijpfx, dijpfy, dijpfz
double precision diipfx, diipfy, diipfz
double precision djjpfx, djjpfy, djjpfz
double precision diipbx, diipby, diipbz
double precision pnd, srfan
double precision pfac1, pfac2, pfac3, unsvol

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: dpdxa, dpdya, dpdza

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(grad(ncelet,3))
allocate(dpdxa(ncelet), dpdya(ncelet), dpdza(ncelet))

! Initialize variables to avoid compiler warnings

pif = 0.d0
pjf = 0.d0

! Memoire


chaine = nomvar(ipp)
cnom   = chaine(1:8)

if(iwarnp.ge.2) then
  if (ischcp.eq.1) then
    WRITE(NFECRA,1000)CNOM,'    CENTRE ',                         &
                                              (1.d0-blencp)*100.d0
  else
    WRITE(NFECRA,1000)CNOM,' 2ND ORDER ',                         &
                                              (1.d0-blencp)*100.d0
  endif
endif

iupwin = 0
if(blencp.eq.0.d0) iupwin = 1

!===============================================================================
! 2.  CALCUL DU BILAN AVEC TECHNIQUE DE RECONSTRUCTION
!===============================================================================

! ======================================================================
! ---> CALCUL DU GRADIENT DE P
! ======================================================================
!    GRAD sert a la fois pour la reconstruction des flux et pour le test
!    de pente. On doit donc le calculer :
!        - quand on a de la diffusion et qu'on reconstruit les flux
!        - quand on a de la convection SOLU
!        - quand on a de la convection, qu'on n'est pas en upwind pur
!          et qu'on reconstruit les flux
!        - quand on a de la convection, qu'on n'est pas en upwind pur
!          et qu'on n'a pas shunte le test de pente

if( (idiffp.ne.0 .and. ircflp.eq.1) .or.                          &
    (iconvp.ne.0 .and. iupwin.eq.0 .and.                          &
       (ischcp.eq.0 .or. ircflp.eq.1 .or. isstpp.eq.0)) ) then

  call grdcel                                                     &
  !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   pvar   , coefap , coefbp ,                                     &
   grad   )

else
  do iel = 1, ncelet
    grad(iel,1) = 0.d0
    grad(iel,2) = 0.d0
    grad(iel,3) = 0.d0
  enddo
endif


! ======================================================================
! ---> CALCUL DU GRADIENT DECENTRE DPDXA, DPDYA, DPDZA POUR TST DE PENTE
! ======================================================================

do iel = 1, ncelet
  dpdxa(iel) = 0.d0
  dpdya(iel) = 0.d0
  dpdza(iel) = 0.d0
enddo

if( iconvp.gt.0.and.iupwin.eq.0.and.isstpp.eq.0 ) then

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    difx = cdgfac(1,ifac) - xyzcen(1,ii)
    dify = cdgfac(2,ifac) - xyzcen(2,ii)
    difz = cdgfac(3,ifac) - xyzcen(3,ii)
    djfx = cdgfac(1,ifac) - xyzcen(1,jj)
    djfy = cdgfac(2,ifac) - xyzcen(2,jj)
    djfz = cdgfac(3,ifac) - xyzcen(3,jj)

    pif = pvar(ii) +difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
    pjf = pvar(jj) +djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

    pfac = pjf
    if( flumas(ifac ).gt.0.d0 ) pfac = pif

    pfac1 = pfac*surfac(1,ifac )
    pfac2 = pfac*surfac(2,ifac )
    pfac3 = pfac*surfac(3,ifac )

    dpdxa(ii) = dpdxa(ii) +pfac1
    dpdya(ii) = dpdya(ii) +pfac2
    dpdza(ii) = dpdza(ii) +pfac3

    dpdxa(jj) = dpdxa(jj) -pfac1
    dpdya(jj) = dpdya(jj) -pfac2
    dpdza(jj) = dpdza(jj) -pfac3

  enddo

  do ifac =1,nfabor
    ii = ifabor(ifac )
    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)
    pfac = inc*coefap(ifac )                                    &
         +coefbp(ifac )*(pvar(ii)+diipbx*grad(ii,1)             &
         +diipby*grad(ii,2)+diipbz*grad(ii,3) )
    dpdxa(ii) = dpdxa(ii) +pfac*surfbo(1,ifac )
    dpdya(ii) = dpdya(ii) +pfac*surfbo(2,ifac )
    dpdza(ii) = dpdza(ii) +pfac*surfbo(3,ifac )
  enddo

  do iel = 1, ncel
    unsvol = 1.d0/volume(iel)
    dpdxa(iel) = dpdxa(iel)*unsvol
    dpdya(iel) = dpdya(iel)*unsvol
    dpdza(iel) = dpdza(iel)*unsvol
  enddo

  ! Synchronization for parallelism or periodicity

  if (irangp.ge.0 .or. iperio.eq.1) then
    call synvec(dpdxa, dpdya, dpdza)
    !==========
  endif

  if (iperot.eq.1.and.ivar.gt.0) then
    call pering(ivarfl(ivar), idimtr, dpdxa, dpdya, dpdza)
    !==========
  endif

endif


! ======================================================================
! ---> ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
! ======================================================================

infac = 0

if(ncelet.gt.ncel) then
  do iel = ncel+1, ncelet
    smbrp(iel) = 0.d0
  enddo
endif


!  --> FLUX UPWIND PUR
!  =====================

if(iupwin.eq.1) then

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    dijpfx = dijpf(1,ifac)
    dijpfy = dijpf(2,ifac)
    dijpfz = dijpf(3,ifac)

    pnd   = pond(ifac)

! ON RECALCULE A CE NIVEAU II' ET JJ'

    diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
    diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
    diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
    djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd  * dijpfx
    djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd  * dijpfy
    djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd  * dijpfz

    dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
    dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
    dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

    pi = pvar(ii)
    pj = pvar(jj)

    pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
    pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

    flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
    fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )

    infac = infac+1

    flux = iconvp*( flui*pi  +fluj*pj  ) + idiffp*viscf(ifac)*( pip -pjp )

    smbrp(ii) = smbrp(ii) - (flux - iconvp*pi*flumas(ifac))
    smbrp(jj) = smbrp(jj) + (flux - iconvp*pj*flumas(ifac))

  enddo


!  --> FLUX SANS TEST DE PENTE
!  ============================

elseif(isstpp.eq.1) then

  iok = 0
  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    dijpfx = dijpf(1,ifac)
    dijpfy = dijpf(2,ifac)
    dijpfz = dijpf(3,ifac)

    pnd   = pond(ifac)

! ON RECALCULE II' ET JJ'

    diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
    diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
    diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
    djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd  * dijpfx
    djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd  * dijpfy
    djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd  * dijpfz

    dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
    dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
    dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

    pi = pvar(ii)
    pj = pvar(jj)

    pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
    pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

    flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
    fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )


!         CENTRE
!        --------

    if (ischcp.eq.1) then

      pif = pnd*pip +(1.d0-pnd)*pjp
      pjf = pif


!         SECOND ORDER
!        --------------

    elseif(ischcp.eq.0) then

      difx = cdgfac(1,ifac) - xyzcen(1,ii)
      dify = cdgfac(2,ifac) - xyzcen(2,ii)
      difz = cdgfac(3,ifac) - xyzcen(3,ii)
      djfx = cdgfac(1,ifac) - xyzcen(1,jj)
      djfy = cdgfac(2,ifac) - xyzcen(2,jj)
      djfz = cdgfac(3,ifac) - xyzcen(3,jj)

!     on laisse la reconstruction de PIF et PJF meme si IRCFLP=0
!     sinon cela revient a faire de l'upwind
      pif = pvar(ii) + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
      pjf = pvar(jj) + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

    else
      write(nfecra,9000)ischcp
      iok = 1
    endif


!        BLENDING
!       ----------

    pif = blencp*pif+(1.d0-blencp)*pvar(ii)
    pjf = blencp*pjf+(1.d0-blencp)*pvar(jj)


!        FLUX
!       ------

    flux = iconvp*( flui*pif +fluj*pjf ) + idiffp*viscf(ifac)*( pip -pjp )


!        ASSEMBLAGE
!       ------------

    smbrp(ii) = smbrp(ii) - (flux - pi*flumas(ifac))
    smbrp(jj) = smbrp(jj) + (flux - pj*flumas(ifac))

  enddo
!        Le call csexit ne doit pas etre dans la boucle, sinon
!         le VPP la devectorise (pour la boucle non vectorisee forcee,
!         ce n'est pas grave, c'est juste pour recopier exactement
!         la precedente)
  if(iok.ne.0) then
    call csexit (1)
  endif




!  --> FLUX AVEC TEST DE PENTE (separe pour vectorisation)
!  =============================

else

  iok = 0
  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    dijpfx = dijpf(1,ifac)
    dijpfy = dijpf(2,ifac)
    dijpfz = dijpf(3,ifac)

    pnd    = pond(ifac)
    srfan  = surfan(ifac)

! ON RECALCULE II' ET JJ'

    diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
    diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
    diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
    djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd  * dijpfx
    djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd  * dijpfy
    djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd  * dijpfz

    dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
    dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
    dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

    pi = pvar(ii)
    pj = pvar(jj)

    pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
    pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

    flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
    fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )


!         TEST DE PENTE
!        ---------------

    testi = dpdxa(ii)*surfac(1,ifac) +dpdya(ii)*surfac(2,ifac)  &
          + dpdza(ii)*surfac(3,ifac)
    testj = dpdxa(jj)*surfac(1,ifac) +dpdya(jj)*surfac(2,ifac)  &
          + dpdza(jj)*surfac(3,ifac)
    testij= dpdxa(ii)*dpdxa(jj)    +dpdya(ii)*dpdya(jj)         &
          + dpdza(ii)*dpdza(jj)

    if( flumas(ifac).gt.0.d0) then
      dcc = grad(ii,1)*surfac(1,ifac) +grad(ii,2)*surfac(2,ifac)    &
          + grad(ii,3)*surfac(3,ifac)
      ddi = testi
      ddj = ( pvar(jj)-pvar(ii) )/dist(ifac) *srfan
    else
      dcc = grad(jj,1)*surfac(1,ifac) +grad(jj,2)*surfac(2,ifac)    &
          + grad(jj,3)*surfac(3,ifac)
      ddi = ( pvar(jj)-pvar(ii) )/dist(ifac) *srfan
      ddj = testj
    endif
    tesqck = dcc**2 -(ddi-ddj)**2


!         UPWIND
!        --------

!MO          IF( (TESTI*TESTJ).LE.0.D0 .OR. TESTIJ.LE.0.D0 ) THEN
    if( tesqck.le.0.d0 .or. testij.le.0.d0 ) then

      pif = pvar(ii)
      pjf = pvar(jj)
      infac = infac+1

    else


!         CENTRE
!        --------

      if (ischcp.eq.1) then

        pif = pnd*pip +(1.d0-pnd)*pjp
        pjf = pif


!         SECOND ORDER
!        --------------

      elseif(ischcp.eq.0) then

        difx = cdgfac(1,ifac) - xyzcen(1,ii)
        dify = cdgfac(2,ifac) - xyzcen(2,ii)
        difz = cdgfac(3,ifac) - xyzcen(3,ii)
        djfx = cdgfac(1,ifac) - xyzcen(1,jj)
        djfy = cdgfac(2,ifac) - xyzcen(2,jj)
        djfz = cdgfac(3,ifac) - xyzcen(3,jj)

!     on laisse la reconstruction de PIF et PJF meme si IRCFLP=0
!     sinon cela revient a faire de l'upwind
        pif = pvar(ii) + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
        pjf = pvar(jj) + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

      else
        write(nfecra,9000)ischcp
        iok = 1
      endif

    endif


!        BLENDING
!       ----------

    pif = blencp*pif+(1.d0-blencp)*pvar(ii)
    pjf = blencp*pjf+(1.d0-blencp)*pvar(jj)


!        FLUX
!       ------

    flux = iconvp*( flui*pif +fluj*pjf ) + idiffp*viscf(ifac)*( pip -pjp )


!        ASSEMBLAGE
!       ------------

    smbrp(ii) = smbrp(ii) - (flux - pi*flumas(ifac))
    smbrp(jj) = smbrp(jj) + (flux - pj*flumas(ifac))

  enddo
!        Le call csexit ne doit pas etre dans la boucle, sinon
!         le VPP la devectorise (pour la boucle non vectorisee forcee,
!         ce n'est pas grave, c'est juste pour recopier exactement
!         la precedente)
  if(iok.ne.0) then
    call csexit (1)
  endif

endif



if(iwarnp.ge.2) then
  write(nfecra,1100)cnom,infac,nfac
endif


! ======================================================================
! ---> ASSEMBLAGE A PARTIR DES FACETTES DE BORD
! ======================================================================

!     Lorsque IIFBRU.GT.0, on ne prend pas en compte le flux convectif
!       sur les faces pour lesquelles on dispose par ailleurs d'un
!       flux de C.L.

!     Lorsque IIFBRU.LE.0, on adopte le schéma standard (i.e. on prend
!       en compte le flux convectitf).

!     Dans les deux cas, si on prend en compte le flux convectif, il
!       pourrait etre utile de vraiment utiliser la condition imposee
!       a la limite, i.e. de ne pas utiliser un schéma upwind (sinon,
!       on ne verra pas certaines C.L.). Actuellement, avec les conditions
!       aux limites proposées, il n'y a que 3 cas pour lesquels on
!       prend en compte le flux convectif ici : paroi, symetrie et
!       sortie supersonique. Dans les deux premiers, le flux est nul.
!       Pour le dernier, un traitement upwind correspond a utiliser
!       effectivement la valeur de bord. On peut donc conserver
!       le code tel quel.

!     On n'impose pas le flux convectif sur les faces pour lesquelles
!       il sera imposé par les C.L.
do ifac = 1, nfabor

  ii = ifabor(ifac)

  diipbx = diipb(1,ifac)
  diipby = diipb(2,ifac)
  diipbz = diipb(3,ifac)

  flui = 0.5d0*( flumab(ifac) +abs(flumab(ifac)) )
  fluj = 0.5d0*( flumab(ifac) -abs(flumab(ifac)) )

  pi = pvar(ii)

  pip = pi +ircflp*(grad(ii,1)*diipbx+grad(ii,2)*diipby+grad(ii,3)*diipbz)

  pfac  = inc*coefap(ifac) +coefbp(ifac)*pip
  pfacd = inc*cofafp(ifac) +cofbfp(ifac)*pip

!            FLUX = ICONVP*( FLUI*PVAR(II) +FLUJ*PFAC )
!     &           + IDIFFP*VISCB(IFAC)*( PIP -PFACD )
  flux = iconvp*( (flui*pi +fluj*pfac )*(1-ifbrus(ifac)) - flumab(ifac)*pi ) &
       + idiffp*viscb(ifac)*pfacd
  smbrp(ii) = smbrp(ii) -flux

enddo

! Free memory
deallocate(grad)
deallocate(dpdxa, dpdya, dpdza)

!--------
! FORMATS
!--------

 1000 format(1X,A8,' : CONVECTION EN ',A11,                             &
                               ' BLENDING A ',F4.0,' % D''UPWIND')
 1100 format(1X,A8,' : ',I10,' FACES UPWIND SUR ',                      &
                               I10,' FACES INTERNES ')
 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS cfbsc2                           ',/,&
'@    =========                                               ',/,&
'@     APPEL DE cfbsc2 POUR ',A8 ,' AVEC ISCHCP = ',I10        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return

end subroutine
