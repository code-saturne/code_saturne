!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine bilsc2 &
!================

 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  , coefap , coefbp , cofafp , cofbfp ,          &
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
! idtvar           ! e  ! <-- ! indicateur du schema temporel                  !
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
! relaxp           ! r  ! <-- ! coefficient de relaxation                      !
! thetap           ! r  ! <-- ! coefficient de ponderation pour le             !
!                  !    !     ! theta-schema (on ne l'utilise pour le          !
!                  !    !     ! moment que pour u,v,w et les scalaire          !
!                  !    !     ! - thetap = 0.5 correspond a un schema          !
!                  !    !     !   totalement centre en temps (mixage           !
!                  !    !     !   entre crank-nicolson et adams-               !
!                  !    !     !   bashforth)                                   !
! pvar (ncelet     ! tr ! <-- ! variable resolue (instant courant)             !
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
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          idtvar
integer          ivar   , iconvp , idiffp , nswrgp , imligp
integer          ircflp , ischcp , isstpp
integer          inc    , imrgra , iccocg
integer          iwarnp , ipp


double precision blencp , epsrgp , climgp, extrap, relaxp , thetap

double precision pvar (ncelet), pvara(ncelet)
double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscf (nfac), viscb (nfabor)
double precision smbrp(ncelet)

! Local variables

character*80     chaine
character*8      cnom
integer          ifac,ii,jj,infac,iel,iupwin, iij, iii, ig, it
integer          idimtr, irpvar
double precision pfac,pfacd,flui,fluj,flux,fluxi,fluxj
double precision difx,dify,difz,djfx,djfy,djfz
double precision pi, pj, pia, pja
double precision pif,pjf,pip,pjp,pir,pjr,pipr,pjpr
double precision pifri,pifrj,pjfri,pjfrj
double precision testi,testj,testij
double precision dpxf,dpyf,dpzf
double precision dcc, ddi, ddj, tesqck
double precision dijpfx, dijpfy, dijpfz
double precision diipfx, diipfy, diipfz
double precision djjpfx, djjpfy, djjpfz
double precision diipbx, diipby, diipbz
double precision pnd, distf, srfan
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
pifri = 0.d0
pifrj = 0.d0
pjfri = 0.d0
pjfrj = 0.d0

! Memoire


chaine = nomvar(ipp)
cnom   = chaine(1:8)

if (iwarnp.ge.2) then
  if (ischcp.eq.1) then
    write(nfecra,1000) cnom, '  CENTERED ', (1.d0-blencp)*100.d0
  else
    write(nfecra,1000) cnom, ' 2ND ORDER ', (1.d0-blencp)*100.d0
  endif
endif

iupwin = 0
if (blencp.eq.0.d0) iupwin = 1

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

if ((idiffp.ne.0 .and. ircflp.eq.1) .or.                          &
    (iconvp.ne.0 .and. iupwin.eq.0 .and.                          &
       (ischcp.eq.0 .or. ircflp.eq.1 .or. isstpp.eq.0))) then

  call grdcel                                                     &
  !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   pvar   , coefap , coefbp ,                                     &
   grad   )

else
  !$omp parallel do
  do iel = 1, ncelet
    grad(iel,1) = 0.d0
    grad(iel,2) = 0.d0
    grad(iel,3) = 0.d0
  enddo
endif


! ======================================================================
! ---> CALCUL DU GRADIENT DECENTRE DPDXA, DPDYA, DPDZA POUR TST DE PENTE
! ======================================================================

!$omp parallel do
do iel = 1, ncelet
  dpdxa(iel) = 0.d0
  dpdya(iel) = 0.d0
  dpdza(iel) = 0.d0
enddo

if (iconvp.gt.0.and.iupwin.eq.0.and.isstpp.eq.0) then

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, difx, dify, difz, djfx, djfy, djfz, &
    !$omp                     pif, pjf, pfac, pfac1, pfac2, pfac3)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)

        difx = cdgfac(1,ifac) - xyzcen(1,ii)
        dify = cdgfac(2,ifac) - xyzcen(2,ii)
        difz = cdgfac(3,ifac) - xyzcen(3,ii)
        djfx = cdgfac(1,ifac) - xyzcen(1,jj)
        djfy = cdgfac(2,ifac) - xyzcen(2,jj)
        djfz = cdgfac(3,ifac) - xyzcen(3,jj)

        pif = pvar(ii) + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
        pjf = pvar(jj) + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

        pfac = pjf
        if (flumas(ifac).gt.0.d0) pfac = pif

        pfac1 = pfac*surfac(1,ifac)
        pfac2 = pfac*surfac(2,ifac)
        pfac3 = pfac*surfac(3,ifac)

        dpdxa(ii) = dpdxa(ii) + pfac1
        dpdya(ii) = dpdya(ii) + pfac2
        dpdza(ii) = dpdza(ii) + pfac3

        dpdxa(jj) = dpdxa(jj) - pfac1
        dpdya(jj) = dpdya(jj) - pfac2
        dpdza(jj) = dpdza(jj) - pfac3

      enddo
    enddo
  enddo

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, diipbx, diipby, diipbz, pfac) &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)
        diipbx = diipb(1,ifac)
        diipby = diipb(2,ifac)
        diipbz = diipb(3,ifac)
        pfac =   inc*coefap(ifac)                                           &
               + coefbp(ifac) * (  pvar(ii)          + diipbx*grad(ii,1)    &
                                 + diipby*grad(ii,2) + diipbz*grad(ii,3))
        dpdxa(ii) = dpdxa(ii) + pfac*surfbo(1,ifac)
        dpdya(ii) = dpdya(ii) + pfac*surfbo(2,ifac)
        dpdza(ii) = dpdza(ii) + pfac*surfbo(3,ifac)

      enddo
    enddo
  enddo

  !$omp parallel do private(unsvol)
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

  if (iperot.eq.1) then
    call pergra(ivar, idimtr, irpvar)
    !==========
    if (idimtr .gt. 0) then
      call pering                                               &
      !==========
      ( idimtr , irpvar , iguper , igrper ,                     &
        dpdxa, dpdya, dpdza,                                    &
        dudxy  , drdxy )
    endif
  endif

endif


! ======================================================================
! ---> Contribution from interior faces
! ======================================================================

infac = 0

if (ncelet.gt.ncel) then
  !$omp parallel do if(ncelet - ncel > thr_n_min)
  do iel = ncel+1, ncelet
    smbrp(iel) = 0.d0
  enddo
endif

! --> Pure upwind flux
! =====================

if (iupwin.eq.1) then

  ! Steady
  if (idtvar.lt.0) then

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,      &
      !$omp                     diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz, &
      !$omp                     dpxf, dpyf, dpzf, pip, pjp, pipr, pjpr,         &
      !$omp                     flui, fluj, pif, pjf, fluxi, fluxj,             &
      !$omp                     pi, pj, pia, pja)                               &
      !$omp             reduction(+:infac)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)
          ! in parallel, face will be counted by one and only one rank
          if (ii.le.ncel) then
            infac = infac+1
          endif

          pi = pvar(ii)
          pj = pvar(jj)
          pia = pvara(ii)
          pja = pvara(jj)

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd = pond(ifac)

          ! Recompute II' and JJ' at this level

          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          ! reconstruction only if IRCFLP = 1
          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          pir = pi/relaxp - (1.d0-relaxp)/relaxp * pia
          pjr = pj/relaxp - (1.d0-relaxp)/relaxp * pja

          pipr = pir                                                   &
               + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjpr = pjr                                                   &
               + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) +abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))

          pif  = pi
          pjf  = pj

          fluxi = iconvp*(flui*pir + fluj*pjf - flumas(ifac)*pi)       &
                + idiffp*viscf(ifac)*(pipr - pjp)
          fluxj = iconvp*(flui*pif + fluj*pjr - flumas(ifac)*pj)       &
                + idiffp*viscf(ifac)*(pip - pjpr)

          smbrp(ii) = smbrp(ii) - fluxi
          smbrp(jj) = smbrp(jj) + fluxj

        enddo
      enddo
    enddo

  ! Unsteady
  else

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,      &
      !$omp                     diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz, &
      !$omp                     dpxf, dpyf, dpzf, pip, pjp, flui, fluj,         &
      !$omp                     pif, pjf, flux, pi, pj)                         &
      !$omp             reduction(+:infac)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)
          ! in parallel, face will be counted by one and only one rank
          if (ii.le.ncel) then
            infac = infac+1
          endif

          pi = pvar(ii)
          pj = pvar(jj)

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd = pond(ifac)

          ! Recompute II' and JJ' at this level

          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) + abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) - abs(flumas(ifac)))

          pif = pi
          pjf = pj

          flux = iconvp*(flui*pif +fluj*pjf) + idiffp*viscf(ifac)*(pip -pjp)

          smbrp(ii) = smbrp(ii) - thetap *(flux - iconvp*flumas(ifac)*pi)
          smbrp(jj) = smbrp(jj) + thetap *(flux - iconvp*flumas(ifac)*pj)

        enddo
      enddo
    enddo

  endif


! --> Flux with no slope test
! ============================

else if (isstpp.eq.1) then

  if (ischcp.lt.0 .or. ischcp.gt.1) then
    write(nfecra,9000) ischcp
    call csexit(1)
  endif

  ! Steady
  if (idtvar.lt.0) then

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,      &
      !$omp                     diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz, &
      !$omp                     dpxf, dpyf, dpzf, pip, pjp, pipr, pjpr, flui,   &
      !$omp                     fluj, pir, pjr, pifri, pjfri, pifrj, pjfrj,     &
      !$omp                     difx, dify, difz, djfx, djfy, djfz,             &
      !$omp                     fluxi, fluxj, pi, pj, pia, pja)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd = pond(ifac)

          pi = pvar(ii)
          pj = pvar(jj)
          pia = pvara(ii)
          pja = pvara(jj)

          ! Recompute II' and JJ' at this level

          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          pir = pi/relaxp - (1.d0 - relaxp)/relaxp*pia
          pjr = pj/relaxp - (1.d0 - relaxp)/relaxp*pja

          pipr = pir                                                    &
               + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjpr = pjr                                                    &
               + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) + abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) - abs(flumas(ifac)))


          ! Centered
          ! --------

          if (ischcp.eq.1) then

            pifri = pnd*pipr +(1.d0-pnd)*pjp
            pjfri = pifri
            pifrj = pnd*pip  +(1.d0-pnd)*pjpr
            pjfrj = pifrj

          ! Second order
          ! ------------

          else ! if (ischcp.eq.0) then

            difx = cdgfac(1,ifac) - xyzcen(1,ii)
            dify = cdgfac(2,ifac) - xyzcen(2,ii)
            difz = cdgfac(3,ifac) - xyzcen(3,ii)
            djfx = cdgfac(1,ifac) - xyzcen(1,jj)
            djfy = cdgfac(2,ifac) - xyzcen(2,jj)
            djfz = cdgfac(3,ifac) - xyzcen(3,jj)

            ! leave reconstruction of PIF and PJF even if IRCFLP=0
            ! otherwise, it is the same as using upwind
            pifri = pir + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
            pifrj = pi + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
            pjfrj = pjr + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)
            pjfri = pj + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

          endif

          ! Blending
          ! --------

          pifri = blencp*pifri+(1.d0-blencp)*pir
          pifrj = blencp*pifrj+(1.d0-blencp)*pi
          pjfri = blencp*pjfri+(1.d0-blencp)*pj
          pjfrj = blencp*pjfrj+(1.d0-blencp)*pjr

          ! Flux
          ! ----

          fluxi =   iconvp*(flui*pifri + fluj*pjfri - flumas(ifac)*pi) &
                  + idiffp*viscf(ifac)*(pipr -pjp)
          fluxj =   iconvp*(flui*pifrj + fluj*pjfrj - flumas(ifac)*pj) &
                  + idiffp*viscf(ifac)*(pip -pjpr)

          ! Assembly
          ! --------

          smbrp(ii) = smbrp(ii) - fluxi
          smbrp(jj) = smbrp(jj) + fluxj

        enddo
      enddo
    enddo

  ! Unsteady
  else

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,      &
      !$omp                     diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz, &
      !$omp                     dpxf, dpyf, dpzf, pip, pjp, flui, fluj, pif,    &
      !$omp                     pjf, difx, dify, difz, djfx, djfy, djfz, flux   &
      !$omp                     pi, pj)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd = pond(ifac)

          pi = pvar(ii)
          pj = pvar(jj)

          ! Recompute II' and JJ' at this level

          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) + abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) - abs(flumas(ifac)))

          ! Centered
          ! --------

          if (ischcp.eq.1) then

            pif = pnd*pip +(1.d0-pnd)*pjp
            pjf = pif

          ! Second order
          ! ------------

          else ! if (ischcp.eq.0) then

            difx = cdgfac(1,ifac) - xyzcen(1,ii)
            dify = cdgfac(2,ifac) - xyzcen(2,ii)
            difz = cdgfac(3,ifac) - xyzcen(3,ii)
            djfx = cdgfac(1,ifac) - xyzcen(1,jj)
            djfy = cdgfac(2,ifac) - xyzcen(2,jj)
            djfz = cdgfac(3,ifac) - xyzcen(3,jj)

            ! leave reconstruction of PIF and PJF even if IRCFLP=0
            ! otherwise, it is the same as using upwind
            pif = pi + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
            pjf = pj + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

          endif

          ! Blending
          ! --------

          pif = blencp*pif+(1.d0-blencp)*pi
          pjf = blencp*pjf+(1.d0-blencp)*pj

          ! Flux
          ! ----

          flux = iconvp*(flui*pif +fluj*pjf) + idiffp*viscf(ifac)*(pip -pjp)

          ! Assembly
          ! --------

          smbrp(ii) = smbrp(ii) - thetap *(flux - iconvp*flumas(ifac)*pi)
          smbrp(jj) = smbrp(jj) + thetap *(flux - iconvp*flumas(ifac)*pj)

        enddo
      enddo
    enddo

  endif

! --> Flux with slope test
! =========================

else

  if (ischcp.lt.0 .or. ischcp.gt.1) then
    write(nfecra,9000) ischcp
    call csexit(1)
  endif

  ! Steady
  if (idtvar.lt.0) then

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,      &
      !$omp                     distf, srfan, diipfx, diipfy, diipfz, djjpfx,   &
      !$omp                     djjpfy, djjpfz, dpxf, dpyf, dpzf, pip, pjp,     &
      !$omp                     pipr, pjpr, flui, fluj, pir, pjr, testi, testj, &
      !$omp                     testij, dcc, ddi, ddj, tesqck, pifri, pjfri,    &
      !$omp                     pifrj, pjfrj, difx, dify, difz, djfx, djfy,     &
      !$omp                     djfz, fluxi, fluxj, pi, pj, pia, pja)           &
      !$omp             reduction(+:infac)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd    = pond(ifac)
          distf  = dist(ifac)
          srfan  = surfan(ifac)

          pi = pvar(ii)
          pj = pvar(jj)
          pia = pvara(ii)
          pja = pvara(jj)

          ! Recompute II' and JJ' at this level
          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          pir = pi/relaxp - (1.d0 - relaxp)/relaxp*pia
          pjr = pj/relaxp - (1.d0 - relaxp)/relaxp*pja

          pipr = pir                                                  &
               + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjpr = pjr                                                  &
               + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) +abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))


          ! Slope test
          ! ----------

          testi =   dpdxa(ii)*surfac(1,ifac) + dpdya(ii)*surfac(2,ifac)  &
                  + dpdza(ii)*surfac(3,ifac)
          testj =   dpdxa(jj)*surfac(1,ifac) + dpdya(jj)*surfac(2,ifac)  &
                  + dpdza(jj)*surfac(3,ifac)
          testij=   dpdxa(ii)*dpdxa(jj)      + dpdya(ii)*dpdya(jj)       &
                  + dpdza(ii)*dpdza(jj)

          if (flumas(ifac).gt.0.d0) then
            dcc =   grad(ii,1)*surfac(1,ifac) +grad(ii,2)*surfac(2,ifac)    &
                  + grad(ii,3)*surfac(3,ifac)
            ddi = testi
            ddj = (pj-pi)/distf *srfan
          else
            dcc =   grad(jj,1)*surfac(1,ifac) +grad(jj,2)*surfac(2,ifac)    &
                  + grad(jj,3)*surfac(3,ifac)
            ddi = (pj-pi)/distf *srfan
            ddj = testj
          endif
          tesqck = dcc**2 -(ddi-ddj)**2

          ! Upwind
          ! ------

          if (tesqck.le.0.d0 .or. testij.le.0.d0) then

            pifri = pir
            pifrj = pi
            pjfri = pj
            pjfrj = pjr
            ! in parallel, face will be counted by one and only one rank
            if (ii.le.ncel) then
              infac = infac+1
            endif

          else

            ! Centered
            ! --------

            if (ischcp.eq.1) then

              pifri = pnd*pipr +(1.d0-pnd)*pjp
              pjfri = pifri
              pifrj = pnd*pip  +(1.d0-pnd)*pjpr
              pjfrj = pifrj

            ! Second order
            ! ------------

            else ! if (ischcp.eq.0) then

              difx = cdgfac(1,ifac) - xyzcen(1,ii)
              dify = cdgfac(2,ifac) - xyzcen(2,ii)
              difz = cdgfac(3,ifac) - xyzcen(3,ii)
              djfx = cdgfac(1,ifac) - xyzcen(1,jj)
              djfy = cdgfac(2,ifac) - xyzcen(2,jj)
              djfz = cdgfac(3,ifac) - xyzcen(3,jj)

              ! leave reconstruction of PIF and PJF even if IRCFLP=0
              ! otherwise, it is the same as using upwind
              pifri = pir + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
              pifrj = pi + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
              pjfrj = pjr + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)
              pjfri = pj + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

            endif

          endif

          ! Blending
          ! --------

          pifri = blencp*pifri+(1.d0-blencp)*pir
          pifrj = blencp*pifrj+(1.d0-blencp)*pi
          pjfri = blencp*pjfri+(1.d0-blencp)*pj
          pjfrj = blencp*pjfrj+(1.d0-blencp)*pjr

          ! Flux
          ! ----

          fluxi =   iconvp*(flui*pifri + fluj*pjfri - flumas(ifac)*pi) &
                  + idiffp*viscf(ifac)*(pipr -pjp)
          fluxj =   iconvp*(flui*pifrj + fluj*pjfrj - flumas(ifac)*pj) &
                  + idiffp*viscf(ifac)*(pip -pjpr)

          ! Assembly
          ! --------

          smbrp(ii) = smbrp(ii) - fluxi
          smbrp(jj) = smbrp(jj) + fluxj

        enddo
      enddo
    enddo

  ! Unsteady
  else

    do ig = 1, ngrpi
       !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,     &
      !$omp                     distf, srfan, diipfx, diipfy, diipfz, djjpfx,   &
      !$omp                     djjpfy, djjpfz, dpxf, dpyf, dpzf, pip, pjp,     &
      !$omp                     flui, fluj, testi, testj, testij, dcc, ddi,     &
      !$omp                     ddj, tesqck, pif, pjf, difx, dify, difz,        &
      !$omp                     djfx, djfy, djfz, flux, pi, pj)                 &
      !$omp             reduction(+:infac)
       do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd    = pond(ifac)
          distf  = dist(ifac)
          srfan  = surfan(ifac)

          pi = pvar(ii)
          pj = pvar(jj)

          ! Recompute II' and JJ' at this level

          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) +abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))

          ! Slope test
          ! ----------

          testi =   dpdxa(ii)*surfac(1,ifac) + dpdya(ii)*surfac(2,ifac)  &
                  + dpdza(ii)*surfac(3,ifac)
          testj =   dpdxa(jj)*surfac(1,ifac) + dpdya(jj)*surfac(2,ifac)  &
                  + dpdza(jj)*surfac(3,ifac)
          testij =   dpdxa(ii)*dpdxa(jj)    + dpdya(ii)*dpdya(jj)        &
                   + dpdza(ii)*dpdza(jj)

          if (flumas(ifac).gt.0.d0) then
            dcc =   grad(ii,1)*surfac(1,ifac) + grad(ii,2)*surfac(2,ifac)    &
                  + grad(ii,3)*surfac(3,ifac)
            ddi = testi
            ddj = (pj-pi)/distf *srfan
          else
            dcc =   grad(jj,1)*surfac(1,ifac) + grad(jj,2)*surfac(2,ifac)    &
                  + grad(jj,3)*surfac(3,ifac)
            ddi = (pj-pi)/distf *srfan
            ddj = testj
          endif
          tesqck = dcc**2 -(ddi-ddj)**2

          ! Upwind
          ! ------

          if (tesqck.le.0.d0 .or. testij.le.0.d0) then

            pif = pi
            pjf = pj
            ! in parallel, face will be counted by one and only one rank
            if (ii.le.ncel) then
              infac = infac+1
            endif

          else

            ! Centered
            ! --------

            if (ischcp.eq.1) then

              pif = pnd*pip +(1.d0-pnd)*pjp
              pjf = pif

            ! Second order
            ! ------------

            else ! if (ischcp.eq.0) then

              difx = cdgfac(1,ifac) - xyzcen(1,ii)
              dify = cdgfac(2,ifac) - xyzcen(2,ii)
              difz = cdgfac(3,ifac) - xyzcen(3,ii)
              djfx = cdgfac(1,ifac) - xyzcen(1,jj)
              djfy = cdgfac(2,ifac) - xyzcen(2,jj)
              djfz = cdgfac(3,ifac) - xyzcen(3,jj)

              ! leave reconstruction of PIF and PJF even if IRCFLP=0
              ! otherwise, it is the same as using upwind
              pif = pi + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
              pjf = pj + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

            endif

          endif

          ! Blending
          ! --------

          pif = blencp*pif+(1.d0-blencp)*pi
          pjf = blencp*pjf+(1.d0-blencp)*pj

          ! Flux
          ! ----

          flux =   iconvp*(flui*pif +fluj*pjf)                        &
                 + idiffp*viscf(ifac)*(pip -pjp)

          ! Assembly
          ! --------

          smbrp(ii) = smbrp(ii) - thetap *(flux - iconvp*flumas(ifac)*pi)
          smbrp(jj) = smbrp(jj) + thetap *(flux - iconvp*flumas(ifac)*pj)

        enddo
      enddo
    enddo

  endif ! idtvar

endif ! iupwin


if (iwarnp.ge.2) then
  if (irangp.ge.0) call parcpt(infac)
  write(nfecra,1100) cnom, infac, nfacgb
endif

! ======================================================================
! ---> Contribution from boundary faces
! ======================================================================

! Steady
if (idtvar.lt.0) then

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, diipbx, diipby, diipbz, flui, fluj,     &
    !$omp                     pir, pipr, pfac, pfacd, flux, pi, pia)            &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        pi = pvar(ii)
        pia = pvara(ii)

        diipbx = diipb(1,ifac)
        diipby = diipb(2,ifac)
        diipbz = diipb(3,ifac)

        ! Remove decentering for coupled faces
        if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
          flui = 0.0d0
          fluj = flumab(ifac)
        else
          flui = 0.5d0*(flumab(ifac) +abs(flumab(ifac)))
          fluj = 0.5d0*(flumab(ifac) -abs(flumab(ifac)))
        endif

        pir  = pi/relaxp - (1.d0-relaxp)/relaxp*pia
        pipr =   pir                                                            &
               + ircflp*(grad(ii,1)*diipbx+grad(ii,2)*diipby+grad(ii,3)*diipbz)

        pfac  = inc*coefap(ifac) +coefbp(ifac)*pipr
        pfacd = inc*cofafp(ifac) +cofbfp(ifac)*pipr

        flux =   iconvp*(flui*pir + fluj*pfac - flumab(ifac)*pi )               &
               + idiffp*viscb(ifac)*(pipr -pfacd)
        smbrp(ii) = smbrp(ii) - flux

      enddo
    enddo
  enddo

! Unsteady
else

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, diipbx, diipby, diipbz, flui, fluj,     &
    !$omp                     pip, pfac, pfacd, flux, pi)                       &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        pi = pvar(ii)

        diipbx = diipb(1,ifac)
        diipby = diipb(2,ifac)
        diipbz = diipb(3,ifac)

        ! Remove decentering for coupled faces
        if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
          flui = 0.0d0
          fluj = flumab(ifac)
        else
          flui = 0.5d0*(flumab(ifac) +abs(flumab(ifac)))
          fluj = 0.5d0*(flumab(ifac) -abs(flumab(ifac)))
        endif

        pip = pi                                                       &
            + ircflp*(grad(ii,1)*diipbx+grad(ii,2)*diipby+grad(ii,3)*diipbz)

        pfac  = inc*coefap(ifac) + coefbp(ifac)*pip
        pfacd = inc*cofafp(ifac) + cofbfp(ifac)*pip

        flux =   iconvp*((flui - flumab(ifac))*pi + fluj*pfac)         &
               + idiffp*viscb(ifac)*(pip -pfacd)
        smbrp(ii) = smbrp(ii) - thetap * flux

      enddo
    enddo
  enddo

endif

! Free memory
deallocate(grad)
deallocate(dpdxa, dpdya, dpdza)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(1X,A8,' : CONVECTION EN ',A11,                             &
                               ' BLENDING A ',F4.0,' % D''UPWIND')
 1100 format(1X,A8,' : ',I10,' FACES UPWIND SUR ',                      &
                               I10,' FACES INTERNES ')
 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS bilsc2                           ',/,&
'@    =========                                               ',/,&
'@     APPEL DE bilsc2 POUR ',A8 ,' AVEC ISCHCP = ',I10        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(1X,A8,' : CONVECTION IN ',A11,                             &
                            ' BLENDING WITH ',F4.0,' % OF UPWIND')
 1100 format(1X,A8,' : ',I10,' FACES WITH UPWIND ON ',                  &
                               I10,' INTERIOR FACES ')
 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN bilsc2                                ',/,&
'@    ========                                                ',/,&
'@     CALL OF bilsc2 FOR ',A8 ,' WITH ISCHCP = ',I10          ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Contact the support.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! End
!----

return

end subroutine
