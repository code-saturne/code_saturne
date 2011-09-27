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

subroutine bilsc4 &
!================

 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg , ivisep ,          &
   ippu   , ippv   , ippw   , iwarnp ,                            &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   vel    , vela   ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscf  , viscb  , secvif , secvib ,          &
   smbr   )

!===============================================================================
! FONCTION :
! ---------

! CALCUL DU BILAN EXPLICITE DE LA VARIABLE VARx VARy VARz (VITESSE)
!                                              --->
! --->   --->      ___ .   --->                --->     -->
! SMBR = SMBR -    \   m   VAR  +Visc   ( grad VAR )  . n
!                  /__  ij    ij     ij             ij    ij
!                   j

! REMARK:
! if ivisep = 1 then we also take MU*TRANSPOSE_GRAD(VAR)
!                               + LAMBDA*TRACE(GRAD(VAR))
!
! where lambda is the secondary viscosity, i.e. usually -2/3*MU

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
! ivisep           ! e  ! <-- ! indicateur = 1 pour la prise en compte         !
!                  !    !     !                div(T Grad(vel))                !
!                  !    !     !                -2/3 Grad(div(vel))             !
!                  !    !     !              0 sinon                           !
! ipp*             ! e  ! <-- ! numero de variable pour post                   !
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
! vel (3,ncelet)   ! tr ! <-- ! vitesse resolue (instant courant)              !
! vela(3,ncelet)   ! tr ! <-- ! vitesse resolue (instant precedent)            !
! coefav           ! tr ! <-- ! tableaux des cond lim pour u, v, w             !
!   (3,nfabor)     !    !     !  sur la normale a la face de bord              !
! coefbv           ! tr ! <-- ! tableaux des cond lim pour u, v, w             !
!   (3,3,nfabor)   !    !     !  sur la normale a la face de bord              !
! cofafv           ! tr ! <-- ! tableaux des cond lim pour le flux de          !
!   (3,nfabor)     !    !     !  diffusion de u, v, w                          !
! cofbfv           ! tr ! <-- ! tableaux des cond lim pour le flux de          !
!   (3,3,nfabor)   !    !     !  diffusion de u, v, w                          !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! viscf (nfac)     ! tr ! <-- ! visc*surface/dist aux faces internes           !
!                  !    !     !  pour second membre                            !
! viscb (nfabor    ! tr ! <-- ! visc*surface/dist aux faces de bord            !
!                  !    !     !  pour second membre                            !
! secvif(nfac)     ! tr ! --- ! secondary viscosity at interior faces          !
! secvib(nfabor)   ! tr ! --- ! secondary viscosity at boundary faces          !
! smbr(3,ncelet)   ! tr ! <-- ! bilan au second membre                         !
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
integer          inc    , imrgra , iccocg , ivisep
integer          iwarnp , ippu   , ippv   , ippw

double precision blencp , epsrgp , climgp, extrap, relaxp , thetap
double precision vel   (3  ,ncelet)
double precision vela  (3  ,ncelet)
double precision coefav(3  ,nfabor)
double precision cofafv(3  ,nfabor)
double precision coefbv(3,3,nfabor)
double precision cofbfv(3,3,nfabor)
double precision flumas(nfac)  , flumab(nfabor)
double precision viscf (nfac)  , viscb (nfabor)
double precision secvif(nfac), secvib(nfabor)
double precision smbr(3,ncelet)


! Local variables

character*80     chaine
character*8      cnom
integer          ifac,ii,jj,infac,iel,iupwin, iok
integer          iiu,iiv,iiw
integer          iitytu
integer          iir11,iir22,iir33
integer          iir12,iir13,iir23
integer          isou, jsou, ityp
logical          ilved
double precision pfac,pfacd,flui,fluj,flux,fluxi,fluxj
double precision vfac(3)
double precision difv(3), djfv(3)
double precision pi , pj, pia, pja
double precision pif,pjf,pip,pjp,pir,pjr,pipr,pjpr
double precision pifr,pjfr,pifri,pifrj,pjfri,pjfrj
double precision testi,testj,testij
double precision dpvf(3)
double precision dcc, ddi, ddj, tesqck
double precision dijpfv(3)
double precision diipfv(3)
double precision djjpfv(3)
double precision diipbv(3)
double precision pnd, distf, srfan
double precision unsvol, visco, grdtrv, tgrdfl, secvis

double precision, dimension(:,:,:), allocatable :: gradv, gradva
double precision, dimension(:), allocatable :: bndcel


!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(gradv(3,3,ncelet))
allocate(gradva(3,3,ncelet))

! Initialize variables to avoid compiler warnings

pif = 0.d0
pjf = 0.d0
pifri = 0.d0
pifrj = 0.d0
pjfri = 0.d0
pjfrj = 0.d0

pi  = 0.d0
pj  = 0.d0
pia = 0.d0
pja = 0.d0

! Memoire

chaine = nomvar(ippu)
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
! ---> CALCUL DU GRADIENT DE VITESSE
! ======================================================================
!    DUDX sert a la fois pour la reconstruction des flux et pour le test
!    de pente. On doit donc le calculer :
!        - quand on a de la diffusion et qu'on reconstruit les flux
!        - quand on a de la convection SOLU
!        - quand on a de la convection, qu'on n'est pas en upwind pur
!          et qu'on reconstruit les flux
!        - quand on a de la convection, qu'on n'est pas en upwind pur
!          et qu'on n'a pas shunte le test de pente

if( (idiffp.ne.0 .and. ircflp.eq.1) .or. ivisep.eq.1 .or.         &
    (iconvp.ne.0 .and. iupwin.eq.0 .and.                          &
    (ischcp.eq.0 .or.  ircflp.eq.1 .or. isstpp.eq.0)) ) then


  ilved = .true.
  iccocg = 0

  call grdvec &
  !==========
( iu     , imrgra , inc    , iccocg ,nswrgp , imligp ,           &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
  ilved  ,                                                       &
  vel    , coefav , coefbv ,                                     &
  gradv )

else
  do iel = 1, ncelet
    do isou =1, 3
      do jsou = 1, 3
        gradv(isou,jsou,iel) = 0.d0
      enddo
    enddo
  enddo
endif


! ======================================================================
! ---> CALCUL DU GRADIENT DECENTRE DPDXA, DPDYA, DPDZA POUR TST DE PENTE
! ======================================================================

do iel = 1, ncelet
  do jsou = 1, 3
    do isou =1, 3
      gradva(isou,jsou,iel) = 0.d0
    enddo
  enddo
enddo

if( iconvp.gt.0.and.iupwin.eq.0.and.isstpp.eq.0 ) then

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    do jsou = 1, 3
      difv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,ii)
      djfv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,jj)
    enddo
    !-----------------
    ! X-Y-Z component, p=u, v, w
    do isou = 1, 3
      pif = vel(isou,ii)
      pjf = vel(isou,jj)
      do jsou = 1, 3
        pif = pif + gradv(isou,jsou,ii)*difv(jsou)
        pjf = pjf + gradv(isou,jsou,jj)*djfv(jsou)
      enddo

      pfac = pjf
      if( flumas(ifac ).gt.0.d0 ) pfac = pif

      ! U gradient
      do jsou = 1, 3
        vfac(jsou) = pfac*surfac(jsou,ifac)

        gradva(isou,jsou,ii) = gradva(isou,jsou,ii) +vfac(jsou)
        gradva(isou,jsou,jj) = gradva(isou,jsou,jj) -vfac(jsou)
      enddo
    enddo

  enddo

  do ifac = 1, nfabor
    ii = ifabor(ifac )

    do jsou = 1, 3
      diipbv(jsou) = diipb(jsou,ifac)
    enddo
    !-----------------
    ! X-Y-Z components, p=u, v, w
    do isou = 1,3
      pfac = inc*coefav(isou,ifac)
      !coefu is a matrix
      do jsou =  1, 3
        pfac = pfac + coefbv(isou,jsou,ifac)*(   vel(jsou,ii)    &
                    + gradv(jsou,1,ii)*diipbv(1)               &
                    + gradv(jsou,2,ii)*diipbv(2)               &
                    + gradv(jsou,3,ii)*diipbv(3)     )
      enddo

      do jsou = 1, 3
        gradva(isou,jsou,ii) = gradva(isou,jsou,ii) +pfac*surfbo(jsou,ifac )
      enddo
    enddo

  enddo


  do iel = 1, ncel
    unsvol = 1.d0/volume(iel)
    do isou = 1, 3
      do jsou = 1, 3
        gradva(isou,jsou,iel) = gradva(isou,jsou,iel)*unsvol
      enddo
    enddo
  enddo

  !     TRAITEMENT DU PARALLELISME, ET DE LA PERIODICITE

  if(irangp.ge.0.or.iperio.eq.1) then
    call syntin (gradva)
    !==========

  endif

endif


! ======================================================================
! ---> ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
! ======================================================================

infac = 0

if(ncelet.gt.ncel) then
  do iel = ncel+1, ncelet
    do isou = 1, 3
      smbr(isou,iel) = 0.d0
    enddo
  enddo
endif


!  --> FLUX UPWIND PUR
!  =====================

if(iupwin.eq.1) then

!     Stationnaire
  if (idtvar.lt.0) then

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
!     en parallele, la face sera comptee d'un cote OU (exclusif) de l'autre
      if (ii.le.ncel) then
        infac = infac+1
      endif

      do jsou = 1, 3
        dijpfv(jsou) = dijpf(jsou,ifac)
      enddo

      pnd   = pond(ifac)

! ON RECALCULE A CE NIVEAU II' ET JJ'
      do jsou = 1, 3
        diipfv(jsou) = cdgfac(jsou,ifac) - (xyzcen(jsou,ii)+                  &
               (1.d0-pnd) * dijpfv(jsou))
        djjpfv(jsou) = cdgfac(jsou,ifac) -  xyzcen(jsou,jj)+                  &
                   pnd  * dijpfv(jsou)
      enddo

      flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
      fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )

!-----------------
! X-Y-Z components, p=u, v, w
      do isou = 1, 3
        do jsou = 1, 3
          dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
        enddo

!     reconstruction uniquement si IRCFLP = 1
        pi  = vel (isou,ii)
        pj  = vel (isou,jj)

        pia = vela(isou,ii)
        pja = vela(isou,jj)

        pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                          +dpvf(2)*diipfv(2)        &
                          +dpvf(3)*diipfv(3))
        pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                          +dpvf(2)*djjpfv(2)        &
                          +dpvf(3)*djjpfv(3))

        pipr = pi /relaxp - (1.d0-relaxp)/relaxp * pia   &
             + ircflp*(dpvf(1)*diipfv(1)                 &
                      +dpvf(2)*diipfv(2)                 &
                      +dpvf(3)*diipfv(3) )
        pjpr = pj /relaxp - (1.d0-relaxp)/relaxp * pja   &
             + ircflp*(dpvf(1)*djjpfv(1)                 &
                      +dpvf(2)*djjpfv(2)                 &
                      +dpvf(3)*djjpfv(3))

        pifr = pi /relaxp - (1.d0-relaxp)/relaxp * pia
        pjfr = pj /relaxp - (1.d0-relaxp)/relaxp * pja

        fluxi = iconvp*( flui*pifr +fluj*pj )                      &
              + idiffp*viscf(ifac)*( pipr -pjp )
        fluxj = iconvp*( flui*pi +fluj*pjfr )                      &
              + idiffp*viscf(ifac)*( pip -pjpr )

        smbr(isou,ii) = smbr(isou,ii) - fluxi
        smbr(isou,jj) = smbr(isou,jj) + fluxj

      enddo

    enddo

!     Instationnaire
  else

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
!     en parallele, la face sera comptee d'un cote OU (exclusif) de l'autre
      if (ii.le.ncel) then
        infac = infac+1
      endif

      do jsou = 1, 3
        dijpfv(jsou) = dijpf(jsou,ifac)
      enddo

      pnd   = pond(ifac)

! ON RECALCULE A CE NIVEAU II' ET JJ'
      do jsou = 1, 3
        diipfv(jsou) = cdgfac(jsou,ifac) - (xyzcen(jsou,ii)+                  &
               (1.d0-pnd) * dijpfv(jsou))
        djjpfv(jsou) = cdgfac(jsou,ifac) -  xyzcen(jsou,jj)+                  &
                   pnd  * dijpfv(jsou)
      enddo

      flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
      fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )

!-----------------
! X-Y-Z components, p=u, v, w
      do isou = 1, 3

        do jsou = 1, 3
          dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
        enddo

        pi = vel(isou,ii)
        pj = vel(isou,jj)

        pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                          +dpvf(2)*diipfv(2)        &
                          +dpvf(3)*diipfv(3))
        pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                          +dpvf(2)*djjpfv(2)        &
                          +dpvf(3)*djjpfv(3))

        flux = iconvp*( flui*pi +fluj*pj )                        &
             + idiffp*viscf(ifac)*( pip -pjp )

        smbr(isou,ii) = smbr(isou,ii) - thetap * flux
        smbr(isou,jj) = smbr(isou,jj) + thetap * flux

      enddo

    enddo

  endif


!  --> FLUX SANS TEST DE PENTE
!  ============================

elseif(isstpp.eq.1) then

!     Stationnaire
  if (idtvar.lt.0) then

    iok = 0
    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      do jsou = 1, 3
        dijpfv(jsou) = dijpf(jsou,ifac)
      enddo

      pnd   = pond(ifac)

! ON RECALCULE A CE NIVEAU II' ET JJ'
      do jsou = 1, 3
        diipfv(jsou) = cdgfac(jsou,ifac) - (xyzcen(jsou,ii)+                  &
               (1.d0-pnd) * dijpfv(jsou))
        djjpfv(jsou) = cdgfac(jsou,ifac) -  xyzcen(jsou,jj)+                  &
                   pnd  * dijpfv(jsou)
      enddo

      flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
      fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )

!  Pour le second ordre on definit IF
      if(ischcp.eq.0) then
        do jsou = 1, 3
          difv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,ii)
          djfv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,jj)
        enddo
      endif


!-----------------
! X-Y-Z components, p=u, v, w
      do isou = 1, 3

        do jsou = 1, 3
          dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
        enddo

        pi = vel (isou,ii)
        pj = vel (isou,jj)

        pia = vela(isou,ii)
        pja = vela(isou,jj)

        pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                          +dpvf(2)*diipfv(2)        &
                          +dpvf(3)*diipfv(3))
        pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                          +dpvf(2)*djjpfv(2)        &
                          +dpvf(3)*djjpfv(3))

        pipr = pi /relaxp - (1.d0-relaxp)/relaxp * pia   &
             + ircflp*(dpvf(1)*diipfv(1)                 &
                      +dpvf(2)*diipfv(2)                 &
                      +dpvf(3)*diipfv(3) )
        pjpr = pj /relaxp - (1.d0-relaxp)/relaxp * pja   &
             + ircflp*(dpvf(1)*djjpfv(1)                 &
                      +dpvf(2)*djjpfv(2)                 &
                      +dpvf(3)*djjpfv(3))

        pir = pi /relaxp - (1.d0 - relaxp)/relaxp* pia
        pjr = pj /relaxp - (1.d0 - relaxp)/relaxp* pja

!         CENTRE
!        --------

        if (ischcp.eq.1) then

          pifri = pnd*pipr +(1.d0-pnd)*pjp
          pjfri = pifri
          pifrj = pnd*pip  +(1.d0-pnd)*pjpr
          pjfrj = pifrj


!         SECOND ORDER
!        --------------

        elseif(ischcp.eq.0) then
! dif* is already defined

!     on laisse la reconstruction de PIF et PJF meme si IRCFLP=0
!     sinon cela revient a faire de l'upwind
          pifri = pir + difv(1)*gradv(isou,1,ii)      &
                      + difv(2)*gradv(isou,2,ii)      &
                      + difv(3)*gradv(isou,3,ii)
          pifrj = pi  + difv(1)*gradv(isou,1,ii)      &
                      + difv(2)*gradv(isou,2,ii)      &
                      + difv(3)*gradv(isou,3,ii)

          pjfrj = pjr + djfv(1)*gradv(isou,1,jj)      &
                      + djfv(2)*gradv(isou,2,jj)      &
                      + djfv(3)*gradv(isou,3,jj)
          pjfri = pj  + djfv(1)*gradv(isou,1,jj)      &
                      + djfv(2)*gradv(isou,2,jj)      &
                      + djfv(3)*gradv(isou,3,jj)
        else
          write(nfecra,9000)ischcp
          iok = 1
        endif


!        BLENDING
!       ----------

        pifri = blencp*pifri+(1.d0-blencp)*pir
        pifrj = blencp*pifrj+(1.d0-blencp)*pif
        pjfri = blencp*pjfri+(1.d0-blencp)*pjf
        pjfrj = blencp*pjfrj+(1.d0-blencp)*pjr


!        FLUX
!       ------

        fluxi = iconvp*( flui*pifri + fluj*pjfri )                  &
               +idiffp*viscf(ifac)*( pipr -pjp )
        fluxj = iconvp*( flui*pifrj +fluj*pjfrj )                   &
               +idiffp*viscf(ifac)*( pip -pjpr )


!        ASSEMBLAGE
!       ------------

        smbr(isou,ii) = smbr(isou,ii) - fluxi
        smbr(isou,jj) = smbr(isou,jj) + fluxj
      enddo
      !end isou

    enddo
!        Position "hors boucle" du CALL CSEXIT pour raisons
!        historiques,pour eviter de devectoriser la boucle
!        -> a conserver si on recherche a optimiser la boucle en
!        vectorisation
    if(iok.ne.0) then
      call csexit (1)
    endif

!     Instationnaire
  else

    iok = 0
    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      do jsou = 1, 3
        dijpfv(jsou) = dijpf(jsou,ifac)
      enddo

      pnd   = pond(ifac)

! ON RECALCULE A CE NIVEAU II' ET JJ'
      do jsou = 1, 3
        diipfv(jsou) = cdgfac(jsou,ifac) - (xyzcen(jsou,ii)+                  &
               (1.d0-pnd) * dijpfv(jsou))
        djjpfv(jsou) = cdgfac(jsou,ifac) -  xyzcen(jsou,jj)+                  &
                   pnd  * dijpfv(jsou)
      enddo

      flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
      fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )

!  Pour le second ordre on definit IF
      if(ischcp.eq.0) then
        do jsou = 1, 3
          difv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,ii)
          djfv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,jj)
        enddo
      endif

!-----------------
! X-Y-Z components, p=u, v, w
      do isou = 1, 3

        do jsou = 1, 3
          dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
        enddo

        pi = vel(isou,ii)
        pj = vel(isou,jj)

        pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                          +dpvf(2)*diipfv(2)        &
                          +dpvf(3)*diipfv(3))
        pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                          +dpvf(2)*djjpfv(2)        &
                          +dpvf(3)*djjpfv(3))

!         CENTRE
!        --------

        if (ischcp.eq.1) then

          pif = pnd*pip +(1.d0-pnd)*pjp
          pjf = pif


!         SECOND ORDER
!        --------------

        elseif(ischcp.eq.0) then
! dif* is already defined

!     on laisse la reconstruction de PIF et PJF meme si IRCFLP=0
!     sinon cela revient a faire de l'upwind
        pif = pi
        pjf = pj
        do jsou = 1, 3
          pif = pif + gradv(isou,jsou,ii)*difv(jsou)
          pjf = pjf + gradv(isou,jsou,jj)*djfv(jsou)
        enddo

        else
          write(nfecra,9000)ischcp
          iok = 1
        endif


!        BLENDING
!       ----------

        pif = blencp*pif+(1.d0-blencp)*pi
        pjf = blencp*pjf+(1.d0-blencp)*pj


!        FLUX
!       ------

        flux = iconvp*( flui*pif +fluj*pjf )                        &
             + idiffp*viscf(ifac)*( pip -pjp )


!        ASSEMBLAGE
!       ------------

        smbr(isou,ii) = smbr(isou,ii) - thetap * flux
        smbr(isou,jj) = smbr(isou,jj) + thetap * flux
      enddo
      !end isou

    enddo
!        Position "hors boucle" du CALL CSEXIT pour raisons
!        historiques,pour eviter de devectoriser la boucle
!        -> a conserver si on recherche a optimiser la boucle en
!        vectorisation
    if(iok.ne.0) then
      call csexit (1)
    endif

  endif




!  --> FLUX AVEC TEST DE PENTE (separe pour vectorisation eventuelle)
!  =============================

else

!     Stationnaire
  if (idtvar.lt.0) then

    iok = 0
    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      do jsou = 1, 3
        dijpfv(jsou) = dijpf(jsou,ifac)
      enddo

      pnd   = pond(ifac)
      distf = dist(ifac)
      srfan = surfan(ifac)
! ON RECALCULE A CE NIVEAU II' ET JJ'
      do jsou = 1, 3
        diipfv(jsou) = cdgfac(jsou,ifac) - (xyzcen(jsou,ii)+                  &
               (1.d0-pnd) * dijpfv(jsou))
        djjpfv(jsou) = cdgfac(jsou,ifac) -  xyzcen(jsou,jj)+                  &
                   pnd  * dijpfv(jsou)
      enddo


      flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
      fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )

!  Pour le second ordre on definit IF
      if(ischcp.eq.0) then
        do jsou = 1, 3
          difv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,ii)
          djfv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,jj)
        enddo
      endif

!-----------------
! X-Y-Z components, p=u, v, w
      do isou = 1, 3

        do jsou = 1, 3
          dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
        enddo

        pi  = vel (isou,ii)
        pj  = vel (isou,jj)

        pia = vela(isou,ii)
        pja = vela(isou,jj)

        pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                          +dpvf(2)*diipfv(2)        &
                          +dpvf(3)*diipfv(3))
        pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                          +dpvf(2)*djjpfv(2)        &
                          +dpvf(3)*djjpfv(3))

        pipr = pi /relaxp - (1.d0-relaxp)/relaxp * pia   &
             + ircflp*(dpvf(1)*diipfv(1)                 &
                      +dpvf(2)*diipfv(2)                 &
                      +dpvf(3)*diipfv(3) )
        pjpr = pj /relaxp - (1.d0-relaxp)/relaxp * pja   &
             + ircflp*(dpvf(1)*djjpfv(1)                 &
                      +dpvf(2)*djjpfv(2)                 &
                      +dpvf(3)*djjpfv(3))


        pir = pi /relaxp - (1.d0 - relaxp)/relaxp*pia
        pjr = pj /relaxp - (1.d0 - relaxp)/relaxp*pja

!         TEST DE PENTE
!        ---------------

        testi = gradva(isou,1,ii)*surfac(1,ifac)                  &
              + gradva(isou,2,ii)*surfac(2,ifac)                  &
              + gradva(isou,3,ii)*surfac(3,ifac)
        testj = gradva(isou,1,jj)*surfac(1,ifac)                  &
              + gradva(isou,2,jj)*surfac(2,ifac)                  &
              + gradva(isou,3,jj)*surfac(3,ifac)
        testij= gradva(isou,1,ii)*gradva(isou,1,jj)             &
              + gradva(isou,2,ii)*gradva(isou,2,jj)             &
              + gradva(isou,3,ii)*gradva(isou,3,jj)

        if( flumas(ifac).gt.0.d0) then
          dcc = gradv(isou,1,ii)*surfac(1,ifac)    &
              + gradv(isou,2,ii)*surfac(2,ifac)    &
              + gradv(isou,3,ii)*surfac(3,ifac)
          ddi = testi
          ddj = ( pj - pi )/distf *srfan
        else
          dcc = gradv(isou,1,jj)*surfac(1,ifac)    &
              + gradv(isou,2,jj)*surfac(2,ifac)    &
              + gradv(isou,3,jj)*surfac(3,ifac)
          ddi = ( pj - pi )/distf *srfan
          ddj = testj
        endif
        tesqck = dcc**2 -(ddi-ddj)**2


!         UPWIND
!        --------

        if( tesqck.le.0.d0 .or. testij.le.0.d0 ) then
          pifri = pir
          pifrj = pi
          pjfri = pj
          pjfrj = pjr
!     en parallele, la face sera comptee d'un cote OU (exclusif) de l'autre
          if (ii.le.ncel) then
            infac = infac+1
          endif

        else


!         CENTRE
!        --------

          if (ischcp.eq.1) then

            pifri = pnd*pipr +(1.d0-pnd)*pjp
            pjfri = pifri
            pifrj = pnd*pip  +(1.d0-pnd)*pjpr
            pjfrj = pifrj


!         SECOND ORDER
!        --------------

          elseif(ischcp.eq.0) then
! difv already defined

!     on laisse la reconstruction de PIF et PJF meme si IRCFLP=0
!     sinon cela revient a faire de l'upwind
            pifri = pir + difv(1)*gradv(isou,1,ii)      &
                        + difv(2)*gradv(isou,2,ii)      &
                        + difv(3)*gradv(isou,3,ii)
            pifrj = pi  + difv(1)*gradv(isou,1,ii)      &
                        + difv(2)*gradv(isou,2,ii)      &
                        + difv(3)*gradv(isou,3,ii)

            pjfrj = pjr + djfv(1)*gradv(isou,1,jj)      &
                        + djfv(2)*gradv(isou,2,jj)      &
                        + djfv(3)*gradv(isou,3,jj)
            pjfri = pj  + djfv(1)*gradv(isou,1,jj)      &
                        + djfv(2)*gradv(isou,2,jj)      &
                        + djfv(3)*gradv(isou,3,jj)
          else
            write(nfecra,9000)ischcp
            iok = 1
          endif

        endif


!        BLENDING
!       ----------

        pifri = blencp*pifri+(1.d0-blencp)*pir
        pifrj = blencp*pifrj+(1.d0-blencp)*pi
        pjfri = blencp*pjfri+(1.d0-blencp)*pj
        pjfrj = blencp*pjfrj+(1.d0-blencp)*pjr


!        FLUX
!       ------

        fluxi = iconvp*( flui*pifri + fluj*pjfri )                  &
               +idiffp*viscf(ifac)*( pipr -pjp )
        fluxj = iconvp*( flui*pifrj +fluj*pjfrj )                   &
               +idiffp*viscf(ifac)*( pip -pjpr )


!        ASSEMBLAGE
!       ------------

        smbr(isou,ii) = smbr(isou,ii) - fluxi
        smbr(isou,jj) = smbr(isou,jj) + fluxj
      enddo
      !end isou

    enddo
!        Position "hors boucle" du CALL CSEXIT pour raisons
!        historiques,pour eviter de devectoriser la boucle
!        -> a conserver si on recherche a optimiser la boucle en
!        vectorisation
    if(iok.ne.0) then
      call csexit (1)
    endif

!     Instationnaire
  else

    iok = 0
    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      do jsou = 1, 3
        dijpfv(jsou) = dijpf(jsou,ifac)
      enddo

      pnd   = pond(ifac)
      distf = dist(ifac)
      srfan = surfan(ifac)

! ON RECALCULE II' ET JJ'
! ON RECALCULE A CE NIVEAU II' ET JJ'
      do jsou = 1, 3
        diipfv(jsou) = cdgfac(jsou,ifac) - (xyzcen(jsou,ii)+                  &
               (1.d0-pnd) * dijpfv(jsou))
        djjpfv(jsou) = cdgfac(jsou,ifac) -  xyzcen(jsou,jj)+                  &
                   pnd  * dijpfv(jsou)
      enddo

      flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
      fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )

!  Pour le second ordre on definit IF
      if(ischcp.eq.0) then
        do jsou = 1, 3
          difv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,ii)
          djfv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,jj)
        enddo
      endif

!-----------------
! X-Y-Z components, p=u, v, w
      do isou = 1, 3

        do jsou = 1, 3
          dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
        enddo

        pi = vel(isou,ii)
        pj = vel(isou,jj)

        pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                          +dpvf(2)*diipfv(2)        &
                          +dpvf(3)*diipfv(3))
        pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                          +dpvf(2)*djjpfv(2)        &
                          +dpvf(3)*djjpfv(3))

!         TEST DE PENTE
!        ---------------

        testi = gradva(isou,1,ii)*surfac(1,ifac)                  &
              + gradva(isou,2,ii)*surfac(2,ifac)                  &
              + gradva(isou,3,ii)*surfac(3,ifac)
        testj = gradva(isou,1,jj)*surfac(1,ifac)                  &
              + gradva(isou,2,jj)*surfac(2,ifac)                  &
              + gradva(isou,3,jj)*surfac(3,ifac)
        testij= gradva(isou,1,ii)*gradva(isou,1,jj)               &
              + gradva(isou,2,ii)*gradva(isou,2,jj)               &
              + gradva(isou,3,ii)*gradva(isou,3,jj)

        if( flumas(ifac).gt.0.d0) then
          dcc = gradv(isou,1,ii)*surfac(1,ifac)    &
              + gradv(isou,2,ii)*surfac(2,ifac)    &
              + gradv(isou,3,ii)*surfac(3,ifac)
          ddi = testi
          ddj = ( pj - pi )/distf *srfan
        else
          dcc = gradv(isou,1,jj)*surfac(1,ifac)    &
              + gradv(isou,2,jj)*surfac(2,ifac)    &
              + gradv(isou,3,jj)*surfac(3,ifac)
          ddi = ( pj - pi )/distf *srfan
          ddj = testj
        endif
        tesqck = dcc**2 -(ddi-ddj)**2


!         UPWIND
!        --------

      if( tesqck.le.0.d0 .or. testij.le.0.d0 ) then
        pif = pi
        pjf = pj
!     en parallele, la face sera comptee d'un cote OU (exclusif) de l'autre
        if (ii.le.ncel) then
          infac = infac+1
        endif

      else


!         CENTRE
!        --------

        if (ischcp.eq.1) then

          pif = pnd*pip +(1.d0-pnd)*pjp
          pjf = pif


!         SECOND ORDER
!        --------------

        elseif(ischcp.eq.0) then
! dif* is already defined
          pif = pi
          pjf = pj
          do jsou = 1, 3
            pif = pif + gradv(isou,jsou,ii)*difv(jsou)
            pjf = pjf + gradv(isou,jsou,jj)*djfv(jsou)
          enddo

!     on laisse la reconstruction de PIF et PJF meme si IRCFLP=0
!     sinon cela revient a faire de l'upwind

        else
          write(nfecra,9000)ischcp
          iok = 1
        endif

      endif


!        BLENDING
!       ----------

        pif = blencp*pif+(1.d0-blencp)*pi
        pjf = blencp*pjf+(1.d0-blencp)*pj


!        FLUX
!       ------

        flux = iconvp*( flui*pif +fluj*pjf )                        &
             + idiffp*viscf(ifac)*( pip -pjp )


!        ASSEMBLAGE
!       ------------

        smbr(isou,ii) = smbr(isou,ii) - thetap * flux
        smbr(isou,jj) = smbr(isou,jj) + thetap * flux
      enddo
      !end isou

    enddo
!        Position "hors boucle" du CALL CSEXIT pour raisons
!        historiques,pour eviter de devectoriser la boucle
!        -> a conserver si on recherche a optimiser la boucle en
!        vectorisation
    if(iok.ne.0) then
      call csexit (1)
    endif

  endif

endif



if(iwarnp.ge.2) then
  if (irangp.ge.0) call parcpt(infac)
  write(nfecra,1100)cnom,infac,nfacgb
endif


! ======================================================================
! ---> ASSEMBLAGE A PARTIR DES FACETTES DE BORD
! ======================================================================

!     Stationnaire
if (idtvar.lt.0) then

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    do jsou = 1, 3
      diipbv(jsou) = diipb(jsou,ifac)
    enddo

    ! On enleve le decentrement pour les faces couplees
    if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
      flui = 0.0d0
      fluj = flumab(ifac)
    else
      flui = 0.5d0*( flumab(ifac) +abs(flumab(ifac)) )
      fluj = 0.5d0*( flumab(ifac) -abs(flumab(ifac)) )
    endif
!-----------------
! X-Y-Z components, p=u, v, w
    do isou = 1, 3
      pfac  = inc*coefav(isou,ifac)
      pfacd = inc*cofafv(isou,ifac)

      !coefu and cofuf are a matrices
      do jsou = 1, 3
        pir  = vel(jsou,ii)/relaxp - (1.d0-relaxp)/relaxp*vela(jsou,ii)

        pipr = pir +ircflp*( gradv(jsou,1,ii)*diipbv(1)           &
                           + gradv(jsou,2,ii)*diipbv(2)           &
                           + gradv(jsou,3,ii)*diipbv(3)    )
        pfac  = pfac  + coefbv(isou,jsou,ifac)*pipr
        pfacd = pfacd + cofbfv(isou,jsou,ifac)*pipr

      enddo

      pir  = vel(isou,ii)/relaxp - (1.d0-relaxp)/relaxp*vela(isou,ii)
      pipr = pir +ircflp*( gradv(isou,1,ii)*diipbv(1)             &
                         + gradv(isou,2,ii)*diipbv(2)             &
                         + gradv(isou,3,ii)*diipbv(3)    )

      flux = iconvp*( flui*pir +fluj*pfac )                       &
           + idiffp*viscb(ifac)*( pipr -pfacd )
      smbr(isou,ii) = smbr(isou,ii) - flux
    enddo
    !end isou

  enddo

!     Instationnaire
else

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    do jsou = 1, 3
      diipbv(jsou) = diipb(jsou,ifac)
    enddo

    ! On enleve le decentrement pour les faces couplees
    if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
      flui = 0.0d0
      fluj = flumab(ifac)
    else
      flui = 0.5d0*( flumab(ifac) +abs(flumab(ifac)) )
      fluj = 0.5d0*( flumab(ifac) -abs(flumab(ifac)) )
    endif

!-----------------
! X-Y-Z components, p=u, v, w
    do isou = 1, 3
      pfac  = inc*coefav(isou,ifac)
      pfacd = inc*cofafv(isou,ifac)

      !coefu and cofuf are a matrices
      do jsou = 1, 3
        pip = vel(jsou,ii) +ircflp*( gradv(jsou,1,ii)*diipbv(1)           &
                                   + gradv(jsou,2,ii)*diipbv(2)           &
                                   + gradv(jsou,3,ii)*diipbv(3)    )
        pfac  = pfac  + coefbv(isou,jsou,ifac)*pip
        pfacd = pfacd + cofbfv(isou,jsou,ifac)*pip
      enddo

      pip = vel(isou,ii) +ircflp*( gradv(isou,1,ii)*diipbv(1)             &
                                 + gradv(isou,2,ii)*diipbv(2)             &
                                 + gradv(isou,3,ii)*diipbv(3)    )

      flux = iconvp*( flui*vel(isou,ii) +fluj*pfac )                      &
           + idiffp*viscb(ifac)*( pip -pfacd )
      smbr(isou,ii) = smbr(isou,ii) - thetap * flux
    enddo
    !end isou

  enddo

endif

!===============================================================================
! 3.  COMPUTATION OF THE TRANSPOSE GRAD(VEL) TERM AND GRAD(-2/3 DIV(VEL))
!===============================================================================

if (ivisep.eq.1) then

  ! We do not know what condition to put in the inlets and the outlets, so we
  ! assume that there is an equilibrium

  ! Allocate a temporary array
  allocate(bndcel(ncelet))

  do iel = 1, ncelet
    bndcel(iel) = 1.d0
  enddo
  do ifac = 1, nfabor
    ityp = itypfb(ifac)
    if (ityp.eq.isolib .or. ityp.eq.ientre) bndcel(ifabor(ifac)) = 0.d0
  enddo

  ! ---> INTERNAL FACES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pnd = pond(ifac)
    secvis = secvif(ifac)
    visco = viscf(ifac)

    grdtrv =        pnd*(gradv(1,1,ii)+gradv(2,2,ii)+gradv(3,3,ii))   &
           + (1.d0-pnd)*(gradv(1,1,jj)+gradv(2,2,jj)+gradv(3,3,jj))

    ! We need to compute trans_grad(u).IJ which is equal to IJ.grad(u)

    do isou = 1, 3

      tgrdfl = dijpf(1,ifac) * (        pnd*gradv(1,isou,ii)         &
                               + (1.d0-pnd)*gradv(1,isou,jj) )       &
             + dijpf(2,ifac) * (        pnd*gradv(2,isou,ii)         &
                               + (1.d0-pnd)*gradv(2,isou,jj) )       &
             + dijpf(3,ifac) * (        pnd*gradv(3,isou,ii)         &
                               + (1.d0-pnd)*gradv(3,isou,jj) )


      flux = visco*tgrdfl + secvis*grdtrv*surfac(isou,ifac)

      smbr(isou,ii) = smbr(isou,ii) + idiffp*flux*bndcel(ii)
      smbr(isou,jj) = smbr(isou,jj) - idiffp*flux*bndcel(jj)

    enddo

  enddo

  ! ---> BOUNDARY FACES
  !      the whole flux term of the stress tensor is already taken into account
  !TODO add the corresponding term in forbr
  !TODO in theory we should take the normal component into account (the
  !tangential one is modeled by the wall law)

  !Free memory
  deallocate(bndcel)

endif

! Free memory
deallocate(gradva)
deallocate(gradv)

!--------
! FORMATS
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
'@ @@ ATTENTION : ARRET DANS bilsc4                           ',/,&
'@    =========                                               ',/,&
'@     APPEL DE bilsc4 POUR ',A8 ,' AVEC ISCHCP = ',I10        ,/,&
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
'@ @@ WARNING: ABORT IN bilsc4                                ',/,&
'@    ========                                                ',/,&
'@     CALL OF bilsc4 FOR ',A8 ,' WITH ISCHCP = ',I10          ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Contact the support.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! FIN
!----

return

end subroutine
