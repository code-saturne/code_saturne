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

!> \file tridim.f90
!> \brief Resolution of incompressible Navier Stokes and scalar transport
!> equations for a time step.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     itrale        ALE iteration number
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     isostd        standard output indicator
!>                              + reference face number
!> \param[in]     dt            time step (per cell)
!> \param[in]     tpucou        velocity-pressure coupling
!> \param[in]     rtpa          calculated variables at cell centers
!>                              (at current and previous time steps)
!> \param[in]     rtp           calculated variables at cell centers
!>                              (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[in]     tslagr        lagrangian returned coupling term
!> \param[in]     coefa         boundary conditions
!> \param[in]     coefb         boundary conditions
!> \param[in]     frcxt         external force generating hydrostatic pressure
!> \param[in]     prhyd         predicted hydrostatic pressure
!______________________________________________________________________________

subroutine tridim &
 ( itrale ,                                                       &
   nvar   , nscal  ,                                              &
   isostd ,                                                       &
   dt     , tpucou , rtpa   , rtp    , propce , propfb ,          &
   tslagr , coefa  , coefb  , frcxt  , prhyd  )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use entsor
use cstphy
use cstnum
use pointe
use albase
use alstru
use alaste
use parall
use period
use ppppar
use ppthch
use ppincl
use cpincl
use coincl
use atincl
use atsoil
use lagpar
use lagdim
use lagran
use vorinc
use ihmpre
use radiat
use cplsat
use ppcpfu
use elincl
use mesh
use field
use cs_f_interfaces

! les " use pp* " ne servent que pour recuperer le pointeur IIZFPP

!===============================================================================

implicit none

! Arguments

integer          itrale
integer          nvar   , nscal

integer          isostd(nfabor+1)

double precision dt(ncelet), tpucou(ncelet,3), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision tslagr(ncelet,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision frcxt(3,ncelet), prhyd(ncelet)

! Local variables

integer          iel   , ifac  , inod  , ivar  , iscal , iappel
integer          ncp   , ncv   , iok
integer          iiptot
integer          nbccou
integer          ntrela

integer          isvhb , isvtb
integer          ii    , jj    , ippcp , ientha, ippcv
integer          ipcrom, ipcroa
integer          iterns, inslst, icvrge
integer          iflmas, iflmab
integer          italim, itrfin, itrfup, ineefl
integer          nbzfmx, nozfmx
integer          isou, ielpdc, ipbrom

double precision cpcst , tditot, tdist2, tdist1, cvcst
double precision xxp0, xyp0, xzp0
double precision relaxk, relaxe, relaxw, relaxn
double precision cosdto, sindto, omgnrm, rrotgb(3,3)
double precision hdls(6)

integer          ipass
data             ipass /0/
save             ipass

integer          infpar
save             infpar

integer, allocatable, dimension(:,:) :: icodcl
integer, allocatable, dimension(:) :: ilzfbr

double precision, allocatable, dimension(:,:) :: uvwk, ximpa, trava
double precision, allocatable, dimension(:,:,:) :: ximpav
double precision, allocatable, dimension(:) :: flmalf, flmalb, xprale
double precision, allocatable, dimension(:,:) :: cofale
double precision, allocatable, dimension(:) :: qcalc
double precision, allocatable, dimension(:,:,:) :: rcodcl
double precision, allocatable, dimension(:) :: hbord, theipb
double precision, allocatable, dimension(:) :: visvdr
double precision, allocatable, dimension(:) :: prdv2f
double precision, dimension(:), pointer :: imasfl, bmasfl

!===============================================================================


!===============================================================================
! 1.  INITIALISATION
!===============================================================================

if(iwarni(iu).ge.1) then
  write(nfecra,1000)
endif


ipass = ipass + 1

! --- Indicateur de stockage d'un scalaire et de son coef
!     d'echange associe.
!     Pour le moment, on stocke uniquement dans le cas couplage SYRTHES.
!     ISVTB donne le numero du scalaire (on suppose qu'il n'y en a qu'un).
!     Dans le cas ou on a un couplage avec le module thermique 1D en paroi,
!     on utilise le meme scalaire que celui qui sert a Syrthes (s'il y a
!     couplage Syrthes), sinon on stocke le scalaire thermique de la phase 1.

call nbcsyr (nbccou)
!==========
isvhb = 0
isvtb = 0
if (nbccou .ge. 1) then
  do iscal = 1, nscal
    if(icpsyr(iscal).eq.1) then
      isvhb = iscal
    endif
  enddo
endif

if ((nfpt1t.gt.0).and.(nbccou.le.0)) then
  isvhb = iscalt
endif

if (iscalt.gt.0) isvtb = iscalt

!     Si la distance a la paroi doit etre mise a jour, on l'initialise a GRAND
!     des maintenant (pour le premier passage dans phyvar en k-omega)
if(ipass.eq.1.and.ineedy.eq.1.and.abs(icdpar).eq.1.and.           &
                                  imajdy.eq.0) then
  do iel = 1, ncel
    dispar(iel) = grand
  enddo
endif

!    Initialisation de pthera (pther est deja initialise ou lu dans
!    un fichier suite pour l'algorithme a masse volumique variable

if (idilat.eq.3) then
  pthera = pther
endif

!===============================================================================
! 2.  AU DEBUT DU CALCUL ON REINITIALISE LA PRESSION
!===============================================================================

! On le fait sur 2 pas de temps, car souvent, le champ de flux de masse
!   initial n'est pas a divergence nulle (CL incluses) et l'obtention
!   d'un flux a divergence nulle coherent avec la contrainte stationnaire
!   peut prendre quelques pas de temps.
! Remarquer que la pression est rattrapee dans l'etape de stokes.
! On ne le fait pas dans le cas de la prise en compte de la pression
!   hydrostatique, ni dans le cas du compressible

if( ntcabs.le.2 .and. isuite.eq.0 .and. (iphydr.eq.0.or.iphydr.eq.2)    &
                .and. ippmod(icompf).lt.0                               &
                .and. idilat .le.1                 ) then

  if(iwarni(ipr).ge.2) then
    write(nfecra,2000) ntcabs
  endif
  iiptot = ipproc(iprtot)
  xxp0   = xyzp0(1)
  xyp0   = xyzp0(2)
  xzp0   = xyzp0(3)
  do iel = 1, ncel
    rtp(iel,ipr) = pred0
    propce(iel,iiptot) = p0                                &
         + ro0*( gx*(xyzcen(1,iel)-xxp0)                   &
         +       gy*(xyzcen(2,iel)-xyp0)                   &
         +       gz*(xyzcen(3,iel)-xzp0) )
  enddo
endif

 2000 format(                                                           &
  ' REINITIALISATION DE LA PRESSION A L''ITERATION ',I10)

!===============================================================================
! 3.  COMMUNICATIONS
!===============================================================================

! ---> On echange ici les variables RTP en debut de pas de temps.
!         Ce n'est pas fait dans caltri pour faciliter la lecture
!         (manipulation des tableaux)
!      Ceci a l'avantage d'echanger egalement ce qui provient de
!         inivar et lecamo (et par rapport a une solution d'echange
!         en fin de pas de temps, on evite une communication)
!      L'inconvenient est que l'on ne dispose pas des halos sitot les
!         RTP calcules (fin de navstv, turbke, turrij, boucle dans scalai)
!      On pourrait penser a echanger aussi les PROPCE et le DT
!         Pour le moment ca ne s'impose pas : on a besoin d'avoir
!         echange RHO dans navstv pour un affichage et les viscosites
!         dans visecv. On fera les transferts localement quand necessaire.
!      Par principe, on suppose que
!         c'est a la charge de celui qui utilise des valeurs voisines de
!           s'assurer qu'elles existent (ss pgm utilisateurs en
!           particulier)
!         seules les RTPA sont echangees de maniere systematique
!           (eviter donc d'utiliser RTP)
!         le calcul du gradient fournit aussi les valeurs sur le halo
!           (utile pour les reconstructions)
!         les seules boucles sur NCELET sont des boucles d'initialisation
!           (on n'a pas a faire de calcul sur les cellules halo, elles
!            elles sont remplies par les routines de communication et
!            on y preleve seulement de l'information)



! ---> Halo synchronization

if (irangp.ge.0 .or. iperio.eq.1) then

  do ivar = 1, nvar
    call synsce (rtp(1,ivar))
    !==========
  enddo

endif

! ---> Periodicity of rotation

if (iperio.eq.1) then

  !  -- Vitesse

  call perrve (rtp(1,iu), rtp(1,iv), rtp(1,iw))
  !==========

  !  -- Reynolds stress tensor

  if (itytur.eq.3) then
    call perrte &
    !==========
  ( rtp(1,ir11), rtp(1,ir12), rtp(1,ir13),           &
    rtp(1,ir12), rtp(1,ir22), rtp(1,ir23),           &
    rtp(1,ir13), rtp(1,ir23), rtp(1,ir33) )
  endif

  !  -- Note for v2f:
  !     v2 (thus phi) is oriented locally, and is handled as a scalar
  !     regarding periodicity of rotation

endif

!===============================================================================
! 4.  POUR IPHYDR ON DOIT COMMUNIQUER FRCXT AU PREMIER PASSAGE
!     (FRCXT SERT DANS TYPECL)
!     SI ICALHY=1, ON COMMUNIQUE AUSSI RHO POUR REMPLIR
!     PROPCE(1,IPPROC(IROMA))
!===============================================================================

if (ipass.eq.1) then

! --- Communication de FRCXT
  if (iphydr.eq.1) then

    if (irangp.ge.0 .or. iperio.eq.1) then
      call synvin (frcxt)
    endif

  endif

! --- Communication de RHO
  if (icalhy.eq.1 .or. idilat.eq.3) then

    ipcrom = ipproc(irom  )
    if (irangp.ge.0 .or. iperio.eq.1) then
      call synsce (propce(1,ipcrom))
      !==========
    endif

  endif

! --- Communication de prhyd
  if (iphydr.eq.2) then

    if (irangp.ge.0 .or. iperio.eq.1) then
      call synsce (prhyd(1))
      !==========
    endif

  endif

endif

!===============================================================================
! 5.  LES VALEURS COURANTES ECRASENT LES VALEURS ANTERIEURES
!===============================================================================

! --- Noter que exceptionnellement, on fait un calcul avec NCELET,
!       pour eviter une nouvelle communication sur RTPA et les autres
!       tableaux du pas de temps precedent

do ivar = 1, nvar
  do iel = 1, ncelet
    rtpa (iel,ivar) = rtp (iel,ivar)
  enddo
enddo

! If required, the density at time step n-1 is updated
if (icalhy.eq.1.or.idilat.gt.1) then
  ipcrom = ipproc(irom  )
  ipcroa = ipproc(iroma )
  do iel = 1, ncelet
    propce(iel,ipcroa) = propce(iel,ipcrom)
  enddo
endif

!===============================================================================
! 6. DANS LE CAS  "zero pas de temps" EN "SUITE" DE CALCUL
!      ON SORT ICI
!===============================================================================
!  on sort avant SCHTMP car sinon a l'ordre 2 en temps la valeur du
!  flux de masse au pas de temps precedent est ecrasee par la valeur
!  au pas de temps actuel et la valeur au pas de temps actuel est
!  remplacee par une extrapolation qui n'a pas lieu d'etre puisque
!  NTCABS n'est pas incremente. Dans le cas INPDT0=1 sans suite, il
!  n'y a pas de probleme puisque tous les flux de masse sont a 0           !
!  Si ITRALE=0, on est a l'iteration d'initialisation de l'ALE,
!  on ne touche pas au flux de masse non plus


if(inpdt0.eq.1.and.isuite.eq.1) goto 200

if (itrale.gt.0) then
  iappel = 1
  call schtmp(nscal, iappel, propce, propfb)
  !==========
endif


!===============================================================================
! 6.  MISE A JOUR DU MAILLAGE POUR UN COUPLAGE ROTOR/STATOR
!===============================================================================

if (imobil.eq.1) then

  ! --- En turbomachine on connaît la valeur exacte de la vitesse de maillage

  omgnrm = sqrt(omegax**2 + omegay**2 + omegaz**2)

  cosdto = cos(ttcmob*omgnrm)
  sindto = sin(ttcmob*omgnrm)

  do ii = 1, 3
    do jj = 1, 3
      rrotgb(ii,jj) = cosdto*irot(ii,jj) + (1.d0 - cosdto)*prot(ii,jj) &
                                         +         sindto *qrot(ii,jj)
    enddo
  enddo

  ! On modifie la géométrie en fonction de la géométrie initiale

  do inod = 1, nnod
    do ii = 1, 3
      xyznod(ii,inod) = 0.d0
      do jj = 1, 3
        xyznod(ii,inod) = xyznod(ii,inod) + rrotgb(ii,jj)*xyzno0(jj,inod)
      enddo
    enddo
  enddo

  call algrma
  !==========

  ! Abort at the end of the current time-step if there is a negative volume
  if (volmin.le.0.d0) ntmabs = ntcabs

endif


!===============================================================================
! 6.  MISE A JOUR DE LA LOCALISATION DES INTERFACES DE COUPLAGE CS/CS
!===============================================================================

! Localisation des interfaces de couplage via la librairie FVM

! On fait cette mise a jour des interfaces de localisation juste apres
! les changements de geometries dus :
!   - soit a la methode ALE (en fin de pas de temps precedent)
!   - soit a un deplacement impose (cf ci-dessus)

if (nbrcpl.gt.0) call cscloc
                 !==========

!===============================================================================
! 7.  CALCUL DES PROPRIETES PHYSIQUES VARIABLES
!      SOIT VARIABLES AU COURS DU TEMPS
!      SOIT VARIABLES LORS D'UNE REPRISE DE CALCUL
!        (VISCOSITES ET MASSE VOLUMIQUE)
!===============================================================================

if(iwarni(iu).ge.1) then
  write(nfecra,1010)
endif

call phyvar                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  )

if (itrale.gt.0) then
  iappel = 2
  call schtmp(nscal, iappel, propce, propfb)
  !==========
endif


! REMPLISSAGE DES COEFS DE PDC
!    ON Y PASSE MEME S'IL N'Y A PAS DE PDC SUR LE PROC COURANT AU CAS OU
!    UN UTILISATEUR DECIDERAIT D'AVOIR UN COEFF DE PDC DEPENDANT DE
!    LA VITESSE MOYENNE OU MAX.


if (ncpdct.gt.0) then

  iappel = 3

  if (iihmpr.eq.1) then
    call uikpdc &
    !==========
  ( iappel, ncelet, ncepdc,             &
    icepdc, ckupdc, rtpa )
  endif

  call uskpdc &
  !==========
( nvar   , nscal  ,                                              &
  ncepdc , iappel ,                                              &
  icepdc , izcpdc ,                                              &
  dt     , rtpa   , rtp    , propce , propfb ,                   &
  ckupdc )

endif


! REMPLISSAGE DES COEFS DE TERME SOURCE DE MASSE

!    ON Y PASSE MEME S'IL N'Y A PAS DE TSM SUR LE PROC COURANT AU CAS OU
!    UN UTILISATEUR DECIDERAIT D'AVOIR UN TSM DEPENDANT DE
!    VALEURS GLOBALES OU MAX.
if(nctsmt.gt.0) then

  !     Mise a zero du tableau de type de TS masse et source
  do ii = 1, ncetsm
    do ivar = 1, nvar
      itypsm(ii,ivar) = 0
      smacel(ii,ivar) = 0.d0
    enddo
  enddo

  iappel = 3
  call  ustsma &
  !===========
( nvar   , nscal  , ncepdc ,                                     &
  ncetsm , iappel ,                                              &
  icepdc ,                                                       &
  icetsm , itypsm , izctsm ,                                     &
  dt     , rtpa   , propce , propfb ,                            &
  ckupdc , smacel )

endif

!===============================================================================
! 8.  CALCUL DU NOMBRE DE COURANT ET DE FOURIER
!     CALCUL DU PAS DE TEMPS SI VARIABLE
!===============================================================================

if(iwarni(iu).ge.1) then
  write(nfecra,1020)
endif

call dttvar &
!==========
 ( nvar   , nscal  ,                                              &
   ncepdc , ncetsm ,                                              &
   iwarni(iu)   ,                                                 &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel )

if (nbaste.gt.0.and.itrale.gt.nalinf) then
  ntrela = ntcabs - ntpabs
  call astpdt(dt, ncelet, ntrela)
  !==========
endif

! Compute the pseudo tensorial time step if needed for the pressure solving
if (idften(ipr).eq.6) then
  do iel = 1, ncel
    dttens(1, iel) = dt(iel)
    dttens(2, iel) = dt(iel)
    dttens(3, iel) = dt(iel)
    dttens(4, iel) = 0.d0
    dttens(5, iel) = 0.d0
    dttens(6, iel) = 0.d0
  enddo
  do ielpdc = 1, ncepdc
    iel = icepdc(ielpdc)

    ! dttens = (1/dt + Kpdc)^-1
    hdls(1) = ckupdc(ielpdc, 1) + 1.d0/dt(iel)
    hdls(2) = ckupdc(ielpdc, 2) + 1.d0/dt(iel)
    hdls(3) = ckupdc(ielpdc, 3) + 1.d0/dt(iel)
    hdls(4) = ckupdc(ielpdc, 4)
    hdls(5) = ckupdc(ielpdc, 5)
    hdls(6) = ckupdc(ielpdc, 6)

    call symmetric_matrix_inverse(dttens(1, iel), hdls)

  enddo
endif


!===============================================================================
!   RECALAGE DE LA PRESSION Pth ET MASSE VOLUMIQUE rho
!     POUR L'AGORITHME A MASSE VOLUMIQUE VARIABLE.
!===============================================================================

if (idilat.eq.3) then
  call pthrbm(nvar, ncetsm, dt, propce, propfb, smacel)
  !==========
endif

!===============================================================================
! 9.  CHARGEMENT ET TRADUCTION DES CONDITIONS AUX LIMITES
!===============================================================================

if(iwarni(iu).ge.1) then
  write(nfecra,1030)
endif

!  -- Methode des vortex en LES :
!     Definition ou modification eventuelle des parametres
!     Mise a jour des vortex

if (ivrtex.eq.1) then

  iappel = 2
  call usvort &
  !==========
 ( nvar   , nscal  ,                                              &
   iappel ,                                                       &
   dt     , rtpa   ,                                              &
   propce , propfb )

!     Verification des donnees entrees par l'utilisateur
!       (au premier passage seulement)
  if (ipass.eq.1) then
    call vorver ( nfabor , iappel )
    !==========
  endif

  if(irangp.le.0) then
    call vortex
    !==========
  endif

! -- Fin de zone Methode des vortex

endif


! --- Methode ALE : debut de boucle d'implicitation du deplacement des
!       structures. ITRFIN=0 indique qu'on a besoin de refaire une iteration
!       pour Syrthes, T1D ou rayonnement.
italim = 1
itrfin = 1
ineefl = 0
if (iale.eq.1 .and. nalimx.gt.1 .and. itrale.gt.nalinf) then
!     On reserve certains tableaux pour permettre le retour a l'etat
!       initial en fin d'iteration ALE
!       - flux de masse
!       - conditions aux limites de gradient de P et U (car on a un appel
!         a GDRCEL pour les non orthogonalites pour calculer les CL reelles)
!         -> n'est peut-etre pas reellement necessaire
!       - la pression initiale (car RTPA est aussi ecrase dans le cas
!         ou NTERUP>1) -> on pourrait optimiser en ne reservant que si
!         necessaire ...
!       Pas la peine de tester les depassements car on passe dans
!       memcli juste apres.
  allocate(flmalf(nfac))
  allocate(flmalb(nfabor))
  allocate(cofale(nfabor,11))
  allocate(xprale(ncelet))
  ineefl = 1

  if (nbccou.gt.0 .or. nfpt1t.gt.0 .or. iirayo.gt.0) itrfin = 0

endif

300 continue


! --- Boucle sur navstv pour couplage vitesse/pression
!     on s'arrete a NTERUP ou quand TOUTES les phases on converge
!     ITRFUP=0 indique qu'on a besoin de refaire une iteration
!     pour Syrthes, T1D ou rayonnement.
itrfup = 1

if (nterup.gt.1.or.isno2t.gt.0) then

  if (.not.allocated(ximpav)) allocate(ximpav(ndim,ndim,ncelet))
  if (.not.allocated(uvwk)) allocate(uvwk(ndim,ncelet))
  if (.not.allocated(trava)) allocate(trava(ndim,ncelet))

  if (nbccou.gt.0 .or. nfpt1t.gt.0 .or. iirayo.gt.0) itrfup = 0

endif
! Compute the number of variable plus the number of turbulent fluxes
nvarcl = nvar
do iscal = 1, nscal
  if (ityturt(iscal).eq.3) nvarcl = nvarcl + 3
enddo

! Allocate temporary arrays for boundary conditions
if (italim .eq. 1) then
  allocate(icodcl(nfabor,nvarcl))
  allocate(rcodcl(nfabor,nvarcl,3))
endif
if (isvhb.gt.0) then
  allocate(hbord(nfabor))
endif
! Boundary value of the thermal scalar in I'
if (iscalt.gt.0) then
  allocate(theipb(nfabor))
endif
if (itytur.eq.4 .and. idries.eq.1) then
  allocate(visvdr(ncelet))
endif

icvrge = 0
inslst = 0
iterns = 1
do while (iterns.le.nterup)

  call precli(nvar, nscal, icodcl, propfb, rcodcl)
  !==========

  !     - Interface Code_Saturne
  !       ======================

  if (iihmpr.eq.1) then

  ! N.B. Zones de face de bord : on utilise provisoirement les zones des
  !    physiques particulieres, meme sans physique particuliere
  !    -> sera modifie lors de la restructuration des zones de bord

    ipbrom = ipprob(irom  )

    call uiclim &
    !==========
  ( ntcabs, nfabor,                                                &
    nozppm, ncharm, ncharb, nclpch,                                &
    iindef, ientre, iesicf, isspcf, iephcf,                        &
    isopcf, iparoi, iparug, isymet, isolib, isca  ,                &
    ipr   , itempk, ienerg, ipbrom,                                &
    iqimp,  icalke, ientat, ientcp, inmoxy, ientox,                &
    ientfu, ientgb, ientgf, iprofm,                                &
    coejou, dpot,   rtpa,   ielcor,                                &
    ipotr,  ipoti,  ipotva, ncelet,                                &
    itypfb, izfppp, icodcl,                                        &
    dtref,  ttcabs, surfbo, cdgfbo,                                &
    qimp,   qimpat, qimpcp, dh,     xintur,                        &
    timpat, timpcp, tkent ,  fment, distch, rcodcl, propfb)

    if (ippmod(iphpar).eq.0) then

    ! ON NE FAIT PAS DE LA PHYSIQUE PARTICULIERE

      nbzfmx = nbzppm
      nozfmx = nozppm
      allocate(ilzfbr(nbzfmx))
      allocate(qcalc(nozfmx))

      call stdtcl &
      !==========
    ( nbzfmx , nozfmx ,                                              &
      iqimp  , icalke , qimp   , dh , xintur,                        &
      itypfb , izfppp , ilzfbr ,                                     &
      propce , propfb ,                                              &
      rcodcl , qcalc  )

      ! Free memory
      deallocate(ilzfbr)
      deallocate(qcalc)

    endif

  endif

  !     - Sous-programme utilisateur
  !       ==========================

  call cs_user_boundary_conditions &
  !===============================
  ( nvar   , nscal  ,                                              &
    icodcl , itrifb , itypfb , izfppp ,                            &
    dt     , rtp    , rtpa   , propce , propfb ,                   &
    rcodcl )

  !     - Interface Code_Saturne
  !       ======================

  if(iihmpr.eq.1) then

    call uiclve &
    !==========
  ( nfabor, nozppm,                                                &
    iindef, ientre, iesicf, iephcf, isspcf, isopcf,                &
    iparoi, iparug, isymet, isolib,                                &
    itypfb, izfppp )

  endif

  ! -- Methode des vortex en L.E.S. :
  !    (Transfert des vortex dans les tableaux RCODCL)

  if (ivrtex.eq.1) then
    call vor2cl(icodcl, itypfb, rcodcl)
    !==========
  endif

  ! --- Couplage code/code entre deux instances (ou plus) de Code_Saturne
  !       On s'occupe ici du couplage via les faces de bord, et de la
  !       transformation de l'information reçue en condition limite.

  if (nbrcpl.gt.0) then

    call cscfbr &
    !==========
  ( nscal  ,                                                       &
    icodcl , itypfb ,                                              &
    dt     , rtp    , propce ,                                     &
    coefa  , coefb  , rcodcl )

  endif

! -- Synthetic Eddy Method en L.E.S. :
!    (Transfert des structures dans les tableaux rcodcl)

    call synthe &
    !==========
  ( nvar   , nscal  ,                                              &
    iu     , iv     , iw     , ipproc(irom)    ,                   &
    ttcabs ,                                                       &
    dt     , rtpa   , rtp    , propce , propfb ,                   &
    coefa  , coefb  , rcodcl )

  ! -- Methode ALE (CL de vitesse de maillage et deplacement aux noeuds)

  if (iale.eq.1) then

    do ii = 1, nnod
      impale(ii) = 0
      disala(1,ii) = depale(1,ii)
      disala(2,ii) = depale(2,ii)
      disala(3,ii) = depale(3,ii)
    enddo

    ! - Interface Code_Saturne
    !   ======================

    if (iihmpr.eq.1) then

      !TODO add ifresf for the free surface

      call uialcl &
      !==========
    ( nfabor, nozppm,                    &
      ibfixe, igliss, ivimpo,            &
      ialtyb, ipnfbr, nodfbr,            &
      impale,                            &
      depale,                            &
      dtref, ttcabs, ntcabs,             &
      iuma, ivma, iwma,                  &
      rcodcl)

    endif

    call usalcl &
    !==========
  ( itrale ,                                                       &
    nvar   , nscal  ,                                              &
    icodcl , itypfb , ialtyb ,                                     &
    impale ,                                                       &
    dt     , rtp    , rtpa   , propce , propfb ,                   &
    rcodcl , xyzno0 , depale )

    !     Au cas ou l'utilisateur aurait touche DEPALE sans mettre IMPALE=1, on
    !       remet le deplacement initial
    do ii  = 1, nnod
      if (impale(ii).eq.0) then
        depale(1,ii) = xyznod(1,ii)-xyzno0(1,ii)
        depale(2,ii) = xyznod(2,ii)-xyzno0(2,ii)
        depale(3,ii) = xyznod(3,ii)-xyzno0(3,ii)
      endif
    enddo

    !     En cas de couplage de structures, on calcule un deplacement predit
    if (nbstru.gt.0.or.nbaste.gt.0) then

      call strpre &
      !==========
    ( itrale , italim , ineefl ,                                     &
      impale ,                                                       &
      rtpa   ,                                                       &
      coefa  , coefb  ,                                              &
      flmalf , flmalb , xprale , cofale , depale )

    endif

  endif

  !     UNE FOIS CERTAINS CODES DE CONDITIONS LIMITES INITIALISES PAR
  !     L'UTILISATEUR, ON PEUT COMPLETER CES CODES PAR LES COUPLAGES
  !     AUX BORDS (TYPE SYRTHES), SAUF SI ON DOIT Y REPASSER ENSUITE
  !     POUR CENTRALISER CE QUI EST RELATIF AU COUPLAGE AVEC SYRTHES
  !     ON POSITIONNE ICI L'APPEL AU COUPLAGE VOLUMIQUE SYRTHES
  !     UTILE POUR BENIFICER DE LA DERNIERE VITESSE CALCULEE SI ON
  !     BOUCLE SUR U/P.
  !     LE COUPLAGE VOLUMIQUE DOIT ETRE APPELE AVANT LE SURFACIQUE
  !     POUR RESPECTER LE SCHEMA DE COMMUNICATION

  if (itrfin.eq.1 .and. itrfup.eq.1) then

    call cpvosy(isvtb, dt, rtp, rtpa, propce, propfb)
    !==========

    call coupbi(nfabor, nvar, nscal, icodcl, rcodcl)
    !==========

    if (nfpt1t.gt.0) then
      call cou1di(nfabor, isvtb, icodcl, rcodcl)
      !==========
    endif

  endif


  if(iirayo.gt.0 .and. itrfin.eq.1 .and. itrfup.eq.1) then

    call raycli &
    !==========
  ( nvar   , nscal  ,                                              &
    icodcl , itypfb ,                                              &
    izfrad ,                                                       &
    dt     , rtp    , rtpa   , propce , propfb , rcodcl ,          &
    coefa  , coefb  )

  endif

  !     ON CALCULE LES COEFFICIENTS ASSOCIES AUX CONDITIONS LIMITES

  call condli &
  !==========
( nvar   , nscal  , iterns ,                                     &
  isvhb  , isvtb  ,                                              &
  icodcl , isostd ,                                              &
  dt     , rtp    , rtpa   , propce , propfb ,                   &
  rcodcl ,                                                       &
  coefa  , coefb  , visvdr ,                                     &
  hbord  , theipb ,                                              &
  frcxt  )

!     ==============================================
!     Appel de l'interface sol-atmosphere
!     ==============================================

  if (ippmod(iatmos).eq.2.and.iatsoil.eq.1.and.nfmodsol.gt.0) then
    ipcrom = ipproc(irom)
    call solvar(rtp(1,isca(iscalt)),rtp(1,isca(itotwt)),rtp(1,ipr),  &
                propce(1,ipcrom)   , dt ,                            &
                rcodcl , rtp)
  endif

  !     UNE FOIS LES COEFFICIENTS CALCULES, ON PEUT EN DEDUIRE PLUS
  !     FACILEMENT (I.E. SANS RECALCULS INUTILES) LES TERMES A
  !     ENVOYER POUR LES COUPLAGES AUX BORDS (TYPE SYRTHES)


  ! On indique si la variable couplee est l'enthalpie
  ientha = 0
  if(isvtb.gt.0) then
    if(iscsth(isvtb).eq.2) then
      ientha = 1
    endif
  endif

  ! Compressible : on indique si la variable couple est l'energie

  if ( ippmod(icompf).ge.0 ) then
    if(isvtb.gt.0) then
      if(iscsth(isvtb).eq.3) then
        ientha = 2
      endif
    endif
  endif

  ! On recupere le Cp de la phase couplee
  if(icp.gt.0) then
    ippcp = ipproc(icp)
    ncp   = ncelet
    cpcst = 0.d0
  else
    ippcp = 1
    ncp   = 1
    cpcst = cp0
  endif

  ! En compressible et si on couple ave l'energie
  ! on recupere de Cv de la phase couplee

  if ( ippmod(icompf).ge.0 .and. ientha .eq. 2 ) then

    if(icv.gt.0) then
      ippcv = ipproc(icv)
      ncv   = ncelet
      cvcst = 0.d0
    else
      ippcv = 1
      ncv   = 1
      cvcst = cv0
    endif
  else
    ippcv = 1
    ncv   = 1
    cvcst = 0.d0
  endif

  ! On envoie le tout vers SYRTHES, en distinguant CP
  !  constant ou variable
  if (itrfin.eq.1 .and. itrfup.eq.1) then

    call coupbo &
    !==========
  ( nvar   , ncp    , ncv    , ientha ,                            &
    dt     , rtp    , rtpa   , propce ,                            &
    cpcst  , propce(1,ippcp) , cvcst  , propce(1,ippcv),           &
    hbord  , theipb )

    if (nfpt1t.gt.0) then
      call cou1do &
      !==========
    ( nvar   , nscal  , ncp    , nfpt1d ,                            &
      ientha , ifpt1d , iclt1d ,                                     &
      tppt1d , tept1d , hept1d ,                                     &
      fept1d , xlmbt1 , rcpt1d , dtpt1d ,                            &
      dt     , rtpa   , propce , propfb ,                            &
      cpcst  , propce(1,ippcp) , hbord  , theipb )
    endif
  endif

  !     ON N'A PLUS BESOIN DE ISVHB OU ISVHT (POUR HBORD ET TBORD)
  !     A PARTIR D'ICI



  !     CALCUL DE LA DISTANCE A LA PAROI
  !       (Nouvel algorithme. L'ancien est dans condli)
  !     En ALE on ne fait ce calcul qu'a la premiere des
  !       sous-iterations d'implicitation ITALIM, car le maillage
  !       n'est pas modifie dans les sous-iterations suivantes

  if (italim.eq.1) then

    if(ineedy.eq.1.and.iwarny.ge.1) then
      call dmtmps(tdist1)
    endif


    ! On ne fait le calcul que s'il y a des parois, 'dispar'  est reserve
    ! et initialise a GRAND avant. S'il n'y a pas de paroi, il restera = GRAND)

    ! Pour le moment, on suppose que l'on peut se contenter de faire
    !  cela au premier passage, sauf avec les maillages mobiles. Attention donc
    !  aux conditions aux limites variables (faces qui deviennent des parois ou
    !  parois qui deviennent autre chose)

    ! Nombre de faces de paroi
    if(ipass.eq.1) then
      if(ineedy.eq.1) then
        infpar = 0
        do ifac = 1, nfabor
          if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
            infpar = infpar+1
          endif
        enddo
        if(irangp.ge.0) then
          call parcpt(infpar)
        endif
      endif
    endif


    !     On calcule la distance a la paroi
    !          si elle doit etre mise a jour
    !       et si on en a besoin,
    !       et si on a choisi ce mode de calcul,
    if( imajdy.eq.0.and.ineedy.eq.1.and.abs(icdpar).eq.1) then

      !     S'il n'y a pas de paroi, on garde l'initialisation a GRAND
      if(infpar.eq.0) then
        imajdy = 1

        !     S'il y a des parois, il faut calculer
      else


        !     On doit conserver la memoire de memcli a cause de 'uetbor'
        !       dans DISTYP (uniquement en LES avec van Driest mais tant pis)

        call distpr(itypfb, dispar)
        !==========

        !     La distance n'a plus a etre mise a jour sauf en ALE
        if (iale.eq.0) imajdy = 1

      endif
    endif

  endif


  !     CALCUL DE L'AMORTISSEMENT DE VAN DRIEST
  !     OU CALCUL DE Y+ POUR LE LAGRANGIEN


  !     On calcule y+ si on en a besoin

  if( (itytur.eq.4.and.idries.eq.1)                 &
       .or. (iilagr.ge.1 .and. iroule.eq.2) ) then

    !       On calcule si on a demande ce mode de calcul
    !               et s'il y a des parois (si pas de paroi, pas de y+)
    if(abs(icdpar).eq.1.and.infpar.gt.0) then

      !     On doit conserver la memoire de memcli a cause de 'uetbor'
      !       dans DISTYP

      call distyp(itypfb, dispar, propce, yplpar)
      !==========

    endif

  endif

  if (itytur.eq.4 .and. idries.eq.1) then

    !     Pas d'amortissement si pas de paroi
    if(infpar.gt.0) then
      call vandri &
      !==========
    ( ndim   , ncelet , ncel   , nfac   , nfabor ,                 &
      itypfb , ifabor , ifapat,                                    &
      xyzcen , cdgfbo , visvdr , yplpar ,                          &
      propce )
    endif

  endif

  if(ineedy.eq.1.and.iwarny.ge.1) then
    call dmtmps(tdist2)
    tditot = tdist2-tdist1
    write(nfecra,4010)tditot
  endif

!===============================================================================
! 10. DANS LE CAS  "zero pas de temps" EN "NON SUITE" DE CALCUL
!      ON SORT ICI
!===============================================================================

  if (inpdt0.eq.1.and.isuite.eq.0) goto 200

  if (iilagr.eq.3) goto 200

!===============================================================================
! 11. RESOLUTION DE LA VITESSE DE MAILLAGE EN ALE
!===============================================================================

  if (iale.eq.1) then

    if (itrale.eq.0 .or. itrale.gt.nalinf) then

      ! otherwise it is done in navstv.f90
      if (itrale.eq.0) then

        call alelav &
        !==========
      ( rtp    , rtpa   , propce , propfb )

      endif


    endif

    if (itrale.eq.0) goto 200

  endif

!===============================================================================
! 11. CALCUL A CHAMP DE VITESSE NON FIGE :
!      ON RESOUT VITESSE ET TURBULENCE
!    ON SUPPOSE QUE TOUTES LES PHASES SONT FIGEES OU AUCUNE
!===============================================================================

  ! En cas de champ de vitesse fige, on ne boucle pas sur U/P
  if (iccvfg.eq.0) then
  !===============

!===============================================================================
! 12. RESOLUTION QUANTITE DE MOUVEMENT ET MASSE
!===============================================================================

    if(iwarni(iu).ge.1) then
      write(nfecra,1040)
    endif


    call field_get_key_int(ivarfl(iu), kimasf, iflmas)
    call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
    call field_get_val_s(iflmas, imasfl)
    call field_get_val_s(iflmab, bmasfl)

    ! Coupled solving of the velocity components

    call navstv &
    !==========
  ( nvar   , nscal  , iterns , icvrge , itrale ,                   &
    isostd ,                                                       &
    dt     , tpucou , rtp    , rtpa   , propce , propfb ,          &
    tslagr , coefa  , coefb  , frcxt  , prhyd  ,                   &
    trava  , ximpav , uvwk   )


    !     Mise a jour de la pression si on utilise un couplage vitesse/pression
    !       par point fixe
    !     En parallele, l'echange est fait au debut de navstv.
    if(nterup.gt.1) then
      do iel = 1, ncel
        rtpa(iel,ipr) = rtp(iel,ipr)
      enddo
    endif

    !     Si c'est la derniere iteration : INSLST = 1
    if((icvrge.eq.1).or.(iterns.eq.nterup)) then

      ! Si on a besoin de refaire une nouvelle iteration pour SYRTHES,
      ! rayonnement, paroi thermique 1D...
      ! ET que l'on est a la derniere iteration en ALE !

      ! ...alors, on remet a zero les indicateurs de convergence
      if (itrfup.eq.0.and.itrfin.eq.1) then
        itrfup = 1
        icvrge = 0
        iterns = iterns - 1

        ! ...sinon, on termine
      else
        inslst = 1
      endif



      !     On teste le flux de masse
      if ((istmpf.eq.0.and.inslst.eq.0) .or. istmpf.ne.0) then
        iappel = 3
        call schtmp(nscal, iappel, propce, propfb)
        !==========
      endif

      if (inslst.eq.1) goto 100

    endif

  endif ! Fin si calcul sur champ de vitesse figee

  iterns = iterns + 1

enddo

100 continue

! Free memory
if (allocated(hbord)) deallocate(hbord)
if (allocated(theipb)) deallocate(theipb)
if (allocated(visvdr)) deallocate(visvdr)

if (nterup.gt.1) then
  deallocate(uvwk, trava)
  deallocate(ximpav)
endif

! Calcul sur champ de vitesse fige SUITE (a cause de la boucle U/P)
if (iccvfg.eq.0) then
!===============

!===============================================================================
! 13.  DEPLACEMENT DES STRUCTURES EN ALE ET TEST DE BOUCLAGE IMPLICITE
!===============================================================================

  if (nbstru.gt.0.or.nbaste.gt.0) then

    call strdep &
    !==========
  ( itrale , italim , itrfin ,                                     &
    nvar   ,                                                       &
    dt     , rtp    , rtpa   ,                                     &
    coefa  , coefb  ,                                              &
    flmalf , flmalb , cofale , xprale )

    !     On boucle eventuellement sur de deplacement des structures
    if (itrfin.ne.-1) then
      italim = italim + 1
      goto 300
    endif

    ! Free memory
    if (allocated(flmalf)) then
      deallocate(flmalf, flmalb)
      deallocate(cofale)
      deallocate(xprale)
    endif

  endif

  !     On ne passe dans SCHTMP que si ISTMPF.EQ.0 (explicite)
  !     On teste le flux de masse de la phase 1 (toutes les phases sont
  !     necessairement traitees de la meme facon, cf. VERINI)
  !     pour conserver
  if (istmpf.eq.0) then
    iappel = 4
    call schtmp(nscal, iappel, propce, propfb)
    !==========
  endif

!===============================================================================
! 14. RESOLUTION TURBULENCE
!===============================================================================

  iok = 0
  if(iwarni(iu).ge.1) then
    if( itytur.eq.2 .or. itytur.eq.3              &
         .or. itytur.eq.5 .or. iturb.eq.60 ) then
      iok = 1
    endif
    if(iok.eq.1) then
      write(nfecra,1050)
    endif
  endif

  ! Si on est en v2f (phi-fbar ou BL-v2/k), on reserve un tableau
  ! de taille NCELET pour eviter de recalculer la production dans RESV2F
  ! (trois appels a GRDCEL)
  if (itytur.eq.5) then
    allocate(prdv2f(ncelet))
  endif

  if( (itytur.eq.2) .or. (itytur.eq.5) ) then

    call turbke &
    !==========
  ( nvar   , nscal  ,                                              &
    ncepdc , ncetsm ,                                              &
    icepdc , icetsm , itypsm ,                                     &
    dt     , rtp    , rtpa   , propce , propfb ,                   &
    tslagr ,                                                       &
    coefa  , coefb  , ckupdc , smacel ,                            &
    prdv2f )

    if( itytur.eq.5 )  then

      call resv2f &
      !==========
    ( nvar   , nscal  ,                                              &
      ncepdc , ncetsm ,                                              &
      icepdc , icetsm , itypsm ,                                     &
      dt     , rtp    , rtpa   , propce , propfb ,                   &
      coefa  , coefb  , ckupdc , smacel ,                            &
      prdv2f )

      ! Free memory
      deallocate(prdv2f)

    endif

    !  RELAXATION DE K ET EPSILON SI IKECOU=0 EN INSTATIONNAIRE
    if (ikecou.eq.0 .and. idtvar.ge.0) then
      relaxk = relaxv(ik)
      relaxe = relaxv(iep)
      do iel = 1,ncel
        rtp(iel,ik) = relaxk*rtp(iel,ik) + (1.d0-relaxk)*rtpa(iel,ik)
        rtp(iel,iep) = relaxe*rtp(iel,iep) + (1.d0-relaxe)*rtpa(iel,iep)
      enddo
    endif

  else if(itytur.eq.3) then

    ! Calcul de Alpha pour l'EBRSM
    if (iturb.eq.32) then

      call resalp &
      !==========
    ( nvar   , nscal  ,                                              &
      dt     , rtp    , rtpa   , propce , propfb ,                   &
      coefa  , coefb  )

    endif

    call turrij &
    !==========
  ( nvar   , nscal  ,                                              &
    ncepdc , ncetsm ,                                              &
    icepdc , icetsm , itypsm ,                                     &
    dt     , rtp    , rtpa   , propce , propfb ,                   &
    tslagr ,                                                       &
    coefa  , coefb  , ckupdc , smacel )

  else if( iturb.eq.60 ) then

    call turbkw &
    !==========
  ( nvar   , nscal  ,                                              &
    ncepdc , ncetsm ,                                              &
    icepdc , icetsm , itypsm ,                                     &
    dt     , rtp    , rtpa   , propce , propfb ,                   &
    tslagr ,                                                       &
    coefa  , coefb  , ckupdc , smacel )

    !  RELAXATION DE K ET OMEGA SI IKECOU=0
    if (ikecou.eq.0 .and. idtvar.ge.0) then
      relaxk = relaxv(ik )
      relaxw = relaxv(iomg)
      do iel = 1,ncel
        rtp(iel,ik)  = relaxk*rtp(iel,ik) +(1.d0-relaxk)*rtpa(iel,ik)
        rtp(iel,iomg) = relaxw*rtp(iel,iomg)+(1.d0-relaxw)*rtpa(iel,iomg)
      enddo
    endif

  else if( iturb.eq.70 ) then

    call turbsa &
    !==========
  ( nvar   , nscal  ,                                              &
    ncepdc , ncetsm ,                                              &
    icepdc , icetsm , itypsm ,                                     &
    dt     , rtp    , rtpa   , propce , propfb ,                   &
    tslagr   ,                                                     &
    coefa  , coefb  , ckupdc , smacel ,                            &
    itypfb )

    !  RELAXATION DE NUSA
    if (idtvar.ge.0) then
      relaxn = relaxv(inusa)
      do iel = 1,ncel
        rtp(iel,inusa) = relaxn*rtp(iel,inusa)+(1.d0-relaxn)*rtpa(iel,inusa)
      enddo
    endif

  endif

endif  ! Fin si calcul sur champ de vitesse fige SUITE


!     Ici on peut liberer les eventuels tableaux SKW et DIVUKW

!===============================================================================
! 15.  RESOLUTION DES SCALAIRES
!===============================================================================


if (nscal.ge.1 .and. iirayo.gt.0) then

  if (iwarni(iu).ge.1 .and. mod(ntcabs,nfreqr).eq.0) then
    write(nfecra,1070)
  endif

  call raydom &
  !==========
 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   izfrad ,                                                       &
   dt     , rtp    , rtpa   , propce , propfb )

endif


if (nscal.ge.1) then

  if(iwarni(iu).ge.1) then
    write(nfecra,1060)
  endif

  call scalai                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   tslagr , coefa  , coefb  )

endif

! Free memory
deallocate(icodcl, rcodcl)

!===============================================================================
! 16.  TRAITEMENT DU FLUX DE MASSE, DE LA VISCOSITE,
!      DE LA MASSE VOLUMIQUE ET DE LA CHALEUR SPECIFIQUE POUR
!      UN THETA SCHEMA
!===============================================================================


iappel = 5
call schtmp(nscal, iappel, propce, propfb)
!==========

!===============================================================================
! 17.  SORTIE DANS LE CAS DE "zero pas de temps" ET INIT ALE
!===============================================================================

 200  continue


!===============================================================================

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  INITIALISATIONS                                            ',/,&
'  ===============                                            ',/)
 1010 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  CALCUL DES GRANDEURS PHYSIQUES                             ',/,&
'  ==============================                             ',/)
 1020 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  CALCUL DU CFL, DU FOURIER ET DU DT VARIABLE                ',/,&
'  ===========================================                ',/)
 1030 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  MISE EN PLACE DES CONDITIONS AUX LIMITES                   ',/,&
'  ========================================                   ',/)
 1040 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  RESOLUTION DES EQUATIONS DE NAVIER-STOKES                  ',/,&
'  =========================================                  ',/)
 1050 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  RESOLUTION DES EQUATIONS DES VARIABLES TURBULENTES         ',/,&
'  ==================================================         ',/)
 1060 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  RESOLUTION DES EQUATIONS SUR L''ENERGIE ET LES SCALAIRES   ',/,&
'  ========================================================   ',/)
 1070 format(/,                                                   &
 '------------------------------------------------------------',/,&
                                                              /,/,&
 ' RESOLUTION DES TRANSFERTS THERMIQUES RADIATIFS             ',/,&
'  ==============================================             ',/)

 4010 format(/,                                                   &
' ** TEMPS POUR LA DISTANCE A LA PAROI : ',E14.5               ,/,&
'    ---------------------------------                        ',/)

#else

 1000 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  INITIALISATIONS                                            ',/,&
'  ===============                                            ',/)
 1010 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  COMPUTATION OF PHYSICAL QUANTITIES                         ',/,&
'  ==================================                         ',/)
 1020 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  COMPUTATION OF CFL, FOURIER AND VARIABLE DT                ',/,&
'  ===========================================                ',/)
 1030 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  SETTING UP THE BOUNDARY CONDITIONS                         ',/,&
'  ==================================                         ',/)
 1040 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  SOLVING NAVIER-STOKES EQUATIONS                            ',/,&
'  ===============================                            ',/)
 1050 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  SOLVING TURBULENT VARIABLES EQUATIONS                      ',/,&
'  =====================================                      ',/)
 1060 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  SOLVING ENERGY AND SCALARS EQUATIONS                       ',/,&
'  ====================================                       ',/)
 1070 format(/,                                                   &
 '------------------------------------------------------------',/,&
                                                              /,/,&
 ' SOLVING THERMAL RADIATIVE TRANSFER                         ',/,&
'  ==================================                         ',/)

 4010 format(/,                                                   &
' ** TIME FOR THE WALL DISTANCE: ',E14.5                       ,/,&
'    ---------------------------                              ',/)

#endif

!----
! FIN
!----

end subroutine
