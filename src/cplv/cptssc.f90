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

subroutine cptssc &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp , ia     ,                            &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt ,                                              &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME CHARBON PULVERISE
!   ON PRECISE LES TERMES SOURCES POUR UN SCALAIRE PP
!   SUR UN PAS DE TEMPS

! ATTENTION : LE TRAITEMENT DES TERMES SOURCES EST DIFFERENT
! ---------   DE CELUI DE USTSSC.F

! ON RESOUT ROVSDT*D(VAR) = SMBRS

! ROVSDT ET SMBRS CONTIENNENT DEJA D'EVENTUELS TERMES SOURCES
!  UTILISATEUR. IL FAUT DONC LES INCREMENTER ET PAS LES
!  ECRASER

! POUR DES QUESTIONS DE STABILITE, ON NE RAJOUTE DANS ROVSDT
!  QUE DES TERMES POSITIFS. IL N'Y A PAS DE CONTRAINTE POUR
!  SMBRS

! DANS LE CAS D'UN TERME SOURCE EN CEXP + CIMP*VAR ON DOIT
! ECRIRE :
!          SMBRS  = SMBRS  + CEXP + CIMP*VAR
!          ROVSDT = ROVSDT + MAX(-CIMP,ZERO)

! ON FOURNIT ICI ROVSDT ET SMBRS (ILS CONTIENNENT RHO*VOLUME)
!    SMBRS en kg variable/s :
!     ex : pour la vitesse            kg m/s2
!          pour les temperatures      kg degres/s
!          pour les enthalpies        Joules/s
!    ROVSDT en kg /s

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! itypfb(nfabor    ! te ! --> ! type des faces de bord                         !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
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
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          itypfb(nfabor)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          izfppp(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smbrs(ncelet), rovsdt(ncelet)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          ivar , ipcrom, iel
integer          numcla , numcha , icla , numtra
integer          ipcgch , ipcgd1 , ipcgd2 , ipcght , ipcsec
integer          ixchcl , ixckcl , iscala
integer          ipcro2 , ipcte1 , ipcte2 , ipcvsl , ipccp
integer          ipcdia , ipcvst
integer          mode, ige
integer          ifac   , ifinra , icoefa , icoefb
integer          ipcx2c , icha , imode , ii
integer          iterch
integer          itermx,nbpauv,nbrich,nbepau,nberic,ipghc2

double precision epsrgp , climgp , extrap , xnuss
double precision aux, rhovst
double precision coefe(ngazem)
double precision t1, t2, hlv , hh2ov
double precision f1mc(ncharm), f2mc(ncharm)
double precision xhdev1 , xhdev2 , xhco , xho2 , gamdv1 , gamdv2

double precision gamhet

double precision xxco,xxo2,xxco2,xxh2o,xco2mx
double precision xkp,xk0p,xkm,xk0m,wplus,wmoins,t0p,t0m
double precision auxp,auxm

double precision anmr,xcot,xo2t,xco2e,xo2e,xcoe,tauchi,tautur
double precision sqh2o , x2

double precision err1mx,err2mx,anmr0

double precision errch,ter1,ddelta,xden

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w11

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet), w6(ncelet))
allocate(w11(ncelet))

! Initialize variables to avoid compiler warnings

ipghc2 = 0

! Memoire

idebia = idbia0
idebra = idbra0

! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
chaine = nomvar(ipprtp(ivar))

! --- Numero des grandeurs physiques (voir usclim)
ipcrom = ipproc(irom)
ipcvst = ipproc(ivisct)

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES POUR LES VARIABLES RELATIVES
!    AUX CLASSES DE PARTICULES
!===============================================================================

! --> Terme source pour la fraction massique de charbon reactif

if ( ivar.ge.isca(ixch(1)) .and. ivar.le.isca(ixch(nclacp)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif
  numcla = ivar-isca(ixch(1))+1
  ipcgch = ipproc(igmdch(numcla))

  do iel = 1, ncel

! ---- Calcul de  W1 = - rho.GMDCH > 0

    w1(iel) = - propce(iel,ipcrom)*propce(iel,ipcgch)             &
               *volume(iel)

! ---- Calcul des parties explicite et implicite du TS

    rovsdt(iel) = rovsdt(iel) + max(w1(iel),zero)
    smbrs(iel) = smbrs(iel) - w1(iel)*rtpa(iel,ivar)

  enddo

endif


! --> Terme source pour la fraction massique de coke

if ( ivar.ge.isca(ixck(1)) .and. ivar.le.isca(ixck(nclacp)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  numcla = ivar-isca(ixck(1))+1
  ipcgch = ipproc(igmdch(numcla))
  ipcgd1 = ipproc(igmdv1(numcla))
  ipcgd2 = ipproc(igmdv2(numcla))
  ixchcl = isca(ixch(numcla))
  ixckcl = isca(ixck(numcla))
  ipcght = ipproc(igmhet(numcla))
  if ( ihtco2 .eq. 1 ) then
    ipghc2 = ipproc(ighco2(numcla))
  endif

  do iel = 1, ncel

! ---- Calcul de W1 = - rho.Xch.GMDCH.Volume > 0

    w1(iel) = - propce(iel,ipcrom)*rtp(iel,ixchcl)                &
                *propce(iel,ipcgch)*volume(iel)

! AE : On prend RTP(IEL,IXCHCL) et pas RTPA(IEL,IXCHCL) afin
!      d'etre conservatif sur la masse

! ---- Calcul de W2 = rho.Xch.(GMDV1+GMDV2)Volume < 0

    w2(iel) = propce(iel,ipcrom)*rtp(iel,ixchcl)                  &
              *(propce(iel,ipcgd1)+propce(iel,ipcgd2))            &
              *volume(iel)

    if ( rtpa(iel,ixckcl) .gt. epsicp ) then

! Reaction C(s) + O2 ---> 0.5CO
! =============================

! ---- Calcul de la partie implicite  > 0 du TS relatif a GMHET

      w3(iel) = - 2.d0/3.d0*propce(iel,ipcrom)*propce(iel,ipcght) &
               /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

! ---- Calcul de la partie explicite < 0 du TS relatif a GMHET

      w4(iel) = propce(iel,ipcrom)*propce(iel,ipcght)             &
                * (rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

    else
      w3(iel) = 0.d0
      w4(iel) = 0.d0
    endif

! ---- Calcul des parties explicite et implicite du TS

    rovsdt(iel) = rovsdt(iel) + max(w3(iel),zero)
    smbrs(iel)  = smbrs(iel)  + w1(iel) + w2(iel) + w4(iel)

  enddo

  if ( ihtco2 .eq. 1 ) then

    do iel = 1, ncel

      if ( rtpa(iel,ixckcl) .gt. epsicp ) then

! Reaction C(s) + CO2 ---> 2CO
! =============================

! ---- Calcul de la partie implicite  > 0 du TS relatif a GMHET

        w5(iel) = - 2.d0/3.d0*propce(iel,ipcrom)                  &
                             *propce(iel,ipghc2)                  &
                 /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

! ---- Calcul de la partie explicite < 0 du TS relatif a GMHET

        w6(iel) = propce(iel,ipcrom)*propce(iel,ipghc2)           &
                 *(rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

      else
        w5(iel) = 0.d0
        w6(iel) = 0.d0
      endif

! ---- Calcul des parties explicite et implicite du TS

      rovsdt(iel) = rovsdt(iel) + max(w5(iel),zero)
      smbrs(iel)  = smbrs(iel)  + w6(iel)

    enddo

  endif

endif


! --> Terme source pour la fraction massique de d'eau

if ( ippmod(icp3pl) .eq. 1 ) then

  if ( ivar.ge.isca(ixwt(1)) .and.                                &
       ivar.le.isca(ixwt(nclacp)) ) then

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif
    numcla = ivar-isca(ixwt(1))+1
    numcha = ichcor(numcla)

    ipcsec = ipproc(igmsec(numcla))
    ipcx2c = ipproc(ix2(numcla))

    do iel = 1, ncel

! ---- Calcul des parties explicite et implicite du TS

     if ( rtpa(iel,ivar).gt. epsicp .and.                         &
          xwatch(numcha).gt. epsicp       ) then
       w1(iel) = propce(iel,ipcrom)*propce(iel,ipcsec)            &
                     *volume(iel)                                 &
                     *(1.d0/propce(iel,ipcx2c))                   &
                     *(1.d0/xwatch(numcha))

        rovsdt(iel) = rovsdt(iel) + max(w1(iel),zero)
        smbrs(iel)  = smbrs(iel)  - w1(iel)*rtpa(iel,ivar)
     endif

    enddo

  endif

endif

! --> Terme source pour l'enthalpie du solide

if ( ivar.ge.isca(ih2(1)) .and. ivar.le.isca(ih2(nclacp)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  numcla = ivar-isca(ih2(1))+1
  numcha = ichcor(numcla)
  ixchcl = isca(ixch(numcla))
  ixckcl = isca(ixck(numcla))
  ipcx2c = ipproc(ix2(numcla))
  ipcro2 = ipproc(irom2(numcla ))
  ipcdia = ipproc(idiam2(numcla))
  ipcte2 = ipproc(itemp2(numcla))
  ipcte1 = ipproc(itemp1)
  ipcght = ipproc(igmhet(numcla))
  if ( ihtco2 .eq. 1 ) then
    ipghc2 = ipproc(ighco2(numcla))
  endif

  ipcgd1 = ipproc(igmdv1(numcla))
  ipcgd2 = ipproc(igmdv2(numcla))
  ipcgch = ipproc(igmdch(numcla))
  ixchcl = isca(ixch(numcla))
  ixckcl = isca(ixck(numcla))

! ---- Calcul preliminaire : calcul de X2(NUMCLA) dans W11
!      ATTENTION tableau a conserver

  do iel = 1, ncel
! Rq : on utilise PROPCE(IEL,IPCX2C) car on veut X2 a l'iteration n
!      (en RTPA)
    w11(iel) = propce(iel,ipcx2c)
  enddo

! ---- Contribution aux bilans explicite et implicite
!        des echanges par diffusion moleculaire
!        6 Lambda Nu / diam**2 / Rho2 * Rho * (T1-T2)

! ------ Calcul de lambda dans W1

  xnuss = 2.d0
  do iel = 1, ncel
    if ( ivisls(ihm).gt.0 ) then
      ipcvsl = ipproc(ivisls(ihm))
      if ( icp.gt.0 ) then
        ipccp   = ipproc(icp)
        w1(iel) = propce(iel,ipcvsl) * propce(iel,ipccp)
      else
        w1(iel) = propce(iel,ipcvsl) * cp0
      endif
    else
      if ( icp.gt.0 ) then
        ipccp   = ipproc(icp)
        w1(iel) = visls0(ihm) * propce(iel,ipccp)
      else
        w1(iel) = visls0(ihm) * cp0
      endif
    endif
  enddo

! ------ Calcul du diametre des particules dans W2
!        On calcule le d20 = (A0.D0**2+(1-A0)*DCK**2)**0.5

  do iel = 1, ncel
    w2(iel) = ( xashch(numcha)*diam20(numcla)**2 +                &
                (1.d0-xashch(numcha))*propce(iel,ipcdia)**2       &
              )**0.5
  enddo

! ------ Contribution aux bilans explicite et implicite de
!        des echanges par diffusion moleculaire

  do iel = 1, ncel
    if ( w11(iel).gt.epsicp ) then
      aux         = 6.d0 * w1(iel) * xnuss / w2(iel)**2           &
                  / propce(iel,ipcro2) * propce(iel,ipcrom)       &
                  * w11(iel)                                      &
                  * volume(iel)
      rhovst      = aux / cp2ch(numcha) /w11(iel)

      smbrs(iel)  = smbrs(iel) - aux *                            &
                  ( propce(iel,ipcte2) - propce(iel,ipcte1) )
      rovsdt(iel) = rovsdt(iel) + max(zero,rhovst)
    endif

  enddo


! ---- Contribution aux bilans explicite et implicite
!        du terme echange d'energie entre les phases :
!        GAMA(dev1) H(mv1,T2)+GAMA(dev2) H(mv2,T2)

  do iel = 1, ncel

!        Gama Dev1 et Gama Dev2

    gamdv1 = propce(iel,ipcrom)*rtp(iel,ixchcl)                   &
            *propce(iel,ipcgd1)

    gamdv2 = propce(iel,ipcrom)*rtp(iel,ixchcl)                   &
            *propce(iel,ipcgd2)

!        H(mv1,T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo

    coefe(ichx1) = a1(numcha)*wmole(ichx1c(numcha))               &
  /( a1(numcha)*wmole(ichx1c(numcha))                             &
    +b1(numcha)*wmole(ico)                                        &
    +c1(numcha)*wmole(ih2o) )
    coefe(ico  ) = b1(numcha)*wmole(ico)                          &
   /( a1(numcha)*wmole(ichx1c(numcha))                            &
     +b1(numcha)*wmole(ico)                                       &
     +c1(numcha)*wmole(ih2o) )
    coefe(ih2o ) = c1(numcha)*wmole(ih2o)                         &
   /( a1(numcha)*wmole(ichx1c(numcha))                            &
     +b1(numcha)*wmole(ico)                                       &
     +c1(numcha)*wmole(ih2o) )

    t2         = propce(iel,ipcte2)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    f1mc(numcha) = 1.d0

    mode      = -1
    call cpthp1                                                   &
    !==========
    ( mode  , xhdev1    , coefe  , f1mc   , f2mc   ,              &
      t2    )

!        H(mv2,T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(ichx2) = a2(numcha)*wmole(ichx2c(numcha))               &
   /( a2(numcha)*wmole(ichx2c(numcha))                            &
     +b2(numcha)*wmole(ico)                                       &
     +c2(numcha)*wmole(ih2o) )
    coefe(ico  ) = b2(numcha)*wmole(ico)                          &
   /( a2(numcha)*wmole(ichx2c(numcha))                            &
     +b2(numcha)*wmole(ico)                                       &
     +c2(numcha)*wmole(ih2o) )
    coefe(ih2o ) = c2(numcha)*wmole(ih2o)                         &
   /( a2(numcha)*wmole(ichx2c(numcha))                            &
     +b2(numcha)*wmole(ico)                                       &
     +c2(numcha)*wmole(ih2o) )

    t2         = propce(iel,ipcte2)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    f2mc(numcha) = 1.d0

    mode      = -1
    call cpthp1                                                   &
    !==========
    ( mode  , xhdev2    , coefe  , f1mc   , f2mc   ,              &
      t2    )

!         Contribution aux bilans explicite et implicite

    smbrs(iel) = smbrs(iel)                                       &
                 + (gamdv1*xhdev1+gamdv2*xhdev2)*volume(iel)

  enddo

! ------ combustion heterogene : C(s) + 02 ---> 0.5 C0
!        GamHET * (28/12 H(CO,T2)-16/12 H(O2,T1) )

  do iel = 1, ncel

!        Calcul de HCO(T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(ico) = 1.d0
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo

    t2        = propce(iel,ipcte2)
    mode      = -1
    call cpthp1                                                   &
    !==========
    ( mode  , xhco    , coefe  , f1mc   , f2mc   ,                &
      t2    )

!        Calcul de HO2(T1)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(io2) = 1.d0
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo

    t1        = propce(iel,ipcte1)
    mode      = -1
    call cpthp1                                                   &
    !==========
    ( mode  , xho2    , coefe  , f1mc   , f2mc   ,                &
      t1    )

!         Contribution aux bilans explicite et implicite

    if ( rtpa(iel,ixckcl) .gt. epsicp ) then

      gamhet = propce(iel,ipcrom)*propce(iel,ipcght)              &
               * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +              &
                2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))      &
                 /(rtpa(iel,ixckcl))**(1.d0/3.d0) )

    else
      gamhet = 0.d0
    endif

   smbrs(iel) = smbrs(iel)                                        &
                 +gamhet                                          &
                  *(28.d0/12.d0*xhco-16.d0/12.d0*xho2)            &
                  *volume(iel)

  enddo

! ------ combustion heterogene : C(s) + C02 ---> 2 C0
!        GamHET * (56/12 H(CO,T2)-44/12 H(O2,T1) )

  if ( ihtco2 .eq. 1 ) then
    do iel = 1, ncel

!        Calcul de HCO(T2)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2        = propce(iel,ipcte2)
      mode      = -1
      call cpthp1                                                 &
      !==========
      ( mode  , xhco    , coefe  , f1mc   , f2mc   ,              &
        t2    )

!        Calcul de HCO2(T1)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico2) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t1        = propce(iel,ipcte1)
      mode      = -1
      call cpthp1                                                 &
      !==========
      ( mode  , xho2    , coefe  , f1mc   , f2mc   ,              &
        t1    )

!         Contribution aux bilans explicite et implicite

      if ( rtpa(iel,ixckcl) .gt. epsicp ) then

        gamhet = propce(iel,ipcrom)*propce(iel,ipghc2)            &
                 * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +            &
                2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))      &
                 /(rtpa(iel,ixckcl))**(1.d0/3.d0) )

      else
        gamhet = 0.d0
      endif

     smbrs(iel) = smbrs(iel)                                      &
                   +gamhet                                        &
                    *(56.d0/12.d0*xhco-44.d0/12.d0*xho2)          &
                    *volume(iel)

    enddo

  endif

!       --> Terme source sur H2 issu du sechage)

  if ( ippmod(icp3pl) .eq. 1 ) then

! ---- Contribution du TS interfacial aux bilans explicite et implicite


    ipcsec = ipproc(igmsec(numcla))
    ipcte2 = ipproc(itemp2(numcla))

    do iel = 1, ncel

!          Calcul de H(H2O) a T2

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ih2o) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2 = propce(iel,ipcte2)
      if ( t2 .gt. 100.d0 ) then
        t2 = 100.d0
      endif
      mode      = -1
      call cpthp1                                                 &
      !==========
      ( mode  , hh2ov    , coefe  , f1mc   , f2mc   ,             &
        t2    )

      hlv = 2.263d+6

!         Contribution aux bilans explicite

      if ( rtpa(iel,isca(ixwt(numcla))).gt. epsicp .and.          &
           xwatch(numcha) .gt. epsicp       ) then

        aux = -propce(iel,ipcrom)*propce(iel,ipcsec)              &
       *(rtp(iel,isca(ixwt(numcla)))/propce(iel,ipcx2c))          &
       *(1.d0                    /xwatch(numcha))                 &
             *hh2ov

      else
        aux = 0.d0
      endif

      smbrs(iel) = smbrs(iel) + aux*volume(iel)

    enddo

  endif

endif

!===============================================================================
! 3. PRISE EN COMPTE DES TERMES SOURCES POUR LES VARIABLES RELATIVES
!    AU MELANGE GAZEUX
!===============================================================================

! --> Terme source pour les matieres volatiles legeres

if ( ivar.ge.isca(if1m(1)) .and. ivar.le.isca(if1m(ncharb)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calcul de GMDEV1 = - SOMME (rho.XCH.GMDV1) > 0  --> W1

  numcha = ivar-isca(if1m(1))+1
  do iel = 1, ncel
    w1(iel) = zero
  enddo
  do icla = 1, nclacp
    ipcgd1 = ipproc(igmdv1(icla))
    ixchcl = isca(ixch(icla))
    if ( ichcor(icla).eq.numcha ) then
      do iel = 1, ncel
        w1(iel) = w1(iel) - propce(iel,ipcrom)*rtp(iel,ixchcl)    &
                * propce(iel,ipcgd1)
      enddo
    endif
  enddo

! ---- Contribution du TS interfacial aux bilans explicite et implicite

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
  enddo

endif


! --> Terme source pour les matieres volatiles lourdes

if ( ivar.ge.isca(if2m(1)) .and. ivar.le.isca(if2m(ncharb)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calcul de GMDEV2 = - SOMME (rho.XCH.GMDV2) >0 --> W1

  numcha = ivar-isca(if2m(1))+1
  do iel = 1, ncel
    w1(iel) = zero
  enddo
  do icla = 1, nclacp
    ipcgd2 = ipproc(igmdv2(icla))
    ixchcl = isca(ixch(icla))
    if ( ichcor(icla).eq.numcha ) then
      do iel = 1, ncel
        w1(iel) = w1(iel) - propce(iel,ipcrom)*rtp(iel,ixchcl)    &
                * propce(iel,ipcgd2)
      enddo
    endif
  enddo

! ---- Contribution du TS interfacial pour le bilan explicite

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
  enddo

endif


! --> Terme source pour le traceur 3 (O2) (C de la comb. het.)

if ( ivar.eq.isca(if3m) ) then

! RQ IMPORTANTE :  On prend les meme TS que pour Xck
!                  afin d'etre conservatif

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calcul prelimimaire pour la partie explicite : W1
!                          pour la partie implicite : W2

  do iel = 1, ncel
    w1(iel) = zero
  enddo

  do icla = 1, nclacp
    ipcght = ipproc(igmhet(icla))
    ixckcl = isca(ixck(icla))
    do iel = 1, ncel
      if ( rtpa(iel,ixckcl) .gt. epsicp ) then
        w1(iel) = w1(iel)                                         &
                 - propce(iel,ipcrom)*propce(iel,ipcght)          &
                 * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +            &
                    2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))  &
                    /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
      endif
    enddo
  enddo

! ---- Contribution du TS interfacial aux bilans explicite et implicite

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
  enddo

endif



! --> Terme source pour le traceur 3 (CO2) (C de la comb. het.)

if ( ihtco2 .eq. 1 ) then
  if ( ivar.eq.isca(if3mc2) ) then

! RQ IMPORTANTE :  On prend les meme TS que pour Xck
!                  afin d'etre conservatif

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

! ---- Calcul prelimimaire pour la partie explicite : W1
!                          pour la partie implicite : W2

    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp
      ipcght = ipproc(ighco2(icla))
      ixckcl = isca(ixck(icla))
      do iel = 1, ncel
        if ( rtpa(iel,ixckcl) .gt. epsicp ) then
          w1(iel) = w1(iel)                                       &
                   - propce(iel,ipcrom)*propce(iel,ipcght)        &
                   * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +          &
                      2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))&
                      /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
        endif
      enddo
    enddo

! ---- Contribution du TS interfacial aux bilans explicite et implicite

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
    enddo

  endif

endif

! --> Terme source pour la variance du traceur 4 (Air)

if ( ivar.eq.isca(if4p2m) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calcul des termes sources explicite et implicite
!      relatif aux echanges interfaciaux entre phases

 numtra = 4
 call cptsvi                                                      &
!!==========
 ( ncelet , ncel   , numtra ,                                     &
   rtp    , propce , volume ,                                     &
   smbrs  , rovsdt ,                                              &
   w1     , w2     ,                                              &
   w3 )

! ---- Calcul des termes sources explicite et implicite
!      relatif aux termes de production et de dissipation

!      Pointeur relatif au scalaire associe
!      (0 si pas de scalaire associe)
 iscala = 0

 call cptsvc                                                      &
!!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  , iscala ,          &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   smbrs  , rovsdt ,                                              &
   ra     )

endif

! --> Terme source pour le traceur 5 (Eau issue du séchage)

if ( ippmod(icp3pl) .eq. 1 ) then

  if ( ivar.eq.isca(if5m) ) then


    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

! ---- Contribution du TS interfacial aux bilans explicite et implicite

    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp

      ipcsec = ipproc(igmsec(icla))
      ipcx2c = ipproc(ix2(icla))
      numcha = ichcor(icla)

      do iel = 1, ncel

        if (  rtpa(iel,isca(ixwt(icla))).gt. epsicp               &
            .and.                                                 &
              xwatch(numcha) .gt. epsicp       ) then

          w1(iel) = w1(iel)                                       &
      + propce(iel,ipcrom)*propce(iel,ipcsec)                     &
       *(rtp(iel,isca(ixwt(icla)))/propce(iel,ipcx2c))            &
       *(1.d0                  /xwatch(numcha))

        endif

      enddo

    enddo

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
    enddo

  endif

endif

! --> Terme source pour CO2

if ( ieqco2 .ge. 1 ) then

  if ( ivar.eq.isca(iyco2) ) then


    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

! ---- Contribution du TS interfacial aux bilans explicite et implicite

! Oxydation du CO
! ===============

!  Dryer Glassman : XK0P en (moles/m3)**(-0.75) s-1

!          XK0P = 1.26D10

!           XK0P = 1.26D7 * (1.1)**(NTCABS)
!           IF ( XK0P .GT. 1.26D10 ) XK0P=1.26D10

!           T0P  = 4807.D0

!  Howard : XK0P en (moles/m3)**(-0.75) s-1

!             XK0P = 4.11D9
!             T0P  = 15090.D0

!  Westbrook & Dryer

    xk0p = 23.256d0
    t0p  = 20096.d0

! Dissociation du CO2 (Trinh Minh Chinh)
! ===================

!          XK0M = 5.D8
!          T0M  = 4807.D0

!          XK0M = 0.D0

!  Westbrook & Dryer

    xk0m = 20.03d0
    t0m  = 20096.d0

    err1mx = 0.d0
    err2mx = 0.d0

! Nombre d'iterations

    itermx = 500

! Nombre de points converges

   nbpauv = 0
   nbepau = 0
   nbrich = 0
   nberic = 0

! Precision pour la convergence

   errch = 1.d-8

   do iel = 1, ncel

     xxco  = propce(iel,ipproc(iym1(ico  )))/wmole(ico)           &
            *propce(iel,ipproc(irom1))
     xxo2  = propce(iel,ipproc(iym1(io2  )))/wmole(io2)           &
             *propce(iel,ipproc(irom1))
     xxco2 = propce(iel,ipproc(iym1(ico2 )))/wmole(ico2)          &
            *propce(iel,ipproc(irom1))
     xxh2o = propce(iel,ipproc(iym1(ih2o )))/wmole(ih2o)          &
            *propce(iel,ipproc(irom1))

     xxco  = max(xxco ,zero)
     xxo2  = max(xxo2 ,zero)
     xxco2 = max(xxco2,zero)
     xxh2o = max(xxh2o,zero)

     xco2mx=xxco2+min(xxco,2.d0*xxo2)

!            XKP  = XK0P*EXP(-T0P/PROPCE(IEL,IPPROC(ITEMP1)))
!            XKM  = XK0M*EXP(-T0M/PROPCE(IEL,IPPROC(ITEMP1)))

     xkp = exp(xk0p-t0p/propce(iel,ipproc(itemp1)))
     xkm = exp(xk0m-t0m/propce(iel,ipproc(itemp1)))

!           Recherche de l'état d'équilibre

!           Recerche itérative sans controle de convergence
!            (pour conserver la parallelisation sur les mailles)
!           sur le nombre de moles de reaction séparant
!           l'etat avant réaction (tel que calculé par Cpcym)
!           de l'état d'équilibre

!           initialisation par l'état transporté

     anmr  = xxco2
     xcot  = xxco + xxco2
     xo2t  = xxo2 + 0.5d0*xxco2
     sqh2o = sqrt(xxh2o)

     iterch = 0
 10        continue

       xco2e = anmr
       xcoe  = xcot-anmr
       xo2e  = xo2t-0.5d0*anmr
       iterch = iterch + 1

       ter1 = xkp*(xo2e**0.25d0)*sqh2o
       if ( abs(ter1) .gt. 1.d-15 ) then
         ddelta = (ter1*xcoe-xkm*xco2e)                           &
                 /(ter1*(1.d0+xcoe/(8.d0*xo2e))+xkm)
       else
         ddelta = 0.d0
       endif

       anmr = anmr + ddelta

!  Clipping
!       Cote Pauvre
       if ( xo2t .gt. 0.5d0*xcot ) then
         anmr = max(zero,min(xcot,anmr))
       else
!       Cote Riche
         anmr = max(zero,min(2.d0*xo2t,anmr))
       endif

! Test de convergence

       if ( abs(ddelta) .gt. errch                                &
           .and. iterch.lt.itermx       ) goto 10

! Calcul du nombre de points non convergé

     if ( iterch.eq.itermx                                        &
         .and. abs(ddelta) .gt. errch ) then
       nberic = nberic +1
     endif
     err1mx = max(err1mx,abs(ddelta))

     xco2e = anmr
     xcoe  = xcot-anmr
     xo2e  = xo2t-0.5d0*anmr

     xden = xkp*sqh2o*(0.5d0*(xo2t+xo2e))**0.25d0
     if ( xden .ne. 0.d0 ) then

       tauchi = 1.d0/xden
       tautur = rtpa(iel,ik)/rtpa(iel,iep)

       x2 = 0.d0
       do icla = 1, nclacp
         x2 = x2 + propce(iel,ipproc(ix2(icla)))
       enddo

       if ( ieqco2 .eq. 1 ) then

!              On transporte CO2

         smbrs(iel)  = smbrs(iel)                                 &
                      +wmole(ico2)/propce(iel,ipproc(irom1))      &
         * (xco2e-xxco2)/(tauchi+tautur)                          &
         * (1.d0-x2)                                              &
         * volume(iel) * propce(iel,ipcrom)

       else if ( ieqco2 .eq. 2 ) then

!              On transporte CO

         smbrs(iel)  = smbrs(iel)                                 &
                      +wmole(ico)/propce(iel,ipproc(irom1))       &
         * (xcoe-xxco)/(tauchi+tautur)                            &
         * (1.d0-x2)                                              &
         * volume(iel) * propce(iel,ipcrom)

       endif

       w1(iel) = volume(iel)*propce(iel,ipcrom)/(tauchi+tautur)
       rovsdt(iel) = rovsdt(iel) +   max(w1(iel),zero)

     else
       rovsdt(iel) = rovsdt(iel) + 0.d0
       smbrs(iel)  = smbrs(iel)  + 0.d0
     endif

   enddo

   if(irangp.ge.0) then
     call parcpt(nberic)
     call parmax(err1mx)
   endif
   write(NFECRA,*) ' Erreur max = ',ERR1MX
   write(NFECRA,*) ' Points non   ',NBERIC

!     Terme source : combustion heterogene par le CO2

   if ( ihtco2 .eq. 1) then

     do iel = 1, ncel

       aux = 0.d0
       do icla = 1,nclacp

         ixckcl = isca(ixck(icla))
         ipghc2 = ipproc(ighco2(icla))

         aux = aux                                                &
              + propce(iel,ipcrom)*propce(iel,ipghc2)             &
               *(rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

       enddo

       rovsdt(iel) = rovsdt(iel) - aux*(wmole(ico2)/0.012)

     enddo

   endif

  endif

endif

! Free memory
deallocate(w1, w2, w3)
deallocate(w4, w5, w6)
deallocate(w11)

!--------
! FORMATS
!--------

 1000 format(' TERMES SOURCES PHYSIQUE PARTICULIERE POUR LA VARIABLE '  &
       ,a8,/)

!----
! FIN
!----

return

end subroutine
