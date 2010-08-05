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

subroutine futssc &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iscal  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   izfppp , idevel , ituser , ia     ,                            &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt ,                                              &
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , w10    , w11    ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME FUEL
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
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iscal            ! i  ! <-- ! scalar number                                  !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! itypfb(nfabor    ! te ! --> ! type des faces de bord                         !
! nfml  ,nprfml    !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
! viscf(nfac)      ! tr ! --- ! tableau de travail    faces internes           !
! viscb(nfabor     ! tr ! --- ! tableau de travail    faces de bord            !
! xam(nfac,2)      ! tr ! --- ! tableau de travail    faces de bord            !
! w1..11(ncelet    ! tr ! --- ! tableau de travail    cellules                 !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
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
use fuincl
use ppincl
use ppcpfu

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp, numtra
integer          nideve , nrdeve , nituse , nrtuse
integer          iscal

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml) , itypfb(nfabor,nphas)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          izfppp(nfabor)
integer          idevel(nideve)
integer          ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smbrs(ncelet), rovsdt(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision xam(nfac,2)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet), w11(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          ivar , ipcrom, iel, iphas , icla , numcla
integer          iexp1 , iexp2 , iexp3 , ifac
integer          iscala
integer          ipcro2 , ipcte1 , ipcte2
integer          ipcdia , ipcvst
integer          imode  , iesp
integer          ipcgev , ipcght , ipchgl
integer          itermx,nbpauv,nbrich,nbepau,nberic

double precision aux, rhovst
double precision t1, h1, t2, h2
double precision xng,rhofol , rom
double precision gameva,fdev,fsd,ftrac,hfov
double precision ho2,hco,xesp(ngazem),xcoke,t2mt1
double precision gmech,gmvap,gmhet

double precision xxco,xxo2,xxco2,xxh2o,xco2mx
double precision xkp,xk0p,xkm,xk0m,wplus,wmoins,t0p,t0m
double precision auxp,auxm, aux1 , aux2 , aux3

double precision xeq,anmr,xcot,xo2t,xco2e,xo2e,xcoe,tauchi,tautur
double precision sqh2o , x2 , wmhcn , wmno

integer          iterch
double precision err1mx,err2mx,anmr0

double precision errch,ter1,ddelta,fn,qpr
double precision auxmax,auxmin
double precision ymoy,volm,volmp,dmp

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
chaine = nomvar(ipprtp(ivar))

! --- Numero de phase associee au scalaire ISCAL
iphas = iphsca(iscal)

! --- Numero des grandeurs physiques (voir usclim)
ipcrom = ipproc(irom(iphas))
ipcvst = ipproc(ivisct(iphas))

! --- Temperature phase gaz

ipcte1 = ipproc(itemp1)

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES POUR LES VARIABLES RELATIVES
!    AUX CLASSES DE PARTICULES
!===============================================================================

! --> Terme source pour l'enthalpie du liquide

if ( ivar .ge. isca(ihlf(1))     .and.                            &
     ivar .le. isca(ihlf(nclafu))      ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  numcla = ivar-isca(ihlf(1))+1

  ipcro2 = ipproc(irom3 (numcla))
  ipcdia = ipproc(idiam3(numcla))
  ipcte2 = ipproc(itemp3(numcla))
  ipcgev = ipproc(igmeva(numcla))
  ipcght = ipproc(igmhtf(numcla))
  ipchgl = ipproc(ih1hlf(numcla))

!       La variable est l'enthalpie du liquide ramenée à
!       la masse de mélange
!       Les flux interfaciaux contribuent à faire varier l'enthalpie
!       du liquide
!       La vapeur emporte son enthalpie
!       flux = PROPCE(IEL,IPPROC(IGMEVA))
!       enthalpie massique reconstituée à partir de EHGAZE(IFOV )
!       à la température de la goutte
!       L'oxydation héterogène comporte un flux "entrant" d'O2
!       un flux sortant de CO
!       le flux net est celui du carbone
!       fluxIN  = 16/12 * PROPCE(IEL,IPPROC(IGMHTF))
!       fluxOUT = 28/12 * PROPCE(IEL,IPPROC(IGMHTF))
!       Enthalpie entrante reconstituée à partir de EHGAZE(IO2 )
!       à la température du gaz environnant
!       Enthalpie sortante reconstituée à partir de EHGAZE(ICO )
!       à la température du grain

  imode = -1
  do iel = 1, ncel

    rom = propce(iel,ipcrom)

    do iesp = 1, ngazem
      xesp(iesp) = zero
    enddo

    xesp(ifov) = 1.d0
    call futhp1(imode,hfov,xesp,propce(iel,ipcte2))

    xesp(ifov) = zero
    xesp(io2)  = 1.d0
    call futhp1(imode,ho2 ,xesp,propce(iel,ipcte1))

    xesp(io2)  = zero
    xesp(ico)  = 1.d0
    call futhp1(imode,hco ,xesp,propce(iel,ipcte2))

    t2mt1 = propce(iel,ipcte2)-propce(iel,ipcte1)

    gmech = -propce(iel,ipchgl)*t2mt1
    gmvap = propce(iel,ipcgev)*hfov*t2mt1
    gmhet = 16.d0/12.d0*propce(iel,ipcght)*ho2                    &
           -28.d0/12.d0*propce(iel,ipcght)*hco

    smbrs(iel) = smbrs(iel) +                                     &
         ( gmech+gmvap+gmhet )*rom*volume(iel)

    rhovst = ( propce(iel,ipchgl)                                 &
              -propce(iel,ipcgev)*hfov )/cp2fol                   &
            *rom*volume(iel)

    rovsdt(iel) = rovsdt(iel) +  max(zero,rhovst)

  enddo

! --> T.S. pour la masse de liquide

elseif ( ivar .ge. isca(iyfol(1))     .and.                       &
         ivar .le. isca(iyfol(nclafu))        ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  numcla = ivar-isca(iyfol(1))+1

  ipcro2 = ipproc(irom3 (numcla))
  ipcdia = ipproc(idiam3(numcla))
  ipcte2 = ipproc(itemp3(numcla))
  ipcgev = ipproc(igmeva(numcla))
  ipcght = ipproc(igmhtf(numcla))
  ipchgl = ipproc(ih1hlf(numcla))

  do iel = 1, ncel

    t2mt1 =  propce(iel,ipcte2)-propce(iel,ipcte1)
    gmvap = -propce(iel,ipcgev)*t2mt1
    gmhet = -propce(iel,ipcght)

    smbrs(iel) = smbrs(iel)                                       &
         - propce(iel,ipcrom)*volume(iel)*(gmvap+gmhet)
!          ROVSDT(IEL) = ROVSDT(IEL) + ZERO

  enddo

! --> T.S. pour le traceur de la vapeur

elseif ( ivar .eq. isca(ifvap) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  do icla = 1, nclafu

    ipcte2 = ipproc(itemp3(icla))
    ipcgev = ipproc(igmeva(icla))

    do iel = 1, ncel

      t2mt1 = propce(iel,ipcte2)-propce(iel,ipcte1)
      gmvap = -propce(iel,ipcgev)*t2mt1

      smbrs(iel) = smbrs(iel)                                     &
                 + gmvap*propce(iel,ipcrom)*volume(iel)
!           ROVSDT(IEL) = ROVSDT(IEL) + ZERO

    enddo

  enddo

! --> T.S. pour le traceur du C ex réaction heterogene

elseif ( ivar .eq. isca(ifhtf) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  do icla = 1, nclafu

    ipcght = ipproc(igmhtf(icla))

    do iel = 1, ncel

      smbrs(iel) = smbrs(iel)                                     &
           -propce(iel,ipcrom)*propce(iel,ipcght)*volume(iel)
!           ROVSDT(IEL) = ROVSDT(IEL) + ZERO

    enddo

  enddo

endif

! --> Terme source pour la variance du traceur 4 (Air)

if ( ivar.eq.isca(if4p2m) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calcul des termes sources explicite et implicite
!      relatif aux echanges interfaciaux entre phases

 numtra = 4
 call futsvi                                                      &
!!==========
 ( ncelet , ncel   , numtra ,                                     &
   rtp    , propce , volume ,                                     &
   smbrs  , rovsdt , w1      )

! ---- Calcul des termes sources explicite et implicite
!      relatif aux termes de production et de dissipation

!      Pointeur relatif au scalaire associe
!      (0 si pas de scalaire associe)

  iscala = 0

  call futsvc                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iscal  , iscala ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   smbrs  , rovsdt ,                                              &
   viscb  ,                                                       &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     ,                                     &
   rdevel , rtuser , ra     )

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

    xk0p = 23.256
    t0p  = 20096.

! Dissociation du CO2 (Trinh Minh Chinh)
! ===================

!          XK0M = 5.D8
!          T0M  = 4807.D0

!          XK0M = 0.D0

!  Westbrook & Dryer

    xk0m = 20.03
    t0m  = 20096.

    err1mx = 0.d0
    err2mx = 0.d0

! Nombre d'iterations

    itermx = 5000

! Nombre de points converges

   nbpauv = 0
   nbepau = 0
   nbrich = 0
   nberic = 0

! Precision pour la convergence

   errch = 1.d-8

   do iel = 1, ncel

     xxco  = propce(iel,ipproc(iym1(ico  )))/wmole(ico)           &
            *propce(iel,ipcrom)
     xxo2  = propce(iel,ipproc(iym1(io2  )))/wmole(io2)           &
            *propce(iel,ipcrom)
     xxco2 = propce(iel,ipproc(iym1(ico2 )))/wmole(ico2)          &
            *propce(iel,ipcrom)
     xxh2o = propce(iel,ipproc(iym1(ih2o )))/wmole(ih2o)          &
            *propce(iel,ipcrom)

     xxco  = max(xxco ,zero)
     xxo2  = max(xxo2 ,zero)
     xxco2 = max(xxco2,zero)
     xxh2o = max(xxh2o,zero)

     xco2mx=xxco2+min(xxco,2.d0*xxo2)

!            XKP  = XK0P*EXP(-T0P/PROPCE(IEL,IPPROC(ITEMP1)))
!            XKM  = XK0M*EXP(-T0M/PROPCE(IEL,IPPROC(ITEMP1)))
!           if ( iel .eq. 36 ) then
     xkp = exp(xk0p-t0p/propce(iel,ipproc(itemp1)))
     xkm = exp(xk0m-t0m/propce(iel,ipproc(itemp1)))
     xeq = exp(xk0m-xk0p)

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

     xo2e = xo2t-0.5d0*anmr

!    Si presence de O2 et H2O

     if ( xo2e .gt. 0.d0 .and. sqh2o .gt. 0.d0 ) then

         if ( iel.eq.36 .or. iel.eq.15151) then
           write(*,*) ' TEST  = ',SQH2O,SQH2O .GT. 0.D0
         endif

       iterch = 0
 10          continue

         xco2e = anmr
         xcoe = xcot-anmr
         xo2e = xo2t-0.5d0*anmr
         iterch = iterch + 1

         if ( iel.eq. 36 .or. iel.eq.15151) then
           write(*,*) ' ITERCH  = ',ITERCH
         endif

         ter1 = xkp*sqh2o

         ddelta = (ter1*xcoe*xo2e-xkm*xco2e*(xo2e**0.75d0))       &
                 /(ter1*(xo2e+xcoe/8.d0)+xkm*(xo2e**0.75d0))


         anmr = anmr + ddelta

!  Clipping
!       Cote Pauvre
        if ( xo2t .gt. 0.5d0*xcot ) then
          anmr = max( zero, min(xcot,anmr))
        else
!       Cote Riche
          anmr = max( zero, min(2.d0*xo2t,anmr))
        endif

! Test de convergence

        if ( abs(ddelta) .gt. errch                               &
             .and. iterch.lt.itermx       ) goto 10

! Calcul du nombre de points non convergé

        if ( iterch.eq.itermx                                     &
            .and. abs(ddelta) .gt. errch ) then

          nberic = nberic +1

        endif
 20          continue
       err1mx = max(err1mx,abs(ddelta))

       xco2e = anmr
       xcoe  = xcot-anmr
       xo2e  = xo2t-0.5d0*anmr

       tauchi = 1.d0/(xkp*sqh2o*(0.5d0*(xo2t+xo2e))**0.25d0)
       tautur = rtpa(iel,ik(iphas))/rtpa(iel,iep(iphas))

       x2 = 0.d0
       do icla = 1,nclafu
         x2 = x2 + rtp(iel,isca(iyfol(icla)))
       enddo

!            On transporte CO2

       smbrs(iel)  = smbrs(iel)                                   &
                   +wmole(ico2)/propce(iel,ipcrom)                &
        * (xco2e-xxco2)/(tauchi+tautur) * volume(iel)             &
        * (1.d0-x2)

       w1(iel) = volume(iel)/(tauchi+tautur)
       rovsdt(iel) = rovsdt(iel) +   max(w1(iel),zero)

     else

!     Pas de CO

       smbrs(iel)  = smbrs(iel) + 0.d0
       rovsdt(iel) = rovsdt(iel)+ 0.d0

     endif

   enddo

   if(irangp.ge.0) then
     call parcpt(nberic)
     call parmax(err1mx)
   endif
   write(NFECRA,*) ' Erreur max = ',ERR1MX
   write(NFECRA,*) ' Points non   ',NBERIC

  endif

endif


! --> Terme source pour HCN et NO : uniquement a partir de la 2eme
!                                   iter

if ( ieqnox .eq. 1 .and. ntcabs .gt. 1) then

  if ( ivar.eq.isca(iyhcn) .or. ivar.eq.isca(iyno) ) then

    iexp1  = ipproc(ighcn1)
    iexp2  = ipproc(ighcn2)
    iexp3  = ipproc(ignoth)

! QPR= %N libéré pendant l'evaporation/taux de matieres volatiles
!          moyen

    qpr = 0.5d0

! YMOY = % vapeur en sorties

    ymoy  = 0.d0
    volm  = 0.d0
    dmp   = 0.d0
    volmp = 0.d0
    do ifac=1,nfabor

      if ( itypfb(ifac,iphas).eq. isolib ) then

        iel = ifabor(ifac)

        ymoy = ymoy + rtpa(iel,isca(ifvap))*volume(iel)
        volm = volm + volume(iel)
        do icla=1,nclafu
          dmp  = dmp  + rtpa(iel,isca(ing(icla)))                 &
                       *volume(iel)*dinifl(icla)
          volmp = volmp + volume(iel)
        enddo

      endif
    enddo

    if ( irangp .ge. 0 ) then
      call parsom(ymoy)
      call parsom(volm)
      call parsom(dmp)
      call parsom(volmp)
    endif

    if ( ymoy .gt. 0.d0 .and. volm .gt. 0.d0 ) then
      ymoy = ymoy/volm
      if ( dmp .gt. 0.d0 .and. volmp .gt.0.d0 ) then
        ymoy = ymoy/(dmp/volmp)
        WRITE(NFECRA,*) ' Taux YMOY = ',YMOY
      else
        write(nfecra,2000)
        ymoy = 0.5d0
      endif
    else
      write(nfecra,2000)
      ymoy = 0.5d0
    endif

! Azote dans le fuel

    fn = 0.015

! Masse molaire

    wmhcn = 0.027d0
    wmno  = 0.030d0

    if ( ivar.eq.isca(iyhcn) ) then

!        Terme source HCN

      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      auxmin = 1.d+20
      auxmax =-1.d+20

      do iel=1,ncel

        aux = volume(iel)*propce(iel,ipcrom)                      &
             *(propce(iel,iexp1)+propce(iel,iexp2))               &
             *rtpa(iel,isca(iyno))                                &
             *propce(iel,ipproc(immel))                           &
             /wmole(io2)

        smbrs(iel)  = smbrs(iel)                                  &
                     -aux*rtpa(iel,ivar)
        rovsdt(iel) = rovsdt(iel) + aux

        gmvap = 0.d0
        gmhet = 0.d0
        do icla=1,nclafu

          ipcgev = ipproc(igmeva(icla))
          ipcght = ipproc(igmhtf(icla))
          ipcte2 = ipproc(itemp3(icla))
          ipcte1 = ipproc(itemp1)

          gmvap = gmvap                                           &
                 + propce(iel,ipcrom)*volume(iel)                 &
                  *propce(iel,ipcgev)                             &
                  *(propce(iel,ipcte2)-propce(iel,ipcte1))

          gmhet = gmhet                                           &
                 +propce(iel,ipcrom)*volume(iel)                  &
                                    *propce(iel,ipcght)

        enddo

        aux = -fn*wmhcn/(wmole(in2)/2.d0)                         &
                 *(qpr*gmvap+(1.d0-qpr*ymoy)/ymoy*gmhet)
        smbrs(iel)  = smbrs(iel) + aux
        auxmin = min(auxmin,aux)
        auxmax = max(auxmax,aux)

      enddo
      write(*,*) ' AUX MIN MAX = ',AUXMIN,AUXMAX

    endif

    if ( ivar.eq.isca(iyno) ) then

!        Terme source NO

      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel=1,ncel

        aux1 = volume(iel)*propce(iel,ipcrom)                     &
              *propce(iel,iexp1)*rtpa(iel,isca(iyhcn))            &
              *propce(iel,ipproc(immel))/wmhcn
        aux2 = volume(iel)*propce(iel,ipcrom)                     &
              *propce(iel,iexp2)*rtpa(iel,isca(iyhcn))            &
              *wmno/wmhcn
        aux3 = volume(iel)*propce(iel,ipcrom)**1.5d0              &
              *propce(iel,iexp3)                                  &
              *(propce(iel,ipproc(iym1(io2))))*0.5d0              &
              * propce(iel,ipproc(iym1(in2)))

        smbrs(iel)  = smbrs(iel) - aux1*rtpa(iel,ivar)            &
                               + aux2 + aux3
        rovsdt(iel) = rovsdt(iel) + aux1
      enddo

    endif

  endif

endif

!--------
! FORMATS
!--------

 1000 format(' TERMES SOURCES PHYSIQUE PARTICULIERE POUR LA VARIABLE '  &
       ,a8,/)
 2000 format(' ATTENTION : LE TAUX DE VAPEUR EN SORTIE',                &
       ' EST NULLE',/,                                      &
       ' ON PREND QPR = 0.5')

!----
! FIN
!----

return

end subroutine
