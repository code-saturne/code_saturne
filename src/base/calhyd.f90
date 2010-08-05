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

subroutine calhyd &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , iphas  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   indhyd ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   fextx  , fexty  , fextz  ,                                     &
   dfextx , dfexty , dfextz ,                                     &
   phydr  , flumas , flumab ,                                     &
   coefa  , coefb  ,                                              &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbr   ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , w10    ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION D'UNE EQUATION DE POISSON SUR LA PRESSION HYDROSTATIQUE
!                DIV( GRAD(P) ) = DIV( F )
!                     ----             -
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
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iphas            ! i  ! <-- ! phase number                                   !
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
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! indhyd           ! e  ! --> ! indicateur de mise a jour de phydr             !
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
! phydr(ncelet)    ! tr ! <-- ! increment de pression hydrostatique            !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor)        !    !     !    faces de bord                               !
! fextx,fexty      ! tr ! <-- ! force exterieure generant la pression          !
! fextz(ncelet)    !    !     !    hydrostatique                               !
! dfextx,dfexty    ! tr ! <-- ! increment de force exterieure                  !
! dfextz(ncelet    !    !     !    generant la pression hydrostatique          !
! viscf(nfac)      ! tr ! --- ! 1*surface/dist aux faces internes              !
! viscb(nfabor     ! tr ! --- ! 1*surface/dist aux faces de bord               !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr  (ncelet    ! tr ! --- ! tableau de travail pour sec mem                !
! w1..10(ncelet    ! tr ! --- ! tableau de travail                             !
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

include "paramx.f90"
include "numvar.f90"
include "entsor.f90"
include "cstnum.f90"
include "optcal.f90"
include "period.f90"
include "parall.f90"
include "mltgrd.f90"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse , iphas


integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          indhyd
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision fextx(ncelet),fexty(ncelet),fextz(ncelet)
double precision dfextx(ncelet),dfexty(ncelet),dfextz(ncelet)
double precision phydr(ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision coefa(nfabor), coefb(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet)
double precision smbr(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

character*80     chaine
integer          lchain
integer          idebia, idebra
integer          iccocg, inc   , init  , isym  , ipol  , isqrt
integer          iel   , ical
integer          ireslp, nswmpr
integer          isweep, niterf, icycle
integer          iphydp
integer          nswrgp, imligp, iwarnp
integer          ipriph
integer          iinvpe
integer          idiffp, iconvp, ndircp
integer          nitmap, imgrp , ncymap, nitmgp
integer          iagmax, nagmax, npstmg
double precision residu, rnorm , rnrmf , rnrmdf
double precision epsrgp, climgp, extrap, epsilp
double precision precre, precab, thetap

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! --- Memoire
idebia = idbia0
idebra = idbra0


! --- Variables
ipriph = ipr(iphas)



! --- Options de resolution
!     Symetrique
!     Preconditionnement diagonal par defaut
isym  = 1
if (iresol(ipriph).eq.-1) then
  ireslp = 0
  ipol   = 0
else
  ireslp = mod(iresol(ipriph),1000)
  ipol   = (iresol(ipriph)-ireslp)/1000
endif

isqrt = 1

!     TEST DE VARIATION DE LA PRESSION HYDROSTATIQUE EN SORTIE


!     on regarde si les terme source ont varie
!     on ne passe dans calhyd que si on a des faces de sortie std
!     la precision pour les tests est a peu pres arbitraire.
precre = sqrt(epzero)
precab = 1.d2*epzero

ical = 0
do iel = 1, ncel
  rnrmf  = fextx(iel)**2+fexty(iel)**2+fextz(iel)**2
  rnrmdf = dfextx(iel)**2+dfexty(iel)**2+dfextz(iel)**2
  if ((rnrmdf.ge.precre*rnrmf).and.(rnrmdf.ge.precab)) then
    ical = 1
  endif
enddo
if (irangp.ge.0) then
  call parcpt (ical)
endif
if (ical.eq.0) then
  do iel = 1,ncel
    phydr(iel) = 0.d0
  enddo
  indhyd = 0
  return
endif

if ( mod(ntcabs,ntlist).eq.0 .or. iwarni(iu(1)) .ge.0 )           &
     write(nfecra,1000)

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'  Calcul de la pression hydrostatique : ',/,               &
'         mise a jour des Dirichlets en sortie (CALHYD)',/)

#else

 1000 format(                                                           &
'  Hydrostatic pressure computation: ',/,                   &
'         updating the Dirichlets at the end (CALHYD)',/)

#endif

indhyd = 1

!===============================================================================
! 2.  PREPARATION DE LA MATRICE DU SYSTEME A RESOUDRE
!===============================================================================

! ---> TERME INSTATIONNAIRE

do iel = 1, ncel
  w1(iel) = 0.d0
enddo

! ---> "VITESSE" DE DIFFUSION FACETTE

do iel = 1, ncel
  w10(iel) = 1.d0
enddo

call viscfa                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse , imvisf ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w10    ,                                                       &
   viscf  , viscb  ,                                              &
   rdevel , rtuser , ra     )


iconvp = 0
idiffp = 1
!  On resout avec des CL de flux nul partout
ndircp = 0

thetap = 1.d0
call matrix                                                       &
!==========
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp ,                                     &
   isym   , nfecra ,                                              &
   thetap ,                                                       &
   ifacel , ifabor ,                                              &
   coefb  , w1     ,                                              &
   flumas , flumab , viscf  , viscb  ,                            &
   dam    , xam    )


!===============================================================================
! 4.  INITIALISATION DU FLUX DE MASSE
!===============================================================================


!     PROJECTION AUX FACES DES TERMES SOURCES
init   = 1
inc    = 0
iccocg = 1
nswrgp = nswrgr(ipriph)
imligp = imligr(ipriph)
iwarnp = iwarni(ipriph)
epsrgp = epsrgr(ipriph)
climgp = climgr(ipriph)

call projts                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dfextx , dfexty , dfextz ,                                     &
   coefb     ,                                                    &
   flumas, flumab ,                                               &
   viscf  , viscb  ,                                              &
   w10    , w10    , w10    ,                                     &
   rdevel , rtuser , ra     )

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
            ifacel,ifabor,flumas,flumab,w7)
call prodsc(ncelet,ncel,isqrt,w7,w7,rnorm)


!===============================================================================
! 5.  PREPARATION DU MULTIGRILLE ALGEBRIQUE
!===============================================================================

if (imgr(ipriph).gt.0) then

!   --- Creation de la hierarchie de maillages

  CHAINE = 'PresHydr'
  iwarnp = iwarni(ipriph)
  iagmax = iagmx0(ipriph)
  nagmax = nagmx0(ipriph)
  npstmg = ncpmgr(ipriph)
  lchain = 8

  call clmlga                                                     &
  !==========
 ( chaine(1:8) ,     lchain ,                                     &
   ncelet , ncel   , nfac   ,                                     &
   isym   , iagmax , nagmax , npstmg , iwarnp ,                   &
   ngrmax , ncegrm ,                                              &
   dam    , xam    )

endif

!===============================================================================
! 6.  BOUCLES SUR LES NON ORTHOGONALITES (RESOLUTION)
!===============================================================================

! --- Nombre de sweeps
nswmpr = nswrsm(ipriph)

! --- Mise a zero des variables
!       RTP(.,IPR) sera l'increment de pression cumule
!       DRTP       sera l'increment d'increment a chaque sweep
!       W7         sera la divergence du flux de masse predit
do iel = 1,ncel
  phydr(iel) = 0.d0
  drtp(iel) = 0.d0
  smbr(iel) = 0.d0
enddo


! --- Boucle de reconstruction : debut
do isweep = 1, nswmpr

! --- Mise a jour du second membre
!     (signe "-" a cause de celui qui est implicitement dans la matrice)
  do iel = 1, ncel
    smbr(iel) = - w7(iel) - smbr(iel)
  enddo

! --- Test de convergence du calcul

  call prodsc(ncelet,ncel,isqrt,smbr,smbr,residu)
  if (iwarni(ipriph).ge.2) then
     CHAINE = 'PresHydr'
     write(nfecra,1400)chaine(1:8),isweep,residu
  endif

!MO IL FAUDRA VERIFIER LA PERTINENCE DU TEST

  if( residu .le. 10.d0*epsrsm(ipriph)*rnorm ) then
!     Si convergence,  sortie

    goto 101

  endif

! --- Resolution implicite sur l'increment d'increment DRTP
  do iel = 1, ncel
    drtp(iel) = 0.d0
  enddo

  CHAINE = 'PresHydr'
  nitmap = nitmax(ipriph)
  imgrp  = imgr  (ipriph)
  ncymap = ncymax(ipriph)
  nitmgp = nitmgf(ipriph)
  iwarnp = iwarni(ipriph)
  epsilp = epsilo(ipriph)
  iinvpe = 1

  call invers                                                     &
  !==========
 ( chaine(1:8)     , idebia , idebra ,                            &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   isym   , ipol   , ireslp , nitmap , imgrp  ,                   &
   ncymap , nitmgp ,                                              &
   iwarnp , nfecra , niterf , icycle , iinvpe ,                   &
   epsilp , rnorm  , residu ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dam    , xam    , smbr   , drtp   ,                            &
   w3     , w4     , w5     , w6     , w8     , w9     ,          &
   rdevel , rtuser , ra     )


  if( isweep.eq.nswmpr ) then
!     Mise a jour de l'increment de pression
    do iel = 1, ncel
      phydr(iel) = phydr(iel) + drtp(iel)
    enddo


  else

! --- Si ce n'est pas le dernier sweep
!       Mise a jour de l'increment de pression et calcul direct de la
!       partie en gradient d'increment de pression du second membre
!       (avec reconstruction)

    do iel = 1, ncel
      phydr(iel) = phydr(iel) + drtp(iel)
    enddo

    iccocg = 1
    init = 1
    inc  = 1
    nswrgp = nswrgr(ipriph)
    imligp = imligr(ipriph)
    iwarnp = iwarni(ipriph)
    epsrgp = epsrgr(ipriph)
    climgp = climgr(ipriph)
    extrap = 0.d0
    iphydp = 1

    call itrgrp                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dfextx , dfexty , dfextz ,                                     &
   phydr  ,                                                       &
   coefa  , coefb  ,                                              &
   viscf  , viscb  ,                                              &
   w10         , w10         , w10         ,                      &
   smbr   ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

  endif

enddo
! --- Boucle de reconstruction : fin

if(iwarni(ipriph).ge.2) then
   CHAINE = 'PresHydr'
   write( nfecra,1600)chaine(1:8),nswmpr
endif

 101  continue


!===============================================================================
! 7.  SUPPRESSION DE LA HIERARCHIE DE MAILLAGES
!===============================================================================

if (imgr(ipriph).gt.0) then
  CHAINE = 'PresHydr'
  lchain = 8
  call dsmlga(chaine(1:8), lchain)
  !==========
endif

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1400 format(1X,A8,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6)
 1600 format(                                                           &
'@                                                            ',/,&
'@ @@ ATTENTION : ',A8 ,' ETAPE DE PRESSION HYDROSTATIQUE     ',/,&
'@    =========                                               ',/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,&
'@                                                            '  )

#else

 1400 format(1X,A8,' : SWEEP = ',I5,' RIGHT HAND SIDE NORM = ',E14.6)
 1600 format(                                                           &
'@                                                            ',/,&
'@ @@ WARNING: ',A8 ,' HYDROSTATIC PRESSURE STEP              ',/,&
'@    ========                                                ',/,&
'@  Maximum number of iterations ',I10   ,' reached           ',/,&
'@                                                            '  )

#endif

!----
! FIN
!----

return

end subroutine
