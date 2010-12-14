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

subroutine fuphyv &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ibrom  , izfppp ,                                              &
   idevel , ituser , ia     ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : COMBUSTION FUEL

! Calcul de RHO du melange


! ATTENTION :
! =========


! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)


!  Il FAUT AVOIR PRECISE ICP(IPHAS) = 1
!     ==================
!    dans usini1 si on souhaite imposer une chaleur specifique
!    CP variable pour la phase IPHAS (sinon: ecrasement memoire).


!  Il FAUT AVOIR PRECISE IVISLS(Numero de scalaire) = 1
!     ==================
!     dans usini1 si on souhaite une diffusivite VISCLS variable
!     pour le scalaire considere (sinon: ecrasement memoire).




! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usini1 :
!             . la masse volumique (initialisee a RO0(IPHAS))
!             . la viscosite       (initialisee a VISCL0(IPHAS))
!      - dans usppiv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a la valeur donnee dans usiniv)

! On peut donner ici les lois de variation aux cellules
!     - de la masse volumique                      ROM    kg/m3
!         (et eventuellememt aux faces de bord     ROMB   kg/m3)
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)


! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)


! Il est conseille de ne garder dans ce sous programme que
!    le strict necessaire.



! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! nphmx            ! e  ! <-- ! nphsmx                                         !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!   (nphmx   )     !    !     !                                                !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
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
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use fuincl
use ppincl
use ppcpfu
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse , nphmx

integer          ibrom(nphmx)
integer          izfppp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ntbfui, ifuwi, ntbfur, ifuwr
integer          ntbwoi, iwori, ntbwor, iworr
integer          ifinia, ifinra
integer          iel, iphas, ipcrom, ipcro2 , ipcte1
integer          izone, ifac , icla
integer          ipbrom, iromf
double precision qtotz
double precision x1sro1, x2sro2, srrom1, uns1pw
double precision x2tot, wmolme, unsro1
double precision x2h2

integer       ipass
data          ipass /0/
save          ipass

!===============================================================================
!===============================================================================
! 0. ON COMPTE LES PASSAGES
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0

! --- Initialisation des tableaux de travail

do iel = 1, ncel
  w1(iel) = zero
  w2(iel) = zero
  w3(iel) = zero
  w4(iel) = zero
  w5(iel) = zero
  w6(iel) = zero
  w7(iel) = zero
  w8(iel) = zero
enddo

!     Pointeur sur masse volumique du gaz aux cellules
iromf = ipproc(irom1)

!===============================================================================
! 2. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE DISPERSEE
!                    VALEURS CELLULES
!                    ----------------
!    FRACTION MASSIQUE DE LIQUIDE
!    DIAMETRE
!    MASSE VOLUMIQUE
!===============================================================================

call fuphy2                                                       &
!==========
 ( ncelet , ncel   ,                                              &
   rtp    , propce )

!===============================================================================
! 3. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE GAZEUSE
!                    VALEURS CELLULES
!                    ----------------
!    TEMPERATURE
!    MASSE VOLUMIQUE
!    CONCENTRATIONS DES ESPECES GAZEUSES
!===============================================================================

! --- Calcul de l'enthalpie du gaz     dans W8 si transport de H2
!                           du melange         si pas de transport de H2
!            de F1M                    dans W2
!            de F2M                    dans W3
!            de F3M                    dans W4
!            de F4M                    dans W5
!            de F3P2M                  dans W6
!            de F4P2M                  dans W7

! ---- W1 = - Somme des X2(i)

do iel = 1, ncel

  w1(iel) = 0.d0
  do icla = 1, nclafu
    w1(iel) = w1(iel) - rtp(iel,isca(iyfol(icla)))
  enddo
  uns1pw  = 1.d0 / ( 1.d0 + w1(iel) )
  w2(iel) =  rtp(iel,isca(ifvap))  * uns1pw
  w4(iel) =  rtp(iel,isca(ifhtf))  * uns1pw
  w6(iel) =  rtp(iel,isca(if4p2m)) * uns1pw

! PPl 09 08 2006
!     Les tableaux de travail contiennent les grandeurs massiques de la
!     phase gaz
!     Les grandeurs variances de F1 et F3 ne sont pas utilisées pour
!     la pdf passant par F4, y placer plutot la variance de F4

!      W6(IEL) = RTP(IEL,ISCA(IF4P2M)) * UNS1PW
! PPl 161205
!       Attention, contrairement au cas charbon, la 2° enthalpie transportée
!       dans le cas du fuel est celle de l'ensemble liquide + vapeur
!       dans les conditions d'ébullition.
!       La reconstitution de l'enthalpie du gaz seul est donc fausse ...

!       X1 * H1 = HM - HLF + FVAP * HrefEvap
!       X2 * H2 = HLF - FVAP * HrefEvap
!       où X1 et X2 sont les fraction s massiques des deux phases
!       H1 l'enthalpie massique de la phase continue
!       H2 l'enthalpie massique de la phase dispersée
!       HM l'enthalpie massique du mélange
!       HLF l'enthalpie du liquide et de la vapeur
!           (ramenée à la masse de mélange)
!       FVAP la fraction massique du traceur attaché à la vapeur
!       HrefVap l'enthalpie massique de la vapeur dans les
!            conditons moyennes d'évaporation soit
!                      0.5*(TEVAP1+Min(Tevap2,Tliqu))

!            TEBMOY = 0.5D0 * ( TEVAP1
!    &              + MIN( PROPCE(IEL,IPPROC(ITEMP3)) , TEVAP2 ) )
!            EH2 = ( RTP(IEL,ISCA(IHLF ))
!    &           - RTP(IEL,ISCA(IFVAP))
!    &           * ( H02FOL + HRFVAP + CP2FOL * (TEBMOY-TREFTH) ) )
!            W8(IEL) = (RTP(IEL,ISCA(IHM))-EH2) * UNS1PW
! PPl 200106 c'était bien beau tant qu'il n'y avait que de la vapeur
!  mais avec la combustion hétérogène c'est le foutoir
!  on repasse (momentanément ?) à une enthalpie de phase

  x2h2 = 0.d0
  do icla = 1, nclafu
    x2h2 = x2h2 + rtp(iel,isca(ihlf(icla)))
  enddo
  w8(iel) = ( rtp(iel,isca(ihm)) - x2h2 )                         &
           *uns1pw

enddo


! --- Gestion memoire
!     Autres tableaux

! ------ Macro tableau d'entiers TBFUI : NTBFUI
!        Macro tableau de reels  TBFUR : NTBFUR
!        Macro tableau d'entiers TBWOI : NTBWOI
!        Macro tableau de reels  TBWOR : NTBWOR

ntbfui = 1
if ( ieqnox .eq. 0 ) then
  ntbfur = 11
else
  ntbfur = 12
endif
ntbwoi = 0
ntbwor = 2

call memfu1                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   ,                                                       &
   ntbfui , ifuwi  ,                                              &
   ntbfur , ifuwr  ,                                              &
   ntbwoi , iwori  ,                                              &
   ntbwor , iworr  ,                                              &
   ifinia , ifinra )

iphas  = 1

call fuphy1                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ncelet , ncel   ,                                              &
   ntbfui , ntbfur , ntbwoi , ntbwor ,                            &
   w2     , w4     , w6     ,                                     &
!         FVAP    FHTF    F4P2M
   w8     ,                                                       &
!         ENTH du gaz
   rtp    , propce  , propce(1,iromf) ,                           &
!                          ---------------- (masse vol. gaz)
   ia(ifuwi) , ra(ifuwr) ,                                        &
   ia(iwori) , ra(iworr)  )

!===============================================================================
! 4. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE DISPERSEE
!                    VALEURS CELLULES
!                    ----------------
!    TEMPERATURE
!===============================================================================

if ( ippmod(icfuel).ge.0 ) then

! --- Transport d'H2

  call futeh2                                                     &
  !==========
 ( ncelet , ncel   , nrtuse ,                                     &
   rtp    , propce , rtuser)

endif


!===============================================================================
! 5. CALCUL DES PROPRIETES PHYSIQUES DU MELANGE
!                    VALEURS CELLULES
!                    ----------------
!    MASSE VOLUMIQUE
!===============================================================================

! --- W2 = - Somme des X2(i)

do iel = 1, ncel
  w2(iel) = zero
enddo

do icla = 1, nclafu
  do iel = 1, ncel
    w2(iel) =  w2(iel)-rtp(iel,isca(iyfol(icla)))
  enddo
enddo

! --- Calcul de Rho du melange : 1/Rho = X1/Rho1 + Somme(X2/Rho2)
!     On sous relaxe quand on a un rho n a disposition, ie
!       a partir du deuxieme passage ou
!       a partir du premier passage si on est en suite de calcul et
!         qu'on a relu la masse volumique dans le fichier suite.

iphas = 1
ipcrom = ipproc(irom(iphas))

if (ipass.gt.1.or.(isuite.eq.1.and.initro(iphas).eq.1)) then
  srrom1 = srrom
else
  srrom1 = 0.d0
endif

do iel = 1, ncel

  x1sro1 = (1.d0+w2(iel)) / propce(iel,iromf)
  x2sro2 = zero
  do icla = 1, nclafu

    ipcro2 = ipproc(irom3(icla))
    propce(iel,ipcro2) = rho0fl

    x2sro2 = x2sro2 + rtp(iel,isca(iyfol(icla)))                  &
                     /propce(iel,ipcro2)

  enddo

! ---- Sous relaxation eventuelle a donner dans ppini1.F

  propce(iel,ipcrom) = srrom1*propce(iel,ipcrom)                  &
                     + (1.d0-srrom1)/(x1sro1+x2sro2)
 enddo

!===============================================================================
! 6. CALCUL DE RHO DU MELANGE

!                      VALEURS FACETTES
!                      ----------------
!===============================================================================

iphas = 1
ibrom(iphas) = 1
ipbrom = ipprob(irom(iphas))
ipcrom = ipproc(irom(iphas))

! ---> Masse volumique au bord pour toutes les facettes
!      Les facettes d'entree seront recalculees.

do ifac = 1, nfabor
  iel = ifabor(ifac)
  propfb(ifac,ipbrom) = propce(iel,ipcrom)
enddo

! ---> Masse volumique au bord pour les facettes d'entree UNIQUEMENT
!     Le test sur IZONE sert pour les reprises de calcul

if ( ipass.gt.1 .or. isuite.eq.1 ) then
  do ifac = 1, nfabor

    izone = izfppp(ifac)
    if(izone.gt.0) then
      if ( ientat(izone).eq.1 .or. ientfl(izone).eq.1 ) then
        qtotz  = qimpfl(izone) + qimpat(izone)
        x2tot  = qimpfl(izone) / qtotz
        x2sro2 = x2tot / rho0fl
        wmolme = (1.d0 + xsi) / (wmole(io2) + xsi * wmole(in2) )
        unsro1 = (wmolme * rr * timpat(izone)) / p0(iphas)
        x1sro1 = (1.d0 - x2tot) * unsro1
        propfb(ifac,ipbrom) = 1.d0 / (x1sro1 + x2sro2)
      endif
    endif

  enddo
endif

!----
! FIN
!----

return
end subroutine
