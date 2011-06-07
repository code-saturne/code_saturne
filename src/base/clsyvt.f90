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

subroutine clsyvt &
!================

 ( nvar   , nscal  ,                                              &
   icodcl ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefu  , rijipb , coefa  , coefb  )

!===============================================================================
! FONCTION :
! --------

! CONDITIONS LIMITES EN SYMETRIE POUR LES VECTEURS ET TENSEURS

! ON SUPPOSE QUE ICODCL(IU) = 4 =>
!                     SYMETRIE POUR LA VITESSE ET RIJ
!  (A PRIORI PEU RESTRICTIF EN MONOPHASIQUE)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
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
! coefu            ! tr ! <-- ! tab de trav pour valeurs en iprime             !
! (nfabor,3   )    !    !     !  des comp de la vitesse au bord                !
! rijipb           ! tr ! <-- ! tab de trav pour valeurs en iprime             !
! (nfabor,6   )    !    !     !  des rij au bord                               !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
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
use pointe
use entsor
use albase
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision coefu(nfabor,ndim), rijipb(nfabor,6)
double precision coefa(nfabor,*), coefb(nfabor,*)

! Local variables

integer          ifac, ii, isou
integer          iclu  , iclv  , iclw
integer          icl11 , icl22 , icl33 , icl12 , icl13 , icl23
integer          icluf , iclvf , iclwf
integer          iclvar
double precision rnx, rny, rnz, rxnn
double precision upx, upy, upz, usn
double precision tx, ty, tz, txn, t2x, t2y, t2z
double precision clsyme
double precision eloglo(3,3), alpha(6,6)
double precision srfbnf, rcodcn

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

icl11 = 0
icl22 = 0
icl33 = 0
icl12 = 0
icl13 = 0
icl23 = 0
iclvar = 0

! --- Memoire

! --- Conditions aux limites
iclu   = iclrtp(iu ,icoef)
iclv   = iclrtp(iv ,icoef)
iclw   = iclrtp(iw ,icoef)
if(itytur.eq.3) then
  icl11  = iclrtp(ir11,icoef)
  icl22  = iclrtp(ir22,icoef)
  icl33  = iclrtp(ir33,icoef)
  icl12  = iclrtp(ir12,icoef)
  icl13  = iclrtp(ir13,icoef)
  icl23  = iclrtp(ir23,icoef)
endif

icluf  = iclrtp(iu ,icoeff)
iclvf  = iclrtp(iv ,icoeff)
iclwf  = iclrtp(iw ,icoeff)


! --- Boucle sur les faces de bord : debut
do ifac = 1, nfabor

! --- Test sur la presence d'une condition de symetrie vitesse : debut
  if( icodcl(ifac,iu).eq.4 ) then

! --- Pour annuler le flux de masse
    isympa(ifac) = 0

! --- Grandeurs geometriques
    srfbnf = surfbn(ifac)

!===============================================================================
! 1. REPERE LOCAL
!      POUR LA VITESSE, SEULE EST NECESSAIRE LA NORMALE
!      POUR RIJ, IL FAUT LE REPERE COMPLET
!===============================================================================

! ---> NORMALE UNITAIRE

    rnx = surfbo(1,ifac)/srfbnf
    rny = surfbo(2,ifac)/srfbnf
    rnz = surfbo(3,ifac)/srfbnf

!     En ALE, on a eventuellement une vitesse de deplacement de la face
!       donc seule la composante normale importe (on continue a determiner
!       TX a partir de la vitesse tangentielle absolue car l'orientation
!       de TX et T2X est sans importance pour les symetries)
    rcodcn = 0.d0
    if (iale.eq.1) then
      rcodcn = rcodcl(ifac,iu,1)*rnx                           &
             + rcodcl(ifac,iv,1)*rny                           &
             + rcodcl(ifac,iw,1)*rnz
    endif

    upx = coefu(ifac,1)
    upy = coefu(ifac,2)
    upz = coefu(ifac,3)

    if (itytur.eq.3) then

! ---> VITESSE TANGENTIELLE RELATIVE

      usn = upx*rnx+upy*rny+upz*rnz
      tx  = upx -usn*rnx
      ty  = upy -usn*rny
      tz  = upz -usn*rnz
      txn = sqrt( tx**2 +ty**2 +tz**2 )

! ---> TANGENTE UNITAIRE

      if( txn.ge.epzero) then

        tx  = tx/txn
        ty  = ty/txn
        tz  = tz/txn

      else

!      SI LA VITESSE EST NULLE, LE VECTEUR T EST NORMAL ET QCQUE

        if(abs(rny).ge.epzero.or.abs(rnz).ge.epzero)then
          rxnn = sqrt(rny**2+rnz**2)
          tx  =  0.d0
          ty  =  rnz/rxnn
          tz  = -rny/rxnn
        elseif(abs(rnx).ge.epzero.or.abs(rnz).ge.epzero)then
          rxnn = sqrt(rnx**2+rnz**2)
          tx  =  rnz/rxnn
          ty  =  0.d0
          tz  = -rnx/rxnn
        else
          write(nfecra,1000)ifac,rnx,rny,rnz
          call csexit (1)
        endif

      endif


! ---> T2 = RN X T (OU X EST LE PRODUIT VECTORIEL)

      t2x = rny*tz - rnz*ty
      t2y = rnz*tx - rnx*tz
      t2z = rnx*ty - rny*tx

!     --> MATRICE ORTHOGONALE DE CHANGEMENT DE BASE ELOGLOij
!         (DE LA BASE LOCALE VERS LA BASE GLOBALE)

!                            |TX  -RNX  T2X|
!                   ELOGLO = |TY  -RNY  T2Y|
!                            |TZ  -RNZ  T2Z|

!         SA TRANSPOSEE ELOGLOt EST SON INVERSE


      eloglo(1,1) =  tx
      eloglo(1,2) = -rnx
      eloglo(1,3) =  t2x
      eloglo(2,1) =  ty
      eloglo(2,2) = -rny
      eloglo(2,3) =  t2y
      eloglo(3,1) =  tz
      eloglo(3,2) = -rnz
      eloglo(3,3) =  t2z

!     --> ON CALCULE ALPHA(6,6)

!       SOIT f LE CENTRE DE LA FACE DE BORD ET
!            I LE CENTRE DE LA CELLULE CORRESPONDANTE

!       EN NOTE RG (RESP RL) INDICE PAR f OU PAR I
!          LE TENSEUR DE REYNOLDS DANS LA BASE GLOBALE (RESP LOCALE)

!       LA MATRICE ALPHA APPLIQUEE AU VECTEUR GLOBAL EN I'
!         (RG11,I'|RG22,I'|RG33,I'|RG12,I'|RG13,I'|RG23,I')t
!         DOIT DONNER LES VALEURS A IMPOSER A LA FACE
!         (RG11,f |RG22,f |RG33,f |RG12,f |RG13,f |RG23,f )t
!         AUX CONDITIONS LIMITES DE DIRICHLET PRES (AJOUTEES ENSUITE)

!       ON LA DEFINIT EN CALCULANT RG,f EN FONCTION DE RG,I' COMME SUIT

!         RG,f = ELOGLO.RL,f.ELOGLOt (PRODUITS MATRICIELS)

!                          | RL,I'(1,1)     B*U*.Uk     C*RL,I'(1,3) |
!           AVEC    RL,f = | B*U*.Uk       RL,I'(2,2)       0        |
!                          | C*RL,I'(1,3)     0         RL,I'(3,3)   |

!                  AVEC    RL,I = ELOGLOt.RG,I'.ELOGLO
!                          B = 0
!                    ET    C = 0 EN PAROI (1 EN SYMETRIE)



!          ON CALCULE EN FAIT   ELOGLO.PROJECTEUR.ELOGLOt


      clsyme=1.d0
      call clca66 ( clsyme , eloglo , alpha )
      !==========

    endif

!===============================================================================
! 2. CONDITIONS SUR LES VITESSES (PARTIELLEMENT IMPLICITES)
!===============================================================================

    coefa(ifac,iclu) = rcodcn*rnx - rnx*(rny*upy+rnz*upz)
    coefb(ifac,iclu) = 1.d0-rnx**2
    coefa(ifac,iclv) = rcodcn*rny - rny*(rnz*upz+rnx*upx)
    coefb(ifac,iclv) = 1.d0-rny**2
    coefa(ifac,iclw) = rcodcn*rnz - rnz*(rnx*upx+rny*upy)
    coefb(ifac,iclw) = 1.d0-rnz**2

!===============================================================================
! 3. CONDITIONS SUR RIJ (PARTIELLEMENT IMPLICITES)
!===============================================================================

    if (itytur.eq.3) then

      do isou = 1, 6

        if(isou.eq.1) iclvar = icl11
        if(isou.eq.2) iclvar = icl22
        if(isou.eq.3) iclvar = icl33
        if(isou.eq.4) iclvar = icl12
        if(isou.eq.5) iclvar = icl13
        if(isou.eq.6) iclvar = icl23

        coefa(ifac,iclvar) = 0.0d0
        coefb(ifac,iclvar) = 0.0d0

      enddo

      do isou = 1,6

        if(isou.eq.1) then
          iclvar = icl11
        else if(isou.eq.2) then
          iclvar = icl22
        else if(isou.eq.3) then
          iclvar = icl33
        else if(isou.eq.4) then
          iclvar = icl12
        else if(isou.eq.5) then
          iclvar = icl13
        else if(isou.eq.6) then
          iclvar = icl23
        endif

!     IMPLICITATION PARTIELLE EVENTUELLE DES CL
        if (iclsyr.eq.1) then
          do ii = 1, 6
            if (ii.ne.isou) then
              coefa(ifac,iclvar) = coefa(ifac,iclvar) +           &
                   alpha(isou,ii) * rijipb(ifac,ii)
            endif
          enddo
          coefb(ifac,iclvar) = alpha(isou,isou)
        else
          do ii = 1, 6
            coefa(ifac,iclvar) = coefa(ifac,iclvar) +             &
                 alpha(isou,ii) * rijipb(ifac,ii)
          enddo
          coefb(ifac,iclvar) = 0.d0
        endif

      enddo

    endif

  endif
! --- Test sur la presence d'une condition de symetrie vitesse : fin

enddo
! ---  Boucle sur les faces de bord : fin

!===============================================================================
! 4. COEFAF et COEFBF BIDONS POUR LES VITESSES
!===============================================================================

  if(iclu.ne.icluf) then
    do ifac = 1, nfabor
      if( icodcl(ifac,iu).eq.4) then
        coefa(ifac,icluf) = coefa(ifac,iclu)
        coefb(ifac,icluf) = coefb(ifac,iclu)
        coefa(ifac,iclvf) = coefa(ifac,iclv)
        coefb(ifac,iclvf) = coefb(ifac,iclv)
        coefa(ifac,iclwf) = coefa(ifac,iclw)
        coefb(ifac,iclwf) = coefb(ifac,iclw)
      endif
    enddo
  endif

!===============================================================================
! 7.  FORMATS
!===============================================================================


#if defined(_CS_LANG_FR)

 1000 format(/,' LA NORMALE A LA FACE DE BORD DE SYMETRIE ',I10,/,&
         ' EST NULLE ; COORDONNEES : ',3E12.5)

#else

 1000 format(/,' THE NORMAL TO THE SYMMETRY BOUNDARY FACE ',I10,/,&
         ' IS NULL; COORDINATES: ',3E12.5)

#endif

!----
! FIN
!----

return
end subroutine
