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

subroutine altycl &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   itypfb , ialtyb , icodcl , impale ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   rcodcl , xyzno0 , depale ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!===============================================================================
! FONCTION :
! --------

! TRAITEMENT DES CODES DE CONDITIONS AUX LIMITES IALTYB POUR L'ALE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! itypfb           ! ia ! <-- ! boundary face types                            !
!  (nfabor, nphas) !    !     !                                                !
! ialtyb(nfabor    ! te ! <-- ! type des faces de bord pour l'ale              !
! icodcl           ! te ! <-- ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
! impale(nnod)     ! te ! <-- ! indicateur de delacement impose                !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! rcodcl           ! tr ! <-- ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2                   !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! depale(nnod,3    ! tr ! <-- ! deplacement aux noeuds                         !
! xyzno0(3,nnod    ! tr ! <-- ! coordonnees noeuds maillage initial            !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
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
use optcal
use cstnum
use cstphy
use entsor
use parall
use pointe
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas

integer          itypfb(nfabor,nphas)
integer          ialtyb(nfabor), icodcl(nfabor,nvar)
integer          impale(nnod)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision rcodcl(nfabor,nvar,3)
double precision depale(nnod,3), xyzno0(3,nnod)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, iel, iphas, iuiph, iviph, iwiph
integer          ii, inod, iecrw, icpt, ierror
double precision ddepx, ddepy, ddepz
double precision srfbnf, rnx, rny, rnz
double precision rcodcx, rcodcy, rcodcz, rcodsn


!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

! Mise a zero des RCODCL non specifies
do ifac = 1, nfabor
  if (rcodcl(ifac,iuma,1).gt.rinfin*0.5d0)                        &
       rcodcl(ifac,iuma,1) = 0.d0
  if (rcodcl(ifac,ivma,1).gt.rinfin*0.5d0)                        &
       rcodcl(ifac,ivma,1) = 0.d0
  if (rcodcl(ifac,iwma,1).gt.rinfin*0.5d0)                        &
       rcodcl(ifac,iwma,1) = 0.d0
enddo

!===============================================================================
! 2.  VERIFICATION DE LA CONSISTANCE DES IALTYB DONNES DANS USALCL
!===============================================================================
!  (valeur 0 autorisee)

ierror = 0
do ifac = 1, nfabor
  if (ialtyb(ifac).ne.0      .and.                                &
      ialtyb(ifac).ne.ibfixe .and.                                &
      ialtyb(ifac).ne.igliss .and.                                &
      ialtyb(ifac).ne.ivimpo ) then
    write(nfecra,1000)ifac,iprfml(ifmfbr(ifac),1),                &
         ialtyb(ifac)
    ierror = ierror + 1
  endif
enddo


!===============================================================================
! 3.  CONVERSION EN RCODCL ICODCL
!===============================================================================

! Calcul ou ecrasement de RCODCL a partir du deplacement impose,
!   si tous les noeuds d'une face sont a deplacement impose

do ifac = 1, nfabor
  iecrw = 0
  ddepx = 0.d0
  ddepy = 0.d0
  ddepz = 0.d0
  icpt  = 0
  do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
    inod = nodfbr(ii)
    if (impale(inod).eq.0) iecrw = iecrw + 1
    icpt = icpt + 1
    ddepx = ddepx + depale(inod,1)+xyzno0(1,inod)-xyznod(1,inod)
    ddepy = ddepy + depale(inod,2)+xyzno0(2,inod)-xyznod(2,inod)
    ddepz = ddepz + depale(inod,3)+xyzno0(3,inod)-xyznod(3,inod)
  enddo
  if (iecrw.eq.0) then
    iel = ifabor(ifac)
    ialtyb(ifac) = ivimpo
    rcodcl(ifac,iuma,1) = ddepx/dt(iel)/icpt
    rcodcl(ifac,ivma,1) = ddepy/dt(iel)/icpt
    rcodcl(ifac,iwma,1) = ddepz/dt(iel)/icpt
  endif
enddo

! Remplissage des autres RCODCL a partir des ITYALB

do ifac = 1, nfabor

  iel = ifabor(ifac)

! --> Faces fixes
!     On force alors les noeuds en question a etre fixes, pour eviter
!       des problemes eventuels aux coins

  if (ialtyb(ifac).eq.ibfixe) then
    icpt = 0
    if (icodcl(ifac,iuma).eq.0) then
      icpt = icpt + 1
      icodcl(ifac,iuma) = 1
      rcodcl(ifac,iuma,1) = 0.d0
    endif
    if (icodcl(ifac,ivma).eq.0) then
      icpt = icpt + 1
      icodcl(ifac,ivma) = 1
      rcodcl(ifac,ivma,1) = 0.d0
    endif
    if (icodcl(ifac,iwma).eq.0) then
      icpt = icpt + 1
      icodcl(ifac,iwma) = 1
      rcodcl(ifac,iwma,1) = 0.d0
    endif
!      Si on a fixe les trois composantes, alors on fixe les noeuds
!        correspondants. Sinon c'est que l'utilisateur a modifie quelque
!        chose     ... on le laisse seul maitre
    if (icpt.eq.3) then
      do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
        inod = nodfbr(ii)
        if (impale(inod).eq.0) then
          depale(inod,1) = 0.d0
          depale(inod,2) = 0.d0
          depale(inod,3) = 0.d0
          impale(inod) = 1
        endif
      enddo
    endif


! --> Faces de glissement
  elseif (ialtyb(ifac).eq.igliss) then

    if (icodcl(ifac,iuma).eq.0) icodcl(ifac,iuma) = 4
    if (icodcl(ifac,ivma).eq.0) icodcl(ifac,ivma) = 4
    if (icodcl(ifac,iwma).eq.0) icodcl(ifac,iwma) = 4

! --> Faces a vitesse imposee
  elseif (ialtyb(ifac).eq.ivimpo) then

    if (icodcl(ifac,iuma).eq.0) icodcl(ifac,iuma) = 1
    if (icodcl(ifac,ivma).eq.0) icodcl(ifac,ivma) = 1
    if (icodcl(ifac,iwma).eq.0) icodcl(ifac,iwma) = 1

  endif

enddo

!===============================================================================
! 4.  VERIFICATION DE LA CONSISTANCE DES ICODCL
!===============================================================================

!     IERROR a ete initialise plus haut
do ifac = 1, nfabor

  if (icodcl(ifac,iuma).ne.1 .and.                                &
      icodcl(ifac,iuma).ne.3 .and.                                &
      icodcl(ifac,iuma).ne.4 ) then
    write(nfecra,2000) ifac,iprfml(ifmfbr(ifac),1),               &
         icodcl(ifac,iuma)
    ierror = ierror + 1
  endif
  if (icodcl(ifac,ivma).ne.1 .and.                                &
      icodcl(ifac,ivma).ne.3 .and.                                &
      icodcl(ifac,ivma).ne.4 ) then
    write(nfecra,2001) ifac,iprfml(ifmfbr(ifac),1),               &
         icodcl(ifac,ivma)
    ierror = ierror + 1
  endif
  if (icodcl(ifac,iwma).ne.1 .and.                                &
      icodcl(ifac,iwma).ne.3 .and.                                &
      icodcl(ifac,iwma).ne.4 ) then
    write(nfecra,2002) ifac,iprfml(ifmfbr(ifac),1),               &
         icodcl(ifac,iwma)
    ierror = ierror + 1
  endif

  if ( ( icodcl(ifac,iuma).eq.4 .or.                              &
         icodcl(ifac,ivma).eq.4 .or.                              &
         icodcl(ifac,iwma).eq.4    ) .and.                        &
       ( icodcl(ifac,iuma).ne.4 .or.                              &
         icodcl(ifac,ivma).ne.4 .or.                              &
         icodcl(ifac,iwma).ne.4    ) ) then
    write(nfecra,3000) ifac,iprfml(ifmfbr(ifac),1),               &
         icodcl(ifac,iuma),icodcl(ifac,ivma),icodcl(ifac,iwma)
    ierror = ierror + 1
  endif

enddo

if (ierror.gt.0) then
  write(nfecra,4000)
  call csexit(1)
  !==========
endif

!===============================================================================
! 5.  VITESSE DE DEFILEMENT POUR LES PAROIS FLUIDES ET SYMETRIES
!===============================================================================

! Pour les symetries on rajoute toujours la vitesse de maillage, car on
!   ne conserve que la vitesse normale
! Pour les parois, on prend la vitesse de maillage si l'utilisateur n'a
!   pas specifie RCODCL, sinon on laisse RCODCL pour la vitesse tangente
!   et on prend la vitesse de maillage pour la composante normale.
! On se base uniquement sur ITYPFB, a l'utilisateur de gere les choses
!   s'il rentre en CL non standards.
iphas = 1
iuiph = iu
iviph = iv
iwiph = iw

do ifac = 1, nfabor

  if (ialtyb(ifac).eq.ivimpo) then

    if ( itypfb(ifac,iphas).eq.isymet ) then
      rcodcl(ifac,iuiph,1) = rcodcl(ifac,iuma,1)
      rcodcl(ifac,iviph,1) = rcodcl(ifac,ivma,1)
      rcodcl(ifac,iwiph,1) = rcodcl(ifac,iwma,1)
    endif

    if ( itypfb(ifac,iphas).eq.iparoi .or.                        &
         itypfb(ifac,iphas).eq.iparug ) then
! Si une des composantes de vitesse de glissement a ete
!    modifiee par l'utilisateur, on ne fixe que la vitesse
!    normale
      if (rcodcl(ifac,iuiph,1).gt.rinfin*0.5d0 .and.              &
           rcodcl(ifac,iviph,1).gt.rinfin*0.5d0 .and.             &
           rcodcl(ifac,iwiph,1).gt.rinfin*0.5d0) then
        rcodcl(ifac,iuiph,1) = rcodcl(ifac,iuma,1)
        rcodcl(ifac,iviph,1) = rcodcl(ifac,ivma,1)
        rcodcl(ifac,iwiph,1) = rcodcl(ifac,iwma,1)
      else
! On met a 0 les composantes de RCODCL non specifiees
        if (rcodcl(ifac,iuiph,1).gt.rinfin*0.5d0)                 &
             rcodcl(ifac,iuiph,1) = 0.d0
        if (rcodcl(ifac,iviph,1).gt.rinfin*0.5d0)                 &
             rcodcl(ifac,iviph,1) = 0.d0
        if (rcodcl(ifac,iwiph,1).gt.rinfin*0.5d0)                 &
             rcodcl(ifac,iwiph,1) = 0.d0

        srfbnf = surfbn(ifac)
        rnx = surfbo(1,ifac)/srfbnf
        rny = surfbo(2,ifac)/srfbnf
        rnz = surfbo(3,ifac)/srfbnf
        rcodcx = rcodcl(ifac,iuiph,1)
        rcodcy = rcodcl(ifac,iviph,1)
        rcodcz = rcodcl(ifac,iwiph,1)
        rcodsn = (rcodcl(ifac,iuma,1)-rcodcx)*rnx                 &
             +   (rcodcl(ifac,ivma,1)-rcodcy)*rny                 &
             +   (rcodcl(ifac,ivma,1)-rcodcz)*rnz
        rcodcl(ifac,iuiph,1) = rcodcx + rcodsn*rnx
        rcodcl(ifac,iviph,1) = rcodcy + rcodsn*rny
        rcodcl(ifac,iwiph,1) = rcodcz + rcodsn*rnz
      endif
    endif

  endif

enddo

!===============================================================================
! FORMATS
!===============================================================================

 1000 format(                                                           &
'@                                                            ',/,&
'@ METHODE ALE                                                ',/,&
'@ TYPE DE BORD NON RECONNU                                   ',/,&
'@   FACE ',I10   ,' ; PROPRIETE 1 :',I10                      ,/,&
'@     ALTYB FACE : ', I10                                     ,/,&
'@                                                            '  )
 2000 format(                                                           &
'@                                                            ',/,&
'@ METHODE ALE                                                ',/,&
'@ TYPE DE CONDITION A LA LIMITE NON RECONNU POUR LA VITESSE  ',/,&
'@   DE MAILLAGE SELON X (IUMA)                               ',/,&
'@   FACE ',I10   ,' ; PROPRIETE 1 :',I10                      ,/,&
'@     ICODCL : ', I10                                         ,/,&
'@                                                            ',/,&
'@ Les seules valeurs autorisees pour ICODCL sont             ',/,&
'@   1 : Dirichlet                                            ',/,&
'@   3 : Neumann                                              ',/,&
'@   4 : Glissement                                           ',/,&
'@                                                            '  )
 2001 format(                                                           &
'@                                                            ',/,&
'@ METHODE ALE                                                ',/,&
'@ TYPE DE CONDITION A LA LIMITE NON RECONNU POUR LA VITESSE  ',/,&
'@   DE MAILLAGE SELON Y (IVMA)                               ',/,&
'@   FACE ',I10   ,' ; PROPRIETE 1 :',I10                      ,/,&
'@     ICODCL : ', I10                                         ,/,&
'@                                                            ',/,&
'@ Les seules valeurs autorisees pour ICODCL sont             ',/,&
'@   1 : Dirichlet                                            ',/,&
'@   3 : Neumann                                              ',/,&
'@   4 : Glissement                                           ',/,&
'@                                                            '  )
 2002 format(                                                           &
'@                                                            ',/,&
'@ METHODE ALE                                                ',/,&
'@ TYPE DE CONDITION A LA LIMITE NON RECONNU POUR LA VITESSE  ',/,&
'@   DE MAILLAGE SELON Z (IWMA)                               ',/,&
'@   FACE ',I10   ,' ; PROPRIETE 1 :',I10                      ,/,&
'@     ICODCL : ', I10                                         ,/,&
'@                                                            ',/,&
'@ Les seules valeurs autorisees pour ICODCL sont             ',/,&
'@   1 : Dirichlet                                            ',/,&
'@   3 : Neumann                                              ',/,&
'@   4 : Glissement                                           ',/,&
'@                                                            '  )
 3000 format(                                                           &
'@                                                            ',/,&
'@ METHODE ALE                                                ',/,&
'@ INCOHERENCE DANS LES TYPES DE CONDITIONS A LA LIMITE       ',/,&
'@   POUR LA VITESSE DE MAILLAGE                              ',/,&
'@   FACE ',I10   ,' ; PROPRIETE 1 :',I10                      ,/,&
'@     ICODCL(.,IUMA) : ', I10                                 ,/,&
'@     ICODCL(.,IVMA) : ', I10                                 ,/,&
'@     ICODCL(.,IWMA) : ', I10                                 ,/,&
'@                                                            ',/,&
'@ Si une composante est traitee en glissement (ICODCL=4),    ',/,&
'@   toutes les composantes doivent etre traitees en          ',/,&
'@   glissement.                                              ',/,&
'@                                                            '  )
 4000 format(                                                           &
'@                                                            ',/,&
'@ METHODE ALE                                                ',/,&
'@ INCOHERENCE DANS LES CONDITIONS AUX LIMITES POUR LA VITESSE',/,&
'@   DE MAILLAGE                                              ',/,&
'@   (cf. message(s) ci-dessus)                               ',/,&
'@                                                            ',/,&
'@ Verifier usalcl.F                                          ',/,&
'@                                                            '  )


return
end subroutine
