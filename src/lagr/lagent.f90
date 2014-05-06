!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine lagent &
!================

 ( lndnod ,                                                       &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itycel , icocel , dlgeo  ,                                     &
   itypfb , itrifb , ifrlag , itepa  ,                            &
   dt     , rtpa   , propce ,                                     &
   ettp   , tepa   , vagaus )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!   Gestion de l'injection des particules dans le domaine de calcul

!     1. initialisation par l'utilisateur via USLAG2
!        des classes de particules et du type d'interaction
!        particule/face de frontiere.

!     2. injection des particules dans le domaine : initialisation
!        des tableau ETTP, ITEPA(IP,JISOR) et TEPA(IP,JRPOI).

!     3. modification des conditions d'injection des particules :
!        retouche des ETTP, ITEPA(IP,JISOR) et TEPA(IP,JRPOI).

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! lndnod           ! e  !  -> ! longueur du tableau icocel                     !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! dlgeo            ! tr ! --> ! tableau contenant les donnees geometriques     !
! (nfabor,ngeol)   !    !     ! pour le sous-modele de depot                   !
! icocel           ! te ! <-- ! connectivite cellules -> faces                 !
! (lndnod)         !    !     !    face de bord si numero negatif              !
! itycel           ! te ! <-- ! connectivite cellules -> faces                 !
! (ncelet+1)       !    !     !    pointeur du tableau icocel                  !
! itypfb(nfabor)   ! ia ! <-- ! boundary face types                            !
! itrifb(nfabor)   ! te ! <-- ! tab d'indirection pour tri des faces           !
! ifrlag           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module lagrangien                     !
! itepa            ! te ! --> ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant prec ou                   !
!                  !    !     !    instant courant si ntcabs = 1)              !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! vagaus           ! tr ! --> ! variables aleatoires gaussiennes               !
!(nbpmax,nvgaus    !    !     !                                                !
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
use entsor
use cstnum
use cstphy
use parall
use period
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use radiat
use ihmpre
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          lndnod
integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itypfb(nfabor) , itrifb(nfabor)
integer          icocel(lndnod) , itycel(ncelet+1)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)

double precision dt(ncelet) , rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision vagaus(nbpmax,*)
double precision dlgeo(nfabor,ngeol)

! Local variables

integer          iel , ifac , ip , nb , nc, ii, ilayer, ifvu
integer          iok , n1 , nd
integer          npt , npar1  , npar2 , mode , idvar
integer          ncmax, nzmax

double precision vn1 , vn2 , vn3 , pis6 , d3
double precision dmasse , rd(1) , aa
double precision xxpart , yypart , zzpart
double precision tvpart , uupart , vvpart , wwpart
double precision ddpart , ttpart
double precision surf   , volp , vitp
double precision dintrf(1)

integer, dimension(3) :: shpe
integer, allocatable, dimension(:) :: iwork
double precision, allocatable, dimension(:) :: surflag
double precision, allocatable, dimension(:,:) :: surlgrg
integer, allocatable, dimension(:) :: ninjrg
integer, allocatable, dimension(:,:,:) :: iusloc
integer, allocatable, dimension(:) :: ilftot

double precision, dimension(:,:), pointer :: vela

double precision unif(1), offset, rapsurf
integer irp, ipart, jj, kk, nfrtot, nlocnew, nbpartall
integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

call field_get_val_prev_v(ivarfl(iu), vela)

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

allocate(ilftot(nflagm))

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Init aberrante pour forcer l'utilisateur a mettre sa valeur

pis6 = pi / 6.d0

do nb = 1, nflagm
  iusncl(nb) = 0
  iusclb(nb) = 0
enddo

do ifac = 1,nfabor
  ifrlag(ifac) = 0
enddo

!    Mise a zero des debits pour chaque zone de bord

do nb = 1,nflagm
  deblag(nb) = 0.d0
enddo

!===============================================================================
! 2. Initialisation utilisateur par classe et par frontiere
!===============================================================================

if (iihmpr.eq.1) then
  call uilag2                                                                  &
  !==========
 ( nfabor, nozppm,                                                             &
   ientrl, isortl,  idepo1, idepo2,                                            &
   idepfa, iencrl,  irebol, isymtl, iphyla,                                    &
   ijnbp,  ijfre,   iclst,  ijuvw,  iuno,    iupt,    ivpt,  iwpt,             &
   ijprpd, ipoit,   idebt,  ijprdp, idpt,    ivdpt,                            &
   iropt,  ijprtp,  itpt,   icpt,   iepsi,                                     &
   nlayer, inuchl,  irawcl, ihpt,   ifrlag, iusncl,  iusclb  )
endif

call uslag2                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   dt     , rtpa   , propce ,                                     &
   ettp   , tepa   )


shpe = shape(iuslag)
ncmax = shpe(1)
nzmax = shpe(2)

allocate(iusloc(ncmax,nzmax,ndlaim))

do nd = 1,ndlaim
  do nb = 1, nzmax
    do nc = 1, ncmax
      iusloc(nc,nb,nd) = iuslag(nc,nb,nd)
    enddo
  enddo
enddo

!===============================================================================
! 3. Controles
!===============================================================================

iok = 0

! --> Les faces doivent toutes appartenir a une zone frontiere

do ifac = 1, nfabor
  if(ifrlag(ifac).le.0 .or. ifrlag(ifac).gt.nflagm) then
    iok = iok + 1
    write(nfecra,1000) ifac,nflagm,ifrlag(ifac)
  endif
enddo

if (iok.gt.0) then
  call csexit (1)
  !==========
endif

! --> On construit une liste des numeros des zones frontieres.

nfrlag = 0
do ifac = 1, nfabor
  ifvu = 0
  do ii = 1, nfrlag
    if (ilflag(ii).eq.ifrlag(ifac)) then
      ifvu = 1
    endif
  enddo
  if(ifvu.eq.0) then
    nfrlag = nfrlag + 1
    if(nfrlag.le.nflagm) then
      ilflag(nfrlag) = ifrlag(ifac)
    else
      write(nfecra,1001) nfrlag
      WRITE(NFECRA,'(I10)') (ILFLAG(II),II=1,NFRLAG)
      call csexit (1)
      !==========
    endif
  endif
enddo


! --> Calculation of the surfaces of the Lagrangian boundary zones

nrangp = irangp

if (irangp.ge.0) then
   nrangp = irangp
   call parmax(nrangp)
   allocate(surflag(nflagm))
   allocate(surlgrg(nflagm, nrangp + 1))
   allocate(ninjrg(nrangp  + 1))

   do kk = 1, nflagm
      surflag(kk) = 0.d0
      do jj = 1, nrangp + 1
         surlgrg(kk,jj) = 0.d0
      enddo
   enddo

   do ii = 1, nfrlag

      surflag(ilflag(ii)) = 0.d0
      do ifac = 1, nfabor
         if (ilflag(ii).eq.ifrlag(ifac)) then
            surflag(ilflag(ii)) = surflag(ilflag(ii)) + surfbn(ifac)
            surlgrg(ilflag(ii), irangp + 1) =                                  &
                 surlgrg(ilflag(ii), irangp + 1) + surfbn(ifac)
         endif
      enddo
   enddo

   do kk = 1, nflagm
      call parsom(surflag(kk))
      do jj = 1, nrangp + 1
         call parsom(surlgrg(kk, jj))
      enddo
   enddo

   if (irangp.eq.0) then
      nfrtot = 0
      jj = 1
      do kk = 1, nflagm
         if (surflag(kk).gt.1.d-15) then
            nfrtot = nfrtot + 1
            ilftot(jj) = kk
            jj = jj + 1
         endif
      enddo
   endif

   call parbci(0, 1, nfrtot)
   call parbci(0, nfrtot, ilftot)

else

   nfrtot = nfrlag
   do ii = 1, nfrlag
      ilftot(ii) = ilflag(ii)
   enddo

endif

! --> Nombre de classes.

do ii = 1,nfrlag
  nb = ilflag(ii)
  if (iusncl(nb).lt.0) then
    iok = iok + 1
    write(nfecra,1010) nb,iusncl(nb)
  endif
enddo

! --> Nombre de particules.

do ii = 1,nfrlag
  nb = ilflag(ii)
  do nc = 1,iusncl(nb)
    if ( iuslag(nc,nb,ijnbp).lt.0 .or.                            &
         iuslag(nc,nb,ijnbp).gt.nbpmax ) then
      iok = iok + 1
      write(nfecra,1020) nc,nb,nbpmax,iuslag(nc,nb,ijnbp)
    endif
  enddo
enddo

!    Verification des classes de particules : Uniquement un warning

if (nbclst.gt.0 ) then
  do ii = 1,nfrlag
    nb = ilflag(ii)
    do nc = 1,iusncl(nb)
      if (iuslag(nc,nb,iclst).le.0      .or.                      &
          iuslag(nc,nb,iclst).gt.nbclst        ) then
        write(nfecra,1021) nc,nb,nbclst,iuslag(nc,nb,iclst)
      endif
    enddo
  enddo
endif

! --> Frequence d'injection.

do ii = 1,nfrlag
  nb = ilflag(ii)
  do nc = 1, iusncl(nb)
    if (iuslag(nc,nb,ijfre).lt.0) then
      iok = iok + 1
      write(nfecra,1030) nb,nc,iuslag(nc,nb,ijfre)
    endif
  enddo
enddo

! --> Conditions au bord.

do ii = 1,nfrlag
  nb = ilflag(ii)
  if ( iusclb(nb).ne.ientrl .and. iusclb(nb).ne.isortl .and.      &
       iusclb(nb).ne.irebol .and. iusclb(nb).ne.idepo1 .and.      &
       iusclb(nb).ne.idepo2 .and. iusclb(nb).ne.isymtl .and.      &
       iusclb(nb).ne.iencrl .and. iusclb(nb).ne.jbord1 .and.      &
       iusclb(nb).ne.jbord2 .and. iusclb(nb).ne.jbord3 .and.      &
       iusclb(nb).ne.jbord4 .and. iusclb(nb).ne.jbord5 .and.      &
       iusclb(nb).ne.idepfa ) then
    iok = iok + 1
    write(nfecra,1040) nb
  endif
enddo

do ii = 1,nfrlag
  nb = ilflag(ii)
  if (iusclb(nb).eq.iencrl .and. iphyla.ne.2) then
    iok = iok + 1
    write(nfecra,1042) nb
  endif
enddo

do ii = 1,nfrlag
  nb = ilflag(ii)
  if (iusclb(nb).eq.iencrl .and. iencra.ne.1) then
    iok = iok + 1
    write(nfecra,1043) nb
  endif
enddo

! --> Type de condition pour le taux de presence.

do ii = 1, nfrlag
  nb = ilflag(ii)
  do nc = 1, iusncl(nb)
    if ( iuslag(nc,nb,ijprpd).lt.1 .or.                           &
         iuslag(nc,nb,ijprpd).gt.2      ) then
      iok = iok + 1
      write(nfecra,1053) nb, nc, iuslag(nc,nb,ijprpd)
    endif
  enddo
enddo

! --> Type de condition pour la vitesse.

do ii = 1, nfrlag
  nb = ilflag(ii)
  do nc = 1, iusncl(nb)
    if ( iuslag(nc,nb,ijuvw).lt.-1 .or.                           &
         iuslag(nc,nb,ijuvw).gt.2      ) then
      iok = iok + 1
      write(nfecra,1050) nb,nc,iuslag(nc,nb,ijuvw)
    endif
  enddo
enddo

! --> Type de condition pour le diametre.

do ii = 1, nfrlag
  nb = ilflag(ii)
  do nc = 1, iusncl(nb)
    if ( iuslag(nc,nb,ijprdp).lt.1 .or.                           &
         iuslag(nc,nb,ijprdp).gt.2      ) then
      iok = iok + 1
      write(nfecra,1051) nb, nc, iuslag(nc,nb,ijprdp)
    endif
  enddo
enddo

! --> Type de condition pour le diametre.

if ( iphyla.eq.1 .and.                                            &
       (itpvar.eq.1 .or. idpvar.eq.1 .or. impvar.eq.1) ) then
  do ii = 1, nfrlag
    nb = ilflag(ii)
    do nc = 1, iusncl(nb)
      if ( iuslag(nc,nb,ijprtp).lt.1 .or.                         &
           iuslag(nc,nb,ijprtp).gt.2      ) then
        iok = iok + 1
        write(nfecra,1052) nb, nc, iuslag(nc,nb,ijprtp)
      endif
    enddo
  enddo
endif

! --> Poids statistiques

do ii = 1, nfrlag
  nb = ilflag(ii)
  do nc = 1, iusncl(nb)
    if (ruslag(nc,nb,ipoit).le.0.d0) then
      iok = iok + 1
      write(nfecra,1055) nb, nc, ruslag(nc,nb,ipoit)
    endif
  enddo
enddo

! --> Debit massique de particule

do ii = 1, nfrlag
  nb = ilflag(ii)
  do nc = 1, iusncl(nb)
    if (ruslag(nc,nb,idebt).lt.0.d0) then
      iok = iok + 1
      write(nfecra,1056) nb, nc, ruslag(nc,nb,idebt)
    endif
  enddo
enddo
do ii = 1, nfrlag
  nb = ilflag(ii)
  do nc = 1, iusncl(nb)
    if (ruslag(nc,nb,idebt).gt.0.d0 .and.                         &
        iuslag(nc,nb,ijnbp).eq.0         ) then
      iok = iok + 1
      write(nfecra,1057) nb, nc, ruslag(nc,nb,idebt),             &
                                 iuslag(nc,nb,ijnbp)
    endif
  enddo
enddo

! --> Proprietes des particules : le diametre, son ecart-type, et rho

do ii = 1, nfrlag
  nb = ilflag(ii)
  do nc = 1, iusncl(nb)
    if (iphyla .ne. 2) then
      if (ruslag(nc,nb,iropt).lt.0.d0 .or.                          &
          ruslag(nc,nb,idpt) .lt.0.d0 .or.                          &
          ruslag(nc,nb,ivdpt).lt.0.d0       ) then
        iok = iok + 1
        write(nfecra,1060) nb, nc,                                  &
                           ruslag(nc,nb,iropt),                     &
                           ruslag(nc,nb,idpt),                      &
                           ruslag(nc,nb,ivdpt)
      endif
    endif
  enddo
enddo

do ii = 1, nfrlag
  nb = ilflag(ii)
  do nc = 1, iusncl(nb)
    if ( ruslag(nc,nb,idpt).lt.3.d0*ruslag(nc,nb,ivdpt) ) then
      iok = iok + 1
      write(nfecra,1065) nb, nc,                                  &
        ruslag(nc,nb,idpt)-3.d0*ruslag(nc,nb,ivdpt)
    endif
  enddo
enddo

! --> Proprietes des particules : Temperature et CP

if (iphyla.eq.1 .and. itpvar.eq.1) then

  do ii = 1, nfrlag
    nb = ilflag(ii)
    do nc = 1, iusncl(nb)
      if (ruslag(nc,nb,icpt)  .lt.0.d0   .or.                     &
          ruslag(nc,nb,itpt)  .lt.tkelvn       ) then
        iok = iok + 1
        write(nfecra,1070)                                        &
        iphyla, itpvar, nb, nc,                                   &
        ruslag(nc,nb,itpt), ruslag(nc,nb,icpt)
      endif
    enddo
  enddo

endif

! --> Proprietes des particules : Emissivite

if (iphyla.eq.1 .and. itpvar.eq.1 .and. iirayo.gt.0) then

  do ii = 1, nfrlag
    nb = ilflag(ii)
    do nc = 1, iusncl(nb)
      if (ruslag(nc,nb,iepsi) .lt.0.d0 .or.                       &
          ruslag(nc,nb,iepsi) .gt.1.d0      ) then
        iok = iok + 1
        write(nfecra,1075)                                        &
        iphyla, itpvar, nb, nc, ruslag(nc,nb,iepsi)
      endif
    enddo
  enddo

endif

! Charbon

if (iphyla.eq.2) then

! --> Numero du charbon

  do ii = 1, nfrlag
    nb = ilflag(ii)
    do nc = 1, iusncl(nb)
      if (iuslag(nc,nb,inuchl).lt.0.d0                            &
          .and.  iuslag(nc,nb,inuchl).gt.ncharb) then
        iok = iok + 1
        write(nfecra,1080) nb, nc, ncharb, iuslag(nc,nb,inuchl)
      endif
    enddo
  enddo

! --> Proprietes des particules de Charbon.

  do ii = 1, nfrlag
    nb = ilflag(ii)

    do nc = 1, iusncl(nb)

      if (iuslag(nc,nb,irawcl).lt.0                            &
          .or.  iuslag(nc,nb,irawcl).gt.1) then

        iok = iok + 1
        write(nfecra,1081) nb, nc, iuslag(nc,nb,irawcl)

      else if (iuslag(nc,nb,irawcl).eq.0                       &
          .and. iuslag(nc,nb,ijprdp).eq.2) then

        write(nfecra,1082) nb, nc, iuslag(nc,nb,irawcl) ,      &
                           iuslag(nc,nb,ijprdp)

      else if (iuslag(nc,nb,irawcl).eq.0                       &
          .and. iuslag(nc,nb,ijprdp).eq.1                      &
          .and. ruslag(nc,nb,ivdpt).gt.0.0d0) then

        iok = iok + 1
        write(nfecra,1083) nb, nc, iuslag(nc,nb,irawcl) ,      &
                           ruslag(nc,nb,ivdpt)
      endif

      do ilayer = 1, nlayer

        if (ruslag(nc,nb,ihpt(ilayer)) .lt. tkelvi) then
          iok = iok + 1
          write(nfecra,1084)                                     &
          iphyla, nb, nc, ilayer,                                &
          ruslag(nc,nb,ihpt(ilayer))
        endif

      enddo

    enddo
  enddo

! --> Proprietes des particules de Charbon.

  do ii = 1, nfrlag

    nb = ilflag(ii)

    do nc = 1, iusncl(nb)

      ! irawcl = 0 --> Composition du charbon definie par l'utilisateur dans uslag2
      ! on verifie les donnes contenues dans le tableau ruslag

      if (iuslag(nc,nb,irawcl) .eq. 0) then

        if (ruslag(nc,nb,iropt)   .lt. 0.d0 .or.                 &
            ruslag(nc,nb,icpt)    .lt. 0.d0 .or.                 &
            ruslag(nc,nb,ifrmwt)  .lt. 0.d0 .or.                 &
            ruslag(nc,nb,ifrmwt)  .gt. 1.d0  ) then
          iok = iok + 1
          write(nfecra,1090)                                     &
          iphyla, nb, nc, ruslag(nc,nb,iropt),                   &
                          ruslag(nc,nb,icpt),                    &
                          ruslag(nc,nb,ifrmwt)
        endif

        do ilayer = 1, nlayer

          if (ruslag(nc,nb,ifrmch(ilayer))  .lt. 0.0d0  .or.     &
              ruslag(nc,nb,ifrmch(ilayer))  .gt. 1.0d0  .or.     &
              ruslag(nc,nb,ifrmck(ilayer))  .lt. 0.0d0  .or.     &
              ruslag(nc,nb,ifrmck(ilayer))  .gt. 1.0d0  .or.     &
              ruslag(nc,nb,irhock0(ilayer)) .lt. 0.0d0  ) then

            iok = iok + 1
            write(nfecra,1091)                                   &
            iphyla, nb, nc, ilayer,                              &
            ruslag(nc,nb,ifrmch(ilayer)),                        &
            ruslag(nc,nb,ifrmck(ilayer)),                        &
            ruslag(nc,nb,irhock0(ilayer))

          endif
        enddo

        if (ruslag(nc,nb,irdck) .lt. 0.0d0    .or.               &
            ruslag(nc,nb,ird0p) .lt. 0.0d0    ) then

          iok = iok + 1
          write(nfecra,1092)                                     &
          iphyla, iuslag(nc,nb,irawcl), nb, nc,                  &
          ruslag(nc,nb,irdck),                                   &
          ruslag(nc,nb,ird0p)

        endif

      ! irawcl = 1 --> Composition du charbon definie dans le fichier XML (DP_FCP)
      ! on verifie les donnes contenues dans le XML

      else if (iuslag(nc,nb,irawcl) .eq. 1 ) then

        if (rho0ch(iuslag(nc,nb,inuchl)).lt. 0.d0 .or.           &
            cp2ch(iuslag(nc,nb,inuchl)) .lt. 0.d0 .or.           &
            xwatch(iuslag(nc,nb,inuchl)).lt. 0.d0 .or.           &
            xwatch(iuslag(nc,nb,inuchl)).gt. 1.d0 .or.           &
            xashch(iuslag(nc,nb,inuchl)).lt. 0.d0 .or.           &
            xashch(iuslag(nc,nb,inuchl)).gt. 1.d0 ) then

          iok = iok + 1
          write(nfecra,1093)                                     &
          iphyla, nb, nc, iuslag(nc,nb,inuchl),                  &
          rho0ch(iuslag(nc,nb,inuchl)),                          &
          cp2ch(iuslag(nc,nb,inuchl)),                           &
          xwatch(iuslag(nc,nb,inuchl)),                          &
          xashch(iuslag(nc,nb,inuchl))

        endif

        if (xwatch(iuslag(nc,nb,inuchl)) +                       &
            xashch(iuslag(nc,nb,inuchl)) .gt. 1.0) then

            iok = iok + 1
            write(nfecra,1094)                                   &
            iphyla, nb, nc, iuslag(nc,nb,inuchl),                &
            xwatch(iuslag(nc,nb,inuchl)),                        &
            xashch(iuslag(nc,nb,inuchl)),                        &
            xwatch(iuslag(nc,nb,inuchl))+xashch(iuslag(nc,nb,inuchl))

        endif

        if (ruslag(nc,nb,iropt)  .ge. 0.0d0  .or.                &
            ruslag(nc,nb,ifrmwt) .ge. 0.0d0  .or.                &
            ruslag(nc,nb,icpt)   .ge. 0.0d0  .or.                &
            ruslag(nc,nb,irdck)  .ge. 0.0d0  .or.                &
            ruslag(nc,nb,ird0p)  .ge. 0.0d0  ) then

          iok = iok + 1
          write(nfecra,1095)                                     &
          iphyla, nb, nc, iuslag(nc,nb,irawcl),-grand,           &
          ruslag(nc,nb,iropt),                                   &
          ruslag(nc,nb,ifrmwt),                                  &
          ruslag(nc,nb,icpt),                                    &
          ruslag(nc,nb,irdck),                                   &
          ruslag(nc,nb,ird0p)

        endif

        do ilayer = 1, nlayer

          if (ruslag(nc,nb,ifrmch(ilayer))  .ge. 0.0d0  .or.     &
              ruslag(nc,nb,ifrmck(ilayer))  .ge. 0.0d0  .or.     &
              ruslag(nc,nb,irhock0(ilayer)) .ge. 0.0d0  ) then

            iok = iok + 1
            write(nfecra,1096)                                   &
            iphyla, nb, nc, ilayer,                              &
            iuslag(nc,nb,irawcl),-grand,                         &
            ruslag(nc,nb,ifrmch(ilayer)),                        &
            ruslag(nc,nb,ifrmck(ilayer)),                        &
            ruslag(nc,nb,irhock0(ilayer))

          endif

        enddo

      endif
    enddo
  enddo
endif

! --> Stop si erreur.

if(iok.gt.0) then
  call csexit (1)
  !==========
endif


!===============================================================================
! 4. Transformation des donnees utilisateur
!===============================================================================

! --> Injection des part 1ere iter seulement si freq d'injection nulle

do ii = 1, nfrtot
  nb = ilftot(ii)
  do nc = 1, iusncl(nb)
    if (iuslag(nc,nb,ijfre).eq.0 .and. iplas.eq.1) then
      iuslag(nc,nb,ijfre) = ntcabs
      iusloc(nc,nb,ijfre) = ntcabs
    endif
    if (iuslag(nc,nb,ijfre).eq.0 .and. iplas.gt.1) then
      iuslag(nc,nb,ijfre) = ntcabs+1
      iusloc(nc,nb,ijfre) = ntcabs+1
    endif
  enddo
enddo

! --> Calcul du nombre de particules a injecter pour cette iteration

nbpnew = 0
dnbpnw = 0.d0

!     Dans le cas ou on a un taux de presence impose dans une zone,
!     on corrige IUSLAG(NC,NB,IJNBP) donne dans USLAG2 qui n'a
!     pas de sens puisque l'on injecte 1 particule par maille

do ii = 1,nfrlag
  nb = ilflag(ii)
!       pour chaque classe :
  do nc = 1, iusncl(nb)
!         si de nouvelles particules doivent entrer :
    if (mod(ntcabs,iuslag(nc,nb,ijfre)).eq.0) then

      if ( iuslag(nc,nb,ijprpd) .eq. 2 ) then
        iuslag(nc,nb,ijnbp)=0
        do ifac = 1,nfabor
          if (ifrlag(ifac).eq.nb) then
            iuslag(nc,nb,ijnbp)=iuslag(nc,nb,ijnbp)+1
          endif
        enddo
      endif
    endif
  enddo
enddo

do ii = 1,nfrtot
  nb = ilftot(ii)
  do nc = 1, iusncl(nb)
    if (mod(ntcabs,iuslag(nc,nb,ijfre)).eq.0) then
      nbpnew = nbpnew + iuslag(nc,nb,ijnbp)
    endif
  enddo
enddo

! --> Limite du nombre de particules a NBPMAX
nbpartall = nbpart
if (irangp.ge.0) then
   call parsom(nbpartall)
endif

if ( (nbpartall+nbpnew).gt.nbpmax ) then
  write(nfecra,3000) nbpartall,nbpnew,nbpmax
  nbpnew = 0
endif

! --> Si pas de new particules alors RETURN
if (nbpnew.eq.0) return

! --> Tirage aleatoire des positions des NBPNEW nouvelles particules
!   au niveau des zones de bord et reperage des cellules correspondantes

!   initialisation du compteur de nouvelles particules

npt = nbpart

! Allocate a work array
allocate(iwork(nbpmax))

! Initialisation du nombre local
! de particules injectées par rang
nlocnew = 0

!     Ensuite, on regarde ou on les met

!     pour chaque zone de bord :
do ii = 1,nfrtot
   nb = ilftot(ii)
   !       pour chaque classe :
   do nc = 1, iusncl(nb)
      !         si de nouvelles particules doivent entrer :
      if (mod(ntcabs,iuslag(nc,nb,ijfre)).eq.0) then

         if ( iuslag(nc,nb,ijprpd) .eq. 1 ) then


          ! Calcul sur le rang 0 du nombre de particules à injecter pour chaque rang
          ! base sur la surface relative de chaque zone d'injection presente sur
          ! chaque rang --> remplissage du tableau ninjrg(nrangp+1)

            if (irangp.eq.0) then

               do irp = 1, nrangp + 1
                  ninjrg(irp) = 0
               enddo

               do ipart = 1, iuslag(nc,nb,ijnbp)

                  call zufall(1, unif)

                  ! blindage
                  unif(1) = unif(1) + 1.d-9

                  irp = 1
                  offset = surlgrg(nb,irp) / surflag(nb)
                  do while (unif(1).gt.offset)
                     irp = irp + 1
                     offset = offset + surlgrg(nb,irp) / surflag(nb)
                  enddo
                  ninjrg(irp) = ninjrg(irp) + 1
               enddo
            endif

            ! Broadcast a tous les rangs
            if (irangp.ge.0) then
               call parbci(0, nrangp + 1, ninjrg)
            endif

            ! Fin du calcul du nombre de particules à injecter

            if (irangp.ge.0) then
               iusloc(nc,nb,ijnbp) = ninjrg(irangp+1)
               nlocnew = nlocnew + ninjrg(irangp+1)
            else
               iusloc(nc,nb,ijnbp) = iuslag(nc,nb,ijnbp)
               nlocnew = nlocnew + iuslag(nc,nb,ijnbp)
            endif

            if  (iusloc(nc,nb,ijnbp).gt.0) then

               call lagnew                                                   &
               !==========
             ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                  &
               npt    , nlocnew ,      iusloc(nc,nb,ijnbp)  ,                &
               nb     ,                                                      &
               ifrlag , itepa(:,jisor)  , iwork  ,                           &
               ettp   )

            endif

         elseif ( iuslag(nc,nb,ijprpd) .eq. 2 ) then

            call lagnpr                                                      &
            !==========
           ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                    &
             npt    , nlocnew ,  iusloc(nc,nb,ijnbp)   ,                     &
             nb     ,                                                        &
             ifrlag , itepa(:,jisor)  , iwork  ,                             &
             ettp   )
         endif

      endif
   enddo
enddo

!-->TEST DE CONTROLE (NE PAS MODIFIER)

if ( (nbpart+nlocnew).ne.npt ) then
  write(nfecra,3010) nlocnew, npt-nbpart
  call csexit (1)
  !==========
endif


!   reinitialisation du compteur de nouvelles particules

npt = nbpart

!     pour chaque zone de bord:
do ii = 1,nfrtot
  nb = ilftot(ii)

!       pour chaque classe :
  do nc = 1, iusncl(nb)

!         si de nouvelles particules doivent entrer :
     if (mod(ntcabs,iuslag(nc,nb,ijfre)).eq.0) then

      do ip = npt+1 , npt + iusloc(nc,nb,ijnbp)
        iel = itepa(ip,jisor)
        ifac = iwork(ip)

        itepa(ip,jord1) = 0

!-->COMPOSANTES DE LA VITESSE DES PARTICULES

!             si composantes de la vitesse imposee :
        if (iuslag(nc,nb,ijuvw).eq.1) then
          ettp(ip,jup) = ruslag(nc,nb,iupt)
          ettp(ip,jvp) = ruslag(nc,nb,ivpt)
          ettp(ip,jwp) = ruslag(nc,nb,iwpt)

!             si norme de la vitesse imposee :
        else if (iuslag(nc,nb,ijuvw).eq.0) then
          aa = -1.d0 / surfbn(ifac)
          vn1 = surfbo(1,ifac) * aa
          vn2 = surfbo(2,ifac) * aa
          vn3 = surfbo(3,ifac) * aa
          ettp(ip,jup) = vn1 * ruslag(nc,nb,iuno)
          ettp(ip,jvp) = vn2 * ruslag(nc,nb,iuno)
          ettp(ip,jwp) = vn3 * ruslag(nc,nb,iuno)

!             si vitesse du fluide vu :
        else if (iuslag(nc,nb,ijuvw).eq.-1) then
          ettp(ip,jup) = vela(1,iel)
          ettp(ip,jvp) = vela(2,iel)
          ettp(ip,jwp) = vela(3,iel)

!             si profil de vitesse impose :
        else if (iuslag(nc,nb,ijuvw).eq.2) then

         idvar = 1
         xxpart = ettp(ip,jxp)
         yypart = ettp(ip,jyp)
         zzpart = ettp(ip,jzp)

         call uslapr                                              &
         !==========
 ( idvar  , iel    , nb     , nc     ,                            &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   xxpart , yypart , zzpart ,                                     &
   tvpart , uupart , vvpart , wwpart , ddpart , ttpart ,          &
   dt     , rtpa   , propce ,                                     &
   ettp   , tepa   )

          ettp(ip,jup) = uupart
          ettp(ip,jvp) = vvpart
          ettp(ip,jwp) = wwpart

        endif

!-->Vitesse du fluide vu

        ettp(ip,juf) = vela(1,iel)
        ettp(ip,jvf) = vela(2,iel)
        ettp(ip,jwf) = vela(3,iel)

!--> TEMPS DE SEJOUR

        tepa(ip,jrtsp) = 0.d0

!--> Diametre

!             si diametre constant imposee :
        if (iuslag(nc,nb,ijprdp).eq.1) then
          if (ruslag(nc,nb,ivdpt) .gt. 0.d0) then
            n1 = 1
            call normalen(n1,rd)
            ettp(ip,jdp) = ruslag(nc,nb,idpt)                     &
                         + rd(1) * ruslag(nc,nb,ivdpt)

!    On verifie qu'on obtient un diametre dans la gamme des 99,7%

            d3 = 3.d0 * ruslag(nc,nb,ivdpt)
            if (ettp(ip,jdp).lt.ruslag(nc,nb,idpt)-d3)            &
              ettp(ip,jdp)= ruslag(nc,nb,idpt)
            if (ettp(ip,jdp).gt.ruslag(nc,nb,idpt)+d3)            &
              ettp(ip,jdp)= ruslag(nc,nb,idpt)
          else
            ettp(ip,jdp) = ruslag(nc,nb,idpt)
          endif

!             si profil pour le diametre  :
        else if (iuslag(nc,nb,ijprdp).eq.2) then

          idvar = 2
          xxpart = ettp(ip,jxp)
          yypart = ettp(ip,jyp)
          zzpart = ettp(ip,jzp)

          call uslapr                                             &
          !==========
 ( idvar  , iel    , nb     , nc     ,                            &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   xxpart , yypart , zzpart ,                                     &
   tvpart , uupart , vvpart , wwpart , ddpart , ttpart ,          &
   dt     , rtpa   , propce ,                                     &
   ettp   , tepa   )

          ettp(ip,jdp) = ddpart

        endif

!--> Autres variables : masse, ... en fonction de la physique

        d3 = ettp(ip,jdp) * ettp(ip,jdp) * ettp(ip,jdp)

        if (nbclst.gt.0) then
          itepa(ip,jclst) = iuslag(nc,nb,iclst)
        endif

        if ( iphyla.eq.0 .or. iphyla.eq.1 ) then

          ettp(ip,jmp) = ruslag(nc,nb,iropt) * pis6 * d3

          if ( iphyla.eq.1 .and. itpvar.eq.1 ) then

!             si Temperature constante imposee :
            if (iuslag(nc,nb,ijprtp).eq.1) then
              ettp(ip,jtp) = ruslag(nc,nb,itpt)
!             si profil pour la temperature :
            else if (iuslag(nc,nb,ijprtp).eq.2) then

              idvar = 3
              xxpart = ettp(ip,jxp)
              yypart = ettp(ip,jyp)
              zzpart = ettp(ip,jzp)

              call uslapr                                         &
              !==========
 ( idvar  , iel    , nb     , nc     ,                            &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   xxpart , yypart , zzpart ,                                     &
   tvpart , uupart , vvpart , wwpart , ddpart , ttpart ,          &
   dt     , rtpa   , propce ,                                     &
   ettp   , tepa   )

              ettp(ip,jtp) = ttpart

            endif

            if ( ippmod(iccoal).ge.0 .or.                         &
                 ippmod(icpl3c).ge.0 .or.                         &
                 ippmod(icfuel).ge.0      ) then

              ettp(ip,jtf) = propce(iel,ipproc(itemp1)) -tkelvi

            else if ( ippmod(icod3p).ge.0 .or.                    &
                      ippmod(icoebu).ge.0 .or.                    &
                      ippmod(ielarc).ge.0 .or.                    &
                      ippmod(ieljou).ge.0      ) then

              ettp(ip,jtf) = propce(iel,ipproc(itemp)) -tkelvi

            else if (itherm.eq.1) then

              if (itpscl.eq.1) then !Kelvin

                ettp(ip,jtf) = rtpa(iel,isca(iscalt)) -tkelvi

              else if (itpscl.eq.2) then ! Celsius

                ettp(ip,jtf) = rtpa(iel,isca(iscalt))

              endif

            else if (itherm.eq.2) then

              mode = 1
              call usthht(mode, rtpa(iel,isca(iscalt)), ettp(ip,jtf))

            endif

            ettp(ip,jcp) = ruslag(nc,nb,icpt)
            tepa(ip,jreps) = ruslag(nc,nb,iepsi)

          endif

        else if ( iphyla.eq.2 ) then

          ! Remplissage de ITEPA
          itepa(ip,jinch)  = iuslag(nc,nb,inuchl)

          ! Remplissage de ETTP
          ettp(ip,jtaux) = 0.0d0 ! non utilise pour iphyla=2
          ettp(ip,jtf) = propce(iel,ipproc(itemp1)) - tkelvi

          do ilayer = 1, nlayer

            ettp(ip,jhp(ilayer))  = ruslag(nc,nb,ihpt(ilayer))

          enddo

          ! user-defined composition (uslag2)
          if (iuslag(nc,nb,irawcl).eq.0) then

            ! Remplissage de ETTP
            ettp(ip,jcp) = ruslag(nc,nb,icpt)
            ettp(ip,jmp) = ruslag(nc,nb,iropt) * pis6 * d3
            ettp(ip,jmwat) = ruslag(nc,nb,ifrmwt) * ettp(ip,jmp)

            do ilayer = 1, nlayer

              ettp(ip,jmch(ilayer)) = ruslag(nc,nb,ifrmch(ilayer)) * ettp(ip,jmp) / float(nlayer)
              ettp(ip,jmck(ilayer)) = ruslag(nc,nb,ifrmck(ilayer)) * ettp(ip,jmp) / float(nlayer)

            enddo

            ! Remplissage de TEPA
            tepa(ip,jrdck) = ruslag(nc,nb,irdck)
            tepa(ip,jrd0p) = ruslag(nc,nb,ird0p)

            do ilayer = 1, nlayer

              tepa(ip,jrhock(ilayer))= ruslag(nc,nb,irhock0(ilayer))

            enddo

          ! composition from DP_FCP
          else if (iuslag(nc,nb,irawcl).eq.1) then

            ! Remplissage de ETTP
            ettp(ip,jcp) = cp2ch(iuslag(nc,nb,inuchl))
            ettp(ip,jmp) = rho0ch(iuslag(nc,nb,inuchl)) * pis6 * d3
            ettp(ip,jmwat) = xwatch(iuslag(nc,nb,inuchl)) * ettp(ip,jmp)

            do ilayer = 1, nlayer

              ettp(ip,jmch(ilayer)) =                                           &
              (1.0d0-xwatch(iuslag(nc,nb,inuchl))-xashch(iuslag(nc,nb,inuchl))) &
                   * ettp(ip,jmp) / float(nlayer)

              ettp(ip,jmck(ilayer)) = 0.0d0

            enddo

            ! Remplissage de TEPA
            tepa(ip,jrdck) = ettp(ip,jdp)
            tepa(ip,jrd0p) = ettp(ip,jdp)

            do ilayer=1,nlayer
              tepa(ip,jrhock(ilayer)) = rho0ch(iuslag(nc,nb,inuchl))
            enddo

         endif
       endif

!--> POIDS STATISTIQUE

       if (iuslag(nc,nb,ijprpd).eq.1) then
          tepa(ip,jrpoi) = ruslag(nc,nb,ipoit)
       else if (iuslag(nc,nb,ijprpd).eq.2) then

          idvar = 0
          xxpart = ettp(ip,jxp)
          yypart = ettp(ip,jyp)
          zzpart = ettp(ip,jzp)

          call uslapr                                             &
          !==========
 ( idvar  , iel    , nb     , nc     ,                            &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   xxpart , yypart , zzpart ,                                     &
   tvpart , uupart , vvpart , wwpart , ddpart , ttpart ,          &
   dt     , rtpa   , propce ,                                     &
   ettp   , tepa   )

          volp = pis6*d3
          surf = sqrt( surfbo(1,ifac)*surfbo(1,ifac)              &
                      +surfbo(2,ifac)*surfbo(2,ifac)              &
                      +surfbo(3,ifac)*surfbo(3,ifac) )
          vitp = sqrt( ettp(ip,jup)*ettp(ip,jup)                  &
                      +ettp(ip,jvp)*ettp(ip,jvp)                  &
                      +ettp(ip,jwp)*ettp(ip,jwp) )
          tepa(ip,jrpoi) =tvpart*(surf*vitp*dtp)/volp

         endif

! Modele de Deposition : Initialisation

         if ( idepst .eq. 1 ) then

           call zufall(1,dintrf(1))

           tepa(ip,jrinpf) = 5.d0 + 15.d0 * dintrf(1)

           tepa(ip,jryplu) = 1000.d0
           itepa(ip,jimark) = -1
           itepa(ip,jdiel)  = 0
           itepa(ip,jdfac)  = 0

         endif

      enddo

      npt = npt + iusloc(nc,nb,ijnbp)

    endif

  enddo
enddo

!-->TEST DE CONTROLE (NE PAS MODIFIER)

if ( (nbpart+nlocnew).ne.npt ) then
  write(nfecra,3010) nlocnew, npt-nbpart
  call csexit (1)
  !==========
endif


!===============================================================================
! 5. MODIFICATION DES POIDS POUR AVOIR LE DEBIT
!===============================================================================

!   reinitialisation du compteur de nouvelles particules

npt = nbpart

!     pour chaque zone de bord :

do ii = 1,nfrlag
   nb = ilflag(ii)

   !         pour chaque classe :

   do nc = 1,iusncl(nb)

      !         si de nouvelles particules sont entrees,
      !         et si on a un debit non nul :

      if ( mod(ntcabs,iuslag(nc,nb,ijfre)).eq.0 .and.               &
           ruslag(nc,nb,idebt) .gt. 0.d0        .and.               &
           iusloc(nc,nb,ijnbp) .gt. 0                 ) then

         if (irangp.ge.0) then
            rapsurf = dble(iusloc(nc,nb,ijnbp)) / iuslag(nc,nb,ijnbp)
         else
            rapsurf = 1.d0
         endif

         dmasse = 0.d0
         do ip = npt+1 , npt + iusloc(nc,nb,ijnbp)
            dmasse = dmasse + ettp(ip,jmp)
         enddo

         !        Calcul des Poids

         if ( dmasse.gt.0.d0 ) then
            do ip = npt+1 , npt+iusloc(nc,nb,ijnbp)
               tepa(ip,jrpoi) = ( ruslag(nc,nb,idebt)*dtp )  * rapsurf &
                    / dmasse
            enddo
         else
            write(nfecra,1057) nb, nc, ruslag(nc,nb,idebt),           &
                 iusloc(nc,nb,ijnbp)
            call csexit (1)
            !==========
         endif


         npt = npt +  iusloc(nc,nb,ijnbp)

      endif

   enddo

enddo


!-->TEST DE CONTROLE (NE PAS MODIFIER)

! FIXME : the following test seems flawed
!if ( (nbpart+nlocnew).ne.npt ) then
!  write(nfecra,3010) nlocnew, npt-nbpart
!  call csexit (1)
!  !==========
!endif



!===============================================================================
! 6. SIMULATION DES VITESSES TURBULENTES FLUIDES INSTANTANEES VUES
!    PAR LES PARTICULES SOLIDES LE LONG DE LEUR TRAJECTOIRE.
!===============================================================================

!   si de nouvelles particules doivent entrer :

npar1 = nbpart+1
npar2 = nbpart+ nlocnew

call lagipn(nbpmax, npar1, npar2, rtpa, vagaus, propce)
!==========

!===============================================================================
! 7. MODIFICATION DES TABLEAUX DE DONNEES PARTICULAIRES
!===============================================================================

call uslain                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nlocnew ,                                                      &
   itypfb , itrifb , itepa  , ifrlag , iwork  ,                   &
   dt     , rtpa   , propce ,                                     &
   ettp   , tepa   , vagaus , icocel , lndnod , itycel , dlgeo,   &
   ncmax  , nzmax  , iusloc )


!===============================================================================
! 7bis. Random id associated with particles (to be initialized later)
!===============================================================================

do npt = npar1,npar2
  call random_number(tepa(npt,jrval))
enddo

!   reinitialisation du compteur de nouvelles particules
npt = nbpart

!     pour chaque zone de bord:
do ii = 1,nfrlag
  nb = ilflag(ii)

!       pour chaque classe :
  do nc = 1, iusncl(nb)

!         si de nouvelles particules doivent entrer :
    if (mod(ntcabs,iusloc(nc,nb,ijfre)).eq.0) then

      do ip = npt+1 , npt+iusloc(nc,nb,ijnbp)

        if (ettp(ip,jdp).lt.0.d0 .and.                            &
            ruslag(nc,nb,ivdpt).gt.0.d0) then
          write(nfecra,4000) ruslag(nc,nb,idpt),                  &
                             ruslag(nc,nb,ivdpt),                 &
                             ettp(ip,jdp)
        endif

      enddo

      npt = npt + iusloc(nc,nb,ijnbp)

    endif

  enddo
enddo

!===============================================================================
! 8. INJECTION "CONTINUE" EVENTUELLE
!===============================================================================

if ( injcon.eq.1 ) then

 write(nfecra,*) "Error : pseudo-continuous injection not implemented "
 call csexit (1)

!  FIXME : Reimplementation of the continuous injection with the
!          the new cs_lagr_tracking.c routines

!
!   reinitialisation du compteur de nouvelles particules

!!$  npt = nbpart
!!$
!!$!       pour chaque zone de bord:
!!$
!!$  do ii = 1,nfrlag
!!$    nb = ilflag(ii)
!!$
!!$!         pour chaque classe :
!!$
!!$    do nc = 1, iusncl(nb)
!!$
!!$!           si de nouvelles particules doivent entrer :
!!$      if ( mod(ntcabs,iuslag(nc,nb,ijfre)).eq.0 ) then
!!$
!!$        call lagnwc                                               &
!!$        !==========
!!$  ( lndnod ,                                                      &
!!$    nbpmax , nvp    , nvp1   , nvep   , nivep  ,                  &
!!$    npt    , nlocnew , iuslag(nc,nb,ijnbp)      ,                  &
!!$    itycel , icocel ,                                             &
!!$    ifrlag , itepa(1,jisor)  , iwork  ,                           &
!!$    ettp   )
!!$
!!$      endif
!!$
!!$    enddo
!!$  enddo

endif

!-->TEST DE CONTROLE (NE PAS MODIFIER)

if ( (nbpart+nlocnew).ne.npt ) then
  write(nfecra,3010) nlocnew, npt-nbpart
  call csexit (1)
  !==========
endif

! Free memory
deallocate(iwork)

!===============================================================================
! 9. CALCUL DE LA MASSE TOTALE INJECTES EN CHAQUE ZONE
!    Attention cette valeur est modifie dans USLABO pour tenir compte
!    des particules qui sortent
!    + calcul du nombres physiques de particules qui rentrent (tenant
!       compte des poids)
!===============================================================================

!   reinitialisation du compteur de nouvelles particules

npt     = nbpart
dnbpnw = 0.d0

!     pour chaque zone de bord :

do ii = 1,nfrlag
  nb = ilflag(ii)
  deblag(nb) = 0.d0

!       pour chaque classe :

  do nc = 1,iusncl(nb)

!        si de nouvelles particules sont entrees,

    if ( mod(ntcabs,iuslag(nc,nb,ijfre)).eq.0 .and.               &
             iusloc(nc,nb,ijfre).gt.0            ) then

      do ip = npt+1 , npt+iusloc(nc,nb,ijnbp)
        deblag(nb) = deblag(nb) + tepa(ip,jrpoi)*ettp(ip,jmp)
        dnbpnw = dnbpnw + tepa(ip,jrpoi)
      enddo

    endif

    npt = npt + iusloc(nc,nb,ijnbp)

  enddo

enddo

!===============================================================================
! 10. NOUVEAU NOMBRE DE PARTICULES TOTAL
!===============================================================================

!     NBPART : NOMBRE DE PARTICULES PRESENTES DANS LE DOMAINE

!     NBPTOT : NOMBRE DE PARTICULES TOTAL INJECTE DANS
!              LE CALCUL DEPUIS LE DEBUT SUITE COMPRISE


nbpart = nbpart + nlocnew
dnbpar = dnbpar + dnbpnw

if (irangp.ge.0) then
  call parcpt(nlocnew)
  !==========
endif

nbptot = nbptot + nlocnew

!===============================================================================

deallocate(iusloc)

if (irangp.ge.0) then
   deallocate(surflag)
   deallocate(surlgrg)
   deallocate(ninjrg)
   deallocate(ilftot)
endif

!--------
! FORMATS
!--------


 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES ',/,&
'@                                                            ',/,&
'@  Le numero de zone associee a la face ',I10   ,' doit etre ',/,&
'@    un entier strictement positif et inferieur ou egal a    ',/,&
'@    NFLAGM = ',I10                                           ,/,&
'@  Ce numero (IFRLAG(IFAC)) vaut ici ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@  Le nombre maximal de zones frontieres qui peuvent etre    ',/,&
'@    definies par l''utilisateur est NFLAGM = ',I10           ,/,&
'@    Il a ete depasse.                                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@  Les zones frontieres NFLAGM premieres zones frontieres    ',/,&
'@    portent ici les numeros suivants :                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES ',/,&
'@                                                            ',/,&
'@  Le nombre de classes de la zone numero ',I10   ,' doit    ',/,&
'@    etre un entier positif ou nul.                          ',/,&
'@  Ce nombre (IUSNCL(NB)  ) vaut ici ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le nombre de particules dans la classe         ',I10       ,/,&
'@                          dans la zone frontiere ',I10       ,/,&
'@  doit etre un entier strictement positif et                ',/,&
'@                      inferieur ou egal a NBPMAX = ',I10     ,/,&
'@                                                            ',/,&
'@  Ce nombre (IUSLAG(NC,NB,IJNBP)) vaut ici ',I10             ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1021 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME A L''EXECUTION DU MODULE LAGRANGIEN',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le numero de groupe statistique de particules             ',/,&
'@     dans la classe         NC= ',I10                        ,/,&
'@     dans la zone frontiere NB= ',I10                        ,/,&
'@  doit etre un entier strictement positif et                ',/,&
'@                      inferieur ou egal a NBCLST = ',I10     ,/,&
'@                                                            ',/,&
'@  Ce nombre (IUSLAG(NC,NB,IJNBP)) vaut ici ',I10             ,/,&
'@                                                            ',/,&
'@  Le calcul continue mais cette classe statistique sera     ',/,&
'@  ignoree                                                   ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  La frequence d''injection des particules doit etre un     ',/,&
'@    entier positif ou nul (nul signifiant que               ',/,&
'@    les particules ne injectees qu''au debut du calcul).    ',/,&
'@                                                            ',/,&
'@  Ce nombre pour la frontiere NB  =',I10                     ,/,&
'@    et      pour la classe    NC  =',I10                     ,/,&
'@    vaut ici IUSLAG (NC,NB,IJFRE) =',I10                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1040 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Les conditions aux bords sont representees par            ',/,&
'@    une variable dont les valeurs sont                      ',/,&
'@    obligatoirement les suivantes :                         ',/,&
'@                                                            ',/,&
'@   = IENTRL  zone d''injection de particules                ',/,&
'@   = ISORTL  sortie du domaine                              ',/,&
'@   = IREBOL  rebond des particules                          ',/,&
'@   = ISYMTL  flux nul pour les particules  (symetrie)       ',/,&
'@   = IDEPO1  deposition definitive                          ',/,&
'@   = IDEPO2  deposition definitive mais la particule reste  ',/,&
'@             en memoire                                     ',/,&
'@   = IENCRL  encrassement (Charbon uniquement IPHYLA = 2)   ',/,&
'@   = JBORD1  interaction particule/frontiere utilisateur    ',/,&
'@   = JBORD2  interaction particule/frontiere utilisateur    ',/,&
'@   = JBORD3  interaction particule/frontiere utilisateur    ',/,&
'@   = JBORD4  interaction particule/frontiere utilisateur    ',/,&
'@   = JBORD5  interaction particule/frontiere utilisateur    ',/,&
'@                                                            ',/,&
'@  Cette valeur pour la frontiere NB = ',I10                  ,/,&
'@     est erronees.                                          ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1042 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  La condition a la limite representee par le type          ',/,&
'@    IENCRL n''est admissible que lorsque les particules     ',/,&
'@    transportees sont des grains de charbon IPHYLA = 2      ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IPHYLA dans USLAG1 et les           ',/,&
'@  conditions aux limites pour la frontiere NB =',I10         ,/,&
'@  dans USLAG2.                                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1043 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  La condition a la limite representee par le type          ',/,&
'@    ENCRAS n''est admissible que lorsque l''option          ',/,&
'@    encrassement est enclenche IENCRA = 1                   ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IENCRA dans USLAG1 et les           ',/,&
'@  conditions aux limites pour la frontiere NB =',I10         ,/,&
'@  dans USLAG2.                                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1050 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le type de condition au bord pour la vitesse est          ',/,&
'@    represente par un entier dont les valeurs sont          ',/,&
'@    obligatoirement les suivantes                           ',/,&
'@   =-1 vitesse fluide imposee                               ',/,&
'@   = 0 vitesse imposee selon la direction normale a la      ',/,&
'@       face de bord et de norme RUSLAG(NB,NC,IUNO)          ',/,&
'@   = 1 vitesse imposee : on donne RUSLAG(NB,NC,IUPT)        ',/,&
'@                                  RUSLAG(NB,NC,IVPT)        ',/,&
'@                                  RUSLAG(NB,NC,IWPT)        ',/,&
'@   = 2 profil de vitesse imposee dans USLAPR                ',/,&
'@                                                            ',/,&
'@  Ce nombre pour la frontiere NB = ',I10                     ,/,&
'@    et      pour la classe    NC = ',I10                     ,/,&
'@    vaut ici IUSLAG(NC,NB,IJUVW) = ',I10                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@  (Si IUSLAG(NC,NB,IJUVW) vaut -2 il n''a pas ete renseigne)',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1051 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le type de condition au bord pour le diametre est         ',/,&
'@    represente par un entier dont les valeurs sont          ',/,&
'@    obligatoirement les suivantes                           ',/,&
'@   = 1 diametre imposee dans la zone :                      ',/,&
'@                  on donne        RUSLAG(NB,NC,IDPT)        ',/,&
'@                                  RUSLAG(NB,NC,IVDPT)       ',/,&
'@   = 2 profil de diametre imposee dans USLAPR               ',/,&
'@                                                            ',/,&
'@  Ce nombre pour la frontiere  NB = ',I10                    ,/,&
'@    et      pour la classe     NC = ',I10                    ,/,&
'@    vaut ici IUSLAG(NC,NB,IJPRDP) = ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@ (Si IUSLAG(NC,NB,IJPRDP) vaut -2 il n''a pas ete renseigne)',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1052 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le type de condition au bord pour la temperature est      ',/,&
'@    represente par un entier dont les valeurs sont          ',/,&
'@    obligatoirement les suivantes                           ',/,&
'@   = 1 temperature imposee dans la zone :                   ',/,&
'@                  on donne        RUSLAG(NB,NC,ITPT))       ',/,&
'@   = 2 profil de temperature imposee dans USLAPR            ',/,&
'@                                                            ',/,&
'@  Ce nombre pour la frontiere  NB = ',I10                    ,/,&
'@    et      pour la classe     NC = ',I10                    ,/,&
'@    vaut ici IUSLAG(NC,NB,IJPRTP) = ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@ (Si IUSLAG(NC,NB,IJPRTP) vaut -2 il n''a pas ete renseigne)',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1053 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le type de condition au bord pour la taux de presence est ',/,&
'@    represente par un entier dont les valeurs sont          ',/,&
'@    obligatoirement les suivantes                           ',/,&
'@   = 1 distribution uniforme                                ',/,&
'@                  on donne        RUSLAG(NB,NC,IPOID))      ',/,&
'@   = 2 profil de taux de presence imposee dans USLAPR       ',/,&
'@                                                            ',/,&
'@  Ce nombre pour la frontiere  NB = ',I10                    ,/,&
'@    et      pour la classe     NC = ',I10                    ,/,&
'@    vaut ici IUSLAG(NC,NB,IJPRPD) = ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@ (Si IUSLAG(NC,NB,IJPRPD) vaut -2 il n''a pas ete renseigne)',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1055 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le poids statistique des particules doit etre un reel     ',/,&
'@    strictement positif.                                    ',/,&
'@                                                            ',/,&
'@  Ce nombre pour la frontiere NB  =',I10                     ,/,&
'@    et      pour la classe    NC  =',I10                     ,/,&
'@    vaut ici RUSLAG (NC,NB,IPOIT) =',E14.5                   ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1056 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le debit massique des particules doit etre un reel        ',/,&
'@    positif ou nul (nul signifiant que le debit n''est pas  ',/,&
'@    pris en compte).                                        ',/,&
'@                                                            ',/,&
'@  Ce nombre pour la frontiere NB  =',I10                     ,/,&
'@    et      pour la classe    NC  =',I10                     ,/,&
'@    vaut ici RUSLAG (NC,NB,IDEBT) =',E14.5                   ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1057 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    ne sont pas physiques.                                  ',/,&
'@  Le debit massique impose                                  ',/,&
'@    vaut ici RUSLAG (NC,NB,IDEBT) =',E14.5                   ,/,&
'@  alors que le nombre de particules injectees est nul       ',/,&
'@             IUSLAG (NC,NB,IJNBP) =',I10                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    ne sont pas physiques :                                 ',/,&
'@    Masse volumique : RUSLAG(NC,NB,IROPT) = ',E14.5          ,/,&
'@    Diametre moyen  : RUSLAG(NC,NB,IDPT)  = ',E14.5          ,/,&
'@    Ecart type      : RUSLAG(NC,NB,IVDPT) = ',E14.5          ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@  Verifier le fichier dp_FCP si l''option Charbon pulverise ',/,&
'@    est activee.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1065 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  L''ecart-type fourni est trop grand par rapport           ',/,&
'@    au diametre moyen pour                                  ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@                                                            ',/,&
'@  Il y a un risque non nul de calcul d''un diametre negatif.',/,&
'@                                                            ',/,&
'@  Theoriquement 99,7% des particules se trouvent entre :    ',/,&
'@    RUSLAG(NC,NB,IDPT)-3*RUSLAG(NC,NB,IVDPT)                ',/,&
'@    RUSLAG(NC,NB,IDPT)+3*RUSLAG(NC,NB,IVDPT)                ',/,&
'@  Pour eviter des diametres aberrants, dans le module       ',/,&
'@    lagrangien, avec un clipping, on impose que 100% des    ',/,&
'@    particules doivent etre dans cet intervalle.            ',/,&
'@                                                            ',/,&
'@  Or on a :                                                 ',/,&
'@    RUSLAG(NC,NB,IDPT)-3*RUSLAG(NC,NB,IVDPT) = ',E14.5       ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@  Verifier le fichier dp_FCP si l''option Charbon pulverise ',/,&
'@    est activee.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1070 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Une equation sur la temperature est associee              ',/,&
'@    aux particules (IPHYLA = ',I10,') :                     ',/,&
'@    ITPVAR = ',I10                                           ,/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    doivent etre renseignees :                              ',/,&
'@    Temperature  : RUSLAG(NC,NB,ITPT) = ',E14.5              ,/,&
'@    Cp           : RUSLAG(NC,NB,ICPT) = ',E14.5              ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1075 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Une equation sur la temperature est associee              ',/,&
'@    aux particules (IPHYLA = ',I10,') :                     ',/,&
'@    ITPVAR = ',I10                                           ,/,&
'@    avec prise em compte des echanges thermiques radiatifs. ',/,&
'@  L''emissivite des particules doit etre renseignee et      ',/,&
'@    comprise entre 0 et 1 (inclus).                         ',/,&
'@                                                            ',/,&
'@  L''emissivite pour la frontiere NB = ',I10                 ,/,&
'@                     et la classe NC = ',I10                 ,/,&
'@    vaut : RUSLAG(NC,NB,IEPSI) = ',E14.5                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1080 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@    L''INDICATEUR SUR LE NUMERO DU CHARBON                  ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@  Le numero du charbon injecte pour                         ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    devrait etre compris entre 1 et NCHARB= ',I10            ,/,&
'@    Le nombre de charbon NCHARB est donne dans dp_FCP.      ',/,&
'@                                                            ',/,&
'@    Il vaut ici  : IUSLAG(NC,NB,INUCHL) = ',I10              ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2 et le fichier dp_FCP.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1081 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le type de condition au bord pour la composition du       ',/,&
'@    charbon est represente par un entier dont les valeurs   ',/,&
'@    sont obligatoirement les suivantes                      ',/,&
'@   = 0 composition définie par l''utilisateur               ',/,&
'@   = 1 composition prise égale à celle du charbon frais     ',/,&
'@                                                            ',/,&
'@  Ce nombre pour la frontiere NB  = ',I10                    ,/,&
'@    et      pour la classe    NC  = ',I10                    ,/,&
'@    vaut ici IUSLAG(NC,NB,IRAWCL) = ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2.                                          ',/,&
'@  (Si IUSLAG(NC,NB,IRAWCL) vaut -2 il n''a pas ete          ',/,&
'@   renseigne)                                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1082 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT DANGEREUSES             ',/,&
'@                                                            ',/,&
'@    IL Y A UN RISQUE DANS LES CONDITIONS AU BORD            ',/,&
'@       DEFINIES (LAGOPT).                                   ',/,&
'@                                                            ',/,&
'@    IRAWCL = 0 correspond à une composition défini par      ',/,&
'@    l''utilisateur                                          ',/,&
'@    IJPRDP = 2 correspond à un diametre défini par          ',/,&
'@    l''utilisateur                                          ',/,&
'@                                                            ',/,&
'@  Ces deux options sont dangereuses car on définit          ',/,&
'@    RUSLAG(NC,NB,IRDCK) et RUSLAG(NC,NB,IRD0P) dans USLAG2  ',/,&
'@    avant meme de connaitre le diametre de la particule     ',/,&
'@    defini dans USLAPR                                      ',/,&
'@                                                            ',/,&
'@    Pour la frontiere NB  = ',I10                            ,/,&
'@    et pour la classe NC  = ',I10                            ,/,&
'@    IUSLAG(NC,NB,IRAWCL)  = ',I10                            ,/,&
'@    IUSLAG(NC,NB,IJPRDP)  = ',I10                            ,/,&
'@                                                            ',/,&
'@  Le calcul continu mais une verification est souhaitable.  ',/,&
'@                                                            ',/,&
'@  Verifier que les expressions du diametre de la particule  ',/,&
'@    dans USLAPR et la donnee de RUSLAG(NC,NB,IRDCK) et      ',/,&
'@    RUSLAG(NC,NB,IRD0P) dans USLAG2 sont bien coherentes    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1083 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@    IL Y A INCOMPATIBILITE ENTRE LES CONDITIONS AU BORD     ',/,&
'@       DEFINIES (LAGOPT).                                   ',/,&
'@                                                            ',/,&
'@    IRAWCL = 0 correspond à une composition défini par      ',/,&
'@    l''utilisateur                                          ',/,&
'@    IVDPT  > 0.0 correspond à un diametre variable          ',/,&
'@                                                            ',/,&
'@  Ces deux options sont incompatibles car on définit        ',/,&
'@    RUSLAG(NC,NB,IRDCK) et RUSLAG(NC,NB,IRD0P) dans USLAG2  ',/,&
'@    avant meme de connaitre le diametre de la particule     ',/,&
'@    tire aleatoirement dans USLAG2                          ',/,&
'@                                                            ',/,&
'@    Pour la frontiere NB  = ',I10                            ,/,&
'@    et pour la classe NC  = ',I10                            ,/,&
'@    IUSLAG(NC,NB,IRAWCL)  = ',I10                            ,/,&
'@    RUSLAG(NC,NB,IVDPT) = ',E14.5                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IUSLAG(NC,NB,IRAWCL) et RUSLAG(NC,NB,IVDPT)      ',/,&
'@      dans USLAG2, on peut avoir:                           ',/,&
'@              IRAWCL = 1  et  IVDPT > 0.0                   ',/,&
'@           ou IRAWCL = 0  et  IVDPP = 0.0                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1084 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  L''option de transport de particules de charbon pulverise ',/,&
'@    est active (IPHYLA = ',I10,')                           ',/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    et la couche ILAYER = ',I10                              ,/,&
'@    doivent etre renseignees :                              ',/,&
'@    Temperature   :                                         ',/,&
'@      RUSLAG(NC,NB,IHPT(ILAYER))    = ',E14.5                ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2 ou le DP_FCP                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1090 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  L''option de transport de particules de charbon pulverise ',/,&
'@    est active (IPHYLA = ',I10,')                           ',/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    doivent etre renseignees :                              ',/,&
'@    Masse volumique :                                       ',/,&
'@      RUSLAG(NC,NB,IROPT)  = ',E14.5                         ,/,&
'@    Cp :                                                    ',/,&
'@      RUSLAG(NC,NB,ICPT)   = ',E14.5                         ,/,&
'@    Fraction massique d'' eau dans la particule :           ',/,&
'@      RUSLAG(NC,NB,IFRMWT) = ',E14.5                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2                                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1091 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  L''option de transport de particules de charbon pulverise ',/,&
'@    est active (IPHYLA = ',I10,')                           ',/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    et la couche ILAYER = ',I10                              ,/,&
'@    doivent etre renseignees :                              ',/,&
'@    Fraction massique de charbon reactif :                  ',/,&
'@      RUSLAG(NC,NB,IFRMCH(ILAYER))  = ',E14.5                ,/,&
'@    Fraction massique de coke :                             ',/,&
'@      RUSLAG(NC,NB,IFRMCK(ILAYER))  = ',E14.5                ,/,&
'@    Masse vol de coke juste après pyrolyse :                ',/,&
'@      RUSLAG(NC,NB,IRHOCK0(ILAYER)) = ',E14.5                ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2 et le fichier dp_FCP.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1092 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  L''option de transport de particules de charbon pulverise ',/,&
'@    est active (IPHYLA = ',I10,') et la composition des     ',/,&
'@    particules est donnée par l''utilisateur                ',/,&
'@    IUSLAG(NC,NB,IRAWCL) = ',I10                             ,/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    doivent etre renseignees :                              ',/,&
'@    Diam de coke : RUSLAG(NC,NB,IRDCK) = ',E14.5             ,/,&
'@    Diam initial : RUSLAG(NC,NB,IRD0P) = ',E14.5             ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2 et le fichier dp_FCP.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1093 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  L''option de transport de particules de charbon pulverise ',/,&
'@    est active (IPHYLA = ',I10,')                           ',/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    doivent etre renseignees :                              ',/,&
'@    Numéro de charbon:                                      ',/,&
'@      ICOAL =         ',I10                                  ,/,&
'@    Masse volumique :                                       ',/,&
'@      RHO0CH(ICOAL) = ',E14.5                                ,/,&
'@    Cp :                                                    ',/,&
'@      CP2CH(ICOAL)  = ',E14.5                                ,/,&
'@    Fraction massique d'' eau dans la particule :           ',/,&
'@      XWATCH(ICOAL) = ',E14.5                                ,/,&
'@    Fraction massique de cendres dans la particule :        ',/,&
'@      XASHCH(ICOAL) = ',E14.5                                ,/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier DP_FCP                                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1094 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  L''option de transport de particules de charbon pulverise ',/,&
'@    est active (IPHYLA = ',I10,')                           ',/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    ne sont pas correctes :                                 ',/,&
'@    Numéro de charbon:                                      ',/,&
'@      ICOAL =         ',I10                                  ,/,&
'@    Fraction massique d'' eau dans la particule :           ',/,&
'@      XWATCH(ICOAL) = ',E14.5                                ,/,&
'@    Fraction massique de cendres dans la particule :        ',/,&
'@      XASHCH(ICOAL) = ',E14.5                                ,/,&
'@  La fraction massique de cendre et d''eau est supérieure   ',/,&
'@    à 1:',E14.5                                              ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier DP_FCP                                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1095 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  L''option de transport de particules de charbon pulverise ',/,&
'@    est active (IPHYLA = ',I10,')                           ',/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    ne sont pas correctes :                                 ',/,&
'@                                                            ',/,&
'@  Le type de condition au limite avec recours au fichier    ',/,&
'@    DP_FCP a été sélectionné:                               ',/,&
'@      IUSLAG(NB,NC,IRAWCL) = ',I10                           ,/,&
'@                                                            ',/,&
'@  Pourtant un des ces grandeur a été initialisée (valeur    ',/,&
'@    différente de ',E14.5,')                                ',/,&
'@    Masse volumique de la particule :                       ',/,&
'@      RUSLAG(NC,NB,IROPT)  = ',E14.5                         ,/,&
'@    Fraction massique d'' eau dans la particule :           ',/,&
'@      RUSLAG(NC,NB,IFRMWT) = ',E14.5                         ,/,&
'@    Cp :                                                    ',/,&
'@      RUSLAG(NC,NB,ICPT)   = ',E14.5                         ,/,&
'@    Diametre de coke :                                      ',/,&
'@      RUSLAG(NC,NB,IRDCK)  = ',E14.5                         ,/,&
'@    Diamètre initial :                                      ',/,&
'@      RUSLAG(NC,NB,IRDCK)  = ',E14.5                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2                                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1096 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN OPTION CHARBON PULVERISE  ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  L''option de transport de particules de charbon pulverise ',/,&
'@    est active (IPHYLA = ',I10,')                           ',/,&
'@  Les proprietes physiques des particules au bord pour      ',/,&
'@    la frontiere NB = ',I10   ,' et la classe NC = ',I10     ,/,&
'@    pour la couche ILAYER = ',I10                            ,/,&
'@    ne sont pas correctes :                                 ',/,&
'@                                                            ',/,&
'@  Le type de condition au limite avec recours au fichier    ',/,&
'@    DP_FCP a été sélectionné:                               ',/,&
'@      IUSLAG(NB,NC,IRAWCL) = ',I10                           ,/,&
'@                                                            ',/,&
'@  Pourtant un des ces grandeur a été initialisée (valeur    ',/,&
'@    différente de ',E14.5,')                                ',/,&
'@    Fraction massique de charbon réactif dans la particule :',/,&
'@      RUSLAG(NC,NB,IFRMCH(ILAYER)) = ',E14.5                 ,/,&
'@    Fraction massique de coke dans la particule :'           ,/,&
'@      RUSLAG(NC,NB,IFRMCK(ILAYER)) = ',E14.5                 ,/,&
'@    Masse volumique initiale du coke :                      ',/,&
'@      RUSLAG(NC,NB,IRHOCK0(ILAYER))  = ',E14.5               ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2                                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE LAGRANGIEN                           ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le nombre de nouvelles particules injectees conduit a un  ',/,&
'@    nombre total de particules superieur au maximum prevu : ',/,&
'@                                                            ',/,&
'@    Nombre de particules courant   : NBPART = ',I10          ,/,&
'@    Nombre de nouvelles particules : NBPNEW = ',I10          ,/,&
'@    Nombre maximal de particules   : NBPMAX = ',I10          ,/,&
'@                                                            ',/,&
'@  Le calcul se poursuit, mais on n''injecte aucune particule',/,&
'@                                                            ',/,&
'@  Ajuster NBPMAX dans USLAG1 et verifier USLAG2.            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le nombre de particules injectees dans le domaine         ',/,&
'@    pour cette iteration Lagrangienne ne correspond pas     ',/,&
'@    a celui specifie dans les conditions aux limites.       ',/,&
'@                                                            ',/,&
'@  Nombre de particules specifie pour l''injection :         ',/,&
'@    NBPNEW = ',I10                                           ,/,&
'@  Nombre de particules effectivement injectees :            ',/,&
'@    NPT-NBPART = ',I10                                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Contacter l''equipe de developpement.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 4000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGENT)                                    ',/,&
'@                                                            ',/,&
'@    LES CONDITIONS AUX LIMITES SONT ERRONEES                ',/,&
'@                                                            ',/,&
'@  Le reconstruction du duametre d''une particule a partir   ',/,&
'@    du diametre moyen et de l''ecart type donne une valeur  ',/,&
'@    de diametre negatif, a cause d''un tirage aleatoire     ',/,&
'@    dans un des "bord de la Gaussienne".                    ',/,&
'@                                                            ',/,&
'@  Diametre moyen : ',E14.5                                   ,/,&
'@  Ecart type : ',E14.5                                       ,/,&
'@  Diametre calcule : ',E14.5                                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Mettre en place une verification du diametre en fonction  ',/,&
'@  des donnees granulometriques dans USLAIN.                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return
end subroutine lagent
