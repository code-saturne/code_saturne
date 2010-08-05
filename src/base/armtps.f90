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

subroutine armtps &
!================

 ( ntcabs , ntmabs )

!===============================================================================

! FONCTION :
! --------
!          ROUTINE PERMETTANT DE STOPPER LE CALCUL PROPREMENT SI
!            LE TEMPS RESTANT ALLOUE AU PROCESS EST INSUFFISANT
!            UNIQUEMENT POUR VPP ET CLUSTER LINUX

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ntcabs           ! e  ! <-- ! numero absolu du pas de temps courant          !
! ntmabs           ! e  ! <-- ! numero absolu du pas de temps final            !
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
include "entsor.f90"
include "parall.f90"

!===============================================================================

! Arguments

integer          ntcabs , ntmabs

! Local variables

integer irangs, lng, itmp(1), itlim
integer imetho
data    imetho /-1/
save    imetho

integer ntcab0
save    ntcab0
double precision trest0, trestp, tcpupr
save             trest0, trestp, tcpupr

double precision tmoy00, titpre, trestc, tmoyit, alpha, titsup
double precision tmarge, aa, bb, cc, tmamin
double precision tcpuco
double precision tresmn, titsmx

!===============================================================================

!     La variable IMETHO indique si l'on a utilise avec succes
!     TREMAIN (1), TCPUMX (2) pour determiner le temps CPU limite,
!     ou s'il ne semble  avoir aucune limite ou que les deux methodes
!     on echoué (0). Avant le premier passage, elle vaut -1.

if (imetho.ne.0) then

!===============================================================================
! 1. AU PREMIER PASSAGE : INITIALISATIONS
!===============================================================================

  if (imetho.eq.-1) then

! ---     Premier passage : on essaie d'abord TREMAI
!           Si l'on obtient pas de temps CPU limite, on essaie
!           TCPUMX, qui se base sur la presence de la variable
!           d'environnement CS_MAXTIME.

    call tremai(trest0, itlim)
    !==========
    if (itlim.eq.1) then
      imetho = 1
    else
      call tcpumx(trest0, itlim)
      !==========
      if (itlim.eq.1) then
        imetho = 2
      endif
    endif

!         Si une des methodes fonctionne et indique une limite :

    if (imetho.ge.0) then

      ntcab0 = ntcabs

! ---       Temps restant et temps CPU a l'iteration courante
!           (qui sera ensuite la precedente)

      trestp = trest0

      call dmtmps(tcpupr)
      !==========

    endif

  else

!===============================================================================
! 2. TEMPS MOYEN PAR ITERATION
!===============================================================================

! --- Temps de l'iteration precedente

    call dmtmps(tcpuco)
    !==========
    titpre = tcpuco-tcpupr


! --- Temps restant courant (temps restant pour le process en cours,
!      le temps des autres processes, compilation etc. non compris)
! --- Temps moyen par iteration depuis le debut

    if (imetho.eq.1) then
      call tremai(trestc, itlim)
      !==========
      tmoy00 = (trest0-trestc)/dble(ntcabs-ntcab0)
    else if (imetho.eq.2) then
!           Ici on utilise le temps alloue initialement
      trestc = max(trest0 - tcpuco,0.d0)
      tmoy00 = tcpuco/dble(ntcabs-ntcab0)
    endif

! --- Estimation du temps par iteration
!       ALPHA -> 0 EST PLUS SUR

    alpha = 0.25d0
    tmoyit = alpha*titpre + (1.d0-alpha)*tmoy00

! --- Temps restant iteration courante (qui sera ensuite la precedente)

    trestp = trestc

! --- Temps CPU iteration courante (qui sera ensuite la precedente)

    tcpupr = tcpuco

!===============================================================================
! 3. TEMPS NECESSAIRE POUR UNE ITERATION DE PLUS
!===============================================================================

! --- Marge pour les sorties ...
!      100 fois une iteration ou 10% du temps alloue au process (-lt)
!        et au moins 50s ou 1% du temps alloue alloue au process (-lt)

!      Soit pour des jobs de
!        moins de    1000 iter     :  10% du temps alloue
!        plus  de    1000 iter et
!          moins de 10000 iter     : 100 fois une iter
!        plus  de   10000 iter     :   1% du temps alloue


    if (tmarus.lt.0.d0) then

      aa = 100.d0
      bb = 0.10d0
      cc = 0.01d0
      tmamin = 50.d0

      tmarge = min(tmoyit*aa,trest0*bb)
      tmarge = max(tmarge,tmamin)
      tmarge = max(tmarge,trest0*cc)

    else

      tmarge = tmarus

    endif


! --- Temps necessaire pour une iteration de plus

    titsup = tmoyit + tmarge

!===============================================================================
! 4. TEST (en parallele, le processeur 0 commande)
!===============================================================================

! On s'arrete si le temps restant minimum est inferieur au
!   temps maximum estime pour une iter de plus.
! Les notions de min et de max sont relatives aux calculs en parallele.
!   En effet, les temps sont estimes separement par les differents
!   processeurs et il faut qu'ils aient tous le temps de finir.
!   Avec cette methode la marge est la meme pour tous.
! Bien noter que UN SEUL processeur doit decider d'arreter le calcul.

    tresmn = trestc
    titsmx = titsup
    if(irangp.ge.0) then
      call parmin(tresmn)
      call parmax(titsmx)
          endif
    if(irangp.lt.0.or.irangp.eq.0) then
      if (tresmn.lt.titsmx) then
        ntmabs = ntcabs
        write(nfecra,1000) ntmabs
      endif
    else
      ntmabs = 0
    endif
! Broadcast
    if(irangp.ge.0) then
      irangs  = 0
      lng     = 1
      itmp(1) = ntmabs
      call parbci(irangs,lng,itmp)
      ntmabs = itmp(1)
    endif

    if (ntcabs.eq.ntmabs) then
      write(nfecra,1100) trestc, titsup, tmoy00, titpre, tmarge
    endif

  endif

endif


!===============================================================================
! 5. FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'===============================================================',&
/,'   ** ARRET PAR MANQUE DE TEMPS ',                             &
/,'      ------------------------- ',                             &
/,'      NOMBRE DE PAS DE TEMPS MAX IMPOSE A NTCABS : ',I10,    /,&
'===============================================================',&
                                                                /)

 1100 format(/,                                                   &
'===============================================================',&
/,'   ** GESTION DU TEMPS RESTANT ',                              &
/,'      ------------------------ ',                              &
/,'      TEMPS RESTANT ALLOUE AU PROCESS          : ',E14.5,      &
/,'      TEMPS ESTIME POUR UNE ITERATION DE PLUS  : ',E14.5,      &
/,'        DUREE MOYENNE D''UNE ITERATION EN TEMPS : ',E14.5,     &
/,'        DUREE DE L''ITERATION PRECEDENTE        : ',E14.5,     &
/,'        MARGE DE SECURITE                      : ',E14.5,   /, &
'===============================================================',&
                                                                /)

#else

 1000 format(/,                                                   &
'===============================================================',&
/,'   ** STOP BECAUSE OF TIME EXCEEDED'                           &
/,'      -----------------------------',                          &
/,'      MAX NUMBER OF TIME STEP SET TO NTCABS: ',I10,          /,&
'===============================================================',&
                                                                /)

 1100 format(/,                                                   &
'===============================================================',&
/,'   ** REMAINING TIME MANAGEMENT ',                             &
/,'      ------------------------- ',                             &
/,'      REMAINING TIME ALLOCATED TO THE PROCESS   : ',E14.5,     &
/,'      ESTIMATED TIME FOR ANOTHER TIME STEP      : ',E14.5,     &
/,'        MEAN TIME FOR A TIME STEP               : ',E14.5,     &
/,'        TIME FOR THE PREVIOUS TIME STEP         : ',E14.5,     &
/,'        SECURITY MARGIN                         : ',E14.5,   /,&
'===============================================================',&
                                                                /)

#endif

end subroutine
