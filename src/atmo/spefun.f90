!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file spefun.f90
!> \brief Atmospheric module - Module for specific math functions

module spefun
!> \defgroup at_spefun
!=============================================================================

use entsor

implicit none

!=============================================================================

contains

  !=============================================================================

  !> \brief  Gamma function
  function tgamma(x)

    double precision, intent(in) :: x
    double precision tgamma

    call csgamma(x, tgamma)

  end function tgamma

  !=============================================================================

  !> \brief Erf function
  function ferf(x)

    double precision, intent(in) :: x
    double precision ferf

    call cserf(x, ferf)

  end function ferf

  !=============================================================================

  !> \brief Hypergeometric function
  !> (see http://mathworld.wolfram.com/hypergeometricfunction.html
  !>   for definition)
  function hypgeo (a, b, c, x)

    double precision, intent(in) :: a, b, c, x
    double precision hypgeo

    ! Local variables

    double precision pp, y1, y2, hyp1, hyp2
    double precision gammaa, gammab, gammac
    double precision gammabma, gammaamb, gammacma, gammacmb
    parameter        (pp=0.1)

    ! Initialization

    gammaa   = tgamma(a)
    gammab   = tgamma(b)
    gammac   = tgamma(c)
    gammabma = tgamma(b-a)
    gammacma = tgamma(c-a)
    gammaamb = tgamma(a-b)
    gammacmb = tgamma(c-b)

    ! Compute hypergeometric function by convergent series for |x|<1

    if (x.ge.-1.+pp) then
      hypgeo = hypser(a, b, c, x)

    else if (x.le.-1-pp) then
      y1     = hypser( a, a+1.-c, a+1.-b, 1.d0/x )
      y2     = hypser( b, b+1.-c, b+1.-a, 1.d0/x )
      hypgeo = (gammac*gammabma*y1*(-1.d0/x)**a)/(gammab*gammacma)              &
             + (gammac*gammaamb*y2*(-1.d0/x)**b)/(gammaa*gammacmb)
    else
      y1     = hypser(a, a+1.-c, a+1.-b, 1.d0/(-1.-pp))
      y2     = hypser(b, b+1.-c, b+1.-a, 1.d0/(-1.-pp))
      hyp1   = (gammac*gammabma*y1*(-1.d0/(-1.-pp))**a)/(gammab*gammacma)       &
             + (gammac*gammaamb*y2*(-1.d0/(-1.-pp))**b)/(gammaa*gammacmb)
      hyp2   = hypser(a, b, c, -1.+pp)
      hypgeo = hyp1 + (x - (-1.-pp))*(hyp2 - hyp1)/(2.d0*pp)

    endif

  end function hypgeo

  !=============================================================================

  !> \brief  Calcul de la fonction hypergeometrique pour |x| < 1 par une serie
  !>     (cf. pour la definition de cette fonction, voir par exemple :
  !>     http://mathworld.wolfram.com/hypergeometricfunction.html )

  function hypser (a, b, c, x)

    use entsor

    implicit none

    double precision, intent(in) :: a, b, c, x
    double precision hypser

    !---------------------------------------------------------------------------

    ! Local variables

    integer          nn
    double precision fac,aa,bb,cc,temp
    integer,parameter :: maxiter=10000
    double precision,parameter :: error=1.d-08

    !---------------------------------------------------------------------------

    if (abs(x).ge.1.) then
      write (nfecra,1120) x
      call csexit(1)
    endif

    fac = 1
    temp = fac
    aa  = a
    bb  = b
    cc  = c

    do nn = 1, maxiter, 1
      fac    = ((aa*bb)/cc)*fac
      fac    = fac*x/nn
      hypser = fac + temp
      if (abs(hypser - temp).le. error) then
        return
      endif
      temp   = hypser
      aa     = aa +1
      bb     = bb +1
      cc     = cc +1
    enddo

    ! Formats
    !--------

#if defined(_CS_LANG_FR)

1120 format (                                                     &
'@                                                            ',/,&
'@ @@ ERREUR : dans la fonction hypser                        ',/,&
'@    ======                                                  ',/,&
'@  Le parametre x doit verifier |x| < 1,  x = ',E12.5         ,/,&
'@                                                            '  )

#else

 1120 format (                                                    &
'@                                                            ',/,&
'@ @@ ERROR: in hypser function                               ',/,&
'@    =====                                                   ',/,&
'@  The x parameter should verify |x| < 1,  x = ', E12.5       ,/,&
'@                                                            '  )

#endif

  end function hypser

  !=============================================================================

  !> \brief Beta function: beta(x,y) = gamma(x)*gamma(y)/gamma(x+y)
  function beta(x, y)

    double precision, intent(in) :: x, y
    double precision beta

    beta = tgamma(x)*tgamma(y)/tgamma(x + y)

  end function beta

  !=============================================================================

end module spefun
