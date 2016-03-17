!***********************************************************************
!
subroutine gauss(nn,mm,aa,xx,bb)
!***************

! solves a linear system with the gauss method

! **********************************************************************

! **********************************************************************

implicit none

integer          ilisti
parameter       (ilisti = 6)

integer          nn, mm
double precision aa(nn,mm), xx(nn), bb(nn)
double precision, allocatable, dimension (:,:) :: ap
double precision, allocatable, dimension (:)   :: bp

integer          ii, jj, kk, ll
integer          ipass, ipivot, ipermu, izero, ipivok
double precision akk, aik, aij, bi, sum, rr, residu

double precision time1, time2

! **********************************************************************

if (nn.ne.mm) then
  print*,'the matrix is not triangular -> stop'
  call csexit(1)
endif


! --- save the matrix which will be modified by the triangularisation

allocate (ap(nn,mm), bp(nn))

print*,"n x m",nn,mm

do jj = 1, nn
  do ii = 1, mm
    ap(ii,jj) = aa(ii,jj)
  enddo
  bp(jj) = bb(jj)
enddo

call cpu_time(time1)

! --- triangulise the matrix

ipass  = 0
ipivot = 0
ipermu = 0

kk = 1
do while (kk.lt.nn)


  ipass = ipass + 1
  if( ipass .gt. (nn) )then
    print*,'ipass > n',ipass,'>',nn,'-> stop'
    print*,'verify the Gauss pivot'
    print*,"gauss diagnostic:"
    print*,"pivoting on ",ipivot,"cells"
    print*,"permuting on ",ipermu,"cells"
    stop
  endif

  ! looking for the greatest coefficient
  ii = kk
  ll = 0
  akk = 0.d0
  do while (ii.le.nn)
    if (abs(ap(ii,kk)).gt.akk) then
      akk = abs(ap(ii,kk))
      ll = ii
    endif
    ii = ii + 1
  enddo

  if (ll.eq.0) then
    print*,"pas de pivot non nul => arret"
    stop
  else
    print*,"pivoting",kk,ll,akk
  endif

  ! pivoting line kk with line ll except if kk = ll
  if (ll.ne.kk) then
    do jj = 1, mm
      aij = ap(kk,jj)
      ap(kk,jj) = ap(ll,jj)
      ap(ll,jj) = aij
    enddo
    bi = bp(kk)
    bp(kk) = bp(ll)
    bp(ll) = bi
  endif

  ! scaling the pivot line
  akk = ap(kk,kk)
  do jj = kk, mm
    ap(kk,jj) = ap(kk,jj)/akk
  enddo
  bp(kk) = bp(kk)/akk

  ! elimination on the last lines
  do ii = kk+1, nn
    aik = ap(ii,kk)
    do jj = kk, mm
      ap(ii,jj) = ap(ii,jj) - aik*ap(kk,jj)
    enddo
    bp(ii) = bp(ii) - aik*bp(kk)
  enddo

  ! call prtmat(nn,mm,ap,bp)

  kk = kk + 1

enddo

! --- solves the triangular linear system

izero = 0
do ii = 1, nn
  if (ap(ii,ii).eq.0.d0) izero = ii
enddo
if (izero.ne.0) then
  write(ilisti,*) "Found zeros on diagonal of triangulised matrix"
  write(ilisti,*) "izero = ","ap(izero,izero) = "
  write(ilisti,*)  izero    , ap(izero,izero)
endif

xx(nn) = bp(nn)/ap(nn,nn)
do ii = nn-1, 1, -1
  sum = 0.d0
  do jj = ii+1, mm
    sum = sum + ap(ii,jj)*xx(jj)
  enddo
  xx(ii) = (bp(ii)-sum)/ap(ii,ii)
enddo

! --- Computes residu
residu = 0.d0
do ii = 1, nn
  rr = bb(ii)
  do jj = 1, nn
    rr = rr - aa(ii,jj)*xx(jj)
  enddo
  residu = residu + rr**2
enddo

call cpu_time(time2)

write(*     ,1000) '*** Gauss algorithm ***'
write(*     ,1002) 'Iterations: ',1
write(*     ,1003) 'Residu    : ',residu
write(*     ,1003) 'CPU time  : ',time2-time1
write(ilisti,1000) '*** Gauss algorithm ***'
write(ilisti,1002) 'Iterations: ',1
write(ilisti,1003) 'Residu    : ',residu
write(ilisti,1003) 'CPU time  : ',time2-time1


deallocate (ap, bp)

return

1000 format(/,a)
1001 format(4x,a12)
1002 format(4x,a12,i14)
1003 format(4x,a12,e14.5)

end subroutine gauss
