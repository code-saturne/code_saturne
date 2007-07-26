          integer i
         double precision x
         do i = 1, 2000
          x=float(i)/100.d0
          write(*,'(6E18.9)')x, 1/0.42*log(x)+5.2,x,
     &                          2/0.42*log(x)+5.2,2*x,2/0.42
         enddo
         end
