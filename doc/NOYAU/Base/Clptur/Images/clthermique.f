          integer i
         double precision x
         do i = 1, 300
          x=float(i)/10.d0
          write(*,'(6E18.9)')x, x*x*x/1000
         enddo
c
c restera a placer deux points pour (a+at)/nu = 1/Pr = 1 
c      dont un a x = 0
c et a placer deux points pour (a+at)/nu = 0.42 x / Prt = 0.42 x  
c      dont un a x = 30 
         end
