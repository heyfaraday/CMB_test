
         implicit real*8 (a-h, o-z)
         dimension a(0:500,-500:500),b(0:500,-500:500)
         dimension f(0:500,0:500),fx(0:500,0:500),fy(0:500,0:500)
         dimension fxx(0:500,0:500),fyy(0:500,0:500),fxy(0:500,0:500)
         dimension dens(0:500,0:500),dfs(0:500,0:500),dfss(0:500,0:500)

         
            pi=4.d0*datan(1.d0)


            open(unit=10,file='f1.dat')
            open(unit=11,file='derivatives1.dat')
            open(unit=12,file='skf1.dat')
            
          print 1000
 1000     format('input number of modes, number of points :',',' )
          read *, Nk,NN
          
          print 1001
 1001     format('input two random numbers :',',' )
          read *, xr1,yr1
          
          
c          xr1=6378836631.d0
c          yr1=77833324111.d0
          xr2=0.d0
          yr2=0.d0
          xr3=0.d0
          yr3=0.d0          
c---------random number generator-----------                    
          ar1=0.d0
          ar2=63308.d0
          ar3=-183326.d0
          br1=86098.d0
          br2=0.d0
          br3=-539608.d0                     
          do 1 k1=0,Nk
             do 2 k2=-Nk,Nk
               do 3 ng=1,2
                 do 4 nrr=1,2
          xr=ar1*xr1+ar2*xr2+ar3*xr3
          yr=br1*yr1+br2*yr2+br3*yr3

          xr=dmod(xr,2147483647.d0)
          yr=dmod(yr,2145483479.d0)
          
          xr3=xr2
          xr2=xr1
          xr1=xr

          yr3=yr2
          yr2=yr1
          yr1=yr
          
          rnd=dmod((xr-yr),2147483647.d0)/2147483647.d0
          rnd=dsqrt(rnd*rnd)
          
          if (nrr.eq.1) then
             rnd1=rnd
               else
             rnd2=rnd
           endif
 4              continue
              grnd =dsqrt(-2.d0*dlog(rnd1))*dcos(2.d0*pi*rnd2)

            if (ng.eq.1) then
               a(k1,k2)=grnd
                else
               b(k1,k2)=grnd
            endif
 3            continue
         
              print *, k1,k2,a(k1,k2),b(k1,k2)
         
 2           continue
 1         continue
c------------------------------------------------------

           a(0,0)=0.d0
           b(0,0)=0.d0

           h=2.d0*pi/dble(NN)
            s=0.d0
           do 5 j1=0,NN-1
              do 6 j2=0,NN-1
                 x=h*dble(j1)
                 y=h*dble(j2)
                 f(j1,j2)=0.d0
                 fx(j1,j2)=0.d0
                 fy(j1,j2)=0.d0
                 fxx(j1,j2)=0.d0
                 fyy(j1,j2)=0.d0
                 fxy(j1,j2)=0.d0
                 fxxx=0.d0
                 fyyy=0.d0
                 fxxy=0.d0
                 fxyy=0.d0
                 do 7 k1=0,Nk
                    do 8 k2=-Nk,Nk
                       arg=x*dble(k1)+y*dble(k2)
                       d1=dble(k1)
                       d2=dble(k2)
               f(j1,j2)=f(j1,j2)+
     *         a(k1,k2)*dcos(arg)+b(k1,k2)*dsin(arg)
               fx(j1,j2)=fx(j1,j2)-
     *         d1*a(k1,k2)*dsin(arg)+d1*b(k1,k2)*dcos(arg)
               fy(j1,j2)=fy(j1,j2)-
     *         d2*a(k1,k2)*dsin(arg)+d2*b(k1,k2)*dcos(arg)
               fxx(j1,j2)=fxx(j1,j2)-
     *         d1*d1*a(k1,k2)*dcos(arg)-d1*d1*b(k1,k2)*dsin(arg)
               fyy(j1,j2)=fyy(j1,j2)-
     *         d2*d2*a(k1,k2)*dcos(arg)-d2*d2*b(k1,k2)*dsin(arg)
               fxy(j1,j2)=fxy(j1,j2)-
     *         d1*d2*a(k1,k2)*dcos(arg)-d1*d2*b(k1,k2)*dsin(arg)
               f111=f111+
     *         d1*d1*d1*a(k1,k2)*dsin(arg)-d1*d1*d1*b(k1,k2)*dcos(arg)
               f222=f222+
     *         d2*d2*d2*a(k1,k2)*dsin(arg)-d2*d2*d2*b(k1,k2)*dcos(arg)
               f112=f112+
     *         d1*d1*d2*a(k1,k2)*dsin(arg)-d1*d1*d2*b(k1,k2)*dcos(arg)
               f122=f122+
     *         d1*d2*d2*a(k1,k2)*dsin(arg)-d1*d2*d2*b(k1,k2)*dcos(arg)
                              
 8                  continue
 7               continue
                 
                 f1=fx(j1,j2)
                 f2=fy(j1,j2)
                 f11=fxx(j1,j2)
                 f22=fyy(j1,j2)
                 f12=fxy(j1,j2)

                 dfs(j1,j2)=f1*f2*(f22-f11)+(f1*f1-f2*f2)*f12
                 
                 dfss(j1,j2)=f2*f2*(f11*f11+f1*f111+f12*f12+f2*f112)+
     *                f1*f1*(f22*f22+f2*f222+f12*f12+f1*f122)-
     *           2.d0*f1*f2*(f12*(f11+f22)+f1*f112+f2*f122)
                 
                 s=s+f(j1,j2)*f(j1,j2)/dble(NN)/dble(NN)

                 dens(j1,j2)=f11+f22
                 
 6            continue
 5          continue

             ss=0.d0
            do 21 j1=0,NN-1
               do 22 j2=0,NN-1
                 x=h*dble(j1)/2.d0/pi
                 y=h*dble(j2)/2.d0/pi
                  f(j1,j2)=f(j1,j2)/dsqrt(s)
                  ss=ss+f(j1,j2)*f(j1,j2)/dble(NN)/dble(NN)
                  write(10,*) x,y,f(j1,j2)
                  write(11,*) fx(j1,j2),fy(j1,j2),
     *                 fxx(j1,j2),fyy(j1,j2),fxy(j1,j2)
                  write(12,*) dens(j1,j2),dfs(j1,j2),dfss(j1,j2) 
 22            continue
 21          continue

             print *, ss

             close(10)
             close(11)
             close(12)

       
          stop
      end
