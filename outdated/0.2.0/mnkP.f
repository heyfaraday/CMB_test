
         implicit real*8 (a-h, o-z)
         dimension f(0:500,0:500),fx(0:500,0:500),fy(0:500,0:500)
         dimension fxx(0:500,0:500),fyy(0:500,0:500),fxy(0:500,0:500)
         dimension xx(4),yy(4),xx1(4),yy1(4),xx2(4),yy2(4)
         dimension dens(0:500,0:500),dfs(0:500,0:500),dfss(0:500,0:500)

         
            pi=4.d0*datan(1.d0)


            open(unit=10,file='P.dat')
            open(unit=11,file='derivativesP.dat')
            open(unit=12,file='skfP.dat')
            
          NN = 254

          h=1.d0/dble(NN)

          fmax=-10000000000000.d0
          fmin=10000000000000.d0
c------------reading---------------------------          
          do 1 i=0,NN-1
             do 2 j=0,NN-1
                read(10,*) f(i,j),x,y
                read(11,*) fx(i,j),fy(i,j),
     *               fxx(i,j),fyy(i,j),fxy(i,j)
c     read(12,*) dens(i,j),dfs(i,j),dfss(i,j)
                 read(12,*) dens(i,j),dfs(i,j),dfss(i,j)
                if (f(i,j).ge.fmax) then
                   fmax=f(i,j)
                  else
                  endif
                  if (f(i,j).le.fmin) then
                   fmin=f(i,j)
                  else
                 endif 
 2           continue
 1        continue

             close(10)
             close(11)
             close(12)
             
             print *, fmin,fmax
c------------------------------------------------

c---------isocontour line------------------------

             open(unit=10,file='levP.dat')

              do 777 l=1,21

                 s=fmin+(fmax-fmin)/10.d0*dble(l)
                 
             dll=0.d0
             
             do 3 i=0,NN-1
                do 4 j=0,NN-1
              x=dble(i)*h
              y=dble(j)*h

               is=0
              
           if (f(i,j).le.s.and.f(i+1,j).gt.s.or.
     *         f(i,j).gt.s.and.f(i+1,j).le.s) then
              is=is+1
              xx(is)=x+(s-f(i,j))/(f(i+1,j)-f(i,j))*h
              yy(is)=y
              write(10,*) xx(is),yy(is)
                else
            endif
                
            if (f(i,j+1).le.s.and.f(i+1,j+1).gt.s.or.
     *           f(i,j+1).gt.s.and.f(i+1,j+1).le.s) then
              is=is+1
              xx(is)=x+(s-f(i,j+1))/(f(i+1,j+1)-f(i,j+1))*h
              yy(is)=y+h
               write(10,*) xx(is),yy(is)
                else
             endif

                
             if (f(i,j).le.s.and.f(i,j+1).gt.s.or.
     *           f(i,j).gt.s.and.f(i,j+1).le.s) then
              is=is+1
              yy(is)=y+(s-f(i,j))/(f(i,j+1)-f(i,j))*h
              xx(is)=x
               write(10,*) xx(is),yy(is)
                else
            endif
                
             if (f(i+1,j).le.s.and.f(i+1,j+1).gt.s.or.
     *           f(i+1,j).gt.s.and.f(i+1,j+1).le.s) then
              is=is+1
              yy(is)=y+(s-f(i+1,j))/(f(i+1,j+1)-f(i+1,j))*h
              xx(is)=x+h
               write(10,*) xx(is),yy(is)
                else
            endif

            if (is.eq.2) then
               dl=dsqrt((xx(2)-xx(1))*(xx(2)-xx(1))+
     *            (yy(2)-yy(1))*(yy(2)-yy(1)))
               dll=dll+dl
            else
               endif
            
 4              continue
 3           continue

               print *, s,dll
             
 777          continue
c----------------------------------------------------
              close(10)

c---------------extrema-----------------------------
              s=0.d0
              dll=0.d0
              
               next=0
               nmax=0
               nmin=0
               nsad=0

               open(unit=11,file='maxP.dat')
               open(unit=12,file='minP.dat')
               open(unit=13,file='sadP.dat')
               open(unit=14,file='skeletP.dat')

               ksk=0
               
               do 5 i=0,NN-1
                do 6 j=0,NN-1
              x=dble(i)*h
              y=dble(j)*h
              
c---------------dfs=0------------------------------
               is=0
              
           if (dfs(i,j).le.s.and.dfs(i+1,j).gt.s.or.
     *         dfs(i,j).gt.s.and.dfs(i+1,j).le.s) then
              is=is+1
              xx(is)=x+(s-dfs(i,j))/(dfs(i+1,j)-dfs(i,j))*h
              yy(is)=y
              write(14,*) xx(is),yy(is)
              ksk=ksk+1
                else
            endif
                
            if (dfs(i,j+1).le.s.and.dfs(i+1,j+1).gt.s.or.
     *          dfs(i,j+1).gt.s.and.dfs(i+1,j+1).le.s) then
              is=is+1
              xx(is)=x+(s-dfs(i,j+1))/(dfs(i+1,j+1)-dfs(i,j+1))*h
              yy(is)=y+h
              write(14,*) xx(is),yy(is)
              ksk=ksk+1
                else
             endif

                
             if (dfs(i,j).le.s.and.dfs(i,j+1).gt.s.or.
     *           dfs(i,j).gt.s.and.dfs(i,j+1).le.s) then
              is=is+1
              yy(is)=y+(s-dfs(i,j))/(dfs(i,j+1)-dfs(i,j))*h
              xx(is)=x
              write(14,*) xx(is),yy(is)
              ksk=ksk+1
                else
            endif
                
             if (dfs(i+1,j).le.s.and.dfs(i+1,j+1).gt.s.or.
     *           dfs(i+1,j).gt.s.and.dfs(i+1,j+1).le.s) then
              is=is+1
              yy(is)=y+(s-dfs(i+1,j))/(dfs(i+1,j+1)-dfs(i+1,j))*h
              xx(is)=x+h
              write(14,*) xx(is),yy(is)
              ksk=ksk+1
                else
            endif

            if (is.eq.2) then
               dl=dsqrt((xx(2)-xx(1))*(xx(2)-xx(1))+
     *            (yy(2)-yy(1))*(yy(2)-yy(1)))
               dll=dll+dl
            else
               endif
c-------------fx=0---------------------------------
               ix=0              
           if (fx(i,j).le.s.and.fx(i+1,j).gt.s.or.
     *         fx(i,j).gt.s.and.fx(i+1,j).le.s) then
              ix=ix+1
              xx1(ix)=x+(s-fx(i,j))/(fx(i+1,j)-fx(i,j))*h
              yy1(ix)=y
                else
            endif
                
            if (fx(i,j+1).le.s.and.fx(i+1,j+1).gt.s.or.
     *           fx(i,j+1).gt.s.and.fx(i+1,j+1).le.s) then
              ix=ix+1
              xx1(ix)=x+(s-fx(i,j+1))/(fx(i+1,j+1)-fx(i,j+1))*h
              yy1(ix)=y+h
                else
             endif
                
             if (fx(i,j).le.s.and.fx(i,j+1).gt.s.or.
     *           fx(i,j).gt.s.and.fx(i,j+1).le.s) then
              ix=ix+1
              yy1(ix)=y+(s-fx(i,j))/(fx(i,j+1)-fx(i,j))*h
              xx1(ix)=x
                else
            endif
                
             if (fx(i+1,j).le.s.and.fx(i+1,j+1).gt.s.or.
     *           fx(i+1,j).gt.s.and.fx(i+1,j+1).le.s) then
              ix=ix+1
              yy1(ix)=y+(s-fx(i+1,j))/(fx(i+1,j+1)-fx(i+1,j))*h
              xx1(ix)=x+h
                else
                endif

c-------------fy=0---------------------------------
                
                 iy=0              
           if (fy(i,j).le.s.and.fy(i+1,j).gt.s.or.
     *         fy(i,j).gt.s.and.fy(i+1,j).le.s) then
              iy=iy+1
              xx2(iy)=x+(s-fy(i,j))/(fy(i+1,j)-fy(i,j))*h
              yy2(iy)=y
                else
            endif
                
            if (fy(i,j+1).le.s.and.fy(i+1,j+1).gt.s.or.
     *           fy(i,j+1).gt.s.and.fy(i+1,j+1).le.s) then
              iy=iy+1
              xx2(iy)=x+(s-fy(i,j+1))/(fy(i+1,j+1)-fy(i,j+1))*h
              yy2(iy)=y+h
                else
             endif
                
             if (fy(i,j).le.s.and.fy(i,j+1).gt.s.or.
     *           fy(i,j).gt.s.and.fy(i,j+1).le.s) then
              iy=iy+1
              yy2(iy)=y+(s-fy(i,j))/(fy(i,j+1)-fy(i,j))*h
              xx2(iy)=x
                else
            endif
                
             if (fy(i+1,j).le.s.and.fy(i+1,j+1).gt.s.or.
     *           fy(i+1,j).gt.s.and.fy(i+1,j+1).le.s) then
              iy=iy+1
              yy2(iy)=y+(s-fy(i+1,j))/(fy(i+1,j+1)-fy(i+1,j))*h
              xx2(iy)=x+h
                else
                endif
              
c---------------------------------------------------------        

            if (ix.eq.2.and.iy.eq.2) then
               
            a1=(yy1(1)-yy1(2))/(xx1(1)-xx1(2))
            b1=yy1(1)-a1*xx1(1)
            a2=(yy2(1)-yy2(2))/(xx2(1)-xx2(2))
            b2=yy2(1)-a2*xx2(1)

            xe=(b2-b1)/(a1-a2)
            ye=a1*xe+b1

            if (xe.ge.x.and.xe.lt.x+h.and.ye.ge.y.and.ye.lt.y+h) then
               next=next+1
               dd=(fxx(i,j)-fyy(i,j))*(fxx(i,j)-fyy(i,j))
     *             +4.d0*fxy(i,j)*fxy(i,j)
               dl1=-fxx(i,j)-fyy(i,j)+dsqrt(dd)
               dl2=-fxx(i,j)-fyy(i,j)-dsqrt(dd)
               if (dl1.le.0.and.dl2.le.0) then
                  nmax=nmax+1
                 write(11,*) xe,ye,f(i, j)
                    else
                    endif
               if (dl1.ge.0.and.dl2.ge.0) then
                  nmin=nmin+1
                  write(12,*) xe,ye,f(i, j)
                    else
                    endif
                    if (dl1.ge.0.and.dl2.le.0.or.
     *                  dl1.le.0.and.dl2.ge.0) then
                   nsad=nsad+1
                  write(13,*) xe,ye 
                    else
                  endif
                    
               
            else
               endif
            
                  else
            endif
            
 6             continue
 5            continue
c----------------------------------------------------

              print *, next,nmax,nmin,nsad
              print *, ksk
              
              close(11)
              close(12)
              close(13)
              close(14)





              
          stop
      end
