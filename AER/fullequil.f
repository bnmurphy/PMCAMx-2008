       subroutine fullequil(acon,aspres,acmet,aakeq,aakhen,awv,
     &  atemp,axsol)
       implicit real*8(a-h,o-z)
       real*4 acon(28),aspres(21),acmet(4),aakeq(17),aakhen(21)
       real*4 awv,atemp,axsol
       real*8 con(28), spres(21), cmet(4), akeq(17), akhen(21)
       real*8 wv,temp,xsol
c
c      CHANGE VARIABLES TO DOUBLE PRECISION TO AVOID LOW PH ERRORS
c
       do k=1,28
       con(k)=acon(k)
       enddo
c
       do k=1,21
       spres(k)=aspres(k)
       akhen(k)=aakhen(k)
       enddo
c
       do k=1,4
       cmet(k)=acmet(k)
       enddo
c
       do k=1,17
       akeq(k)=aakeq(k)
       enddo
c
       wv=awv
       temp=atemp
c
c      WE FIND THE INITIAL INTERVAL [aa,bb] FOR THE BISECTION METHOD
c      NEW VERSION (31/10/87)
       x=10.0d0**(-14)
       call electro(x,fa,con,spres,cmet,akeq,akhen,wv,temp)
       aa=x
       do 1035 i=-14,1
       x=10.0d0**i
       call electro(x,f,con,spres,cmet,akeq,akhen,wv,temp)
       if (f*fa .ge. 0.0d0) then
       aa=x
       fa=f
       else
       bb=x
       go to 1040
       end if
1035   continue
1037   format(' MISTAKE IN FULLEQUIL')
1040   continue
c
c      BISECTION METHOD FOR THE SOLUTION OF THE EQUATION
c      rtol : RELATIVE TOLERANCE
       ntry=0
       rtol=0.00001d0
1050   error= dabs(bb-aa)/aa
       ntry=ntry+1
       if (error .le. rtol) then
       xsol=(aa+bb)/2.0d0
       axsol=xsol                        ! single precision
       return
       end if
       xm=(aa+bb)/2.0d0
       call electro(xm,fm,con,spres,cmet,akeq,akhen,wv,temp)
       if (fa*fm .gt.  0.0d0) then
       aa=xm
       fa=fm
       else
       bb=xm
c       fb=fm
       end if
       go to 1050
       end
