CKF  subroutine decisions(aerosol,gas,total,chtype)

      subroutine decisions(aerosol,gas,total,chtype,p,temp)
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'dropcom.inc'
      real*4 gas(ngas_aq), aerosol(nsect,naers), total, newso4
      integer chtype
            
c       WHAT DEGREE OF SIZE RESOLUTION?
      if (ires .eq. 1) then
        chtype = 1
        goto 150
      else if (ires .eq. 2) then
        goto 140
      else
        goto 80
      endif
     
c     Bulk calculations instead of size-resolved are executed for
c     those initial conditions in which the bulk model produces 
c     sulfate within 5% of the sized-resolved prediction at
c     t = 120 minutes
c     conditions for activation diam = 0.7, T = 283
c      
c     ad = alkaline dust.  This is alkaline dust concentration in ug/m3
c    
80    dust = 0
c      do isect = 5,10 ! for sections 0.1 to 10 um
      do isect = 4,10 ! for sections 0.04 to 40 um
      dust = dust + aerosol(isect,kcru)
      enddo
      ad = caratio*dust

      if (gas(ngh2o2) .ge. gas(ngso2)+1.0e-3) then
      chtype = 1
      goto 150
      else
      goto 81
      endif
      
81    if(gas(ngso2) .gt. 12e-3 .OR. gas(ngn) .gt. 15.e-3) then
      goto 140
      else
      goto 82
      endif
      
82    if (gas(ngn) .gt. gas(nga)) then
      goto 139
      else 
      goto 87
      endif 
      
87    if (gas(nga) .ge. 1.e-3+gas(ngn)) then
      if (gas(ngh2o2) .ge. 0.95*gas(ngso2)) then
      chtype = 1
      goto 150
      endif
      else
      goto 90
      endif
      
90    if (total .lt. 0.1) then
      goto 140
      else 
      goto 95
      endif

95    if (gas(nga) .le. 1e-3 .AND. ad .ge. 5.0) then 
      chtype = 1
      goto 150
      else if (gas(nga) .le. 5.e-3 .AND. ad .ge. 10.) then
      chtype = 1   
      goto 150
c new
      else if (gas(ngn).le.1.e-3.AND.gas(nga).ge.2.e-3+gas(ngn)) then
      if (gas(ngso2) .le. 7.e-3) then
      chtype = 1
      goto 150 
      else if (ad .ge. 2.) then
      chtype = 1
      goto 150 
      endif
       else if (gas(ngn).le.3.e-3.AND.gas(nga).ge.4.e-3+gas(ngn)) then
      chtype = 1
      goto 150
      else if (gas(ngn).le.7.e-3.AND.gas(nga).ge.3.e-3+gas(ngn)) then
      if (gas(ngso2) .le. 5.e-3) then
      chtype = 1
      goto 150 
      else if (ad .ge. 4. .AND. gas(ngso2) .le. 9.e-3) then
      chtype = 1
      goto 150 
      endif
      else if (ad .ge. 3 .AND. gas(nga) .le. 3.e-3) then
      if (gas(ngso2) .le. 4.e-3) then
      chtype = 1
      goto 150
      endif
      else if (ad .ge. 5 .AND. gas(ngso2) .le. 5.e-3) then
      if (gas(nga) .le. 7.e-3) then
      chtype = 1
      goto 150
      endif
c endnew        
      else if (gas(nga) .ge. 2.e-3+gas(ngn)) then
      if (gas(ngso2) .le. 5.e-3) then
      chtype = 1
      goto 150
      endif
      else if (gas(nga) .ge. 4.e-3+gas(ngn)) then
      if (gas(ngso2) .le. 10.e-3) then
      chtype = 1
      goto 150
      endif
      else if (ad .ge. 2.0 .AND. gas(nga) .le. 10.e-3) then
      if (gas(ngh2o2) .ge. gas(ngso2)) then
      chtype = 1
      goto 150
      endif
      else if (gas(nga) .le. 1.e-3 .AND. gas(ngso2) .ge. 3.e-3) then
      if (gas(ngh2o2) .ge. gas(ngso2)) then
      chtype = 1
      goto 150
      endif
      else
      goto 100
      endif

100	  if (total .ge. 0.3) then
          goto 110
	  else
	  goto 140
	  endif
	  
110       if (gas(nga) .ge. 5.e-3+gas(ngn)) then
          if (gas(ngso2) .le. 10.e-3) then
          chtype = 1
	  goto 150
	  endif
c new	  
 	  else if (gas(ngn).le.1.e-3.AND.gas(nga).ge.2.e-3+gas(ngn)) then
          chtype = 1
	  goto 150 
	  else if (gas(ngn).le.3.e-3.AND.gas(nga).ge.4.e-3+gas(ngn)) then
          chtype = 1
	  goto 150 
	  else if (gas(ngn).le.7.e-3.AND.gas(nga).ge.3.e-3+gas(ngn)) then
          chtype = 1
	  goto 150
	  else if (ad .ge. 3. .AND. gas(nga) .le. 10.e-3) then
	  if (gas(ngso2) .le. 5.e-3) then
	  chtype = 1
	  goto 150
	  endif
	  else if (ad .ge. 5. .AND. gas(nga) .le. 7.e-3) then
	  if (gas(ngso2) .le. 8.e-3) then
	  chtype = 1
	  goto 150
	  endif
c endnew	  
	  else if (gas(ngso2) .ge. 1.5e-3) then
	  if (gas(ngh2o2) .ge. gas(ngso2)) then
	  chtype = 1
	  goto 150 
	  endif
	  else if (gas(nga) .le. 12.e-3 .AND. ad .ge. 10) then
	  chtype = 1
	  goto 150
	  else if (gas(nga) .le. 1.e-3 .AND. ad .ge. 4.0) then
	  if (gas(ngso2) .le. 10.e-3) then
	  chtype = 1
	  goto 150 
	  endif
	  else if (gas(nga) .le. 5.e-3 .AND. ad .ge. 6.0) then
	  if (gas(ngso2) .le. 10.e-3) then
	  chtype = 1
	  goto 150
	  endif 
	  else if (gas(nga) .le. 7.0e-3 .AND. ad .ge. 8) then
	  if (gas(ngso2) .le. 10.e-3) then
	  chtype = 1
	  goto 150 
	  endif
	  else
	  goto 115 
	  endif	
	  
115	  if (total .ge. 0.5) then
          goto 120
	  else
	  goto 140
	  endif  
	  	  	  	  	  
120       if (gas(ngh2o2) .ge. 0.9*gas(ngso2)) then
	  chtype = 1 
	  goto 150
	  else if (gas(nga) .le. 1.e-3 .AND. ad .ge. 5.0) then
	  if (gas(ngso2) .le. 10.e-3) then
	  chtype = 1
	  goto 150
	  endif
	  else
	  goto 140
	  endif

139       if (total .ge. 0.1 .AND. gas(ngso2) .ge. 5.e-3) then
          if (gas(ngso2) .le. gas(ngh2o2)) then
	  chtype = 1
	  goto 150
	  endif
	  else if (total .ge. 0.3 .AND. gas(ngso2) .ge. 3.e-3) then
	  if (gas(ngso2) .le. gas(ngh2o2)) then
	  chtype = 1
	  goto 150
	  endif
c new
          else if (ad .ge. 5. .AND. gas(ngh2o2) .ge. gas(ngso2)) then
	  if (total .ge. 0.5) then 
	  chtype = 1
	  goto 150
	  else if (total .ge. 0.1 .AND. gas(ngn) .le. 2.e-3+gas(nga)) then
	  chtype = 1
	  goto 150
	  endif
c end new	  
	  else
	  goto 140
	  endif
	  
140       chtype = 2

CKF  At low SO2 concentrations, in order to avoid numerical difficulties,
CKF  transfer all SO2 to SO4=

150       if (gas(ngso2) .le. minso2) then

          chtype = 1
	  if (gas(ngh2o2) .ge. gas(ngso2)) then
	  gas(ngh2o2) = gas(ngh2o2) - gas(ngso2)
	  else
	  gas(ngh2o2) = 0.0
	  endif
	  fso2 = 1000*p*64./8.314e-2/temp
	  newso4 = gas(ngso2)*fso2*96./64.

	  gas(ngso2) = 0.0

	  do i=1,nsect
	  aerosol(i,na4) = aerosol(i,na4)+fdist(i)*newso4
	  enddo

	  endif
CKF
	  return
          end
