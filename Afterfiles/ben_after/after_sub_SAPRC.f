      subroutine  after_sub_SAPRC(indir,filename,idate,ndate,nhr,nlay,cspc,mconc,
     &                            ncomp, nspc, UOMOC, ROMOC)
c######################################################################
c   BNM
c     This is a subroutine for processing .avrg data to be compared to 
c        ambient observations.  Program processes PMCAMx results
c        Values are NOT averaged over an entire day
c
c     Ben Murphy wrote this post-processor 5-7-08
c                 modified 9-24-08 to work with universal ambient code
c       outputs Gas, PM 2.5, and NOT(PM 10) for each hour
c
c######################################################################

c
c Define constants and arrays
c

      parameter(ncol = 97)
      parameter(nrow = 90)
      parameter(ifin  = 10)
      parameter(ifout = 30)

      integer             ndate, idate, nlay, nspc
      real*4         sconh(nhr,ncol,nrow,ncomp)         !Conc for output (hourly)
      real*4         mconc(ndate,31,ncol,nrow,ncomp)    !Conc for output (hourly)
      character*7    cspc(ncomp)                        !Output species name vector
      character*198  filename
      character*198  indir
      real                UOMOC, ROMOC
      
      parameter(gasfac = 1000)

      character*198 fname
      character*4   name(10),note(60),mspec(10,nspc)
      character*10  spcname(nspc)    !Input species name string

      real*4 maxh(nspc)
      real*4 btime,etime,rdum,xorg,yorg,delx,dely
      real*4 tbtime
      real*4 conh(nspc,nhr,ncol,nrow,nlay)  !Conc read from input file
      real*4 date1
      
      integer*4 ione,nspec,ibdate,iedate,idum,izero,nx,ny,nz
      integer*4 itbdate,col,row
      integer*4 ihr,msec,icut,line
      integer*4 iday,ix,iy, icol, irow

c
c   Open Input File
c
        fname = indir(1:INDEX(indir,' ')-1)//filename(1:INDEX(filename,' ')-1)
        write(6,*)'INPUT FILE: ',fname
        open(unit=ifin, file=fname, form='UNFORMATTED', status='old')
     
c
c   Read Input Header
c
      read(ifin,End=990)name,note,ione,nspec,ibdate,btime,iedate,etime
      write(6,*)(name(i)(1:1),i=1,10)
      write(6,*)(note(i)(1:1),i=1,60)
      write(6,*),nspec,ibdate,btime,iedate,etime
      If(nspec .ne. nspc) then
          print *,'nspec(=',nspec,') is not equal to nspc(=',nspc,')'
          stop
      end if
      read(ifin,END=990)rdum,rdum,iutm,xorg,yorg,delx,dely,nx,ny,nz
     &                  ,idum,idum,rdum,rdum,rdum
      write(6,*)iutm,xorg,yorg,delx,dely
      write(6,*)nx,ny,nz
      read(ifin, END=990)izero,izero,nx,ny
      write(6,*)izero,nx,ny
      if(nx.ne.ncol)stop'nx is not equal to ncol'
      if(ny.ne.nrow)stop'ny is not equal to nrow'
c      if(nz.ne.nlay)stop'nz is not equal to nlay'

      read(ifin,END=990)((mspec(i,j),i=1,10),j=1,nspec)
      do isp = 1,nspc
          write (spcname(isp),*), (mspec(n,isp)(1:1),n=1,10)
      end do

c
c   Read Conc'ns into memory and close input file
c
      print *,'\nGetting Hourly Conc Data'
      do ihr = 1,nhr
          read(ifin)ibdate,btime,iedate,etime
          write(6,*)'ihr = ',ihr-1
          do isp = 1,nspec
            do ilay = 1,1
               read(ifin,END=990)ione,(mspec(n,isp),n=1,10),
     &                  ((conh(isp,ihr,i,j,ilay),i=1,nx),j=1,ny)
            end do  !ilay
          end do  !isp
      end do  !ihr

      close(ifin)  !Close input file

c
c   Initialize Output Species Matrix
c
          print *,'Initializing Output Concentration Arrays...'
          do m=1,ncomp
           do i = 1,nx
              do j = 1,ny
                 do k = 1,nhr
                    sconh(k,i,j,m)=0.0
                 end do
              end do
           end do
          end do

c
c   Sum up concentrations
c
      print *,'Assigning Internal Species to Output Species and Summing Aggregates...'
      print *,'\nSpecies to be tallied: ',(cspc(i),i=1,ncomp)
      do isp = 1,nspec
        
          isund = INDEX(spcname(isp),'_')
          if (isund.eq.0) then !Gas
c                  ***** Gases *****
            do kspc = 1,ncomp  !Loop through output species

              if (spcname(isp)(1:INDEX(spcname(isp),' ')-1).eq.
     &                  cspc(kspc)(1:INDEX(cspc(kspc),' ')-1)) then
                do i = 1,ncol
                  do j = 1,nrow
                    do khr = 1,nhr
                            sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc) 
     &                            + conh(isp,khr,i,j,1)*gasfac
                    end do
                  end do
                end do
              elseif (spcname(isp)(1:3).eq.cspc(kspc)(1:3).and.
     &               (spcname(isp)(1:3).eq.'CPO'.or.
     &                spcname(isp)(1:3).eq.'COO'.or.
     &                spcname(isp)(1:3).eq.'CBS'.or.
     &                spcname(isp)(1:3).eq.'CNS'.or.
     &                spcname(isp)(1:3).eq.'CAS')    ) then   !Sum up all volatility bins
                do i = 1,ncol
                  do j = 1,nrow
                    do khr = 1,nhr
                            sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc)
     &                            + conh(isp,khr,i,j,1)
                    end do
                  end do
                end do
              elseif (cspc(kspc).eq.'NOx'.and.
     &                (spcname(isp).eq.'NO'.or.spcname(isp).eq.'NO2')) then
                do i = 1,ncol
                  do j = 1,nrow
                    do khr = 1,nhr
                    sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc)
     &                        + conh(isp,khr,i,j,1)*gasfac
                    end do
                  end do
                end do
              elseif ((cspc(kspc).eq.'TotNH3'.or.cspc(kspc).eq.'NH3g')
     &                        .and.spcname(isp).eq.'NH3') then
                do i = 1,ncol
                  do j = 1,nrow
                    do khr = 1,nhr
                    sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc)
     &                        + conh(isp,khr,i,j,1)*(1.0133*100000*17)/(8.314*300)  !ug/m3
                    end do
                  end do
                end do
              elseif ((cspc(kspc).eq.'TotNO3'.or.cspc(kspc).eq.'HNO3g')
     &                        .and.spcname(isp).eq.'HNO3') then
                do i = 1,ncol
                  do j = 1,nrow
                    do khr = 1,nhr
                            sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc)
     &                            + conh(isp,khr,i,j,1)*(1.0133*100000*63)/(8.314*300)  !ug/m3
                    end do
                  end do
                end do
              end if
            end do
              
           else   !Aerosol

c                 ***** Aerosols *****
c                Doesn't do PM 10 yet!!!
            read (spcname(isp)(isund+1:isund+3),*), bin
            if (bin.le.6) then  !PM 2.5
               do kspc = 1,ncomp !Loop through output species
                   isblk = INDEX(cspc(kspc),' ')-1
 
c                ** Add up all of the non-aggregate species **
                       if (spcname(isp)(1:4).eq.cspc(kspc)(1:4)
     &                        .and.isblk.eq.4 ) then
                     do i = 1,ncol
                        do j = 1,nrow
                           do khr = 1,nhr
                                   sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc)
     &                                  + conh(isp,khr,i,j,1)
                           end do
                        end do
                     end do
                  elseif (spcname(isp)(1:3).eq.cspc(kspc)(1:3)
     &                            .and.isblk.eq.3) then
                     do i = 1,ncol
                        do j = 1,nrow
                           do khr = 1,nhr
                                   sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc)
     &                                  + conh(isp,khr,i,j,1)
                           end do
                        end do
                     end do
                  elseif (spcname(isp)(1:2).eq.'NA'.and.cspc(kspc)(1:2).eq.'NA') then
                     do i = 1,ncol
                        do j = 1,nrow
                           do khr = 1,nhr
                                   sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc)
     &                                  + conh(isp,khr,i,j,1)
                           end do
                        end do
                     end do
                  end if

c                ** Add up aggregate species
                  if (cspc(kspc).eq.'TotNH3' .and.
     &                      spcname(isp)(1:4).eq.'PNH4') then
                     do i = 1,ncol
                        do j = 1,nrow
                          do khr = 1,nhr
                                  sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc)
     &                                  + conh(isp,khr,i,j,1)
                          end do
                        end do
                     end do
                  elseif (cspc(kspc).eq.'TotNO3' .and.
     &                          spcname(isp)(1:4).eq.'PNO3') then
                     do i = 1,ncol
                        do j = 1,nrow
                          do khr = 1,nhr
                                  sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc)
     &                                    + conh(isp,khr,i,j,1)
                          end do
                        end do
                     end do
                  elseif ((cspc(kspc).eq.'TOTOC')
     &                        .and.(spcname(isp)(1:3).eq.'ABS'.or.
     &                              spcname(isp)(1:3).eq.'AOO'.or.
     &                              spcname(isp)(1:3).eq.'AAS'.or.
     &                              spcname(isp)(1:3).eq.'ANS'.or.
     &                              spcname(isp)(1:3).eq.'POC')    ) then
                     do i = 1,ncol
                        do j = 1,nrow
                          do khr = 1,nhr
                             sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc) + 
     &                                conh(isp,khr,i,j,1)  /ROMOC
                          end do
                        end do
                     end do
                  elseif ((cspc(kspc).eq.'TOTOC')
     &                        .and.(spcname(isp)(1:3).eq.'APO') ) then
                     do i = 1,ncol
                        do j = 1,nrow
                          do khr = 1,nhr
                           sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc) + 
     &                                conh(isp,khr,i,j,1)  /UOMOC
                          end do
                        end do
                     end do
                  elseif (cspc(kspc).eq.'SOA' .and.
     &                          (spcname(isp)(1:3).eq.'ABS'.or.
     &                           spcname(isp)(1:3).eq.'ANS'.or.
     &                           spcname(isp)(1:3).eq.'AAS')    ) then
                     do i = 1,ncol
                        do j = 1,nrow
                          do khr = 1,nhr
                           sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc) +
     &                                conh(isp,khr,i,j,1)
                          end do
                        end do
                     end do
                  elseif (cspc(kspc).eq.'POA' .and.
     &                          (spcname(isp)(1:3).eq.'APO'.or.
     &                           spcname(isp)(1:3).eq.'AOO')    ) then
                     do i = 1,ncol
                        do j = 1,nrow
                          do khr = 1,nhr
                           sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc) +
     &                                conh(isp,khr,i,j,1)
                          end do
                        end do
                     end do
                  elseif (cspc(kspc).eq.'OOA' .and.
     &                          (spcname(isp)(1:3).eq.'ABS'.or.
     &                           spcname(isp)(1:3).eq.'AAS'.or.
     &                           spcname(isp)(1:3).eq.'ANS'.or.
     &                           spcname(isp)(1:3).eq.'POC'.or.
     &                           spcname(isp)(1:3).eq.'AOO')    ) then
                     do i = 1,ncol
                        do j = 1,nrow
                          do khr = 1,nhr
                           sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc) +
     &                                conh(isp,khr,i,j,1)
                          end do
                        end do
                     end do
                  elseif (cspc(kspc).eq.'HOA' .and.
     &                          spcname(isp)(1:3).eq.'APO') then
                     do i = 1,ncol
                        do j = 1,nrow
                          do khr = 1,nhr
                           sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc) +
     &                                conh(isp,khr,i,j,1)
                          end do
                        end do
                     end do
                  elseif (cspc(kspc).eq.'TOTOM' .and.
     &                          (spcname(isp)(1:3).eq.'ABS'.or.
     &                           spcname(isp)(1:3).eq.'AAS'.or.
     &                           spcname(isp)(1:3).eq.'ANS'.or.
     &                           spcname(isp)(1:3).eq.'POC'.or.
     &                           spcname(isp)(1:3).eq.'APO'.or.
     &                           spcname(isp)(1:3).eq.'AOO')    ) then

          do i = 1,ncol
                        do j = 1,nrow
                          do khr = 1,nhr
                           sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc) +
     &                                conh(isp,khr,i,j,1)
                          end do
                        end do
                     end do
                  elseif (cspc(kspc).eq.'PM25' .and.
     &             (spcname(isp)(1:3).eq.'APO' .or.spcname(isp)(1:3).eq.'AOO' .or.
     &              spcname(isp)(1:3).eq.'POC' .or.spcname(isp)(1:3).eq.'AAS' .or.
     &              spcname(isp)(1:3).eq.'ABS' .or.spcname(isp)(1:3).eq.'PEC' .or.
     &              spcname(isp)(1:4).eq.'PSO4'.or.spcname(isp)(1:4).eq.'PNH4'.or.
     &              spcname(isp)(1:4).eq.'PNO3'.or.spcname(isp)(1:3).eq.'PCL' .or.
     &              spcname(isp)(1:4).eq.'CRST'.or.spcname(isp)(1:3).eq.'ANS' .or.
     &              spcname(isp)(1:2).eq.'NA') ) then
                     do i = 1,ncol
                        do j = 1,nrow
                          do khr = 1,nhr
                           sconh(khr,i,j,kspc) = sconh(khr,i,j,kspc) +
     &                                conh(isp,khr,i,j,1)
                          end do
                        end do
                     end do
                   
                  end if  ! end aggregates
               end do
            end if  !PM 2.5
                
            end if  !Gas or Aerosols
      end do            !isp

c
c   Assign sconh --> mconc
c
      do ihr = 1,nhr
          do icol = 1,ncol
            do irow = 1,nrow
                do ispec = 1,ncomp
                    mconc(idate,ihr,icol,irow,ispec) = sconh(ihr,icol,irow,ispec)
                enddo
            enddo
          enddo
      enddo
c      print *,'sconh: ',(sconh(i,65,51,3),i=1,nhr)
c      print *,'In subroutine: Date=',idate,'  Species=',cspc(3),'  mconc=',(mconc(idate,i,65,51,3),i=1,nhr)

c
c  Calculate Aggregate Quantities
c    Species Index MUST match up with vector <cspc>
c

c        do khr = 1,nhr
c        do i = 1,nx
c        do j = 1,ny
        
c        NH4 = PNH4 + NH3
c              sconh(32,khr,i,j) = sconh(29,khr,i,j) + sconh(14,khr,i,j)*(1.0133*10**2*17/(8.314*300))
c        NO3 = PNO3 + HNO3        
c        sconh(33,khr,i,j) = sconh(30,khr,i,j) + sconh(7,khr,i,j)*(1.0133*10**2*63/(8.314*300))
c        CIS
c          sconh(62,khr,i,j) = sconh(62,khr,i,j)*(1.0133*10**2*250/(8.314*300))
c        CR1
c          sconh(63,khr,i,j) = sconh(63,khr,i,j)*(1.0133*10**2*250/(8.314*300))
c        POA = APO + AOO + POC      
c          sconh(35,khr,i,j) = sconh(19,khr,i,j) + sconh(20,khr,i,j) + sconh(23,khr,i,j)
c        SOA = ABS + AAS        
c          sconh(34,khr,i,j) = sconh(21,khr,i,j) + sconh(22,khr,i,j)
c        OC = POA + SOA
c          sconh(36,khr,i,j) = sconh(35,khr,i,j) + sconh(34,khr,i,j)
c        PM2.5 = OC + PEC + CRST + PCL + NA + PH2O + PNH4 + PNO3 + PSO4        
c          sconh(37,khr,i,j) = sconh(36,khr,i,j) + sconh(24,khr,i,j) + sconh(25,khr,i,j) + 
c    &                              sconh(27,khr,i,j) + sconh(28,khr,i,j) + sconh(26,khr,i,j) + 
c     &                              sconh(29,khr,i,j) + sconh(30,khr,i,j) + sconh(31,khr,i,j)

c        end do
c        end do
c        end do

c
c  Done with Day -> Start Next One
c
       print *,'Done with Day \n\n'

 990   end
