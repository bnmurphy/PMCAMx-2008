      PROGRAM after_main
c#######################################################################
c#
c#   Afterwrap performs all post-processing tasks on a new set of output
c#        by calling several subroutines already written
c#
c#        Calls:  after_saprc     daily avg conc'ns used in maps
c#                after_di_saprc        diurnal profile for specified grid cell
c#
c#            Obsolete:
c#                cmp_IMPROVE        match up IMPROVE measurements with 
c#                                 daily avg predictions
c#                cmp_AIRS        match AIRS measurements
c#                cmp_PAQS        match PAQS diurnal measurements with
c#                                 Pittsburgh predictions
c#                Processing ambient comparisons has been moved to the 
c#                  uni_ambient program!!!
c#
c#######################################################################
c#
c#    Written:  Ben Murphy 05-28-08
c#
c#######################################################################

      implicit none

      integer ndate, nlay, ncol, nrow, mspc, nhr, nspc
      real    UOMOC, ROMOC

      parameter (ndate = 17)   ! Number of simulation days
      parameter (nlay  = 1)   ! Layers of interest for diurnal data
      parameter (ncol  = 97)  !Number of columns
      parameter (nrow  = 90)  !Number of rows
      parameter (mspc  = 115)
      parameter (nhr   = 24)
      parameter (nspc  = 459) !CB4=158, SAPRC=459, SAPRC with NT-SOA=525

c      parameter (UOMOC = 1.4)  ! Urban OM:OC ratio
      parameter (UOMOC = 1.4)  ! Urban OM:OC ratio
      parameter (ROMOC = 2.0)  ! Rural OM:OC ratio

      integer tasks(3)

      character*199 indir, outdir, indir1     ! Input/Output Directories
      character*199 fname, fname1
      character*199 indirtmp                 ! Dummy Directory
      character*199 scenario              ! Run title
      character*199 amb_in              ! Ambient data location
      character*3   cdate(ndate)      ! Julian date array
      character*2   clay(nlay)        ! Character array for layers
      character*7   cspc(mspc)


      real*4  mconc(ndate,31,ncol,nrow,mspc)

      integer ilay, idate, io
      integer i, j, k

c
c   Set Task List
c
c           /daily avg map data,
c             diurnal profiles for average species,
c             diurnal maps,
c      data tasks /1,1,1/

c
c   Set Dates
c
      data cdate /'193','194','195','196','197','198','199','200',
     &            '201','202','203','204','205','206','207','208',
     &            '209'/
c      data cdate /'193','194','195','196','197','198','199','200','201',
c     &                   '202','203','204','205','206','207'/
c      data cdate /'203','204','205','206','207','208','209'/
c      data cdate /'195','196','197','198','199','200','201',
c     &                  '202','203','204','205','206','207','208','209'/
c       data cdate /'357','358','359','360','361','362','363','364','365'/
c       data cdate /'172','173','174','175','176','177','178','179','180','181'/
c       data cdate /'001','002','003'/

c      data cdate /'274','275','276','277','278','279','280','281',
c     &                    '282','283','284','285','286','287','288','289',
c     &                    '290' /

c      data cdate /'001','002','003','004','005','006','007'/  !,'008','009','010',
c     &            '011','012','013','014','015','016','017','018','019','020',
c     &            '021','022','023','024','025','026','027','028','029','030'/

c      data cdate /'182','183','184','185','186','187','188','189','190',
c     &            '191','192','193','194','195','196','197','198','199',
c     &            '200','201','202','203','204','205','206','207','208',
c     &            '209','210','211','212'/


c      data clay /'01','02','03','04','05','06','07','08','09','10',
c     &                 '11','12','13','14'/
      data clay /'01'/

c
c   Set Output Species
c
c      data cspc /'NO','NO2','O3','PAN','PAN2','MPAN','PBZN','NPHE',
c     &           'RNO3','CRES','DCB2','DCB3','HNO4','BALD','HONO',
c     &           'XN','HCHO','CCHO','RCHO','BACL','PROD','DCB1',
c     &           'PHEN','ISOP','ISPD','MVK','METH','MGLY','GLY',
c     &           'TERP','BPIN','LIMO','MONO','SESQ','HNO3','HO2H',
c     &           'HC2H','CO2H','CO3H','RC2H','RC3H','ACET','MEK',
c     &           'MEOH','COOH','ROOH','CO','ETHE','ALK1','ALK2',
c     &           'ALK3','ALK4','ALK5','ARO1','ARO2','OLE1','OLE2',
c     &           'NXOY','SO2','SULF','NH3','CPO','COO','CBS','CAS',
c     &           'APO','AOO','ABS','AAS','POC','PEC','CRST','PH2O',
c     &           'PCL','NA','PNH4','PNO3','PSO4',
c     &           'AAS1','AAS2','AAS3','AAS4','ABS1','ABS2','ABS3','ABS4',
c     &           'TotNH3','TotNO3','SOA','POA','TotOM','PM25','NOx',
c     &           'OOA','HOA'/

c      data cspc /'HCHO','CCHO','RCHO','LIMO','ETHE'/

c      mspc = 54
c      data cspc /'SO2', 'NH3','O3','OLE1','OLE2','ARO1','ARO2','CPO', 'COO', 'CAS', 'CBS',
c     &           'CBS1','CBS2','CBS3','CBS4','CBS5','CAS1','CAS2','CAS3','CAS4','CAS5',
c     &           'CNS','ANS','ABS1','ABS2','ABS3','ABS4','ABS5','AAS1','AAS2','AAS3','AAS4','AAS5',
c     &           'APO', 'AOO', 'ABS', 'AAS', 'POC', 'PEC','CRST','PH2O',
c     &           'PCL', 'NA',  'PNH4','PNO3','PSO4',
c     &           'TotNH3','TotNO3',   'SOA', 'POA', 'TOTOM','TOTOC','PM25','NOx'/


c     New OA SPECIATION WITH 8 BINS AND NT-SOA;   mspc = 128, mspc reduced = 116
c      data cspc /'NO',  'NO2', 'O3', 'TOTOM', 'TOTOC', 
cc     &                 'ISOP','ISPD','MGLY','GLY' ,'RCHO','HCHO','CCHO','LIMO',
c     &           'TERP','BPIN','SESQ','HNO3', !'HO2H',
c     &           'CO',  'ALK1','ALK2', !'ETHE',
c     &           'ALK3','ALK4','ALK5','ARO1','ARO2','OLE1','OLE2',
c     &           'SO2', 'NH3', 'CPO', 'COO', 'CAS', 'CBS','CNS',
c     &           'CPO1','CPO2','CPO3','CPO4','CPO5','CPO6','CPO7','CPO8',
c     &           'COO1','COO2','COO3','COO4','COO5','COO6','COO7','COO8',
c     &                 'CBS1','CBS2','CBS3','CBS4','CBS5',
c     &                 'CAS1','CAS2','CAS3','CAS4','CAS5',
c     &           'CNS1','CNS2','CNS3','CNS4','CNS5','CNS6','CNS7','CNS8',
c     &           'APO', 'AOO', 'ABS', 'AAS', 'ANS', 'POC', 'PEC','CRST','PH2O',
c     &           'PCL', 'NA',  'PNH4','PNO3','PSO4',
c     &           'AAS1','AAS2','AAS3','AAS4','AAS5',
c     &           'ABS1','ABS2','ABS3','ABS4','ABS5',
c     &           'TotNH3','TotNO3',   'SOA', 'POA','PM25','NOx',
c     &           'OOA', 'HOA', 'APO1','APO2','APO3','APO4',
c     &                         'APO5','APO6','APO7','APO8',
c     &                 'AOO1','AOO2','AOO3','AOO4','AOO5','AOO6','AOO7','AOO8',
c     &                 'ANS1','ANS2','ANS3','ANS4','ANS5','ANS6','ANS7','ANS8'/

c        mspc = 115
      data cspc /'NO',  'NO2', 'O3',  'RNO3','HNO4','HONO',
     &                 'ISOP','ISPD','MGLY','GLY' ,'RCHO','HCHO','CCHO','LIMO',
     &           'TERP','BPIN','SESQ','HNO3','HO2H',
     &           'CO',  'ETHE','ALK1','ALK2',
     &           'ALK3','ALK4','ALK5','ARO1','ARO2','OLE1','OLE2',
     &           'SO2', 'NH3', 'CPO', 'COO', 'CAS', 'CBS',
     &           'CPO1','CPO2','CPO3','CPO4','CPO5','CPO6','CPO7','CPO8','CPO9','CPO10',
     &           'COO1','COO2','COO3','COO4','COO5','COO6','COO7','COO8','COO9','COO10',
     &           'CBS1','CBS2','CBS3','CBS4','CAS1','CAS2','CAS3','CAS4',
     &           'APO', 'AOO', 'ABS', 'AAS', 'POC', 'PEC','CRST','PH2O',
     &           'PCL', 'NA',  'PNH4','PNO3','PSO4',
     &           'AAS1','AAS2','AAS3','AAS4','ABS1','ABS2','ABS3','ABS4',
     &           'TotNH3','TotNO3',   'SOA', 'POA', 'TOTOM','PM25','NOx','TOTOC',
     &           'OOA', 'HOA', 'APO1','APO2','APO3','APO4','APO5','APO6','APO7','APO8','APO9','APO10',
     &           'AOO1','AOO2','AOO3','AOO4','AOO5','AOO6','AOO7','AOO8','AOO9','AOO10'/

c      data cspc /'APO1','APO2','APO3','APO4','APO5','APO6','APO7','APO8','APO9','APO10',
c     &                 'AOO1','AOO2','AOO3','AOO4','AOO5','AOO6','AOO7','AOO8','AOO9','AOO10'/

c      data cspc /'HNO3','PNO3','HNO3g','NH3','PNH4','NH3g'/

c
c        CBM-IV Species
c        mspc = 27
c      data cspc /'NO',  'NO2', 'O3', 'TOL', 'OLE', 'OLE2','CG1', 'CG2','CG3','CG4',
c     &                 'SO2', 'NH3', 'HCL','SOA1','SOA2','SOA3','SOA4','POC','PEC','HNO3',
c     &                 'PSO4','PNO3','PNH4','TotOM','SOA','PTOC1','PTOC0'/


c
c   Set User-Defined Task List
c
      print *,'Enter 1 to output daily averaged concentrations (daily maps): '
      read(*,*),tasks(1)
      print *,'Enter 1 to output timeseries at specific sites (sites): '
      read(*,*),tasks(2)
      print *,'Enter 1 to output diurnal averaged concentrations (diurnal maps): '
      read (*,*), tasks(3)


c
c   Input/Output Structure
c
      print *,'Enter the folder name (eg. SAPRC.July/): '
      read(*,*), fname
      print *,'Enter the scenario name (eg. semivol.fxEmisBC.060308): '
      read(*,*), scenario
      indir    = '/home/mday/Output/' !//fname(1:INDEX(fname,' ')-1)//'/'
      outdir   = '/home/mday/processed_output/' !//fname(1:INDEX(fname,' ')-1)//'/'
c      scenario = 'semivol.fxEmisBC.060308'

      amb_in = '/home/bnmurphy/Research/Ambient/'

c
c   Begin Looping Over Layers
c
      do ilay = 1,nlay

c
c Iterate Days (input files)
c
         do idate=1,ndate

c
c   Set Input File
c
            if (nlay.gt.1) then
                    fname1='/4rpos.baseE.'//cdate(idate)//'.'//
     &              scenario(1:INDEX(scenario,' ')-1)//'.avrg'//clay(ilay)
            else
                    fname1='/4rpos.baseE.'//cdate(idate)//'.'//
     &              scenario(1:INDEX(scenario,' ')-1)//'.avrg'//clay(ilay)
c     &              scenario(1:INDEX(scenario,' ')-1)//'.avrg'
            endif
            indir1 = indir(1:INDEX(indir,' ')-1)//
     &          scenario(1:INDEX(scenario,' ')-1)
            write(6,*)'INPUT FILE: ',indir1(1:INDEX(indir1,' ')-1),fname1

c   Clear Temporary Variables
c
        do io = 1,mspc
            do i = 1,ncol
                do j = 1,nrow
                    do k = 1,nhr
                        mconc(idate,k,i,j,io) = 0.0
                    end do
                end do
            end do
        end do

c
c   Calling After Subroutine
c
        call after_sub_SAPRC(indir1, fname1, idate, ndate, nhr, nlay, cspc, 
     &             mconc, mspc, nspc, UOMOC, ROMOC)

c
c  Done with Day -> Start Next One
c
        print *,'Done with Day ',cdate(idate),'\n\n'
      end do !Day -> idate iterates

c
c   Call Subroutines
c
      if (tasks(1).eq.1) then
        print *,'TALLYING DAILY AVERAGE DATA'
        call after_daily_map(outdir,scenario,nlay,clay,ilay,ndate,cdate,mspc,cspc,nspc,mconc)
      endif
      if (tasks(2).eq.1) then
        print *,'TALLYING DIURNAL PROFILES FOR AVERAGE SPECIES AT SPECIFIC SITES'
        call after_timeseries(outdir,scenario,nlay,clay,ilay,ndate,cdate,mspc,cspc,nspc,mconc)
      end if
      if (tasks(3).eq.1) then
        print *,'Calculating Diurnal Profile for Domain Map'
        call after_diurnal_map(outdir,scenario,nlay,clay,ilay,ndate,cdate,mspc,cspc,nspc,mconc)
      end if

c    Processing comparisons to Ambient data has been moved to the uni_ambient program!!!

      end do   ! Iterate to next layer

      end
