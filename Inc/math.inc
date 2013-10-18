c     include file for SVODE parameters and non-changing values
c     input to hybrid.f

      integer itol,iopt,mf,worki,lrw1,liw1,lrw2,liw2
      real*4 tola,tolr,workr
      integer numfunc, mode, nprint, maxfev, ml, mu, lr,ldfjac
      real factor, epsfcn, xtol 
      
c     FOR SVODE  
      parameter (itol = 1)                 ! Type of error control
      parameter (tola = 5.e-4)             ! Absolute tolerance
      parameter (tolr = 1.e-5)             ! Relative tolerance
      parameter (iopt = 1)
      parameter (mf = 22)   
      parameter (worki = 100000)           ! maximum # steps allowed
      parameter (workr = 300.0)            ! maximum time step
c  for bulk
      parameter (lrw1 = 22+9*meqn1+2*meqn1**2)
      parameter (liw1 = 30+meqn1) 
c  for sr
      parameter (lrw2 = 22+9*meqn2+2*meqn2**2)
      parameter (liw2 = 30+meqn2)                                                  
c
c   WHERE
c      itol: 4=use arrays for tolerances
c      tola: absolute tolerance in ug/m3
c      tolr: relative tolerance
c      itask: 1 for normal computation of output values of Y at t = TOUT.
c      istate: Integer flag (input and output).  Set ISTATE = 1.
c      iopt: 0 to indicate no optional input used.
c      rwork: Real work array of length at least..
c             20 + 16*NEQ                      for MF = 10,
c             22 +  9*NEQ + 2*NEQ**2           for MF = 21 or 22,
c             22 + 11*NEQ + (3*ML + 2*MU)*NEQ  for MF = 24 or 25.
c      lrw: Declared length of RWORK (in user's DIMENSION statement).
c      iwork: Integer work array of length at least..
c             30        for MF = 10,
c             30 + NEQ  for MF = 21, 22, 24, or 25.
c          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
c          and upper half-bandwidths ML,MU.
c      liw: Declared length of IWORK (in user's DIMENSION statement).
c      mf: Method flag.  Standard values are..
c          10 for nonstiff (Adams) method, no Jacobian used.
c          21 for stiff (BDF) method, user-supplied full Jacobian.
c          22 for stiff method, internally generated full Jacobian.
c          24 for stiff method, user-supplied banded Jacobian.
c          25 for stiff method, internally generated banded Jacobian.
c      iopt: 1 = some optional parameters used 
c           HERE:  workr: rwork(6) (max absolute step size allowed -
c                                    default value is infinite.)
c                  worki: iwork(6) (Maximum number of (internally defined) 
c                                    steps allowed during one call to the 
c				    solver. The default value is 500.)

c      FOR HYBRID.F

ckf       parameter (numfunc = 7)  
       parameter (numfunc = 6)
       parameter (xtol = 0.1e0**3)
       parameter (maxfev = 300*(numfunc+1) )
       parameter (ml = numfunc - 1, mu = numfunc -1)
       parameter (epsfcn = 0.0e0, factor = 100., mode = 2)
       parameter (nprint = 0)
       parameter (lr = numfunc*(numfunc+1)/2, ldfjac = numfunc)  
c
c      numfunc : number of functions and variables 
c      xtol : termination occurs when the rel error  between two consecutive 
c             iterates is at most xtol
c      maxfev : termination occurs when the number of calls to fcn is at least maxfev
c      ml     : specifies the number of subdiagonals within the band of the 
c               jacobian matrix.  if the jacobian is not banded, set ml to at 
c               least n -1.
c      mu     : specifies the number of superdiagonals within the band of the 
c               jacobian matrix.  if the jacobian is not banded, set mu to at 
c               least n -1.
c      epsfcn : used in determining a suitable step length for the 
c               forward-difference approximation 
c      factor : used in determining the initial step bound
c      mode   : if 1, the variables will be scaled internally; if 2, the 
c               scaling is specified by the input diag.
c      nprint : input variable that enables controlled
c               printing of iterates if it is positive. in this case,
c               fcn is called with iflag = 0 at the beginning of the first
c               iteration and every nprint iterations thereafter and
c               immediately prior to return, with x and fvec available
c               for printing. if nprint is not positive, no special calls
c               of fcn with iflag = 0 are made. 
