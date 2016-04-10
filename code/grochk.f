

      integer function grochk(thrsinc,dayhrs,tave)
      use parm

      implicit none


!!...Argument declarations
      real  ::   tave, dayhrs
      logical :: thrsinc
!!...This function determines whether or not there is potential growth
!!...and if it is the first month of the growing season.

!!...grochk = 1:  Month of greenup in deciduous forests
!!...grochk = 0:  No greenup

!!...Notes:
!!...  * Assumes 30 daily time periods for deciduous greenup
!!...  * This routine assumes greenup will occur before day length
!!...     starts decreasing! -mdh 6/7/01

      integer :: numGreenUpPeriods, idf, j
      logical :: startd

      data     startd /.FALSE./
      data     numGreenUpPeriods /0/

      save     startd, numGreenUpPeriods
      

      grochk = 0
      j = ihru
      idf = idplt(j)
      

!!...If this is the first year for this block, reset STARTD.
      if ((curyr == 1) .and. (.not. startd)) then
        startd = .FALSE.
        numGreenUpPeriods = 0
      endif

!!...If it is spring and the temperature is high enough and
!!...you haven't already reapportioned carbon to leaves...
!!...Add number of hours in day to the check for when leaf out occurs, the
!!...value used here is extracted from the prephenology.c file from the BGC
!!...model code, cak - 06/10/02

      if ((thrsinc) .and. (.not. startd) .and.
     &    (tave .gt. TMPLFSf(idf)) .and. (dayhrs .gt. 10.917)) then
        grochk = 1
        startd = .TRUE.
        decidgrow = .TRUE.
      elseif ((.not. thrsinc) .or. (.not. decidgrow)) then
        startd = .FALSE.
        numGreenUpPeriods = 0
      endif
      
!!...New block - mdh 6/6/01
      if (startd) then
        numGreenUpPeriods = numGreenUpPeriods + 1
      endif
      if (numGreenUpPeriods.ge.1 .and. numGreenUpPeriods.le.30) then
        grochk = 1    
      endif             
     
      return
      end
