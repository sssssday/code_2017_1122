
      subroutine froota (a2drat,h2ogef)
      use parm

      implicit none


!!... Argument declarations
     
      real    a2drat(3), h2ogef 

!!... This function determines the fraction of production going to fine
!!... roots in crops, grasses, and woody plants based based on water
!!... and nutrient availability.
!!...
!!... CALLED FROM:  cropDynC
!!...               treeDynC
!!...
!!... a2drat(iel) - the ratio of available mineral to mineral demand by
!!...               the plant
!!... h2ogef      - the effect of water on root to shoot ratio,
!!...               computed in POTCRP or POTFOR
!!... h2oeff      - effect of water stress on fraction of root carbon
!!... ntreff      - effect of nutrient stress in fraction of root carbon
!!... rtsh        - root to shoot ratio
!!...
!!... For trees only:
!!...   tfrtcn(1) - maximum fraction of C allocated to fine roots under
!!...               maximum nutrient stress
!!...   tfrtcn(2) - minimum fraction of C allocated to fine roots with no
!!...               nutrient stress
!!...   tfrtcw(1) - maximum fraction of C allocated to fine roots under
!!...               maximum water stress
!!...   tfrtcw(2) - minimum fraction of C allocated to fine roots with no
!!...               water stress
!!...
!!... For grasses/crops only (crop.100):
!!...   frtcindx  - (0) Use Great Plains eqn,
!!...               (1) perennial plant,
!!...               (2) annual plant, 
!!...               (3) perennial, use growing degree day implementation
!!...               (4) non-grain filling annual plant, growing degree day
!!...                   implementation, dynamic carbon allocation
!!...               (5) grain filling annual plant, growing degree day
!!...                   implementation, dynamic carbon allocation
!!...               (6) grain filling annual plant that requires a
!!...                   vernalization period (i.e. winter wheat), growing
!!...                   degree day implementation, dynamic carbon allocation
!!...   frtc(1)   - fraction of C allocated to roots at planting, with no
!!...               water or nutrient stress, used when FRTCINDX = 2, 4, 5, or
!!...               6
!!...   frtc(2)   - fraction of C allocated to roots at time FRTC(3), with no
!!...               water or nutrient stress, used when FRTCINDX = 2, 4, 5, or
!!...               6
!!...   frtc(3)   - time after planting (days with soil temperature greater
!!...               than RTDTMP) at which the FRTC(2) value is reached, used
!!...               when FRTCINDX = 2, 4, 5, or 6
!!...   frtc(4)   - maximum increase in the fraction of C going to the roots
!!...               due to water stress, used when FRTCINDX = 2, 4, 5, or 6
!!...   frtc(5)   - maximum increase in the fraction of C going to the roots
!!...               due to nutrient stress, used when FRTCINDX = 2, 4, 5, or 6
!!...   cfrtcn(1) - maximum fraction of C allocated to roots under maximum
!!...               water stress, used when FRTCINDX = 1 or 3
!!...   cfrtcn(2) - minimum fraction of C allocated to roots with no water
!!...               stress, used when FRTCINDX = 1 or 3
!!...   cfrtcw(1) - maximum fraction of C allocated to roots under maximum
!!...               nutrient stress, used when FRTCINDX = 1 or 3
!!...   cfrtcw(2) - minimum fraction of C allocated to roots with no nutrient
!!...               stress, used when FRTCINDX = 1 or 3

!!... Function declarations
      real line, ramp
      external line, ramp

!!... Local Variables
      integer :: iel, j, idf
      real    temp
      real    h2oeff, ntreff
      real    rtsh
      !!integer :: nelem
      !!real :: co2crs
      !!real :: ffroota
      
      j = 0
      idf = 0
          
      j = ihru
      idf = idplt(j)
      
     !! nelem = 1 !! only consider N limitation. change later. 2, means consider both N and P, 3 means N,P,S


      do 5 iel = 1, nelem
        if ((a2drat(iel) .lt. 0.0) .or. (a2drat(iel) .gt. 1.0)) then
          write(*,*) 'Error in froota, a2drat value out of bounds'
          STOP
        endif
5     continue

    

     
!!..... Effect of water limitation - inverse of growth curve
!!..... Compute water stress limitation on trees the same way
!!..... as we are computing water stress limitation for crops
!!..... and grasses cak - 09/12/03
c        h2oeff = frfrac(1) + ((frfrac(2) - frfrac(1)) * (1.0 - h2ogef))
        h2oeff = line(h2ogef, 0.0, TFRTCWf(j,1), 1.0, TFRTCWf(j,2))
!!..... Effect of nutrient limitation
        ntreff = 0.0
        do 10 iel = 1, nelem
          temp = line(a2drat(iel), 0.0, TFRTCNf(j,1), 1.0, TFRTCNf(j,2))
          ntreff = max(temp, ntreff) 
10      continue
!!..... Compute fraction of C to go to fine roots
        ffroota = max(h2oeff, ntreff)
        ffroota = min(ffroota, 0.99)
!!..... Change fraction of C to go to fine roots by effect of CO2
        ffroota = ffroota * co2crs



        
      if (ffroota .lt. 0.0) then
          ffroota = 0.0
      endif
      if (ffroota .gt. 1.0) then
      ffroota = 1.0
      endif 
      print*, "haha"
      
      return 
      end
