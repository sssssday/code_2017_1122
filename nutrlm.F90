

      subroutine nutrlm(elimit, tnelem, nparts, cfrac, eavail, maxec,   &
                       minec, maxeci, mineci, snfxmx, cprodl, eprodl,   &
                       plantNfix, afert, taufert, ctob)
     
      use parm

      implicit none


!!...Argument declarations
      integer tnelem, nparts
      real    cfrac(nparts), cprodl, ctob, elimit,      &
             eprodl(MAXIEL), eavail(MAXIEL), maxec(MAXIEL),             &
             minec(MAXIEL), maxeci(FPARTS-1,MAXIEL),                    &
             mineci(FPARTS-1,MAXIEL), plantNfix, snfxmx,                &
             afert(MAXIEL), taufert

!!...Nutrient limitation for plants is based on demand

!!...Local Variables
!!...NOTE:  Local variables cannot have adjustable array size.  ECFOR
!!...       is set to the largest array size which may occur.
      integer   iel, ipart,j,idf
      real      cpbe(MAXIEL), demand, ecfor(FPARTS-1,MAXIEL),           &
               maxNfix, totale, sum

!!...Definitions of Local variables
!!...  cpbe    - Carbon production limited by element
!!...  demand  - demand based on maximum E/C ratio
!!...  ecfor   - Actual E/C ratio by part 
!!...  maxNfix - N fixation
!!...  totale  - E available
!!...  ctob    - Compute weighted average carbon to biomass conversion factor

!!...Automatic fertilizer explanation
!!...  aufert 0.0 to 1.0
!!...    automatic fertilizer maybe applied to remove some nutrient
!!...    stress so that relyld is at least the value of aufert


      
      j = 0
      idf = 0

      
      j = ihru
      idf = idplt(j)
      
      if ((tnelem .le. 0) .or. (tnelem .gt. 3)) then
        write(*,*) 'Error in nutrlm, nelem out of bounds = ', tnelem
        write(*,*) 'Check <site>.100 file'
        STOP
      endif
      if (nparts .gt. FPARTS-1) then
        write(*,*) 'Error in nutrlm, nparts out of bounds = ', nparts
        STOP
      endif

!!...Compute production limitation

!!...P and S first so that potential N fixation accounts for P and S
!!...limitation on yield, AKM 18/8/2000
      if (tnelem .ge. 2) then

        do 10 iel =  tnelem, 2, -1
        
       

!! ...... DEMAND based on the maximum E/C ratio.
          demand = cprodl * maxec(iel)
          totale = eavail(iel)



!! ...... New E/C ratios by part based on E available.
          if (totale .gt. demand) then
            do 20 ipart = 1, nparts
              ecfor(ipart,iel) = maxeci(ipart,iel)
20          continue

          else
            if (demand .eq. 0.0) then
              write(*,*) 'Error in nutrlm, demand = 0.0'
              STOP
            endif
            do 30 ipart = 1, nparts
              ecfor(ipart,iel) = mineci(ipart,iel) +                &
               (maxeci(ipart,iel) - mineci(ipart,iel)) *            &
               totale / demand
30          continue
          endif

!! ...... Initialize local variables to zero
          cpbe(iel) = 0.0

!! ...... Total potential production with nutrient limitation
          do 40 ipart = 1, nparts
!! ........ Calculate the average nutrient content
            if (ipart .lt. 3) then
              cpbe(iel) = cpbe(iel) +                                   &
                         ((cfrac(ipart) * ecfor(ipart,iel)) / 2.5)          
            else
              cpbe(iel) = cpbe(iel) +                                   &
                         ((cfrac(ipart) * ecfor(ipart,iel)) / 2.0)
            endif
40        continue
!! ...... Calculate average E/C
          !cpbe(iel) was average nutrient content weighted by (biomass to carbon conversion factor)
          !ctob is weighted average (carbon to biomass conversion factor)
          cpbe(iel) = cpbe(iel) * ctob                          
          !cpbe(iel) is weighted avrage nutreint content as a fraction of C or average E/C                                                    
          
          
          if (cpbe(iel) .eq. 0.0) then
            write(*,*) 'Error in nutrlm, cpbe(iel) = 0.0'
            STOP
          endif
!! ...... Calculate potential production for the element with nutrient limitation
          cpbe(iel) = totale / cpbe(iel)
          
          
!! ...... Automatic fertilization, AKM 18/8/2000
          if ((taufert .gt. 0).and.(cpbe(iel)/cprodl .lt. taufert)) then
            cpbe(iel) = cprodl * taufert
!! ........ Mathcad equations
            sum = 0.0
            do ipart = 1, nparts
              sum = sum + cfrac(ipart) *                                    &
                   (maxeci(ipart,iel)-mineci(ipart,iel))
            end do
            afert(iel) = max((-taufert * cprodl * demand * minec(iel) /         &
                            (taufert * cprodl * sum - demand)) -                &
                             totale, 0.0)
          endif
10      continue
      endif

!!...Do the N

!!...Initialize fixation to 0
      maxNfix = 0.0
!!...DEMAND based on the maximum E/C ratio.
      demand = cprodl * maxec(N)
!!...N FIXATION
      maxNfix = fixn
      totale = eavail(N) + maxNfix


!!...New E/C ratios by part based on E available.
      if (totale .gt. demand) then
        do 120 ipart = 1, nparts
          ecfor(ipart,N) = maxeci(ipart,N)
120     continue

      else
        if (demand .eq. 0.0) then
          write(*,*) 'Error in nutrlm, demand = 0.0'
          STOP
        endif
        do 130 ipart = 1, nparts
          ecfor(ipart,N) = mineci(ipart,N) +            &
           (maxeci(ipart,N) - mineci(ipart,N)) *        &
           (totale / demand)
130     continue
      endif

!!...Initialize local variables to zero
      cpbe(N) = 0.0

!!...Total potential production with nutrient limitation
      do 140 ipart = 1, nparts
!! ..... Calculate the average nutrient content
        if (ipart .lt. 3) then
          cpbe(N) = cpbe(N) + ((cfrac(ipart) * ecfor(ipart,N)) / 2.5)
        else
          cpbe(N) = cpbe(N) + ((cfrac(ipart) * ecfor(ipart,N)) / 2.0)
        endif
140   continue
!!...Calculate average E/C
      cpbe(N) = cpbe(N) * ctob
      if (cpbe(N) .eq. 0.0) then
        write(*,*) 'Error in nutrlm, cpbe(N) = 0.0'
        STOP
      endif
!!...Calculate potential production for the element with nutrient limitation
      cpbe(N) = totale / cpbe(N)

!!...Put automatic fertilization here when necessary
!!...Automatic fertilization, AKM 18/8/2000
      if ((taufert .gt. 0) .and. (cpbe(N)/cprodl .lt. taufert)) then
        cpbe(N) = cprodl * taufert
!! ..... Mathcad equations
        sum = 0.0
        do ipart = 1, nparts
          sum = sum + cfrac(ipart) *                                &
               (maxeci(ipart,N) - mineci(ipart,N))
        end do
        afert(N) = max((-taufert * cprodl * demand * minec(N) /     &
                      (taufert * cprodl * sum - demand)) -          &
                       totale, 0.0)
      endif

!!...Compute the limiting element
      elimit = 0.0
      do 150 iel = 1, tnelem
        if (cprodl .gt. cpbe(iel)) then 
          cprodl = cpbe(iel)
          elimit = real(iel)
        endif
150   continue



!!...Recompute EPRODL
      do 170 iel = 1, tnelem
        eprodl(iel) = 0.0

!! ..... Total potential production with nutrient limitation
        do 160 ipart = 1, nparts
          eup(ipart,iel,j) = cprodl * cfrac(ipart) * ecfor(ipart,iel)
          eprodl(iel) = eprodl(iel) + eup(ipart,iel,j)
160     continue
170   continue

!!...Check to make sure the total flow won't exceed what's available,
!!...mdh - 10/09/02
      if (eprodl(N) - (eavail(N) + maxNfix) .gt. 0.0) then
        if (eprodl(N) - (eavail(N) + maxNfix) .gt. 0.05) then
!!          call message(' Warning: nutrlm - ')
!!          call message('   total N flow exceeds available N')
!!          write(*,*) 'eprodl(N)-(eavail(N)+maxNfix) = ',
!!     &                eprodl(N) - (eavail(N) + maxNfix)
        endif
!! ..... Prevent precision error - mdh 10/8/02
        eprodl(N) = eavail(N) + maxNfix
      endif
!!...Compute N fixation which actually occurs
      plantNfix = max(eprodl(N) - eavail(N), 0.0)

      return
      end
