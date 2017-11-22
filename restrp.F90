

      subroutine restrp(elimit, tnelem, cerat, nparts, cfrac,      &
                       potenc, rimpct, storage, snfxmx, cprodl,             &
                       eprodl, uptake,  plantNfix, relyld)
     
       use parm

      implicit none

!!... Argument declarations
!!    need to check how to add one dimension to these arrang              
      integer   tnelem, nparts
      real    cerat(2,FPARTS,MAXIEL),                                              &
               cfrac(FPARTS), cprodl, elimit, eprodl(MAXIEL),                   &
               plantNfix, potenc, relyld, rimpct, snfxmx,                   &
               storage(1:MAXIEL), uptake(4,MAXIEL)
     
     

!!... Restrict the actual production based on C/E ratios.  Calculate
!!... minimum, and maximum whole plant nutrient concentrations.

!!... Local variables
!!... NOTE:  Local variables cannot have adjustable array size.  MINECI
!!...        and MAXECI are set to the largest array size which may occur.
      integer   iel, ipart, j, idf
      real      afert(MAXIEL), ctob, eavail(MAXIEL), maxec(MAXIEL),         &
               maxeci(FPARTS-1,MAXIEL), minec(MAXIEL),                      &
               mineci(FPARTS-1,MAXIEL), ustorg(MAXIEL)

     
      real, dimension (2) :: favail 
     
      j = 0
      idf = 0
      
 

      
      j = ihru
      idf = idplt(j)
      
      favail = 0.0
      iel = 1
     

!!... Definitions of Local variables
!!...   ctob - weighted average carbon to biomass conversion factor

      if ((tnelem .le. 0) .or. (tnelem .gt. 3)) then
        write(*,*) 'Error in restrp, nelem out of bounds = ', tnelem
        write(*,*) 'Check <site>.100 file'
        STOP
      endif
      if (nparts .gt. FPARTS-1) then
        write(*,*) 'Error in restrp, nparts out of bounds = ', nparts
        STOP
      endif

!!... Reset variables to zero
      cprodl = 0.0
      do 30 iel = 1, tnelem
        eprodl(iel) = 0.0
        afert(iel) = 0.0
        do 20 ipart = 1, nparts
          eup(ipart,iel,j) = 0.0
20      continue
30    continue

!!... There is no production if one of the mineral elements is not
!!... available.
      do 40 iel = 1, tnelem
!!        if ((availm(iel) .le. 1E-10) .and. (snfxmx .eq. 0.0) .and.
        if ((availm(j,iel) .le. 1E-4) .and.                                    &
            (SNFXMX2f(idf) .eq. 0.0) .and.                                      &!N FIXATION:maxNfix = snfxmx * cprodl
           (aufert .eq. 0.0)) then
          goto 999
        endif
40    continue

!!... Initialize cprodl
      cprodl = potenc
      
      favail(1) = SFAVAIL2f(idf) !!Calculate soil available nutrients, based on a maximum fraction (favail)
        
      !! calculation of favail(2) based on the first soil layer        
      favail(2) = max(0.2,                                                  &
                      min(0.2 + (sol_no3(1,j)/10. + sol_nh3(1,j)/10. )*     &
                   (0.4 - 0.2) / 2.0, 0.4))

!!... Calculate soil available nutrients, based on a maximum fraction
!!... (favail) and the impact of root biomass (rimpct), adding storage.
      do 45 iel = 1, tnelem   !! check what if there are three elimit here?
        eavail(iel) = (availm(j,iel) * favail(iel) * rimpct) +             &
                      forstg(iel,j)
45    continue

!!... Compute weighted average carbon to biomass conversion factor
      ctob = 0.0
      do 70 ipart = 1, nparts
        if (ipart .lt. 3) then
          ctob = ctob + (cfrac(ipart) * 2.5)
        else
          ctob = ctob + (cfrac(ipart) * 2.0)
        endif
70    continue
      !print *, "aaa"
!!... Calculate average E/C of whole plant (crop, grass, or tree)
      do 50 iel = 1, tnelem
        minec(iel) = 0.0
        maxec(iel) = 0.0
        do 60 ipart = 1, nparts
          if (cerat(IMAX,ipart,iel) .eq. 0.0) then
            write(*,*) 'Error is restrp, cerat(IMAX,ipart,iel) = 0.0'
             
            STOP
          endif
          if (cerat(IMIN,ipart,iel) .eq. 0.0) then
            write(*,*) 'Error is restrp, cerat(IMIN,ipart,iel) = 0.0'
            STOP
          endif
          mineci(ipart,iel) = 1.0 / cerat(IMAX,ipart,iel)
          maxeci(ipart,iel) = 1.0 / cerat(IMIN,ipart,iel)
          minec(iel) = minec(iel) + cfrac(ipart) * mineci(ipart,iel)
60      continue
!!..... Calculate average nutrient content based on roots and shoots only,
!!..... cak - 07/25/02
        maxec(iel) = maxec(iel) +                   &
                    ((cfrac(FROOT) * maxeci(FROOT,iel)) / 2.5)
        maxec(iel) = maxec(iel) +                   &
                    (((1.0 - cfrac(FROOT)) * maxeci(LEAF,iel)) / 2.5)
!!..... Calculate average E/C
        maxec(iel) = maxec(iel) * ctob
50    continue

!!... Compute the limitation
      call nutrlm(elimit, tnelem, nparts, cfrac, eavail, maxec, minec,          &
                 maxeci, mineci, snfxmx, cprodl, eprodl,               &
                 plantNfix, afert, aufert, ctob)

!!... Calculate relative yield - reimplemented by AKM 18/8/2000
      relyld = cprodl/potenc

!!... Calculate uptakes from all sources: storage, soil, plantNfix, and
!!... automatic fertilizer
      do 200 iel = 1, tnelem
        ustorg(iel) = min(storage(iel), eprodl(iel))
!!..... If storage pool contains all needed for uptake
        if (eprodl(iel) .le. ustorg(iel)) then
          uptake(ESTOR,iel) = eprodl(iel)
          uptake(ESOIL,iel) = 0.0
!!..... Otherwise, extra necessary from the soil pool
        elseif (eprodl(iel) .gt. ustorg(iel)) then
          uptake(ESTOR,iel) = storage(iel)
!!          uptake(ESOIL,iel) = eprodl(iel) - storage(iel)
!!....... subtract plantNfix -mdh 3/8/99
!!....... soil N uptake should not include monthly symbiotic N
!!....... fixation amount
          if (iel .eq. N) then
!!            uptake(ESOIL,iel)=min(eprodl(iel)-storage(iel)-plantNfix,
!!     &                            eavail(iel)-storage(iel)-plantNfix)
!!          else
!!            uptake(ESOIL,iel)=min(eprodl(iel)-storage(iel),
!!     &                            eavail(iel)-storage(iel))
            uptake(ESOIL,iel)=eprodl(iel)-storage(iel)-plantNfix
          else
            uptake(ESOIL,iel)=eprodl(iel)-storage(iel)
          endif
        endif
!!..... If (eprodl .gt. (uptake(ESTOR,iel) + uptake(ESOIL,iel) + plantNfix)
!!..... then pull from afert; afert should be > 0 only in this case, so
!!..... don't actually test eprodl, AKM 18/8/2000
        uptake(EFERT,iel) = afert(iel)
200   continue
!!... N fixation uptake was computed in the limitation routines
      uptake(ENFIX,N) = plantNfix

999   continue

      return
      end
