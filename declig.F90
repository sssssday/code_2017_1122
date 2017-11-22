
     
      subroutine declig(flag, aminrl,ligcon,lyr,tnelem,ps1co2,rnew,rsplig,tcflow,tcstva)

      use parm
      
      implicit none
      
      integer j, idf
      
      integer   lyr, tnelem,  flag  !! flag=1: wood1 (branch) pools; flag = 2:wood2 (large wood) pools; flag = 3:wood3 (coarse root) pools; 
        
      real  :: ligcon , ps1co2(2), rsplig, tcflow, tcstva,netmnr(MAXIEL) !,som1e(MAXLYR,MAXIEL)
      real  :: aminrl(MAXIEL), rnew(MAXIEL,2), minerl(MAXLYR,MAXIEL)!, som2e(MAXLYR,MAXIEL)!, elstva(MACIEL), 
               !cstatv(nlr,ISOS),  gromin(MAXIEL),                &
               !resp(ISOS), som1ci(2,ISOS),           &
               !som2ci(2,ISOS), som2e(2,MAXIEL)
     
      real      cflowtosom1,cflowtosom2
      
      double precision newminrl

!! ...... Decompose stuff containing lignin (structural and wood).
!! ...... Call the compartment to be decomposed 'Box A'.
!! ...... C/E ratios of new material are computed once at the beginning of
!! ...... the simulation.

!! ...... Modified to flow to either soil or surface SOM2

!! ...... Input:
!! ......   aminrl(iel) - iel=1,3 available mineral N, P, S.
!! ......                 the sum of the positive layers of minerl 
!! ......                 before uptake by plants
!! ......   ligcon      - lignin content of Box A
!! ......   lyr         - layer (1=surface, 2=soil)
!! ......   nelem       - number of elements (besides C) to be simulated
!! ......   nlr         - total number of layers modeled for Box A;
!! ......                 (=2 for structural, metabolic, som1, som2;
!! ......                  =1 for som3, wood compartments)
!! ......   ps1co2(lyr) - controls amount of co2 loss when structural
!! ......                 decomposes to som1, subscripted for surface 
!! ......                 and soil layer (1)=surface, (2)=soil
!! ......   rnew(iel,1) - C/E ratio of new material being added to som1
!! ......                 (iel=1,3)
!! ......   rnew(iel,2) - C/E ratio of new material being added to som2
!! ......                 (iel=1,3)
!! ......   rsplig      - fraction of lignin flow lost to respiration
!! ......                 (rsplig was formerly named psco2l)
!! ......   tcflow      - total C flow out of Box A
!! ......   tcstva      - total C (unlabeled + labeled) in layer 'lyr' of
!! ......                 Box A.  For components with only 1 layer,
!! ......                 tcstva will be dimensioned (1).
!! ...... Transput:
!! ......   csrsnk(iso)     - C source/sink (iso=1,2)
!! ......   cstatv(lyr,iso) - C state variables for Box A;
!! ......                     cstatv(lyr,iso) represents C in layer 'lyr',
!! ......                     isotope 'iso' where iso=1 for unlabeled C
!! ......                     and iso=2 for labeled C.  For components with
!! ......                     only 1 layer, the 1st dimension of cstatv will
!! ......                     be 1.
!! ......   elstva(lyr,iel) - N, P, and S in Box A by layer and element
!! ......                     (iel=1,3)
!! ......   gromin(iel)     - gross mineralization (iel=1,3)
!! ......   minerl(1,iel)   - labile N, P, or S in layer 1 (iel=1,3).
!! ......   netmnr(lyr,iel) - net mineralization for layer lyr (N, P, or S)
!! ......                     For components with only 1 layer, the 1st
!! ......                     dimension of netmnr will be 1 (iel=1,3).
!! ......   resp(iso)       - output variable; annual accumulator for C flows
!! ......                     associated with microbial respiration (iso=1,2)
!! ......   som1ci(lyr,iso) - C by isotope in soil organic matter with fast
!! ......                     turnover rate (g/m2)
!! ......                     (iso=1,2) (1)=unlabeled; (2)=labeled
!! ......   som1e(lyr,iel)  - N, P, and S in soil organic matter with fast
!! ......                     turnover rate (g/m2)
!! ......                     (iel=1,3) (1)=N     (2)=P    (3)=S
!! ......   som2ci(lyr,iso) - C by isotope in soil organic matter with
!! ......                     intermediate turnover rate (g/m2)
!! ......                     (iso=1,2) (1) unlabeled; (2) labeled
!! ......   som2e(lyr,iel)  - N, P, and S in soil organic matter with
!! ......                     intermediate turnover rate (g/m2)
!! ......                     (iel=1,3) (1)=N     (2)=P     (3)=S

!! ...... Function declarations
      logical   candec
      external  candec

!! ...... Local variables
      integer   iel, actlyr, llyr, sol_lyr,nlr
      !!! integer   idf
      real      accum, co2los, mnrflo, tosom1, tosom2,rnew1(MAXIEL), rnew2(MAXIEL)
     
      real      boxa_ele(flag, MAXIEL)                      !! element in the 3 woody debris pools
      real, dimension (10):: tosom1_wood3,tosom2_wood3
       
     !! integer lyrr 
     
      j = 0
      idf = 0
            
      j = ihru
      idf = idplt(j)
      

      accum = 0.0
      
      !!lyrr = 0
      rnew1 = 0.      
      boxa_ele = 0.
      
      tosom1_wood3 = 0.
      tosom2_wood3 = 0.
      
      
      if (flag == 1) then  
      
            !!Store the C/E ratios for what goes to SOM1 in
            !!an array of correct size to be passed to CANDEC
            do 5 iel = 1, MAXIEL
                rnew1(iel) = rnew(iel,SRFC)
5           continue

            !!See if Box A can decompose to som1.
            !!If it can go to som1, it will also go to som2.
            !!If it can't go to som1, it can't decompose at all.

            do 6 iel = 1, MAXIEL
                boxa_ele(flag,iel) = woode1(iel,j)
6           continue      

            !!aminrl already in g/m2
            !!If Box A can decompose
            if (candec(tnelem,flag,aminrl,woodc1(j),boxa_ele,rnew1)) then

                !!Decompose Box A to SOM2
                !! -----------------------
                !!Gross C flow to som2
                tosom2 = tcflow * ligcon   !ligcon      - lignin content of Box A

                !!Respiration associated with decomposition to som2
                co2los = tosom2 * rsplig

                call respir(flag,co2los,lyr) !this function relates to isotope

                !!Net C flow to SOM2
                tosom2 = tosom2 - co2los
                !cflowtosom2 = tosom2        !! not accounted for yet, need to with Xuesong

                !!Partition and schedule C flows by isotope
                woodc1(j) =  woodc1(j) - tosom2       
       
                !!Zhang added
                sol_HSC(lyr,j) = sol_HSC(lyr,j)  + tosom2   * 10
                do 10 iel = 1, tnelem 
                
                    !call esched(tosom2,tcstva(lyr),rnew(iel,2), &
                    !            boxa_ele(lyr,iel),som2e(lyr,iel), &
                    !            minerl(lyr,iel))
                    !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    if (iel == 1) then
                        call esched(tosom2,woodc1(j),rnew(iel,2), &
                                    boxa_ele(lyr,iel),sol_HSN(lyr,j), &
                                    sol_nh3(lyr,j))
                        !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    end if
                    if (iel == 2) then
                        call esched(tosom2,woodc1(j),rnew(iel,2), &
                                    boxa_ele(lyr,iel),sol_orgp(lyr,j), &
                                    sol_solp(lyr,j))
                        !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    end if
                    if (iel .eq. N) then
                        newminrl = newminrl + mnrflo
                    endif
10              continue

                !!Decompose Box A to SOM1
                !! -----------------------

                !!Gross C flow to som1
                tosom1 = tcflow - tosom2 - co2los

                !!Respiration associated with decomposition to som1 
                co2los = tosom1 * ps1co2(lyr)

                call respir(flag,co2los,lyr)  !this function relates to isotope
     
     

                !!Net C flow to SOM1
                tosom1 = tosom1 - co2los
                cflowtosom1 = tosom1   !! not accounted for yet

                !!Partition and schedule C flows by isotope
                woodc1(j) =  woodc1(j) - tosom1
          
                !!Compute and schedule N, P, and S flows.

                !!Update mineralization accumulators.
                    
                !!Zhang added
                sol_BMC(sol_lyr,j) = sol_BMC(sol_lyr,j)  + tosom1*10 
                do 20 iel = 1, tnelem 
                    !call esched(tosom1,tcstva(lyr),rnew(iel,1), &
                    !            boxa_ele(lyr,iel),som2e(lyr,iel), &
                    !            minerl(lyr,iel))
                    !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    if (iel == 1) then
                        call esched(tosom2,woodc1(j),rnew(iel,1), &
                                    boxa_ele(lyr,iel),sol_BMN(lyr,j), &
                                    sol_nh3(lyr,j))
                        !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    end if
                    if (iel == 2) then
                        call esched(tosom2,woodc1(j),rnew(iel,1), &
                                    boxa_ele(lyr,iel),sol_orgp(lyr,j), &
                                    sol_solp(lyr,j))
                        !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    end if
                    if (iel .eq. N) then
                        newminrl = newminrl + mnrflo
                    endif
20              continue   
            end if
      end if       
       
       
       
      if (flag == 2) then  
      
            !!Store the C/E ratios for what goes to SOM1 in
            !!an array of correct size to be passed to CANDEC
            do 15 iel = 1, MAXIEL
                rnew1(iel) = rnew(iel,SRFC)
15           continue

            !!See if Box A can decompose to som1.
            !!If it can go to som1, it will also go to som2.
            !!If it can't go to som1, it can't decompose at all.

            do 16 iel = 1, MAXIEL
                boxa_ele(flag,iel) = woode2(iel,j)
16           continue      

            !!aminrl already in g/m2
            !!If Box A can decompose
            if (candec(tnelem,flag,aminrl,woodc2(j),boxa_ele,rnew1)) then

                !!Decompose Box A to SOM2
                !! -----------------------
                !!Gross C flow to som2
                tosom2 = tcflow * ligcon

                !!Respiration associated with decomposition to som2
                co2los = tosom2 * rsplig

                call respir(flag,co2los,lyr)

                !!Net C flow to SOM2
                tosom2 = tosom2 - co2los
                !cflowtosom2 = tosom2        !! not accounted for yet, need to with Xuesong

                !!Partition and schedule C flows by isotope
                woodc2(j) =  woodc2(j) - tosom2       
       
                !!Zhang added
                sol_HSC(lyr,j) = sol_HSC(lyr,j)  + tosom2  *10
                do 110 iel = 1, tnelem 
                    !call esched(tosom2,tcstva(lyr),rnew(iel,2), &
                    !            boxa_ele(lyr,iel),som2e(lyr,iel), &
                    !            minerl(lyr,iel))
                    !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    if (iel == 1) then
                        call esched(tosom2,woodc2(j),rnew(iel,2), &
                                    boxa_ele(lyr,iel),sol_HSN(lyr,j), &
                                    sol_nh3(lyr,j))
                        !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    end if
                    if (iel == 2) then
                        call esched(tosom2,woodc2(j),rnew(iel,2), &
                                    boxa_ele(lyr,iel),sol_orgp(lyr,j), &
                                    sol_solp(lyr,j))
                        !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    end if
                    if (iel .eq. N) then
                        newminrl = newminrl + mnrflo
                    endif
110              continue

                !!Decompose Box A to SOM1
                !! -----------------------

                !!Gross C flow to som1
                tosom1 = tcflow - tosom2 - co2los

                !!Respiration associated with decomposition to som1 
                co2los = tosom1 * ps1co2(lyr)

                call respir(flag,co2los,lyr)
     
     

                !!Net C flow to SOM1
                tosom1 = tosom1 - co2los
                cflowtosom1 = tosom1   !! not accounted for yet

                !!Partition and schedule C flows by isotope
                woodc2(j) =  woodc2(j) - tosom1
          
                !!Compute and schedule N, P, and S flows.

                !!Update mineralization accumulators.
                    
                !!Zhang added
                sol_BMC(sol_lyr,j) = sol_BMC(sol_lyr,j)  + tosom1  *10
                do 120 iel = 1, tnelem 
                    !call esched(tosom1,tcstva(lyr),rnew(iel,1), &
                    !            boxa_ele(lyr,iel),som2e(lyr,iel), &
                    !            minerl(lyr,iel))
                    !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    if (iel == 1) then
                        call esched(tosom1,woodc2(j),rnew(iel,1), &
                                    boxa_ele(lyr,iel),sol_BMN(lyr,j), &
                                    sol_nh3(lyr,j))
                        !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    end if
                    if (iel == 2) then
                        call esched(tosom1,woodc2(j),rnew(iel,1), &
                                    boxa_ele(lyr,iel),sol_orgp(lyr,j), &
                                    sol_solp(lyr,j))
                        !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                    end if
                    if (iel .eq. N) then
                        newminrl = newminrl + mnrflo
                    endif
120              continue   
            end if
      end if          
       
       
      if (flag == 3) then  

            !!actlyr (actual layer) is used when coarse roots are decomposed.
            !!In this case lyr = 2 and nlr = 1.  An error results when
            !!tcstva(lyr) is accessed because tcstva is dimensioned with nlr.
            !!Set actlyr to lyr and lyr to 1 for the special case.

            !!llyr (temporary variable) is used instead of 'lyr' so that it can
            !!be modified by the special case of SOIL=2 and nlr=1.

            Do lyr = 2, sol_nly(j)

                llyr = lyr
                actlyr = lyr
                
                if ((lyr .eq. 2) .and. (nlr .eq. 1)) llyr = 1
          
                !!Store the C/E ratios for what goes to SOM1 in
                !!an array of correct size to be passed to CANDEC
                do 7 iel = 1, MAXIEL
                    rnew2(iel) = rnew(iel,SOIL)
    7           continue

                !!See if Box A can decompose to som1.
                !!If it can go to som1, it will also go to som2.
                !!If it can't go to som1, it can't decompose at all.

                do 8 iel = 1, MAXIEL
                    !!need to uncomment
                    boxa_ele(flag,iel) = woode3(iel,lyr,j)
    8           continue

                !!aminrl already in g/m2
                !!If Box A can decompose
                !!need to uncomment
                !if (candec(tnelem,flag,aminrl,woodc(j,flag,lyr),boxa_ele,rnew2)) then
                if (candec(tnelem,flag,aminrl,woodc3(lyr,j),boxa_ele,rnew2)) then
                    !!Decompose Box A to SOM2
                    !! -----------------------
                    !!Gross C flow to som2
                    tosom2 = tcflow * ligcon

                    !!Respiration associated with decomposition to som2
                    co2los = tosom2 * rsplig

                    call respir(flag,co2los,lyr)
                    !respir(flag,co2los,nlr,lyr)
                    !!Net C flow to SOM2
                    tosom2 = tosom2 - co2los
                    !cflowtosom2 = tosom2        !! not accounted for yet, need to with Xuesong

                    !!need to uncomment
                    !woodc(j,flag,lyr) =  woodc(j,flag,lyr) - tosom2       
                    woodc3(lyr,j) =  woodc3(lyr,j) - tosom2  
                    
                    !!Zhang added
                    sol_HSC(lyr,j) = sol_HSC(lyr,j)  + tosom2  *10
                    do 30 iel = 1, tnelem 

                        if (iel == N) then
                            call esched(tosom2,woodc3(lyr,j),rnew(iel,2), &
                                        boxa_ele(lyr,iel),sol_HSN(lyr,j), &
                                        sol_nh3(lyr,j))
                            !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                        end if
                        if (iel == P) then
                            call esched(tosom2,woodc3(lyr,j),rnew(iel,2), &
                                        boxa_ele(lyr,iel),sol_orgp(lyr,j), &
                                        sol_solp(lyr,j))
                            !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                        end if
                        
                        if (iel .eq. N) then
                            newminrl = newminrl + mnrflo
                        endif
    30              continue

                    !!Decompose Box A to SOM1
                    !! -----------------------

                    !!Gross C flow to som1
                    tosom1 = tcflow - tosom2 - co2los

                    !!Respiration associated with decomposition to som1 
                    co2los = tosom1 * ps1co2(lyr)

                    call respir(flag,co2los,lyr)     
         
                    !!Net C flow to SOM1
                    tosom1 = tosom1 - co2los
                    cflowtosom1 = tosom1   !! not accounted for yet

                    !!Partition and schedule C flows by isotope
                    !!need to uncomment
                    !woodc(j,flag,lyr) =  woodc(j,flag,lyr) - tosom1
                    woodc3(lyr,j) =  woodc3(lyr,j) - tosom1
                    !!Compute and schedule N, P, and S flows.

                    !!Update mineralization accumulators.
                        
                    !!Zhang added
                    sol_BMC(sol_lyr,j) = sol_BMC(sol_lyr,j)  + tosom1  *10
                    do 40 iel = 1, tnelem 
                        if (iel == N) then
                            call esched(tosom1,woodc3(lyr,j),rnew(iel,1), &
                                        boxa_ele(lyr,iel),sol_BMN(lyr,j), &
                                        sol_nh3(lyr,j))
                            !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                        end if
                        if (iel == P) then
                            call esched(tosom1,woodc3(lyr,j),rnew(iel,1), &
                                        boxa_ele(lyr,iel),sol_orgp(lyr,j), &
                                        sol_solp(lyr,j))
                            !call mnracc(mnrflo,gromin(iel),netmnr(lyr,iel))
                        end if
                                            
                        if (iel .eq. N) then
                            newminrl = newminrl + mnrflo
                        endif
    40              continue   

                end if         
       
        END DO

      endif


      return
      end