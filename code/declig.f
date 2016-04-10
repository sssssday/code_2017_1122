
     
      subroutine declig(flag, aminrl,ligcon,lyr,tnelem,nlr,ps1co2,rnew,
     &                  rsplig,tcflow)
   !!  ,
   !!   &                  gromin,minerl,netmnr,resp,som1ci,som1e,
   !!  &                  som2ci,som2e,cflowtosom1,cflowtosom2,newminrl)





      use parm
      
      implicit none
      
      integer j, idf
      

      integer   lyr, tnelem, nlr, flag
      real      aminrl(MAXIEL), ligcon, ps1co2(2), rnew(MAXIEL,2),
     &          rsplig, tcflow, tcstva(nlr), 
     &          cstatv(nlr,ISOS), elstva, gromin(MAXIEL),
     &          minerl(MAXLYR,MAXIEL), 
     &          netmnr(nlr,MAXIEL),resp(ISOS), som1ci(2,ISOS),
     &          som1e(2,MAXIEL), som2ci(2,ISOS), som2e(2,MAXIEL)
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
      integer   iel, actlyr, llyr
      !!! integer   idf
      real      accum, co2los, mnrflo, tosom1, tosom2,
     &          rnew1(MAXIEL)
       
     !! integer lyrr 
     
      j = 0
      idf = 0
            
      j = ihru
      idf = idplt(j)
      

      accum = 0.0
      
      !!lyrr = 0
      rnew1 = 0.
      

!! ...... actlyr (actual layer) is used when coarse roots are decomposed.
!! ...... In this case lyr = 2 and nlr = 1.  An error results when
!! ...... tcstva(lyr) is accessed because tcstva is dimensioned with nlr.
!! ...... Set actlyr to lyr and lyr to 1 for the special case.

!! ...... llyr (temporary variable) is used instead of 'lyr' so that it can
!! ...... be modified by the special case of SOIL=2 and nlr=1.

      llyr = lyr
      actlyr = lyr
      if ((lyr .eq. 2) .and. (nlr .eq. 1)) llyr = 1

!! ...... Store the C/E ratios for what goes to SOM1 in
!! ...... an array of correct size to be passed to CANDEC
      do 5 iel = 1, MAXIEL
        rnew1(iel) = rnew(iel,llyr)
5     continue

!! ...... See if Box A can decompose to som1.
!! ...... If it can go to som1, it will also go to som2.
!! ...... If it can't go to som1, it can't decompose at all.

!! ...... If Box A can decompose

      if (candec(tnelem,aminrl,woodc(j,flag),woode(j,flag,1),nlr,
     &           llyr,rnew1)) then

!! ........         Decompose Box A to SOM2
!! ........         -----------------------
!! ........ Gross C flow to som2
        tosom2 = tcflow * ligcon

!! ........ Respiration associated with decomposition to som2
        co2los = tosom2 * rsplig

        call respir(flag, co2los,nlr,llyr)

!! ........ Net C flow to SOM2
        tosom2 = tosom2 - co2los
        cflowtosom2 = tosom2        !! not accounted for yet

!! ........ Partition and schedule C flows by isotope

       woodc(j,flag) =  woodc(j,flag) - tosom2
        !!XXXXXXXXX sink not accounted yet
!!        call csched(tosom2,cstatv(llyr,LABELD),tcstva(llyr),
!!     &              cstatv(llyr,UNLABL),som2ci(actlyr,UNLABL),
!!    &              cstatv(llyr,LABELD),som2ci(actlyr,LABELD),
!!     &              1.0,accum)

!! ........ Compute and schedule N, P, and S flows.

!! ........ Update mineralization accumulators.
        do 10 iel = 1, tnelem
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  not account fall yet
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX          call esched(tosom2,tcstva(llyr),rnew(iel,2),
 !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    &                elstva(llyr,iel),som2e(actlyr,iel),
 !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    &                minerl(1,iel),mnrflo)
 !!         call mnracc(mnrflo,gromin(iel),netmnr(llyr,iel))
!! .......... newminrl should be updated only for nitrogen
!! .......... akm via cak 07/31/01
 !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX         if (iel .eq. N) then
  !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX          newminrl = newminrl + mnrflo
  !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        endif
10      continue

!! ........         Decompose Box A to SOM1
!! ........         -----------------------

!! ........ Gross C flow to som1
        tosom1 = tcflow - tosom2 - co2los

!! ........ Respiration associated with decomposition to som1 
        co2los = tosom1 * ps1co2(llyr)

        call respir(flag,co2los,nlr,llyr)
     
     

!! ........ Net C flow to SOM1
        tosom1 = tosom1 - co2los
        cflowtosom1 = tosom1   !! not accounted for yet

!! ........ Partition and schedule C flows by isotope

          woodc(j,flag) =  woodc(j,flag) - tosom1
          !!XXXXXXXXX sink not accounted yet
 !!       call csched(tosom1,cstatv(llyr,LABELD),tcstva(llyr),
 !!    &              cstatv(llyr,UNLABL),som1ci(actlyr,UNLABL),
 !!    &              cstatv(llyr,LABELD),som1ci(actlyr,LABELD),
 !!    &              1.0,accum)

!! ........ Compute and schedule N, P, and S flows.

!! ........ Update mineralization accumulators.
        do 20 iel = 1, tnelem
!! not accounted for yet
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX          call esched(tosom1,tcstva(llyr),rnew(iel,1),
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     &                elstva(llyr,iel),som1e(actlyr,iel),
 !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    &                minerl(1,iel),mnrflo)
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX          call mnracc(mnrflo,gromin(iel),netmnr(llyr,iel))
!! .......... newminrl should be updated only for nitrogen
!! .......... akm via cak 07/31/01
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX          if (iel .eq. N) then
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX            newminrl = newminrl + mnrflo
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX          endif
20      continue

      endif

      return
      end