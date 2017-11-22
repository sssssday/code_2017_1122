



      subroutine respir(flag,co2los,lyr)

      
      use parm
      implicit none


!! ......Argument declarations
      integer   nlr, lyr,flag 
      real      co2los, tcstva, cstatv(ISOS), resp(ISOS), estatv
     
      double precision newminrl

!! ......Compute flows associated with microbial respiration.

!! ......Input:
!! ......  co2los      = CO2 loss associated with decomposition
!! ......  nlr         = total number of layers modeled for Box A;
!! ......                (=2 for structural, metabolic, som1, som2;
!! ......                 =1 for som3, wood compartments)
!! ......  lyr         = layer for which respiration is being computed
!! ......  tcstva(nlr) = total C (unlabeled + labeled) in layer 'lyr' of 
!! ......                Box A.  For components with only 1 layer,
!! ......                tcstva will be dimensioned (1).
!! ......  time        = passed in /zztim/
!! ...
!! ......Transput:
!! ......  cstatv          = C state variables for Box A;
!! ......                    cstatv(lyr,iso) represents C in layer 'lyr',
!! ......                    isotope 'iso' where iso=1 for unlabeled C
!! ......                    and iso=2 for labeled C.  For components with
!! ......                    only 1 layer, the 1st dimension of cstatv will
!! ......                    be 1.
!! ......  csrsnk          = C source/sink 
!! ......  resp            = output variable
!! ......  estatv          = N, P, and S variables for Box A.  In 'estatv(lyr,iel)',
!! ......                    lyr represents layer and iel indicates N, P, or S.  For
!! ......                    components with only 1 layer, the 1st dimension of
!! ......                    estatv will be 1.
!! ......  minerl(1,iel)   = labile N, P, or S in layer 1.
!! ......  gromin(iel)     = gross mineralization -- stored in /comput/
!! ......  netmnr(lyr,iel) = net mineralization for layer lyr (N, P, or S)
!! ......                    For components with only 1 layer, the 1st
!! ......                    dimension of netmnr will be 1.

!! ......Local variables
      integer j, idf
      integer   iel
      real      mnrflo

!! ......C flow from cstatv to CO2
           
      j = 0
      idf = 0
            
      j = ihru
      idf = idplt(j)


      csrsnk (j) = csrsnk (j) + co2los
      rspc_d(j) = rspc_d(j) + co2los*10.  !*10 to conter from g/m2 to kg/ha
      if (flag == 1) then
        woodc1(j) =  woodc1(j) - co2los
      end if
      if (flag == 2) then
        woodc2(j) =  woodc2(j) - co2los
      end if
      if (flag == 3) then
        woodc3(lyr,j) =  woodc3(lyr,j) - co2los
      end if            
   !!   if (labtyp .eq. 2) then
   !!     call csched(co2los,cstatv(lyr,LABELD),tcstva(lyr),
   !!  &              cstatv(lyr,UNLABL),csrsnk(UNLABL),
   !!  &              cstatv(lyr,LABELD),csrsnk(LABELD),
   !!  &              dresp,resp)
   !!   else
   !!!!     call csched(co2los,cstatv(lyr,LABELD),tcstva(lyr),
   !!  &              cstatv(lyr,UNLABL),csrsnk(UNLABL),
   !!  &              cstatv(lyr,LABELD),csrsnk(LABELD),
   !!  &              1.0,resp)
   !!   endif
 
!! ......Mineralization associated with respiration
      do 10 iel = 1, nelem
        if (flag == 1) then
            mnrflo = co2los*woode1(iel,j)/woodc1(j)
            woode1(iel,j) = woode1(iel,j) - mnrflo
        end if
        if (flag == 2) then
            mnrflo = co2los*woode2(iel,j)/woodc2(j)
            woode2(iel,j) = woode2(iel,j) - mnrflo
        end if        

        if (flag == 3) then
            mnrflo = co2los*woode3(iel,lyr,j)/woodc3(lyr,j)
            woode3(iel,lyr,j) = woode3(iel,lyr,j) - mnrflo
        end if  
        
        if (iel == N) then
          sol_no3 (lyr,j) = sol_no3 (lyr,j) +  mnrflo  &
         * 10. /(sol_no3 (lyr,j) + sol_nh3 (lyr,j) )
        
          sol_nh3 (lyr,j) = sol_nh3 (lyr,j) +  mnrflo * 10. / &
         (sol_no3 (lyr,j) + sol_nh3 (lyr,j) )
     
        else 
          sol_solp(lyr,j) = sol_solp(lyr,j) + mnrflo * 10.
        
        !!zXXXXXX add P code here
       
        end if
        
        
       !! call flow(estatv(lyr,iel),minerl(SRFC,iel),time,mnrflo)

!! ..... Update gross mineralization
       !! gromin(iel) = gromin(iel) + mnrflo

!!..... Update net mineralization
!!        netmnr(lyr,iel) = netmnr(lyr,iel) + mnrflo
!! ..... newminrl should be updated only for nitrogen
!! ..... akm via cak 07/31/01
!!        if (iel .eq. N) then
!!          newminrl = newminrl + mnrflo
 !!       endif
      
10    continue
 
      return
      end
