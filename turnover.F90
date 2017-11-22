      subroutine turnover



      use parm
      !!implicit none
      real:: avc, maxrate, tmprate
      real::retivelai, systemc, systemn, totalc, totaln
      real:: acvegc
      real :: tlfc, tlfn, tbrchc, tbrchn, tcrtc, tcrtn, tfrtjc, tfrtjn
      real :: tfrtmc, tfrtmn
      real :: tseedc, tseedn, tlargwc, tlargwn
      real :: ltlfc, ltlfn, ltbrchc, ltbrchn, ltcrtc, ltcrtn, ltfrtc,ltfrtn, ltlargwc, ltlargwn
      real :: ltseedc, ltseedn
      real:: lai, tlai
      integer :: j, idf, yday
      real:: manurec, manuren, manurecn
      real:: tolargwc, tolargwn,rnlargwn
      real:: frtolitn,lftolitn,lftolab
      real:: phen_yday, phen_day
      real:: hvtolitterc,hvtolittern,pharvestc,pharvestn
      real:: ltharvestc, ltharvestn, templfc

      real:: pothw_n, avnhw_n
      real:: day_stage

      !!real, dimension (20):: prdis
      real :: BLG1, BLG2, BLG3,  CLG, sf
      real :: sol_min_n,  resnew_n, resnew_ne
      real :: LMF, LSF, LSLF, LSNF,LMNF 
      real :: no3_t, nh4_t, volatilziation, lmn_t, lsn_t,lmc_t,lsc_t
      real :: cal_phen
      integer::get_stage
      real :: cleafmax_pot,toheartwoodc, returnsapn, toheartwoodn
      real :: n11, n22, n33, n44, n55, n66, n77, n88, n99, n1010
      
      
      orgc_f = 0.
      BLG1 = 0.
      BLG2 = 0.
      BLG3 = 0.
      CLG = 0.
      sf = 0.
      sol_min_n = 0.
      resnew = 0.
      resnew_n = 0.
      resnew_ne = 0.
      LMF = 0.
      LSF = 0.
      LSLF = 0.
      LSNF = 0.
      LMNF = 0.


      j = 0
      j = ihru
      idf =idplt(j)
      maxrate=0.


        tlfc = leafc(j)
        tlfn = leafn(j)
        tbrchc = brchc(j)
        tbrchn = brchn(j)
        tlargwc = largwc(j)
        tlargwn = largwn(j)
        tcrtc = csrootc(j)
        tcrtn = csrootn(j)
        tfrtjc = frootjc(j)
        tfrtmc = frootmc(j)
       !! tfrtn = frootn(j)
        tfrtjn = frootjn(j)
        tfrtmn = frootmn(j)
        tseedc = seedc(j)
        tseedn = seedn(j)

        

       ltlfc=0.
       ltlfn=0.
       ltbrchc = 0.
       ltbrchn = 0.
       ltlargwc = 0.
       ltlargwn = 0.
       ltcrtc = 0.
       ltcrtn = 0.
       ltfrtc = 0.
       ltfrtn = 0.
       ltseedc = 0.
       ltseedn = 0.

        manuren = 0.
        manurec = 0.
        manurecn = 15.0

        tolargwn = 0.
        rnlargwn = 0.
        lflttolab = 1.0

        hvtolitterc = 0.
        hvtolittern = 0.

        pharvestc = 0.
        pharvestn = 0.

        ltharvestc = 0.
        ltharvestn = 0.
        cleafmax_pot = 0.
        toheartwoodc = 0.
        returnsapn = 0.
        toheartwoodn= 0.
        frtolitn = 0.
        lftolitn = 0.

    
     if(acvegc > inicdf(idf)) then

            !!add by zhang
            !!===================
            if (cswat == 2) then
                rsdc_d(j) = rsdc_d(j) +  ltlfc !! change later
                !rsdc_d(j) = rsdc_d(j) + resnew*0.42
            end if
            !!add by zhang
            !!===================

            !!insert new biomss by zhang
            !!=============================
            if (cswat == 2) then
            
                !!===============adding leaf and seeds
	          !!all the lignin from STD is assigned to LSL, 
	            !!add STDL calculation


	          !if (k == 1) then
		        sf = 0.05
	          !else
		        !sf = 0.1
	          !end if	

               !kg/ha  
	          sol_min_n = 0.	
	          sol_min_n = (sol_no3(1,j)+sol_nh3(1,j))
	          	          
	          resnew = ((ltlfc+ltseedc)*10)/0.45 !resnew !bm_dieoff(idplt(j)) * bio_ms(j) 
	          resnew_n = (lftolitn+ltseedn)*10 !bm_dieoff(idplt(j)) * plantn(j)   	    
        	    resnew_ne = resnew_n + sf * sol_min_n
                
                !update no3 and nh3 in soil
                sol_no3(1,j) = sol_no3(1,j) * (1-sf)
                sol_nh3(1,j) = sol_nh3(1,j) * (1-sf)

        	     !Not sure 1000 should be here or not!
        	    !RLN = 1000*(resnew * CLG/(resnew_n+1.E-5))
        	    CLG = 0.15
        	    RLN = (resnew * CLG/(resnew_n+1.E-5))
        	    RLR = MIN(.8, resnew * CLG/1000/(resnew/1000+1.E-5))
        	    
        	    LMF = 0.85 - 0.018 * RLN
        	    if (LMF <0.01) then
        	        LMF = 0.01
        	    else
        	        if (LMF >0.7) then
        	            LMF = 0.7
        	        end if
        	    end if
	          !if ((resnew * CLG/(resnew_n+1.E-5)) < 47.22) then
		        !    LMF = 0.85 - 0.018 * (resnew * CLG/(resnew_n+1.E-5))
	          !else
		        !    LMF = 0.
	          !end if 	

	          LSF =  1 - LMF  
        	  
	          sol_LM(1,j) = sol_LM(1,j) + LMF * resnew
	          sol_LS(1,j) = sol_LS(1,j) + LSF * resnew
        	    

                
	          !here a simplified assumption of 0.5 LSL
	          !LSLF = 0.0
	          !LSLF = CLG          
	          
	          sol_LSL(1,j) = sol_LSL(1,j) + RLR* resnew	          
	          sol_LSC(1,j) = sol_LSC(1,j) + LSF*resnew*0.45 
	          
	          sol_LSLC(1,j) = sol_LSLC(1,j) + RLR*resnew*0.45 
	          sol_LSLNC(1,j) = sol_LSC(1,j) - sol_LSLC(1,j)              
                
                !X3 = MIN(X6,0.42*LSF * resnew/150) 
                
	          if (resnew_ne >= (LSF * (resnew*0.45 ) /150)) then
		         sol_LSN(1,j) = sol_LSN(1,j) + LSF * (resnew*0.45 ) / 150
		         sol_LMN(1,j) = sol_LMN(1,j) + resnew_ne - (LSF * (resnew*0.45 ) / 150) + 1.E-25
	          else
		         sol_LSN(1,j) = sol_LSN(1,j) + resnew_ne
		         sol_LMN(1,j) = sol_LMN(1,j) + 1.E-25
	          end if	
        	
	          !LSNF = sol_LSN(1,j)/(sol_LS(1,j)+1.E-5)	
        	  
	          sol_LMC(1,j) = sol_LMC(1,j) + LMF * (resnew*0.45 )	
	          !LMNF = sol_LMN(1,j)/(sol_LM(1,j) + 1.E-5)           
                

            
            !!insert new biomss by zhang
            !!===========================


            sol_rsd(1,j) = sol_rsd(1,j) + resnew
            sol_rsd(1,j) = Max(sol_rsd(1,j),0.)
            sol_fon(1,j) = sol_fon(1,j) + resnew_ne
            sol_fop(1,j) = sol_fop(1,j) + resnew_ne/20
            !!===============adding leaf and seeds


                 !!===============adding stem
	          !!all the lignin from STD is assigned to LSL, 
	            !!add STDL calculation


	          !if (k == 1) then
		        sf = 0.05
	          !else
		        !sf = 0.1
	          !end if	

               !kg/ha  
	          sol_min_n = 0.	
	          sol_min_n = (sol_no3(1,j)+sol_nh3(1,j))
          	          
	          	          
	          resnew = ((ltbrchc+ltlargwc)*10)/0.45 !resnew !bm_dieoff(idplt(j)) * bio_ms(j) 
	          resnew_n = (ltbrchn+ltlargwn)*10 !bm_dieoff(idplt(j)) * plantn(j)   	    
        	    resnew_ne = resnew_n + sf * sol_min_n
                
                !update no3 and nh3 in soil
                sol_no3(1,j) = sol_no3(1,j) * (1-sf)
                sol_nh3(1,j) = sol_nh3(1,j) * (1-sf)

        	     !Not sure 1000 should be here or not!
        	    !RLN = 1000*(resnew * CLG/(resnew_n+1.E-5))
        	    CLG = 0.30
        	    RLN = (resnew * CLG/(resnew_n+1.E-5))
        	    RLR = MIN(.8, resnew * CLG/1000/(resnew/1000+1.E-5))
        	    
        	    LMF = 0.85 - 0.018 * RLN
        	    if (LMF <0.01) then
        	        LMF = 0.01
        	    else
        	        if (LMF >0.7) then
        	            LMF = 0.7
        	        end if
        	    end if
	          !if ((resnew * CLG/(resnew_n+1.E-5)) < 47.22) then
		        !    LMF = 0.85 - 0.018 * (resnew * CLG/(resnew_n+1.E-5))
	          !else
		        !    LMF = 0.
	          !end if 	

	          LSF =  1 - LMF  
        	  
	          sol_LM(1,j) = sol_LM(1,j) + LMF * resnew
	          sol_LS(1,j) = sol_LS(1,j) + LSF * resnew
        	    

                
	          !here a simplified assumption of 0.5 LSL
	          !LSLF = 0.0
	          !LSLF = CLG          
	          
	          sol_LSL(1,j) = sol_LSL(1,j) + RLR* resnew	          
	          sol_LSC(1,j) = sol_LSC(1,j) + LSF*resnew*0.45 
	          
	          sol_LSLC(1,j) = sol_LSLC(1,j) + RLR*resnew*0.45 
	          sol_LSLNC(1,j) = sol_LSC(1,j) - sol_LSLC(1,j)              
                
                !X3 = MIN(X6,0.42*LSF * resnew/150) 
                
	          if (resnew_ne >= (LSF * (resnew*0.45 ) /150)) then
		         sol_LSN(1,j) = sol_LSN(1,j) + LSF * (resnew*0.45 ) / 150
		         sol_LMN(1,j) = sol_LMN(1,j) + resnew_ne - (LSF * (resnew*0.45 ) / 150) + 1.E-25
	          else
		         sol_LSN(1,j) = sol_LSN(1,j) + resnew_ne
		         sol_LMN(1,j) = sol_LMN(1,j) + 1.E-25
	          end if	
        	
	          !LSNF = sol_LSN(1,j)/(sol_LS(1,j)+1.E-5)	
        	  
	          sol_LMC(1,j) = sol_LMC(1,j) + LMF * (resnew*0.45 )	
	          !LMNF = sol_LMN(1,j)/(sol_LM(1,j) + 1.E-5)           
                

           
            !!insert new biomss by zhang
            !!===========================


            sol_rsd(1,j) = sol_rsd(1,j) + resnew
            sol_rsd(1,j) = Max(sol_rsd(1,j),0.)
            sol_fon(1,j) = sol_fon(1,j) + resnew_ne
            sol_fop(1,j) = sol_fop(1,j) + resnew_ne/20
            !!===============adding stem         


        !!===============adding roots
	          !!all the lignin from STD is assigned to LSL, 
	            !!add STDL calculation
	            total_root = 0.
!!	      do ly = 2, sol_nly(j)
!!       
!!              if (ly<sol_nly(j)) then 
!!                    prdis(ly) = 0.5*(exp(-raf(idf)*sol_z(ly-1,j)/1000.)+exp(-rbf(idf)*sol_z(ly-1,j)/1000.)-exp(-raf(idf)*sol_z(ly,j)/1000.)- exp(-rbf(idf)*sol_z(ly,j)/1000.))
!!               else
!!                    prdis(ly) = 0.5*(exp(-raf(idf)*sol_z(ly-1,j)/1000.)+exp(-rbf(idf)*sol_z(ly-1,j)/1000.))               
!!               end if 
               
!!             rdis(ly) = prdis(ly)
!!               total_root = total_root + rdis(ly)
!!               !write (*,*) rdis(ly)
!!           end do   
          do ly = 2, sol_nly(j)
       
!!               if (ly<sol_nly(j)) then 
!!                    prdis(ly) = 0.5*(exp(-raf(idf)*sol_z(ly-1,j)/1000.)+exp(-rbf(idf)*sol_z(ly-1,j)/1000.)-exp(-raf(idf)*sol_z(ly,j)/1000.)- exp(-rbf(idf)*sol_z(ly,j)/1000.))
!!               else
!!                    prdis(ly) = 0.5*(exp(-raf(idf)*sol_z(ly-1,j)/1000.)+exp(-rbf(idf)*sol_z(ly-1,j)/1000.))               
!!               end if 
               
!!               rdis(ly) = prdis(ly)/total_root
 
 	              !if (k == 1) then
		            sf = 0.05
	              !else
		            !sf = 0.1
	              !end if	

                   !kg/ha  
	              sol_min_n = 0.	
	              sol_min_n = (sol_no3(ly,j)+sol_nh3(ly,j))
   	          	          
	              resnew = ((ltcrtc+ltfrtc)*10)/0.45 * rdis(ly,j) !resnew !bm_dieoff(idplt(j)) * bio_ms(j) 
	              resnew_n = (ltcrtn+frtolintn)*10 * rdis(ly,j) !bm_dieoff(idplt(j)) * plantn(j)   	    
        	        resnew_ne = resnew_n + sf * sol_min_n
                    
                    !update no3 and nh3 in soil
                    sol_no3(ly,j) = sol_no3(ly,j) * (1-sf)
                    sol_nh3(ly,j) = sol_nh3(ly,j) * (1-sf)

        	         !Not sure 1000 should be here or not!
        	        !RLN = 1000*(resnew * CLG/(resnew_n+1.E-5))
        	        CLG = 0.20
        	        RLN = (resnew * CLG/(resnew_n+1.E-5))
        	        RLR = MIN(.8, resnew * CLG/1000/(resnew/1000+1.E-5))
            	    
        	        LMF = 0.85 - 0.018 * RLN
        	        if (LMF <0.01) then
        	            LMF = 0.01
        	        else
        	            if (LMF >0.7) then
        	                LMF = 0.7
        	            end if
        	        end if
	              !if ((resnew * CLG/(resnew_n+1.E-5)) < 47.22) then
		            !    LMF = 0.85 - 0.018 * (resnew * CLG/(resnew_n+1.E-5))
	              !else
		            !    LMF = 0.
	              !end if 	

	              LSF =  1 - LMF  
            	  
	              sol_LM(ly,j) = sol_LM(ly,j) + LMF * resnew
	              sol_LS(ly,j) = sol_LS(ly,j) + LSF * resnew
            	    

                    
	              !here a simplified assumption of 0.5 LSL
	              !LSLF = 0.0
	              !LSLF = CLG          
    	          
	              sol_LSL(ly,j) = sol_LSL(ly,j) + RLR* resnew	          
	              sol_LSC(ly,j) = sol_LSC(ly,j) + LSF*resnew*0.45 
    	          
	              sol_LSLC(ly,j) = sol_LSLC(ly,j) + RLR*resnew*0.45 
	              sol_LSLNC(ly,j) = sol_LSC(ly,j) - sol_LSLC(ly,j)              
                    
                    !X3 = MIN(X6,0.42*LSF * resnew/150) 
                    
	              if (resnew_ne >= (LSF * (resnew*0.45 ) /150)) then
		             sol_LSN(ly,j) = sol_LSN(ly,j) + LSF * (resnew*0.45 ) / 150
		             sol_LMN(ly,j) = sol_LMN(ly,j) + resnew_ne - (LSF * (resnew*0.45 ) / 150) + 1.E-25
	              else
		             sol_LSN(ly,j) = sol_LSN(ly,j) + resnew_ne
		             sol_LMN(ly,j) = sol_LMN(ly,j) + 1.E-25
	              end if	
            	
	              !LSNF = sol_LSN(ly,j)/(sol_LS(ly,j)+1.E-5)	
            	  
	              sol_LMC(ly,j) = sol_LMC(ly,j) + LMF * (resnew*0.45 )	
	              !LMNF = sol_LMN(1,j)/(sol_LM(1,j) + 1.E-5)           

                !!insert new biomss by zhang
                !!===========================


                sol_rsd(ly,j) = sol_rsd(ly,j) + resnew
                sol_rsd(ly,j) = Max(sol_rsd(ly,j),0.)
                sol_fon(ly,j) = sol_fon(ly,j) + resnew_ne
                sol_fop(ly,j) = sol_fop(ly,j) + resnew_ne/20
                !!===============adding roots 
            end do

            end if

        else

                lai = leafc(j) * avg_slaf(idf)

                if(day_stage == 2 .and. phen_yday > 0.) then

                tlai = lai*(phen_day / phen_yday)
                end if

            if(lai>=tlai .and. lai > 0. .and. day_stage == 2) then
                        ltlfc = (lai-tlai)/avg_slaf(idf)
            end if

                if(tlfc>0.) ltlfn = (ltlfc/tlfc) * tlfn;

               leafc(j)= leafc(j)-ltlfc
               leafn(j)= leafn(j)-ltlfn

     

            !!insert new biomss by zhang
            !!=============================
            if (cswat == 2) then
            
                !!===============adding leaf and seeds
	          !!all the lignin from STD is assigned to LSL, 
	            !!add STDL calculation


	          !if (k == 1) then
		        sf = 0.05
	          !else
		        !sf = 0.1
	          !end if	

               !kg/ha  
	          sol_min_n = 0.	
	          sol_min_n = (sol_no3(1,j)+sol_nh3(1,j))
	          	          
	          resnew = ((ltlfc)*10)/0.45 !resnew !bm_dieoff(idplt(j)) * bio_ms(j) 
	          resnew_n = (ltlfn)*10 !bm_dieoff(idplt(j)) * plantn(j)   	    
        	    resnew_ne = resnew_n + sf * sol_min_n
                
                !update no3 and nh3 in soil
                sol_no3(1,j) = sol_no3(1,j) * (1-sf)
                sol_nh3(1,j) = sol_nh3(1,j) * (1-sf)

        	     !Not sure 1000 should be here or not!
        	    !RLN = 1000*(resnew * CLG/(resnew_n+1.E-5))
        	    CLG = 0.15
        	    RLN = (resnew * CLG/(resnew_n+1.E-5))
        	    RLR = MIN(.8, resnew * CLG/1000/(resnew/1000+1.E-5))
        	    
        	    LMF = 0.85 - 0.018 * RLN
        	    if (LMF <0.01) then
        	        LMF = 0.01
        	    else
        	        if (LMF >0.7) then
        	            LMF = 0.7
        	        end if
        	    end if
	          !if ((resnew * CLG/(resnew_n+1.E-5)) < 47.22) then
		        !    LMF = 0.85 - 0.018 * (resnew * CLG/(resnew_n+1.E-5))
	          !else
		        !    LMF = 0.
	          !end if 	

	          LSF =  1 - LMF  
        	  
	          sol_LM(1,j) = sol_LM(1,j) + LMF * resnew
	          sol_LS(1,j) = sol_LS(1,j) + LSF * resnew
        	    

                
	          !here a simplified assumption of 0.5 LSL
	          !LSLF = 0.0
	          !LSLF = CLG          
	          
	          sol_LSL(1,j) = sol_LSL(1,j) + RLR* resnew	          
	          sol_LSC(1,j) = sol_LSC(1,j) + LSF*resnew*0.45 
	          
	          sol_LSLC(1,j) = sol_LSLC(1,j) + RLR*resnew*0.45 
	          sol_LSLNC(1,j) = sol_LSC(1,j) - sol_LSLC(1,j)              
                
                !X3 = MIN(X6,0.42*LSF * resnew/150) 
                
	          if (resnew_ne >= (LSF * (resnew*0.45 ) /150)) then
		         sol_LSN(1,j) = sol_LSN(1,j) + LSF * (resnew*0.45 ) / 150
		         sol_LMN(1,j) = sol_LMN(1,j) + resnew_ne - (LSF * (resnew*0.45 ) / 150) + 1.E-25
	          else
		         sol_LSN(1,j) = sol_LSN(1,j) + resnew_ne
		         sol_LMN(1,j) = sol_LMN(1,j) + 1.E-25
	          end if	
        	
	          !LSNF = sol_LSN(1,j)/(sol_LS(1,j)+1.E-5)	
        	  
	          sol_LMC(1,j) = sol_LMC(1,j) + LMF * (resnew*0.45 )	
	          !LMNF = sol_LMN(1,j)/(sol_LM(1,j) + 1.E-5)           
                

            
            !!insert new biomss by zhang
            !!===========================


            sol_rsd(1,j) = sol_rsd(1,j) + resnew
            sol_rsd(1,j) = Max(sol_rsd(1,j),0.)
            sol_fon(1,j) = sol_fon(1,j) + resnew_ne
            sol_fop(1,j) = sol_fop(1,j) + resnew_ne/20
            !!===============adding leaf and seeds
            end if




        end if

            no3_t = 0.
            nh4_t = 0.
            volatilziation = 0.
            lmn_t = 0.
            lsn_t = 0.
            lmc_t = 0.
            lsc_t = 0.
            !deni = 0.
            !leaching = 0.
            !surfqn = 0.
             do ly = 1, sol_nly(j)
                no3_t = no3_t + sol_no3(ly,j) 
                nh4_t = nh4_t + sol_nh3(ly,j)
                lmn_t = lmn_t + sol_LMN(ly,j)
                lsn_t = lsn_t + sol_LSN(ly,j)
                lmc_t = lmc_t + sol_LMC(ly,j)
                lsc_t = lsc_t + sol_LSC(ly,j)                
            end do
       
       !!cal_temp(10)= ltlfc + ltcrtc + ltfrtc + ltseedc + ltbrchc + ltlargwc
    
        cal_temp(10)= lftolitn + ltcrtn + frtolitn + ltseedn + ltbrchn + ltlargwn
               
       
        
       return
       end

