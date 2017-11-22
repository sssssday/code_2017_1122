
      subroutine partit(cpart,recres,lyr,frlign)


      use parm
      
      implicit none
      
      integer j, idf
      

!! ...... Argument declarations
      integer lyr   !SURF or SOIL
      real    cpart, recres(MAXIEL), &!! cdonor(ISOS), edonor(MAXIEL), 
             frlign, friso

!! ...... Partition residue from compartments cdonor and edonor
!! ...... into layer lyr of structural and metabolic.
!! ...... cpart is the amount of carbon in the residue.
!! ...... recres contains the n/c, p/c, and s/c ratios in the residue.
!! ...... frlign is the fraction of the incoming material which is lignin.
!! ...... friso is the fraction of cpart which is labeled.

!! ...... Local variables
      integer   iel 
      integer   slyr
      real      accum, caddm, cadds, dirabs(MAXIEL),            &
               eaddm, eadds, epart(MAXIEL),                     &
               fligst, frmet, frn, rcetot, rlnres, delin,       &
               dellig, c13c12r, c13frac, c13lig, c13nlig,       &
               c13struc, c12struc, c13met, c12met 
      real      namt
      real      no3nh4
      double precision frac_nh4, frac_no3
      character subname*10
      real      damr(10,3)   !damr(1,1):fraction of surface N absorbed by residue; Range: 0.0 to 1.0
                            !damr(1,2):fraction of surface P absorbed by residue; Range: 0.0 to 1.0
                            !damr(1,3):fraction of surface S absorbed by residue; Range: 0.0 to 1.0
                            !damr(2,1):fraction of soil N absorbed by residue;    Range: 0.0 to 1.0
                            !damr(2,2):fraction of soil P absorbed by residue            Range: 0.0 to 1.0
                            !damr(2,3):fraction of soil S absorbed by residue            Range: 0.0 to 1.0
      real      pabres      !amount of residue which will give maximum direct absorption of N (gC/m2)
      real      damrmn (3)  !If C/E ratio is too low, transfer just enough to make C/E of residue = damrmn
                            !damrmn(1)       minimum C/N ratio allowed in residue after direct absorption            Range: ? to ?
                            !damrmn(2)       minimum C/P ratio allowed in residue after direct absorption            Range: ? to ?
                            !damrmn(3)       minimum C/S ratio allowed in residue after direct absorption            Range: ? to ?
                            
      !!real      frac_nh4suf, frac_no3suf
      !!real      frac_nh4sol, frac_no3sol
      real      spl(2)      !spl(1)          intercept parameter for metabolic (vs. structural) split;            Range: 0.0 to 1.0
                            !spl(2)          slope parameter for metabolic split (fraction metabolic is a function of lignin to N ratio)            Range: 0.0 to 1.0
      
   !!  fake variables defined here to explain the patition process to xuesong. Real varialbes need to be find from xuesong's soc process
   
         
      real metcis  !! metbolic carbon pool
      real strcis  !! structural carbon pool
      real RCESTR(3)  !! carbon to elementy ratio in structural pool
                        !rcestr(1)       C/N ratio for structural material (fixed parameter value)            Range: ? to ?
                        !rcestr(2)       C/P ratio for structural material (fixed parameter value)            Range: ? to ?
                        !rcestr(3)       C/S ratio for structural material (fixed parameter value)            Range: ? to ?
      real struce(2,3)  !! structural element pool
      real metabe(2,3)  !! metabolic element pool
      
      j = 0
      idf = 0
      
      damr = 0.
      slyr = 2
      
      j = ihru
      idf = idplt(j)
      
      !rcestr = 1. !! give fake values to this array, SWAT should have have this variable. Need to check with Xuesong 
      RCESTR(1) = 200.0
      RCESTR(2) = 500.0
      RCESTR(3) = 500.0


      
      no3nh4 = sol_no3(lyr,j) + sol_nh3(lyr,j)  !! daycent only consider impacts of nitrogen
      
      
      
      frac_nh4 = sol_nh3(lyr,j) / no3nh4
      frac_no3 = 1 - frac_nh4

      accum = 0.0
      
      damr(1,1) = 0.00000               !damr(1,1):fraction of surface N absorbed by residue; Range: 0.0 to 1.0   
      damr(1,2) = 0.00000               !damr(1,2):fraction of surface P absorbed by residue; Range: 0.0 to 1.0    
      damr(1,3) = 0.01000               !damr(1,3):fraction of surface S absorbed by residue; Range: 0.0 to 1.0
      
      
      do 5 slyr = 2, 10
          damr(slyr,1) = 0.02000               !damr(2,1):fraction of soil N absorbed by residue;    Range: 0.0 to 1.0  
          damr(slyr,2) = 0.02000               !damr(2,2):fraction of soil P absorbed by residue            Range: 0.0 to 1.0
          damr(slyr,3) = 0.04000               !damr(2,3):fraction of soil S absorbed by residue            Range: 0.0 to 1.0
          
     
5     continue
      
      damrmn (1) = 15.00000             !damrmn(1)       minimum C/N ratio allowed in residue after direct absorption            Range: ? to ?     
      damrmn (2) = 150.00000            !damrmn(2)       minimum C/P ratio allowed in residue after direct absorption            Range: ? to ? 
      damrmn (3) = 150.00000            !damrmn(3)       minimum C/S ratio allowed in residue after direct absorption            Range: ? to ?
       
                            
                            
      spl(1) = 0.85000          !! 'SPL(1)'    
      spl(2) = 0.01300           !!'SPL(2)'     
          
      pabres = 100.             !pabres          amount of residue which will give maximum direct absorption of N (gC/m2)        Range: ? to ?
     

      if (cpart .lt. 1.e-07) then
        goto 999
      endif


!! ...... For each mineral element...
      do 10 iel = 1, 2 !! only consider N here

!! ........ Compute amount of elesent in residue.
        epart(iel) = cpart * recres(iel)        !recres contains the n/c, p/c, and s/c ratios in the residue

!! ........ Direct absorption of mineral element by residue
!! ........ (mineral will be transferred to donor compartment
!! ........ and then partitioned into structural and metabolic
!! ........ using flow routines.)

!! ........ If minerl(SRFC,iel) is negative then dirabs = zero.
        if (no3nh4.lt. 0.) then
          dirabs(iel) = 0.0                                     !Direct absorption
        else

          dirabs(iel)=damr(lyr,iel)*(0.1*no3nh4)*         &   !(0.1*no3nh4) convert kg/ha to g/m2 !damr fraction of surface and soil element absorbed by residue
                     amax1(cpart/pabres,1.)                   !pabres          amount of residue which will give maximum direct absorption of N (gC/m2)        Range: ? to ?
        
        endif

!! ........ If C/E ratio is too low, transfer just enough to make
!! ........ C/E of residue = damrmn
        if (epart(iel)+dirabs(iel) .le. 0.0) then
          rcetot = 0.0
        else
          rcetot = cpart/(epart(iel)+dirabs(iel))
        endif

        if (rcetot .lt. damrmn(iel)) then
          dirabs(iel) = cpart/damrmn(iel) - epart(iel)
        endif
        if (dirabs(iel) .lt. 0.) then
          dirabs(iel) = 0.
        endif
        
        if (iel .eq. N) then
            namt = -1.0*dirabs(iel)
            sol_no3 (lyr,j) = sol_no3 (lyr,j) - 10.*namt * frac_no3 
            sol_nh3 (lyr,j) = sol_nh3 (lyr,j) - 10.*namt * frac_nh4     
        endif
      
10    continue

!! ...... Partition carbon into structural and metabolic fraction of
!! ...... residue (including direct absorption) which is nitrogen
      frn = (epart(1)+dirabs(1)) / (cpart*2.5)  !g/m2

!! ...... Lignin/nitrogen ratio of residue
      rlnres = frlign/frn

!! ...... Carbon added to metabolic
!! ...... Compute the fraction of carbon that goes to metabolic.
      frmet = spl(INTCPT)-spl(SLOPE)*rlnres

!! ...... Make sure the fraction of residue which is lignin isn't
!! ...... greater than the fraction which goes to structural.  -rm 12/91
      if (frlign .gt. (1.0 - frmet)) then
        frmet = (1.0 - frlign)
      endif

!! ...... Make sure at least 20% goes to metabolic
      if (frmet .lt. 0.20) then
        frmet = .20
      endif

!! ...... Compute amounts to flow
      caddm = cpart * frmet !! carbon added to metabolic pool
      if (caddm .lt. 0) then
        caddm = 0.0
      endif
      cadds = cpart-caddm   !! carbon added to structural pool
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX STOP here because xuesong may have done the particition
!! ...... Adjust lignin content of structural.
!! ...... fligst is the fraction of incoming structural residue
!! ...... which is lignin; restricting it to a maximum of .8
      fligst = frlign/(cadds/cpart)

!! ...... Changed allowable maximum from .8 to .6 -rm 5/92
!! ...... Changed maximum fraction from .6 to 1.0  -lh 1/93
      if (fligst .gt. 1.0) then
         fligst = 1.0
      endif

! ... Adjust lignin
!!XXXXX      call adjlig(strcis,fligst,cadds,strlig(lyr))   !! need to check with Xuesong if this process is needed, the key is if SWAT has a lignin fraction for the structural pool


            j = ihru
            sol_LS(lyr,ihru) = sol_LS(lyr,ihru) + cadds * 10 /0.42            !*10 changes g/m2 to kg/ha
            sol_LM(lyr,ihru) = sol_LM(lyr,ihru) + caddm * 10 /0.42  
            sol_LSC(lyr,ihru) = sol_LSC(lyr,ihru) + cadds * 10                ! *10 Carbon added to structural !factor 10 changes g/m2 to kg/ha
            sol_LSL(lyr,ihru) = sol_LSL(lyr,ihru) + cadds * 10 *fligst
            sol_LMC(lyr,ihru) = sol_LMC(lyr,ihru) + caddm * 10                !Carbon added to metabolic
            
        !!  sol_HSC(:,:)    : mass of C present in slow humus (kg ha-1)
        !!  sol_HSN(:,:)    : mass of N present in slow humus (kg ha-1)
        !!  sol_HPC(:,:)    : mass of C present in passive humus (kg ha-1)
        !!  sol_HPN(:,:)    : mass of N present in passive humus (kg ha-1)
        !!  sol_LM(:,:)     : mass of metabolic litter (kg ha-1)
        !!  sol_LMC(:,:     : mass of C in metabolic litter (kg ha-1)
        !!  sol_LMN(:,:)    : mass of N in metabolic litter (kg ha-1)
        !!  sol_LS(:,:)     : mass of structural litter (kg ha-1)
        !!  sol_LSC(:,:)    : mass of C in structural litter (kg ha-1)
        !!  sol_LSL(:,:)    : mass of lignin in structural litter (kg ha-1)
        !!  sol_LSN(:,:)    : mass of N in structural litter (kg ha-1)

  
       

!! ...... Adjust lignin
    !1  call adjlig(strucc(lyr),fligst,cadds,strlig(lyr))

!! ...... Partition mineral elements into structural and metabolic
      do 20 iel = 1, nelem
    !! ........ Flow into structural
            !!!!eadds = cadds/rcestr(iel)
           !! call flow(edonor(iel),struce(lyr,iel),time,eadds)
           !!  edonor(iel) = edonor(iel) - eadds
             !!!!struce(lyr,iel) = struce(lyr,iel) + eadds
    !! ........ Flow into metabolic
            !!!!eaddm = epart(iel)+dirabs(iel)-eadds
           !! call flow(edonor(iel),metabe(lyr,iel),time,eaddm)
            !! edonor(iel) = edonor(iel) - eaddm
            !!!!metabe(lyr,iel) = metabe(lyr,iel) + eaddm
        if (iel == 1) then
            eadds = cadds/rcestr(iel)
            sol_LSN(lyr,j) = sol_LSN(lyr,j) + eadds * 10 !factor 10 changes g/m2 to kg/ha
            eaddm = epart(iel)+dirabs(iel)-eadds
            sol_LMN(lyr,j) = sol_LMN(lyr,j) + eaddm  * 10 !factor 10 changes g/m2 to kg/ha
        end if
        if (iel == 2) then
            sol_fop(lyr,j) = sol_fop(lyr,j) + cpart/rcestr(iel) * 10 !factor 10 changes g/m2 to kg/ha
        end if        
 20    continue

       
      sol_rsd(lyr,j)= sol_LS(lyr,j)+sol_LM(lyr,j)            
      sol_fon(lyr,j) = sol_LMN(lyr,j) + sol_LSN(lyr,j) 



      !if cpart is too small, skip
 999   continue

      return
      end
