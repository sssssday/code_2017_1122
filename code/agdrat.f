


      real function agdrat (aminrl,varat,iel)

      use parm
      implicit none

!! ......Argument declarations
      integer    iel
      real       aminrl(3), varat(3,3)

!! ......AboveGround Decomposition RATio computation.
!! ......Determine C/E of new material entering 'Box B'.
!! ......'E' represents N, P, or S.

!! ......The C/E ratios for structural and wood can be computed once; 
!! ......they then remain fixed throughout the run.  The ratios for 
!! ......metabolic and som may vary and must be recomputed each time
!! ......step.

!! ......Input:
!! ......  aminrl = mineral E in top soil layer
!! ......  varat : varat(1,iel) = maximum C/E ratio for material
!! ......                         entering 'Box B'
!! ......          varat(2,iel) = minimum C/E ratio for material
!! ......          varat(3,iel) = amount of E present when minimum
!! ......                         ratio applies
!! ......  iel  = index to element for which ratio is being computed;
!! ......         1 for N, 2 for P, 3 for S

!! ......Output:
!! ......  agdrat = C/E ratio of new material where E is N, P, or S
!! ......           depending on the value of iel

!! ......Ratio depends on available E 

      if (aminrl(iel) .le. 0.0) then

c ..... Set ratio to maximum allowed
        agdrat = varat(1,iel)

      elseif (aminrl(iel) .gt. varat(3,iel)) then

c ..... Set ratio to minimum allowed
        agdrat = varat(2,iel)

      else

c ..... aminrl(iel) > 0 and <= varat(3,iel)
        agdrat = (1.-aminrl(iel)/varat(3,iel))*
     &           (varat(1,iel)-varat(2,iel))+varat(2,iel)
c ..... where: varat(1,iel) = maximum C/E ratio for material
c .....                       entering 'Box B'
c .....        varat(2,iel) = minimum C/E ratio for material
c .....        varat(3,iel) = amount of E present when minimum
c .....                       ratio applies

      endif

      return
      end