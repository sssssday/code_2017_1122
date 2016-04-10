



      real function bgdrat (aminrl,varat,iel)

      use parm
      implicit none

!! ...... Argument declarations
      integer   iel
      real      aminrl(3), varat(3,3)

!! ...... BelowGround Decomposition RATio computation.

!! ...... Determine C/E of new material entering 'Box B'.

!! ...... Input:
!! ......   aminrl = mineral E in top soil layer
!! ......   varat : varat(1,iel) = maximum C/E ratio for material
!! ......                          entering 'Box B'
!! ......           varat(2,iel) = minimum C/E ratio for material
!! ......           varat(3,iel) = amount of E present when minimum
!! ......                          ratio applies
!! ......   iel = index to element for which ratio is being computed;
!! ......         1 for N, 2 for P, 3 for S

!! ...... Output:
!! ......   bgdrat = C/E ratio of new material where E is N, P, or S
!! ......            depending on the value of iel

!! ...... Ratio depends on available E 

      if (aminrl(iel) .le. 0) then

!! ........ Set ratio to maximum allowed
        bgdrat = varat(1,iel)

      elseif (aminrl(iel) .gt. varat(3,iel)) then

!! ........ Set ratio to minimum allowed
        bgdrat = varat(2,iel)

      else

!! ........ aminrl(iel) > 0 and <= varat(3,iel)
        bgdrat = (1.-aminrl(iel)/varat(3,iel))*
     &           (varat(1,iel)-varat(2,iel))+varat(2,iel)
!! ........ where: varat(1,iel) = maximum C/E ratio for material
!! ........                       entering 'Box B'
!! ........        varat(2,iel) = minimum C/E ratio for material
!! ........        varat(3,iel) = amount of E present when minimum
!! ........                       ratio applies

      endif

      return
      end
