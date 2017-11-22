
      real function line(x, x1, y1, x2, y2)
      use parm

      implicit none

!! ... Argument declarations
      real:: x, x1, y1, x2, y2

!! ... This function is the generic equation of a line from
!! ... two points.
!!  ... Given 2 known points and a new X, calculate Y.
!!  ... slope = (y2 - y1) / (x2 - x1)
!!  ... y = slope * (x - x2) + y2
!!  ...


      line = (y2 - y1) / (x2 - x1) * (x - x2) + y2

      return
      end
