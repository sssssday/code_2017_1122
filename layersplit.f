      subroutine layersplit(dep_new)

      use parm
      !!implicit none
      integer nly,ltn,j  !! ltn is changed from n in the original verison by qichun
	integer :: flag
	real, intent(in):: dep_new
	nly = sol_nly(ihru)

!!    create a septic layer
!! changed all sol_zmx(ihru) in subroutine to dep_new 1/27/09 gsm 
      flag = 0
      
      do j = 2, nly 
        xx = 0.
        xx = Abs(dep_new - sol_z(j,ihru))
        !! if values are within 10 mm of one another, reset boundary
        if (xx < 10.) then
          sol_z(j,ihru) = dep_new
          exit
        end if

        !! set a soil layer at dep_new and adjust all lower layers
        if (sol_z(j,ihru) > dep_new) then
          flag = 1
          sol_nly(ihru) = sol_nly(ihru) + 1
          nly = nly + 1
          jj = 0
          jj = j + 1
          do ltn = nly, jj, -1
            sol_z(ltn,ihru) = sol_z(ltn-1,ihru)
            sol_bd(ltn,ihru) = sol_bd(ltn-1,ihru)
            sol_awc(ltn,ihru) = sol_awc(ltn-1,ihru)
            sol_k(ltn,ihru) = sol_k(ltn-1,ihru)
            sol_cbn(ltn,ihru) = sol_cbn(ltn-1,ihru)
            sol_rock(ltn,ihru) = sol_rock(ltn-1,ihru) !!! Armen 13 Jan 2008 MJW rev 490
            sol_clay(ltn,ihru) = sol_clay(ltn-1,ihru)
            sol_sand(ltn,ihru) = sol_sand(ltn-1,ihru) !!! Claire 2 Dec 2009 MJW rev 490
            sol_silt(ltn,ihru) = sol_silt(ltn-1,ihru) !!! Claire 2 Dec 2009 MJW rev 490
            sol_ec(ltn,ihru) = sol_ec(ltn-1,ihru)
            sol_no3(ltn,ihru) = sol_no3(ltn-1,ihru)
            sol_orgn(ltn,ihru) = sol_orgn(ltn-1,ihru)
            sol_orgp(ltn,ihru) = sol_orgp(ltn-1,ihru)
            sol_solp(ltn,ihru) = sol_solp(ltn-1,ihru)
            sol_mc(ltn,ihru) = sol_mc(ltn-1,ihru)
            sol_mn(ltn,ihru) = sol_mn(ltn-1,ihru)
            sol_mp(ltn,ihru) = sol_mp(ltn-1,ihru)
		    sol_n(ltn,ihru) = sol_n(ltn-1,ihru)
		    sol_ph(ltn,ihru) = sol_ph(ltn-1,ihru) !! MJW rev 490
		    sol_cal(ltn,ihru) = sol_cal(ltn-1,ihru) !! MJW rev 490

          end do
          sol_z(j,ihru) = dep_new
        end if
        if (flag == 1) exit
      end do
	
	iseptic = j 
      end             