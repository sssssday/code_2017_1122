      integer function get_stage (tday,tjth)
        implicit none
        integer, intent(in):: tday, tjth
        integer :: ytday
        real :: tday_phe, ytday_phe
        real :: stage
        real :: cal_phen
        !!real :: get_stage

        if(tday==0) then
           ytday = 364
        else
           ytday = tday-1
        endif

        tday_phe= cal_phen(tday,tjth)
        ytday_phe = cal_phen(ytday, tjth)

        if(tday_phe == ytday_phe) then
           stage = 0
        else if(tday_phe > ytday_phe) then
          stage = 1

        else
          stage = 2

        end if

        get_stage = stage
        
        return

        end function
