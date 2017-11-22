      real function daylenth(dayth,j)
!! qichun this function has been modified to calculate day legth of the previous day      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates the daylength, distribution of 
!!    radiation throughout the day and maximum radiation for day

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    iida        |julian date   |day being simulated (current julian date)
!!    ipet        |none          |code for potential ET method
!!                               |0 Priestley-Taylor method
!!                               |1 Penman/Monteith method
!!                               |2 Hargreaves method
!!    j           |none          |HRU number
!!    latcos(:)   |none          |Cos(Latitude) for HRU
!!    latsin(:)   |none          |Sin(Latitude) for HRU
!!    npcp(:)     |none          |prior day category
!!                               |1 dry day
!!                               |2 wet day
!!    subp(:)     |mm H2O        |precipitation for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dayl(:)     |hours         |day length
!!    frad(:,:)   |none          |fraction of solar radiation occuring during 
!!                               |hour in day in HRU
!!    hru_rmx(:)  |MJ/m^2        |maximum possible radiation for the day in HRU
!!    npcp(:)     |none          |prior day category
!!                               |1 dry day
!!                               |2 wet day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ch          |none          |(-Tan(sd)*Tan(lat))
!!    cosrho(:)   |none          |Cos(zenith angle for hour)
!!    dd          |none          |relative distance of the earth from the sun
!!    h           |none          |Acos(-Tan(sd)*Tan(lat))
!!    ii          |none          |counter
!!    sd          |radians       |solar declination: latitude at which the sun
!!                               |is directly overhead at noon
!!    totrho      |none          |sum of cosrho values for all hours of day
!!    yc          |none          |Cos(sd)*Cos(lat)
!!    ys          |none          |Sin(sd)*Sin(lat)
!!    w           |none          |hour angle
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sin, Tan, Acos, Cos, Real

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer, intent (in) :: j, dayth
      integer :: ii
      real :: sd, ch, h, ys, yc, dd, w, cosrho(nstep), totrho


      

      !! Calculate Daylength !!
      !! calculate solar declination: equation 2.1.2 in SWAT manual
      sd = 0.
      sd = Asin(.4 * Sin((Real(dayth) - 82.) / 58.09))  !!365/2pi = 58.09

      !! calculate the relative distance of the earth from the sun
      !! the eccentricity of the orbit
      !! equation 2.1.1 in SWAT manual
      dd = 0.
      dd = 1.0 + 0.033 * Cos(Real(dayth) / 58.09)

      !!daylength = 2 * Acos(-Tan(sd) * Tan(lat)) / omega
      !!where the angular velocity of the earth's rotation, omega, is equal
      !! to 15 deg/hr or 0.2618 rad/hr and 2/0.2618 = 7.6374
      !! equation 2.1.6 in SWAT manual

      ch = 0.
      h = 0.
      ch = -latsin(hru_sub(j)) * Tan(sd) / latcos(hru_sub(j))
      if (ch > 1.) then    !! ch will be >= 1. if latitude exceeds +/- 66.5 deg in winter
        h = 0.
      elseif (ch >= -1.) then
        h = Acos(ch)
      else
        h = 3.1416         !! latitude exceeds +/- 66.5 deg in summer
      endif
      daylenth = 7.6394 * h
          
      !! Calculate Potential (maximum) Radiation !!
      !! equation 2.2.7 in SWAT manual

      

      return
      end function