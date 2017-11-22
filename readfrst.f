      subroutine readfrst

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine reads input parameters from the forest
!!    database (forest.dat). Added by Qichun Yang 2015.8.12

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name      |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  ~ ~ ~ ~
!!    mcrdb     |none             |maximum number of crops/landcover in 
!!                                |database file (crop.dat)
!!    rsdco     |none             |residue decomposition coefficient
!!                                |The fraction of residue which will decompose
!!                                |in a day assuming optimal moisture,
!!                                |temperature, C:N ratio, and C:P ratio
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name       |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    alai_min(:)|m**2/m**2        |minimum LAI during winter dormant period
!!    bio_e(:)   |(kg/ha)/(MJ/m**2)|biomass-energy ratio
!!                                 |The potential (unstressed) growth rate per
!!                                 |unit of intercepted photosynthetically
!!                                 |active radiation.
!!    bio_leaf(:)|none             |fraction of leaf/needle biomass that drops 
!!                                 |during dormancy (for trees only)
!!    bio_n1(:)  |none             |1st shape parameter for plant N uptake 
!!                                 |equation
!!    bio_n2(:)  |none             |2nd shape parameter for plant N uptake 
!!                                 |equation
!!    bio_p1(:)  |none             |1st shape parameter for plant P uptake 
!!                                 |equation
!!    bio_p2(:)  |none             |2st shape parameter for plant P uptake 
!!                                 |equation
!!    blai(:)    |none             |maximum (potential) leaf area index
!!    bm_dieoff(:) fraction        |fraction above ground biomass that dies
!!                                 |   off at dormancy
!!    chtmx(:)   |m                |maximum canopy height
!!    cnyld(:)   |kg N/kg yield    |fraction of nitrogen in yield
!!    cpnm(:)    |NA               |four character code to represent crop name
!!    cpyld(:)   |kg P/kg yield    |fraction of phosphorus in yield
!!    cvm(:)     |none             |natural log of USLE_C
!!    dlai(:)    |none             |fraction of growing season when leaf
!!                                 |area declines
!!    gsi(:)     |m/s              |maximum stomatal conductance
!!    hvsti(:)   |(kg/ha)/(kg/ha)  |harvest index: crop yield/aboveground 
!!                                 |biomass
!!    idc(:)     |none             |crop/landcover category:
!!               |                 |1 warm season annual legume
!!               |                 |2 cold season annual legume
!!               |                 |3 perennial legume
!!               |                 |4 warm season annual
!!               |                 |5 cold season annual
!!               |                 |6 perennial
!!               |                 |7 trees
!!    leaf1(:)   |none             |1st shape parameter for leaf area
!!                                 |development equation.
!!    leaf2(:)   |none             |2nd shape parameter for leaf area
!!                                 |development equation.
!!    pltnfr(1,:)|kg N/kg biomass  |nitrogen uptake parameter #1: normal
!!                                 |fraction of N in crop biomass at emergence
!!    pltnfr(2,:)|kg N/kg biomass  |nitrogen uptake parameter #2: normal
!!                                 |fraction of N in crop biomass at 0.5 
!!                                 |maturity 
!!    pltnfr(3,:)|kg N/kg biomass  |nitrogen uptake parameter #3: normal
!!                                 |fraction of N in crop biomass at maturity
!!    pltpfr(1,:)|kg P/kg biomass  |phosphorus uptake parameter #1: normal 
!!                                 |fraction of P in crop biomass at emergence
!!    pltpfr(2,:)|kg P/kg biomass  |phosphorus uptake parameter #2: normal
!!                                 |fraction of P in crop biomass at 0.5 
!!                                 |maturity
!!    pltpfr(3,:)|kg P/kg biomass  |phosphorus uptake parameter #3: normal
!!                                 |fraction of P in crop biomass at maturity
!!    rdmx(:)    |m                |maximum root depth
!!    rsdco_pl(:)|none             |plant residue decomposition coefficient. The
!!                                 |fraction of residue which will decompose in
!!                                 |a day assuming optimal moisture,
!!                                 |temperature, C:N ratio, and C:P ratio
!!    rsr1c      |                 |initial root to shoot ratio at the beg of growing season
!!    rsr2c      |                 |root to shoot ratio at the end of the growing season
!!    t_base(:)  |deg C            |minimum temperature for plant growth
!!    t_opt(:)   |deg C            |optimal temperature for plant growth
!!    vpd2(:)    |(m/s)*(1/kPa)    |rate of decline in stomatal conductance per
!!                                 |unit increase in vapor pressure deficit 
!!    wac21(:)   |none             |1st shape parameter for radiation use
!!                                 |efficiency equation.
!!    wac22(:)   |none             |2nd shape parameter for radiation use
!!                                 |efficiency equation.
!!    wavp(:)    |none             |Rate of decline in radiation use efficiency
!!                                 |as a function of vapor pressure deficit
!!    wsyf(:)    |(kg/ha)/(kg/ha)  |Value of harvest index between 0 and HVSTI 
!!                                 |which represents the lowest value expected 
!!                                 |due to water stress
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name      |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    alaimin   |m**2/m**2        |minimum leaf area index for winter dormant 
!!                                |period
!!    b1        |none             |variable to hold calculation results
!!    b2        |none             |variable to hold calculation results
!!    b3        |none             |variable to hold calculation results
!!    bioehi    |(kg/ha)/(MJ/m**2)|biomass-energy ratio when plant is in
!!                                |an environment with CO2 level equal to
!!                                |the value of CO2HI. This biomass-energy
!!                                |ratio is used to set the 2nd point on the
!!                                |radiation use efficiency curve
!!    bioleaf   |none             |fraction of biomass accumulated each year
!!                                |that is leaf/needle
!!    c1        |none             |variable to hold calculation results
!!    co2hi     |uL CO2/L air     |CO2 concetration higher than the ambient
!!                                |corresponding to the 2nd point on radiation
!!                                |use efficiency curve
!!    frgmax    |none             |fraction of maximum stomatal conductance
!!                                |that is achieved at the vapor pressure
!!                                |deficit defined by VPDFR
!!    frgrw1    |none             |fraction of the growing season corresponding
!!                                |to the 1st point on optimal leaf area 
!!                                |development curve
!!    frgrw2    |none             |fraction of the growing season corresponding
!!                                |to the 2nd point on optimal leaf area 
!!                                |development curve
!!    eof       |none             |end of file flag (=-1 of eof, else =0)
!!    ic        |none             |landuse/landcover array storage number
!!                                |when a land cover is assigned in the 
!!                                |.mgt file, the variables for the land
!!                                |cover are accessed by the array number.
!!                                |Landuse/landcover numbers (ICNUM) in 
!!                                |crop.dat need to be assigned consecutively
!!                                |to ensure that the crop number used by the
!!                                |user is the same as the array storage number
!!    icnum     |none             |crop/landcover number. Reference number only.
!!    laimx1    |none             |fraction of maximum leaf area index 
!!                                |corresponding to the 1st point on optimal 
!!                                |leaf area development curve
!!    laimx2    |none             |fraction of maximum leaf area index 
!!                                |corresponding to the 2nd point on optimal 
!!                                |leaf area development curve
!!    usle_c    |none             |minimum value of the USLE C factor for water
!!                                |erosion
!!    vpdfr     |kPa              |vapor pressure deficit at which FRGMAX is
!!                                |valid
!!    xx        |none             |dummy variable to hold IDC expressed as a
!!                                |real number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    NInt, Int, Log, ascrv 

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
!!     added by Qichun Yang 2015.8.12
!!     part1 physiology parameter
!!        float salb;                   //short wave albedo
!!        float slalb;                 // soil albedo
!!        float slfg;                  // ratio of ground heat flux to ground
!!        float slrb;                  // soil boundary resistance
!!
!!
!!        float omga;                    //(index) foliage clumping index
!!        float o3af;                  //ozone effect on photosynthesis
!!        float hvtp;                //hearvest part   0-1    Yield/biomass
!!
!!        float inisn;        // initial soil mineral N (NH4 + NO3) content (gN/m2); added by ZC 6/20/06
!!        float smaxnh4;          //maximum soil NH4 gN/m2
!!        float smaxno3;       //maximum soil NO3
!!
!!        float inic;                    //seedC
!!        float inin;                    //seedN
!!
!!        float kltow0;                  //a function of the leaf biomass to the live wood biomass
!!        float kltow1;                  //a function of the leaf biomass to the live wood biomass
!!                                        //y = k0 * exp(k1*Csap)
!!        float kw_maxcarbon;             //for control the max lai
!!
!!        float maxlai;                  //maximum lai
!!        float rinilfc;                 //fractional of potmaximum leafc to start
!!                                        //
!!        int maxage;                     //maximum age for woody plant. added by CGS 4/22/2008
!!
!!        float alloc_root0;               //initial allocation rate to root
!!        float alloc_stem0;               //initial allocation rate to stem
!!        float adjStem_to_Root;
!!        float alloc_coarseroot;           //percentage of allocated root carbon to coarseroot
!!        float alloc_fineroot;             //percentage of allocated root carbon to fineroot
!!        float alloc_reproduction;         //reporduction part from allocatable npp   0-1
!!        float alloc_K_leaf;               //allocation coefficient to leaf as L = exp(-K*LAI)

!!        float roota;            //parameters for root vertical distributions  Zeng et al., 2001
!!        float rootb;

!!        float max_Nup;          //maximum n uptake speed gN/day
!!        float knuphalf;         //the half-saturation coefficient for nitrogen uptake (gN/m2)
!!       float nup_tmin;        // = -5.0; // minmum N uptake temperature; unit: 0C
!!        float nup_tmax;        // = 30.0; // max N uptake temperature; unit: 0C
!!        float nup_topt;        // = 15.0; // reference temperature at which nupTempFact = 1; unit: 0C
!!        float knroot;           // parameter controlling the effect of fineroot on nitrogen uptake
!!
!!        float anfix;            //gN/m2/yr annual nitrogen fixation rate
!!        float buffer_NO3;       // Dissolved NO3 buffer capacity in soils
!!        float buffer_NH4;       // Dissolved NH4 buffer capacity in soils
!!        float buffer_DIC;       // Dissolved DIC buffer capacity in soils
!!        float buffer_POM;       // Particulate organic matter buffer capacity in soils
!!
!!        float den_kcn;          //coeffficient to estimate pot denitrification rate
!!        float maxnit;           //the potential nitrification rate under optimal conditions
!!        float kgasnitri;      //proportion of N intermediates resulting in N2O emissions from nitrification
!!        float nh3VolitCoeff; // NH3 volatization coefficient
!!        float knsorb_half;
!!        float knsorb_max; // coefficients to calculate the SIN available for loss (avn4Loss).
!!        float knleak;  // (0 to 1.0) DON & DIN loss rate (except for leaching)
!!
!!        //float respcoeff;       //coefficient in respiration equations gC/gN/day
!!        float rq10;            //respiration Q10
!!        float mr_leaf;         //(gC/gC/day) mr coefficient for leaf
!!        float mr_stem;         //(gC/gC/day) mr coefficient for stem
!!        float mr_csroot;        //(gC/gC/day) mr coefficient for coarse root
!!        float mr_froot;        //(gC/gC/day) mr coefficient for fineroot
!!        float mr_rep;          //(gC/gC/day) mr coefficient for reproducts

!!        //float avg_proj_sla;             //(m2/gC) canopy average proj. SLA
!!        float avg_sla;                  //(m2/gc) specific leaf area (SLA)
!!        float lai_ratio;                //(DIM) ratio of (all-sided LA / one-sided LA)
!!        float ext_coef;                 //(DIM) canopy light extinction coefficient
!!        float ext_n;                    //(0.5) nitrogen distribution coe�cient
!!        float max_canopyw_coef;         //(mmH2O/mm PPT/LAI) water interception coefficient
!!                                         //LAI: all-side LAI
!!        float max_canopyw_daymm;        //(mmH2O) maximum interception per day on whole plant

!!        float fcn_leafmin;      // minimum leaf C:N ratio
!!        float fcn_heartwood;
!!        float fcn_sapwood;
!!        float fcn_coarseroot;
!!        float fcn_fineroot;
!!        float fcn_reproduction;

!!        float leaflitterCNratio;        //
!!        float finerootlitterCNratio;    //
!!        float ftcsroot;                 //(yr) coarse root turn over rate
!!        float ftfroot;
!!        float ftleaf;                  //(yr) leaf turnover rate in growing season
!!        float ftsapwood;               //(yr) sapwood or branch turnover rate
!!        float ftprod;
!!        float ftsapwood_heartwood;      //fraction of turnovered sapwood to heart wood, for nonwoody biomes, it is zero

!!        float fharvest;           //part of residues after crop harvesting to be removed from ecosystem (through fire etc.)!!        int irrigation_frequency;       //(days) to irrigate
!!        float proddecay_Nreturn;  //part of nitrogen of decaded PROD1 to SOM (labile) for crop ecosystems

!!        //methane parameters
!!        float CH4ProMaxRate;            //maximum rate of ch4 production gC/m3/day
!!        float CH4OxidairMaxRate;        //maximum rate of air ch4 oxidization gC/m3/day
!!        float CH4OxidsoilMaxRate;        //maximum rate of soil ch4 oxidization gC/m3/day
!!        float CH4OxidsoilMaxRateSatw;    //maximum rate of soil ch4 oxidization gC/m3/day under saterated water
!!        float Kmch4pro;                 //concentration gC/m3
!!        float Kmch4airMaxrate;          //ppm
!!        float Kmch4soilMaxrate;         //gC/m3
!!        float Kmch4soilMaxrateSatw;     //gC/m3 under conditions of saterated water

!!        //fire parameters
!!        float fireEmission_CO2;        //gCO2 /kg dry matter
!!        float fireEmission_CO;         //gCO / kg dry matter
!!        float fireEmission_CH4;        //gCH4 /kg dry matter
!!        float fireEmission_NMHC;       //gNMHC /kg dry matter

!!        float fireEmission_NOy;        //gNOy/kg d.m. total nitrogen oxides
!!        float fireEmission_NH3;        //gNH3/kg d.m.
!!        float fireEmission_N2O;        //gN2O/kg d.m.

!!        float CC_leaf;      //%  combustion completeness of leaf
!!        float CC_stem;      //%  combustion completeness of stem
!!        float CC_root;      //%  combustion completeness of root
!!        float CC_rep;       //%  combustion completeness of reproduct organ
!!        float CC_litter;    //%  combustion completeness of litter
!!        float CC_CWD;       //%  combustion completeness of CWD

!!        float mort_leaf;        // postfire mortality of leaf
!!        float mort_stem;        // postfire mortality of stem
!!        float mort_root;        // postfire mortality of root
!!        float mort_rep;         // postfire mortality of rep



!!        float moist_extinct;

!!        float usle_cmin;               //minimum value of the USLE C factor for water
!!                                       //erosion

!!        // parameters for VOC emission
!!       float Eps_isop;                 // (ug C/g dry foliar biomass/day) emission capacity for isoprene
!!        float Eps_mono;                 // (ug C/g dry foliar biomass/day) emission capacity for monoterpens
!!        float Eps_ovoc;                 // (ug C/g dry foliar biomass/day) emission capacity for other VOC
!!        float Eps_orvoc;                // (ug C/g dry foliar biomass/day) emission capacity for other reactive VOC
!!        float Eps_co;
!!


      use parm
      implicit none
      integer :: ic, eof, icnum, yrsmat, idtype
      character (len=100):: paraname
   !! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!     part1 photosynthesis variables
      real :: salbd, slalb, slfg, omga 
      real :: inicd, avg_sla, ext_n,
     &        fcn_leafmin,
     &        ra, rb,
     &        ko25, kc25, vcmx25, ako, akc, avcmx, vpd_open, vpd_close, 
     &        toptd,
     &        gmax, gmin

      real :: chtmax, tbase, matyrs
      
!!      part2 variables from daycent

        integer :: DECID
        real :: PRDX2
        real :: PPDF1
        real :: PPDF2
        real :: PPDF3
        real :: PPDF4
        real :: CERFOR111
        real :: CERFOR112
        real :: CERFOR113
        real :: CERFOR121
        real :: CERFOR122
        real :: CERFOR123
        real :: CERFOR131
        real :: CERFOR132
        real :: CERFOR133
        real :: CERFOR141
        real :: CERFOR142
        real :: CERFOR143
        real :: CERFOR151
        real :: CERFOR152
        real :: CERFOR153
        real :: CERFOR211
        real :: CERFOR212
        real :: CERFOR213
        real :: CERFOR221
        real :: CERFOR222
        real :: CERFOR223
        real :: CERFOR231
        real :: CERFOR232
        real :: CERFOR233
        real :: CERFOR241
        real :: CERFOR242
        real :: CERFOR243
        real :: CERFOR251
        real :: CERFOR252
        real :: CERFOR253
        real :: CERFOR311
        real :: CERFOR312
        real :: CERFOR313
        real :: CERFOR321
        real :: CERFOR322
        real :: CERFOR323
        real :: CERFOR331
        real :: CERFOR332
        real :: CERFOR333
        real :: CERFOR341
        real :: CERFOR342
        real :: CERFOR343
        real :: CERFOR351
        real :: CERFOR352
        real :: CERFOR353
        real :: DECW1
        real :: DECW2
        real :: DECW3
        real :: FCFRAC11
        real :: FCFRAC21
        real :: FCFRAC31
        real :: FCFRAC41
        real :: FCFRAC51
        real :: FCFRAC12
        real :: FCFRAC22
        real :: FCFRAC32
        real :: FCFRAC42
        real :: FCFRAC52
        real :: TFRTCN1
        real :: TFRTCN2
        real :: TFRTCW1
        real :: TFRTCW2
        real :: LEAFDR1
        real :: LEAFDR2
        real :: LEAFDR3
        real :: LEAFDR4
        real :: LEAFDR5
        real :: LEAFDR6
        real :: LEAFDR7
        real :: LEAFDR8
        real :: LEAFDR9
        real :: LEAFDR10
        real :: LEAFDR11
        real :: LEAFDR12
        real :: BTOLAI
        real :: KLAI
        real :: LAITOP
        real :: MAXLAI
        real :: MAXLDR
        real :: FORRTF1
        real :: FORRTF2
        real :: FORRTF3
        real :: SAPK
        real :: SWOLD
        real :: WDLIG1
        real :: WDLIG2
        real :: WDLIG3
        real :: WDLIG4
        real :: WDLIG5
        real :: WDLIG6
        real :: WOODDR1
        real :: WOODDR2
        real :: WOODDR3
        real :: WOODDR4
        real :: WOODDR5
        real :: WOODDR6
        real :: WRDSRFC
        real :: WMRTFRAC
        real :: SNFXMX2
        real :: DEL13C
        real :: CO2IPR2
        real :: CO2ITR2
        real :: CO2ICE211
        real :: CO2ICE212
        real :: CO2ICE213
        real :: CO2ICE221
        real :: CO2ICE222
        real :: CO2ICE223
        real :: CO2IRS2
        real :: BASFC2
        real :: BASFCT
        real :: SITPOT
        real :: MAXNP
        real :: FKMRSPMX1
        real :: FKMRSPMX2
        real :: FKMRSPMX3
        real :: FKMRSPMX4
        real :: FKMRSPMX5
        real :: FKMRSPMX6
        real :: FMRSPLAI1
        real :: FMRSPLAI2
        real :: FMRSPLAI3
        real :: FMRSPLAI4
        real :: FMRSPLAI5
        real :: FMRSPLAI6
        real :: FGRESP1
        real :: FGRESP2
        real :: FGRESP3
        real :: FGRESP4
        real :: FGRESP5
        real :: FGRESP6
        real :: NO3PREF2
        integer :: TLAYPG
        real :: TMIX
        real :: TMPLFF
        real :: TMPLFS
        integer :: FURGDYS
        real :: FLSGRES
        real :: TMXTURN
        real :: WSCOEFF21
        real :: WSCOEFF22
        real :: NPP2CS2
        real :: SFAVAIL2

   

 !! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      character (len=4) :: cname

      eof = 0
      read (333,5333,iostat=eof) paraname



      do
!!
!!      initialize locals in loop
!! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    part1 physiology parameter

             salbd = 0.0
             slalb = 0.0 
             slfg = 0.0 
             omga = 0.0 
             inicd = 0.0
             avg_sla = 0.0 
             ext_n = 0.0
             fcn_leafmin = 0.0
             ra = 0.0 
             rb = 0.0
             ko25 = 0.0 
             kc25 = 0.0 
             vcmx25 = 0.0 
             ako = 0.0 
             akc = 0.0 
             avcmx = 0.0 
             vpd_open = 0.0 
             vpd_close = 0.0 
             toptd = 0.0
             gmax = 0.0 
             gmin = 0.0
             chtmax = 0.0 
             tbase = 0.0
             matyrs=0.0
!!  part2 variables from daycent

        DECID = 0 
        PRDX2 = 0.0
        PPDF1 = 0.0
        PPDF2 = 0.0
        PPDF3 = 0.0
        PPDF4 = 0.0
        CERFOR111 = 0.0
        CERFOR112 = 0.0
        CERFOR113 = 0.0
        CERFOR121 = 0.0
        CERFOR122 = 0.0
        CERFOR123 = 0.0
        CERFOR131 = 0.0
        CERFOR132 = 0.0
        CERFOR133 = 0.0
        CERFOR141 = 0.0
        CERFOR142 = 0.0
        CERFOR143 = 0.0
        CERFOR151 = 0.0
        CERFOR152 = 0.0
        CERFOR153 = 0.0
        CERFOR211 = 0.0
        CERFOR212 = 0.0
        CERFOR213 = 0.0
        CERFOR221 = 0.0
        CERFOR222 = 0.0
        CERFOR223 = 0.0
        CERFOR231 = 0.0
        CERFOR232 = 0.0
        CERFOR233 = 0.0
        CERFOR241 = 0.0
        CERFOR242 = 0.0
        CERFOR243 = 0.0
        CERFOR251 = 0.0
        CERFOR252 = 0.0
        CERFOR253 = 0.0
        CERFOR311 = 0.0
        CERFOR312 = 0.0
        CERFOR313 = 0.0
        CERFOR321 = 0.0
        CERFOR322 = 0.0
        CERFOR323 = 0.0
        CERFOR331 = 0.0
        CERFOR332 = 0.0
        CERFOR333 = 0.0
        CERFOR341 = 0.0
        CERFOR342 = 0.0
        CERFOR343 = 0.0
        CERFOR351 = 0.0
        CERFOR352 = 0.0
        CERFOR353 = 0.0
        DECW1 = 0.0
        DECW2 = 0.0
        DECW3 = 0.0
        FCFRAC11 = 0.0
        FCFRAC21 = 0.0
        FCFRAC31 = 0.0
        FCFRAC41 = 0.0
        FCFRAC51 = 0.0
        FCFRAC12 = 0.0
        FCFRAC22 = 0.0
        FCFRAC32 = 0.0
        FCFRAC42 = 0.0
        FCFRAC52 = 0.0
        TFRTCN1 = 0.0
        TFRTCN2 = 0.0
        TFRTCW1 = 0.0
        TFRTCW2 = 0.0
        LEAFDR1 = 0.0
        LEAFDR2 = 0.0
        LEAFDR3 = 0.0
        LEAFDR4 = 0.0
        LEAFDR5 = 0.0
        LEAFDR6 = 0.0
        LEAFDR7 = 0.0
        LEAFDR8 = 0.0
        LEAFDR9 = 0.0
        LEAFDR10 = 0.0
        LEAFDR11 = 0.0
        LEAFDR12 = 0.0
        BTOLAI = 0.0
        KLAI = 0.0
        LAITOP = 0.0
        MAXLAI = 0.0
        MAXLDR = 0.0
        FORRTF1 = 0.0
        FORRTF2 = 0.0
        FORRTF3 = 0.0
        SAPK = 0.0
        SWOLD = 0.0
        WDLIG1 = 0.0
        WDLIG2 = 0.0
        WDLIG3 = 0.0
        WDLIG4 = 0.0
        WDLIG5 = 0.0
        WDLIG6 = 0.0
        WOODDR1 = 0.0
        WOODDR2 = 0.0
        WOODDR3 = 0.0
        WOODDR4 = 0.0
        WOODDR5 = 0.0
        WOODDR6 = 0.0
        WRDSRFC = 0.0
        WMRTFRAC = 0.0
        SNFXMX2 = 0.0
        DEL13C = 0.0
        CO2IPR2 = 0.0
        CO2ITR2 = 0.0
        CO2ICE211 = 0.0
        CO2ICE212 = 0.0
        CO2ICE213 = 0.0
        CO2ICE221 = 0.0
        CO2ICE222 = 0.0
        CO2ICE223 = 0.0
        CO2IRS2 = 0.0
        BASFC2 = 0.0
        BASFCT = 0.0
        SITPOT = 0.0
        MAXNP = 0.0
        FKMRSPMX1 = 0.0
        FKMRSPMX2 = 0.0
        FKMRSPMX3 = 0.0
        FKMRSPMX4 = 0.0
        FKMRSPMX5 = 0.0
        FKMRSPMX6 = 0.0
        FMRSPLAI1 = 0.0
        FMRSPLAI2 = 0.0
        FMRSPLAI3 = 0.0
        FMRSPLAI4 = 0.0
        FMRSPLAI5 = 0.0
        FMRSPLAI6 = 0.0
        FGRESP1 = 0.0
        FGRESP2 = 0.0
        FGRESP3 = 0.0
        FGRESP4 = 0.0
        FGRESP5 = 0.0
        FGRESP6 = 0.0
        NO3PREF2 = 0.0
        TLAYPG = 0
        TMIX = 0.0
        TMPLFF = 0.0
        TMPLFS = 0.0
        FURGDYS = 0
        FLSGRES = 0.0
        TMXTURN = 0.0
        WSCOEFF21 = 0.0
        WSCOEFF22 = 0.0
        NPP2CS2 = 0.0
        SFAVAIL2 = 0.0
!! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        read (333,*,iostat=eof) ic, cname, idtype
         if (eof < 0) exit

        read (333,*,iostat=eof) salbd, slalb, slfg, omga, 
     &  inicd, avg_sla, ext_n, fcn_leafmin,ra, rb,
     &  ko25, kc25, vcmx25, ako, akc, avcmx, vpd_open, vpd_close, 
     &  toptd,
     &  gmax,gmin, chtmax, tbase, matyrs,DECID, PRDX2,
     &  PPDF1,
     &  PPDF2,
     &  PPDF3,
     &  PPDF4,
     &  CERFOR111,
     &  CERFOR112,
     &  CERFOR113,
     &  CERFOR121,
     &  CERFOR122,
     &  CERFOR123,
     &  CERFOR131,
     &  CERFOR132,
     &  CERFOR133,
     &  CERFOR141,
     &  CERFOR142,
     &  CERFOR143,
     &  CERFOR151,
     &  CERFOR152,
     &  CERFOR153,
     &  CERFOR211,
     &  CERFOR212,
     &  CERFOR213,
     &  CERFOR221,
     &  CERFOR222,
     &  CERFOR223,
     &  CERFOR231,
     &  CERFOR232,
     &  CERFOR233,
     &  CERFOR241,
     &  CERFOR242,
     &  CERFOR243,
     &  CERFOR251,
     &  CERFOR252,
     &  CERFOR253,
     &  CERFOR311,
     &  CERFOR312,
     &  CERFOR313,
     &  CERFOR321,
     &  CERFOR322,
     &  CERFOR323,
     &  CERFOR331,
     &  CERFOR332,
     &  CERFOR333,
     &  CERFOR341,
     &  CERFOR342,
     &  CERFOR343,
     &  CERFOR351,
     &  CERFOR352,
     &  CERFOR353,
     &  DECW1,
     &  DECW2,
     &  DECW3,
     &  FCFRAC11,
     &  FCFRAC21,
     &  FCFRAC31,
     &  FCFRAC41,
     &  FCFRAC51,
     &  FCFRAC12,
     &  FCFRAC22,
     &  FCFRAC32,
     &  FCFRAC42,
     &  FCFRAC52,
     &  TFRTCN1,
     &  TFRTCN2,
     &  TFRTCW1,
     &  TFRTCW2,
     &  LEAFDR1,
     &  LEAFDR2,
     &  LEAFDR3,
     &  LEAFDR4,
     &  LEAFDR5,
     &  LEAFDR6,
     &  LEAFDR7,
     &  LEAFDR8,
     &  LEAFDR9,
     &  LEAFDR10,
     &  LEAFDR11,
     &  LEAFDR12,
     &  BTOLAI,
     &  KLAI,
     &  LAITOP,
     &  MAXLAI,
     &  MAXLDR,
     &  FORRTF1,
     &  FORRTF2,
     &  FORRTF3,
     &  SAPK,
     &  SWOLD,
     &  WDLIG1,
     &  WDLIG2,
     &  WDLIG3,
     &  WDLIG4,
     &  WDLIG5,
     &  WDLIG6,
     &  WOODDR1,
     &  WOODDR2,
     &  WOODDR3,
     &  WOODDR4,
     &  WOODDR5,
     &  WOODDR6,
     &  WRDSRFC,
     &  WMRTFRAC,
     &  SNFXMX2,
     &  DEL13C,
     &  CO2IPR2,
     &  CO2ITR2,
     &  CO2ICE211,
     &  CO2ICE212,
     &  CO2ICE213,
     &  CO2ICE221,
     &  CO2ICE222,
     &  CO2ICE223,
     &  CO2IRS2,
     &  BASFC2,
     &  BASFCT,
     &  SITPOT,
     &  MAXNP,
     &  FKMRSPMX1,
     &  FKMRSPMX2,
     &  FKMRSPMX3,
     &  FKMRSPMX4,
     &  FKMRSPMX5,
     &  FKMRSPMX6,
     &  FMRSPLAI1,
     &  FMRSPLAI2,
     &  FMRSPLAI3,
     &  FMRSPLAI4,
     &  FMRSPLAI5,
     &  FMRSPLAI6,
     &  FGRESP1,
     &  FGRESP2,
     &  FGRESP3,
     &  FGRESP4,
     &  FGRESP5,
     &  FGRESP6,
     &  NO3PREF2,
     &  TLAYPG,
     &  TMIX,
     &  TMPLFF,
     &  TMPLFS,
     &  FURGDYS,
     &  FLSGRES,
     &  TMXTURN,
     &  WSCOEFF21,
     &  WSCOEFF22,
     &  NPP2CS2,
     &  SFAVAIL2



      if (eof < 0) exit


!! read local varibles to incoming arrays. added by Qichun Yang 2015. 8.12
!! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!      part1 photosynthesis variables

             cpnmfrst(ic) = cname
             idcfrst(ic) = idtype
             salbdf(ic) = salbd
             slalbf(ic) = slalb
             slfgf(ic) = slfg
             omgaf(ic) = omga 
             inicdf(ic) = inicd
             avg_slaf(ic) = avg_sla 
             ext_nf(ic) = ext_n
             fcn_leafminf(ic) = fcn_leafmin
             raf(ic) = ra 
             rbf(ic) = rb
             ko25f(ic) = ko25 
             kc25f(ic) = kc25 
             vcmx25f(ic) = vcmx25 
             akof(ic) = ako 
             akcf(ic) = akc 
             avcmxf(ic) = avcmx 
             vpd_openf(ic) = vpd_open 
             vpd_closef(ic) = vpd_close 
             toptdf(ic) = toptd
             gmaxf(ic) = gmax 
             gminf(ic) = gmin
             chtmaxf(ic) = chtmax
             tbasef(ic) = tbase 
             matyrsf(ic) = matyrs
             
  !!      part2 daycent variables
  
          DECIDf(ic) = DECID
        PRDX2f(ic) = PRDX2
        PPDFf(ic,1) = PPDF1
		PPDFf(ic,2) = PPDF2
		PPDFf(ic,3) = PPDF3
		PPDFf(ic,4) = PPDF4
        CERFORf(ic,1,1,1) = CERFOR111
		CERFORf(ic,1,1,2) = CERFOR112
		CERFORf(ic,1,1,3) = CERFOR113
		CERFORf(ic,1,2,1) = CERFOR121
		CERFORf(ic,1,2,2) = CERFOR122
		CERFORf(ic,1,2,3) = CERFOR123
		CERFORf(ic,1,3,1) = CERFOR131
		CERFORf(ic,1,3,2) = CERFOR132
		CERFORf(ic,1,3,3) = CERFOR133
		CERFORf(ic,1,4,1) = CERFOR141
		CERFORf(ic,1,4,2) = CERFOR142
		CERFORf(ic,1,4,3) = CERFOR143
		CERFORf(ic,1,5,1) = CERFOR151
		CERFORf(ic,1,5,2) = CERFOR152
		CERFORf(ic,1,5,3) = CERFOR153
		CERFORf(ic,2,1,1) = CERFOR211
		CERFORf(ic,2,1,2) = CERFOR212
		CERFORf(ic,2,1,3) = CERFOR213
		CERFORf(ic,2,2,1) = CERFOR221
		CERFORf(ic,2,2,2) = CERFOR222
		CERFORf(ic,2,2,3) = CERFOR223
		CERFORf(ic,2,3,1) = CERFOR231
		CERFORf(ic,2,3,2) = CERFOR232
		CERFORf(ic,2,3,3) = CERFOR233
		CERFORf(ic,2,4,1) = CERFOR241
		CERFORf(ic,2,4,2) = CERFOR242
		CERFORf(ic,2,4,3) = CERFOR243
		CERFORf(ic,2,5,1) = CERFOR251
		CERFORf(ic,2,5,2) = CERFOR252
		CERFORf(ic,2,5,3) = CERFOR253
		CERFORf(ic,3,1,1) = CERFOR311
		CERFORf(ic,3,1,2) = CERFOR312
		CERFORf(ic,3,1,3) = CERFOR313
		CERFORf(ic,3,2,1) = CERFOR321
		CERFORf(ic,3,2,2) = CERFOR322
		CERFORf(ic,3,2,3) = CERFOR323
		CERFORf(ic,3,3,1) = CERFOR331
		CERFORf(ic,3,3,2) = CERFOR332
		CERFORf(ic,3,3,3) = CERFOR333
		CERFORf(ic,3,4,1) = CERFOR341
		CERFORf(ic,3,4,2) = CERFOR342
		CERFORf(ic,3,4,3) = CERFOR343
		CERFORf(ic,3,5,1) = CERFOR351
		CERFORf(ic,3,5,2) = CERFOR352
		CERFORf(ic,3,5,3) = CERFOR353		
		DECWf(ic,1) = DECW1
		DECWf(ic,2) = DECW2
		DECWf(ic,3) = DECW3
        FCFRACf(ic,1,1) = FCFRAC11
        FCFRACf(ic,2,1) = FCFRAC21
		FCFRACf(ic,3,1) = FCFRAC31
		FCFRACf(ic,4,1) = FCFRAC41
		FCFRACf(ic,5,1) = FCFRAC51
		FCFRACf(ic,1,2) = FCFRAC12
		FCFRACf(ic,2,2) = FCFRAC22
		FCFRACf(ic,3,2) = FCFRAC32
		FCFRACf(ic,4,2) = FCFRAC42
		FCFRACf(ic,5,2) = FCFRAC52
        TFRTCNf(ic,1) = TFRTCN1
		TFRTCNf(ic,2) = TFRTCN2
		TFRTCWf(ic,1) = TFRTCW1
		TFRTCWf(ic,2) = TFRTCW2
		LEAFDRf(ic,1) = LEAFDR1
		LEAFDRf(ic,2) = LEAFDR2
		LEAFDRf(ic,3) = LEAFDR3
		LEAFDRf(ic,4) = LEAFDR4
		LEAFDRf(ic,5) = LEAFDR5
		LEAFDRf(ic,6) = LEAFDR6
		LEAFDRf(ic,7) = LEAFDR7
		LEAFDRf(ic,8) = LEAFDR8
		LEAFDRf(ic,9) = LEAFDR9
		LEAFDRf(ic,10) =LEAFDR10 
		LEAFDRf(ic,11) =LEAFDR11 
		LEAFDRf(ic,12) =LEAFDR12 
		BTOLAIf(ic) = BTOLAI
        KLAIf(ic) = KLAI
        LAITOPf(ic) = LAITOP
        MAXLAIf(ic) = MAXLAI
        MAXLDRf(ic) = MAXLDR
        FORRTFf(ic,1) = FORRTF1
		FORRTFf(ic,2) = FORRTF2
		FORRTFf(ic,3) = FORRTF3
        SAPKf(ic) = SAPK
        SWOLDf(ic) = SWOLD
        WDLIGf(ic,1) = WDLIG1
		WDLIGf(ic,2) = WDLIG2
		WDLIGf(ic,3) = WDLIG3
		WDLIGf(ic,4) = WDLIG4
		WDLIGf(ic,5) = WDLIG5
		WDLIGf(ic,6) = WDLIG6
        WOODDRf(ic,1) = WOODDR1
		WOODDRf(ic,2) = WOODDR2
		WOODDRf(ic,3) = WOODDR3
		WOODDRf(ic,4) = WOODDR4
		WOODDRf(ic,5) = WOODDR5
		WOODDRf(ic,6) = WOODDR6
        WRDSRFCf(ic) = WRDSRFC
        WMRTFRACf(ic) = WMRTFRAC 
        SNFXMX2f(ic) = SNFXMX2
        DEL13Cf(ic) = DEL13C
        CO2IPR2f(ic) = CO2IPR2
        CO2ITR2f(ic) = CO2ITR2
        CO2ICE2f(ic,1,1) = CO2ICE211
		CO2ICE2f(ic,1,2) = CO2ICE212
		CO2ICE2f(ic,1,3) = CO2ICE213
		CO2ICE2f(ic,2,1) = CO2ICE221
		CO2ICE2f(ic,2,2) = CO2ICE222
		CO2ICE2f(ic,2,3) = CO2ICE223
        CO2IRS2f(ic) = CO2IRS2
        BASFC2f(ic) = BASFC2
        BASFCTf(ic) = BASFCT
        SITPOTf(ic) = SITPOT
        MAXNPf(ic) = MAXNP
        FKMRSPMXf(ic,1) = FKMRSPMX1
		FKMRSPMXf(ic,2) = FKMRSPMX2
		FKMRSPMXf(ic,3) = FKMRSPMX3
		FKMRSPMXf(ic,4) = FKMRSPMX4
		FKMRSPMXf(ic,5) = FKMRSPMX5
		FKMRSPMXf(ic,6) = FKMRSPMX6
	    FMRSPLAIf(ic,1) = FMRSPLAI1
		FMRSPLAIf(ic,2) = FMRSPLAI2
		FMRSPLAIf(ic,3) = FMRSPLAI3
		FMRSPLAIf(ic,4) = FMRSPLAI4
		FMRSPLAIf(ic,5) = FMRSPLAI5
		FMRSPLAIf(ic,6) = FMRSPLAI6
        FGRESPf(ic,1) = FGRESP1
		FGRESPf(ic,2) = FGRESP2
		FGRESPf(ic,3) = FGRESP3
		FGRESPf(ic,4) = FGRESP4
		FGRESPf(ic,5) = FGRESP5
		FGRESPf(ic,6) = FGRESP6
        NO3PREF2f(ic) = NO3PREF2
        TLAYPGf(ic) = TLAYPG
        TMIXf(ic) = TMIX
        TMPLFFf(ic) = TMPLFF
        TMPLFSf(ic) = TMPLFS
        FURGDYSf(ic) = FURGDYS
        FLSGRESf(ic) = FLSGRES
        TMXTURNf(ic) = TMXTURN
        WSCOEFF2f(ic,1) = WSCOEFF21
		WSCOEFF2f(ic,2) = WSCOEFF22
        NPP2CS2f(ic) = NPP2CS2
        SFAVAIL2f(ic) = SFAVAIL2 
  
  
       

!! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  

      end do
        
      close (333)
      return
 5333 format (a80)
      end
