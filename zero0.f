      subroutine zero0

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine initializes the values for some of the arrays 

      use parm
      
      co_p = .1
      
      ifirstatmo = 1
      mo_atmo = 0
      
      ires_nut = 0
!!    apex command initialize
      idapa = 0
      iypa = 0
      flodaya = 0.
      seddaya = 0.
      orgndaya = 0.
      orgpdaya = 0.
      no3daya = 0.
      minpdaya = 0.
!!    apex command initialize

      lat_orgn = 0.
      lat_orgp = 0.
      cf = 0.0
      cfh = 0.0
      cfdec = 0.0
      pnd_d50 = 0.0
      isnow = 0
      imgt = 0
      iwtr = 0


!!    initialize mike van liew variables from .bsn file
      anion_excl_bsn = 0.
      ch_onco_bsn = 0.
      ch_opco_bsn = 0.
      hlife_ngw_bsn = 0.
      rcn_sub_bsn = 0.
      bc1_bsn = 0.
      bc2_bsn = 0.
      bc3_bsn = 0.
      bc4_bsn = 0.
      decr_min = 0.
!!    initialize mike van liew variables from .bsn file

      
      wtab = 0.8
!    Drainmod tile equations  01/2006 
      cumeira = 0.
      cumei = 0.
      cumrai = 0.
      cumrt = 0.
      ranrns_hru = 20.
!    Drainmod tile equations  01/2006
      afrt_surface = 0.
      alai_min = 0.
      lai_init = 0.
      aird = 0.
      alpha_bf = 0.
      alpha_bfe = 0.
      ammonian = 0.
      amp_r = 0.
      ap_ef = 0.
      atmo_day = 0.
      auto_eff = 0.
      auto_nyr = 0.
      auto_napp = 0.
      auto_nstrs = -1.
      auto_wstr = 0.
      bactkddb = 0.
      bactlp_plt = 0.
      bactlpdb = 0.
      bactp_plt = 0.
      bactpdb = 0.
      bc1 = 0.
      bc2 = 0.
      bc3 = 0.
      bc4 = 0.
      bio_e = 0.
      bio_hv = 0.
      bio_init = 0.
      bio_leaf = 0.
      bio_min = 0.
      bio_ms = 0.
      bio_n1 = 0.
      bio_n2 = 0.
      bio_p1 = 0.
      bio_p2 = 0.
      bio_targ = 0.
      blai = 0.
      bio_eat = 0.
      bio_trmp = 0.
      brt = 0.
      bss = 0.
      canmx = 0.
      canstor = 0.
      cbodcnst = 0.
      cbodmon = 0.
      cbodyr = 0.
      cfrt_id = 0
      cfrt_kg = 0.
      ch_d = 0.
      ch_di = 0.
      ch_k = 0.
      ch_li = 0.
      ch_onco = 0.
      ch_opco = 0.
      ch_orgn = 0.
      ch_orgp = 0.
      ch_n = 0.
      ch_s = 0.
      ch_si = 0.
      ch_w = 0.
      ch_wi = 0.
      ch_wdr = 0.

!!    Initialization by balaji
      ch_bnk_bd = 0.
      ch_bed_bd = 0.
      ch_bnk_kd = 0.
      ch_bed_kd = 0.
      tc_bnk = 0.
      tc_bed = 0.
      ch_bnk_d50 = 0.
      ch_bed_d50 = 0.
      ch_bed_san = 0.
      ch_bed_sil = 0.
      ch_bed_cla = 0.
      ch_bed_gra = 0.
      ch_bnk_san = 0.
      ch_bnk_sil = 0.
      ch_bnk_cla = 0.
      ch_bnk_gra = 0.
      ch_eqn = 0 
      ch_erodmo = 0.      !CB 12/2/09
      ch_cov1 = 0.
      ch_cov2 = 0.
      chlacnst = 0.
      chlamon = 0.
      chlayr = 0.
      chlora = 0.
      chpst_conc = 0.
      chpst_koc = 0.
      chpst_mix = 0.
      chpst_rea = 0.
      chpst_rsp = 0.
      chpst_stl = 0.
      chpst_vol = 0.
      chside = 0.
      chtmx = 0.
      cncoef_sub = 0.
      cn1 = 0.
      cn2 = 0.
      cn3 = 0.
      cnop = 0.
      cnyld = 0.
      co2 = 0.
!    Drainmod tile equations  01/2006 
      conk = 0.
!    Drainmod tile equations  01/2006
      conv_wt = 0.
      cpst_id = 0
      cpnm = ""
      cpyld = 0.
      curbden = 0.
      curyr = 0
      curyr_mat = 0
      igrotree = 0
      cvm = 0.
      daylmn = 0.
      daru_km = 0.
!    Drainmod tile equations  01/2006 
      dc = 0.
      drain_co_bsn = 0.
!    Drainmod tile equations  01/2006
      ddrain = 0.
      ddrain_bsn = 0.
      decay_f = 0.
      decay_s = 0.
      deepirr = 0.
      deepst = 0.
      delay = 0.
      depimp_bsn = 0.
      dep_imp = 0.
      deptil = 0.
      det_san = 0.
      det_sil = 0.
      det_cla = 0.
      det_sag = 0.
      det_lag = 0.
      dewpt = 0.
      dirtmx = 0.
      dis_stream = 0.
      disolvp = 0.
      disoxcnst = 0.
      disoxmon = 0.
      disoxyr = 0.
      divmax = 0.
      dlai = 0.
      dormhr = 0.
      dorm_hr = -1.
      dorm_flag = 0
      driftco = 0.
!    Drainmod tile equations  01/2006 
      dtwt = 600.
      dz = 0.
!    Drainmod tile equations  01/2006
      drydep_no3 = 0.
      drydep_nh4 = 0.
      eo_30d = 0.
      effmix = 0.
      elevb = 0.
      elevb_fr = 0.
      epco = 0.
      esco = 0.
      evpnd = 0.
      evpot = 0.
      evrsv = 0.
      evwet = 0.
      fcst_reg = 0
      fcstaao = 0.
      fcstcnt = 0
      fertnm = ""
      ffc = 0.
      ffcst = 0
      filterw = 0.
      fimp = 0.
      fld_fr = 0.
      flomon = 0.
      flowmin = 0.
      fminn = 0.
      fminp = 0.
      fnh3n = 0.
      forgn = 0.
      forgp = 0.
      fpcp_stat = 0.
      fpr_w = 0.
      frac_harvk = 0.
      frt_kg = 0.
      frt_surface = 0.
      ftmpmn = 0.
      ftmpmx = 0.
      ftmpstdmn = 0.
      ftmpstdmx = 0.
      gdrain = 0.
      gdrain_bsn = 0.
!    Drainmod tile equations  01/2006 
      gee = 0.
!    Drainmod tile equations  01/2006
      gsi = 0.
      gw_delaye = 0.
      gw_q = 0.
      gw_revap = 0.
      gw_spyld = 0.
      gwht = 0.
      gwq_ru = 0.
!     latq_ru = 0.
!      surfq_ru = 0.
!      infl_ru = 0.
      gwqmn = 0.
!    Drainmod tile equations  01/2006 
      hdrain = 0.
      hdrain_bsn = 0.
!    Drainmod tile equations  01/2006
      hi_ovr = 0.
      hi_targ = 0.
      hlife_f = 0.
      hlife_s = 0.
      hqd = 0.
      hqdsave = 0.
      hru_dafr = 0.
      hru_fr = 0.
      hru_ha = 0.
      hru_km = 0.
      hru_rufr = 0.
      hrugis = 0
      hrupest = 0
      hrupsta = 0.
      hrupstm = 0.
      hrupsty = 0.
      hru_rufr = 0.
      huminc = 0.
      hvsti = 0.
      hyd_dakm = 0.
      iafrttyp = 0
      iatmodep = 0
      icalen = 0
      icfrt = 0
      icodes = 0
      icpst = 0
      icr = 0
      icrmx = 0
      iday_fert = 0
      idc = 0
      idop = 0
      idorm = 0
      idplt = 0
      idplrot = 1
      idtill = 0
      ihv_gbm = 0
      wstrs_id = 0
      ifirstpcp = 1
      ifirsthr = 1
      ifirsta = 1 
      ifirstr = 1 
      ifirstt = 1 
      ifld = 0
      iflod1 = 0
      iflod2 = 0
      ifrt_freq = 0
      ifrttyp = 0
      irelh = 1
      manure_id = 0
      no_lup = 1
      igro = 0
      igrz = 0
      ihouts = 0
      ils2 = 0
      ils2flag = 0
      ils_nofig = 0
      inum1s = 0
      inum2s = 0
      inum3s = 0
      inum4s = 0
      inum5s = 0
      inum6s = 0
      inum7s = 0
      inum8s = 1
      iop = 0
      ioper = 1
      iopera = 1
      iopday = 0
      iopyr = 0
      ipdhru = 0
      ipdvab = 0
      ipdvar = 0
      ipdvas = 0
      ipest = 0
 !!     ipot = 0
      ireg = 1
      irgage = 0
      irip = 0
      iroutunit = 0
      irn = 0
      irramt = 0.
      irreff = 1.
      irrefm = 1.
      irrsalt = 0.
      irrsc = 0
      irrno = 0
      irr_sca = 0
      irr_noa = 0
      irrsq = 0
      irrno = 0
      irrsc = 0
      irr_sc = 0
      irr_no = 0
      irr_sq = 0
      irr_asq = 0
      irr_sca = 0
      irr_noa = 0
      itb = 0 
      itemp = 0
      itgage = 0
      iurban = 0
      ivar_orig = 0
      kirr = " "
      laiday = 0.
!    Drainmod tile equations  02/2009 
      latksatf = 0.
      latksatf_bsn = 0.
!    Drainmod tile equations  02/2009
      lat_sed = 0.
      lat_ttime = 0.
      latcos = 0.
      latsin = 0.
      leaf1 = 0.
      leaf2 = 0.
      mcr = 1
      mcrhru = 0
      mgtop = 0
      mgt1iop = 0
      mgt2iop = 0
      mgt3iop = 0
      mgt4op = 0.0
      mgt5op = 0.0
      mgt6op = 0.0
      mgt7op = 0.0
      mgt8op = 0.0
      mgt9op = 0.0
      mgt10iop = 0
      mo_transb = 0
      mo_transe = 0
      ncrops = 0 
      ncut = 1
      ndcfrt = 0
      fert_days = 0
      grz_days = 0
      nair = 1
      irr_mx = 0.
      latno3 = 0.
      nicr = 0
      ndmo = 0
      ndtarg = 0
      newrti = 0.
      nitraten = 0.
      nitriten = 0.
      nmgt = 0
      nope = 0
      nopmx = 0
      npcp = 1
      npno = 0
      nro = 1
      nrot = 0
      ntil = 1
      orig_alai = 0.
      orig_bioms = 0.
      orig_deepst = 0.
      orig_igro = 0
      organicn = 0.
      organicp = 0.
      orgn_con = 0.
      orgp_con = 0.
      ov_n = 0.
!    Drainmod tile equations  01/2006 
      pc = 0.
      pc_bsn = 0.
!    Drainmod tile equations  01/2006
      phubase = 0.
      pltnfr = 0.
      pltpfr = 0.
      pot_seep = 0.
!----------------------- !Moriasi 4/8/2014  
      prf = 0.
      prf_bsn = 0. 
      spcon_bsn = 0.
      spexp_bsn = 0.
      r2adj = 0.
      r2adj_bsn = 1.
!----------------------- !Moriasi 4/8/2014       
!! drainmod tile equations   06/2006
      ranrns = 0.
!! drainmod tile equations   06/2006
      qird = 0.
      rammo_sub = 0.
      rch_cbod = 0.
      rch_dox = 0.
      rchrg_n = 0.
      rcn_sub = 0.
      rcn_sub_bsn = 0.
      reccnstps = '             '
      recmonps = '             '
      rammo_mo = 0.
      rcn_mo = 0.
      drydep_nh4_mo = 0.
      drydep_no3_mo = 0.
      rammo_d = 0.
      rcn_d = 0.
      drydep_nh4_d = 0.
      drydep_no3_d = 0.
!! routing 5/3/2010 gsm per jga
      idum = 0
      mhyd1 = 0
      irtun = 0
      
!    Drainmod tile equations  01/2006 
      sdrain = 0.
      sdrain_bsn = 0.
      sstmaxd = 0.
      sstmaxd_bsn = 0.
 !New water table depth parameters D. Moriasi 4/8/2014
      sol_swpwt = 0.
      sol_stpwt = 0.
      vwt = 0.
      wat_tbl = 0.         
!    Drainmod tile equations  01/2006
      rsr1 = 0.
      rsr2 = 0.
      rsr1 = 0.
      rsr2 = 0.
      sed_con = 0.
      sepcrk = 0.
      sq_rto = 0.
      sol_clay = 0. 
!    Drainmod tile equations  01/2006 
      stmaxd = 0.
      stmaxd_bsn = 0.
!    Drainmod tile equations  01/2006      
      sol_ec = 0.
      sol_sand = 0.
      sol_silt = 0.
      sol_clay = 0.
!!   added for Srini in output.mgt nitrogen and phosphorus nutrients per JGA by gsm 9/8/2011
      sol_sumn03 = 0.
      sol_sumsolp = 0.
      strsw = 1.
      strsw_sum = 0.
      strstmp_sum = 0.
      strsn_sum = 0.
      strsp_sum = 0.
      strsa_sum = 0.
      stsol_rd = 0.
      soln_con = 0.
      solp_con = 0.
      subed = 0
      sub_elev = 0.
      subfr_nowtr = 0.
      sub_lat = 0.
      sub_latq = 0.
      sub_tileq = 0.
      sub_latno3 = 0.
      sub_smtmp = 0.
      sub_tileno3 = 0.
      sub_gwq_d = 0.
      alpha_bf_d = 0.
      gw_qdeep = 0.
      subgis = 0
      tb_adj = 0.
      tdrain = 0.
      tdrain_bsn = 0.
      tile_no3 = 0.
      tileq = 0.
      tile_ttime = 0.
      uh = 0.
      vfsratio = 0.
      vfscon = 0.
      vfsi = 0.
      vfsch = 0.
      wt_shall = 0.
      wshd_aamon = 0.
!      wshddayo = 0.
      yr_skip = 0


          !! initialize variables, added by Qichun Yang
      !! need to change the initial values of carbon pools based on the initial values in the parameter table

      idcfrst = 0.
      cpnmfrst = 'aaaa'
!!         part1 photosynthesis
             salbdf = 0.
             slalbf = 0.
             slfgf = 0.
             omgaf = 0.
             inicdf = 0.
             avg_slaf = 0. 
             ext_nf = 0.
             fcn_leafminf = 0.
             raf = 0. 
             rbf = 0.
             ko25f = 0.
             kc25f = 0. 
             vcmx25f = 0. 
             akof = 0. 
             akcf = 0. 
             avcmxf = 0. 
             vpd_openf = 0.
             vpd_closef = 0.
             toptdf = 0.
             gmaxf = 0. 
             gminf = 0.
             chtmaxf = 0.
             tbasef = 0.
!!         part2 daycent variables

         wood1tosom1 = 0.
         wood1tosom2 = 0.
         wood2tosom1 = 0.
         wood2tosom2 = 0.
         wood3tosom1 = 0.
         wood3tosom2 = 0.
         
         tswatlyr = 0
         rdis = 0.
         
        DECIDf = 0
        PRDX2f = 0.
        PPDFf = 0.
        CERFORf = 0.
        DECWf = 0.
        FCFRACf = 0.
        TFRTCNf = 0.
        TFRTCWf = 0.
        LEAFDRf = 0.
        BTOLAIf = 0.
        KLAIf = 0.
        LAITOPf = 0.
        MAXLAIf = 0.
        MAXLDRf = 0.
        FORRTFf = 0.
        SAPKf = 0.
        SWOLDf = 0.
        WDLIGf = 0.
        WOODDRf = 0.
        WRDSRFCf = 0.
        WMRTFRACf = 0.
        SNFXMX2f = 0.
        DEL13Cf = 0.
        CO2IPR2f = 0.
        CO2ITR2f = 0.
        CO2ICE2f = 0.
        CO2ICE2f = 0.
        CO2ICE2f = 0.
        CO2ICE2f = 0.
        CO2ICE2f = 0.
        CO2ICE2f = 0.
        CO2IRS2f = 0.
        BASFC2f = 0.
        BASFCTf = 0.
        SITPOTf = 0.
        MAXNPf = 0.
        FKMRSPMXf = 0.
        FMRSPLAIf = 0.
        FGRESPf = 0.
        NO3PREF2f = 0.
        TLAYPGf = 0
        TMIXf = 0.
        TMPLFFf = 0.
        TMPLFSf = 0.
        FURGDYSf = 0
        FLSGRESf = 0.
        TMXTURNf = 0.
        WSCOEFF2f = 0.
        NPP2CS2f = 0.
        SFAVAIL2f = 0.
        BASFC2f = 0.
        hrsinc = .FALSE.
        decidgrow  = .FALSE.
        ffroota = 0.
        
       
!!      part3, biomass pools
        leafc = 0.
        leafn = 0.
        leafp = 0.
        brchc = 0.
        brchn = 0.
        brchp = 0.
        largwc = 0.
        largwn = 0.
        largwp = 0.
        frootjc = 0.
        frootjn = 0.
        frootjp = 0.
        frootmc = 0.
        frootmn = 0.
        frootmp = 0.
        frootp = 0.
        csrootc = 0.
        csrootn = 0.
        csrootp = 0.
        seedc = 0.
        seedn = 0.
        seedp = 0.
        
        eleaf = 0.
        ebrch = 0. 
        elargw = 0.
        efrootj = 0.
        efrootm = 0.
        ecsroot = 0.
        
        csrsnk = 0.
        esrsnk = 0.
        carbostg = 0.
        forstg = 0.
        availm = 0.
        co2cce = 0.
        tree_a2drat = 0.
        woode = 0.
        woodc = 0.
!!      part4 carbon fluxes        
        flai = 0.
        f_gpp = 0.
        f_pgpp = 0.
        f_npp = 0.
        f_mr = 0.
        f_gr = 0.
        ab_litter = 0.
        bl_litter = 0.
        f_moist = 0.
        tavgdate = 0.
        pwther = 0.
!! assign values to the constant vairables from daycent  

!! qichun varibles from daycent
!!  auto matic fertlizer
      aufert = 0. 
      co2crs = 0.
      nelem = 2

!!... CRPSYS, FORSYS, SAVSYS is an enumeration for the system
      CRPSYS = 1
      FORSYS = 2
      SAVSYS = 3

!!... MAXIEL is the total # of elements available
!!...   N = Nitrogen, P = Phosphorus, S = Sulphur
      MAXIEL = 3
      N = 1
      P = 2
      S = 3

!!... MAXLYR is the maximum number of layers
      MAXLYR = 10

!!... SWMAXLYR is the maximum number of layers in the soil water model
      SWMAXLYR = 21

!!... MONTHS is the number of months
      MONTHS = 12

!!... SRFC, SOIL are enumerations for surface, soil
      SRFC = 1
      SOIL = 2

!!... CPARTS is the number of parts in the grassland/crop system:
!!...   ABOVE is an enumeration for aboveground
!!...   BELOW is an enumeration for belowground
!!...   BELOWJ is an enumeration for juvenile belowground
!!...   BELOWM is an enumeration for mature belowground
      CPARTS = 3
      ABOVE = 1
      BELOW = 2
      BELOWJ = 2
      BELOWM = 3

!!... FPARTS is the number of parts in the forest system:
!!...   LEAF   =  leaf forest part
!!...   FROOT =  fine root forest part
!!...   FROOTJ =  juvenile fine root forest part
!!...   FROOTM =  mature fine root forest part
!!...   FBRCH  =  fine branch forest part
!!...   LWOOD  =  large wood forest part
!!...   CROOT  =  coarse root forest part
      FPARTS = 6
      LEAF = 1
      FROOT = 2
      FROOTJ = 2
      FBRCH = 3
      LWOOD = 4
      CROOT = 5
      FROOTM = 6

!!... NEWFOR, OLDFOR are the new and old forests array pointers
      NEWFOR = 1
      OLDFOR = 2

!!... IMIN, IMAX, IVAL are enumerations for indexing arrays
      IMIN = 1
      IMAX = 2
      IVAL = 3

!!... INTCPT, SLOPE are the intercept and slope array pointers
      INTCPT = 1
      SLOPE = 2

!!... UNLABL, LABELD are the unlabeled, labeled array pointers
      UNLABL = 1
      LABELD = 2

!!... ISOS is the total number of isotopes (unlabeld, labeled)
      ISOS = 2

!!... ESTOR, ESOIL, ENFIX, EFERT are enumerations used only in
!!... restrp.f, growth.f, trees.f
      ESTOR = 1
      ESOIL = 2
      ENFIX = 3
      EFERT = 4

!!... Constant values
!!... Change value of PEEDEE constant used in del13c computations,
!!... cak - 12/19/01
c      parameter (PEEDEE=0.0112372)
      PEEDEE=0.01119490

!!... Fraction of the total carbon that is C14 carbon
c      parameter (FRAC_C14=0.00103)
      FRAC_C14=0.00000103     
      DPI = 3.141592653589793
      growday = 0
      pligst(1) = 3.0
      pligst(2) = 3.0

      return
      end
