MODULE usrdef_nam
   !!======================================================================
   !!                       ***  MODULE  usrdef_nam  ***
   !!
   !!                      ===  BASIN configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2017-10  (J. Chanut)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce, ONLY: nimpp , njmpp            ! i- & j-indices of the local domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam, merc_proj, usr_def_nam_cfg   ! called by nemogcm.F90

   !                              !!* namusr_def namelist *!!
   REAL(wp), PUBLIC ::   rn_e1_deg      =     1     ! Resolution in degrees of longitude (Mercator grid)
   REAL(wp), PUBLIC ::   rn_phi_min     =   -70.    ! Latitude of the south frontier (T point) [degrees] (approximative)
   REAL(wp), PUBLIC ::   rn_phi_max     =    70.    ! Latitude of the north frontier (T point) [degrees] (approximative)
   REAL(wp), PUBLIC ::   rn_lam_min     =     0.    ! Longitude of the west frontier (T point) [degrees] (approximative)
   REAL(wp), PUBLIC ::   rn_lam_max     =    50.    ! Longitude of the east frontier (T point) [degrees] (approximative)
   INTEGER , PUBLIC ::   nn_jeq_s       =    100    ! Number of grid points from southern boundary to equator
   INTEGER , PUBLIC ::   nn_iglo        =    42     ! Number of grid points along i direction
   INTEGER , PUBLIC ::   nn_jglo        =    62     ! Number of grid points along j direction
   INTEGER , PUBLIC ::   nn_k           =    42     ! Number of grid points along k direction
   INTEGER , PUBLIC ::   nn_forcingtype =     0     ! Surface forcing type
   REAL(wp), PUBLIC ::   rn_ztau0       =    0.1    ! Magnitude of the wind-forcing [Pa]
   LOGICAL , PUBLIC ::   ln_diu_cyc     =  .TRUE.   ! Use diurnal cycle for qsr or not
   LOGICAL , PUBLIC ::   ln_ann_cyc     =  .TRUE.   ! Use an annual cycle or not
   LOGICAL , PUBLIC ::   ln_emp_field   =  .FALSE.  ! EmP is read from file and replaces saltflx from S-restoring
   LOGICAL , PUBLIC ::   ln_qns_field   =  .FALSE.  ! Qtot is read from file and replaces heatflx from T-restoring
   REAL(wp), PUBLIC ::   rn_emp_prop    =     1.    ! Proportionality factor to apply on the EMP
   REAL(wp), PUBLIC ::   rn_trp         =   -40.    ! Retroaction term on T*, must be negative  [W/m2/K]
   REAL(wp), PUBLIC ::   rn_srp         =     0.0   ! Restoring term on S*, must be negative    [kg/m2/s]
   REAL(wp), PUBLIC ::   rn_sstar_s     =     35.   ! Salinity restoring value at the southern boundary
   REAL(wp), PUBLIC ::   rn_sstar_n     =     35.1  ! Salinity restoring value at the northern boundary 
   REAL(wp), PUBLIC ::   rn_sstar_eq   =     37.25 ! Salinity restoring maximum without exponential
   REAL(wp), PUBLIC ::   rn_tstar_s     =    -0.5   ! Retroaction value on T at the southern boundary
   REAL(wp), PUBLIC ::   rn_tstar_n     =     5.0   ! Retroaction value on T at the northern boundary
   REAL(wp), PUBLIC ::   rn_tstar_eq    =     27.   ! Retroaction value on T at the equator
   LOGICAL , PUBLIC ::   ln_qsr         =  .TRUE.   ! Solar heat or not
   INTEGER , PUBLIC ::   nn_botcase     =     0     ! bottom definition (0:flat, 1:bowl in cosh, 2: bowl in 1-x**4)
   REAL(wp), PUBLIC ::   rn_H           =  4000.    ! Maximum depth of the basin or asymptotical depth of the basin (depending on the basin shape)
   REAL(wp), PUBLIC ::   rn_hborder     =  2000.    ! Depth of the borders of the basin
   REAL(wp), PUBLIC ::   rn_distLam     =     3.    ! Typical length scale of the slope in longitude
   REAL(wp), PUBLIC ::   rn_dzmin       =     10.   ! minimum value of e3 at the surface   [m]
   REAL(wp), PUBLIC ::   rn_kth         =     35.   ! position of the inflexion point
   REAL(wp), PUBLIC ::   rn_hco         =   1000.   ! layer thickness with z-coordinate [m]
   REAL(wp), PUBLIC ::   rn_acr         =    10.5   ! slope of the tanh
   INTEGER , PUBLIC ::   nn_initcase    =     0     ! initial condition case (0=rest+uniform, 1=rest+stratification)
   LOGICAL , PUBLIC ::   ln_Iperio      = .TRUE.    ! periodicity in i
   REAL(wp), PUBLIC ::   rn_cha_min     =   -60.    ! chanel width [degrees]
   REAL(wp), PUBLIC ::   rn_cha_max     =   -50.    ! midpoint latitude of the chanel [degrees] (approximate)
   REAL(wp), PUBLIC ::   rn_slp_cha     =    1.5    ! Slope of the boundaries in the channel
   LOGICAL , PUBLIC ::   ln_zco_nam     = .FALSE.   ! z               vertical coordinate
   LOGICAL , PUBLIC ::   ln_zps_nam     = .FALSE.   ! z-partial steps vertical coordinate
   LOGICAL , PUBLIC ::   ln_sco_nam     = .TRUE.    ! s               vertical coordinate
   INTEGER , PUBLIC ::   nn_ztype       =     0     ! type of vertical grid spacing (0: uniform) (z-coordinate or s pure coordinate)
   LOGICAL , PUBLIC ::   ln_mid_ridge   = .FALSE.   ! Including a Mid-Atlantic Ridge
   REAL(wp), PUBLIC ::   rn_mr_depth    =  2000.    ! Depth of the Mid-Atlantic ridge
   REAL(wp), PUBLIC ::   rn_mr_width    =     8.    ! Width of the Mid-Atlantic ridge
   REAL(wp), PUBLIC ::   rn_mr_lat_s    =   -50.    ! southern edge of the Mid-Atlantic ridge [degrees]
   REAL(wp), PUBLIC ::   rn_mr_lat_n    =    55.    ! northern edge of the Mid-Atlantic ridge [degrees]
   INTEGER , PUBLIC ::   nn_mr_edge     =      1    ! shape of the southern/northern edge
   LOGICAL , PUBLIC ::   ln_drake_sill  = .FALSE.   ! Including circular Drake-Sill
   REAL(wp), PUBLIC ::   rn_ds_depth   =  3000.     ! Depth of the circular Drake-Sill
   REAL(wp), PUBLIC ::   rn_ds_width   =     2.     ! Width of the circular Drake-Sill 
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 10074 2018-08-28 16:15:49Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_nam( cd_cfg, kk_cfg, kpi, kpj, kpk, ldIperio, ldJperio, ldNFold, cdNFtype )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read user defined namelist and define the domain size
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here BASIN configuration
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(out) ::   cd_cfg               ! configuration name
      INTEGER         , INTENT(out) ::   kk_cfg               ! configuration resolution
      INTEGER         , INTENT(out) ::   kpi, kpj, kpk        ! global domain sizes
      LOGICAL         , INTENT(out) ::   ldIperio, ldJperio   ! i- and j- periodicity
      LOGICAL         , INTENT(out) ::   ldNFold              ! North pole folding
      CHARACTER(len=1), INTENT(out) ::   cdNFtype             ! Folding type: T or F
      !
      INTEGER ::   ios, ii               ! Local integer
      REAL(wp)::   zh, ziglo             ! Local scalars
      ! REAL(wp)::   zarg_min, zarg_max, zjeq_min, zjeq_max, ijeq_max, rn_iglo
      INTEGER ::   nn_jeq_n
      !!
      NAMELIST/namusr_def/ rn_e1_deg, rn_phi_min, rn_phi_max, rn_lam_min         &
         &                 , rn_lam_max, nn_k, rn_emp_prop, rn_ztau0             &
         &                 , nn_botcase, nn_initcase, nn_forcingtype             &
         &                 , ln_Iperio, rn_cha_min, rn_cha_max, rn_slp_cha       &
         &                 , ln_zco_nam, ln_zps_nam, ln_sco_nam                  &
         &                 , nn_ztype, rn_H, rn_hborder, rn_distLam              &
         &                 , ln_mid_ridge, ln_drake_sill, ln_ann_cyc             &
         &                 , ln_qns_field, ln_emp_field                          &
         &                 , rn_trp, rn_srp, ln_qsr, ln_diu_cyc                  &
         &                 , rn_sstar_s, rn_sstar_n, rn_sstar_eq                 &
         &                 , rn_tstar_s, rn_tstar_n, rn_tstar_eq                 &
         &                 , rn_dzmin, rn_kth, rn_hco, rn_acr,  nn_mr_edge       &
         &                 , rn_mr_depth, rn_mr_width, rn_mr_lat_s               &
         &                 , rn_mr_lat_n, rn_ds_depth, rn_ds_width
      !!----------------------------------------------------------------------
      !
      ii = 1

      !REWIND( numnam_cfg )          ! Namelist namusr_def (exist in namelist_cfg only)
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      !
      WRITE( numond, namusr_def )
      !
      cd_cfg = 'DINO'             ! name & resolution (not used)
      kk_cfg = rn_e1_deg
      !
      nn_jeq_n = merc_proj(rn_phi_max, rn_e1_deg)
      nn_jeq_s = merc_proj(rn_phi_min, rn_e1_deg)
      IF(lwp) WRITE(numout,*) '          Index of the equator (north from merc_proj) on the MERCATOR grid:', nn_jeq_n
      IF(lwp) WRITE(numout,*) '          Index of the equator (south from merc_proj) on the MERCATOR grid:', nn_jeq_s
      nn_jglo = (nn_jeq_s - nn_jeq_n) + 1
      ! Number of gridpoints in i-direction
      ziglo = (rn_lam_max - rn_lam_min) / rn_e1_deg
      nn_iglo = FLOOR(ziglo)
      IF( ABS( REAL( nn_iglo, wp ) - ziglo ) > 0.5 )   nn_iglo = nn_iglo + 1
      ! To conserve volume across resolutions in i-direction
      nn_iglo = nn_iglo + 2
      ! Global domain size
      kpi = nn_iglo
      kpj = nn_jglo
      kpk = nn_k
      !
      ldIperio = ln_Iperio   ;   ldJperio = .FALSE.   ! DINO configuration : with periodic channel 
      ldNFold  = .FALSE.     ;   cdNFtype = '-'
      
      !
      !kperio = nn_perio   ! 0: closed basin, 8: south symmetrical
      !IF( kperio == 8 )   rn_phi0 = -rn_e1_deg   ! manually set the longitude of the global first (southern) T point
      !                             ! control print
      WRITE(numout,*) '   '                                                                                                  ;   ii = ii + 1
      WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'                           ;   ii = ii + 1
      WRITE(numout,*) '~~~~~~~~~~~ '                                                                                         ;   ii = ii + 1
      WRITE(numout,*) '   Namelist namusr_def : DINO test case'                                                             ;   ii = ii + 1
      WRITE(numout,*) '   Resolution in degrees of longitude (Mercator grid)      rn_e1_deg      = ', rn_e1_deg, 'degrees'   ;   ii = ii + 1
      WRITE(numout,*) '   Latitude of the south frontier (T point) [degrees]      rn_phi0_min    = ', rn_phi_min, 'degrees'  ;   ii = ii + 1
      WRITE(numout,*) '   Latitude of the south frontier (T point) [degrees]      rn_phi_max     = ', rn_phi_max, 'degrees'  ;   ii = ii + 1
      WRITE(numout,*) '   Latitude of the south frontier (T point) [degrees]      rn_lam_min     = ', rn_lam_min, 'degrees'  ;   ii = ii + 1
      WRITE(numout,*) '   Latitude of the south frontier (T point) [degrees]      rn_lam_max     = ', rn_lam_max, 'degrees'  ;   ii = ii + 1
      WRITE(numout,*) '   Number of grid points along i direction                 nn_piglo       = ', kpi                    ;   ii = ii + 1
      WRITE(numout,*) '   Number of grid points along j direction                 nn_pjglo       = ', kpj                    ;   ii = ii + 1
      WRITE(numout,*) '   Number of grid points along k direction                 nn_k           = ', kpk                    ;   ii = ii + 1
      WRITE(numout,*) '   Use surface forcing like Gyre (1) or remade one (0)     nn_forcingtype = ', nn_forcingtype         ;   ii = ii + 1
      WRITE(numout,*) '   Magnitude of wind forcing                               rn_ztau0       = ', rn_ztau0               ;   ii = ii + 1
      WRITE(numout,*) '   Use annual cycle or not                                 ln_ann_cyc     = ', ln_ann_cyc             ;   ii = ii + 1
      WRITE(numout,*) '   Proportionality factor applied on EMP                   rn_emp_prop    = ', rn_emp_prop            ;   ii = ii + 1
      WRITE(numout,*) '   Bottom definition                                       nn_botcase     = ', nn_botcase             ;   ii = ii + 1
      WRITE(numout,*) '      (0:flat, 1:bowl in cosh, 2: bowl in 1-x**4)'                                                    ;   ii = ii + 1
      WRITE(numout,*) '   Maximum / asymptotical depth of the basin               rn_H           = ', rn_H, 'm'              ;   ii = ii + 1
      WRITE(numout,*) '   Depth of the bathymetry at the coast                    rn_hborder     = ', rn_hborder, 'm'        ;   ii = ii + 1
      WRITE(numout,*) '   Typical length scale of the slope in longitude          rn_distLam     = ', rn_distLam, 'degrees'
      WRITE(numout,*) '   Initial condition case                                  nn_initcase    = ', nn_initcase            ;   ii = ii + 1
      WRITE(numout,*) '      (0:rest and constant T and S, '                                                                 ;   ii = ii + 1
      WRITE(numout,*) '       1: rest and stratification)'                                                                   ;   ii = ii + 1
      WRITE(numout,*) '      (0:closed, 8:south symmetrical)'                                                                ;   ii = ii + 1
      WRITE(numout,*) '   Include mid-atlantic ridge                              ln_mid_ridge   = ', ln_mid_ridge           ;   ii = ii + 1
      WRITE(numout,*) '   Include circular drake sill                             ln_drake_sill  = ', ln_drake_sill          ;   ii = ii + 1
   END SUBROUTINE usr_def_nam


   SUBROUTINE usr_def_nam_cfg( )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read user defined namelist when reading from domain configuration file (ln_read_cfg=T)
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here DINO configuration
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      !
      INTEGER ::   ios, ii               ! Local integer
      !!
      NAMELIST/namusr_def/ rn_e1_deg, rn_phi_min, rn_phi_max, rn_lam_min         &
         &                 , rn_lam_max, nn_k, rn_emp_prop, rn_ztau0             &
         &                 , nn_botcase, nn_initcase, nn_forcingtype             &
         &                 , ln_Iperio, rn_cha_min, rn_cha_max, rn_slp_cha       &
         &                 , ln_zco_nam, ln_zps_nam, ln_sco_nam                  &
         &                 , nn_ztype, rn_H, rn_hborder, rn_distLam              &
         &                 , ln_mid_ridge, ln_drake_sill, ln_ann_cyc             &
         &                 , ln_qns_field, ln_emp_field                          &
         &                 , rn_trp, rn_srp, ln_qsr, ln_diu_cyc                  &
         &                 , rn_sstar_s, rn_sstar_n, rn_sstar_eq                 &
         &                 , rn_tstar_s, rn_tstar_n, rn_tstar_eq                 &
         &                 , rn_dzmin, rn_kth, rn_hco, rn_acr,  nn_mr_edge       &
         &                 , rn_mr_depth, rn_mr_width, rn_mr_lat_s               &
         &                 , rn_mr_lat_n, rn_ds_depth, rn_ds_width
      !!----------------------------------------------------------------------
      !
      ii = 1

      !REWIND( numnam_cfg )          ! Namelist namusr_def (exist in namelist_cfg only)
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      WRITE( numond, namusr_def )
      !
   END SUBROUTINE usr_def_nam_cfg


   FUNCTION merc_proj( pphi, pres ) RESULT( kphi )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE  ***
      !!
      !! ** Purpose :   Compute the number of gridpoints to the equator from a given latitude (MERCATOR).
      !!
      !! ** Method  :   Find number of gridpoints between the equator and the (approximate) latitude
      !!  pphi for the mercator projection. This is useful to find locations from the physical domain
      !!  on the numerical grid. We can also ensure that the equator is always on a U/T point.
      !!                
      !!                
      !!                
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )             ::   pphi     ! approximate latitude   [degrees]
      REAL(wp), INTENT(in   )             ::   pres     ! resolution of the grid [meters]
      INTEGER                             ::   kphi     ! number of gridpoints between equator and (approximately) pphi
      !!----------------------------------------------------------------------
      INTEGER  ::   jk
      REAL(wp) ::   zarg, zjeq, zsur     ! Computed parameters
      !!----------------------------------------------------------------------
      !
      zarg = rpi / 4. - rpi / 180. * pphi / 2.
      zjeq = ABS( 180./rpi * LOG( COS( zarg ) / SIN( zarg ) ) / pres )
      IF(  pphi > 0 )  zjeq = -zjeq
      zjeq =  zjeq + 1._wp          ! Fortran indexing starts at 1
      kphi = NINT( zjeq )
      !
   END FUNCTION merc_proj
   !!======================================================================
END MODULE usrdef_nam
