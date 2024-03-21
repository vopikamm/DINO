MODULE step_oce
   !!======================================================================
   !!                       ***  MODULE step_oce  ***
   !! Ocean time-stepping : module used in both initialisation phase and time stepping
   !!                                     (i.e. nemo_init and stp or stp_MLF routines)
   !!======================================================================
   !! History :   3.3  !  2010-08  (C. Ethe)  Original code - reorganisation of the initial phase
   !!             3.7  !  2014-01  (G. Madec) LDF simplication 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE domtile

   USE daymod          ! calendar                         (day     routine)

   USE sbc_oce         ! surface boundary condition: ocean
   USE sbcmod          ! surface boundary condition       (sbc     routine)
   USE sbcrnf          ! surface boundary condition: runoff variables
   USE sbccpl          ! surface boundary condition: coupled formulation (call send at end of step)
   USE sbcapr          ! surface boundary condition: atmospheric pressure
   USE sbcwave         ! Wave intialisation
   USE tide_mod        ! tides

   USE bdy_oce  , ONLY : ln_bdy
   USE bdydta          ! open boundary condition data     (bdy_dta routine)
   USE bdytra          ! bdy cond. for tracers            (bdy_tra routine)
   USE bdydyn3d        ! bdy cond. for baroclinic vel.  (bdy_dyn3d routine)

   USE isf_oce         ! ice shelf boundary condition
   USE isfstp          ! ice shelf boundary condition     (isf_stp routine)

   USE sshwzv          ! vertical velocity and ssh        (ssh_nxt routine)
   !                                                      (ssh_swp routine)
   !                                                      (wzv     routine)
   USE domvvl          ! variable vertical scale factors  (dom_vvl_sf_nxt routine)
   !                                                      (dom_vvl_sf_swp routine)
   
   USE divhor          ! horizontal divergence            (div_hor routine)
   USE dynadv          ! advection                        (dyn_adv routine)
   USE dynvor          ! vorticity term                   (dyn_vor routine)
   USE dynhpg          ! hydrostatic pressure grad.       (dyn_hpg routine)
   USE dynldf          ! lateral momentum diffusion       (dyn_ldf routine)
   USE zanna_bolton    ! Zanna & Bolton parameterization  (ZB_2020 routine)
   USE dynzdf          ! vertical diffusion               (dyn_zdf routine)
   USE dynspg          ! surface pressure gradient        (dyn_spg routine)
   USE dynatf          ! time-filtering                   (dyn_atf routine)
   USE dyndmp          ! current damping                  (dyn_dmp routine)

   USE traqsr          ! solar radiation penetration      (tra_qsr routine)
   USE traisf          ! ice shelf                        (tra_isf routine)
   USE trasbc          ! surface boundary condition       (tra_sbc routine)
   USE trabbc          ! bottom boundary condition        (tra_bbc routine)
   USE trabbl          ! bottom boundary layer            (tra_bbl routine)
   USE tradmp          ! internal damping                 (tra_dmp routine)
   USE traadv          ! advection scheme control     (tra_adv_ctl routine)
   USE traldf          ! lateral mixing                   (tra_ldf routine)
   USE trazdf          ! vertical mixing                  (tra_zdf routine)
   USE traatf          ! time filtering                   (tra_atf routine)
   USE tranpc          ! non-penetrative convection       (tra_npc routine)

   USE eosbn2          ! equation of state                (eos_bn2 routine)

   USE stopar          ! Stochastic parametrization       (sto_par routine)
   USE stopts 

   USE ldfslp          ! iso-neutral slopes               (ldf_slp routine)
   USE ldfdyn          ! lateral eddy viscosity coef.     (ldf_dyn routine)
   USE ldftra          ! lateral eddy diffusive coef.     (ldf_tra routine)

   USE zdf_oce         ! ocean vertical physics variables
   USE zdfphy          ! vertical physics manager      (zdf_phy_init routine)
   USE zdfdrg   , ONLY : ln_drgimp   ! implicit top/bottom friction
   USE zdfosm   , ONLY : osm_rst, dyn_osm, tra_osm      ! OSMOSIS routines used in step.F90
   USE zdfmfc          ! Mass FLux Convection routine used in step.F90

   USE diu_layers      ! diurnal SST bulk and coolskin routines
   USE sbc_oce         ! surface fluxes  
   
   USE zpshde          ! partial step: hor. derivative     (zps_hde routine)

   USE diawri          ! Standard run outputs             (dia_wri routine)
   USE diaptr          ! poleward transports              (dia_ptr routine)
   USE diadct          ! sections transports              (dia_dct routine)
   USE diaar5          ! AR5 diagnosics                   (dia_ar5 routine)
   USE diahth          ! thermocline depth                (dia_hth routine)
   USE diahsb          ! heat, salt and volume budgets    (dia_hsb routine)
   USE diacfl          ! CFL diagnostics                  (dia_cfl routine)
   USE diaobs          ! Observation operator             (dia_obs routine)
   USE diadetide       ! Weights computation for daily detiding of model diagnostics
   USE diamlr          ! IOM context management for multiple-linear-regression analysis
   USE flo_oce         ! floats variables
   USE floats          ! floats computation               (flo_stp routine)

   USE crsfld          ! Standard output on coarse grid   (crs_fld routine)

   USE asminc          ! assimilation increments      (tra_asm_inc routine)
   !                                                   (dyn_asm_inc routine)
   USE asmbkg          ! writing out state trajectory
   USE stpctl          ! time stepping control            (stp_ctl routine)
   USE restart         ! ocean restart                    (rst_wri routine)
   USE prtctl          ! Print control                    (prt_ctl routine)

   USE in_out_manager  ! I/O manager
   USE iom             !
   USE lbclnk
   USE timing          ! Timing

#if defined key_xios
   USE xios            ! I/O server
#endif
#if defined key_agrif
   USE agrif_oce_sponge ! Momemtum and tracers sponges
   USE agrif_all_update ! Main update driver
   USE agrif_oce_update
#endif
#if defined key_top
   USE trcstp, ONLY : trc_stp    ! passive tracer time-stepping      (trc_stp routine)
#endif
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: step_oce.F90 15023 2021-06-18 14:35:25Z gsamson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
END MODULE step_oce
