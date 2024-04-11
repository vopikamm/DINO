MODULE step
   !!======================================================================
   !!                       ***  MODULE step  ***
   !! Time-stepping   : manager of the ocean, tracer and ice time stepping
   !!======================================================================
   !! History :  OPA  !  1991-03  (G. Madec)  Original code
   !!             -   !  1991-11  (G. Madec)
   !!             -   !  1992-06  (M. Imbard)  add a first output record
   !!             -   !  1996-04  (G. Madec)  introduction of dynspg
   !!             -   !  1996-04  (M.A. Foujols)  introduction of passive tracer
   !!            8.0  !  1997-06  (G. Madec)  new architecture of call
   !!            8.2  !  1997-06  (G. Madec, M. Imbard, G. Roullet)  free surface
   !!             -   !  1999-02  (G. Madec, N. Grima)  hpg implicit
   !!             -   !  2000-07  (J-M Molines, M. Imbard)  Open Bondary Conditions
   !!   NEMO     1.0  !  2002-06  (G. Madec)  free form, suppress macro-tasking
   !!             -   !  2004-08  (C. Talandier) New trends organization
   !!             -   !  2005-01  (C. Ethe) Add the KPP closure scheme
   !!             -   !  2005-11  (G. Madec)  Reorganisation of tra and dyn calls
   !!             -   !  2006-01  (L. Debreu, C. Mazauric)  Agrif implementation
   !!             -   !  2006-07  (S. Masson)  restart using iom
   !!            3.2  !  2009-02  (G. Madec, R. Benshila)  reintroduicing z*-coordinate
   !!             -   !  2009-06  (S. Masson, G. Madec)  TKE restart compatible with key_cpl
   !!            3.3  !  2010-05  (K. Mogensen, A. Weaver, M. Martin, D. Lea) Assimilation interface
   !!             -   !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
   !!            3.4  !  2011-04  (G. Madec, C. Ethe) Merge of dtatem and dtasal
   !!            3.6  !  2012-07  (J. Simeon, G. Madec. C. Ethe)  Online coarsening of outputs
   !!            3.6  !  2014-04  (F. Roquet, G. Madec) New equations of state
   !!            3.6  !  2014-10  (E. Clementi, P. Oddo) Add Qiao vertical mixing in case of waves
   !!            3.7  !  2014-10  (G. Madec)  LDF simplication
   !!             -   !  2014-12  (G. Madec) remove KPP scheme
   !!             -   !  2015-11  (J. Chanut) free surface simplification (remove filtered free surface)
   !!            4.0  !  2017-05  (G. Madec)  introduction of the vertical physics manager (zdfphy)
   !!            4.1  !  2019-08  (A. Coward, D. Storkey) rewrite in preparation for new timestepping scheme
   !!----------------------------------------------------------------------

#if defined key_qco   ||   defined key_linssh
   !!----------------------------------------------------------------------
   !!   'key_qco'      EMPTY MODULE      Quasi-Eulerian vertical coordinate
   !!                                OR
   !!   'key_linssh    EMPTY MODULE       Fixed in time vertical coordinate
   !!----------------------------------------------------------------------
#else
   !!----------------------------------------------------------------------
   !!   stp             : OCE system time-stepping
   !!----------------------------------------------------------------------
   USE step_oce         ! time stepping definition modules

   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp   ! called by nemogcm.F90

   !                                          !**  time level indices  **!
   INTEGER, PUBLIC ::   Nbb, Nnn, Naa, Nrhs   !: used by nemo_init

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: step.F90 15398 2021-10-19 08:49:42Z timgraham $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_agrif
   RECURSIVE SUBROUTINE stp( )
      INTEGER             ::   kstp   ! ocean time-step index
#else
   SUBROUTINE stp( kstp )
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
#endif
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp  ***
      !!
      !! ** Purpose : - Time stepping of OCE  (momentum and active tracer eqs.)
      !!              - Time stepping of SI3 (dynamic and thermodynamic eqs.)
      !!              - Time stepping of TRC  (passive tracer eqs.)
      !!
      !! ** Method  : -1- Update forcings and data
      !!              -2- Update ocean physics
      !!              -3- Compute the t and s trends
      !!              -4- Update t and s
      !!              -5- Compute the momentum trends
      !!              -6- Update the horizontal velocity
      !!              -7- Compute the diagnostics variables (rd,N2, hdiv,w)
      !!              -8- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk, jtile   ! dummy loop indice
      !! ---------------------------------------------------------------------
#if defined key_agrif
      IF( nstop > 0 ) RETURN   ! avoid to go further if an error was detected during previous time step (child grid)
      kstp = nit000 + Agrif_Nb_Step()
      Kbb_a = Nbb; Kmm_a = Nnn; Krhs_a = Nrhs   ! agrif_oce module copies of time level indices
      IF( lk_agrif_debug ) THEN
         IF( Agrif_Root() .and. lwp)   WRITE(*,*) '---'
         IF(lwp)   WRITE(*,*) 'Grid Number', Agrif_Fixed(),' time step ', kstp, 'int tstep', Agrif_NbStepint()
      ENDIF
      IF( kstp == nit000 + 1 )   lk_agrif_fstep = .FALSE.
# if defined key_xios
      IF( Agrif_Nbstepint() == 0 )   CALL iom_swap( cxios_context )
# endif
#endif
      !
      IF( ln_timing )   CALL timing_start('stp')
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! model timestep
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      IF( l_1st_euler ) THEN     ! start or restart with Euler 1st time-step
         rDt   = rn_Dt   
         r1_Dt = 1._wp / rDt
      ENDIF
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! update I/O and calendar
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      IF( kstp == nit000 ) THEN                       ! initialize IOM context (must be done after nemo_init for AGRIF+XIOS+OASIS)
                             CALL iom_init( cxios_context, ld_closedef=.FALSE. )   ! for model grid (including possible AGRIF zoom)
         IF( lk_diamlr   )   CALL dia_mlr_iom_init    ! with additional setup for multiple-linear-regression analysis
                             CALL iom_init_closedef
         IF( ln_crs      )   CALL iom_init( TRIM(cxios_context)//"_crs" )  ! for coarse grid
      ENDIF
      IF( kstp == nitrst .AND. lwxios ) THEN
                             CALL iom_swap(                     cw_ocerst_cxt )
                             CALL iom_init_closedef(            cw_ocerst_cxt )
                             CALL iom_setkt( kstp - nit000 + 1, cw_ocerst_cxt )
#if defined key_top
                             CALL iom_swap(                     cw_toprst_cxt )
                             CALL iom_init_closedef(            cw_toprst_cxt )
                             CALL iom_setkt( kstp - nit000 + 1, cw_toprst_cxt )
#endif
      ENDIF
      IF( kstp + nn_fsbc - 1 == nitrst .AND. lwxios ) THEN
#if defined key_si3
                             CALL iom_swap(                     cw_icerst_cxt )
                             CALL iom_init_closedef(            cw_icerst_cxt )
                             CALL iom_setkt( kstp - nit000 + 1, cw_icerst_cxt )
#endif
         IF( ln_abl      ) THEN
                             CALL iom_swap(                     cw_ablrst_cxt )
                             CALL iom_init_closedef(            cw_ablrst_cxt )
                             CALL iom_setkt( kstp - nit000 + 1, cw_ablrst_cxt )
         ENDIF
      ENDIF
      IF( kstp /= nit000 )   CALL day( kstp )         ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp - nit000 + 1,      cxios_context          )   ! tell IOM we are at time step kstp
      IF( ln_crs         )   CALL iom_setkt( kstp - nit000 + 1, TRIM(cxios_context)//"_crs" )   ! tell IOM we are at time step kstp

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update external forcing (tides, open boundaries, ice shelf interaction and surface boundary condition (including sea-ice)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_tide    )   CALL tide_update( kstp )                     ! update tide potential
      IF( ln_apr_dyn )   CALL sbc_apr ( kstp )                        ! atmospheric pressure (NB: call before bdy_dta which needs ssh_ib)
      IF( ln_bdy     )   CALL bdy_dta ( kstp, Nnn )                   ! update dynamic & tracer data at open boundaries
      IF( ln_isf     )   CALL isf_stp ( kstp, Nnn )
                         CALL sbc     ( kstp, Nbb, Nnn )              ! Sea Boundary Condition (including sea-ice)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update stochastic parameters and random T/S fluctuations
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF( ln_sto_eos )   CALL sto_par( kstp )                         ! Stochastic parameters
      IF( ln_sto_eos )   CALL sto_pts( ts(:,:,:,:,Nnn)  )             ! Random T/S fluctuations

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Ocean physics update
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  THERMODYNAMICS
                         CALL eos_rab( ts(:,:,:,:,Nbb), rab_b, Nnn )       ! before local thermal/haline expension ratio at T-points
                         CALL eos_rab( ts(:,:,:,:,Nnn), rab_n, Nnn )       ! now    local thermal/haline expension ratio at T-points
                         CALL bn2    ( ts(:,:,:,:,Nbb), rab_b, rn2b, Nnn ) ! before Brunt-Vaisala frequency
                         CALL bn2    ( ts(:,:,:,:,Nnn), rab_n, rn2, Nnn  ) ! now    Brunt-Vaisala frequency

      !  VERTICAL PHYSICS
      ! lbc_lnk needed for zdf_sh2 when using nn_hls = 2, moved here to allow tiling in zdf_phy
      IF( nn_hls == 2 .AND. l_zdfsh2 ) CALL lbc_lnk( 'stp', avm_k, 'W', 1.0_wp )

      IF( ln_tile ) CALL dom_tile_start         ! [tiling] ZDF tiling loop
      DO jtile = 1, nijtile
         IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )

                         CALL zdf_phy( kstp, Nbb, Nnn, Nrhs )   ! vertical physics update (top/bot drag, avt, avs, avm + MLD)
      END DO
      IF( ln_tile ) CALL dom_tile_stop

      !  LATERAL  PHYSICS
      !
      IF( ln_zps .OR. l_ldfslp ) CALL eos( ts(:,:,:,:,Nbb), rhd, gdept_0(:,:,:) )               ! before in situ density

      IF( ln_zps .AND. .NOT. ln_isfcav)                                    &
            &            CALL zps_hde    ( kstp, jpts, ts(:,:,:,:,Nbb), gtsu, gtsv,  &  ! Partial steps: before horizontal gradient
            &                                          rhd, gru , grv    )       ! of t, s, rd at the last ocean level

      IF( ln_zps .AND.       ln_isfcav)                                                &
            &            CALL zps_hde_isf( kstp, jpts, ts(:,:,:,:,Nbb), gtsu, gtsv, gtui, gtvi,  &  ! Partial steps for top cell (ISF)
            &                                          rhd, gru , grv , grui, grvi   )       ! of t, s, rd at the first ocean level

      IF( l_ldfslp ) THEN                             ! slope of lateral mixing
         IF( ln_traldf_triad ) THEN
                         CALL ldf_slp_triad( kstp, Nbb, Nnn )             ! before slope for triad operator
         ELSE
                         CALL ldf_slp     ( kstp, rhd, rn2b, Nbb, Nnn )   ! before slope for standard operator
         ENDIF
      ENDIF
      !                                                                        ! eddy diffusivity coeff.
      IF( l_ldftra_time .OR. l_ldfeiv_time )   CALL ldf_tra( kstp, Nbb, Nnn )  !       and/or eiv coeff.
      IF( l_ldfdyn_time                    )   CALL ldf_dyn( kstp, Nbb )       ! eddy viscosity coeff.

      !  BBL coefficients
      !
      IF( ln_trabbl )   CALL bbl( kstp, nit000, Nbb, Nnn )   ! BBL diffusion coefficients and transports

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  Ocean dynamics : hdiv, ssh, e3, u, v, w
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                         CALL ssh_nxt       ( kstp, Nbb, Nnn, ssh, Naa )   ! after ssh (includes call to div_hor)
      IF( .NOT.ln_linssh )   &
                       & CALL dom_vvl_sf_nxt( kstp, Nbb, Nnn,      Naa )   ! after vertical scale factors
                         CALL wzv           ( kstp, Nbb, Nnn, Naa, ww  )   ! now cross-level velocity
      IF( ln_zad_Aimp )  CALL wAimp         ( kstp,      Nnn           )   ! Adaptive-implicit vertical advection partitioning
                         CALL eos    ( ts(:,:,:,:,Nnn), rhd, rhop, gdept(:,:,:,Nnn) )  ! now in situ density for hpg computation


                         uu(:,:,:,Nrhs) = 0._wp            ! set dynamics trends to zero
                         vv(:,:,:,Nrhs) = 0._wp

      IF( ln_tile ) CALL dom_tile_start         ! [tiling] DYN tiling loop (1)
      DO jtile = 1, nijtile
         IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )

         IF(  lk_asminc .AND. ln_asmiau .AND. ln_dyninc )   &
                  &         CALL dyn_asm_inc   ( kstp, Nbb, Nnn, uu, vv, Nrhs )  ! apply dynamics assimilation increment
         IF( ln_bkgwri )    CALL asm_bkg_wri( kstp, Nnn )     ! output background fields
         IF( ln_bdy     )   CALL bdy_dyn3d_dmp ( kstp, Nbb,      uu, vv, Nrhs )  ! bdy damping trends
#if defined key_agrif
      END DO
      IF( ln_tile ) CALL dom_tile_stop

      IF(.NOT. Agrif_Root())  &
               &         CALL Agrif_Sponge_dyn        ! momentum sponge

      IF( ln_tile ) CALL dom_tile_start         ! [tiling] DYN tiling loop (1, continued)
      DO jtile = 1, nijtile
         IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )
#endif
                                 CALL dyn_adv( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! advection (VF or FF)	==> RHS
                                 CALL dyn_vor( kstp,      Nnn      , uu, vv, Nrhs )  ! vorticity           	==> RHS
                                 CALL dyn_ldf( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! lateral mixing
         IF( ln_zanna_bolton )   CALL ZB_apply( kstp,Nbb      , uu, vv, Nrhs )  ! Zanna & Bolton 2020 parameterization                   
         IF( ln_zdfosm  )        CALL dyn_osm( kstp,      Nnn      , uu, vv, Nrhs )  ! OSMOSIS non-local velocity fluxes ==> RHS
                                 CALL dyn_hpg( kstp,      Nnn      , uu, vv, Nrhs )  ! horizontal gradient of Hydrostatic pressure
      END DO
      IF( ln_tile ) CALL dom_tile_stop

                            CALL dyn_spg( kstp, Nbb, Nnn, Nrhs, uu, vv, ssh, uu_b, vv_b, Naa )  ! surface pressure gradient

                                                      ! With split-explicit free surface, since now transports have been updated and ssh(:,:,Nrhs) as well
      IF( ln_dynspg_ts ) THEN                         ! vertical scale factors and vertical velocity need to be updated
         IF( ln_tile ) CALL dom_tile_start      ! [tiling] DYN tiling loop (2- div_hor only)
         DO jtile = 1, nijtile
            IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )

                             CALL div_hor       ( kstp, Nbb, Nnn )               ! Horizontal divergence  (2nd call in time-split case)
         END DO
         IF( ln_tile ) CALL dom_tile_stop

         IF(.NOT. ln_linssh) CALL dom_vvl_sf_nxt( kstp, Nbb, Nnn, Naa, kcall=2 )  ! after vertical scale factors (update depth average component)
      ENDIF

      IF( ln_tile ) CALL dom_tile_start         ! [tiling] DYN tiling loop (3- dyn_zdf only)
      DO jtile = 1, nijtile
         IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )

                               CALL dyn_zdf    ( kstp, Nbb, Nnn, Nrhs, uu, vv, Naa  )  ! vertical diffusion
      END DO
      IF( ln_tile ) CALL dom_tile_stop

      IF( ln_dynspg_ts ) THEN                                                       ! vertical scale factors and vertical velocity need to be updated
                            CALL wzv        ( kstp, Nbb, Nnn, Naa, ww )             ! Nnn cross-level velocity
         IF( ln_zad_Aimp )  CALL wAimp      ( kstp,      Nnn )                      ! Adaptive-implicit vertical advection partitioning
      ENDIF


      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! cool skin
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF ( ln_diurnal )  CALL diurnal_layers( kstp )

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! diagnostics and outputs
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_floats  )   CALL flo_stp   ( kstp, Nbb, Nnn )      ! drifting Floats
      IF( ln_diacfl  )   CALL dia_cfl   ( kstp,      Nnn )      ! Courant number diagnostics
                         CALL dia_hth   ( kstp,      Nnn )      ! Thermocline depth (20 degres isotherm depth)
      IF( ln_diadct  )   CALL dia_dct   ( kstp,      Nnn )      ! Transports
                         CALL dia_ar5   ( kstp,      Nnn )      ! ar5 diag
                         CALL dia_ptr   ( kstp,      Nnn )      ! Poleward adv/ldf TRansports diagnostics
                         CALL dia_wri   ( kstp,      Nnn )      ! ocean model: outputs
      IF( ln_crs     )   CALL crs_fld   ( kstp,      Nnn )      ! ocean model: online field coarsening & output
      IF( lk_diadetide ) CALL dia_detide( kstp )                ! Weights computation for daily detiding of model diagnostics
      IF( lk_diamlr  )   CALL dia_mlr                           ! Update time used in multiple-linear-regression analysis

#if defined key_top
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Passive Tracer Model
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL trc_stp    ( kstp, Nbb, Nnn, Nrhs, Naa )           ! time-stepping
#endif

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Active tracers
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         ts(:,:,:,:,Nrhs) = 0._wp         ! set tracer trends to zero

      IF( ln_tile ) CALL dom_tile_start         ! [tiling] TRA tiling loop (1)
      DO jtile = 1, nijtile
         IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )

         IF(  lk_asminc .AND. ln_asmiau .AND. &
            & ln_trainc )   CALL tra_asm_inc( kstp, Nbb, Nnn, ts, Nrhs )  ! apply tracer assimilation increment
                            CALL tra_sbc    ( kstp,      Nnn, ts, Nrhs )  ! surface boundary condition
         IF( ln_traqsr  )   CALL tra_qsr    ( kstp,      Nnn, ts, Nrhs )  ! penetrative solar radiation qsr
         IF( ln_isf     )   CALL tra_isf    ( kstp,      Nnn, ts, Nrhs )  ! ice shelf heat flux
         IF( ln_trabbc  )   CALL tra_bbc    ( kstp,      Nnn, ts, Nrhs )  ! bottom heat flux
         IF( ln_trabbl  )   CALL tra_bbl    ( kstp, Nbb, Nnn, ts, Nrhs )  ! advective (and/or diffusive) bottom boundary layer scheme
         IF( ln_tradmp  )   CALL tra_dmp    ( kstp, Nbb, Nnn, ts, Nrhs )  ! internal damping trends
         IF( ln_bdy     )   CALL bdy_tra_dmp( kstp, Nbb,      ts, Nrhs )  ! bdy damping trends
      END DO
      IF( ln_tile ) CALL dom_tile_stop

#if defined key_agrif
      IF(.NOT. Agrif_Root() )   THEN
                            CALL Agrif_Sponge_tra        ! tracers sponge
      ENDIF
#endif

      ! TEMP: [tiling] Separate loop over tile domains (due to tra_adv workarounds for tiling)
      IF( ln_tile ) CALL dom_tile_start         ! [tiling] TRA tiling loop (2)
      DO jtile = 1, nijtile
         IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )

                            CALL tra_adv    ( kstp, Nbb, Nnn, ts, Nrhs )  ! hor. + vert. advection	==> RHS
         IF( ln_zdfmfc  )   CALL tra_mfc    ( kstp, Nbb,      ts, Nrhs )  ! Mass Flux Convection
         IF( ln_zdfosm  ) THEN
                            CALL tra_osm    ( kstp,      Nnn, ts, Nrhs )  ! OSMOSIS non-local tracer fluxes ==> RHS
            IF( lrst_oce )  CALL osm_rst    ( kstp,      Nnn, 'WRITE'  )  ! write OSMOSIS outputs + ww (so must do here) to restarts
         ENDIF
                            CALL tra_ldf    ( kstp, Nbb, Nnn, ts, Nrhs )  ! lateral mixing

                            CALL tra_zdf    ( kstp, Nbb, Nnn, Nrhs, ts, Naa  )  ! vertical mixing and after tracer fields
         IF( ln_zdfnpc  )   CALL tra_npc    ( kstp,      Nnn, Nrhs, ts, Naa  )  ! update after fields by non-penetrative convection
      END DO
      IF( ln_tile ) CALL dom_tile_stop

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Set boundary conditions, time filter and swap time levels
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!jc1: For agrif, it would be much better to finalize tracers/momentum here (e.g. bdy conditions) and move the swap
!!    (and time filtering) after Agrif update. Then restart would be done after and would contain updated fields.
!!    If so:
!!    (i) no need to call agrif update at initialization time
!!    (ii) no need to update "before" fields
!!
!!    Apart from creating new tra_swp/dyn_swp routines, this however:
!!    (i) makes boundary conditions at initialization time computed from updated fields which is not the case between
!!    two restarts => restartability issue. One can circumvent this, maybe, by assuming "interface separation",
!!    e.g. a shift of the feedback interface inside child domain.
!!    (ii) requires that all restart outputs of updated variables by agrif (e.g. passive tracers/tke/barotropic arrays) are done at the same
!!    place.
!!
!!jc2: dynnxt must be the latest call. e3t(:,:,:,Nbb) are indeed updated in that routine
                         CALL tra_atf       ( kstp, Nbb, Nnn, Naa, ts )                      ! time filtering of "now" tracer arrays
                         CALL dyn_atf       ( kstp, Nbb, Nnn, Naa, uu, vv, e3t, e3u, e3v  )  ! time filtering of "now" velocities and scale factors
                         CALL ssh_atf       ( kstp, Nbb, Nnn, Naa, ssh )                     ! time filtering of "now" sea surface height
      !
      ! Swap time levels
      Nrhs = Nbb
      Nbb = Nnn
      Nnn = Naa
      Naa = Nrhs
      !
      IF(.NOT.ln_linssh) CALL dom_vvl_sf_update( kstp, Nbb, Nnn, Naa )  ! recompute vertical scale factors
      !
      IF( ln_diahsb  )   CALL dia_hsb       ( kstp, Nbb, Nnn )  ! - ML - global conservation diagnostics

!!gm : This does not only concern the dynamics ==>>> add a new title
!!gm2: why ouput restart before AGRIF update?
!!
!!jc: That would be better, but see comment above
!!
      IF( lrst_oce   )   CALL rst_write    ( kstp, Nbb, Nnn )   ! write output ocean restart file
      IF( ln_sto_eos )   CALL sto_rst_write( kstp )   ! write restart file for stochastic parameters

#if defined key_agrif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! AGRIF recursive integration
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         Kbb_a = Nbb; Kmm_a = Nnn; Krhs_a = Nrhs      ! agrif_oce module copies of time level indices
                         CALL Agrif_Integrate_ChildGrids( stp )       ! allows to finish all the Child Grids before updating

#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL stp_ctl      ( kstp, Nnn )

#if defined key_agrif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! AGRIF update
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( Agrif_NbStepint() == 0 .AND. nstop == 0 )   &
         &               CALL Agrif_update_all( )                  ! Update all components

#endif
      IF( ln_diaobs .AND. nstop == 0 )   &
         &               CALL dia_obs( kstp, Nnn )  ! obs-minus-model (assimilation) diags (after dynamics update)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! File manipulation at the end of the first time step
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( kstp == nit000 ) THEN                          ! 1st time step only
                                        CALL iom_close( numror )   ! close input  ocean restart file
         IF( lrxios )                   CALL iom_context_finalize( cr_ocerst_cxt )
         IF(lwm)                        CALL FLUSH    ( numond )   ! flush output namelist oce
         IF(lwm .AND. numoni /= -1 )    CALL FLUSH    ( numoni )   ! flush output namelist ice (if exist)
      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Coupled mode
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( lk_oasis .AND. nstop == 0 )   CALL sbc_cpl_snd( kstp, Nbb, Nnn )     ! coupled mode : field exchanges
      !
#if defined key_xios
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Finalize contextes if end of simulation or error detected
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( kstp == nitend .OR. nstop > 0 ) THEN
                      CALL iom_context_finalize(      cxios_context          ) ! needed for XIOS+AGRIF
         IF( ln_crs ) CALL iom_context_finalize( trim(cxios_context)//"_crs" ) !
      ENDIF
#endif
      !
      IF( l_1st_euler ) THEN         ! recover Leap-frog timestep
         rDt   = 2._wp * rn_Dt
         r1_Dt = 1._wp / rDt
         l_1st_euler = .FALSE.
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('stp')
      !
   END SUBROUTINE stp
   !
#endif
   !!======================================================================
END MODULE step
