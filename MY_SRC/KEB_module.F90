MODULE KEB_module
   !!==============================================================================
   !!                       ***  MODULE  KEB_module  ***
   !!==============================================================================
   !! History : NEMO
   !!            3.6  ! 2019-05  (P.Perezhogin) original code for Cartesian coordinates
   !!                 (published in P. Perezhogin, RJNAMM, 2020)
   !!            3.6  ! 2019-10  (P. Perezhogin) curvilinear coordinates were added
   !!            3.6  ! 2021-06  (P. Perezhogin) preparing the code to Jones Colin
   !!----------------------------------------------------------------------

   USE oce                    ! ocean dynamics and tracers
   USE dom_oce                ! ocean space and time domain
   USE lib_mpp                ! MPP library
   USE lbclnk                 ! ocean lateral boundary conditions (or mpp link)
   USE iom                    ! input-output manager
   USE storng                 ! generation of gaussian noise
   USE in_out_manager
   USE bdy_par, only: lk_bdy  ! for Warning

   USE KEB_operators          ! spatial operators for KEB
   USE KEB_testing            ! check spatial operators
   
   IMPLICIT NONE   
   PRIVATE
   
   ! subroutines
   PUBLIC KEB_init               ! called by nemogcm.F90      | init KEB parameters and allocate arrays
   PUBLIC KEB_rst                ! called by step.F90         | writes restart file
   PUBLIC KEB_apply              ! called by step.F90         | apply KEB tendencies
   
   ! logical keys
   PUBLIC KEB_on
   PUBLIC KEB_test

   PUBLIC Ediss_check            ! defined in dynldf_bilap.F90 | direct computing of Ediss for testing
   PUBLIC Ediss                  ! defined in dynldf_bilap.F90 | get Ediss from horizontal visocosity operator
   
   !!------------- parameters -------------!!
   LOGICAL  :: KEB_on               ! key for KEB
   LOGICAL  :: KEB_test             ! test KEB operators
   LOGICAL  :: KEB_negvisc          ! key for negative viscosity KEB
   LOGICAL  :: KEB_AR1              ! key for stochastic AR1 in time KEB

   LOGICAL  :: tke_adv, tke_diff    ! keys for advection and diffusion of TKE
   LOGICAL  :: dirichlet_TKE        ! true, if zero Dirichlet B.C., false if zero Neuman
   REAL(wp) :: nu_TKE               ! diffusion coef. for TKE, m^2/s
   LOGICAL  :: tke_filter           ! filtering in z direction to avoid numerical instability

   LOGICAL  :: constant_cdiss       ! true = constant cdiss, false = cdiss varies according to Klower
   REAL(wp) :: c_diss               ! constant backscatter rate, i.e. amount of backscattered energy, typically from 0 to 1
   REAL(wp) :: R_diss               ! constant for Klower formula
   INTEGER  :: ndiss                ! the number of filters to smooth energy tendency Ediss
   LOGICAL  :: dirichlet_filter     ! true, if zero Dirichlet B.C., false if zero Neuman


   REAL(wp) :: c_back               ! nondimensional parameter for determing neg. visc. from TKE, typical value given by Jansen 2015 is 0.4*sqrt(2.)
   INTEGER  :: nback                ! the number of filters to smooth negative viscosity backscatter tendency

   INTEGER  :: nstoch               ! the number of filters applied to stochastic streamfunction
   REAL(wp) :: T_decorr             ! decorrelation time in seconds, for KEB_AR1

   REAL(wp) :: Ediss_check          ! compare flux and advective forms of energy dissipation

   !! ------------ prognostic variables --------------!!
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: TKE               ! subgrid kinetic energy in T points
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: nu2t              ! predicted negative viscosity in T points
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: rhs_adv, rhs_diff ! RHS in TKE equation for advection and diffusion
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: Ediss, Eback ! Dissipation and backscatter of energy
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: Esource           ! Source of TKE (cdiss * Ediss)
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: local_cdiss       ! local value of cdiss

   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: rhsu, rhsv        ! negvisc KEB tendency, modifies momentum equation
   
   REAL(wp), ALLOCATABLE, DIMENSION (:,:)   :: phi               ! 2d random field
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: psi, psib         ! 3d random streamfunction, at now and before time steps
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: fx, fy, fxb, fyb  ! tendency of stochastic KEB, at now and before time steps

   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: ffmask            ! fmask with 0 on the boundary, without rn_shlat

   ! a posteriori correction for AR1 parameterization   
   REAL(wp) :: Estoch_mean  ! AR1 backscattered energy before correction, time-mean
   REAL(wp) :: Esource_mean ! Energy to be backscattered, time-mean
   
CONTAINS 

   SUBROUTINE KEB_init ( )
   
      INTEGER(KIND=8) :: zseed1, zseed2, zseed3, zseed4
      
      INTEGER :: ji, jj, jk
      
      NAMELIST/KEB_prms/ KEB_on, KEB_test, constant_cdiss, c_diss, R_diss, ndiss, dirichlet_filter, &
      KEB_negvisc, tke_adv, tke_diff, dirichlet_TKE, nu_TKE, tke_filter, c_back, nback, KEB_AR1, nstoch, T_decorr
      
      REWIND( numnam_cfg )
      READ  ( numnam_cfg, KEB_prms )
      
      ! if KEB is turned off, turn off all tendencies
      IF (not(KEB_on)) THEN
         KEB_negvisc = .false.
         KEB_AR1 = .false.
      END IF

      ! if all parameterizations off, turn off KEB
      IF (not(KEB_negvisc .or. KEB_AR1)) KEB_on = .false.

      IF (KEB_negvisc .and. KEB_AR1) CALL ctl_stop('KEB error: use only one parameterization')

      ! nback must be even and positive
      nback = nback / 2
      nback = max(nback * 2, 0)

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'KEB_init: kinetic energy backscatter parameterizations - initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
         IF (lk_bdy) WRITE(numout,*) 'KEB Warning: lk_bdy = .true. Check tmask in KEB_module.f90 before use, because lk_bdy key modifies it'
         IF (lk_vvl) WRITE(numout,*) 'KEB Warning: lk_vvl = .true. Check scale factors e3t_0, e3u_0, e3v_0, e3f_0 in KEB files' 

         WRITE(numout,*) '   Namelist KEB_prms : set KEB parameters'
         WRITE(numout,*) '      KEB is working                    KEB_on =', KEB_on
         WRITE(numout,*) '      testing KEB operators           KEB_test =', KEB_test
         
         WRITE(numout,*) ''

         WRITE(numout,*) '      constant cdiss            constant_cdiss =', constant_cdiss
         WRITE(numout,*) '      constant backscatter rate         c_diss =', REAL(c_diss,4)
         WRITE(numout,*) '      threshold Rossby number (Klower)  R_diss =', REAL(R_diss,4)
         WRITE(numout,*) '      number of filters for Ediss        ndiss =', ndiss
         WRITE(numout,*) '      filter b.c. for Ediss   dirichlet_filter =', dirichlet_filter

         WRITE(numout,*) ''

         WRITE(numout,*) '      negative viscosity KEB       KEB_negvisc =', KEB_negvisc
         WRITE(numout,*) '      advection of TKE                 tke_adv =', tke_adv
         WRITE(numout,*) '      diffusion of TKE                tke_diff =', tke_diff
         WRITE(numout,*) '      dirichlet B.C. for TKE     dirichlet_TKE =', dirichlet_TKE
         WRITE(numout,*) '      diffusion coefficient of TKE      nu_TKE =', REAL(nu_TKE,4)
         WRITE(numout,*) '      filter TKE in z direction     tke_filter =', tke_filter
         WRITE(numout,*) '      neg. visc. ~ c_back sqrt(TKE)     c_back =', REAL(c_back,4)
         WRITE(numout,*) '      number of filters for negvsc tend. nback =', nback                  

         WRITE(numout,*) ''
         
         WRITE(numout,*) '      autoregressive KEB               KEB_AR1 =', KEB_AR1
         WRITE(numout,*) '      number of filters for Psi         nstoch =', nstoch
         WRITE(numout,*) '      decor. time in sec. for AR1     T_decorr =', REAL(T_decorr,4)
      END IF
      
      IF (not(KEB_on)) THEN
         return
      END IF

      ALLOCATE(TKE (jpi,jpj,jpk))
      ALLOCATE(nu2t(jpi,jpj,jpk))
      ALLOCATE(rhs_adv (jpi,jpj,jpk)) 
      ALLOCATE(rhs_diff(jpi,jpj,jpk))
      ALLOCATE(Ediss(jpi,jpj,jpk))
      ALLOCATE(Eback(jpi,jpj,jpk))
      ALLOCATE(Esource(jpi,jpj,jpk))
      ALLOCATE(local_cdiss(jpi,jpj,jpk))
      
      ALLOCATE(rhsu(jpi,jpj,jpk)) 
      ALLOCATE(rhsv(jpi,jpj,jpk))
      
      ALLOCATE(phi(jpi,jpj)) 
      ALLOCATE(psi(jpi,jpj,jpk))
      ALLOCATE(psib(jpi,jpj,jpk))
      ALLOCATE(fx(jpi,jpj,jpk)) 
      ALLOCATE(fy(jpi,jpj,jpk))
      ALLOCATE(fxb(jpi,jpj,jpk)) 
      ALLOCATE(fyb(jpi,jpj,jpk))

      ALLOCATE(ffmask(jpi,jpj,jpk))

      ! without rn_shlat
      ffmask = 0._wp
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               ffmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
                                * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
            END DO
         END DO
      END DO
      CALL lbc_lnk( ffmask, 'F', 1._wp )

      ! init subgrid TKE with zeros
      TKE     = 0._wp
      nu2t    = 0._wp;
      rhs_adv = 0._wp; rhs_diff = 0._wp
      Ediss   = 0._wp; Eback = 0._wp
      Esource = 0._wp
      local_cdiss = 0._wp
      
      rhsu    = 0._wp; rhsv = 0._wp
      
      phi = 0._wp; 
      psi = 0._wp; psib = 0._wp
      fx = 0._wp; fxb = 0._wp
      fy = 0._wp; fyb = 0._wp
      
      ! typical values of energy exchange per unit volume in [m2/s3]
      ! needed to smoothly start a posteriori correction
      Estoch_mean  = 1.e-9
      Esource_mean = 1.e-9

      ! init seed for stochastic KEB
      CALL kiss_reset()                                 ! set default seeds
      CALL kiss_state( zseed1, zseed2, zseed3, zseed4 ) ! get default seeds
      zseed1 = zseed1 + int(mpprank,8)                  ! unique seed for each rank
      CALL kiss_seed ( zseed1, zseed2, zseed3, zseed4 ) 

      CALL KEB_rst_read()
   END SUBROUTINE KEB_init

   SUBROUTINE KEB_rst_read()
      INTEGER :: id1, id2, id3, id4, id5
      
      IF ( ln_rstart ) THEN
         id1 = iom_varid( numror, 'KEB_tke', ldstop = .FALSE.)
         IF (id1 > 0) THEN
            CALL iom_get( numror, jpdom_autoglo, 'KEB_tke', TKE)
            WRITE(numout,*) '~~~~~~~~~~~~'
            IF (lwp) WRITE(numout,*) '      KEB: TKE set from restart file'
            WRITE(numout,*) '~~~~~~~~~~~~'
         ELSE
            WRITE(numout,*) '~~~~~~~~~~~~'
            IF (lwp) WRITE(numout,*) '      KEB: TKE set to zero'
            WRITE(numout,*) '~~~~~~~~~~~~'
         END IF

         id1 = iom_varid( numror, 'KEB_fx', ldstop = .FALSE.)
         id2 = iom_varid( numror, 'KEB_fy', ldstop = .FALSE.)
         id3 = iom_varid( numror, 'KEB_psi', ldstop = .FALSE.)
         id4 = iom_varid( numror, 'KEB_source', ldstop = .FALSE.)
         id5 = iom_varid( numror, 'KEB_stoch', ldstop = .FALSE.)
         IF (id1 > 0 .and. id2 > 0 .and. id3 > 0 .and. id4 > 0 .and. id5 > 0) THEN
            CALL iom_get( numror, jpdom_autoglo, 'KEB_fx', fx)
            CALL iom_get( numror, jpdom_autoglo, 'KEB_fy', fy)
            CALL iom_get( numror, jpdom_autoglo, 'KEB_psi', psi)
            CALL iom_get( numror, 'KEB_source', Esource_mean)
            CALL iom_get( numror, 'KEB_stoch', Estoch_mean)
            WRITE(numout,*) '~~~~~~~~~~~~'
            IF (lwp) WRITE(numout,*) '      KEB: fx, fy, psi, Estoch_mean, Esource_mean set from restart file'
            WRITE(numout,*) '~~~~~~~~~~~~'
         ELSE
            WRITE(numout,*) '~~~~~~~~~~~~'
            IF (lwp) WRITE(numout,*) '      KEB: fx, fy, psi, Estoch_mean, Esource_mean set to zero'
            WRITE(numout,*) '~~~~~~~~~~~~'
         END IF
      END IF

      ! no need for MPI exchange because iom_get do it

   END SUBROUTINE KEB_rst_read

   SUBROUTINE KEB_rst ( kt )
      INTEGER :: kt ! ocean time step

      IF (lwp) WRITE(numout,*) '---- KEB-rst ----'
      CALL iom_rstput( kt, nitrst, numrow, 'KEB_tke', TKE)
      CALL iom_rstput( kt, nitrst, numrow, 'KEB_fx', fx)
      CALL iom_rstput( kt, nitrst, numrow, 'KEB_fy', fy)
      CALL iom_rstput( kt, nitrst, numrow, 'KEB_psi', psi)
      CALL iom_rstput( kt, nitrst, numrow, 'KEB_source', Esource_mean)
      CALL iom_rstput( kt, nitrst, numrow, 'KEB_stoch', Estoch_mean)

   END SUBROUTINE KEB_rst

   SUBROUTINE KEB_apply( kt, Kbb, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE KEB_apply  ***   
      !! ** Purpose :   Multiplies Ediss by cdiss and applies KEB tendences
      !!                Additionally, statistics are collected
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt               ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kbb, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv 
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zurhs, zvrhs

      IF (constant_cdiss) THEN
         local_cdiss = c_diss
      ELSE
         CALL Klower_cdiss( un, vn, Ediss, R_diss, local_cdiss, ffmask )
      END IF
      Esource = Ediss * local_cdiss ! Ediss is computed in dynldf_bilap.F90

      CALL lbc_lnk(Esource, 'T', 1.)
      CALL filter_laplace_T3D_ntimes(Esource, Esource, ndiss, dirichlet_filter)
      
      ! produce backscatter given Esource field
      IF (KEB_negvisc) CALL KEB_negvisc_tendency(kt, Kbb, puu(:,:,:,Kbb), pvv(:,:,:,Kbb), zurhs, zvrhs)
      IF (KEB_AR1    ) CALL KEB_AR1_tendency(kt, Kbb, puu(:,:,:,Kbb), pvv(:,:,:,Kbb), zurhs, zvrhs)
      
      CALL KEB_statistics( )

      puu(:,:,:,Krhs) = puu(:,:,:,Krhs) + zurhs
      pvv(:,:,:,Krhs) = pvv(:,:,:,Krhs) + zvrhs
      !
   END SUBROUTINE KEB_apply
   
   SUBROUTINE KEB_negvisc_tendency(kt, Kbb, puu, pvv, kebu, kebv)      
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE KEB_negvisc_tendency  ***   
      !! ** Purpose :   Computes negative viscosity momentum tendency
      !!                Computes Eback, updates TKE
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step index
      INTEGER, INTENT(in) ::   Kbb  ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  puu, pvv         ! before velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  kebu, kebv       ! u-, v-component of KEB tendency

      REAL(wp), DIMENSION(jpi,jpj,jpk) :: rot       ! u-, v-component of KEB tendency 
      
      INTEGER  ::   ji, jj, jk, n           ! dummy loop indices
      
      ! compute negative viscosity and rotation at now time step
      DO jk = 1, jpkm1
         nu2t(:,:,jk) = sqrt(2._wp) * c_back * sqrt(e12t(:,:) * max(TKE(:,:,jk),0._wp)) ! lbc_lnk is already applied to TKE

         DO_2D( nn_hls-1, nn_hls, nn_hls-1, nn_hls )
            !                                      ! (warning: computed for ji-1,jj-1); dk: see dynldf_lap_blp
            rot(ji-1,jj-1,jk) = r1_e1e2f(ji-1,jj-1) * ffmask(ji-1, jj-1,jk)       &
               &     * (  e2v(ji  ,jj-1) * pvv(ji  ,jj-1,jk) - e2v(ji-1,jj-1) * pvv(ji-1,jj-1,jk)  &
               &        - e1u(ji-1,jj  ) * puu(ji-1,jj  ,jk) + e1u(ji-1,jj-1) * puu(ji-1,jj-1,jk)  )
         END_2D
      END DO

      CALL KEB_ldf_lap( rot, nu2t, nback, Eback, kebu, kebv, ffmask ) 

      CALL iom_put('negviscx', kebu)
      CALL iom_put('negviscy', kebv)
      
      ! update TKE
      IF (tke_adv)   CALL upwind_advection ( TKE, un, vn, wn, rhs_adv, .true. ) ! .true. = free surface b.c.
      IF (tke_diff)  CALL laplace_T3D( TKE, rhs_diff, nu_TKE, dirichlet_TKE ) 

      ! IF (KEB_test) CALL test_negvisc_KEB( TKE, rhs_adv, rhs_diff, Esource, Eback, Ediss, &
      !                                      local_cdiss, ffmask, rhsu, rhsv, Ediss_check )

      ! update subgrid energy to after time step
      TKE = TKE + tmask * (rhs_adv + rhs_diff + Esource - Eback) * rdt

      IF (tke_filter) CALL z_filter(TKE, TKE)
      
      CALL lbc_lnk( TKE, 'T', 1.)
      
   END SUBROUTINE KEB_negvisc_tendency
   
   SUBROUTINE KEB_AR1_tendency(kt, Kbb, puu, pvv, kebu, kebv)   
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE KEB_AR1_tendency  ***   
      !! ** Purpose :   Computes stochastic streamfunction given Esource,
      !!                Computes Eback, updates momentum tendency
      !! ** Method  :   Stochastic streamfunction construction:
      !!                1) generate 2d white-noise field (phi)
      !!                2) weighting with energy source and spatial filtering
      !!                3) time-filtering with AR1 process 
      !!                4) compute curl
      !!                5) a posteriori correction of energy input
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step index
      INTEGER, INTENT(in) ::   Kbb  ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  puu, pvv         ! before velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  kebu, kebv       ! u-, v-component of KEB tendency

      REAL(wp) :: amp_increase
      
      INTEGER  ::   ji, jj, jk, n           ! dummy loop indices

      ! swap arrays
      fxb = fx
      fyb = fy
      psib = psi
      
      ! N(0,1) 2d noise
      CALL gauss_white_noise_2d( phi )

      ! weight with energy source and correct amplitude
      CALL compute_psi( Esource, phi, nstoch, T_decorr, psi, ffmask )

      ! correlation in time, doesn't change variance
      psi = psib * (1._wp - rdt / T_decorr) + sqrt(rdt / T_decorr * (2._wp - rdt / T_decorr)) * psi

      ! convert streamfunction to velocity tendency
      CALL horizontal_curl( psi, fx, fy )

      ! IF (KEB_test) CALL test_AR1_KEB( Esource, Ediss, local_cdiss, ffmask, fx, fy, Ediss_check )   

      ! a posteriori correction of energy input
      CALL compute_Eback_AR1( fx, fy, puu, pvv, Eback )
      amp_increase = aposteriori_correction( Esource, Eback, Esource_mean, Estoch_mean)

      fx = fx * amp_increase
      fy = fy * amp_increase
      Eback = Eback * amp_increase

      ! modified LF scheme (see NEMO book); dk: don't know if this still makes sense in 4.2
      kebu =  (fx + fxb) * 0.5_wp
      kebv =  (fy + fyb) * 0.5_wp   

      CALL iom_put('SKEBamp', amp_increase)
      CALL iom_put('SKEBpsi', psi)
      CALL iom_put('SKEBpsi_s', psi(:,:,1))
      CALL iom_put('SKEBfx' , fx)
      CALL iom_put('SKEBfy' , fy)

   END SUBROUTINE KEB_AR1_tendency
   
   SUBROUTINE KEB_statistics( )
      
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: z3d
      REAL(wp), SAVE :: Eback_sum = 0._wp, Ediss_sum = 0._wp
      REAL(wp) :: Eback_m, TKE_m
      REAL(wp) :: min_TKE, time_TKE
      REAL(wp) :: average_cdiss

      CALL put_fields('nu2t', nu2t)
      CALL put_fields('TKE', TKE)
      CALL put_fields('Ediss', Ediss)
      CALL put_fields('Eback', Eback)
      CALL put_fields('Esource', Esource)
      CALL put_fields('local_cdiss', local_cdiss)

      IF (iom_use('min_TKE')) THEN
         min_TKE = min_xyz(TKE, tmask)
         CALL iom_put('min_TKE', min_TKE)
      END IF

      IF (iom_use('average_cdiss')) THEN
         Eback_sum = Eback_sum + average_xyz(Eback,e1t,e2t,e3t_0,tmask)
         Ediss_sum = Ediss_sum + average_xyz(Ediss,e1t,e2t,e3t_0,tmask)
         average_cdiss = 0._wp
         IF (Ediss_sum > 1.e-16) THEN
            average_cdiss = Eback_sum / Ediss_sum
         END IF
         CALL iom_put('average_cdiss', average_cdiss)
      END IF

      ! time scale of returning energy in days
      IF (iom_use('time_TKE')) THEN
         Eback_m = average_xyz(Eback,e1t,e2t,e3t_0,tmask)
         TKE_m   = average_xyz(TKE,e1t,e2t,e3t_0,tmask)
         time_TKE = 0._wp
         IF (Eback_m > 1.e-16) THEN
            time_TKE = TKE_m / Eback_m / 86400._wp
         END IF
         CALL iom_put('time_TKE', time_TKE)
      END IF

      ! IF (iom_use('est_TKE') .OR. iom_use('est_TKE_s') .OR. &
      !     iom_use('est_TKE_z') .OR. iom_use('est_TKE_0d')) THEN
      !    CALL estimate_TKE(un, vn, rotn, z3d, ffmask)
      !    CALL put_fields('est_TKE', z3d)
      END IF

   END SUBROUTINE KEB_statistics

   SUBROUTINE put_fields(cdname, pfield3d)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE put_fields  ***   
      !! ** Purpose :   save 3d, surface, depth-averaged and xyz-averaged
      !!                fields
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)                , INTENT(in) :: cdname
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: pfield3d
      
      REAL(wp), DIMENSION(jpi,jpj) :: z2d
      REAL(wp) :: z0d
      CHARACTER(LEN=20) :: cdname_s, cdname_z, cdname_0d

      cdname_s  = cdname//'_s'
      cdname_z  = cdname//'_z'
      cdname_0d = cdname//'_0d'

      CALL iom_put(cdname, pfield3d)
      CALL iom_put(trim(cdname_s), pfield3d(:,:,1))
      IF (iom_use(trim(cdname_0d))) THEN
         z0d = average_xyz(pfield3d, e1t, e2t, e3t_0, tmask)
         CALL iom_put(trim(cdname_0d), z0d)
      END IF
      IF (iom_use(trim(cdname_z))) THEN
         CALL average_z(pfield3d, z2d, e3t_0, tmask)
         CALL iom_put(trim(cdname_z), z2d)
      END IF

   END SUBROUTINE put_fields
   
END MODULE KEB_module