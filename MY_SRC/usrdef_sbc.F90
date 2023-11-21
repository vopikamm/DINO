MODULE usrdef_sbc
   !!======================================================================
   !!                       ***  MODULE  usrdef_sbc  ***
   !! 
   !!                      ===  BASIN configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 2017-11  (J.Chanut)  user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_sbc    : user defined surface bounday conditions
   !!----------------------------------------------------------------------
   USE oce, ONLY : ts                                             ! ocean dynamics and tracers
   USE dom_oce, ONLY:                                             ! ocean space and time domain
   USE sbc_oce, ONLY: utau, vtau, taum, wndm, emp, sfx, qns, qsr  ! Surface boundary condition: ocean fields
   USE phycst                                                     ! physical constants
   USE sbcdcy, ONLY: sbc_dcy, nday_qsr                            ! surface boundary condition: diurnal cycle
   !
   USE usrdef_nam, ONLY : nn_forcingtype, rn_ztau0, rn_emp_prop, ln_ann_cyc, ln_diu_cyc,  &
               &          rn_trp, rn_srp, ln_qsr, rn_phi_max, rn_phi_min, rn_sstar_s,     &
               &          rn_sstar_n, rn_sstar_eq, rn_tstar_s, rn_tstar_n, rn_tstar_eq,   &
               &          ln_emp_field
   !
   USE in_out_manager  ! I/O manager
   USE fldread         ! read input fields
   USE iom             ! 
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined) 

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce      ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau  ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx  ! routine called by icestp.F90 for ice thermo

   INTEGER , PARAMETER ::   jp_emp  = 1            ! index of evaporation-precipation file
   INTEGER , PARAMETER ::   jpfld   = 1            ! maximum number of files to read
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf    ! structure of input fields (file informations, fields read)
   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztstar, zqsr_dayMean, zsstar   !: ztstar used in the heat forcing, zqsr_dayMean is the dayly averaged solar heat flux, zsstar for sfx

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_sbc.F90 10074 2018-08-28 16:15:49Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt, Kbb )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usrdef_sbc_oce  ***
      !!              
      !! ** Purpose :   provide at each time-step the GYRE surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   analytical configuration.
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt      ! ocean time step
      INTEGER, INTENT(in) ::   Kbb  ! ocean time index
      !
      INTEGER  ::   ji, jj
      REAL(wp) ::   z1_2L 
      REAL(wp) ::   zts_eq, ztstar_n, ztstar_s, zdeltaT
      REAL(wp) ::   zdeltaemp, zemp_mean, zF0, zconv, zaemp, zb, zdeltaemp2, zaemp2
      REAL(wp), DIMENSION(:), ALLOCATABLE  ::  znds_wnd_phi, znds_emp_phi, znds_tmp_phi, znds_slt_phi    ! Latitude of nodes    [degrees]
      REAL(wp), DIMENSION(:), ALLOCATABLE  ::  znds_wnd_val, znds_emp_val, znds_tmp_val, znds_slt_val    ! Values of nodes      [Pa]
      REAL(wp) ::   zdeltatau, zatau
      REAL(wp) ::   zf1, zf2, zf3, zweight2, zweight3
      REAL(wp) ::   za1, za2, za3
      REAL(wp) ::   zphi01, zphi02, zphi03
      REAL(wp) ::   z1_d1, zd2, z1_d2, zd3, z1_d3
      REAL(wp) ::   z1_s
      REAL(wp), DIMENSION(jpi,jpj) ::   zSurf
      REAL(wp) ::   zcos_sais1, zcos_sais2
      !!---------------------------------------------------------------------
      !
      SELECT CASE( nn_forcingtype )
      CASE(0)
         CALL ctl_stop( 'usrdef_sbc_oce : option not available anymore' )
      CASE(1)
         CALL ctl_stop( 'usrdef_sbc_oce : option not available anymore' )
      CASE(2)
         ! just a zonal wind stress, going eastward
         IF( kt .EQ. nit000 .AND. lwp ) THEN
            WRITE(numout,*)'usrdef_sbc_oce : analytical surface fluxes'
            WRITE(numout,*) 'Using just a zonal wind stress, goind eastward'
         ENDIF
         utau(:,:) = -0.1_wp
         vtau(:,:) = 0._wp
         taum(:,:) = 0._wp
         wndm(:,:) = 0._wp
         !
         emp (:,:) = 0._wp
         sfx (:,:) = 0._wp
         qns (:,:) = 0._wp
         qsr (:,:) = 0._wp
      CASE(3)
         ! no forcing
         IF( kt .EQ. nit000 .AND. lwp ) THEN
            WRITE(numout,*)'usrdef_sbc_oce : analytical surface fluxes'
            WRITE(numout,*) 'no forcing'
         ENDIF
         utau(:,:) = 0._wp
         vtau(:,:) = 0._wp
         taum(:,:) = 0._wp
         wndm(:,:) = 0._wp
         !
         emp (:,:) = 0._wp
         sfx (:,:) = 0._wp
         qns (:,:) = 0._wp
         qsr (:,:) = 0._wp
      CASE(4)
         ! Forcing close to Wolfe and Cessi 2014, JPO
         IF( kt == nit000 .AND. lwp ) THEN
            WRITE(numout,*) 'usrdef_sbc_oce : analytical surface fluxes'
            WRITE(numout,*) '~~~~~~~~~~~~~~'
            WRITE(numout,*) '   Forcing close to Wolfe and Cessi 2014, JPO'
            WRITE(numout,*) '      Zonal annual wind'
            WRITE(numout,*) '      Zonal annual E-P'
            IF( ln_qsr )   THEN
               WRITE(numout,*) '      Solar heat flux'
               IF( ln_ann_cyc )   THEN
                  WRITE(numout,*) '         Annual cycle of solar heat flux'
               ELSE
                  WRITE(numout,*) '         Annual average of solar heat flux (December = June)'
               ENDIF
               IF( ln_diu_cyc )   THEN
                  WRITE(numout,*) '         Diurnal cycle of solar heat flux'
               ELSE
                  WRITE(numout,*) '         Daily average of solar heat flux (day = night)'
               ENDIF
            ELSE
               WRITE(numout,*) '      No solar heat flux'
            ENDIF
            WRITE(numout,*) '      Zonal annual T* (heat forcing proportional to (SST - T*)'
            IF( rn_srp /= 0._wp )   THEN
               WRITE(numout,*) '      Zonal annual S* (salt flux proportional to (SSS - S*)'
            ELSE
               WRITE(numout,*) '      No salt restoring'
            ENDIF
         ENDIF
         ! Initialization of parameters
         !
         ! Computation of the day of the year (from Gyre)
         CALL compute_day_of_year( kt, zcos_sais1, zcos_sais2, ln_ann_cyc)
         !
         ! Zonal wind profile as is Marques et al. (2022)
         znds_wnd_phi    = (/rn_phi_min, -45._wp, -15._wp, 0._wp, 15._wp, 45._wp, rn_phi_max /)
         znds_wnd_val    = (/0._wp, 0.2_wp, -0.1_wp, -0.02_wp, -0.1_wp, 0.1_wp, 0._wp /)
         !
         ! Temperature restoring profile
         ! Chosen with the southern boundary always colder than the northern boundary
         ! Seasonnal cycle on T* coming from zcos_sais2 with extrema in July/January
         ! zts_eq         =  27._wp       ! [deg C] Temperature at the equator
         ! zdts_n         =  8._wp        ! [deg C] seasonal temperature difference in the north
         ! zdts_s         =   1._wp       ! [deg C] seasonal temperature difference in the south
         ! znds_tmp_phi    = (/rn_phi_min, 10. * zcos_sais2, rn_phi_max /)
         ! znds_tmp_val    = (/-0.5 * zdts_s * (1 + zcos_sais2) , zts_eq, 0.5 * zdts_n * (1 + zcos_sais2)/)
         !
         ! Evaporation - Precipitation
         ! znds_emp_phi    = (/rn_phi_min, -50._wp, -20._wp, 0._wp, 20._wp, 50._wp, rn_phi_max /)
         ! znds_emp_val    = (/-0.00001_wp, -0.00002_wp, 0.000035_wp, -0.000025_wp, 0.000035_wp, -0.00002_wp, -0.00001_wp /)
         !
         ! znds_slt_phi    = (/rn_phi_min, -40._wp, 0._wp, 40._wp, rn_phi_max /)
         ! znds_slt_val    = (/35.401_wp, 34.505_wp, 36.09_wp, 34.931_wp, 35.401_wp /)         
         !
         IF( kt == nit000 ) THEN
            ALLOCATE( ztstar(jpi,jpj) )   ! Allocation of ztstar
            ALLOCATE( zqsr_dayMean(jpi,jpj) )   ! Allocation of zqsr_dayMean
            IF( rn_srp /= 0._wp )   THEN
               ALLOCATE( zsstar(jpi,jpj) )   ! Allocation of zsstar
               DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                  ! Munday et al.
                  IF (gphit(ji, jj) <= 0) THEN
                     zsstar(ji,jj) = rn_sstar_s                                                                               &
                        & + (rn_sstar_eq - rn_sstar_s) * (1 + COS( 2 * rpi * gphit(ji,jj) / ( rn_phi_max - rn_phi_min) )) / 2 &
                        & - 1.25_wp * EXP( - gphit(ji,jj) ** 2 / 7.5_wp ** 2 )
                  ELSE
                     zsstar(ji,jj) = rn_sstar_n                                                                               &
                        & + (rn_sstar_eq - rn_sstar_n) * (1 + COS( 2 * rpi * gphit(ji,jj) / ( rn_phi_max - rn_phi_min) )) / 2 &
                        & - 1.25_wp * EXP( - gphit(ji,jj) ** 2 / 7.5_wp ** 2 )
                  ENDIF
                  ! Caneill
                  !zsstar(ji,jj) = 37.12_wp * EXP( - gphit(ji,jj)**2 / 260._wp**2 ) - 1.1_wp * EXP( - gphit(ji,jj)**2 / 7.5_wp**2 )
                  ! S-curve interpolation
                  !zsstar(ji,jj)  = znl_cbc(znds_slt_phi, znds_slt_val, gphit(ji,jj))
               END_2D
            ENDIF
         ENDIF
         !
         vtau(:,:) = 0._wp   ! no meridional wind stress
         wndm(:,:) = 0._wp   ! no use of 10 m wind speed
         !
         ! Time dependant forcing:
         !
         ! SALT FLUX
         IF( rn_srp /= 0._wp )   THEN
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               sfx(ji,jj) = ( rn_srp * ( ts(ji,jj,1,jp_sal, Kbb) - zsstar(ji,jj) ) ) * tmask(ji,jj,1)   ! Restoring salt flux
            END_2D
         ELSE
            sfx (:,:) = 0._wp   ! no salt flux
         ENDIF
         !
         ! Seasonnal cycle on T* coming from zcos_sais2
         !
         ztstar_s    = rn_tstar_s - 0.5_wp * zcos_sais2
         ztstar_n    = rn_tstar_n + 3.0_wp * zcos_sais2
         !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            ! Marques et al. (2022)
            utau(ji,jj)    = znl_cbc(znds_wnd_phi, znds_wnd_val, gphiu(ji,jj))
            taum(ji,jj)    = ABS( utau(ji,jj) )
            IF( utau(ji,jj) > 0 )   taum(ji,jj) = taum(ji,jj) * 1.3_wp   ! Boost in westerlies for TKE
            !
            ! EMP inspired from Wolfe and Cessi 2014, JPO and IPSL climate model output
            ! emp(ji,jj)     = rn_emp_prop * znl_cbc(znds_emp_phi, znds_emp_val, gphit(ji,jj))
            !
            ! T* inspired from Wolfe and Cessi 2014, JPO and IPSL climate model output
            ! ztstar(ji,jj)  = znl_cbc(znds_tmp_phi, znds_tmp_val, gphit(ji,jj))
            !
            ! T* inspired from Munday et al. (2012)
            IF ( gphit(ji, jj) <= 0 ) THEN
               ztstar(ji,jj) = ztstar_s                                                                                    &
                  & + (rn_tstar_eq - ztstar_s) * SIN( rpi * ( gphit(ji,jj) + rn_phi_max ) /  ( rn_phi_max - rn_phi_min) )
            ELSE
               ztstar(ji,jj) = ztstar_n                                                                                    &
                  & + (rn_tstar_eq - ztstar_n) * SIN( rpi * ( gphit(ji,jj) + rn_phi_max ) /  ( rn_phi_max - rn_phi_min) )
            ENDIF
         END_2D
         !
         IF( ln_emp_field ) THEN
            Call emp_flx( kt )
         ELSE
            emp(:,:) = 0._wp
         ENDIF
         !
         ! CALL remove_emp_mean()
         !
         ! Q SOLAR (from Gyre)
         ! see https://www.desmos.com/calculator/87duqiuxsf
         IF( ln_qsr )   THEN
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               zqsr_dayMean(ji,jj) = MAX(230._wp * COS( rpi * (gphit(ji,jj) - 23.5 * zcos_sais1 ) / ( 180._wp ) ) * tmask(ji,jj,1), 0._wp)
            END_2D
            CALL compute_diurn_cycle( kt, zqsr_dayMean, ln_diu_cyc )   ! Adding diurnal cycle if needed
         ELSE
            zqsr_dayMean(:,:) = 0._wp
            qsr(:,:) = 0._wp
         ENDIF
         !
         ! QNS
         ! take (SST - T*) into account, heat content of emp, remove qsr
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            qns(ji,jj) = (  rn_trp * ( ts(ji,jj,1,jp_tem, Kbb) - ztstar(ji,jj) ) &
                 &      - emp(ji,jj) * ts(ji,jj,1,jp_tem, Kbb) * rcp           &
                 &      - zqsr_dayMean(ji,jj)                         ) * tmask(ji,jj,1)
         END_2D
      CASE(5)
         ! Forcing inspired from Wolfe and Cessi 2014, JPO
         IF( kt == nit000 .AND. lwp ) THEN
            WRITE(numout,*) 'usrdef_sbc_oce : analytical surface fluxes'
            WRITE(numout,*) '~~~~~~~~~~~~~~'
            WRITE(numout,*) '   Forcing inspired from *Wolfe and Cessi 2014, JPO*, and from the *GYRE configuration*'
            WRITE(numout,*) '      Zonal annual wind'
            WRITE(numout,*) '      Zonal annual E-P'
            IF( ln_qsr )   THEN
               WRITE(numout,*) '      Solar heat flux'
               IF( ln_ann_cyc )   THEN
                  WRITE(numout,*) '         Annual cycle of solar heat flux'
               ELSE
                  WRITE(numout,*) '         Annual average of solar heat flux (December = June)'
               ENDIF
               IF( ln_diu_cyc )   THEN
                  WRITE(numout,*) '         Diurnal cycle of solar heat flux'
               ELSE
                  WRITE(numout,*) '         Daily average of solar heat flux (day = night)'
               ENDIF
            ELSE
               WRITE(numout,*) '      No solar heat flux'
            ENDIF
            WRITE(numout,*) '      Zonal T* (heat forcing proportional to (SST - T*)'
            IF( ln_ann_cyc )   THEN
               WRITE(numout,*) '         Annual cycle for T*'
            ELSE
               WRITE(numout,*) '         Annual average for T* (December = June)'
            ENDIF
            IF( rn_srp /= 0._wp )   THEN
               WRITE(numout,*) '      Zonal annual S* (salt flux proportional to (SSS - S*)'
            ELSE
               WRITE(numout,*) '      No salt restoring'
            ENDIF
         ENDIF
         ! Initialization of parameters
         ! Computation of the day of the year (from Gyre)
         CALL compute_day_of_year( kt, zcos_sais1, zcos_sais2, ln_ann_cyc)
         ! 
         ! zL = 132                                                                                                                                                                                                                                 ! [degrees] Approximative meridional extend of the basin
         ! Wind stress
         !ztau0         =   0.1_wp    ! [Pa]
         znds_wnd_phi    = (/rn_phi_min, -45._wp, -15._wp, 0._wp, 15._wp, 45._wp, rn_phi_max /)
         znds_wnd_val    = (/0._wp, 0.2_wp, -0.1_wp, -0.02_wp, -0.1_wp, 0.1_wp, 0._wp /)
         ! zatau         =   0.8_wp    ! [no unit]
         ! zdeltatau     =   5.77_wp   ! [deg North]
         ! T star and qns
         !ztrp          = -40._wp     ! [W/m2/K] retroaction term on heat fluxes 
         zts_eq        =  25._wp     ! [deg C] Temperature at the equator
         !zdts_n         =   1._wp     ! [deg C] Temperature in the north
         !zdts_s         =   0._wp     ! [deg C] Temperature in the south
         z1_2L         =   1._wp / (2._wp * rn_phi_max)
         zdeltaT       = 2           ! [deg C] half difference of temperature during winter and summer in the north (magnitude of the cos) !!rc TODO set in namelist
         ! zdeltaT_s     = 2           ! [deg C] half difference of temperature during winter and summer in the north (magnitude of the cos) !!rc TODO set in namelist !!dk necessary?
         ! EMP
         zconv         =   1._wp / ( 86400._wp)   ! convertion factor: 1 mm/day => 1/(3600*24) mm/s
         !!rc TODO put a1, a2 and a3 in namelist
         za1 = -3.24_wp              ! [mm/day] Set the amplitude of EMP at the equator
         za2 = 4.15_wp               ! [mm/day] Set the amplitude of EMP at mid latitude
         za3 = -1.59_wp              ! [mm/day] Set the amplitude of EMP at the northern part !!dk TODO is there an asymmetry in the emp forcing as well?
         zphi01 = 0._wp              ! [deg North]
         zphi02 = 20._wp             ! [deg North]
         zphi03 = 50._wp             ! [deg North]
         z1_d1 = 1._wp / 8._wp       ! [1 / deg North]
         zd2 = 30._wp                ! [deg North]
         z1_d2 = 1._wp / zd2         ! [1 / deg North]
         zd3 = 40._wp                ! [deg North]
         z1_d3 = 1._wp / zd3         ! [1 / deg North]
         z1_s = 1._wp / 10._wp       ! streching of the tanh function (i.e. smoothness of the filter)
         za1 = za1 * zconv           ! [mm/s] after  conversion
         za2 = za2 * zconv           ! [mm/s] after  conversion
         za3 = za3 * zconv           ! [mm/s] after  conversion
         !
         IF( kt == nit000 ) THEN
            ALLOCATE( ztstar(jpi,jpj) )   ! Allocation of ztstar
            ALLOCATE( zqsr_dayMean(jpi,jpj) )   ! Allocation of zqsr_dayMean
            IF( rn_srp /= 0._wp )   THEN
               ALLOCATE( zsstar(jpi,jpj) )   ! Allocation of zsstar
               ! See https://www.desmos.com/calculator/qrapqqrbfa
               DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                  zsstar(ji,jj) = 37.12_wp * EXP( - gphit(ji,jj)**2 / 260._wp**2 ) - 1.1_wp * EXP( - gphit(ji,jj)**2 / 7.5_wp**2 )
               END_2D
            ENDIF
         ENDIF
         ! necessary to compute at each time step because seasonnal variation of ztstar and solar heat flux
         vtau(:,:) = 0._wp   ! no meridional wind stress
         wndm(:,:) = 0._wp   ! no use of 10 m wind speed
         !
         ! SALT FLUX
         IF( rn_srp /= 0._wp )   THEN
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               sfx(ji,jj) = ( rn_srp * ( ts(ji,jj,1,jp_sal, Kbb) - zsstar(ji,jj) ) ) * tmask(ji,jj,1)   ! Restoring salt flux
            END_2D
         ELSE
            sfx (:,:) = 0._wp   ! no salt flux
         ENDIF
         !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            utau(ji,jj)    = znl_cbc(znds_wnd_phi, znds_wnd_val, gphiu(ji,jj))
            !utau(ji,jj) = rn_ztau0 * ( -COS((3 * rpi * gphit(ji,jj))/(2 * rn_phi_max)) + zatau * EXP(-gphit(ji,jj)**2/zdeltatau**2) )
            taum(ji,jj) = ABS( utau(ji,jj) )
            IF( utau(ji,jj) > 0 )   taum(ji,jj) = taum(ji,jj) * 1.3_wp   ! Boost in westerlies for TKE
            !
            ! EMP
            ! See https://www.desmos.com/calculator/v0vpbcc81h
            ! weights
            zweight2 = 0.5 * ( TANH((gphit(ji,jj) - zphi02 + zd2 * 0.5) * z1_s) - TANH((gphit(ji,jj) - zphi02 - zd2 * 0.5) * z1_s) )
            zweight3 = 0.5 * ( TANH((gphit(ji,jj) - zphi03 + zd3 * 0.5) * z1_s) - TANH((gphit(ji,jj) - zphi03 - zd3 * 0.5) * z1_s) )
            ! each component
            zf1 = za1 * EXP(-(gphit(ji,jj) - zphi01)**2 * z1_d1**2               )
            zf2 = za2 * SIN( (gphit(ji,jj) - zphi02)    * z1_d2 * rpi + 0.5 * rpi) * zweight2
            zf3 = za3 * SIN( (gphit(ji,jj) - zphi03)    * z1_d3 * rpi + 0.5 * rpi) * zweight3
            ! total
            emp (ji,jj) = rn_emp_prop * ( zf1 + zf2 + zf3 )
            ! The mean is removed later on (no simple analytical function)
            !
            ! T*
            ! See https://www.desmos.com/calculator/zij8tgy5yr
            ! Seasonnal cycle on T* coming from zcos_sais2
         END_2D
         !
         ! emp(:,:) = rn_emp_prop * emp(:,:)   ! taking the proportionality factor into account
         CALL remove_emp_mean()
         !
         ! Q SOLAR (flux from GYRE)
         ! see https://www.desmos.com/calculator/87duqiuxsf
         IF( ln_qsr )   THEN
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
              zqsr_dayMean(ji,jj) = 230._wp * COS( rpi * (gphit(ji,jj) - 23.5 * zcos_sais1 ) / ( 0.9_wp * 180._wp ) ) * tmask(ji,jj,1)
            END_2D
            CALL compute_diurn_cycle( kt, zqsr_dayMean, ln_diu_cyc )   ! Adding diurnal cycle if needed
         ELSE
            zqsr_dayMean(:,:) = 0._wp
            qsr(:,:) = 0._wp
         ENDIF
         !
         ! QNS
         ! take (SST - T*) into account, heat content of emp, remove zqsr_dayMean
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            qns(ji,jj) = (  rn_trp * ( ts(ji,jj,1,jp_tem, Kbb) - ztstar(ji,jj) ) &
                 &      - emp(ji,jj) * ts(ji,jj,1,jp_tem, Kbb) * rcp           &
                 &      - zqsr_dayMean(ji,jj)                         ) * tmask(ji,jj,1)
         END_2D
      END SELECT
     ! We call lbc_lnk to take the boundaries into account (especially the equator symmetrical condition)
     CALL lbc_lnk( 'usrdef_sbc', taum, 'T',  1. )
     CALL lbc_lnk( 'usrdef_sbc', wndm, 'T',  1. )
     CALL lbc_lnk( 'usrdef_sbc', utau, 'U', -1. )
     CALL lbc_lnk( 'usrdef_sbc', vtau, 'V', -1. )
     !
     CALL lbc_lnk( 'usrdef_sbc', emp , 'T',  1. )
     CALL lbc_lnk( 'usrdef_sbc', sfx , 'T',  1. )
     CALL lbc_lnk( 'usrdef_sbc', qns , 'T',  1. )
     CALL lbc_lnk( 'usrdef_sbc', qsr , 'T',  1. )
   END SUBROUTINE usrdef_sbc_oce

   
   SUBROUTINE usrdef_sbc_ice_tau( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_tau

   
   SUBROUTINE usrdef_sbc_ice_flx( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_flx

    
   SUBROUTINE remove_emp_mean()   ! !!dk this is not used for now  
     !!---------------------------------------------------------------------
     !!                    ***  ROUTINE remove_emp_mean  ***
     !!              
     !! ** Purpose :   Remove the average on emp (leading to a balanced flux)
     !!
     !! ** Method  :   emp = emp - MEAN(emp)
     !!                where MEAN is an average taking the surface of the cells into account
     !!
     !! ** Action  :   set emp
     !!---------------------------------------------------------------------
     REAL(wp), DIMENSION(jpi,jpj) ::   zSurf
     REAL(wp) ::   zemp_mean
     !!---------------------------------------------------------------------
     !
     ! Setting global E-P to 0
        ! 0: closed
        ! 1: cyclic east-west
        ! 2:cyclic north-south
        ! 7: cyclic east-west and north-south
        !
        ! Computing the surface of the cells
      zSurf( :     , :     ) = 0._wp                                         ! masking the halo + boundary points
      zSurf(2:jpi - 1,2:jpj - 1) = e1t(2:jpi - 1,2:jpj - 1) * e2t(2:jpi - 1,2:jpj - 1)   ! surface of the cells
      zemp_mean = glob_sum( 'usrdef_sbc', emp(:,:) * zSurf(:,:) ) / glob_sum( 'usrdef_sbc', zSurf(:,:) )
     !
      emp(:,:) = emp(:,:) - zemp_mean                     ! freshwater flux (=0 in domain average)
   END SUBROUTINE remove_emp_mean

   
   SUBROUTINE compute_day_of_year( kt, pcos_sais1, pcos_sais2, ll_ann_cyc )
     !!---------------------------------------------------------------------
     !!                    ***  SUBROUTINE compute_day_of_year  ***
     !!              
     !! ** Purpose :   Computation of the day of the year (from Gyre) as a cosine
     !!                pcos_sais1 is for heat flux, and is min the 21th of December, max the 21th of June
     !!                pcos_sais2 is for T*       , and is min the 21th of January , max the 21th of July
     !!
     !! ** Method  :   
     !!
     !! ** Action  : 
     !!---------------------------------------------------------------------
     INTEGER , INTENT(in   ) ::   kt      ! ocean time step
     REAL(wp), INTENT(  out) ::   pcos_sais1, pcos_sais2   ! cosine of the day of year (1 is for solar heat flux, 2 is for T* cycle)
     LOGICAL , INTENT(in   ), optional ::   ll_ann_cyc    ! if .false., the cos are set to zero. Default behaviour is true
     !
     LOGICAL  ::   ld_compute   ! local variable to know if we need to compute the cosines or set them to 0
     ! Variables to get the day of the year
     INTEGER  ::   zyear0                 ! initial year 
     INTEGER  ::   zmonth0                ! initial month
     INTEGER  ::   zday0                  ! initial day
     INTEGER  ::   zday_year0             ! initial day since january 1st
     REAL(wp) ::   ztime                  ! time in hour
     REAL(wp) ::   ztimemax1, ztimemin1   ! 21th June, and 21th decem. if date0 = 1st january
     REAL(wp) ::   ztimemax2, ztimemin2   ! 21th June, and 21th decem. if date0 = 1st january
     REAL(wp) ::   zyydd                 ! number of days in one year
     !!---------------------------------------------------------------------
     !
      IF( PRESENT(ll_ann_cyc) )   THEN
        ld_compute = ll_ann_cyc
      ELSE
        ld_compute = .TRUE.
      ENDIF
     !
      IF( ld_compute )   THEN
        zyydd = REAL(nyear_len(1),wp)
        zyear0     =   ndate0 / 10000._wp                                ! initial year
        zmonth0    = ( ndate0 - zyear0 * 10000._wp ) / 100._wp           ! initial month
        zday0      =   ndate0 - zyear0 * 10000._wp - zmonth0 * 100._wp   ! initial day betwen 1 and 30
        zday_year0 = ( zmonth0 - 1._wp ) * 30._wp + zday0                ! initial day betwen 1 and 360
        !
        ! current day (in hours) since january the 1st of the current year
        ztime = REAL( kt ) * rn_dt / (rmmss * rhhmm)   &       !  total incrementation (in hours)
             &      - (nyear  - 1) * rjjhh * zyydd              !  minus years since beginning of experiment (in hours)

        ztimemax1 = ((5.*30.)+21.)* 24.                      ! 21th june     at 24h in hours
        ztimemin1 = ztimemax1 + rjjhh * zyydd / 2            ! 21th december        in hours
        ztimemax2 = ((6.*30.)+21.)* 24.                      ! 21th july     at 24h in hours
        ztimemin2 = ztimemax2 - rjjhh * zyydd / 2            ! 21th january         in hours
        !                                                    ! NB: rjjhh * zyydd / 4 = one seasonal cycle in hours
        !
        ! 1/2 period between 21th June and 21th December and between 21th July and 21th January (1 for solar heat flux, 2 for T*)
        pcos_sais1 = COS( (ztime - ztimemax1) / (ztimemin1 - ztimemax1) * rpi )
        pcos_sais2 = COS( (ztime - ztimemax2) / (ztimemax2 - ztimemin2) * rpi )
      ELSE
        pcos_sais1 = 0._wp
        pcos_sais2 = 0._wp
      ENDIF
   END SUBROUTINE compute_day_of_year


   SUBROUTINE compute_diurn_cycle( kt, pqsr_dayMean, ll_diu_cyc )
     !!---------------------------------------------------------------------
     !!                    ***  SUBROUTINE compute_diurn_cycle  ***
     !!              
     !! ** Purpose :   Set the value of qsr.
     !!                If ll_diu_cyc is .true. or is not present, use the diurnal cycle.
     !!                If ll_diu_cyc is .false. use the daily mean.
     !!
     !! ** Method  :   
     !!
     !! ** Action  : 
     !!---------------------------------------------------------------------
     INTEGER ,                     INTENT(in   ) ::   kt      ! ocean time step
     REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pqsr_dayMean
     LOGICAL , optional          , INTENT(in   ) ::   ll_diu_cyc    ! if .false. use the mean value, if .true. use the diurnal cycle.
     !
     INTEGER  ::   ji, jj                 ! dummy loop indices
     REAL(wp), DIMENSION(jpi,jpj)                ::   zdiu_cyc ! diurnal cycle
     LOGICAL  ::   ld_compute   ! local variable to know if we need to compute the diurnal cycle
     !!---------------------------------------------------------------------
     !
     ! Adding diurnal cycle if needed
     IF( PRESENT(ll_diu_cyc) )   THEN
        ld_compute = ll_diu_cyc
     ELSE
        ld_compute = .TRUE.
     ENDIF
     IF( ld_compute )   THEN
         IF(  kt == nit000 )   nday_qsr = -1
         zdiu_cyc = sbc_dcy( pqsr_dayMean(:,:)) !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            qsr(ji,jj) = zdiu_cyc(ji,jj) * tmask(ji,jj,1)
         END_2D
     ELSE
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         qsr(ji,jj) =  pqsr_dayMean(ji,jj)
      END_2D
     ENDIF
   END SUBROUTINE compute_diurn_cycle

   FUNCTION znl_cbc( pnodes_phi, pnodes_val, pPhi ) RESULT( pprofile )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE  ***
      !!
      !! ** Purpose :   Fit a cubic zonal profile to a given set of lat/value pairs at the profile nodes.
      !!
      !! ** Method  :   see TODO'
      !!
      !!     
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in    ) ::   pnodes_phi   ! Latitude of nodes     [degrees]
      REAL(wp), DIMENSION(:), INTENT(in    ) ::   pnodes_val   ! Values of nodes       [unit]
      REAL(wp),               INTENT(in    ) ::   pPhi         ! Latitude at u-point   [degrees]
      !
      INTEGER                                ::   jn           ! dummy index
      INTEGER                                ::   ks, kn, kmin ! southern/northern node
      REAL(wp), DIMENSION(:), ALLOCATABLE    ::   zdphi        ! difference to node
      REAL(wp)                               ::   zs           ! saw-function
      REAL(wp)                               ::   zscurve      ! cubic s-curve between nodes
      REAL(wp)                               ::   pprofile     ! Zonal wind stress  [Pa]
      !
      ALLOCATE(zdphi(SIZE(pnodes_phi)))
      DO jn=1, SIZE(pnodes_phi)
         zdphi(jn) = pnodes_phi(jn) - pPhi
      END DO
      !
      kmin = MINLOC( ABS( zdphi ), DIM=1 )       ! nearest node to latitude
      !
      IF( zdphi(kmin)<=0 ) THEN
         ks = kmin
         kn = kmin+1
      ELSE
         kn = kmin
         ks = kmin-1
      ENDIF
      !
      zs = MIN( 1._wp, MAX( 0._wp, ( pPhi - pnodes_phi(ks) ) / ( pnodes_phi(kn) - pnodes_phi(ks) ) ) )
      pprofile = pnodes_val(ks) + ( pnodes_val(kn) - pnodes_val(ks) ) * ( 3 - 2 * zs ) * zs ** 2
      !!----------------------------------------------------------------------

   END FUNCTION znl_cbc

   SUBROUTINE emp_flx( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_flx  ***
      !!
      !! ** Purpose :   provide at each time step the Evaporation - Precipitation
      !!
      !! ** Method  : - READ net upward freshwater (evapo - precip) emp   (kg/m2/s)
      !!                   salt flux                              sfx   (pss*dh*rho/dt => g/m2/s)
      !!
      !!      CAUTION :  - never mask the surface stress fields
      !!
      !! ** Action  :   update at each time-step
      !!              - emp         upward mass flux (evap. - precip.)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj, jf            ! dummy indices
      INTEGER  ::   ierror                ! return error code
      INTEGER  ::   ios                   ! Local integer output status for namelist read
      !
      CHARACTER(len=100) ::  cn_dir                ! Root directory for location of flx files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i     ! array of namelist information structures
      TYPE(FLD_N) ::                     sn_emp    ! informations about the fields to be read
      NAMELIST/namsbc_flx/ cn_dir, sn_emp
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN                ! First call kt=nit000
         ! set file information
!         READ  ( numnam_ref, namsbc_flx, IOSTAT = ios, ERR = 901)
!901      IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_flx in reference namelist' )

         READ  ( numnam_cfg, namsbc_flx, IOSTAT = ios, ERR = 902 )
902      IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_flx in configuration namelist' )
         IF(lwm) WRITE ( numond, namsbc_flx )
         !                                         
         slf_i(jp_emp ) = sn_emp                   ! store namelist information in an array
         !
         ALLOCATE( sf(jpfld), STAT=ierror )        ! set sf structure
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'sbc_flx: unable to allocate sf structure' )   ;   RETURN
         ENDIF
         DO ji= 1, jpfld
            ALLOCATE( sf(ji)%fnow(jpi,jpj,1) )
            IF( slf_i(ji)%ln_tint ) ALLOCATE( sf(ji)%fdta(jpi,jpj,1,2) )
         END DO
         !                                         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir, 'emp_flx', 'flux formulation for ocean surface boundary condition', 'namsbc_flx' )
         !
      ENDIF

      CALL fld_read( kt, nn_fsbc, sf )                            ! input fields provided at the current time-step

      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN                        ! update ocean fluxes at each SBC frequency
      !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )                  ! set the ocean fluxes from read fields
            emp (ji,jj) =   sf(jp_emp )%fnow(ji,jj,1) * tmask(ji,jj,1)
         END_2D
         !
      ENDIF
      !
   END SUBROUTINE emp_flx

   SUBROUTINE test_compute_day_of_year()
		!!---------------------------------------------------------------------
     	!!                    ***  SUBROUTINE test_compute_day_of_year  ***
     	!!              
     	!! ** Purpose :   Testing SUBROUTINE 'compute_day_of_year' from 'usrdef_sbc.F90'.
     	!!
     	!!---------------------------------------------------------------------
		INTEGER 						::   kt       						! ocean time step
		REAL(wp) 					::   zcos_sais1, zcos_sais2   ! cosine of the day of year (1 is for solar heat flux, 2 is for T* cycle)
		REAL(wp), DIMENSION(13) ::   zcos1, zcos2   				! array of above
		INTEGER 					   ::   ksteps, inum 			   ! steps per day

		DO kt = 1, 6
			ksteps = INT((kt - 1) * 24 * 3600 / rn_dt)
			CALL compute_day_of_year( ksteps, zcos_sais1, zcos_sais2, .true.)
			zcos1(kt) = zcos_sais1
			zcos2(kt) = zcos_sais2
		END DO

		CALL iom_open( 'Cos_day_of_year', inum, ldwrt = .TRUE.)
      CALL iom_rstput( 0, 0, inum, 'solar_heat'	, zcos1 )
      CALL iom_rstput( 0, 0, inum, 't_star'		, zcos2 )	
      CALL iom_close( inum )

   END SUBROUTINE test_compute_day_of_year

   !!======================================================================
END MODULE usrdef_sbc
