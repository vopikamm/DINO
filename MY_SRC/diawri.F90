MODULE diawri
   !!======================================================================
   !!                     ***  MODULE  diawri  ***
   !! Ocean diagnostics :  write ocean output files
   !!=====================================================================
   !! History :  OPA  ! 1991-03  (M.-A. Foujols)  Original code
   !!            4.0  ! 1991-11  (G. Madec)
   !!                 ! 1992-06  (M. Imbard)  correction restart file
   !!                 ! 1992-07  (M. Imbard)  split into diawri and rstwri
   !!                 ! 1993-03  (M. Imbard)  suppress writibm
   !!                 ! 1998-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
   !!                 ! 1999-02  (E. Guilyardi)  name of netCDF files + variables
   !!            8.2  ! 2000-06  (M. Imbard)  Original code (diabort.F)
   !!   NEMO     1.0  ! 2002-06  (A.Bozec, E. Durand)  Original code (diainit.F)
   !!             -   ! 2002-09  (G. Madec)  F90: Free form and module
   !!             -   ! 2002-12  (G. Madec)  merge of diabort and diainit, F90
   !!                 ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2008-11  (B. Lemaire) creation from old diawri
   !!            3.7  ! 2014-01  (G. Madec) remove eddy induced velocity from no-IOM output
   !!                 !                     change name of output variables in dia_wri_state
   !!            4.0  ! 2020-10  (A. Nasser, S. Techene) add diagnostic for SWE
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_wri       : create the standart output files
   !!   dia_wri_state : create an output NetCDF file for a single instantaeous ocean state and forcing fields
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers 
   USE isf_oce
   USE isfcpl
   USE abl            ! abl variables in case ln_abl = .true.
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE dianam         ! build name of file (routine)
   USE diahth         ! thermocline diagnostics
   USE dynadv   , ONLY: ln_dynadv_vec
   USE icb_oce        ! Icebergs
   USE icbdia         ! Iceberg budgets
   USE ldftra         ! lateral physics: eddy diffusivity coef.
   USE ldfdyn         ! lateral physics: eddy viscosity   coef.
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE sbcssr         ! restoring term toward SST/SSS climatology
   USE sbcwave        ! wave parameters
   USE wet_dry        ! wetting and drying
   USE zdf_oce        ! ocean vertical physics
   USE zdfdrg         ! ocean vertical physics: top/bottom friction
   USE zdfmxl         ! mixed layer
   USE zdfosm         ! mixed layer
   !
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager ! I/O manager
   USE dia25h         ! 25h Mean output
   USE iom            ! 
   USE ioipsl         ! 

#if defined key_si3
   USE ice 
   USE icewri 
#endif
   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary
   USE diu_bulk        ! diurnal warm layer
   USE diu_coolskin    ! Cool skin

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_wri                 ! routines called by step.F90
   PUBLIC   dia_wri_state
   PUBLIC   dia_wri_alloc           ! Called by nemogcm module
#if ! defined key_xios   
   PUBLIC   dia_wri_alloc_abl       ! Called by sbcabl  module (if ln_abl = .true.)
#endif
   INTEGER ::   nid_T, nz_T, nh_T, ndim_T, ndim_hT   ! grid_T file
   INTEGER ::          nb_T              , ndim_bT   ! grid_T file
   INTEGER ::   nid_U, nz_U, nh_U, ndim_U, ndim_hU   ! grid_U file
   INTEGER ::   nid_V, nz_V, nh_V, ndim_V, ndim_hV   ! grid_V file
   INTEGER ::   nid_W, nz_W, nh_W                    ! grid_W file
   INTEGER ::   nid_A, nz_A, nh_A, ndim_A, ndim_hA   ! grid_ABL file   
   INTEGER ::   ndex(1)                              ! ???
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hT, ndex_hU, ndex_hV
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hA, ndex_A ! ABL
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_T, ndex_U, ndex_V
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_bT

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: diawri.F90 15141 2021-07-23 14:20:12Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_xios
   !!----------------------------------------------------------------------
   !!   'key_xios'                                        use IOM library
   !!----------------------------------------------------------------------
   INTEGER FUNCTION dia_wri_alloc()
      !
      dia_wri_alloc = 0
      !
   END FUNCTION dia_wri_alloc

   
   SUBROUTINE dia_wri( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri  ***
      !!                   
      !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
      !!      NETCDF format is used by default 
      !!
      !! ** Method  :  use iom_put
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kmm     ! ocean time level index
      !!
      INTEGER ::   ji, jj, jk       ! dummy loop indices
      INTEGER ::   ikbot            ! local integer
      REAL(wp)::   zztmp , zztmpx   ! local scalar
      REAL(wp)::   zztmp2, zztmpy   !   -      -
      REAL(wp)::   ze3
      REAL(wp), DIMENSION(A2D(     0))     ::   z2d   ! 2D workspace
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   z3d   ! 3D workspace
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   z3d_2 ! 2nd 3D workspace
      !!----------------------------------------------------------------------
      ! 
      IF( ln_timing )   CALL timing_start('dia_wri')
      ! 
      ! Output the initial state and forcings
      IF( ninist == 1 ) THEN                       
         CALL dia_wri_state( Kmm, 'output.init' )
         ninist = 0
      ENDIF

      ! initialize arrays
      z2d(:,:)   = 0._wp
      z3d(:,:,:) = 0._wp
      
      ! Output of initial vertical scale factor
      CALL iom_put("e3t_0", e3t_0(:,:,:) )
      CALL iom_put("e3u_0", e3u_0(:,:,:) )
      CALL iom_put("e3v_0", e3v_0(:,:,:) )
      CALL iom_put("e3f_0", e3f_0(:,:,:) )
      !
      IF ( iom_use("tpt_dep") ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) = gdept(ji,jj,jk,Kmm)
         END_3D
         CALL iom_put( "tpt_dep", z3d )
      ENDIF

      ! --- vertical scale factors --- !
      IF ( iom_use("e3t") .OR. iom_use("e3tdef") ) THEN  ! time-varying e3t
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) =  e3t(ji,jj,jk,Kmm)
         END_3D
         CALL iom_put( "e3t", z3d )
         IF ( iom_use("e3tdef") ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk )
               z3d(ji,jj,jk) = ( ( z3d(ji,jj,jk) - e3t_0(ji,jj,jk) ) / e3t_0(ji,jj,jk) * 100._wp * tmask(ji,jj,jk) ) ** 2
            END_3D
            CALL iom_put( "e3tdef", z3d ) 
         ENDIF
      ENDIF 
      IF ( iom_use("e3u") ) THEN                         ! time-varying e3u
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) =  e3u(ji,jj,jk,Kmm)
         END_3D 
         CALL iom_put( "e3u" , z3d )
      ENDIF
      IF ( iom_use("e3v") ) THEN                         ! time-varying e3v
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) =  e3v(ji,jj,jk,Kmm)
         END_3D
         CALL iom_put( "e3v" , z3d )
      ENDIF
      IF ( iom_use("e3w") ) THEN                         ! time-varying e3w
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) =  e3w(ji,jj,jk,Kmm)
         END_3D
         CALL iom_put( "e3w" , z3d )
      ENDIF
      IF ( iom_use("e3f") ) THEN                         ! time-varying e3f caution here at Kaa
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) =  e3f(ji,jj,jk)
         END_3D
         CALL iom_put( "e3f" , z3d )
      ENDIF

      IF ( iom_use("ssh") ) THEN
         IF( ll_wd ) THEN                                ! sea surface height (brought back to the reference used for wetting and drying)
            CALL iom_put( "ssh" , (ssh(:,:,Kmm)+ssh_ref)*ssmask(:,:) )
         ELSE
            CALL iom_put( "ssh" ,  ssh(:,:,Kmm) )        ! sea surface height
         ENDIF
      ENDIF

      IF( iom_use("wetdep") )    CALL iom_put( "wetdep" , ht_0(:,:) + ssh(:,:,Kmm) )   ! wet depth
         
#if defined key_qco
      IF( iom_use("ht") )   CALL iom_put( "ht" , ht(:,:)     )   ! water column at t-point
      IF( iom_use("hu") )   CALL iom_put( "hu" , hu(:,:,Kmm) )   ! water column at u-point
      IF( iom_use("hv") )   CALL iom_put( "hv" , hv(:,:,Kmm) )   ! water column at v-point
      IF( iom_use("hf") )   CALL iom_put( "hf" , hf_0(:,:)*( 1._wp + r3f(:,:) ) )   ! water column at f-point (caution here at Naa)
#endif

      ! --- tracers T&S --- !      
      CALL iom_put( "toce", ts(:,:,:,jp_tem,Kmm) )    ! 3D temperature
      CALL iom_put(  "sst", ts(:,:,1,jp_tem,Kmm) )    ! surface temperature

      IF ( iom_use("sbt") ) THEN
         DO_2D( 0, 0, 0, 0 )
            ikbot = mbkt(ji,jj)
            z2d(ji,jj) = ts(ji,jj,ikbot,jp_tem,Kmm)
         END_2D
         CALL iom_put( "sbt", z2d )                ! bottom temperature
      ENDIF
      
      CALL iom_put( "soce", ts(:,:,:,jp_sal,Kmm) )    ! 3D salinity
      CALL iom_put(  "sss", ts(:,:,1,jp_sal,Kmm) )    ! surface salinity
      IF ( iom_use("sbs") ) THEN
         DO_2D( 0, 0, 0, 0 )
            ikbot = mbkt(ji,jj)
            z2d(ji,jj) = ts(ji,jj,ikbot,jp_sal,Kmm)
         END_2D
         CALL iom_put( "sbs", z2d )                ! bottom salinity
      ENDIF

      IF( .NOT.lk_SWE )   CALL iom_put( "rhop", rhop(:,:,:) )          ! 3D potential density (sigma0)

      ! --- momentum --- !
      IF ( iom_use("taubot") ) THEN                ! bottom stress
         zztmp = rho0 * 0.25_wp
         z2d(:,:) = 0._wp
         DO_2D( 0, 0, 0, 0 )
            zztmp2 = (  ( rCdU_bot(ji+1,jj)+rCdU_bot(ji  ,jj) ) * uu(ji  ,jj,mbku(ji  ,jj),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji  ,jj)+rCdU_bot(ji-1,jj) ) * uu(ji-1,jj,mbku(ji-1,jj),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj  ) ) * vv(ji,jj  ,mbkv(ji,jj  ),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji,jj  )+rCdU_bot(ji,jj-1) ) * vv(ji,jj-1,mbkv(ji,jj-1),Kmm)  )**2
            z2d(ji,jj) = zztmp * SQRT( zztmp2 ) * tmask(ji,jj,1) 
            !
         END_2D
         CALL iom_put( "taubot", z2d )           
      ENDIF
         
      CALL iom_put( "uoce", uu(:,:,:,Kmm) )            ! 3D i-current
      CALL iom_put(  "ssu", uu(:,:,1,Kmm) )            ! surface i-current
      IF ( iom_use("sbu") ) THEN
         DO_2D( 0, 0, 0, 0 )
            ikbot = mbku(ji,jj)
            z2d(ji,jj) = uu(ji,jj,ikbot,Kmm)
         END_2D
         CALL iom_put( "sbu", z2d )                ! bottom i-current
      ENDIF
      
      CALL iom_put( "voce", vv(:,:,:,Kmm) )            ! 3D j-current
      CALL iom_put(  "ssv", vv(:,:,1,Kmm) )            ! surface j-current
      IF ( iom_use("sbv") ) THEN
         DO_2D( 0, 0, 0, 0 )
            ikbot = mbkv(ji,jj)
            z2d(ji,jj) = vv(ji,jj,ikbot,Kmm)
         END_2D
         CALL iom_put( "sbv", z2d )                ! bottom j-current
      ENDIF

      !                                            ! vertical velocity
      IF( ln_zad_Aimp ) THEN
         IF( iom_use('woce') ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk )
               z3d(ji,jj,jk) = ww(ji,jj,jk) + wi(ji,jj,jk)
            END_3D
            CALL iom_put( "woce", z3d )   ! explicit plus implicit parts
         ENDIF
      ELSE
         CALL iom_put( "woce", ww )
      ENDIF

      IF( iom_use('w_masstr') .OR. iom_use('w_masstr2') ) THEN   ! vertical mass transport & its square value
         !                     ! Caution: in the VVL case, it only correponds to the baroclinic mass transport.
         IF( ln_zad_Aimp ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk )
               z3d(ji,jj,jk) = rho0 * e1e2t(ji,jj) * ( ww(ji,jj,jk) + wi(ji,jj,jk) )
            END_3D
         ELSE
            DO_3D( 0, 0, 0, 0, 1, jpk )
               z3d(ji,jj,jk) = rho0 * e1e2t(ji,jj) * ww(ji,jj,jk)
            END_3D
         ENDIF
         CALL iom_put( "w_masstr" , z3d )  
         IF( iom_use('w_masstr2') )   CALL iom_put( "w_masstr2", z3d * z3d )
      ENDIF

      CALL iom_put( "avt" , avt )                  ! T vert. eddy diff. coef.
      CALL iom_put( "avs" , avs )                  ! S vert. eddy diff. coef.
      CALL iom_put( "avm" , avm )                  ! T vert. eddy visc. coef.

      IF( iom_use('logavt') )   CALL iom_put( "logavt", LOG( MAX( 1.e-20_wp, avt(:,:,:) ) ) )
      IF( iom_use('logavs') )   CALL iom_put( "logavs", LOG( MAX( 1.e-20_wp, avs(:,:,:) ) ) )

      IF ( iom_use("sssgrad") .OR. iom_use("sssgrad2") ) THEN
         DO_2D( 0, 0, 0, 0 )                       ! sss gradient
            zztmp  = ts(ji,jj,1,jp_sal,Kmm)
            zztmpx = (ts(ji+1,jj,1,jp_sal,Kmm) - zztmp) * r1_e1u(ji,jj) + (zztmp - ts(ji-1,jj  ,1,jp_sal,Kmm)) * r1_e1u(ji-1,jj)
            zztmpy = (ts(ji,jj+1,1,jp_sal,Kmm) - zztmp) * r1_e2v(ji,jj) + (zztmp - ts(ji  ,jj-1,1,jp_sal,Kmm)) * r1_e2v(ji,jj-1)
            z2d(ji,jj) = 0.25_wp * ( zztmpx * zztmpx + zztmpy * zztmpy )   &
               &                 * umask(ji,jj,1) * umask(ji-1,jj,1) * vmask(ji,jj,1) * vmask(ji,jj-1,1)
         END_2D
         CALL iom_put( "sssgrad2",  z2d )          ! square of module of sss gradient
         IF ( iom_use("sssgrad") ) THEN
            DO_2D( 0, 0, 0, 0 )
               z2d(ji,jj) = SQRT( z2d(ji,jj) )
            END_2D
            CALL iom_put( "sssgrad",  z2d )        ! module of sss gradient
         ENDIF
      ENDIF
         
      IF ( iom_use("sstgrad") .OR. iom_use("sstgrad2") ) THEN
         DO_2D( 0, 0, 0, 0 )                       ! sst gradient
            zztmp  = ts(ji,jj,1,jp_tem,Kmm)
            zztmpx = ( ts(ji+1,jj,1,jp_tem,Kmm) - zztmp ) * r1_e1u(ji,jj) + ( zztmp - ts(ji-1,jj  ,1,jp_tem,Kmm) ) * r1_e1u(ji-1,jj)
            zztmpy = ( ts(ji,jj+1,1,jp_tem,Kmm) - zztmp ) * r1_e2v(ji,jj) + ( zztmp - ts(ji  ,jj-1,1,jp_tem,Kmm) ) * r1_e2v(ji,jj-1)
            z2d(ji,jj) = 0.25_wp * ( zztmpx * zztmpx + zztmpy * zztmpy )   &
               &                 * umask(ji,jj,1) * umask(ji-1,jj,1) * vmask(ji,jj,1) * vmask(ji,jj-1,1)
         END_2D
         CALL iom_put( "sstgrad2",  z2d )          ! square of module of sst gradient
         IF ( iom_use("sstgrad") ) THEN
            DO_2D( 0, 0, 0, 0 )
               z2d(ji,jj) = SQRT( z2d(ji,jj) )
            END_2D
            CALL iom_put( "sstgrad",  z2d )        ! module of sst gradient
         ENDIF
      ENDIF
         
      ! heat and salt contents
      IF( iom_use("heatc") ) THEN
         z2d(:,:)  = 0._wp 
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            z2d(ji,jj) = z2d(ji,jj) + e3t(ji,jj,jk,Kmm) * ts(ji,jj,jk,jp_tem,Kmm) * tmask(ji,jj,jk)
         END_3D
         CALL iom_put( "heatc", rho0_rcp * z2d )   ! vertically integrated heat content (J/m2)
      ENDIF

      IF( iom_use("saltc") ) THEN
         z2d(:,:)  = 0._wp 
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            z2d(ji,jj) = z2d(ji,jj) + e3t(ji,jj,jk,Kmm) * ts(ji,jj,jk,jp_sal,Kmm) * tmask(ji,jj,jk)
         END_3D
         CALL iom_put( "saltc", rho0 * z2d )       ! vertically integrated salt content (PSU*kg/m2)
      ENDIF
      !
      IF( iom_use("salt2c") ) THEN
         z2d(:,:)  = 0._wp 
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            z2d(ji,jj) = z2d(ji,jj) + e3t(ji,jj,jk,Kmm) * ts(ji,jj,jk,jp_sal,Kmm) * ts(ji,jj,jk,jp_sal,Kmm) * tmask(ji,jj,jk)
         END_3D
         CALL iom_put( "salt2c", rho0 * z2d )      ! vertically integrated square of salt content (PSU2*kg/m2)
      ENDIF
      !
      IF ( iom_use("ke") .OR. iom_use("ke_int") ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zztmpx = uu(ji-1,jj  ,jk,Kmm) + uu(ji,jj,jk,Kmm)
            zztmpy = vv(ji  ,jj-1,jk,Kmm) + vv(ji,jj,jk,Kmm)
            z3d(ji,jj,jk) = 0.25_wp * ( zztmpx*zztmpx + zztmpy*zztmpy )
         END_3D
         CALL iom_put( "ke", z3d )                 ! kinetic energy

         z2d(:,:)  = 0._wp 
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            z2d(ji,jj) = z2d(ji,jj) + e3t(ji,jj,jk,Kmm) * z3d(ji,jj,jk) * e1e2t(ji,jj) * tmask(ji,jj,jk)
         END_3D
         CALL iom_put( "ke_int", z2d )             ! vertically integrated kinetic energy
      ENDIF
      !
      IF ( iom_use("sKE") ) THEN                   ! surface kinetic energy at T point
         z2d(:,:) = 0._wp
         DO_2D( 0, 0, 0, 0 )
            z2d(ji,jj) = 0.25_wp * ( uu(ji  ,jj,1,Kmm) * uu(ji  ,jj,1,Kmm) * e1e2u(ji  ,jj) * e3u(ji  ,jj,1,Kmm)  &
               &                   + uu(ji-1,jj,1,Kmm) * uu(ji-1,jj,1,Kmm) * e1e2u(ji-1,jj) * e3u(ji-1,jj,1,Kmm)  &
               &                   + vv(ji,jj  ,1,Kmm) * vv(ji,jj  ,1,Kmm) * e1e2v(ji,jj  ) * e3v(ji,jj  ,1,Kmm)  & 
               &                   + vv(ji,jj-1,1,Kmm) * vv(ji,jj-1,1,Kmm) * e1e2v(ji,jj-1) * e3v(ji,jj-1,1,Kmm)  )  &
               &                 * r1_e1e2t(ji,jj) / e3t(ji,jj,1,Kmm) * ssmask(ji,jj)
         END_2D
         IF ( iom_use("sKE" ) )  CALL iom_put( "sKE" , z2d )   
      ENDIF
      !    
      IF ( iom_use("ssKEf") ) THEN                 ! surface kinetic energy at F point
         z2d(:,:) = 0._wp                          ! CAUTION : only valid in SWE, not with bathymetry
         DO_2D( 0, 0, 0, 0 )
            z2d(ji,jj) = 0.25_wp * ( uu(ji,jj  ,1,Kmm) * uu(ji,jj  ,1,Kmm) * e1e2u(ji,jj  ) * e3u(ji,jj  ,1,Kmm)  &
               &                   + uu(ji,jj+1,1,Kmm) * uu(ji,jj+1,1,Kmm) * e1e2u(ji,jj+1) * e3u(ji,jj+1,1,Kmm)  &
               &                   + vv(ji  ,jj,1,Kmm) * vv(ji,jj  ,1,Kmm) * e1e2v(ji  ,jj) * e3v(ji  ,jj,1,Kmm)  & 
               &                   + vv(ji+1,jj,1,Kmm) * vv(ji+1,jj,1,Kmm) * e1e2v(ji+1,jj) * e3v(ji+1,jj,1,Kmm)  )  &
               &                 * r1_e1e2f(ji,jj) / e3f(ji,jj,1) * ssfmask(ji,jj)
         END_2D
         CALL iom_put( "ssKEf", z2d )                     
      ENDIF
      !
      CALL iom_put( "hdiv", hdiv )                 ! Horizontal divergence
      !
      IF( iom_use("u_masstr") .OR. iom_use("u_masstr_vint") .OR. iom_use("u_heattr") .OR. iom_use("u_salttr") ) THEN
         
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) = rho0 * uu(ji,jj,jk,Kmm) * e2u(ji,jj) * e3u(ji,jj,jk,Kmm) * umask(ji,jj,jk)
         END_3D
         CALL iom_put( "u_masstr"     , z3d )      ! mass transport in i-direction
         
         IF( iom_use("u_masstr_vint") ) THEN
            z2d(:,:) = 0._wp 
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk)
            END_3D
            CALL iom_put( "u_masstr_vint", z2d )   ! mass transport in i-direction vertical sum
         ENDIF
         IF( iom_use("u_heattr") ) THEN
            z2d(:,:) = 0._wp 
            zztmp = 0.5_wp * rcp
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               z2d(ji,jj) = z2d(ji,jj) + zztmp * z3d(ji,jj,jk) * ( ts(ji,jj,jk,jp_tem,Kmm) + ts(ji+1,jj,jk,jp_tem,Kmm) )
            END_3D
            CALL iom_put( "u_heattr", z2d )        ! heat transport in i-direction
         ENDIF
         IF( iom_use("u_salttr") ) THEN
            z2d(:,:) = 0._wp 
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               z2d(ji,jj) = z2d(ji,jj) +   0.5 * z3d(ji,jj,jk) * ( ts(ji,jj,jk,jp_sal,Kmm) + ts(ji+1,jj,jk,jp_sal,Kmm) )
            END_3D
            CALL iom_put( "u_salttr", z2d )        ! heat transport in i-direction
         ENDIF
         
      ENDIF
      
      IF( iom_use("v_masstr") .OR. iom_use("v_heattr") .OR. iom_use("v_salttr") ) THEN
         
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) = rho0 * vv(ji,jj,jk,Kmm) * e1v(ji,jj) * e3v(ji,jj,jk,Kmm) * vmask(ji,jj,jk)
         END_3D
         CALL iom_put( "v_masstr", z3d )           ! mass transport in j-direction
         
         IF( iom_use("v_heattr") ) THEN
            z2d(:,:) = 0._wp
            zztmp = 0.5_wp * rcp
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               z2d(ji,jj) = z2d(ji,jj) + zztmp * z3d(ji,jj,jk) * ( ts(ji,jj,jk,jp_tem,Kmm) + ts(ji,jj+1,jk,jp_tem,Kmm) )
            END_3D
            CALL iom_put( "v_heattr", z2d )        !  heat transport in j-direction
         ENDIF
         IF( iom_use("v_salttr") ) THEN
            z2d(:,:) = 0._wp 
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               z2d(ji,jj) = z2d(ji,jj) +   0.5 * z3d(ji,jj,jk) * ( ts(ji,jj,jk,jp_sal,Kmm) + ts(ji,jj+1,jk,jp_sal,Kmm) )
            END_3D
            CALL iom_put( "v_salttr", z2d )        !  heat transport in j-direction
         ENDIF

      ENDIF

      IF( iom_use("v_voltr") .OR. iom_use("v_heattr_3D") .OR. iom_use("v_salttr_3D") ) THEN
         
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) = vv(ji,jj,jk,Kmm) * e1v(ji,jj) * e3v(ji,jj,jk,Kmm) * vmask(ji,jj,jk)
         END_3D
         CALL iom_put( "v_voltr", z3d )                 ! volume transport in j-direction
         
         IF( iom_use("v_heattr_3D") ) THEN
            zztmp = 0.5_wp * rho0_rcp
            z3d_2(:,:,:) = 0._wp
            DO_3D( 0, 0, 0, 0, 1, jpk )
               z3d_2(ji,jj,jk) = zztmp * z3d(ji,jj,jk) * ( ts(ji,jj,jk,jp_tem,Kmm) + ts(ji,jj+1,jk,jp_tem,Kmm) )
            END_3D
            CALL iom_put( "v_heattr_3D", z3d_2 )        !  heat transport in j-direction
         ENDIF

         IF( iom_use("v_salttr_3D") ) THEN
            z3d_2(:,:,:) = 0._wp 
            DO_3D( 0, 0, 0, 0, 1, jpk )
               z3d_2(ji,jj,jk) = 0.5_wp * z3d(ji,jj,jk) * ( ts(ji,jj,jk,jp_sal,Kmm) + ts(ji,jj+1,jk,jp_sal,Kmm) )
            END_3D
            CALL iom_put( "v_salttr_3D", z3d_2 )        !  salt transport in j-direction
         ENDIF

      ENDIF

      IF( iom_use("tosmint") ) THEN
         z2d(:,:) = 0._wp
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            z2d(ji,jj) = z2d(ji,jj) + rho0 * e3t(ji,jj,jk,Kmm) * ts(ji,jj,jk,jp_tem,Kmm)
         END_3D
         CALL iom_put( "tosmint", z2d )            ! Vertical integral of temperature
      ENDIF
      IF( iom_use("somint") ) THEN
         z2d(:,:) = 0._wp
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            z2d(ji,jj) = z2d(ji,jj) + rho0 * e3t(ji,jj,jk,Kmm) * ts(ji,jj,jk,jp_sal,Kmm)
         END_3D
         CALL iom_put( "somint", z2d )             ! Vertical integral of salinity
      ENDIF

      CALL iom_put( "bn2", rn2 )                   ! Brunt-Vaisala buoyancy frequency (N^2)
      
      IF (ln_dia25h)   CALL dia_25h( kt, Kmm )     ! 25h averaging
      
      ! Output of surface vorticity terms
      !
      CALL iom_put( "ssplavor", ff_f )             ! planetary vorticity ( f )
      !
      IF ( iom_use("ssrelvor")    .OR. iom_use("ssEns")    .OR.   &
         & iom_use("ssrelpotvor") .OR. iom_use("ssabspotvor") ) THEN
         !
         z2d(:,:) = 0._wp 
         DO_2D( 0, 0, 0, 0 )
            z2d(ji,jj) = (   e2v(ji+1,jj  ) * vv(ji+1,jj  ,1,Kmm) - e2v(ji,jj) * vv(ji,jj,1,Kmm)    &
            &              - e1u(ji  ,jj+1) * uu(ji  ,jj+1,1,Kmm) + e1u(ji,jj) * uu(ji,jj,1,Kmm)  ) * r1_e1e2f(ji,jj)
         END_2D
         CALL iom_put( "ssrelvor", z2d )           ! relative vorticity ( zeta ) 
         !
         IF ( iom_use("ssEns") .OR. iom_use("ssrelpotvor") .OR. iom_use("ssabspotvor") ) THEN
            DO_2D( 0, 0, 0, 0 )  
               ze3 = (  e3t(ji,jj+1,1,Kmm) * e1e2t(ji,jj+1) + e3t(ji+1,jj+1,1,Kmm) * e1e2t(ji+1,jj+1)    &
                  &    + e3t(ji,jj  ,1,Kmm) * e1e2t(ji,jj  ) + e3t(ji+1,jj  ,1,Kmm) * e1e2t(ji+1,jj  )  ) * r1_e1e2f(ji,jj)
               IF( ze3 /= 0._wp ) THEN   ;   ze3 = 4._wp / ze3
               ELSE                      ;   ze3 = 0._wp
               ENDIF
               z2d(ji,jj) = ze3 * z2d(ji,jj) 
            END_2D
            CALL iom_put( "ssrelpotvor", z2d )     ! relative potential vorticity (zeta/h)
            !
            IF ( iom_use("ssEns") .OR. iom_use("ssabspotvor") ) THEN
               DO_2D( 0, 0, 0, 0 )
                  ze3 = (  e3t(ji,jj+1,1,Kmm) * e1e2t(ji,jj+1) + e3t(ji+1,jj+1,1,Kmm) * e1e2t(ji+1,jj+1)    &
                     &    + e3t(ji,jj  ,1,Kmm) * e1e2t(ji,jj  ) + e3t(ji+1,jj  ,1,Kmm) * e1e2t(ji+1,jj  )  ) * r1_e1e2f(ji,jj)
                  IF( ze3 /= 0._wp ) THEN   ;   ze3 = 4._wp / ze3
                  ELSE                      ;   ze3 = 0._wp
                  ENDIF
                  z2d(ji,jj) = ze3 * ff_f(ji,jj) + z2d(ji,jj) 
               END_2D
               CALL iom_put( "ssabspotvor", z2d )  ! absolute potential vorticity ( q )
               !
               IF ( iom_use("ssEns") ) THEN
                  DO_2D( 0, 0, 0, 0 )  
                     z2d(ji,jj) = 0.5_wp * z2d(ji,jj) * z2d(ji,jj) 
                  END_2D
                  CALL iom_put( "ssEns", z2d )     ! potential enstrophy ( 1/2*q2 )
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF( ln_timing )   CALL timing_stop('dia_wri')
      !
   END SUBROUTINE dia_wri

#else
   !!----------------------------------------------------------------------
   !!   Default option                                  use IOIPSL  library
   !!----------------------------------------------------------------------

   INTEGER FUNCTION dia_wri_alloc()
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ierr
      !!----------------------------------------------------------------------
      IF( nn_write == -1 ) THEN
         dia_wri_alloc = 0
      ELSE    
         ierr = 0
         ALLOCATE( ndex_hT(jpi*jpj) , ndex_T(jpi*jpj*jpk) ,     &
            &      ndex_hU(jpi*jpj) , ndex_U(jpi*jpj*jpk) ,     &
            &      ndex_hV(jpi*jpj) , ndex_V(jpi*jpj*jpk) , STAT=ierr(1) )
         !
         dia_wri_alloc = MAXVAL(ierr)
         CALL mpp_sum( 'diawri', dia_wri_alloc )
         !
      ENDIF
      !
   END FUNCTION dia_wri_alloc
 
   INTEGER FUNCTION dia_wri_alloc_abl()
      !!----------------------------------------------------------------------
	  ALLOCATE(   ndex_hA(jpi*jpj), ndex_A (jpi*jpj*jpkam1), STAT=dia_wri_alloc_abl)
      CALL mpp_sum( 'diawri', dia_wri_alloc_abl )
      !
   END FUNCTION dia_wri_alloc_abl

   
   SUBROUTINE dia_wri( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri  ***
      !!                   
      !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
      !!      NETCDF format is used by default 
      !!
      !! ** Method  :   At the beginning of the first time step (nit000), 
      !!      define all the NETCDF files and fields
      !!      At each time step call histdef to compute the mean if ncessary
      !!      Each nn_write time step, output the instantaneous or mean fields
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kmm  ! ocean time level index
      !
      LOGICAL ::   ll_print = .FALSE.                        ! =T print and flush numout
      CHARACTER (len=40) ::   clhstnam, clop, clmx           ! local names
      INTEGER  ::   inum = 11                                ! temporary logical unit
      INTEGER  ::   ji, jj, jk                               ! dummy loop indices
      INTEGER  ::   ierr                                     ! error code return from allocation
      INTEGER  ::   iimi, iima, ipk, it, itmod, ijmi, ijma   ! local integers
      INTEGER  ::   ipka                                     ! ABL
      INTEGER  ::   jn, ierror                               ! local integers
      REAL(dp) ::   zsto, zout, zmax, zjulian                ! local scalars
      !
      REAL(wp), DIMENSION(jpi,jpj    ) :: z2d     ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: z3d     ! 3D workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zw3d_abl   ! ABL 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( ninist == 1 ) THEN     !==  Output the initial state and forcings  ==!
         CALL dia_wri_state( Kmm, 'output.init' )
         ninist = 0
      ENDIF
      !
      IF( nn_write == -1 )   RETURN   ! we will never do any output
      ! 
      IF( ln_timing )   CALL timing_start('dia_wri')
      !
      ! 0. Initialisation
      ! -----------------

      ll_print = .FALSE.                  ! local variable for debugging
      ll_print = ll_print .AND. lwp

      ! Define frequency of output and means
      clop = "x"         ! no use of the mask value (require less cpu time and otherwise the model crashes)
#if defined key_diainstant
      zsto = nn_write * rn_Dt
      clop = "inst("//TRIM(clop)//")"
#else
      zsto=rn_Dt
      clop = "ave("//TRIM(clop)//")"
#endif
      zout = nn_write * rn_Dt
      zmax = ( nitend - nit000 + 1 ) * rn_Dt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = Nis0   ;   iima = Nie0
      ijmi = Njs0   ;   ijma = Nje0
      ipk = jpk
      IF(ln_abl) ipka = jpkam1

      ! define time axis
      it = kt
      itmod = kt - nit000 + 1

      ! 1. Define NETCDF files and fields at beginning of first time step
      ! -----------------------------------------------------------------

      IF( kt == nit000 ) THEN

         ! Define the NETCDF files (one per grid)

         ! Compute julian date from starting date of the run
         CALL ymds2ju( nyear, nmonth, nday, rn_Dt, zjulian )
         zjulian = zjulian - adatrj   !   set calendar origin to the beginning of the experiment
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'Date 0 used :', nit000, ' YEAR ', nyear,   &
            &                    ' MONTH ', nmonth, ' DAY ', nday, 'Julian day : ', zjulian
         IF(lwp)WRITE(numout,*) ' indexes of zoom = ', iimi, iima, ijmi, ijma,   &
                                 ' limit storage in depth = ', ipk

         ! WRITE root name in date.file for use by postpro
         IF(lwp) THEN
            CALL dia_nam( clhstnam, nn_write,' ' )
            CALL ctl_opn( inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
            WRITE(inum,*) clhstnam
            CLOSE(inum)
         ENDIF

         ! Define the T grid FILE ( nid_T )

         CALL dia_nam( clhstnam, nn_write, 'grid_T' )
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam    ! filename
         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           &  ! Horizontal grid: glamt and gphit
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, rn_Dt, nh_T, nid_T, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_T, "deptht", "Vertical T levels",      &  ! Vertical grid: gdept
            &           "m", ipk, gdept_1d, nz_T, "down" )
         !                                                            ! Index of ocean points
         CALL wheneq( jpi*jpj*ipk, tmask, 1, 1., ndex_T , ndim_T  )      ! volume
         CALL wheneq( jpi*jpj    , tmask, 1, 1., ndex_hT, ndim_hT )      ! surface
         !
         IF( ln_icebergs ) THEN
            !
            !! allocation cant go in dia_wri_alloc because ln_icebergs is only set after 
            !! that routine is called from nemogcm, so do it here immediately before its needed
            ALLOCATE( ndex_bT(jpi*jpj*nclasses), STAT=ierror )
            CALL mpp_sum( 'diawri', ierror )
            IF( ierror /= 0 ) THEN
               CALL ctl_stop('dia_wri: failed to allocate iceberg diagnostic array')
               RETURN
            ENDIF
            !
            !! iceberg vertical coordinate is class number
            CALL histvert( nid_T, "class", "Iceberg class",      &  ! Vertical grid: class
               &           "number", nclasses, class_num, nb_T )
            !
            !! each class just needs the surface index pattern
            ndim_bT = 3
            DO jn = 1,nclasses
               ndex_bT((jn-1)*jpi*jpj+1:jn*jpi*jpj) = ndex_hT(1:jpi*jpj)
            ENDDO
            !
         ENDIF

         ! Define the U grid FILE ( nid_U )

         CALL dia_nam( clhstnam, nn_write, 'grid_U' )
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam    ! filename
         CALL histbeg( clhstnam, jpi, glamu, jpj, gphiu,           &  ! Horizontal grid: glamu and gphiu
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, rn_Dt, nh_U, nid_U, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_U, "depthu", "Vertical U levels",      &  ! Vertical grid: gdept
            &           "m", ipk, gdept_1d, nz_U, "down" )
         !                                                            ! Index of ocean points
         CALL wheneq( jpi*jpj*ipk, umask, 1, 1., ndex_U , ndim_U  )      ! volume
         CALL wheneq( jpi*jpj    , umask, 1, 1., ndex_hU, ndim_hU )      ! surface

         ! Define the V grid FILE ( nid_V )

         CALL dia_nam( clhstnam, nn_write, 'grid_V' )                   ! filename
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam
         CALL histbeg( clhstnam, jpi, glamv, jpj, gphiv,           &  ! Horizontal grid: glamv and gphiv
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, rn_Dt, nh_V, nid_V, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_V, "depthv", "Vertical V levels",      &  ! Vertical grid : gdept
            &          "m", ipk, gdept_1d, nz_V, "down" )
         !                                                            ! Index of ocean points
         CALL wheneq( jpi*jpj*ipk, vmask, 1, 1., ndex_V , ndim_V  )      ! volume
         CALL wheneq( jpi*jpj    , vmask, 1, 1., ndex_hV, ndim_hV )      ! surface

         ! Define the W grid FILE ( nid_W )

         CALL dia_nam( clhstnam, nn_write, 'grid_W' )                   ! filename
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam
         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           &  ! Horizontal grid: glamt and gphit
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
            &          nit000-1, zjulian, rn_Dt, nh_W, nid_W, domain_id=nidom, snc4chunks=snc4set )
         CALL histvert( nid_W, "depthw", "Vertical W levels",      &  ! Vertical grid: gdepw
            &          "m", ipk, gdepw_1d, nz_W, "down" )

         IF( ln_abl ) THEN 
         ! Define the ABL grid FILE ( nid_A )
            CALL dia_nam( clhstnam, nn_write, 'grid_ABL' )
            IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam    ! filename
            CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,           &  ! Horizontal grid: glamt and gphit
               &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,       &
               &          nit000-1, zjulian, rn_Dt, nh_A, nid_A, domain_id=nidom, snc4chunks=snc4set )
            CALL histvert( nid_A, "ght_abl", "Vertical T levels",      &  ! Vertical grid: gdept
               &           "m", ipka, ght_abl(2:jpka), nz_A, "up" )
            !                                                            ! Index of ocean points
			ALLOCATE( zw3d_abl(jpi,jpj,ipka) ) 
			zw3d_abl(:,:,:) = 1._wp 
			CALL wheneq( jpi*jpj*ipka, zw3d_abl, 1, 1., ndex_A , ndim_A  )      ! volume
            CALL wheneq( jpi*jpj     , zw3d_abl, 1, 1., ndex_hA, ndim_hA )      ! surface
			DEALLOCATE(zw3d_abl)
         ENDIF
         !

         ! Declare all the output fields as NETCDF variables

         !                                                                                      !!! nid_T : 3D
         CALL histdef( nid_T, "votemper", "Temperature"                        , "C"      ,   &  ! tn
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         CALL histdef( nid_T, "vosaline", "Salinity"                           , "PSU"    ,   &  ! sn
            &          jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         IF(  .NOT.ln_linssh  ) THEN
            CALL histdef( nid_T, "vovvle3t", "Level thickness"                    , "m"      ,&  ! e3t n
            &             jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
            CALL histdef( nid_T, "vovvldep", "T point depth"                      , "m"      ,&  ! e3t n
            &             jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
            CALL histdef( nid_T, "vovvldef", "Squared level deformation"          , "%^2"    ,&  ! e3t n
            &             jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout )
         ENDIF
         !                                                                                      !!! nid_T : 2D
         CALL histdef( nid_T, "sosstsst", "Sea Surface temperature"            , "C"      ,   &  ! sst
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosaline", "Sea Surface Salinity"               , "PSU"    ,   &  ! sss
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sossheig", "Sea Surface Height"                 , "m"      ,   &  ! ssh
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowaflup", "Net Upward Water Flux"              , "Kg/m2/s",   &  ! (emp-rnf)
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sorunoff", "River runoffs"                      , "Kg/m2/s",   &  ! runoffs
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sosfldow", "downward salt flux"                 , "PSU/m2/s",  &  ! sfx
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         IF(  ln_linssh  ) THEN
            CALL histdef( nid_T, "sosst_cd", "Concentration/Dilution term on temperature"     &  ! emp * ts(:,:,1,jp_tem,Kmm)
            &                                                                  , "KgC/m2/s",  &  ! sosst_cd
            &             jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "sosss_cd", "Concentration/Dilution term on salinity"        &  ! emp * ts(:,:,1,jp_sal,Kmm)
            &                                                                  , "KgPSU/m2/s",&  ! sosss_cd
            &             jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         ENDIF
         CALL histdef( nid_T, "sohefldo", "Net Downward Heat Flux"             , "W/m2"   ,   &  ! qns + qsr
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "soshfldo", "Shortwave Radiation"                , "W/m2"   ,   &  ! qsr
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         IF( ALLOCATED(hmld) ) THEN   ! zdf_mxl not called by SWE
            CALL histdef( nid_T, "somixhgt", "Turbocline Depth"                   , "m"      ,   &  ! hmld
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "somxl010", "Mixed Layer Depth 0.01"             , "m"      ,   &  ! hmlp
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         ENDIF
         CALL histdef( nid_T, "soicecov", "Ice fraction"                       , "[0,1]"  ,   &  ! fr_i
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sowindsp", "wind speed at 10m"                  , "m/s"    ,   &  ! wndm
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         !
         IF( ln_abl ) THEN
            CALL histdef( nid_A, "t_abl", "Potential Temperature"     , "K"        ,       &  ! t_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout )
            CALL histdef( nid_A, "q_abl", "Humidity"                  , "kg/kg"    ,       &  ! q_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout ) 
            CALL histdef( nid_A, "u_abl", "Atmospheric U-wind   "     , "m/s"        ,     &  ! u_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout )
            CALL histdef( nid_A, "v_abl", "Atmospheric V-wind   "     , "m/s"    ,         &  ! v_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout ) 
            CALL histdef( nid_A, "tke_abl", "Atmospheric TKE   "     , "m2/s2"    ,        &  ! tke_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout ) 
            CALL histdef( nid_A, "avm_abl", "Atmospheric turbulent viscosity", "m2/s"   ,  &  ! avm_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout ) 
            CALL histdef( nid_A, "avt_abl", "Atmospheric turbulent diffusivity", "m2/s2",  &  ! avt_abl
               &          jpi, jpj, nh_A, ipka, 1, ipka, nz_A, 32, clop, zsto, zout ) 
            CALL histdef( nid_A, "pblh", "Atmospheric boundary layer height "  , "m",      &  ! pblh
               &          jpi, jpj, nh_A,  1  , 1, 1   , -99 , 32, clop, zsto, zout )		 			   
#if defined key_si3
            CALL histdef( nid_A, "oce_frac", "Fraction of open ocean"  , " ",      &  ! ato_i
               &          jpi, jpj, nh_A,  1  , 1, 1   , -99 , 32, clop, zsto, zout )
#endif
            CALL histend( nid_A, snc4chunks=snc4set )
         ENDIF
         !
         IF( ln_icebergs ) THEN
            CALL histdef( nid_T, "calving"             , "calving mass input"                       , "kg/s"   , &
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "calving_heat"        , "calving heat flux"                        , "XXXX"   , &
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "berg_floating_melt"  , "Melt rate of icebergs + bits"             , "kg/m2/s", &
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "berg_stored_ice"     , "Accumulated ice mass by class"            , "kg"     , &
               &          jpi, jpj, nh_T, nclasses  , 1, nclasses  , nb_T , 32, clop, zsto, zout )
            IF( ln_bergdia ) THEN
               CALL histdef( nid_T, "berg_melt"           , "Melt rate of icebergs"                    , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_buoy_melt"      , "Buoyancy component of iceberg melt rate"  , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_eros_melt"      , "Erosion component of iceberg melt rate"   , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_conv_melt"      , "Convective component of iceberg melt rate", "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_virtual_area"   , "Virtual coverage by icebergs"             , "m2"     , &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "bits_src"           , "Mass source of bergy bits"                , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "bits_melt"          , "Melt rate of bergy bits"                  , "kg/m2/s", &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "bits_mass"          , "Bergy bit density field"                  , "kg/m2"  , &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_mass"           , "Iceberg density field"                    , "kg/m2"  , &
                  &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
               CALL histdef( nid_T, "berg_real_calving"   , "Calving into iceberg class"               , "kg/s"   , &
                  &          jpi, jpj, nh_T, nclasses  , 1, nclasses  , nb_T , 32, clop, zsto, zout )
            ENDIF
         ENDIF

         IF( ln_ssr ) THEN
            CALL histdef( nid_T, "sohefldp", "Surface Heat Flux: Damping"         , "W/m2"   ,   &  ! qrp
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "sowafldp", "Surface Water Flux: Damping"        , "Kg/m2/s",   &  ! erp
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
            CALL histdef( nid_T, "sosafldp", "Surface salt flux: damping"         , "Kg/m2/s",   &  ! erp * sn
               &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         ENDIF
       
         clmx ="l_max(only(x))"    ! max index on a period
!         CALL histdef( nid_T, "sobowlin", "Bowl Index"                         , "W-point",   &  ! bowl INDEX 
!            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clmx, zsto, zout )
#if defined key_diahth
         CALL histdef( nid_T, "sothedep", "Thermocline Depth"                  , "m"      ,   & ! hth
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "so20chgt", "Depth of 20C isotherm"              , "m"      ,   & ! hd20
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "so28chgt", "Depth of 28C isotherm"              , "m"      ,   & ! hd28
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
         CALL histdef( nid_T, "sohtc300", "Heat content 300 m"                 , "J/m2"   ,   & ! htc3
            &          jpi, jpj, nh_T, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
#endif

         CALL histend( nid_T, snc4chunks=snc4set )

         !                                                                                      !!! nid_U : 3D
         CALL histdef( nid_U, "vozocrtx", "Zonal Current"                      , "m/s"    ,   &  ! uu(:,:,:,Kmm)
            &          jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout )
         IF( ln_wave .AND. ln_sdw) THEN
            CALL histdef( nid_U, "sdzocrtx", "Stokes Drift Zonal Current"         , "m/s"    ,   &  ! usd
               &          jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout )
         ENDIF
         !                                                                                      !!! nid_U : 2D
         CALL histdef( nid_U, "sozotaux", "Wind Stress along i-axis"           , "N/m2"   ,   &  ! utau
            &          jpi, jpj, nh_U, 1  , 1, 1  , - 99, 32, clop, zsto, zout )

         CALL histend( nid_U, snc4chunks=snc4set )

         !                                                                                      !!! nid_V : 3D
         CALL histdef( nid_V, "vomecrty", "Meridional Current"                 , "m/s"    ,   &  ! vv(:,:,:,Kmm)
            &          jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout )
         IF( ln_wave .AND. ln_sdw) THEN
            CALL histdef( nid_V, "sdmecrty", "Stokes Drift Meridional Current"    , "m/s"    ,   &  ! vsd
               &          jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout )
         ENDIF
         !                                                                                      !!! nid_V : 2D
         CALL histdef( nid_V, "sometauy", "Wind Stress along j-axis"           , "N/m2"   ,   &  ! vtau
            &          jpi, jpj, nh_V, 1  , 1, 1  , - 99, 32, clop, zsto, zout )

         CALL histend( nid_V, snc4chunks=snc4set )

         !                                                                                      !!! nid_W : 3D
         CALL histdef( nid_W, "vovecrtz", "Vertical Velocity"                  , "m/s"    ,   &  ! ww
            &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         CALL histdef( nid_W, "votkeavt", "Vertical Eddy Diffusivity"          , "m2/s"   ,   &  ! avt
            &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         CALL histdef( nid_W, "votkeavm", "Vertical Eddy Viscosity"             , "m2/s"  ,   &  ! avm
            &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )

         IF( ln_zdfddm ) THEN
            CALL histdef( nid_W,"voddmavs","Salt Vertical Eddy Diffusivity"    , "m2/s"   ,   &  ! avs
               &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         ENDIF
         
         IF( ln_wave .AND. ln_sdw) THEN
            CALL histdef( nid_W, "sdvecrtz", "Stokes Drift Vertical Current"   , "m/s"    ,   &  ! wsd
               &          jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout )
         ENDIF
         !                                                                                      !!! nid_W : 2D
         CALL histend( nid_W, snc4chunks=snc4set )

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization'
         IF(ll_print) CALL FLUSH(numout )

      ENDIF

      ! 2. Start writing data
      ! ---------------------

      ! ndex(1) est utilise ssi l'avant dernier argument est different de 
      ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
      ! donne le nombre d'elements, et ndex la liste des indices a sortir

      IF( lwp .AND. MOD( itmod, nn_write ) == 0 ) THEN 
         WRITE(numout,*) 'dia_wri : write model outputs in NetCDF files at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~ '
      ENDIF

      IF( .NOT.ln_linssh ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) = ts(ji,jj,jk,jp_tem,Kmm) * e3t(ji,jj,jk,Kmm)
         END_3D
         CALL histwrite( nid_T, "votemper", it, z3d, ndim_T , ndex_T  )   ! heat content
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) = ts(ji,jj,jk,jp_sal,Kmm) * e3t(ji,jj,jk,Kmm)
         END_3D
         CALL histwrite( nid_T, "vosaline", it, z3d, ndim_T , ndex_T  )   ! salt content
         DO_2D( 0, 0, 0, 0 )
            z2d(ji,jj   ) = ts(ji,jj, 1,jp_tem,Kmm) * e3t(ji,jj, 1,Kmm)
         END_2D
         CALL histwrite( nid_T, "sosstsst", it, z2d, ndim_hT, ndex_hT )   ! sea surface heat content
         DO_2D( 0, 0, 0, 0 )
            z2d(ji,jj   ) = ts(ji,jj, 1,jp_sal,Kmm) * e3t(ji,jj, 1,Kmm)
         END_2D
         CALL histwrite( nid_T, "sosaline", it, z2d, ndim_hT, ndex_hT )   ! sea surface salinity content
      ELSE
         CALL histwrite( nid_T, "votemper", it, ts(:,:,:,jp_tem,Kmm) , ndim_T , ndex_T  )   ! temperature
         CALL histwrite( nid_T, "vosaline", it, ts(:,:,:,jp_sal,Kmm) , ndim_T , ndex_T  )   ! salinity
         CALL histwrite( nid_T, "sosstsst", it, ts(:,:,1,jp_tem,Kmm) , ndim_hT, ndex_hT )   ! sea surface temperature
         CALL histwrite( nid_T, "sosaline", it, ts(:,:,1,jp_sal,Kmm) , ndim_hT, ndex_hT )   ! sea surface salinity
      ENDIF
      IF( .NOT.ln_linssh ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
           z3d(ji,jj,jk) = e3t(ji,jj,jk,Kmm)     ! 3D workspace for qco substitution
         END_3D
         CALL histwrite( nid_T, "vovvle3t", it, z3d        , ndim_T , ndex_T  )   ! level thickness
         DO_3D( 0, 0, 0, 0, 1, jpk )
           z3d(ji,jj,jk) = gdept(ji,jj,jk,Kmm)   ! 3D workspace for qco substitution
         END_3D
         CALL histwrite( nid_T, "vovvldep", it, z3d        , ndim_T , ndex_T  )   ! t-point depth 
         DO_3D( 0, 0, 0, 0, 1, jpk )
            z3d(ji,jj,jk) = ( ( e3t(ji,jj,jk,Kmm) - e3t_0(ji,jj,jk) ) / e3t_0(ji,jj,jk) * 100._wp * tmask(ji,jj,jk) ) ** 2
         END_3D         
         CALL histwrite( nid_T, "vovvldef", it, z3d        , ndim_T , ndex_T  )   ! level thickness deformation
      ENDIF
      CALL histwrite( nid_T, "sossheig", it, ssh(:,:,Kmm)  , ndim_hT, ndex_hT )   ! sea surface height
      DO_2D( 0, 0, 0, 0 )
         z2d(ji,jj) = emp(ji,jj) - rnf(ji,jj)
      END_2D
      CALL histwrite( nid_T, "sowaflup", it, z2d           , ndim_hT, ndex_hT )   ! upward water flux
      CALL histwrite( nid_T, "sorunoff", it, rnf           , ndim_hT, ndex_hT )   ! river runoffs
      CALL histwrite( nid_T, "sosfldow", it, sfx           , ndim_hT, ndex_hT )   ! downward salt flux 
                                                                                  ! (includes virtual salt flux beneath ice 
                                                                                  ! in linear free surface case)
      IF( ln_linssh ) THEN
         DO_2D( 0, 0, 0, 0 )
            z2d(ji,jj) = emp (ji,jj) * ts(ji,jj,1,jp_tem,Kmm)
         END_2D
         CALL histwrite( nid_T, "sosst_cd", it, z2d, ndim_hT, ndex_hT )          ! c/d term on sst
         DO_2D( 0, 0, 0, 0 )
            z2d(ji,jj) = emp (ji,jj) * ts(ji,jj,1,jp_sal,Kmm)
         END_2D
         CALL histwrite( nid_T, "sosss_cd", it, z2d, ndim_hT, ndex_hT )          ! c/d term on sss
      ENDIF
      DO_2D( 0, 0, 0, 0 )
         z2d(ji,jj) = qsr(ji,jj) + qns(ji,jj)
      END_2D
      CALL histwrite( nid_T, "sohefldo", it, z2d           , ndim_hT, ndex_hT )   ! total heat flux
      CALL histwrite( nid_T, "soshfldo", it, qsr           , ndim_hT, ndex_hT )   ! solar heat flux
      IF( ALLOCATED(hmld) ) THEN   ! zdf_mxl not called by SWE
         CALL histwrite( nid_T, "somixhgt", it, hmld          , ndim_hT, ndex_hT )   ! turbocline depth
         CALL histwrite( nid_T, "somxl010", it, hmlp          , ndim_hT, ndex_hT )   ! mixed layer depth
      ENDIF
      CALL histwrite( nid_T, "soicecov", it, fr_i          , ndim_hT, ndex_hT )   ! ice fraction   
      CALL histwrite( nid_T, "sowindsp", it, wndm          , ndim_hT, ndex_hT )   ! wind speed   
      !
      IF( ln_abl ) THEN 
         ALLOCATE( zw3d_abl(jpi,jpj,jpka) )
         IF( ln_mskland )   THEN 
            DO jk=1,jpka
               zw3d_abl(:,:,jk) = tmask(:,:,1)
            END DO       
         ELSE
            zw3d_abl(:,:,:) = 1._wp		 
         ENDIF			
         CALL histwrite( nid_A,  "pblh"   , it, pblh(:,:)                  *zw3d_abl(:,:,1     ), ndim_hA, ndex_hA )   ! pblh 
         CALL histwrite( nid_A,  "u_abl"  , it, u_abl   (:,:,2:jpka,nt_n  )*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! u_abl
         CALL histwrite( nid_A,  "v_abl"  , it, v_abl   (:,:,2:jpka,nt_n  )*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! v_abl
         CALL histwrite( nid_A,  "t_abl"  , it, tq_abl  (:,:,2:jpka,nt_n,1)*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! t_abl
         CALL histwrite( nid_A,  "q_abl"  , it, tq_abl  (:,:,2:jpka,nt_n,2)*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! q_abl		 
         CALL histwrite( nid_A,  "tke_abl", it, tke_abl (:,:,2:jpka,nt_n  )*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! tke_abl
         CALL histwrite( nid_A,  "avm_abl", it, avm_abl (:,:,2:jpka       )*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! avm_abl
         CALL histwrite( nid_A,  "avt_abl", it, avt_abl (:,:,2:jpka       )*zw3d_abl(:,:,2:jpka), ndim_A , ndex_A  )   ! avt_abl	
#if defined key_si3
         CALL histwrite( nid_A,  "oce_frac"   , it, ato_i(:,:)                                  , ndim_hA, ndex_hA )   ! ato_i
#endif
         DEALLOCATE(zw3d_abl)
      ENDIF
      !
      IF( ln_icebergs ) THEN
         !
         CALL histwrite( nid_T, "calving"             , it, berg_grid%calving      , ndim_hT, ndex_hT )  
         CALL histwrite( nid_T, "calving_heat"        , it, berg_grid%calving_hflx , ndim_hT, ndex_hT )         
         CALL histwrite( nid_T, "berg_floating_melt"  , it, berg_grid%floating_melt, ndim_hT, ndex_hT )  
         !
         CALL histwrite( nid_T, "berg_stored_ice"     , it, berg_grid%stored_ice   , ndim_bT, ndex_bT )
         !
         IF( ln_bergdia ) THEN
            CALL histwrite( nid_T, "berg_melt"           , it, berg_melt        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "berg_buoy_melt"      , it, buoy_melt        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "berg_eros_melt"      , it, eros_melt        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "berg_conv_melt"      , it, conv_melt        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "berg_virtual_area"   , it, virtual_area     , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "bits_src"            , it, bits_src         , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "bits_melt"           , it, bits_melt        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "bits_mass"           , it, bits_mass        , ndim_hT, ndex_hT   )  
            CALL histwrite( nid_T, "berg_mass"           , it, berg_mass        , ndim_hT, ndex_hT   )  
            !
            CALL histwrite( nid_T, "berg_real_calving"   , it, real_calving     , ndim_bT, ndex_bT   )
         ENDIF
      ENDIF

      IF( ln_ssr ) THEN
         CALL histwrite( nid_T, "sohefldp", it, qrp           , ndim_hT, ndex_hT )   ! heat flux damping
         CALL histwrite( nid_T, "sowafldp", it, erp           , ndim_hT, ndex_hT )   ! freshwater flux damping
         DO_2D( 0, 0, 0, 0 )
            z2d(ji,jj) = erp(ji,jj) * ts(ji,jj,1,jp_sal,Kmm) * tmask(ji,jj,1)
         END_2D
         CALL histwrite( nid_T, "sosafldp", it, z2d           , ndim_hT, ndex_hT )   ! salt flux damping
      ENDIF
!      zw2d(:,:) = FLOAT( nmln(:,:) ) * tmask(:,:,1)
!      CALL histwrite( nid_T, "sobowlin", it, zw2d          , ndim_hT, ndex_hT )   ! ???

#if defined key_diahth
      CALL histwrite( nid_T, "sothedep", it, hth           , ndim_hT, ndex_hT )   ! depth of the thermocline
      CALL histwrite( nid_T, "so20chgt", it, hd20          , ndim_hT, ndex_hT )   ! depth of the 20 isotherm
      CALL histwrite( nid_T, "so28chgt", it, hd28          , ndim_hT, ndex_hT )   ! depth of the 28 isotherm
      CALL histwrite( nid_T, "sohtc300", it, htc3          , ndim_hT, ndex_hT )   ! first 300m heaat content
#endif

      CALL histwrite( nid_U, "vozocrtx", it, uu(:,:,:,Kmm) , ndim_U , ndex_U )    ! i-current
      CALL histwrite( nid_U, "sozotaux", it, utau          , ndim_hU, ndex_hU )   ! i-wind stress

      CALL histwrite( nid_V, "vomecrty", it, vv(:,:,:,Kmm) , ndim_V , ndex_V  )   ! j-current
      CALL histwrite( nid_V, "sometauy", it, vtau          , ndim_hV, ndex_hV )   ! j-wind stress

      IF( ln_zad_Aimp ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
           z3d(ji,jj,jk) = ww(ji,jj,jk) + wi(ji,jj,jk)
         END_3D         
         CALL histwrite( nid_W, "vovecrtz", it, z3d         , ndim_T, ndex_T )    ! vert. current
      ELSE
         CALL histwrite( nid_W, "vovecrtz", it, ww          , ndim_T, ndex_T )    ! vert. current
      ENDIF
      CALL histwrite( nid_W, "votkeavt", it, avt            , ndim_T, ndex_T )    ! T vert. eddy diff. coef.
      CALL histwrite( nid_W, "votkeavm", it, avm            , ndim_T, ndex_T )    ! T vert. eddy visc. coef.
      IF( ln_zdfddm ) THEN
         CALL histwrite( nid_W, "voddmavs", it, avs         , ndim_T, ndex_T )    ! S vert. eddy diff. coef.
      ENDIF

      IF( ln_wave .AND. ln_sdw ) THEN
         CALL histwrite( nid_U, "sdzocrtx", it, usd         , ndim_U , ndex_U )    ! i-StokesDrift-current
         CALL histwrite( nid_V, "sdmecrty", it, vsd         , ndim_V , ndex_V )    ! j-StokesDrift-current
         CALL histwrite( nid_W, "sdvecrtz", it, wsd         , ndim_T , ndex_T )    ! StokesDrift vert. current
      ENDIF

      ! 3. Close all files
      ! ---------------------------------------
      IF( kt == nitend ) THEN
         CALL histclo( nid_T )
         CALL histclo( nid_U )
         CALL histclo( nid_V )
         CALL histclo( nid_W )
         IF(ln_abl) CALL histclo( nid_A )
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('dia_wri')
      !
   END SUBROUTINE dia_wri
#endif

   SUBROUTINE dia_wri_state( Kmm, cdfile_name )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_state  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ocean state and forcing fields.
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! ** Method  :   NetCDF files using ioipsl
      !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
      !!      File 'output.abort.nc' is created in case of abnormal job end
      !!----------------------------------------------------------------------
      INTEGER           , INTENT( in ) ::   Kmm              ! time level index
      CHARACTER (len=* ), INTENT( in ) ::   cdfile_name      ! name of the file created
      !!
      INTEGER ::   ji, jj, jk       ! dummy loop indices
      INTEGER ::   inum
      REAL(wp), DIMENSION(jpi,jpj)     :: z2d      
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: z3d      
      !!----------------------------------------------------------------------
      ! 
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
         WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
         WRITE(numout,*) '                and named :', cdfile_name, '...nc'
      ENDIF 
      !
      CALL iom_open( TRIM(cdfile_name), inum, ldwrt = .TRUE. )
      !
      CALL iom_rstput( 0, 0, inum, 'votemper', ts(:,:,:,jp_tem,Kmm) )    ! now temperature
      CALL iom_rstput( 0, 0, inum, 'vosaline', ts(:,:,:,jp_sal,Kmm) )    ! now salinity
      CALL iom_rstput( 0, 0, inum, 'sossheig', ssh(:,:,Kmm)         )    ! sea surface height
      CALL iom_rstput( 0, 0, inum, 'vozocrtx', uu(:,:,:,Kmm)        )    ! now i-velocity
      CALL iom_rstput( 0, 0, inum, 'vomecrty', vv(:,:,:,Kmm)        )    ! now j-velocity
      IF( ln_zad_Aimp ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
           z3d(ji,jj,jk) = ww(ji,jj,jk) + wi(ji,jj,jk)
         END_3D
         CALL iom_rstput( 0, 0, inum, 'vovecrtz', z3d            )    ! now k-velocity
      ELSE
         CALL iom_rstput( 0, 0, inum, 'vovecrtz', ww             )    ! now k-velocity
      ENDIF
      CALL iom_rstput( 0, 0, inum, 'risfdep', risfdep            )
      CALL iom_rstput( 0, 0, inum, 'ht'     , ht(:,:)            )    ! now water column height
      !
      IF ( ln_isf ) THEN
         IF (ln_isfcav_mlt) THEN
            CALL iom_rstput( 0, 0, inum, 'fwfisf_cav', fwfisf_cav          )
            CALL iom_rstput( 0, 0, inum, 'rhisf_cav_tbl', rhisf_tbl_cav    )
            CALL iom_rstput( 0, 0, inum, 'rfrac_cav_tbl', rfrac_tbl_cav    )
            CALL iom_rstput( 0, 0, inum, 'misfkb_cav', REAL(misfkb_cav,wp) )
            CALL iom_rstput( 0, 0, inum, 'misfkt_cav', REAL(misfkt_cav,wp) )
            CALL iom_rstput( 0, 0, inum, 'mskisf_cav', REAL(mskisf_cav,wp), ktype = jp_i1 )
         END IF
         IF (ln_isfpar_mlt) THEN
            CALL iom_rstput( 0, 0, inum, 'isfmsk_par', REAL(mskisf_par,wp) )
            CALL iom_rstput( 0, 0, inum, 'fwfisf_par', fwfisf_par          )
            CALL iom_rstput( 0, 0, inum, 'rhisf_par_tbl', rhisf_tbl_par    )
            CALL iom_rstput( 0, 0, inum, 'rfrac_par_tbl', rfrac_tbl_par    )
            CALL iom_rstput( 0, 0, inum, 'misfkb_par', REAL(misfkb_par,wp) )
            CALL iom_rstput( 0, 0, inum, 'misfkt_par', REAL(misfkt_par,wp) )
            CALL iom_rstput( 0, 0, inum, 'mskisf_par', REAL(mskisf_par,wp), ktype = jp_i1 )
         END IF
      END IF
      !
      IF( ALLOCATED(ahtu) ) THEN
         CALL iom_rstput( 0, 0, inum,  'ahtu', ahtu              )    ! aht at u-point
         CALL iom_rstput( 0, 0, inum,  'ahtv', ahtv              )    ! aht at v-point
      ENDIF
      IF( ALLOCATED(ahmt) ) THEN
         CALL iom_rstput( 0, 0, inum,  'ahmt', ahmt              )    ! ahmt at u-point
         CALL iom_rstput( 0, 0, inum,  'ahmf', ahmf              )    ! ahmf at v-point
      ENDIF
      DO_2D( 0, 0, 0, 0 )
         z2d(ji,jj) = emp(ji,jj) - rnf(ji,jj)
      END_2D
      CALL iom_rstput( 0, 0, inum, 'sowaflup', z2d               )    ! freshwater budget
      DO_2D( 0, 0, 0, 0 )
         z2d(ji,jj) = qsr(ji,jj) + qns(ji,jj)
      END_2D
      CALL iom_rstput( 0, 0, inum, 'sohefldo', z2d               )    ! total heat flux
      CALL iom_rstput( 0, 0, inum, 'soshfldo', qsr               )    ! solar heat flux
      CALL iom_rstput( 0, 0, inum, 'soicecov', fr_i              )    ! ice fraction
      CALL iom_rstput( 0, 0, inum, 'sozotaux', utau              )    ! i-wind stress
      CALL iom_rstput( 0, 0, inum, 'sometauy', vtau              )    ! j-wind stress
      IF(  .NOT.ln_linssh  ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
           z3d(ji,jj,jk) = gdept(ji,jj,jk,Kmm)   ! 3D workspace for qco substitution
         END_3D
         CALL iom_rstput( 0, 0, inum, 'vovvldep', z3d            )    !  T-cell depth 
         DO_3D( 0, 0, 0, 0, 1, jpk )
           z3d(ji,jj,jk) = e3t(ji,jj,jk,Kmm)     ! 3D workspace for qco substitution
         END_3D
         CALL iom_rstput( 0, 0, inum, 'vovvle3t', z3d            )    !  T-cell thickness  
      END IF
      IF( ln_wave .AND. ln_sdw ) THEN
         CALL iom_rstput( 0, 0, inum, 'sdzocrtx', usd            )    ! now StokesDrift i-velocity
         CALL iom_rstput( 0, 0, inum, 'sdmecrty', vsd            )    ! now StokesDrift j-velocity
         CALL iom_rstput( 0, 0, inum, 'sdvecrtz', wsd            )    ! now StokesDrift k-velocity
      ENDIF
      IF ( ln_abl ) THEN
         CALL iom_rstput ( 0, 0, inum, "uz1_abl",   u_abl(:,:,2,nt_a  ) )   ! now first level i-wind
         CALL iom_rstput ( 0, 0, inum, "vz1_abl",   v_abl(:,:,2,nt_a  ) )   ! now first level j-wind
         CALL iom_rstput ( 0, 0, inum, "tz1_abl",  tq_abl(:,:,2,nt_a,1) )   ! now first level temperature
         CALL iom_rstput ( 0, 0, inum, "qz1_abl",  tq_abl(:,:,2,nt_a,2) )   ! now first level humidity
      ENDIF
      IF( ln_zdfosm ) THEN
         CALL iom_rstput( 0, 0, inum, 'hbl', hbl*tmask(:,:,1)  )      ! now boundary-layer depth
         CALL iom_rstput( 0, 0, inum, 'hml', hml*tmask(:,:,1)  )      ! now mixed-layer depth
         CALL iom_rstput( 0, 0, inum, 'avt_k', avt_k*wmask     )      ! w-level diffusion
         CALL iom_rstput( 0, 0, inum, 'avm_k', avm_k*wmask     )      ! now w-level viscosity
         CALL iom_rstput( 0, 0, inum, 'ghamt', ghamt*wmask     )      ! non-local t forcing
         CALL iom_rstput( 0, 0, inum, 'ghams', ghams*wmask     )      ! non-local s forcing
         CALL iom_rstput( 0, 0, inum, 'ghamu', ghamu*umask     )      ! non-local u forcing
         CALL iom_rstput( 0, 0, inum, 'ghamv', ghamv*vmask     )      ! non-local v forcing
         IF( ln_osm_mle ) THEN
            CALL iom_rstput( 0, 0, inum, 'hmle', hmle*tmask(:,:,1)  ) ! now transition-layer depth
         END IF
      ENDIF
      !
      CALL iom_close( inum )
      ! 
#if defined key_si3
      IF( nn_ice == 2 ) THEN   ! condition needed in case agrif + ice-model but no-ice in child grid
         CALL iom_open( TRIM(cdfile_name)//'_ice', inum, ldwrt = .TRUE., kdlev = jpl, cdcomp = 'ICE' )
         CALL ice_wri_state( inum )
         CALL iom_close( inum )
      ENDIF
      !
#endif
   END SUBROUTINE dia_wri_state

   !!======================================================================
END MODULE diawri
