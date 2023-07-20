MODULE ldfslp
   !!======================================================================
   !!                       ***  MODULE  ldfslp  ***
   !! Ocean physics: slopes of neutral surfaces
   !!======================================================================
   !! History :  OPA  ! 1994-12  (G. Madec, M. Imbard)  Original code
   !!            8.0  ! 1997-06  (G. Madec)  optimization, lbc
   !!            8.1  ! 1999-10  (A. Jouzeau)  NEW profile in the mixed layer
   !!   NEMO     1.0  ! 2002-10  (G. Madec)  Free form, F90
   !!             -   ! 2005-10  (A. Beckmann)  correction for s-coordinates
   !!            3.3  ! 2010-10  (G. Nurser, C. Harris, G. Madec)  add Griffies operator
   !!             -   ! 2010-11  (F. Dupond, G. Madec)  bug correction in slopes just below the ML
   !!            3.7  ! 2013-12  (F. Lemarie, G. Madec)  add limiter on triad slopes
   !!            4.2  ! 2023-01  (G. Madec)  add s-coordinate capability
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_slp       : calculates the slopes of neutral surface   (Madec operator)
   !!   ldf_slp_triad : calculates the triads of isoneutral slopes (Griffies operator)
   !!   ldf_slp_mxl   : calculates the slopes at the base of the mixed layer (Madec operator)
   !!   ldf_slp_init  : initialization of the slopes computation
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE isf_oce        ! ice shelf
   USE dom_oce        ! ocean space and time domain
   !USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE phycst         ! physical constants
   USE zdfmxl         ! mixed layer depth
   USE eosbn2         ! equation of states
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distribued memory computing library
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_slp         ! routine called by step.F90
   PUBLIC   ldf_slp_triad   ! routine called by step.F90
   PUBLIC   ldf_slp_init    ! routine called by nemogcm.F90

   LOGICAL , PUBLIC ::   l_ldfslp = .FALSE.     !: slopes flag

   LOGICAL , PUBLIC ::   ln_traldf_iso   = .TRUE.       !: iso-neutral direction                           (nam_traldf namelist)
   LOGICAL , PUBLIC ::   ln_traldf_triad = .FALSE.      !: griffies triad scheme                           (nam_traldf namelist)
   LOGICAL , PUBLIC ::   ln_dynldf_iso                  !: iso-neutral direction                           (nam_dynldf namelist)

   LOGICAL , PUBLIC ::   ln_triad_iso    = .FALSE.      !: pure horizontal mixing in ML                    (nam_traldf namelist)
   LOGICAL , PUBLIC ::   ln_botmix_triad = .FALSE.      !: mixing on bottom                                (nam_traldf namelist)
   REAL(wp), PUBLIC ::   rn_sw_triad     = 1._wp        !: =1 switching triads ; =0 all four triads used   (nam_traldf namelist)
   REAL(wp), PUBLIC ::   rn_slpmax       = 0.01_wp      !: slope limit                                     (nam_traldf namelist)

   LOGICAL , PUBLIC ::   ln_ldfeiv           !: eddy induced velocity flag

   LOGICAL , PUBLIC ::   l_grad_zps = .FALSE.           !: special treatment for Horz Tgradients w partial steps (triad operator)
   
   !                                                     !! Classic operator (Madec)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   uslp, wslpi          !: i_slope at U- and W-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   vslp, wslpj          !: j-slope at V- and W-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   eipsi_uw, eipsi_vw   !: eiv psi at uw- and vw-points
   !                                                     !! triad operator (Griffies)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   wslp2                !: wslp**2 from Griffies quarter cells
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) ::   triadi_g, triadj_g   !: skew flux  slopes relative to geopotentials
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) ::   triadi  , triadj     !: isoneutral slopes relative to model-coordinate
   !                                                     !! both operators
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   ah_wslp2             !: ah * slope^2 at w-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   akz                  !: stabilizing vertical diffusivity
   
   !                                                     !! Madec operator
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   omlmask           ! mask of the surface mixed layer at T-pt   
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   uslp_hml, wslpi_hml, uwslp_hml   ! i_slope / hml at U- and W-points just below the mixed layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   vslp_hml, wslpj_hml, vwslp_hml   ! i_slope / hml at U- and W-points just below the mixed layer

   REAL(wp) ::   repsln = 1.e-25_wp       ! tiny value used as minium of di(rho), dj(rho) and dk(rho)

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: ldfslp.F90 15062 2021-06-28 11:19:48Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_slp( kt, prd, pn2, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE ldf_slp  ***
      !!
      !! ** Purpose :   Compute the slopes of neutral surface (slope of isopycnal
      !!              surfaces referenced locally) (ln_traldf_iso=T).
      !!
      !! ** Method  :   The slope in the i-direction is computed at U- and
      !!      W-points (uslp, wslpi) and the slope in the j-direction is
      !!      computed at V- and W-points (vslp, wslpj).
      !!      They are bounded by 1/100 over the whole ocean, and within the
      !!      surface layer they are bounded by the distance to the surface
      !!      ( slope<= depth/l  where l is the length scale of horizontal
      !!      diffusion (here, aht=2000m2/s ==> l=20km with a typical velocity
      !!      of 10cm/s)
      !!        A horizontal shapiro filter is applied to the slopes
      !!        ln_sco=T, s-coordinate, add to the previously computed slopes
      !!      the slope of the model level surface.
      !!        macro-tasked on horizontal slab (jk-loop)  (2, jpk-1)
      !!      [slopes already set to zero at level 1, and to zero or the ocean
      !!      bottom slope (ln_sco=T) at level jpk in inildf]
      !!
      !! ** Action : - uslp, wslpi, and vslp, wslpj, the i- and  j-slopes
      !!               of now neutral surfaces at u-, w- and v- w-points, resp.
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in)                   ::   kt    ! ocean time-step index
      INTEGER , INTENT(in)                   ::   Kbb, Kmm   ! ocean time level indices
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   prd   ! in situ density
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   pn2   ! Brunt-Vaisala frequency (locally ref.)
      !!
      INTEGER  ::   ji , jj , jk    ! dummy loop indices
      INTEGER  ::   ii0, ii1        ! temporary integer
      INTEGER  ::   ij0, ij1        ! temporary integer
      REAL(wp) ::   zeps, z1_16, zcofw ! local scalars
      !
      REAL(wp) ::   zrho0_g, zrho0_2g, zzdzR, zzdzRi, zzdzrj, zzdiRw, zzdjRw        
      REAL(wp) ::   zsi_raw, zsi_coord, zsi_g_raw, zsi_g_lim, ze3_e1
      REAL(wp) ::   zsj_raw, zsj_coord, zsj_g_raw, zsj_g_lim, ze3_e2
      !
      REAL(wp) ::   zci, zfi, zau, zbu, zai, zbi   !   -      -
      REAL(wp) ::   zcj, zfj, zav, zbv, zaj, zbj   !   -      -
      REAL(wp) ::   zck, zfk, zcku, zckv, zbw      !   -      -
      REAL(wp) ::   zdepu, zdepuw, zdepw                  !   -      -
      REAL(wp) ::   zdepv, zdepvw                  !   -      -
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zgru, zwz, zdzr
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zgrv, zww
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ldf_slp')
      !
      zeps   =  1.e-20_wp           !==   Local constant initialization   ==!
      z1_16  =  1.0_wp / 16._wp
      zrho0_g  = rho0 / grav
      zrho0_2g = rho0 / ( 2._wp * grav ) 
      !
      zww(:,:,:) = 0._wp
      zwz(:,:,:) = 0._wp
      !
      DO_3D( 1, 0, 1, 0, 1, jpk )   !==   i- & j-gradient of density   ==!
         zgru(ji,jj,jk) = ( prd(ji+1,jj  ,jk) - prd(ji,jj,jk) ) * umask(ji,jj,jk)
         zgrv(ji,jj,jk) = ( prd(ji  ,jj+1,jk) - prd(ji,jj,jk) ) * vmask(ji,jj,jk)
      END_3D
      IF( ln_zps ) THEN                           ! partial steps correction at the bottom ocean level
         DO_2D( 1, 0, 1, 0 )
            zgru(ji,jj,mbku(ji,jj)) = gru(ji,jj)
            zgrv(ji,jj,mbkv(ji,jj)) = grv(ji,jj)
         END_2D
      ENDIF
      IF( ln_zps .AND. ln_isfcav ) THEN           ! partial steps correction at the bottom ocean level
         DO_2D( 1, 0, 1, 0 )
            IF( miku(ji,jj) > 1 )   zgru(ji,jj,miku(ji,jj)) = grui(ji,jj) 
            IF( mikv(ji,jj) > 1 )   zgrv(ji,jj,mikv(ji,jj)) = grvi(ji,jj)
         END_2D
      ENDIF
      !
      !                               !==   Local vertical density gradient at T-point   == !   (evaluated from N^2)
!!gm new

      zdzr(:,:,1) = zeps   ;   zdzr(:,:,jpk) = zeps
      DO jk = 2, jpkm1
         zdzr(:,:,jk) = - MAX( zeps , zrho0_2g * ( pn2(:,:,jk) + pn2(:,:,jk+1) ) ) 
      END DO

      !
      !                             !==   Slopes just below the mixed layer   ==!
      CALL ldf_slp_mxl( prd, pn2, zgru, zgrv, zdzr, Kmm )        ! output: uslp_hml, vslp_hml, wslpi_hml, wslpj_hml
      !                                                          !                             eipsi_uw , eipsi_vw

      ! ===========================
      ! I.  slopes at u and v point      | uslp = d/di( prd ) / d/dz( prd )
      ! ===========================      | vslp = d/dj( prd ) / d/dz( prd )
      !
      DO_3D( 0, 0, 0, 0, 2, jpkm1 )        !* Slopes at u and v points
         !
         !                                      ! raw slopes
         zsi_raw = zgru(ji,jj,jk) * r1_e1u(ji,jj) * 2._wp / ( zdzr(ji,jj,jk) + zdzr(ji+1,jj  ,jk) )
         zsj_raw = zgrv(ji,jj,jk) * r1_e2v(ji,jj) * 2._wp / ( zdzr(ji,jj,jk) + zdzr(ji  ,jj+1,jk) )
         !
         !                                      ! coordinate slopes reference to z=0
         zsi_coord = ( gdept_z0(ji+1,jj,jk,Kbb) - gdept_z0(ji,jj,jk,Kbb) ) * r1_e1u(ji,jj)
         zsj_coord = ( gdept_z0(ji,jj+1,jk,Kbb) - gdept_z0(ji,jj,jk,Kbb) ) * r1_e2v(ji,jj)
         !
         !                                      ! slopes referenced to geopot surfaces
         zsi_g_raw = zsi_raw - zsi_coord
         zsj_g_raw = zsj_raw - zsj_coord
         !                                      ! slope limiter
         !                                            ! limit required in bilap case and usefull in lap case
         ze3_e1    = e3u(ji,jj,jk,Kbb) * r1_e1u(ji,jj)
         ze3_e2    = e3v(ji,jj,jk,Kbb) * r1_e2v(ji,jj)
         !                                            ! NB: hard coded factor 5 (may be a namelist parameter...)
         zsi_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e1, ABS( zsi_g_raw ) ), zsi_g_raw )
         zsj_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e2, ABS( zsj_g_raw ) ), zsj_g_raw )
         !                                            ! Mixed Layer slope linear flattening
         zfi = MAX( omlmask(ji,jj,jk), omlmask(ji+1,jj  ,jk) )
         zfj = MAX( omlmask(ji,jj,jk), omlmask(ji  ,jj+1,jk) )
         !                                            ! thickness of water column between z=0 and level k at u/v point
#if defined key_isf
         !                                            ! ISF : depth of u/v points referenced to the base of iceshelves or z=0
         zdepu = 0.5_wp * ( (gdept_z0(ji,jj,jk,Kbb) + gdept_z0(ji+1,jj,jk,Kbb)) - MAX( risfdep(ji,jj) , risfdep(ji+1,jj) ) * 2._wp )
         zdepv = 0.5_wp * ( (gdept_z0(ji,jj,jk,Kbb) + gdept_z0(ji,jj+1,jk,Kbb)) - MAX( risfdep(ji,jj) , risfdep(ji,jj+1) ) * 2._wp )
#else
         !                                            ! no ISF: depth of u/v points referenced to z=0
         zdepu = 0.5_wp * (  gdept_z0(ji,jj,jk,Kbb) + gdept_z0(ji+1,jj,jk,Kbb)  )
         zdepv = 0.5_wp * (  gdept_z0(ji,jj,jk,Kbb) + gdept_z0(ji,jj+1,jk,Kbb)  )
#endif
         !
         zwz(ji,jj,jk) = ( ( 1._wp - zfi) * zsi_g_lim  +  zfi * zdepu * uslp_hml(ji,jj) ) * umask(ji,jj,jk)
         zww(ji,jj,jk) = ( ( 1._wp - zfj) * zsj_g_lim  +  zfj * zdepv * vslp_hml(ji,jj) ) * vmask(ji,jj,jk)
         !
      END_3D
      !
      CALL lbc_lnk( 'ldfslp', zwz, 'U', -1.0_wp,  zww, 'V', -1.0_wp )      ! lateral boundary conditions
      !
      !                                    !* horizontal Shapiro filter AND add the coordinate slopes
      DO jk = 2, jpkm1
         !
         DO_2D( 0, 0, 0, 0 )                 ! rows jj=2 and =jpjm1 only
            !                                      ! coordinate slopes reference to z=0
            zsi_coord = ( gdept_z0(ji+1,jj,jk,Kbb) - gdept_z0(ji,jj,jk,Kbb) ) * r1_e1u(ji,jj)
            zsj_coord = ( gdept_z0(ji,jj+1,jk,Kbb) - gdept_z0(ji,jj,jk,Kbb) ) * r1_e2v(ji,jj)
            !                                      ! decrease along coastal boundaries
            zcku = ( umask(ji,jj+1,jk) + umask(ji,jj-1,jk  ) )     &
               & * ( umask(ji,jj  ,jk) + umask(ji,jj  ,jk+1) ) * 0.25_wp
            zckv = ( vmask(ji+1,jj,jk) + vmask(ji-1,jj,jk  ) )     &
               & * ( vmask(ji  ,jj,jk) + vmask(ji  ,jj,jk+1) ) * 0.25_wp
            !
            uslp(ji,jj,jk) = zcku * z1_16 * (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)      &
               &                              +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)      &
               &                              + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)      &
               &                              +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )    &
               &                              + 4.*  zwz(ji  ,jj  ,jk)                        ) + zsi_coord
            vslp(ji,jj,jk) = zckv * z1_16 * (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)      &
               &                              +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)      &
               &                              + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)      &
               &                              +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )    &
               &                              + 4.*  zww(ji,jj    ,jk)                        ) + zsj_coord
         END_2D
         !
      END DO


      ! ===========================
      ! II.  slopes at w point           | wslpi = mik( d/di( prd ) ) / d/dz( prd )
      ! ===========================      | wslpj = mjk( d/dj( prd ) ) / d/dz( prd )
      !
      DO_3D( 0, 0, 0, 0, 2, jpkm1 )
         !                                  !* Local vertical density gradient evaluated from N^2
         zzdzR = - MAX( zeps , zrho0_g * pn2(ji,jj,jk) )
         !
         !                                  !* Slopes at w point
         !                                        ! i- & j-slopes at w-points
         zci = MAX(  umask(ji-1,jj,jk  ) + umask(ji,jj,jk  )           &
            &      + umask(ji-1,jj,jk-1) + umask(ji,jj,jk-1) , zeps  ) * zzdzR * e1t(ji,jj)
         zcj = MAX(  vmask(ji,jj-1,jk  ) + vmask(ji,jj,jk-1)           &
            &      + vmask(ji,jj-1,jk-1) + vmask(ji,jj,jk  ) , zeps  ) * zzdzR * e2t(ji,jj)
         !
         zsi_raw = (  zgru(ji-1,jj,jk  ) +  zgru(ji,jj,jk-1)           &
            &       + zgru(ji-1,jj,jk-1) +  zgru(ji,jj,jk  )   ) / zci * wmask (ji,jj,jk)
         zsj_raw = (  zgrv(ji,jj-1,jk  ) +  zgrv(ji,jj,jk-1)           &
            &       + zgrv(ji,jj-1,jk-1) +  zgrv(ji,jj,jk  )   ) / zcj * wmask (ji,jj,jk)
         !
         !                                  !* coordinate slopes reference to z=0 at w-point
!!gm ISF case ??  think about it !
         zsi_coord = 0.5_wp * ( gdepw_z0(ji+1,jj,jk,Kbb) - gdepw_z0(ji-1,jj,jk,Kbb) ) * r1_e1t(ji,jj)
         zsj_coord = 0.5_wp * ( gdepw_z0(ji,jj+1,jk,Kbb) - gdepw_z0(ji,jj-1,jk,Kbb) ) * r1_e2t(ji,jj)
         !
         !                                  !* slopes referenced to geopot surfaces
         zsi_g_raw = zsi_raw - zsi_coord
         zsj_g_raw = zsj_raw - zsj_coord
         !                                      ! slope limiter
         !                                            ! limit required in bilap case and usefull in lap case
         ze3_e1    = e3w(ji,jj,jk,Kbb) * r1_e1u(ji,jj)
         ze3_e2    = e3w(ji,jj,jk,Kbb) * r1_e2v(ji,jj)
         !                                            ! NB: hard coded factor 5 (may be a namelist parameter...)
         zsi_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e1, ABS( zsi_g_raw ) ), zsi_g_raw )
         zsj_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e2, ABS( zsj_g_raw ) ), zsj_g_raw )
         !                                            ! Mixed Layer slope linear flattening
         zfk = MAX( omlmask(ji,jj,jk), omlmask(ji,jj,jk-1) )   ! zfk=1 in the ML otherwise zfk=0
#if defined key_isf
         zdepw = ( gdepw(ji,jj,jk,Kbb) - gdepw(ji,jj,mikt(ji,jj),Kbb) ) )
#else
         zdepw =   gdepw_z0(ji,jj,jk,Kbb)
#endif
         zwz(ji,jj,jk) = (  zsi_g_lim * ( 1._wp - zfk ) + zdepw * wslpi_hml(ji,jj) * zfk  ) * wmask(ji,jj,jk)
         zww(ji,jj,jk) = (  zsj_g_lim * ( 1._wp - zfk ) + zdepw * wslpj_hml(ji,jj) * zfk  ) * wmask(ji,jj,jk)
         !
      END_3D
      !
      CALL lbc_lnk( 'ldfslp', zwz, 'T', -1.0_wp,  zww, 'T', -1.0_wp )      ! lateral boundary conditions
      !
      !                                           !* horizontal Shapiro filter AND add the coordinate slopes
      DO jk = 2, jpkm1
         !
         DO_2D( 0, 0, 0, 0 )                      ! rows jj=2 and =jpjm1 only
            !                                         ! coordinate slopes reference to z=0 at w-point
            zsi_coord = 0.5_wp * ( gdepw_z0(ji+1,jj,jk,Kbb) - gdepw_z0(ji-1,jj,jk,Kbb) ) * r1_e1t(ji,jj)
            zsj_coord = 0.5_wp * ( gdepw_z0(ji,jj+1,jk,Kbb) - gdepw_z0(ji,jj-1,jk,Kbb) ) * r1_e2t(ji,jj)
            !                                         ! decrease in vicinity of topography
            zck =   ( umask(ji,jj,jk) + umask(ji-1,jj,jk) )   &
               &  * ( vmask(ji,jj,jk) + vmask(ji,jj-1,jk) ) * 0.25
            zcofw = z1_16 * zck * wmask(ji,jj,jk)
            wslpi(ji,jj,jk) = (         zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)     &
                 &               +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)     &
                 &               + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)     &
                 &               +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )   &
                 &               + 4.*  zwz(ji  ,jj  ,jk)                         ) * zcofw  + zsi_coord

            wslpj(ji,jj,jk) = (         zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)     &
                 &               +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)     &
                 &               + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)     &
                 &               +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )   &
                 &               + 4.*  zww(ji  ,jj  ,jk)                         ) * zcofw  + zsj_coord
         END_2D
         !
      END DO


      ! =================================
      ! III.  Eddy Induced streamfunction at uw/vw points     (computed from slopes referenced to z=0  WITHOUT horizontal filter)
      ! =================================----------------     | eipsi_uw = mk[zdxR] / mi[d/dzR]
      !                                                       ! eipsi_vw = mk[zdyR] / mj[d/dzR]
      IF( ln_ldfeiv ) THEN
         !
         DO_3D( 0, 0, 0, 0, 2, jpkm1 )
            !
            !                                  !* Local vertical density gradient evaluated from N^2
            zzdzRi    = - MAX( zeps , zrho0_g * ( pn2(ji,jj,jk) + pn2(ji+1,jj,jk) ) )
            zzdzRj    = - MAX( zeps , zrho0_g * ( pn2(ji,jj,jk) + pn2(ji,jj+1,jk) ) )
            !
            !                                  !* local iso-level density gradient
            zzdiRw    =  ( zgru(ji+1,jj,jk) + zgru(ji,jj,jk) ) * r1_e1u(ji,jj) * wumask(ji,jj,jk)
            zzdjRw    =  ( zgrv(ji,jj+1,jk) + zgrv(ji,jj,jk) ) * r1_e2v(ji,jj) * wvmask(ji,jj,jk)
            !
            !                                  !* raw slopes
            zsi_raw   = zzdiRw / zzdzRi
            zsj_raw   = zzdjRw / zzdzRj
            !
            !                                  !* coordinate slopes at w-level reference to z=0
            zsi_coord = ( gdepw_z0(ji+1,jj,jk,Kbb) - gdepw_z0(ji,jj,jk,Kbb) ) * r1_e1u(ji,jj)
            zsj_coord = ( gdepw_z0(ji,jj+1,jk,Kbb) - gdepw_z0(ji,jj,jk,Kbb) ) * r1_e2v(ji,jj)
            !
            !                                  !* slopes at uw/vw-points referenced to geopot surfaces
            zsi_g_raw = zsi_raw - zsi_coord
            zsj_g_raw = zsj_raw - zsj_coord
            !
            !                                  !* slope limiter
            !                                         ! limit required in bilap case and usefull in lap case
            ze3_e1    = e3w(ji,jj,jk,Kbb) * r1_e1u(ji,jj)
            ze3_e2    = e3w(ji,jj,jk,Kbb) * r1_e2v(ji,jj)
            !                                         ! NB: hard coded factor 5 (may be a namelist parameter...)
            zsi_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e1, ABS( zsi_g_raw ) ), zsi_g_raw )
            zsj_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e2, ABS( zsj_g_raw ) ), zsj_g_raw )
            !                                         !
            !                                         ! Mixed Layer slope linear flattening
            zfi = MAX( omlmask(ji+1,jj  ,jk-1), omlmask(ji,jj,jk-1) )   ! zfi/j =1 in the ML
            zfj = MAX( omlmask(ji  ,jj+1,jk-1), omlmask(ji,jj,jk-1) )   !       =0 at the ML base and below
            !
            !                                            ! thickness of water column between z=0 and level k at uw/vw point
#if defined key_isf
            !                                            ! ISF : depth of u/v points referenced to the base of iceshelves or z=0
            zdepuw = 0.5_wp * (  (   gdepw_z0(ji,jj,jk,Kbb) + gdepw_z0(ji+1,jj,jk,Kbb)  )  &
               &               - MAX( risfdep(ji,jj)        ,  risfdep(ji+1,jj) ) * 2._wp  )
            zdepvw = 0.5_wp * (  (   gdepw_z0(ji,jj,jk,Kbb) + gdepw_z0(ji,jj+1,jk,Kbb) )   &
               &               - MAX( risfdep(ji,jj)        ,  risfdep(ji,jj+1) ) * 2._wp  )   
#else
            !                                            ! no ISF: depth of u/v points referenced to z=0
            zdepuw = 0.5_wp * ( gdepw_z0(ji,jj,jk,Kbb) + gdepw_z0(ji+1,jj  ,jk,Kbb) )
            zdepvw = 0.5_wp * ( gdepw_z0(ji,jj,jk,Kbb) + gdepw_z0(ji  ,jj+1,jk,Kbb) )
#endif
            !                                            ! Eddy Induced streamfunction (without * aei : done in traadv) 
            eipsi_uw(ji,jj,jk) = (  (1._wp - zfi) * zsi_g_lim + zfi * zdepuw * uwslp_hml(ji,jj)  ) * wumask(ji,jj,jk)
            eipsi_vw(ji,jj,jk) = (  (1._wp - zfj) * zsj_g_lim + zfj * zdepvw * vwslp_hml(ji,jj)  ) * wvmask(ji,jj,jk)
            !                                            ! Eddy Induced streamfunction (limited slope * aei) 
!             eipsi_uw(ji,jj,jk) = (  (1._wp - zfi) * zsi_g_lim + zfi * zdepuw * uwslp_hml(ji,jj)  ) * aeiu(ji,jj,1) * wumask(ji,jj,jk)
!             eipsi_vw(ji,jj,jk) = (  (1._wp - zfj) * zsj_g_lim + zfj * zdepvw * vwslp_hml(ji,jj)  ) * aeiv(ji,jj,1) * wvmask(ji,jj,jk)
            !
         END_3D
         !
      ENDIF


      ! IV. Lateral boundary conditions
      ! ===============================
      IF( ln_ldfeiv ) THEN
         CALL lbc_lnk( 'ldfslp', uslp    , 'U', -1.0_wp , vslp    , 'V', -1.0_wp , &
            &                    wslpi   , 'W', -1.0_wp , wslpj   , 'W', -1.0_wp , &
            &                    eipsi_uw, 'U', -1.0_wp,  eipsi_vw, 'V', -1.0_wp   )
      ELSE
         CALL lbc_lnk( 'ldfslp', uslp    , 'U', -1.0_wp , vslp    , 'V', -1.0_wp , &
            &                    wslpi   , 'W', -1.0_wp , wslpj   , 'W', -1.0_wp   )
      ENDIF

      IF(sn_cfctl%l_prtctl) THEN
         CALL prt_ctl(tab3d_1=uslp    , clinfo1=' slp  - u : '    , tab3d_2=vslp    , clinfo2=' v     : ')
         CALL prt_ctl(tab3d_1=wslpi   , clinfo1=' slp  - wi: '    , tab3d_2=wslpj   , clinfo2=' wj    : ')
         CALL prt_ctl(tab3d_1=eipsi_uw, clinfo1=' slp  - psi_uw: ', tab3d_2=eipsi_uw, clinfo2=' psi_vw: ')
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('ldf_slp')
      !
   END SUBROUTINE ldf_slp


   SUBROUTINE ldf_slp_triad ( kt, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE ldf_slp_triad  ***
      !!
      !! ** Purpose :   Compute the squared slopes of neutral surfaces (slope
      !!      of iso-pycnal surfaces referenced locally) (ln_traldf_triad=T)
      !!      at W-points using the Griffies quarter-cells.
      !!
      !! ** Method  :   calculates alpha and beta at T-points
      !!
      !! ** Action : - triadi_g, triadj_g   T-pts i- and j-slope triads relative to geopot. (used for eiv)
      !!             - triadi , triadj    T-pts i- and j-slope triads relative to model-coordinate
      !!             - wslp2              squared slope of neutral surfaces at w-points.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt             ! ocean time-step index
      INTEGER , INTENT(in)  ::   Kbb, Kmm       ! ocean time level indices
      !!
      INTEGER  ::   ji, jj, jk, jl, ip, jp, kp  ! dummy loop indices
      INTEGER  ::   iku, ikv                    ! local integer
      REAL(wp) ::   zfacti, zfactj              ! local scalars
      REAL(wp) ::   znot_thru_surface           ! local scalars
      REAL(wp) ::   zdit, zdis, zdkt, zbu, zbti, zisw
      REAL(wp) ::   zdjt, zdjs, zdks, zbv, zbtj, zjsw
      REAL(wp) ::   zdxrho_raw, zti_coord, zti_raw, zti_lim, zti_g_raw, zti_g_lim
      REAL(wp) ::   zdyrho_raw, ztj_coord, ztj_raw, ztj_lim, ztj_g_raw, ztj_g_lim
      REAL(wp) ::   zdzrho_raw
      REAL(wp) ::   zbeta0, ze3_e1, ze3_e2
      REAL(wp), DIMENSION(jpi,jpj)     ::   z1_mlbw
      REAL(wp), DIMENSION(jpi,jpj,jpk,0:1) ::   zdxrho , zdyrho, zdzrho     ! Horizontal and vertical density gradients
      REAL(wp), DIMENSION(jpi,jpj,0:1,0:1) ::   zti_mlb, ztj_mlb            ! for Griffies operator only
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ldf_slp_triad')
      !
      !--------------------------------!
      !  Some preliminary calculation  !
      !--------------------------------!
      !
      DO jl = 0, 1                            !==  unmasked before density i- j-, k-gradients  ==!
         !
         ip = jl   ;   jp = jl                ! guaranteed nonzero gradients ( absolute value larger than repsln)
         DO_3D( nn_hls, nn_hls-1, nn_hls, nn_hls-1, 1, jpkm1 )        ! done each pair of triad ! NB: not masked ==>  a minimum value is set
            zdit = ( ts(ji+1,jj,jk,jp_tem,Kbb) - ts(ji,jj,jk,jp_tem,Kbb) )    ! i-gradient of T & S at u-point
            zdis = ( ts(ji+1,jj,jk,jp_sal,Kbb) - ts(ji,jj,jk,jp_sal,Kbb) )
            zdjt = ( ts(ji,jj+1,jk,jp_tem,Kbb) - ts(ji,jj,jk,jp_tem,Kbb) )    ! j-gradient of T & S at v-point
            zdjs = ( ts(ji,jj+1,jk,jp_sal,Kbb) - ts(ji,jj,jk,jp_sal,Kbb) )
            zdxrho_raw = ( - rab_b(ji+ip,jj   ,jk,jp_tem) * zdit + rab_b(ji+ip,jj   ,jk,jp_sal) * zdis ) * r1_e1u(ji,jj)
            zdyrho_raw = ( - rab_b(ji   ,jj+jp,jk,jp_tem) * zdjt + rab_b(ji   ,jj+jp,jk,jp_sal) * zdjs ) * r1_e2v(ji,jj)
            zdxrho(ji+ip,jj   ,jk,1-ip) = SIGN(  MAX( repsln, ABS( zdxrho_raw ) ), zdxrho_raw  )   ! keep the sign
            zdyrho(ji   ,jj+jp,jk,1-jp) = SIGN(  MAX( repsln, ABS( zdyrho_raw ) ), zdyrho_raw  )
         END_3D
         !
         IF( ln_zps .AND. l_grad_zps ) THEN     ! partial steps: correction of i- & j-grad on bottom
            DO_2D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 )
               iku  = mbku(ji,jj)          ;   ikv  = mbkv(ji,jj)             ! last ocean level (u- & v-points)
               zdit = gtsu(ji,jj,jp_tem)   ;   zdjt = gtsv(ji,jj,jp_tem)      ! i- & j-gradient of Temperature
               zdis = gtsu(ji,jj,jp_sal)   ;   zdjs = gtsv(ji,jj,jp_sal)      ! i- & j-gradient of Salinity
               zdxrho_raw = ( - rab_b(ji+ip,jj   ,iku,jp_tem) * zdit + rab_b(ji+ip,jj   ,iku,jp_sal) * zdis ) * r1_e1u(ji,jj)
               zdyrho_raw = ( - rab_b(ji   ,jj+jp,ikv,jp_tem) * zdjt + rab_b(ji   ,jj+jp,ikv,jp_sal) * zdjs ) * r1_e2v(ji,jj)
               zdxrho(ji+ip,jj   ,iku,1-ip) = SIGN( MAX( repsln, ABS( zdxrho_raw ) ), zdxrho_raw )   ! keep the sign
               zdyrho(ji   ,jj+jp,ikv,1-jp) = SIGN( MAX( repsln, ABS( zdyrho_raw ) ), zdyrho_raw )
            END_2D
         ENDIF
         !
      END DO

      DO kp = 0, 1                            !==  unmasked before density i- j-, k-gradients  ==!
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpkm1 )      ! done each pair of triad ! NB: not masked ==>  a minimum value is set
            IF( jk+kp > 1 ) THEN              ! k-gradient of T & S a jk+kp
               zdkt = ( ts(ji,jj,jk+kp-1,jp_tem,Kbb) - ts(ji,jj,jk+kp,jp_tem,Kbb) )
               zdks = ( ts(ji,jj,jk+kp-1,jp_sal,Kbb) - ts(ji,jj,jk+kp,jp_sal,Kbb) )
            ELSE
               zdkt = 0._wp                                             ! 1st level gradient set to zero
               zdks = 0._wp
            ENDIF
            zdzrho_raw = ( - rab_b(ji,jj,jk   ,jp_tem) * zdkt & 
                       &   + rab_b(ji,jj,jk   ,jp_sal) * zdks &
                       & ) / e3w(ji,jj,jk+kp,Kmm)  
            zdzrho(ji,jj,jk,kp) = - MIN( - repsln , zdzrho_raw )    ! force zdzrho >= repsln
         END_3D
      END DO
      !
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )                   !== Reciprocal depth of the w-point below ML base  ==!
         jk = MIN( nmln(ji,jj), mbkt(ji,jj) ) + 1     ! MIN in case ML depth is the ocean depth
         z1_mlbw(ji,jj) = 1._wp / gdepw(ji,jj,jk,Kmm)
      END_2D
      !
      !                                       !==  intialisations to zero  ==!
      !
      wslp2  (:,:,:)     = 0._wp              ! wslp2 will be cumulated 3D field set to zero
      triadi_g(:,:,1,:,:) = 0._wp   ;   triadi_g(:,:,jpk,:,:) = 0._wp   ! set surface and bottom slope to zero
      triadj_g(:,:,1,:,:) = 0._wp   ;   triadj_g(:,:,jpk,:,:) = 0._wp
      !!gm _iso set to zero missing
      triadi  (:,:,1,:,:) = 0._wp   ;   triadj  (:,:,jpk,:,:) = 0._wp   ! set surface and bottom slope to zero
      triadj  (:,:,1,:,:) = 0._wp   ;   triadj  (:,:,jpk,:,:) = 0._wp

      !-------------------------------------!
      !  Triads just below the Mixed Layer  !
      !-------------------------------------!
      !
      DO jl = 0, 1                            ! calculate slope of the 4 triads immediately ONE level below mixed-layer base
         DO kp = 0, 1                         ! with only the slope-max limit   and   MASKED
            DO_2D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 )
               ip = jl   ;   jp = jl
               !
               jk = nmln(ji+ip,jj) + 1
               IF( jk > mbkt(ji+ip,jj) ) THEN   ! ML reaches bottom
                  zti_mlb(ji+ip,jj   ,1-ip,kp) = 0._wp
               ELSE                             
                  ! Add s-coordinate slope at t-points (do this by *subtracting* gradient of depth)
                  zti_g_raw = (  zdxrho(ji+ip,jj,jk-kp,1-ip) / zdzrho(ji+ip,jj,jk-kp,kp)      &
                     &          - ( gdept(ji+1,jj,jk-kp,Kmm) - gdept(ji,jj,jk-kp,Kmm) ) * r1_e1u(ji,jj)  ) * umask(ji,jj,jk)
                  ze3_e1    =  e3w(ji+ip,jj,jk-kp,Kmm) * r1_e1u(ji,jj) 
                  zti_mlb(ji+ip,jj   ,1-ip,kp) = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e1  , ABS( zti_g_raw ) ), zti_g_raw )
               ENDIF
               !
               jk = nmln(ji,jj+jp) + 1
               IF( jk >  mbkt(ji,jj+jp) ) THEN  !ML reaches bottom
                  ztj_mlb(ji   ,jj+jp,1-jp,kp) = 0._wp
               ELSE
                  ztj_g_raw = (  zdyrho(ji,jj+jp,jk-kp,1-jp) / zdzrho(ji,jj+jp,jk-kp,kp)      &
                     &      - ( gdept(ji,jj+1,jk-kp,Kmm) - gdept(ji,jj,jk-kp,Kmm) ) / e2v(ji,jj)  ) * vmask(ji,jj,jk)
                  ze3_e2    =  e3w(ji,jj+jp,jk-kp,Kmm) / e2v(ji,jj)
                  ztj_mlb(ji   ,jj+jp,1-jp,kp) = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e2  , ABS( ztj_g_raw ) ), ztj_g_raw )
               ENDIF
            END_2D
         END DO
      END DO

      !-------------------------------------!
      !  Triads with surface limits         !
      !-------------------------------------!
      !
      DO kp = 0, 1                            ! k-index of triads
         DO jl = 0, 1
            ip = jl   ;   jp = jl             ! i- and j-indices of triads (i-k and j-k planes)
            DO jk = 1, jpkm1
               ! Must mask contribution to slope from dz/dx at constant s for triads jk=1,kp=0 that poke up though ocean surface
               znot_thru_surface = REAL( 1-1/(jk+kp), wp )  !jk+kp=1,=0.; otherwise=1.0
               DO_2D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 )
                  !
                  ! Calculate slope relative to geopotentials used for GM skew fluxes
                  ! Add s-coordinate slope at t-points (do this by *subtracting* gradient of depth)
                  ! Limit by slope *relative to geopotentials* by rn_slpmax, and mask by psi-point
                  ! masked by umask taken at the level of dz(rho)
                  !
                  ! raw slopes: unmasked unbounded slopes (relative to geopotential (zti_g) and model surface (zti)
                  !
                  zti_raw   = zdxrho(ji+ip,jj   ,jk,1-ip) / zdzrho(ji+ip,jj   ,jk,kp)                   ! unmasked
                  ztj_raw   = zdyrho(ji   ,jj+jp,jk,1-jp) / zdzrho(ji   ,jj+jp,jk,kp)
                  !
                  ! Must mask contribution to slope for triad jk=1,kp=0 that poke up though ocean surface
                  zti_coord = znot_thru_surface * ( gdept(ji+1,jj  ,jk,Kmm) - gdept(ji,jj,jk,Kmm) ) * r1_e1u(ji,jj)
                  ztj_coord = znot_thru_surface * ( gdept(ji  ,jj+1,jk,Kmm) - gdept(ji,jj,jk,Kmm) ) * r1_e2v(ji,jj)     ! unmasked
                  zti_g_raw = zti_raw - zti_coord      ! ref to geopot surfaces
                  ztj_g_raw = ztj_raw - ztj_coord
                  ! additional limit required in bilaplacian case and usefull in laplacian case
                  ze3_e1    = e3w(ji+ip,jj   ,jk+kp,Kmm) * r1_e1u(ji,jj)
                  ze3_e2    = e3w(ji   ,jj+jp,jk+kp,Kmm) * r1_e2v(ji,jj)
                  ! NB: hard coded factor 5 (may be a namelist parameter...)
                  zti_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e1, ABS( zti_g_raw ) ), zti_g_raw )
                  ztj_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e2, ABS( ztj_g_raw ) ), ztj_g_raw )
                  !
                  ! Below  ML use limited zti_g as is & mask
                  ! Inside ML replace by linearly reducing sx_mlb towards surface & mask
                  !
                  zfacti = REAL( 1 - 1/(1 + (jk+kp-1)/nmln(ji+ip,jj)), wp )  ! k index of uppermost point(s) of triad is jk+kp-1
                  zfactj = REAL( 1 - 1/(1 + (jk+kp-1)/nmln(ji,jj+jp)), wp )  ! must be .ge. nmln(ji,jj) for zfact=1
                  !                                                          !                   otherwise  zfact=0
                  zti_g_lim =          ( zfacti   * zti_g_lim                       &
                     &      + ( 1._wp - zfacti ) * zti_mlb(ji+ip,jj,1-ip,kp)   &
                     &                           * gdepw(ji+ip,jj,jk+kp,Kmm) * z1_mlbw(ji+ip,jj) ) * umask(ji,jj,jk+kp)
                  ztj_g_lim =          ( zfactj   * ztj_g_lim                       &
                     &      + ( 1._wp - zfactj ) * ztj_mlb(ji,jj+jp,1-jp,kp)   &
                     &                           * gdepw(ji,jj+jp,jk+kp,Kmm) * z1_mlbw(ji,jj+jp) ) * vmask(ji,jj,jk+kp)
                  !
                  triadi_g(ji+ip,jj   ,jk,1-ip,kp) = zti_g_lim
                  triadj_g(ji   ,jj+jp,jk,1-jp,kp) = ztj_g_lim
                  !
                  ! Get coefficients of isoneutral diffusion tensor
                  ! 1. Utilise gradients *relative* to s-coordinate, so add t-point slopes (*subtract* depth gradients)
                  ! 2. We require that isoneutral diffusion  gives no vertical buoyancy flux
                  !     i.e. 33 term = (real slope* 31, 13 terms)
                  ! To do this, retain limited sx**2  in vertical flux, but divide by real slope for 13/31 terms
                  ! Equivalent to tapering A_iso = sx_limited**2/(real slope)**2
                  !
                  zti_lim  = ( zti_g_lim + zti_coord ) * umask(ji,jj,jk+kp)    ! remove coordinate slope => relative to coordinate surfaces
                  ztj_lim  = ( ztj_g_lim + ztj_coord ) * vmask(ji,jj,jk+kp)
                  !
                  IF( ln_triad_iso ) THEN
                     zti_raw = zti_lim*zti_lim / zti_raw
                     ztj_raw = ztj_lim*ztj_lim / ztj_raw
                     zti_raw = SIGN( MIN( ABS(zti_lim), ABS( zti_raw ) ), zti_raw )
                     ztj_raw = SIGN( MIN( ABS(ztj_lim), ABS( ztj_raw ) ), ztj_raw )
                     zti_lim = zfacti * zti_lim + ( 1._wp - zfacti ) * zti_raw
                     ztj_lim = zfactj * ztj_lim + ( 1._wp - zfactj ) * ztj_raw
                  ENDIF
                  !                                      ! switching triad scheme 
                  zisw = (1._wp - rn_sw_triad ) + rn_sw_triad    &
                     &            * 2._wp * ABS( 0.5_wp - kp - ( 0.5_wp - ip ) * SIGN( 1._wp , zdxrho(ji+ip,jj,jk,1-ip) )  )
                  zjsw = (1._wp - rn_sw_triad ) + rn_sw_triad    &
                     &            * 2._wp * ABS( 0.5_wp - kp - ( 0.5_wp - jp ) * SIGN( 1._wp , zdyrho(ji,jj+jp,jk,1-jp) )  )
                  !
                  triadi(ji+ip,jj   ,jk,1-ip,kp) = zti_lim * zisw
                  triadj(ji   ,jj+jp,jk,1-jp,kp) = ztj_lim * zjsw
                  !
                  zbu  = e1e2u(ji   ,jj   ) * e3u(ji   ,jj   ,jk   ,Kmm)
                  zbv  = e1e2v(ji   ,jj   ) * e3v(ji   ,jj   ,jk   ,Kmm)
                  zbti = e1e2t(ji+ip,jj   ) * e3w(ji+ip,jj   ,jk+kp,Kmm)
                  zbtj = e1e2t(ji   ,jj+jp) * e3w(ji   ,jj+jp,jk+kp,Kmm)
                  !
                  wslp2(ji+ip,jj,jk+kp) = wslp2(ji+ip,jj,jk+kp) + 0.25_wp * zbu / zbti * zti_g_lim*zti_g_lim      ! masked
                  wslp2(ji,jj+jp,jk+kp) = wslp2(ji,jj+jp,jk+kp) + 0.25_wp * zbv / zbtj * ztj_g_lim*ztj_g_lim
               END_2D
            END DO
         END DO
      END DO
      !
      wslp2(:,:,1) = 0._wp                ! force the surface wslp to zero

      CALL lbc_lnk( 'ldfslp', wslp2, 'W', 1.0_wp )      ! lateral boundary confition on wslp2 only   ==>>> gm : necessary ? to be checked
      !
      IF( ln_timing )   CALL timing_stop('ldf_slp_triad')
      !
   END SUBROUTINE ldf_slp_triad


   SUBROUTINE ldf_slp_mxl( prd, pn2, p_gru, p_grv, p_dzr, Kbb )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_slp_mxl  ***
      !!
      !! ** Purpose :   Compute the slopes of iso-neutral surface just below
      !!              the mixed layer.
      !!
      !! ** Method  :   The slope in the i-direction is computed at u- & w-points
      !!              (uslp_hml, wslpiml) and the slope in the j-direction is computed
      !!              at v- and w-points (vslp_hml, wslpjml) with the same bounds as
      !!              in ldf_slp.
      !!
      !! ** Action  :   uslp_hml, wslpiml :  i- &  j-slopes of neutral surfaces
      !!                vslp_hml, wslpjml    just below the mixed layer
      !!                omlmask         :  mixed layer mask
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   prd            ! in situ density
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   pn2            ! Brunt-Vaisala frequency (locally ref.)
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   p_gru, p_grv   ! i- & j-gradient of density (u- & v-pts)
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   p_dzr          ! z-gradient of density      (T-point)
      INTEGER , INTENT(in)                   ::   Kbb            ! ocean time level indices
      !!
      INTEGER  ::   ji , jj , jk                   ! dummy loop indices
      INTEGER  ::   iku, ikv, ik, ikm1             ! local integers
      REAL(wp) ::   zci, zfi, zau, zbu, zai, zbi   !   -      -
      REAL(wp) ::   zcj, zfj, zav, zbv, zaj, zbj   !   -      -
      REAL(wp) ::   zck, zfk,      zbw             !   -      -
      !
      REAL(wp) ::   zeps                           ! local scalars
      REAL(wp) ::   zrho0_g, zrho0_2g, zzdzR, zzdzRi, zzdzrj, zzdiRw, zzdjRw
      REAL(wp) ::   zhml_w, zhml_uw, zhml_vw
      REAL(wp) ::   zdepu, zdepv, zdepuw, zdepvw      
      REAL(wp) ::   zsi_raw, zsi_coord, zsi_g_raw, zsi_g_lim, ze3_e1
      REAL(wp) ::   zsj_raw, zsj_coord, zsj_g_raw, zsj_g_lim, ze3_e2
      !!----------------------------------------------------------------------
      !
      zeps   =  1.e-20_wp        !==   Local constant initialization   ==!
      zrho0_g  = rho0 / grav
      zrho0_2g = rho0 / ( 2._wp * grav ) 
      !
      uslp_hml (1,:) = 0._wp      ;      uslp_hml (jpi,:) = 0._wp
      vslp_hml (1,:) = 0._wp      ;      vslp_hml (jpi,:) = 0._wp
      wslpi_hml(1,:) = 0._wp      ;      wslpi_hml(jpi,:) = 0._wp
      wslpj_hml(1,:) = 0._wp      ;      wslpj_hml(jpi,:) = 0._wp
      !
      !                          !==  mixed layer mask  ==!
      !
      DO_3D( 1, 1, 1, 1, 1, jpk )                  ! =1 inside the mixed layer, =0 otherwise
         ik = nmln(ji,jj) - 1
         IF( jk <= ik ) THEN   ;   omlmask(ji,jj,jk) = 1._wp
         ELSE                  ;   omlmask(ji,jj,jk) = 0._wp
         ENDIF
      END_3D


      !                          !==  Slopes just below the bottom of mixed layer  ==!
      ! --------------------------------------------------------------
      ! The slope are computed as in the 3D case.
      ! A key point here is the definition of the mixed layer at u- and v-points.
      ! It is assumed to be the maximum of the two neighbouring T-point mixed layer depth.
      ! Otherwise, a n2 value inside the mixed layer can be involved in the computation
      ! of the slope, resulting in a too steep diagnosed slope and thus a spurious eddy
      ! induce velocity field near the base of the mixed layer.
      !-----------------------------------------------------------------------
      !
      DO_2D( 0, 0, 0, 0 )         !==   Slope/MLD at u- & v-points ==!   NB: slope just below the Mixed Layer : jk=mi,j[nmln]
         !
         !                                       !* Mixed layer level at u- & v-point
         iku = MIN(  MAX(  2, nmln(ji,jj) , nmln(ji+1,jj) ) , jpkm1  )
         ikv = MIN(  MAX(  2, nmln(ji,jj) , nmln(ji,jj+1) ) , jpkm1  )

         !                                      !* raw slopes
         zsi_raw = p_gru(ji,jj,iku) * r1_e1u(ji,jj) * 2._wp / ( p_dzr(ji,jj,iku) + p_dzr(ji+1,jj  ,iku) )
         zsj_raw = p_grv(ji,jj,ikv) * r1_e2v(ji,jj) * 2._wp / ( p_dzr(ji,jj,ikv) + p_dzr(ji  ,jj+1,ikv) )
         !
         !                                      !* coordinate slopes reference to z=0
         zsi_coord = ( gdept_z0(ji+1,jj,iku,Kbb) - gdept_z0(ji,jj,iku,Kbb) ) * r1_e1u(ji,jj)
         zsj_coord = ( gdept_z0(ji,jj+1,ikv,Kbb) - gdept_z0(ji,jj,ikv,Kbb) ) * r1_e2v(ji,jj)
         !
         !                                      !* slopes referenced to geopot surfaces
         zsi_g_raw = zsi_raw - zsi_coord
         zsj_g_raw = zsj_raw - zsj_coord
         !                                      !* slope limiter
         !                                            
         ze3_e1    = e3u(ji,jj,iku,Kbb) * r1_e1u(ji,jj)     ! limit required in bilap case and usefull in lap case
         ze3_e2    = e3v(ji,jj,ikv,Kbb) * r1_e2v(ji,jj)
         !                                                  ! NB: hard coded factor 5 (may be a namelist parameter...)
         zsi_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e1, ABS( zsi_g_raw ) ), zsi_g_raw )
         zsj_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e2, ABS( zsj_g_raw ) ), zsj_g_raw )
         !
         !                                      !* slope / MLD_z0   MLD reference to z=0 (and based of ISF)
#if defined key_isf
         !                                            ! ISF : depth of u/v points referenced to the base of iceshelves or z=0
         zdepu = 0.5_wp * ( (gdepw_z0(ji,jj,iku,Kbb) + gdepw_z0(ji+1,jj,iku,Kbb)) - MAX( risfdep(ji,jj),risfdep(ji+1,jj) ) * 2._wp )
         zdepv = 0.5_wp * ( (gdepw_z0(ji,jj,ikv,Kbb) + gdepw_z0(ji,jj+1,ikv,Kbb)) - MAX( risfdep(ji,jj),risfdep(ji,jj+1) ) * 2._wp )  
#else
         !                                            ! no ISF: depth of u/v points referenced to z=0
         zdepu = 0.5_wp * (  gdepw_z0(ji,jj,iku,Kbb) + gdepw_z0(ji+1,jj,iku,Kbb)  )
         zdepv = 0.5_wp * (  gdepw_z0(ji,jj,ikv,Kbb) + gdepw_z0(ji,jj+1,ikv,Kbb)  )
#endif
         uslp_hml(ji,jj) = zsi_g_lim / zdepu * umask(ji,jj,iku)   ! mask needed as zsi_coord unmasked
         vslp_hml(ji,jj) = zsj_g_lim / zdepv * vmask(ji,jj,ikv)
      END_2D
      !
      !
      DO_2D( 0, 0, 0, 0 )         !==  i- & j-slopes/MLD at w-points  ==!   NB: slope just below the Mixed Layer : jk=mi,j[nmln]+1
         !
         ik   = MIN( nmln(ji,jj) + 1, jpk )
         ikm1 = MAX( 1, ik-1 )
         !                                !* Local vertical density gradient evaluated from N^2
         zzdzR = - MAX( zeps , zrho0_g * pn2(ji,jj,ik) )
         !
         !                                !* i- & j-slopes at w point
         zci = MAX(  umask(ji-1,jj,ik  ) + umask(ji,jj,ik  )           &
            &      + umask(ji-1,jj,ikm1) + umask(ji,jj,ikm1) , zeps  ) * zzdzR * e1t(ji,jj)
         zcj = MAX(  vmask(ji,jj-1,ik  ) + vmask(ji,jj,ikm1)           &
            &      + vmask(ji,jj-1,ikm1) + vmask(ji,jj,ik  ) , zeps  ) * zzdzR * e2t(ji,jj)
         !
         zsi_raw = (  p_gru(ji-1,jj,ik  ) +  p_gru(ji,jj,ikm1)           &
            &       + p_gru(ji-1,jj,ikm1) +  p_gru(ji,jj,ik  )   ) / zci * wmask (ji,jj,ik)
         zsj_raw = (  p_grv(ji,jj-1,ik  ) +  p_grv(ji,jj,ikm1)           &
            &       + p_grv(ji,jj-1,ikm1) +  p_grv(ji,jj,ik  )   ) / zcj * wmask (ji,jj,ik)
         !
         !                                !* coordinate slopes reference to z=0 at w-point
!!gm ISF case ??  think about it !  probably NOT OK  
         zsi_coord = 0.5_wp * ( gdepw_z0(ji+1,jj  ,ik,Kbb) - gdepw_z0(ji-1,jj  ,ik,Kbb) ) * r1_e1t(ji,jj)
         zsj_coord = 0.5_wp * ( gdepw_z0(ji  ,jj+1,ik,Kbb) - gdepw_z0(ji  ,jj-1,ik,Kbb) ) * r1_e2t(ji,jj)
         !
         !                                !* slopes referenced to geopot surfaces
         zsi_g_raw = zsi_raw - zsi_coord
         zsj_g_raw = zsj_raw - zsj_coord
         !                                      ! slope limiter
         !                                            ! limit required in bilap case and usefull in lap case
         ze3_e1    = e3w(ji,jj,ik,Kbb) * r1_e1u(ji,jj)
         ze3_e2    = e3w(ji,jj,ik,Kbb) * r1_e2v(ji,jj)
         !                                            ! NB: hard coded factor 5 (may be a namelist parameter...)
         zsi_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e1, ABS( zsi_g_raw ) ), zsi_g_raw )
         zsj_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e2, ABS( zsj_g_raw ) ), zsj_g_raw )
         !
         !                                !* i- & j-slope / MLD_z0      MLD reference to z=0 (and based of ISF)
#if defined key_ISF
         zhml_w =  MAX( 10._wp , hmlp(ji,jj)-risfdep(ji,jj)-ssh(ji,jj,Kbb) )
#else
         zhml_w =  MAX( 10._wp , hmlp(ji,jj)-ssh(ji,jj,Kbb) )
#endif
         wslpi_hml(ji,jj) = zsi_g_lim / zhml_w * umask(ji,jj,iku)
         wslpj_hml(ji,jj) = zsi_g_lim / zhml_w * vmask(ji,jj,ikv)
         !
      END_2D




      ! =================================
      ! III.  Eddy Induced streamfunction at uw/vw points     (computed from slopes referenced to z=0  WITHOUT horizontal filter)
      ! =================================----------------     | eipsi_uw = mk[zdxR] / mi[d/dzR]
      !                                                       ! eipsi_vw = mk[zdyR] / mj[d/dzR]
      IF( ln_ldfeiv ) THEN
         !
         DO_2D( 0, 0, 0, 0 )
            !
            !
            ik   = MIN( nmln(ji,jj) + 1, jpk )
            ikm1 = MAX( 1, ik-1 )
            !                                  !* Local vertical density gradient evaluated from N^2
            zzdzRi    = - MAX( zeps , zrho0_g * ( pn2(ji,jj,ik) + pn2(ji+1,jj,ik) ) )
            zzdzRj    = - MAX( zeps , zrho0_g * ( pn2(ji,jj,ik) + pn2(ji,jj+1,ik) ) )
            !
            !                                  !* local iso-level density gradient
            zzdiRw    =  ( p_gru(ji+1,jj,ik) + p_gru(ji,jj,ik) ) * r1_e1u(ji,jj) * wumask(ji,jj,ik)
            zzdjRw    =  ( p_grv(ji,jj+1,ik) + p_grv(ji,jj,ik) ) * r1_e2v(ji,jj) * wvmask(ji,jj,ik)
            !
            !                                  !* raw slopes
            zsi_raw   = zzdiRw / zzdzRi
            zsj_raw   = zzdjRw / zzdzRj
            !
            !                                  !* coordinate slopes at w-level reference to z=0
            zsi_coord = ( gdepw_z0(ji+1,jj,ik,Kbb) - gdepw_z0(ji,jj,ik,Kbb) ) * r1_e1u(ji,jj)
            zsj_coord = ( gdepw_z0(ji,jj+1,ik,Kbb) - gdepw_z0(ji,jj,ik,Kbb) ) * r1_e2v(ji,jj)
            !
            !                                  !* slopes at uw/vw-points referenced to geopot surfaces
            zsi_g_raw = zsi_raw - zsi_coord
            zsj_g_raw = zsj_raw - zsj_coord
            !
            !                                  !* slope limiter
            !                                         ! limit required in bilap case and usefull in lap case
            ze3_e1    = e3w(ji,jj,ik,Kbb) * r1_e1u(ji,jj)
            ze3_e2    = e3w(ji,jj,ik,Kbb) * r1_e2v(ji,jj)
            !                                         ! NB: hard coded factor 5 (may be a namelist parameter...)
            zsi_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e1, ABS( zsi_g_raw ) ), zsi_g_raw )
            zsj_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e2, ABS( zsj_g_raw ) ), zsj_g_raw )
            !                                         !
            !                                         ! Mixed Layer slope linear flattening
            zfi = MAX( omlmask(ji+1,jj  ,ik-1), omlmask(ji,jj,ikm1) )   ! zfi/j =1 in the ML
            zfj = MAX( omlmask(ji  ,jj+1,ik-1), omlmask(ji,jj,ikm1) )   !       =0 at the ML base and below
            !
            !                                            ! thickness of water column between z=0 and level k at uw/vw point
#if defined key_isf
            !                                            ! ISF : depth of u/v points referenced to the base of iceshelves or z=0
            zdepuw = 0.5_wp * (  (   gdepw_z0(ji,jj,ik,Kbb) + gdepw_z0(ji+1,jj,ik,Kbb)  )  &
               &               - MAX( risfdep(ji,jj)        ,  risfdep(ji+1,jj) ) * 2._wp  )
            zdepvw = 0.5_wp * (  (   gdepw_z0(ji,jj,ik,Kbb) + gdepw_z0(ji,jj+1,ik,Kbb) )   &
               &               - MAX( risfdep(ji,jj)        ,  risfdep(ji,jj+1) ) * 2._wp  )   
#else
            !                                            ! no ISF: depth of u/v points referenced to z=0
            zdepuw = 0.5_wp * ( gdepw_z0(ji,jj,ikm1,Kbb) + gdepw_z0(ji+1,jj  ,ikm1,Kbb) )
            zdepvw = 0.5_wp * ( gdepw_z0(ji,jj,ikm1,Kbb) + gdepw_z0(ji  ,jj+1,ikm1,Kbb) )
#endif
            !                                !* uw- & vw-slope / MLD_z0      MLD reference to z=0 (and based of ISF  ! NOT YET)
            uwslp_hml(ji,jj) = zsi_g_lim / zdepuw * wumask(ji,jj,ik)
            vwslp_hml(ji,jj) = zsj_g_lim / zdepvw * wvmask(ji,jj,ik)
            !
         END_2D
         !
         CALL lbc_lnk( 'ldfslp', uslp_hml , 'U', -1.0_wp, vslp_hml , 'V', -1.0_wp,   &
            &                    wslpi_hml, 'W', -1.0_wp, wslpj_hml, 'W', -1.0_wp,   &
            &                    uwslp_hml, 'U', -1.0_wp, vwslp_hml, 'V', -1.0_wp )
         !
      ELSE
      !
      !!gm this lbc_lnk should be useless....
         CALL lbc_lnk( 'ldfslp', uslp_hml , 'U', -1.0_wp, vslp_hml , 'V', -1.0_wp,   &
            &                    wslpi_hml, 'W', -1.0_wp, wslpj_hml, 'W', -1.0_wp )
      !
      ENDIF
   END SUBROUTINE ldf_slp_mxl


   SUBROUTINE ldf_slp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_slp_init  ***
      !!
      !! ** Purpose :   Initialization for the isopycnal slopes computation
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      INTEGER ::   ierr         ! local integer
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_slp_init : direction of lateral mixing'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      ALLOCATE( ah_wslp2(jpi,jpj,jpk) , akz(jpi,jpj,jpk) , STAT=ierr )
      IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'ldf_slp_init : unable to allocate ah_slp2 or akz' )
      !
      IF( ln_traldf_triad ) THEN        ! Griffies operator : triad of slopes
         IF(lwp) WRITE(numout,*) '   ==>>>   triad) operator (Griffies)'
         ALLOCATE( triadi_g(jpi,jpj,jpk,0:1,0:1) , triadj_g(jpi,jpj,jpk,0:1,0:1) ,     &
            &      triadi  (jpi,jpj,jpk,0:1,0:1) , triadj  (jpi,jpj,jpk,0:1,0:1) ,     &
            &      wslp2   (jpi,jpj,jpk)                                         , STAT=ierr )
         IF( ierr > 0      )   CALL ctl_stop( 'STOP', 'ldf_slp_init : unable to allocate Griffies operator slope' )
         IF( ln_dynldf_iso )   CALL ctl_stop( 'ldf_slp_init: Griffies operator on momentum not supported' )
         !
      ELSE                             ! Madec operator : slopes at u-, v-, and w-points
         IF(lwp) WRITE(numout,*) '   ==>>>   iso operator (Madec)'
         ALLOCATE( omlmask(jpi,jpj,jpk) ,                                                                        &
            &      uslp(jpi,jpj,jpk) , uslp_hml(jpi,jpj) , wslpi(jpi,jpj,jpk) , wslpi_hml(jpi,jpj) ,     &
            &      vslp(jpi,jpj,jpk) , vslp_hml(jpi,jpj) , wslpj(jpi,jpj,jpk) , wslpj_hml(jpi,jpj) , STAT=ierr )
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'ldf_slp_init : unable to allocate Madec operator slope ' )

         ! Direction of lateral diffusion (tracers and/or momentum)
         ! ------------------------------
         uslp (:,:,:) = 0._wp   ;   uslp_hml (:,:) = 0._wp      ! set the slope to zero (even in s-coordinates)
         vslp (:,:,:) = 0._wp   ;   vslp_hml (:,:) = 0._wp
         wslpi(:,:,:) = 0._wp   ;   wslpi_hml(:,:) = 0._wp
         wslpj(:,:,:) = 0._wp   ;   wslpj_hml(:,:) = 0._wp
         !
         IF( ln_ldfeiv ) THEN
            ALLOCATE( eipsi_uw(jpi,jpj,jpk), uwslp_hml(jpi,jpj) , eipsi_vw(jpi,jpj,jpk) , vwslp_hml(jpi,jpj) )
            eipsi_uw(:,:,:) = 0._wp   ;   uwslp_hml(:,:) = 0._wp   
            eipsi_vw(:,:,:) = 0._wp   ;   vwslp_hml(:,:) = 0._wp
         ENDIF
         !
         !!gm I no longer understand this.....
!!gm         IF( (ln_traldf_hor .OR. ln_dynldf_hor) .AND. .NOT. (.NOT.ln_linssh .AND. ln_rstart) ) THEN
!            IF(lwp)   WRITE(numout,*) '          Horizontal mixing in s-coordinate: slope = slope of s-surfaces'
!
!            ! geopotential diffusion in s-coordinates on tracers and/or momentum
!            ! The slopes of s-surfaces are computed once (no call to ldfslp in step)
!            ! The slopes for momentum diffusion are i- or j- averaged of those on tracers
!
!            ! set the slope of diffusion to the slope of s-surfaces
!            !      ( c a u t i o n : minus sign as dep has positive value )
!            DO jk = 1, jpk
!               DO jj = 2, jpjm1
!                  DO ji = 2, jpim1   ! vector opt.
!                     uslp (ji,jj,jk) = - ( gdept(ji+1,jj,jk,Kmm) - gdept(ji ,jj ,jk,Kmm) ) * r1_e1u(ji,jj) * umask(ji,jj,jk)
!                     vslp (ji,jj,jk) = - ( gdept(ji,jj+1,jk,Kmm) - gdept(ji ,jj ,jk,Kmm) ) * r1_e2v(ji,jj) * vmask(ji,jj,jk)
!                     wslpi(ji,jj,jk) = - ( gdepw(ji+1,jj,jk,Kmm) - gdepw(ji-1,jj,jk,Kmm) ) * r1_e1t(ji,jj) * wmask(ji,jj,jk) * 0.5
!                     wslpj(ji,jj,jk) = - ( gdepw(ji,jj+1,jk,Kmm) - gdepw(ji,jj-1,jk,Kmm) ) * r1_e2t(ji,jj) * wmask(ji,jj,jk) * 0.5
!                  END DO
!               END DO
!            END DO
!            CALL lbc_lnk( 'ldfslp', uslp , 'U', -1. ; CALL lbc_lnk( 'ldfslp', vslp , 'V', -1.,  wslpi, 'W', -1.,  wslpj, 'W', -1. )
!!gm         ENDIF
      ENDIF
      !
   END SUBROUTINE ldf_slp_init

   !!======================================================================
END MODULE ldfslp
