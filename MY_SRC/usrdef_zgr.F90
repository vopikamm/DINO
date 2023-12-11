MODULE usrdef_zgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_zgr  ***
   !!
   !!                      ===  BASIN configuration  ===
   !!
   !! User defined : vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2017-11  (J. Chanut)  Original code
   !!                 ! 2019-05  (R.Caneill and G.Madec) Adaptation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system
   !!      zgr_z      : reference 1D z-coordinate 
   !!      zgr_top_bot: ocean top and bottom level indices
   !!      zgr_zco    : 3D vertical coordinate in pure z-coordinate case
   !!---------------------------------------------------------------------
   USE oce    !, ONLY:            ! ocean variables
   USE dom_oce!, ONLY:        ! ocean domain
   USE phycst         ! physical constants
   !
   USE usrdef_nam, ONLY: nn_botcase, rn_distlam, rn_H, rn_hborder,            &
                     &   ln_sco_nam, ln_zco_nam, ln_zps_nam, nn_ztype,        &
                     &   ln_Iperio, rn_cha_min, rn_cha_max, rn_e1_deg,        &
                     &   nn_jeq_s, merc_proj, ln_mid_ridge, ln_drake_sill,  &
                     &   rn_dzmin, rn_hco, rn_kth, rn_acr, rn_slp_cha,        &
                     &   rn_ds_depth, rn_ds_width, rn_mr_depth, rn_mr_width,  &
                     &   rn_mr_lat_n, rn_mr_lat_s, nn_mr_edge, rn_phi_min     
   !
   !!rc USE depth_e3       ! depth <=> e3
   USE zgr_lib    ! tools for vertical coordinate
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library !!rc
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr        ! called by domzgr.F90

  !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS             

   SUBROUTINE usr_def_zgr( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
      &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,    &   ! 1D reference vertical coordinate
      &                    pdept , pdepw ,                             &   ! 3D t & w-points depth
      &                    pe3t  , pe3u  , pe3v   , pe3f ,             &   ! vertical scale factors
      &                    pe3w  , pe3uw , pe3vw         ,             &   !     -      -      -
      &                    k_top , k_bot    )                              ! top & bottom ocean level
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(out) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags
      LOGICAL                   , INTENT(out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pdept, pdepw                ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3w , pe3uw, pe3vw         ! i-scale factors 
      INTEGER , DIMENSION(:,:)  , INTENT(out) ::   k_top, k_bot                ! first & last ocean level
      !
      INTEGER  ::   ji, jj                                     ! dummy loop index
      REAL(wp) ::   zdzmin, zkth, zh_co, zacr                  ! Local scalars for the coordinate stretching
      REAL(wp) ::   zHmax

      REAL(wp), DIMENSION(jpi,jpj)           ::   zbathy, z2d   ! bathymetry, 2D variable
      REAL(wp) ::   zmidPhi            ! latitude of the central channel
      INTEGER  ::   nn_cha_min         ! index of the southern edge of the channel
      INTEGER  ::   nn_cha_max         ! index of the northern edge of the channel
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : BASIN configuration'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      ! Vertical coordinate type
      ld_zco      = ln_zco_nam
      ld_zps      = ln_zps_nam
      ld_sco      = ln_sco_nam
      ld_isfcav   = .FALSE.   ! no ice-shelves
      !
      zHmax = rn_H   ! Maximum depth of the basin

      !nn_cha_loc = NINT( ABS( rn_phi_min - rn_cha_loc )/ rn_e1_deg )

      CALL zgr_bat( zHmax, zbathy )   ! User creation of bathymetry
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      !
      !   === z-coordinate   ===   !
      IF( ld_zco ) THEN
         !
         CALL zgr_z1D( zHmax, nn_ztype, pdept_1d, pdepw_1d )   ! Reference z-coordinate system
         !
         !                                                       ! z-coordinate (3D arrays) from the 1D z-coord.
         CALL zgr_zco( pdept_1d, pdepw_1d,                   &   ! <==>  1D reference vertical coordinate
            &          pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
            &          pdept   , pdepw   ,                   &   !       3D t & w-points depth
            &          pe3t    , pe3u    , pe3v    , pe3f,   &   !       vertical scale factors
            &          pe3w    , pe3uw   , pe3vw           )     !           -      -      -
         !
         CALL zgr_msk_top_bot( pdept_1d, zbathy, k_top, k_bot )                 ! masked top and bottom ocean t-level indices
         !
      ENDIF
      !
      !   === z-coordinate with partial steps   ==   !
      !IF ( ld_zps ) THEN      !==  zps-coordinate  ==!   (partial bottom-steps)
      !   !
      !   CALL zgr_z1D( zHmax, nn_ztype, pdept_1d, pdepw_1d )   ! Reference z-coordinate system
      !   !
      !   ze3min = 0.1_wp * rn_dz
      !   IF(lwp) WRITE(numout,*) '   minimum thickness of the partial cells = 10 % of e3 = ', ze3min
      !   !
      !   !
      !   !                                !* bottom ocean compute from the depth of grid-points
      !   k_bot(:,:) = jpkm1
      !   DO jk = jpkm1, 1, -1
      !      WHERE( zht(:,:) < pdepw_1d(jk) + ze3min )   k_bot(:,:) = jk-1
      !   END DO
      !   !
      !   !                                !* vertical coordinate system
      !   DO jk = 1, jpk                      ! initialization to the reference z-coordinate
      !      pdept(:,:,jk) = pdept_1d(jk)
      !      pdepw(:,:,jk) = pdepw_1d(jk)
      !      pe3t (:,:,jk) = pe3t_1d (jk)
      !      pe3u (:,:,jk) = pe3t_1d (jk)
      !      pe3v (:,:,jk) = pe3t_1d (jk)
      !      pe3f (:,:,jk) = pe3t_1d (jk)
      !      pe3w (:,:,jk) = pe3w_1d (jk)
      !      pe3uw(:,:,jk) = pe3w_1d (jk)
      !      pe3vw(:,:,jk) = pe3w_1d (jk)
      !   END DO
      !   DO jj = 1, jpj                      ! bottom scale factors and depth at T- and W-points
      !      DO ji = 1, jpi
      !         ik = k_bot(ji,jj)
      !            pdepw(ji,jj,ik+1) = MIN( zht(ji,jj) , pdepw_1d(ik+1) )
      !            pe3t (ji,jj,ik  ) = pdepw(ji,jj,ik+1) - pdepw(ji,jj,ik)
      !            pe3t (ji,jj,ik+1) = pe3t (ji,jj,ik  ) 
      !            !
      !            pdept(ji,jj,ik  ) = pdepw(ji,jj,ik  ) + pe3t (ji,jj,ik  ) * 0.5_wp
      !            pdept(ji,jj,ik+1) = pdepw(ji,jj,ik+1) + pe3t (ji,jj,ik+1) * 0.5_wp
      !            pe3w (ji,jj,ik+1) = pdept(ji,jj,ik+1) - pdept(ji,jj,ik)              ! = pe3t (ji,jj,ik  )
      !      END DO
      !   END DO         
      !   !                                   ! bottom scale factors and depth at  U-, V-, UW and VW-points
      !   !                                   ! usually Computed as the minimum of neighbooring scale factors
      !   !pe3u (:,:,:) = pe3t(:,:,:)          ! HERE OVERFLOW configuration : 
      !   !pe3v (:,:,:) = pe3t(:,:,:)          !    e3 increases with i-index and identical with j-index
      !   !pe3f (:,:,:) = pe3t(:,:,:)          !    so e3 minimum of (i,i+1) points is (i) point
      !   !pe3uw(:,:,:) = pe3w(:,:,:)          !    in j-direction e3v=e3t and e3f=e3v
      !   !pe3vw(:,:,:) = pe3w(:,:,:)          !    ==>>  no need of lbc_lnk calls
      !   
      !   DO jk = 1,jpk                         ! Computed as the minimum of neighbooring scale factors
      !      DO jj = 1, jpj - 1
      !         DO ji = 1, jpi - 1   ! vector opt.
      !            pe3u (ji,jj,jk) = MIN( pe3t(ji,jj,jk), pe3t(ji+1,jj,jk) )
      !            pe3v (ji,jj,jk) = MIN( pe3t(ji,jj,jk), pe3t(ji,jj+1,jk) )
      !            pe3f (ji,jj,jk) = MIN( pe3t(ji,jj,jk), pe3t(ji+1,jj,jk), pe3t(ji,jj+1,jk), pe3t(ji+1,jj+1,jk) )
      !            pe3uw(ji,jj,jk) = MIN( pe3w(ji,jj,jk), pe3w(ji+1,jj,jk) )
      !            pe3vw(ji,jj,jk) = MIN( pe3w(ji,jj,jk), pe3w(ji,jj+1,jk) )
      !         END DO
      !      END DO
      !   END DO
      !   !
      !   CALL zgr_msk_top_bot( pdept_1d, zbathy, k_top, k_bot )
      !   !
      !ENDIF               ! masked top and bottom ocean t-level indices
   !
!   ===   s-coordinate   ===   !
      IF( ld_sco ) THEN
         !
         zdzmin  = rn_dzmin
         zkth    = rn_kth
         zh_co   = rn_hco
         zacr    = rn_acr
         !
         CALL zgr_sco_mi96( pdept_1d, pdepw_1d,                   &   ! ==>>  1D reference vertical coordinate
            &               zbathy  , zHmax   ,                   &   ! <<==  bathymetry
            &               zdzmin  , zkth    , zh_co , zacr,     &   ! <<==  parameters for the tanh stretching
            &               pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
            &               pdept   , pdepw   ,                   &   !       3D t & w-points depth
            &               pe3t    , pe3u    , pe3v    , pe3f,   &   !       vertical scale factors
            &               pe3w    , pe3uw   , pe3vw           )     !           -      -      -
         !
         ! Slope tests
         CALL zgr_test_slopes( zbathy, pdept, pe3u, pe3v )
         !
         !                                           ! masked top and bottom ocean t-level indices
         k_top(:,:) = 1    
         k_bot(:,:) = jpkm1
         !
      ENDIF
   
      zbathy(:,:) = 1._wp                         ! Use zbathy as land-sea mask to compute k_top and k_bot
      !
      IF (ln_Iperio) THEN
         nn_cha_min  = nn_jeq_s - merc_proj(rn_cha_min, rn_e1_deg) + 1
         nn_cha_max  = nn_jeq_s - merc_proj(rn_cha_max, rn_e1_deg) + 1
         DO jj = 1, jpj
            IF (mjg(jj) <= nn_cha_min) THEN
               zbathy(  mi0(     1+nn_hls):mi1(     1+nn_hls),jj) = 0._wp   ! first column of inner global domain at 0
               zbathy(  mi0(jpiglo-nn_hls):mi1(jpiglo-nn_hls),jj) = 0._wp   ! last  column of inner global domain at 0
            ENDIF
            IF (mjg(jj) >= nn_cha_max) THEN
               zbathy(  mi0(     1+nn_hls):mi1(     1+nn_hls),jj) = 0._wp   ! first column of inner global domain at 0
               zbathy(  mi0(jpiglo-nn_hls):mi1(jpiglo-nn_hls),jj) = 0._wp   ! last  column of inner global domain at 0
            ENDIF
         END DO
      ENDIF

      CALL lbc_lnk( 'usr_def_zgr', zbathy, 'T', 1._wp )

      k_top(:,:) = k_top(:,:) * INT( zbathy(:,:) )
      k_bot(:,:) = k_bot(:,:) * INT( zbathy(:,:) )
      !
   END SUBROUTINE usr_def_zgr

   
   SUBROUTINE zgr_bat( pH, pbathy )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_bat  ***
      !!
      !! ** Purpose :  Create analyticaly the bathymetry at T points, according to the nn_botcase option
      !!
      !! ** Method  :  Depending on nn_botcase will create different types of bathymetry.
      !!               The lon/lat coordinates are normalized between -1 and 1, and an analytical function is applied
      !!                  If nn_botcase == 0: flat bottom at a depth pH
      !!                  If nn_botcase == 1: bowl shape with a cosh
      !!                  If nn_botcase == 2: bowl shape with a (1-x**4)
      !!----------------------------------------------------------------------
      REAL(wp)                , INTENT(inout) ::   pH       ! ocean depth maximum   [m]
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pbathy   ! Ocean bathymetry      [m]
      !
      INTEGER  ::   ji, jj, jhl, inum             ! dummy loop indices
      REAL(wp) ::   zmaxlam, zmaxphi, zminlam, zminphi         ! local scalars (horizontal extent)
      REAL(wp) ::   zcha_min, zcha_max, zslp_cha               ! local scalars (channel)
      REAL(wp) ::   zmr_depth, zmr_width, zds_depth, zds_width ! local scalars (MR and DS)
      REAL(wp) ::   zmr_lat_s, zmr_lat_n, zcha_width, zgauss   ! local scalars (MR and DS)
      REAL(wp) ::   z1_H, z1_dLam, z1_dPhi, zxnorm, zynorm     ! local scalars
      REAL(wp) ::   zdistPhi, ztaper, zmidLam, zmidPhi, zrad   ! local scalars
      REAL(wp) ::   zdistLam, zwidth, zratio, zds_taper, zds   ! local scalars (boundary)
      REAL(wp) ::   zx, zy, zy_cha                            ! local normalised bathymetry
      REAL(wp), DIMENSION(jpi,jpj)     ::   zlamt, zphit       ! local horizontal coordinate arrays
      REAL(wp), DIMENSION(jpi,jpj)     ::   z2d
      !!----------------------------------------------------------------------
      !
      IF(lwp)   THEN
         WRITE(numout,*)
         WRITE(numout,*)    '    zgr_bat : defines the bathymetry at T points.'
         WRITE(numout,*)    '    ~~~~~~~'
         WRITE(numout,*)    '       BASIN case : flat bottom or bowl shape'
         WRITE(numout,*)    '       nn_botcase = ', nn_botcase
      ENDIF
      !!----------------------------------------------------------------------
      !Extend the horizontal coordinates beyond inner points for bathymetry computation
      zlamt(:,:) = glamt(:,:)
      zphit(:,:) = gphit(:,:)
      !
      IF( mig(jpi) == jpiglo ) THEN
         zlamt(jpi,:) = glamt(jpi-1,:) + 1./rad * ASIN ( TANH( rn_e1_deg * rad) )
      ENDIF
      !
      IF( mjg(jpj) == jpjglo ) THEN
         zphit(:,jpj) = gphit(:,jpj-1) + rn_e1_deg
      ENDIF
      !
      SELECT CASE(nn_botcase)
      !
      CASE(0)   ! flat bottom
         IF(lwp)   THEN
            WRITE(numout,*) '          i.e. flat bottom'
            WRITE(numout,*) '          Depth of the bathymetry   pH   = ', pH
         ENDIF
         pbathy(:,:) = pH
      CASE(1)   ! new bowl shape with exponentials
         IF(lwp)   THEN
            WRITE(numout,*) '          i.e. bowl shape in cosh(x)'
            WRITE(numout,*) '          see https://www.desmos.com/calculator/onczibqzb4'
            WRITE(numout,*) '          Asymptotical depth of the exponentials              pHasymp   = ', pH,      ' m'
            WRITE(numout,*) '          Depth of the bathymetry at the coast             rn_hborder   = ', rn_hborder, ' m'
            WRITE(numout,*) '          Typical length scale of the slope in longitude   rn_distLam   = ', rn_distLam,   ' degrees'
            WRITE(numout,*) '                                                                     '
            WRITE(numout,*) '             |    ^                                   ^              '
            WRITE(numout,*) '             | hborder                                |              '
            WRITE(numout,*) '             |    v                                   |              '
            WRITE(numout,*) '             |*                                       |              '
            WRITE(numout,*) '             | *                                     pH              '
            WRITE(numout,*) '             |  *                                     |              '
            WRITE(numout,*) '             |   *                                    |              '
            WRITE(numout,*) '             |    \**                                 |              '
            WRITE(numout,*) '             |     \ ***                              |              '
            WRITE(numout,*) '             |      \   *****                         v              '
            WRITE(numout,*) '             |       \       *********************************       '
            WRITE(numout,*) '              < zdist>                                               '
         ENDIF
         !
         CALL zgr_get_boundaries( zmaxlam, zminlam, zmaxphi, zminphi )
         zdistLam = rn_distLam
         zdistPhi =  COS( rad * zmaxphi ) * rn_distLam
         zwidth   = ( zmaxlam - zminlam )
         zslp_cha = rn_slp_cha
         zcha_min = rn_cha_min
         zcha_max = rn_cha_max
         zcha_width = ( zcha_max - zcha_min )
         !
         IF( ln_Iperio )   THEN
            ! Flattening the bathymetry in the channel
            If( zcha_min <= rn_phi_min ) THEN 
               ! Southern boundary is already southern boundary of the channel 
               zcha_min = zcha_max - zwidth     ! dk: TODO cleaner with separate routines left/right bathymetry 
            ENDIF
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                  zy_cha   = exp_bathy( zphit(ji,jj), zcha_min, zcha_max,  &
                     &            zwidth, zdistLam, zcha_width / 2 )                
                  zx       = exp_bathy( zlamt(ji,jj), zminlam, zmaxlam, zwidth,     &
                     &            zdistLam, zcha_width / 2 )                        &
                     &     * (1._wp - zy_cha) + zy_cha
                  zy       = exp_bathy( zphit(ji,jj), zminphi, zmaxphi, zwidth,     &
                     &            zdistPhi, zcha_width / 2 )
                  pbathy(ji,jj) = zx * zy * ( pH - rn_hborder ) + rn_hborder
            END_2D
         ELSE
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )               
                  zx       = exp_bathy( zlamt(ji,jj), zminlam, zmaxlam, zwidth,     &
                     &            zdistLam, zcha_width / 2 )                        
                  zy       = exp_bathy( zphit(ji,jj), zminphi, zmaxphi, zwidth,     &
                     &            zdistPhi, zcha_width / 2 )
                  pbathy(ji,jj) = zx * zy * ( pH - rn_hborder ) + rn_hborder
            END_2D
         ENDIF
         !
      CASE(2)   ! bowl shape 
         IF(lwp)   THEN
            WRITE(numout,*) '          i.e. bowl shape in (1-x**4)'
         ENDIF
         CALL zgr_get_boundaries( zmaxlam, zminlam, zmaxphi, zminphi )
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            ! Normalization between -1 and 1
            zx = 2._wp*(zlamt(ji,jj) - zminlam)/(zmaxlam - zminlam) - 1._wp
            zy = 2._wp*(zphit(ji,jj) - zminphi)/(zmaxphi - zminphi) - 1._wp
            pbathy(ji,jj) =  0.5_wp * pH * ( 1 + (1-zx**4)*(1-zy**4) )
         END_2D
      CASE DEFAULT
         CALL ctl_stop( 'zgr_bat: The chosen case for the bathymetry does not exist (nn_botcase).' )
      !
      END SELECT
      !
      ! Check that zHmax is really the maximum on the inner domain
      !
      z2d(:,:) = pbathy(:,:)
      CALL lbc_lnk( 'usr_def_zgr', z2d, 'T', 1._wp )   ! put the land boundaries to 0
      pH = MAXVAL(z2d(:,:))
      IF( lk_mpp )   CALL mpp_max( 'zgr_bat', pH )
      IF(lwp)   WRITE(numout,*) '      Effective Hmax   = ', pH, ' m      while rn_H   = ', rn_H, ' m'
      !
      ! Mid-Atlantic ridge
      !
      zmr_depth = rn_mr_depth
      zmr_width  = rn_mr_width
      zmr_lat_s  = rn_mr_lat_s
      zmr_lat_n  = rn_mr_lat_n
      !
      IF (ln_mid_ridge) THEN
         zmidLam = (zmaxlam + zminlam) / 2
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            IF (zphit(ji, jj) >= zmr_lat_s .AND. zphit(ji, jj) <= zmr_lat_n) THEN
               pbathy(ji,jj) = gauss_bathy( zmr_width, zmidLam, zlamt(ji, jj), zmr_depth, pbathy(ji, jj))
            ENDIF
         END_2D
         !
         SELECT CASE(nn_mr_edge)
         !
         CASE(0)     ! Edge in cosh shape
         !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               IF (zphit(ji, jj) >= zmr_lat_s - zcha_width / 2 .AND. zphit(ji, jj) <= rn_mr_lat_s) THEN
                  ztaper = cosh_bathy(rn_slp_cha, zmr_lat_s - zcha_width, zmr_lat_s, zphit(ji,jj))
                  pbathy(ji,jj) = pH * ztaper + pbathy(ji,jj) * ( 1._wp - ztaper )
               ELSE IF (zphit(ji, jj) >= zmr_lat_n .AND. zphit(ji, jj) <= rn_mr_lat_n + zcha_width / 2) THEN
                  ztaper = cosh_bathy(rn_slp_cha, zmr_lat_n, zmr_lat_n + zcha_width, zphit(ji,jj))
                  pbathy(ji,jj) = pH * ztaper + pbathy(ji,jj) * ( 1._wp - ztaper )
               ENDIF
         END_2D
         !
         CASE(1)     ! Edge in gaussian shape
         !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               IF (zphit(ji, jj) <= rn_mr_lat_s) THEN
                  zgauss = gauss_bathy(zmr_width, zmidLam, zlamt(ji, jj), zmr_depth, pH) !dk:TODO since ph, not pbathy maybe slight differences?
                  pbathy(ji,jj) = gauss_bathy( 2 * zslp_cha, rn_mr_lat_s, zphit(ji, jj), zgauss, pbathy(ji, jj))
               ELSE IF (zphit(ji, jj) >= rn_mr_lat_n) THEN
                  zgauss = gauss_bathy(zmr_width, zmidLam, zlamt(ji, jj), zmr_depth, pH)
                  pbathy(ji,jj) = gauss_bathy( 2 * zslp_cha, rn_mr_lat_n, zphit(ji, jj), zgauss, pbathy(ji, jj))
               ENDIF
         END_2D
         !
         CASE DEFAULT
            CALL ctl_stop( 'zgr_bat: The chosen case for the mid-atlantic ridge does not exist (nn_mr_edge).' )
         !
         END SELECT
         !
      ENDIF
      !
      ! Drake-Sill
      !
      zds_depth = rn_ds_depth
      zds_width = rn_ds_width
      !
      IF (ln_drake_sill) THEN
         zmidPhi = (rn_cha_max + rn_cha_min) / 2
         zrad    = (rn_cha_max - rn_cha_min) / 2
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            zds_taper  = smooth_step(zlamt(ji, jj), zminlam, zminlam + zds_width)
            zds        = gauss_ring( zds_width, zminlam, zmidphi, zrad, zlamt(ji, ji), zphit(ji, jj), zds_depth, pbathy(ji, jj) )
            pbathy(ji, jj) = zds_taper * zds + (1._wp - zds_taper) * pbathy(ji, jj)
         END_2D
      ENDIF
   END SUBROUTINE zgr_bat


   SUBROUTINE zgr_get_boundaries( pmaxlam, pminlam, pmaxphi, pminphi , cdpoint)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_get_boundaries  ***
      !!
      !! ** Purpose :   Compute the maximum and minimum values of latitude and longitude in the global inner domain
      !!
      !! ** Method  :   We take only the inner domain for the min/max of phi and lambda
      !!                   max phi is taken on gphiv( 2:jpim1 , 1:jpj - 1 ) so the last unmasked V point represents the boundary of the bathymetry
      !!                   max lam -    -   -  glamu( 1:jpim1 , 2:jpj - 1 ) -  -   last   -      U   -       -       -   -       -   -      -
      !!                   min lam -    -   -  glamu( 1:jpim1 , 2:jpj - 1 ) -  -  first   -      U   -       -       -   -       -   -      -
      !!                   min phi -    -   -  gphiv( 2:jpim1 , 1:jpj - 1 ) -  -  first   -      V   -       -       -   -       -   -      -
      !!                If jperio==8 (symmetrical at the equator), returns pminphi=-pmaxphi
      !!
      !!                If cdpoint is present (can be any character), we do the calculations on the inner T points of the domain
      !!----------------------------------------------------------------------
      CHARACTER, OPTIONAL, INTENT(in   ) ::   cdpoint
      REAL(wp)           , INTENT(  out) ::   pmaxlam, pminlam, pmaxphi, pminphi
      !
      REAL(wp), DIMENSION(4) :: zmaxmin
      !!----------------------------------------------------------------------
      IF( PRESENT(cdpoint) ) THEN   ! We take the min and max values at T points
         pmaxlam    = MAXVAL( glamt( 2:jpi-1 , 2:jpj-1 ) )
         pminlam    = MINVAL( glamt( 2:jpi-1 , 2:jpj-1 ) )
         pmaxphi    = MAXVAL( gphit( 2:jpi-1 , 2:jpj-1 ) )
         pminphi    = MINVAL( gphit( 2:jpi-1 , 2:jpj-1 ) )
      ELSE                        ! We take min and max values at the coast (i.e. U points for longitude, V points for latitude)
         pmaxlam    = MAXVAL( glamu( 1:jpi-1 , 2:jpj-1 ) )
         pminlam    = MINVAL( glamu( 1:jpi-1 , 2:jpj-1 ) )
         pmaxphi    = MAXVAL( gphiv( 2:jpi-1 , 1:jpj-1 ) )
         pminphi    = MINVAL( gphiv( 2:jpi-1 , 1:jpj-1 ) )
         !
      ENDIF
      IF( lk_mpp ) THEN
         ! max and min over the global domain
         ! the array is just used for the mpp communications (to reduce the time cost)
         !
         ! We use here the propertie for the min: MIN(a,b) = -MAX(-a,-b) if a,b >= 0
         ! we know that  90+minphi > 0,
         !              360+minlam > 0
         zmaxmin(1) = pmaxphi
         zmaxmin(2) = pmaxlam
         zmaxmin(3) = - (  90 + pminphi )
         zmaxmin(4) = - ( 360 + pminlam )
         CALL mpp_max( 'usrdef_zgr', zmaxmin )
         pmaxphi =   zmaxmin(1)
         pmaxlam =   zmaxmin(2)
         pminphi = - zmaxmin(3) -  90
         pminlam = - zmaxmin(4) - 360
         !
      ENDIF
      !
      IF( lwp ) WRITE(numout,*)
      IF( lwp ) WRITE(numout,*) '      zgr_get_boundaries : '
      IF( lwp ) WRITE(numout,*) '      ~~~~~~~~~~~~~~~~~~ '
      IF( lwp ) WRITE(numout,*) '         pminlam = ', pminlam
      IF( lwp ) WRITE(numout,*) '         pmaxlam = ', pmaxlam
      IF( lwp ) WRITE(numout,*) '         pminphi = ', pminphi
      IF( lwp ) WRITE(numout,*) '         pmaxphi = ', pmaxphi
   END SUBROUTINE zgr_get_boundaries

       
   SUBROUTINE zgr_z1D( pHmax, pztype, pdept_1d, pdepw_1d )   ! 1D reference vertical coordinate
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z  ***
      !!
      !! ** Purpose :   set the 1D depth of model levels and the resulting 
      !!              vertical scale factors.
      !!
      !! ** Method  :   1D z-coordinate system (use in all type of coordinate)
      !!       The depth of model levels is set from dep(k), an analytical function:
      !!                   w-level: depw_1d  = dep(k)
      !!                   t-level: dept_1d  = dep(k+0.5)
      !!       The scale factors are the discrete derivative of the depth:
      !!                   e3w_1d(jk) = dk[ dept_1d ] 
      !!                   e3t_1d(jk) = dk[ depw_1d ]
      !!           with at top and bottom :
      !!                   e3w_1d( 1 ) = 2 * ( dept_1d( 1 ) - depw_1d( 1 ) )
      !!                   e3t_1d(jpk) = 2 * ( dept_1d(jpk) - depw_1d(jpk) )
      !!       The depth are then re-computed from the sum of e3. This ensures 
      !!    that depths are identical when reading domain configuration file. 
      !!    Indeed, only e3. are saved in this file, depth are compute by a call
      !!    to the e3_to_depth subroutine.
      !!
      !!       Here the Madec & Imbard (1996) function is used.
      !!
      !! ** Action  : - pdept_1d, pdepw_1d : depth of T- and W-point (m)
      !!              - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!
      !! Reference : Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
      !!             Madec and Imbard, 1996, Clim. Dyn.
      !!----------------------------------------------------------------------
      REAL(wp)               , INTENT(in   ) ::   pHmax                ! ocean depth maximum   [m]
      INTEGER                , INTENT(in   ) ::   pztype               ! type of z grid (0: constant)
      REAL(wp), DIMENSION(:) , INTENT(  out) ::   pdept_1d, pdepw_1d   ! 1D grid-point depth   [m]
      !
      INTEGER  ::   jk       ! dummy loop indices
      REAL(wp) ::   zd       ! local scalar
      REAL(wp) ::   zt, zw   ! local scalars
      REAL(wp) ::   zsur, za0, za1, zkth, zacr   ! Values for the Madec & Imbard (1996) function  
      !!----------------------------------------------------------------------
      !
      zd = pHmax / REAL(jpkm1,wp)
      !
      IF(lwp) THEN            ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates '
         WRITE(numout,*) '    ~~~~~~~'
      ENDIF
      !
      ! 1D Reference z-coordinate    (using Madec & Imbard 1996 function)   !!rc TODO set non uniform spacing cf ORCA
      ! -------------------------
      !
      SELECT CASE( pztype )
      CASE( 0 )   ! Uniform vertical grid
         IF(lwp) THEN
            WRITE(numout,*) '       BASIN case : uniform vertical grid :'
            WRITE(numout,*) '                     with thickness = ', zd
         ENDIF
         pdepw_1d(1) = 0._wp
         pdept_1d(1) = 0.5_wp * zd
         ! 
         DO jk = 2, jpk          ! depth at T and W-points
            pdepw_1d(jk) = pdepw_1d(jk-1) + zd 
            pdept_1d(jk) = pdept_1d(jk-1) + zd 
         END DO
      CASE( 1 )   ! Non uniform spacing
         IF(lwp)   WRITE(numout,*) '       BASIN case : non uniform vertical grid'
         !CALL ctl_stop( 'zgr_z1D: The chosen case (pztype = 1) has not been implemeted yet' )
         ! taken from GYRE
         !!----------------------------------------------------------------------
         !
         ! Set parameters of z(k) function
         ! -------------------------------
         zsur = -2033.194295283385_wp       
         za0  =   155.8325369664153_wp 
         za1  =   146.3615918601890_wp
         zkth =    17.28520372419791_wp   
         zacr =     5.0_wp       
         !
         IF(lwp) THEN            ! Parameter print
            WRITE(numout,*)
            WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates '
            WRITE(numout,*) '    ~~~~~~~'
            WRITE(numout,*) '       GYRE case : MI96 function with the following coefficients :'
            WRITE(numout,*) '                 zsur = ', zsur
            WRITE(numout,*) '                 za0  = ', za0
            WRITE(numout,*) '                 za1  = ', za1
            WRITE(numout,*) '                 zkth = ', zkth
            WRITE(numout,*) '                 zacr = ', zacr
         ENDIF
   
         !
         ! 1D Reference z-coordinate    (using Madec & Imbard 1996 function)
         ! -------------------------
         !
         DO jk = 1, jpk          ! depth at T and W-points
            zw = REAL( jk , wp )
            zt = REAL( jk , wp ) + 0.5_wp
            pdepw_1d(jk) = ( zsur + za0 * zw + za1 * zacr *  LOG( COSH( (zw-zkth) / zacr ) )  )
            pdept_1d(jk) = ( zsur + za0 * zt + za1 * zacr *  LOG( COSH( (zt-zkth) / zacr ) )  )
         END DO
         !
   
      CASE DEFAULT
         CALL ctl_stop( 'zgr_z1D: The chosen case for the vertical 1D grid does not exist' )
      END SELECT
      !
   END SUBROUTINE zgr_z1D


   SUBROUTINE zgr_msk_top_bot( pdept_1d , pbathy, k_top , k_bot )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_msk_top_bot  ***
      !!
      !! ** Purpose :   set the masked top and bottom ocean t-levels for z-coordinates
      !!
      !! ** Method  :   BASIN case = closed or south symmetrical box ocean without ocean cavities
      !!                      k_top = 1                                       except along boundaries according to jperio
      !!                      k_bot = {first point under the ocean bottom}    except along boundaries according to jperio
      !!
      !! ** Action  : - k_top : first wet ocean level index
      !!              - k_bot : last  wet ocean level index
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)  , INTENT(in   ) ::   pdept_1d
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pbathy          ! bathymetry   [m]
      INTEGER , DIMENSION(:,:), INTENT(  out) ::   k_top , k_bot   ! first & last wet ocean level
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d      ! local 2D variable
      INTEGER :: jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_top_bot : defines the top and bottom wet ocean levels for z-coordinates.'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '       BASIN case : = closed or south symmetrical box ocean without ocean cavities'
      IF(lwp) WRITE(numout,*) '          k_top = 1                                       except along boundaries according to jperio'
      IF(lwp) WRITE(numout,*) '          k_bot = {first point under the ocean bottom}    except along boundaries according to jperio'
      !
      !
      ! bottom ocean mask computed from the depth of grid-points
      k_bot(:,:) = jpkm1   ! initialisation of k_bot to its maximum authorized value
      DO jk = 1, jpkm1
         WHERE( pdept_1d(jk) < pbathy(:,:) .AND. pbathy(:,:) <= pdept_1d(jk+1) )   z2d(:,:) = REAL(jk,wp)
      END DO
      !
      CALL lbc_lnk( 'usrdef_zgr', z2d, 'T', 1. )           ! set surrounding land point to zero (depending on jperio value)
      !
      k_bot(:,:) = INT( z2d(:,:) )
      k_top(:,:) = MIN( 1 , k_bot(:,:) )     ! = 1 over the ocean point, =0 elsewhere
      
      WRITE(numout,*) '    |     ', k_bot(5,2), '    |     ' , k_bot(5,3), '    |     ', k_bot(5,4), '    |     ' , k_bot(5,5), '    |     '
      WRITE(numout,*) '    |     ', k_bot(4,2), '    |     ' , k_bot(4,3), '    |     ', k_bot(4,4), '    |     ' , k_bot(4,5), '    |     '
      WRITE(numout,*) '    |     ', k_bot(3,2), '    |     ' , k_bot(3,3), '    |     ', k_bot(3,4), '    |     ' , k_bot(3,5), '    |     '
      WRITE(numout,*) '    |     ', k_bot(2,2), '    |     ' , k_bot(2,3), '    |     ', k_bot(2,4), '    |     ' , k_bot(2,5), '    |     '
      !
   END SUBROUTINE zgr_msk_top_bot


   SUBROUTINE zgr_test_slopes( pbathy, pdept, pe3u, pe3v )
!!rc TODO implement multi processors
!!rc TODO remove masked U and V points
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_test_slopes  ***
      !!
      !! ** Purpose :   TODO
      !!
      !! ** Method  :   TODO
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)  , INTENT(in) ::   pbathy
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   pdept       ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   pe3u , pe3v  ! vertical scale factors  [m]
      !
      INTEGER  ::   inum, ji, jj, jk, i_criteria_not_respected_i, i_criteria_not_respected_j
      REAL(WP) ::   zsi_bat, zsj_bat, zsmax_bat   ! vars for checking the slopes (    bathymetry       relative to the    horizontal          )
      REAL(WP), DIMENSION(2)                 ::   zsi_bat_loc, zsj_bat_loc      ! locs for checking the slopes (    bathymetry       relative to the    horizontal          )
      REAL(WP) ::   zsi_dia, zsj_dia, zsmin_dia   ! vars for checking the slopes (diagonal of the grid relative to the slope of the grid      )
      REAL(WP) ::   zsi_gri, zsj_gri, zsmax_gri   ! vars for checking the slopes (      grid           relative to a horizontal stratification)
      REAL(wp), DIMENSION(jpkm1)             ::   zsi1d_dia, zsj1d_dia, zsi1d_gri, zsj1d_gri
      REAL(wp), DIMENSION(jpi-2,jpj-2,jpkm1) ::   zglo_ri, zglo_rj    ! global ratio of slopes of the grid
      REAL(wp), DIMENSION(jpi,jpj)           ::   z2d   ! 2D variable
      !!----------------------------------------------------------------------
      ! maximum slope of bathymetry
      ! We take here only the inner values
      zsi_bat   = MAXVAL(        ABS( (pbathy(3:jpi   , 2:jpj - 1)      - pbathy(2:jpi - 1 , 2:jpj - 1))      * r1_e1u(2:jpi - 1 , 2:jpj - 1) ) )
      zsi_bat_loc   = MAXLOC(        ABS( (pbathy(3:jpi   , 2:jpj - 1)      - pbathy(2:jpi - 1 , 2:jpj - 1))      * r1_e1u(2:jpi - 1 , 2:jpj - 1) ) )

      zsj_bat   = MAXVAL(        ABS( (pbathy(2:jpi - 1 , 3:jpj  )      - pbathy(2:jpi - 1 , 2:jpj - 1))      * r1_e2v(2:jpi - 1 , 2:jpj - 1) ) )
      zsj_bat_loc   = MAXLOC(        ABS( (pbathy(2:jpi - 1 , 3:jpj  )      - pbathy(2:jpi - 1 , 2:jpj - 1))      * r1_e2v(2:jpi - 1 , 2:jpj - 1) ) )
      !
      DO jk = 1, jpkm1
         ! Maximum slope of the grid relative to a horizontal stratification 
         zsi1d_gri(jk) = MAXVAL( ABS( (pdept( 3:jpi   , 2:jpj - 1 , jk) - pdept( 2:jpi - 1 , 2:jpj - 1 , jk)) * r1_e1u(2:jpi - 1 , 2:jpj - 1) ) )
         zsj1d_gri(jk) = MAXVAL( ABS( (pdept( 2:jpi - 1 , 3:jpj   , jk) - pdept( 2:jpi - 1 , 2:jpj - 1 , jk)) * r1_e2v(2:jpi - 1 , 2:jpj - 1) ) )
         ! minimum slopes in the diagonals of the grid cells
         zsi1d_dia(jk) = MINVAL( pe3u(2:jpi - 1 , 2:jpj - 1 , jk) * r1_e1u(2:jpi - 1 , 2:jpj - 1) )
         zsj1d_dia(jk) = MINVAL( pe3v(2:jpi - 1 , 2:jpj - 1 , jk) * r1_e2v(2:jpi - 1 , 2:jpj - 1) )
         ! global ratio of slope (must be <= 1)
         !!rc WRITE(numout,*) 'usrdef_zgr DEBUG   jpi, jpj', jpi, jpj
         
         zglo_ri(:,:,jk) = ABS( (pdept( 3:jpi   , 2:jpj - 1 , jk) - pdept( 2:jpi - 1 , 2:jpj - 1 , jk)) * r1_e1u(2:jpi - 1 , 2:jpj - 1) ) * ( e1u(2:jpi - 1 , 2:jpj - 1) / pe3u(2:jpi - 1 , 2:jpj - 1 , jk))
         zglo_rj(:,:,jk) = ABS( (pdept( 2:jpi - 1 , 3:jpj   , jk) - pdept( 2:jpi - 1 , 2:jpj - 1 , jk)) * r1_e2v(2:jpi - 1 , 2:jpj - 1) ) * ( e2v(2:jpi - 1 , 2:jpj - 1) / pe3v(2:jpi - 1 , 2:jpj - 1 , jk))
         !
         !!rc zglo_ri(:,:,jk) = ABS( (pdept( 2:jpi   , 1:jpj - 1 , jk) - pdept( 1:jpi - 1 , 1:jpj - 1 , jk)) * r1_e1u(1:jpi - 1 , 1:jpj - 1) ) * ( e1u(1:jpi - 1 , 1:jpj - 1) / pe3u(1:jpi - 1 , 1:jpj - 1 , jk))
         !!rc zglo_rj(:,:,jk) = ABS( (pdept( 1:jpi - 1 , 2:jpj   , jk) - pdept( 1:jpi - 1 , 1:jpj - 1 , jk)) * r1_e2v(1:jpi - 1 , 1:jpj - 1) ) * ( e2v(1:jpi - 1 , 1:jpj - 1) / pe3v(1:jpi - 1 , 1:jpj - 1 , jk))
      END DO
      !
      zsmax_bat = MAX(          zsi_bat ,          zsj_bat  )
      zsmax_gri = MAX( MAXVAL(zsi1d_gri), MAXVAL(zsj1d_gri) )
      zsmin_dia = MIN( MINVAL(zsi1d_dia), MINVAL(zsj1d_dia) )
      !
      i_criteria_not_respected_i = SUM( MERGE( 1,0,SUM(MERGE(1,0,zglo_ri > 3), 3) >= 1 ) ) ! get 1 if zglo_ri > 3 somewhere in the column, 0 else
      i_criteria_not_respected_j = SUM( MERGE( 1,0,SUM(MERGE(1,0,zglo_rj > 3), 3) >= 1 ) )
      !
      ! Open and create a netcdf file
      CALL iom_open( 'slope_ratio', inum, ldwrt = .TRUE.)
      z2d(:,:) = 0._wp
      z2d(2:jpi - 1,2:jpj - 1) = SUM(MERGE(1,0,zglo_ri > 3), 3)
      CALL iom_rstput( 0, 0, inum, 'ratio_slopei', z2d)
      z2d(2:jpi - 1,2:jpj - 1) = SUM(MERGE(1,0,zglo_rj > 3), 3)
      CALL iom_rstput( 0, 0, inum, 'ratio_slopej', z2d)
      CALL iom_close( inum )
      !
      ! If mpp, we compute the min / max over the whole domain
      IF( lk_mpp ) THEN
         CALL mpp_max( 'usr_def_zgr', zsmax_bat )
         CALL mpp_max( 'usr_def_zgr', zsmax_gri )
         CALL mpp_min( 'usr_def_zgr', zsmin_dia )
      ENDIF
      !
      IF( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' usr_def_zgr : Diag for the slopes (used in lateral diffusion)'
         WRITE(numout,*) ' ~~~~~~~~~~~'
         WRITE(numout,*) '    Maximum slope of the bathymetry                        = ', zsmax_bat * 100, ' per cent'
         WRITE(numout,*) '    Location of maximum slope of the bathymetry along i    = ', zsi_bat_loc * 100, ' per cent'
         WRITE(numout,*) '    Location of maximum slope of the bathymetry along j    = ', zsj_bat_loc * 100, ' per cent'
         WRITE(numout,*) '    Maximum slope of the grid                              = ', zsmax_gri * 100, ' per cent'
         WRITE(numout,*) '    Minimum value of e3/e1 (i.e. slope of the diagonals)   = ', zsmin_dia * 100    , ' per cent'
         WRITE(numout,*) '    Corresponding maximum accepted slope                   = ', zsmin_dia * 100 * 3, ' per cent, with a factor 3'
         WRITE(numout,*) '    Corresponding maximum accepted slope                   = ', zsmin_dia * 100 * 5, ' per cent, with a factor 5'
         WRITE(numout,*) '    Maximum value of ratio (slope grid relative to horiz) / (slope diagonal relative to slope of the grid)'
         WRITE(numout,*) '                                                           = ', MAXVAL(zglo_ri), ' must be <= to 3 or 5 depending on the factor TODO use a var + mpp_max'
         WRITE(numout,*) '    Number of i point where the criteria is not respected  = ', i_criteria_not_respected_i
         WRITE(numout,*) '    Number of j point where the criteria is not respected  = ', i_criteria_not_respected_j
         WRITE(numout,*) '    Number tot of none masked horizontal points            = ', (jpj-2)*(jpi-2)
      ENDIF
      !
   END SUBROUTINE zgr_test_slopes

   FUNCTION cosh_bathy( pdistLam, pminLam, pmaxLam, pLam) RESULT( pcosh )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE  ***
      !!
      !! ** Purpose :   Compute a bowl-shaped bathymetry.
      !!
      !! ** Method  :   see https://www.desmos.com/calculator/onczibqzb4'
      !!
      !!                    <---------------------pwid---------------------->
      !!                   |    ^                   ^                        | 
      !!                   | hborder                |                        | 
      !!                   |    v                   |                        | 
      !!                   |*                       |                       *| 
      !!                   | *                     pH                      * | 
      !!                   |  *                     |                     *  | 
      !!                   |   *                    |                    *   | 
      !!                   |    \**                 |                 **/    | 
      !!                   |     \ ***              |              *** /     | 
      !!                   |      \   *****         v         *****   /      | 
      !!                   |       \       *******************       /       | 
      !!                    < zdist>
      !!     
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )             ::   pdistLam    ! Length scale of the slope               [degrees]
      REAL(wp), INTENT(in   )             ::   pmaxLam     ! Maximum longitude/latitude              [degrees]
      REAL(wp), INTENT(in   )             ::   pminLam     ! Minimum longitude/latitude              [degrees]
      REAL(wp), INTENT(in   )             ::   pLam        ! Value of longitude/latitude
      REAL(wp)                            ::   pcosh       ! Resulting tapering factor
      !!----------------------------------------------------------------------
      
      REAL(wp) ::   znorm     ! Computed parameters
      !!----------------------------------------------------------------------
      !
      znorm    = 1._wp + EXP( ( pminLam - pmaxLam ) / pdistLam )
      !
      pcosh   = 1._wp - 1._wp / znorm * (  EXP(   ( pLam - pmaxLam ) / pdistLam )   &
         &                                + EXP( - ( pLam - pminLam ) / pdistLam ) )
      !
   END FUNCTION cosh_bathy

   FUNCTION exp_bathy( plam,      plam_lft,      plam_rgt,   pwidth,    &
      &                  pdist_lam, pdist_taper                            ) RESULT( pexp )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE  ***
      !!
      !! ** Purpose :   Compute a bowl-shaped bathymetry tapered to reach the bottom.

      !!
      !! ** Method  :   see https://www.desmos.com/calculator/onczibqzb4'
      !!
      !!                    <---------------------pwid---------------------->
      !!                   |    ^                   ^                        | 
      !!                   | hborder                |                        | 
      !!                   |    v                   |                        | 
      !!                   |*                       |                       *| 
      !!                   | *                     pH                      * | 
      !!                   |  *                     |                     *  | 
      !!                   |   *                    |                    *   | 
      !!                   |    \**                 |                 **/    | 
      !!                   |     \ ***              |              *** /     | 
      !!                   |      \   *****         |         *****   /      | 
      !!                   |       \       *******  v    ******      /       |
      !!                   |<    ztaper_dist     >******            /
      !!                    < zdist> \                             /        
      !!     
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )    ::   plam         ! Value of longitude/latitude                      [degrees]
      REAL(wp), INTENT(in   )    ::   plam_lft     ! Location of the left boundary                    [degrees]
      REAL(wp), INTENT(in   )    ::   plam_rgt     ! Location of the right boundary                   [degrees]
      REAL(wp), INTENT(in   )    ::   pwidth       ! Width of the basin                               [degrees]
      REAL(wp), INTENT(in   )    ::   pdist_lam    ! Length scale of the slope                        [degrees]
      REAL(wp), INTENT(in   )    ::   pdist_taper  ! Distance where the bathymetry meets the bottom   [degrees]
      REAL(wp)                   ::   pexp         ! Resulting tapered bathymetry
      !!----------------------------------------------------------------------
      REAL(wp) ::   zstep                          ! Smooth step function for tampering  
      REAL(wp) ::   znorm                          ! Normalisation
      !!----------------------------------------------------------------------
      !
      IF(plam < plam_lft) THEN
         pexp = 0.0
      ELSE IF(plam >= plam_lft .AND. plam <= plam_lft + pdist_taper) THEN
         zstep = smooth_step(plam, plam_lft, plam_lft + pdist_taper)
         znorm    = 1 + EXP(-( pwidth) / pdist_lam) 
         pexp     = (1 - (EXP( - (plam - plam_lft) / pdist_lam)) / znorm )    &
            &     * (1 - zstep) + zstep
      ELSE IF(plam > plam_lft + pdist_taper .AND. plam < plam_rgt - pdist_taper ) THEN
         pexp = 1.0
      ELSE IF(plam >= plam_rgt - pdist_taper .AND. plam <= plam_rgt) THEN
         zstep    = 1 - smooth_step(plam, plam_rgt - pdist_taper, plam_rgt)
         znorm    = 1 + EXP(-( pwidth) / pdist_lam) 
         pexp     = (1 - (EXP( (plam - plam_rgt) / pdist_lam)) / znorm )    &
            &     * (1 - zstep) + zstep
      ELSE
         pexp = 0.0
      ENDIF
      !
   END FUNCTION exp_bathy

   FUNCTION smooth_step( plam, plam_lft, plam_rgt) RESULT( pstep )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE  ***
      !!
      !! ** Purpose :   Compute a smooth step function for tapering
      !! 
      !!             |                                ***|*****---> 1 
      !!             |                           *****   | 
      !!             |                       ****        | 
      !!             |                     **            | 
      !!             |                    *              | 
      !!             |                   *               | 
      !!             |                 **                | 
      !!             |              ***                  | 
      !!             |         *****                     | 
      !!             |   ******                          |
      !!  0 <---*****|***                                |
      !!             v                                   v
      !!          plam_lft                            plam_rgt
      !!     
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )    ::   plam         ! Value of longitude/latitude                      [degrees]
      REAL(wp), INTENT(in   )    ::   plam_lft     ! Location of the left boundary                    [degrees]
      REAL(wp), INTENT(in   )    ::   plam_rgt     ! Location of the right boundary                   [degrees]
      REAL(wp)                   ::   pstep        ! Smooth step function
      !!----------------------------------------------------------------------
      !
      IF(plam < plam_lft) THEN
         pstep = 0.0
      ELSE IF(plam >= plam_lft .AND. plam <= plam_rgt) THEN
         pstep =     6 * ( (plam - plam_lft) / (plam_rgt - plam_lft) ) ** 5    &
            &     - 15 * ( (plam - plam_lft) / (plam_rgt - plam_lft) ) ** 4    &
            &     + 10 * ( (plam - plam_lft) / (plam_rgt - plam_lft) ) ** 3
      ELSE
         pstep = 1.0
      ENDIF
      !
   END FUNCTION smooth_step



   FUNCTION gauss_bathy( pdistLam, pmidLam, pLam, pdep_top, pdep_bot) RESULT( pgauss )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE  ***
      !!
      !! ** Purpose :   Compute a gaussian bump.
      !!
      !! ** Method  :   see https://www.desmos.com/calculator/dgcewoh3bx'
      !!
      !!     
      !!----------------------------------------------------------------------
      REAL(wp)    , INTENT(in    )           ::   pdistLam     ! Length scale of the slope               [degrees]
      REAL(wp)    , INTENT(in    )           ::   pmidLam      ! longitude/latitude of the bump          [degrees]
      REAL(wp)    , INTENT(in    )           ::   pLam         ! Value of longitude/latitude             [degrees]
      REAL(wp)    , INTENT(in    )           ::   pdep_top     ! depth of the gaussian bump              [meters]
      REAL(wp)    , INTENT(in    )           ::   pdep_bot     ! depth of the basin                      [meters]
      REAL(wp)                               ::   pgauss       ! Resulting bathymetry depth              [meters]
      !!----------------------------------------------------------------------
      !
      IF (pdep_bot >= pdep_top) THEN
         pgauss   = (pdep_top - pdep_bot) * EXP( - ( ( pLam - pmidLam ) / pdistLam )**2 ) + pdep_bot
      ELSE
         pgauss   = pdep_bot
      ENDIF
      !
   END FUNCTION gauss_bathy

   FUNCTION gauss_ring( pdistLam, plam0, pphi0, prad, pLam, pPhi, pdep_top, pdep_bot) RESULT( pring )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE  ***
      !!
      !! ** Purpose :   Compute a ring of gaussian bumps.
      !!
      !! ** Method  :   
      !!
      !!     
      !!----------------------------------------------------------------------
      REAL(wp)    , INTENT(in    )           ::   pdistLam     ! Length scale of the slope               [degrees]
      REAL(wp)    , INTENT(in    )           ::   plam0        ! Longitude of the center                 [degrees]
      REAL(wp)    , INTENT(in    )           ::   pphi0        ! Latitude of the center                  [degrees]
      REAL(wp)    , INTENT(in    )           ::   prad         ! radius of the circle                    [degrees]
      REAL(wp)    , INTENT(in    )           ::   pLam         ! Longitude                               [degrees]
      REAL(wp)    , INTENT(in    )           ::   pPhi         ! Latitude                                [degrees]
      REAL(wp)    , INTENT(in    )           ::   pdep_top     ! depth of the gaussian bump              [meters]
      REAL(wp)    , INTENT(in    )           ::   pdep_bot     ! depth of the basin                      [meters]
      REAL(wp)                               ::   pring        ! Resulting bathymetry depth              [meters]
      !!----------------------------------------------------------------------
      !
      REAL(wp)                               ::   zx, zy, zexp   ! local scalars
      zx = pLam - plam0
      zy = pPhi - pphi0

      zexp = (- zx**2 - zy**2 + 2 * prad * SQRT( zx**2 + zy**2 ) - prad**2) / pdistLam**2
      IF (pdep_bot >= pdep_top) THEN
         pring   = (pdep_top - pdep_bot) * EXP( zexp ) + pdep_bot
      ELSE
         pring   = pdep_bot
      ENDIF
      !
   END FUNCTION gauss_ring
       
   !!======================================================================
END MODULE usrdef_zgr
