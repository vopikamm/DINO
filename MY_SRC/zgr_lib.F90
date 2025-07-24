MODULE zgr_lib
   !!==============================================================================
   !!                       ***  MODULE usrlib_zgr   ***
   !! Ocean domain : library routines used to build a vertical coordinate 
   !!==============================================================================
   !! History :  
   !!            4.0  ! 2019-05  ( R. Caneil, G. Madec)  Original code 
   !!            4.2  ! 2025-03  ( D. Kamm )             Adapt for DINO hybrid coordinates
   !....
   !!             -   ! 2005-10  (A. Beckmann)  modifications for hybrid s-ccordinates & new stretching function
   !!            3.0  ! 2008-06  (G. Madec)  insertion of domzgr_zps.h90 & coding style
   !!            3.3  ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!            3.4  ! 2012-08  (J. Siddorn) added Siddorn and Furner stretching function
   !!            3.6  ! 2014-11  (P. Mathiot and C. Harris) add ice shelf capabilitye  
   !!            3.?  ! 2015-11  (H. Liu) Modifications for Wetting/Drying
   !....
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zgr_bot_level: deepest ocean level for t-, u, and v-points
   !!   zgr_z        : reference z-coordinate 
   !!   zgr_zps      : z-coordinate with partial steps
   !!   zgr_sco      : s-coordinate
   !!   fssig        : tanh stretch function
   !!   fssig1       : Song and Haidvogel 1994 stretch function
   !!   fgamma       : Siddorn and Furner 2012 stretching function
   !!---------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zgr_zco      : z-coordinate
   !!
   !!   depth_to_e3    : use the depth of t- and w-points to calculate e3t & e3w
   !!                    (generic interface for 1D and 3D fields)
   !!   e3_to_depth    : use e3t & e3w to calculate the depth of t- and w-points
   !!                    (generic interface for 1D and 3D fields)
   !!  e3tw_to_other_e3: set e3u, e3v, e3f from e3t and e3uw, e3vw, e3w from e3w
   !!---------------------------------------------------------------------
   USE oce               ! ocean variables
   USE dom_oce           ! ocean domain
   !
   USE in_out_manager    ! I/O manager
   USE lbclnk            ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp           ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE
  
   INTERFACE depth_to_e3
      MODULE PROCEDURE depth_to_e3_1d, depth_to_e3_3d
   END INTERFACE

   INTERFACE e3_to_depth
      MODULE PROCEDURE e3_to_depth_1d, e3_to_depth_3d
   END INTERFACE
   PUBLIC   zgr_zco            ! called by usrdef_zgr
   PUBLIC   zgr_sco_mi96       ! called by usrdef_zgr
   PUBLIC   depth_to_e3        ! called by usrdef_zgr
   PUBLIC   e3_to_depth        ! called by usrdef_zgr and domzgr.F90
      
   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: depth_e3.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS             

   SUBROUTINE zgr_zco( pdept_1d, pdepw_1d,                   &   ! <==>  1D reference vertical coordinate
      &                pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
      &                pdept   , pdepw   ,                   &   ! ==>>  3D t & w-points depth
      &                pe3t    , pe3u    , pe3v    , pe3f,   &   ! ==>>  vertical scale factors
      &                pe3w    , pe3uw   , pe3vw           )     ! ==>>      -      -      -
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE zgr_zco  ***
      !!
      !! ** Purpose :   User defined the full step vertical coordinate  (zco)
      !!              from a 1D t- and w-points depth
      !!                recompute the depth as the sum of scale factor in 
      !!              order to preserve the truncation between  pdepw and ht_0
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(inout) ::   pdept_1d, pdepw_1d          ! 1D t- and w-depth          [m]
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pe3t_1d, pe3w_1d            ! 1D t and w-scale factors   [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! 3D t and w-depth           [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         ! at t, u, v, f, w, uw, vw points
      !
      INTEGER  ::   jk   ! do-loop index
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'zgr_lib: zgr_zco : define full step vertical coord. system from 1D z-coordinate'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~'
      !
      DO jk = 1, jpk          !==  3D depth = 1D depth  ==!
         pdept(:,:,jk) = pdept_1d(jk) 
         pdepw(:,:,jk) = pdepw_1d(jk)
      END DO
      !                       !==  t- and w- scale factors from depth  ==!
      CALL depth_to_e3( pdept   , pdepw   , pe3t   , pe3w    )
      CALL depth_to_e3( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d )
      !
      !                       !==  update t- and w-depths from SUM(e3)  <== needed to get the same last digit in ht_0 calculation
      CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d ) 
      CALL e3_to_depth( pe3t   ,  pe3w  , pdept   , pdepw    ) 
      !
      !                       !==  other scale factors  ==!
      pe3u (:,:,:) = pe3t(:,:,:)          ! t-level 
      pe3v (:,:,:) = pe3t(:,:,:)
      pe3f (:,:,:) = pe3t(:,:,:)
      pe3uw(:,:,:) = pe3w(:,:,:)          ! w-level 
      pe3vw(:,:,:) = pe3w(:,:,:)
      !
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level   pdept_1d   pdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
   END SUBROUTINE zgr_zco

   SUBROUTINE zgr_sco_mi96( pdept_1d, pdepw_1d,                   &   ! ==>>  1D reference vertical coordinate
      &                     pbathy  , pHmax   ,                   &   ! <<==  bathymetry
      &                     pdzmin  , pkth    , ph_co , pacr,     &   ! <<==  parameters for the tanh stretching
      &                     pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
      &                     pdept   , pdepw   ,                   &   !       3D t & w-points depth
      &                     pe3t    , pe3u    , pe3v    , pe3f,   &   !       vertical scale factors
      &                     pe3w    , pe3uw   , pe3vw           )     !           -      -      -
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_sco  ***
      !!                     
      !! ** Purpose :   TODO
      !!
      !! ** Method  :   TODO
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pdept_1d, pdepw_1d          ! 1D t- and w-depth                    [m]
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pbathy                      ! bathymetry                           [m]
      REAL(wp)                  , INTENT(in   ) ::   pHmax                       ! maximum depth   [m]
      REAL(wp)                  , INTENT(in   ) ::   pdzmin                      ! minimum value of e3 at the surface   [m]
      REAL(wp)                  , INTENT(in   ) ::   pkth                        ! position of the inflexion point
      REAL(wp)                  , INTENT(in   ) ::   ph_co                       ! layer thickness with z-coordinate [m]
      REAL(wp)                  , INTENT(in   ) ::   pacr                        ! slope of the tanh
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pe3t_1d, pe3w_1d            ! 1D t and w-scale factors (computed at 4000m)  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! 3D t and w-depth           [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         ! at t, u, v, f, w, uw, vw points
      !
      INTEGER  ::   ji, jj, jk                     ! do-loop index                           
      INTEGER  ::   kkconst                        ! index until z-coordinates are used
      REAL(wp), DIMENSION(2,jpk)   ::   zlev_dep   ! depth of the levels
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'zgr_lib: zgr_sco : define full vertical s-coord. system using the 2d bathymetry field'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~'
      !
      !!----------------------------------------------------------------------
      ! 1D profile
      zlev_dep = mi96_1d( pHmax, pdzmin, 0._wp, pkth, 0, pacr, ldprint=.TRUE. )
      pdepw_1d(:) = zlev_dep(1,:)
      pdept_1d(:) = zlev_dep(2,:)
      !
      !                       !==  t- and w- scale factors from depth  ==!
      CALL depth_to_e3( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d )
      !
      ! nearest index and its depth value to the depth where s-coordinates should start
      kkconst = MINLOC(abs(pdepw_1d - ph_co), 1)
      !
      !!----------------------------------------------------------------------
      ! 3D profile
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         zlev_dep = mi96_1d( pbathy(ji,jj), pe3w_1d(kkconst), pdepw_1d(kkconst), pkth, kkconst-1, pacr )
         pdepw(ji,jj,:) = pdepw_1d(:)
         pdept(ji,jj,:) = pdept_1d(:)
         !
         DO jk = kkconst+1, jpk
            pdepw(ji,jj,jk) = zlev_dep(1,jk)
            pdept(ji,jj,jk) = zlev_dep(2,jk)
         END DO
      END_2D
      !
      !                       !==  t- and w- scale factors from depth  ==!
      CALL depth_to_e3( pdept   , pdepw   , pe3t   , pe3w    )
      !                       !==  update t- and w-depths from SUM(e3)  <== needed to get the same last digit in ht_0 calculation
      CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d ) 
      !
      CALL e3_to_depth( pe3t   ,  pe3w  , pdept   , pdepw    ) 
      !
      !                       !==  other scale factors  ==!
      CALL e3tw_to_other_e3( pe3t , pe3w ,         &   ! <<==  e3 at t- and w-levels
         &                   pe3u , pe3v , pe3f,   &   ! ==>>  e3 at u- ,v- ,f-points
         &                   pe3uw, pe3vw        )     !         and uw-,vw-points
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
      !
   END SUBROUTINE zgr_sco_mi96 

   SUBROUTINE e3tw_to_other_e3( pe3t , pe3w ,         &   ! <<==  e3 at t- and w-levels
      &                         pe3u , pe3v , pe3f,   &   ! ==>>  e3 at u- ,v- ,f-points
      &                         pe3uw, pe3vw        )     !         and uw-,vw-points
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE e3tw_to_other_e3  ***
      !!                     
      !! ** Purpose :   set the vertical scale factors used in the computation
      !!
      !! ** Method  :   s-coordinate : scale factors are simple t-level (w-level) 
      !!              averaging from e3t (e3w) neighbours.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pe3t, pe3w           ! t- and w-scale factors   [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3u , pe3v , pe3f   ! vertical scale factors    [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3uw, pe3vw  ! at t, u, v, f, w, uw, vw points
      REAL(wp), DIMENSION(jpi,jpj,jpk)          ::   ze3fw  ! at fw points
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop argument
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '   e3tw_to_other_e3:   set e3u ,e3v , e3f from e3t'
      IF(lwp) WRITE(numout,*) '   ~~~~~~~~~~~~~~~~~   and e3uw,e3vw, e3w from e3w'
      !
      !                       !* set values except on last row or/and column
      !
      DO_3D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 , 1 , jpk)
         pe3u (ji,jj,jk) = 0.50_wp * ( pe3t(ji,jj,jk) + pe3t(ji+1,jj,jk) )   
         pe3uw(ji,jj,jk) = 0.50_wp * ( pe3w(ji,jj,jk) + pe3w(ji+1,jj,jk) )
         pe3v (ji,jj,jk) = 0.50_wp * ( pe3t(ji,jj,jk) + pe3t(ji,jj+1,jk) )
         pe3vw(ji,jj,jk) = 0.50_wp * ( pe3w(ji,jj,jk) + pe3w(ji,jj+1,jk) )
         pe3f (ji,jj,jk) = 0.25_wp * (  pe3t(ji  ,jj,jk) + pe3t(ji  ,jj+1,jk)  &
            &                         + pe3t(ji+1,jj,jk) + pe3t(ji+1,jj+1,jk)  )
         ze3fw(ji,jj,jk) = 0.25_wp * (  pe3w(ji  ,jj,jk) + pe3w(ji  ,jj+1,jk)  &
            &                         + pe3w(ji+1,jj,jk) + pe3w(ji+1,jj+1,jk)  )
      END_3D
      !
      !*  Apply lateral boundary condition   (use of optional argument cd_mpp)
      ! 
      IF( mig(1, nn_hls) == 1 ) THEN   ! first column of the local domain is the global domain one (which is closed here)
         !                     ! set first inner v-, vw-, t-point to corresponding u-, uw-, f-point  
         pe3v (1+nn_hls,:,:)       = pe3f (1+nn_hls,:,:)
         pe3vw(1+nn_hls,:,:)       = ze3fw(1+nn_hls,:,:)
         pe3t (1+nn_hls,1:jpj,:)   = pe3u (1+nn_hls,1:jpj,:)
         pe3w (1+nn_hls,1:jpj,:)   = pe3uw(1+nn_hls,1:jpj,:)
         !
      ENDIF
      !
      IF( mjg(1, nn_hls) == 1 ) THEN   ! last row of the local domain is the global domain one (which is closed here)
         !                            ! Extend inner domain value on the last row
         pe3u (:,1+nn_hls,:) = pe3f (:,1+nn_hls,:)
         pe3uw(:,1+nn_hls,:) = ze3fw(:,1+nn_hls,:)
         pe3t (:,1+nn_hls,:) = pe3v (:,1+nn_hls,:)
         pe3w (:,1+nn_hls,:) = pe3vw(:,1+nn_hls,:)
      ENDIF
      !
      CALL lbc_lnk( 'e3tw_to_other_e3', pe3u , 'U', 1., pe3v , 'V', 1., pe3f , 'F', 1.,     &
            &                           pe3uw, 'U', 1., pe3uw, 'V', 1., pe3t , 'T', 1.,     &
            &                           pe3w , 'T', 1., kfillmode = jpfillcopy )
      !
   END SUBROUTINE e3tw_to_other_e3


   FUNCTION mi96_1d( pH, pdzmin, ph_co, pkth, kkconst, pacr, ldprint ) RESULT( zlev_dep )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!                
      !!                
      !!                
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )             ::   pH       ! depth of the bottom   [meters]
      REAL(wp), INTENT(in   )             ::   pdzmin   ! value of e3 at the connection point   [meters]
      REAL(wp), INTENT(in   )             ::   ph_co    ! depth of the connection between z- and s-coordinates   [meters]
      REAL(wp), INTENT(in   )             ::   pkth     ! position of the inflexion point
      INTEGER , INTENT(in   )             ::   kkconst  ! number of levels with z-coordinates [meters]
      REAL(wp), INTENT(in   )             ::   pacr     ! slope of the tanh
      LOGICAL , INTENT(in   ), OPTIONAL   ::   ldprint  ! Whether to print infos about the arguments (default .false.)
      REAL(wp), DIMENSION(2,jpk)          ::   zlev_dep ! depth of the levels
      !!----------------------------------------------------------------------
      INTEGER  ::   jk
      REAL(wp) ::   za0, za1, zsur     ! Computed parameters for the tanh
      REAL(wp) ::   zw, zt             ! local scalar
      REAL(wp) ::   zh_co              ! local depth of the connection between z- and s-coordinates   [meters]
      !
      !IF( kkconst == 0 ) THEN
      !   zh_co = 0._wp
      !ELSE
      !   zh_co = ph_co
      !ENDIF
      !
      za1  = (  pdzmin - (pH - ph_co) / REAL(jpkm1-kkconst,wp)  )                                                      &
         & / ( TANH((1-pkth)/pacr) - pacr / REAL(jpkm1-kkconst) * (  LOG( COSH( (jpk - kkconst - pkth) / pacr) )      &
         &                                                         - LOG( COSH( ( 1            - pkth) / pacr) )  )  )
      za0  = pdzmin  - za1 *             TANH( (1-pkth) / pacr )
      zsur =   - za0 - za1 * pacr * LOG( COSH( (1-pkth) / pacr )  )
      !
      IF( lwp .AND. PRESENT(ldprint) ) THEN                         ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    mi96_id : Vertical coordinate ... TODO'
         WRITE(numout,*) '    ~~~~~~~'
         WRITE(numout,*) '       Depth of the bottom                                   pH    = ', pH
         WRITE(numout,*) '       Value of e3 in the pure z-coordinate'
         WRITE(numout,*) '          (= value of e3 at the first s-coordinate level)    pdzmin   = ', pdzmin
         WRITE(numout,*) '       Position of the inflexion point                       pkth     = ', pkth
         WRITE(numout,*) '       Number of pure z levels                               kkconst  = ', kkconst
         WRITE(numout,*) '       Slope of the tanh                                     pacr     = ', pacr
      ENDIF
      !
      DO jk = kkconst+1, jpk          ! depth at T and W-points
         zw = REAL( jk-kkconst, wp )
         zt = REAL( jk-kkconst, wp ) + 0.5_wp
         !
         zlev_dep(1,jk) = zsur + za0 * zw + za1 * pacr *  LOG( COSH( (zw-pkth) / pacr ) ) + ph_co
         zlev_dep(2,jk) = zsur + za0 * zt + za1 * pacr *  LOG( COSH( (zt-pkth) / pacr ) ) + ph_co
      END DO
      !
   END FUNCTION mi96_1d
   
   SUBROUTINE depth_to_e3_1d( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE depth_to_e3_1d  ***
      !!
      !! ** Purpose :   compute e3t & e3w scale factors from t- & w-depths of model levels
      !!
      !! ** Method  :   The scale factors are given by the discrete derivative 
      !!              of the depth:
      !!                               e3w(jk) = dk[ dept_1d ] 
      !!                               e3t(jk) = dk[ depw_1d ]
      !!              with, at top and bottom :
      !!                      e3w( 1 ) = 2 * ( dept( 1 ) - depw( 1 ) )
      !!                      e3t(jpk) = 2 * ( dept(jpk) - depw(jpk) )   
      !!
      !! ** Action  : - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in   ) ::   pdept_1d, pdepw_1d   ! depths          [m]
      REAL(wp), DIMENSION(:), INTENT(  out) ::   pe3t_1d , pe3w_1d    ! e3.=dk[depth]   [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      ! use pdep. at w- and t-points to compute e3. (e3. = dk[depth])
      !
      pe3w_1d( 1 ) = 2._wp * ( pdept_1d(1) - pdepw_1d(1) ) 
      DO jk = 1, jpkm1
         pe3w_1d(jk+1) = pdept_1d(jk+1) - pdept_1d(jk) 
         pe3t_1d(jk  ) = pdepw_1d(jk+1) - pdepw_1d(jk) 
      END DO
      pe3t_1d(jpk) = 2._wp * ( pdept_1d(jpk) - pdepw_1d(jpk) )
      !
   END SUBROUTINE depth_to_e3_1d
   
      
   SUBROUTINE depth_to_e3_3d( pdept_3d, pdepw_3d, pe3t_3d, pe3w_3d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE depth_to_e3_3d  ***
      !!
      !! ** Purpose :   compute e3t & e3w scale factors from t- & w-depths of model levels
      !!
      !! ** Method  :   The scale factors are given by the discrete derivative 
      !!              of the depth:
      !!                               e3w(jk) = dk[ dept_1d ] 
      !!                               e3t(jk) = dk[ depw_1d ]
      !!              with, at top and bottom :
      !!                      e3w( 1 ) = 2 * ( dept( 1 ) - depw( 1 ) )
      !!                      e3t(jpk) = 2 * ( dept(jpk) - depw(jpk) )   
      !!
      !! ** Action  : - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdept_3d, pdepw_3d   ! depth           [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t_3d , pe3w_3d    ! e3.=dk[depth]   [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------      
      pe3w_3d(:,:, 1 ) = 2._wp * ( pdept_3d(:,:,1) - pdepw_3d(:,:,1) ) 
      DO jk = 1, jpkm1
         pe3w_3d(:,:,jk+1) = pdept_3d(:,:,jk+1) - pdept_3d(:,:,jk) 
         pe3t_3d(:,:,jk  ) = pdepw_3d(:,:,jk+1) - pdepw_3d(:,:,jk) 
      END DO
      pe3t_3d(:,:,jpk) = 2._wp * ( pdept_3d(:,:,jpk) - pdepw_3d(:,:,jpk) )   
      !
   END SUBROUTINE depth_to_e3_3d


   SUBROUTINE e3_to_depth_1d( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE e3_to_depth_1d  ***
      !!
      !! ** Purpose :   compute t- & w-depths of model levels from e3t & e3w scale factors
      !!
      !! ** Method  :   The t- & w-depth are given by the summation of e3w & e3t, resp. 
      !!
      !! ** Action  : - pe3t_1d, pe3w_1d : scale factor of t- and w-point (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in   ) ::   pe3t_1d , pe3w_1d    ! vert. scale factors   [m]
      REAL(wp), DIMENSION(:), INTENT(  out) ::   pdept_1d, pdepw_1d   ! depth = SUM( e3 )     [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      pdepw_1d(1) = 0.0_wp
      pdept_1d(1) = 0.5_wp * pe3w_1d(1)
      DO jk = 2, jpk
         pdepw_1d(jk) = pdepw_1d(jk-1) + pe3t_1d(jk-1) 
         pdept_1d(jk) = pdept_1d(jk-1) + pe3w_1d(jk  ) 
      END DO
      !
   END SUBROUTINE e3_to_depth_1d
   
      
   SUBROUTINE e3_to_depth_3d( pe3t_3d, pe3w_3d, pdept_3d, pdepw_3d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE e3_to_depth_3d  ***
      !!
      !! ** Purpose :   compute t- & w-depths of model levels from e3t & e3w scale factors
      !!
      !! ** Method  :   The t- & w-depth are given by the summation of e3w & e3t, resp. 
      !!
      !! ** Action  : - pe3t_1d, pe3w_1d : scale factor of t- and w-point (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pe3t_3d , pe3w_3d    ! vert. scale factors   [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept_3d, pdepw_3d   ! depth = SUM( e3 )     [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------      
      !
      pdepw_3d(:,:,1) = 0.0_wp
      pdept_3d(:,:,1) = 0.5_wp * pe3w_3d(:,:,1)
      DO jk = 2, jpk
         pdepw_3d(:,:,jk) = pdepw_3d(:,:,jk-1) + pe3t_3d(:,:,jk-1) 
         pdept_3d(:,:,jk) = pdept_3d(:,:,jk-1) + pe3w_3d(:,:,jk  ) 
      END DO
      !
   END SUBROUTINE e3_to_depth_3d

   !!======================================================================
END MODULE zgr_lib
