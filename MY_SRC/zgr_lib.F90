MODULE zgr_lib
   !!==============================================================================
   !!                       ***  MODULE usrlib_zgr   ***
   !! Ocean domain : library routines used to build a vertical coordinate 
   !!==============================================================================
   !! History :  4.0  ! 2019-05  ( R. Caneil, G. Madec)  Original code 
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

   !SUBROUTINE zgr_zps( pdept_1d, pdepw_1d,                   &   ! <==>  1D reference vertical coordinate
   !   &                pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
   !   &                pdept   , pdepw   ,                   &   ! ==>>  3D t & w-points depth
   !   &                pe3t    , pe3u    , pe3v    , pe3f,   &   ! ==>>  vertical scale factors
   !   &                pe3w    , pe3uw   , pe3vw           )     ! ==>>      
   !   !!----------------------------------------------------------------------
   !   !!                  ***  ROUTINE zgr_zps  ***
   !   !!                     
   !   !! ** Purpose :   the depth and vertical scale factor in partial step
   !   !!              reference z-coordinate case
   !   !!
   !   !! ** Method  :   Partial steps : computes the 3D vertical scale factors
   !   !!      of T-, U-, V-, W-, UW-, VW and F-points that are associated with
   !   !!      a partial step representation of bottom topography.
   !   !!
   !   !!        The reference depth of model levels is defined from an analytical
   !   !!      function the derivative of which gives the reference vertical
   !   !!      scale factors.
   !   !!        From  depth and scale factors reference, we compute there new value
   !   !!      with partial steps  on 3d arrays ( i, j, k ).
   !   !!
   !   !!              w-level:  pdepw(i,j,k)  =  pdep(k)
   !   !!                       e3w(i,j,k) = dk( pdep)(k)     = e3(i,j,k)
   !   !!              t-level:  pdept(i,j,k)  =  pdep(k+0.5)
   !   !!                       e3t(i,j,k) = dk( pdep)(k+0.5) = e3(i,j,k+0.5)
   !   !!
   !   !!        With the help of the bathymetric file ( bathymetry_depth_ORCA_R2.nc),
   !   !!      we find the mbathy index of the depth at each grid point.
   !   !!      This leads us to three cases:
   !   !!
   !   !!              - bathy = 0 => mbathy = 0
   !   !!              - 1 < mbathy < jpkm1    
   !   !!              - bathy >  pdepw(jpk) => mbathy = jpkm1  
   !   !!
   !   !!        Then, for each case, we find the new depth at t- and w- levels
   !   !!      and the new vertical scale factors at t-, u-, v-, w-, uw-, vw- 
   !   !!      and f-points.
   !   !! 
   !   !!        This routine is given as an example, it must be modified
   !   !!      following the user s desiderata. nevertheless, the output as
   !   !!      well as the way to compute the model levels and scale factors
   !   !!      must be respected in order to insure second order accuracy
   !   !!      schemes.
   !   !!
   !   !!         c a u t i o n :  pdept_1d,  pdepw_1d and e3._1d are positives
   !   !!         - - - - - - -    pdept,  pdepw and e3. are positives
   !   !!      
   !   !!  Reference :   Pacanowsky & Gnanadesikan 1997, Mon. Wea. Rev., 126, 3248-3270.
   !   !!----------------------------------------------------------------------
   !   REAL(wp), DIMENSION(:)    , INTENT(inout) ::   pdept_1d, pdepw_1d          ! 1D t- and w-depth          [m]
   !   REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pe3t_1d, pe3w_1d            ! 1D t and w-scale factors   [m]
   !   REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! 3D t and w-depth           [m]
   !   REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors     [m]
   !   REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         ! at t, u, v, f, w, uw, vw points
!
   !   INTEGER  ::   ji, jj, jk       ! dummy loop indices
   !   INTEGER  ::   ik, it, ikb, ikt ! temporary integers
   !   REAL(wp) ::   ze3tp , ze3wp    ! Last ocean level thickness at T- and W-points
   !   REAL(wp) ::   zdepwp, zdepth   ! Ajusted ocean depth to avoid too small e3t
   !   REAL(wp) ::   zdiff            ! temporary scalar
   !   REAL(wp) ::   zmax             ! temporary scalar
   !   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  zprt
   !   !!---------------------------------------------------------------------
   !   !
   !   ALLOCATE( zprt(jpi,jpj,jpk) )
   !   !
   !   IF(lwp) WRITE(numout,*)
   !   IF(lwp) WRITE(numout,*) '    zgr_zps : z-coordinate with partial steps'
   !   IF(lwp) WRITE(numout,*) '    ~~~~~~~ '
   !   IF(lwp) WRITE(numout,*) '              mbathy is recomputed : bathy_level file is NOT used'
!
   !   ! compute position of the ice shelf grounding line
   !   ! set bathy and isfdraft to 0 where grounded
   !   IF ( ln_isfcav ) CALL zgr_isf_zspace
!
   !   ! bathymetry in level (from bathy_meter)
   !   ! ===================
   !   zmax =  pdepw_1d(jpk) + pe3t_1d(jpk)        ! maximum depth (i.e. the last ocean level thickness <= 2*pe3t_1d(jpkm1) )
   !   bathy(:,:) = MIN( zmax ,  bathy(:,:) )    ! bounded value of bathy (min already set at the end of zgr_bat)
   !   WHERE( bathy(:,:) == 0._wp )   ;   mbathy(:,:) = 0       ! land  : set mbathy to 0
   !   ELSE WHERE                     ;   mbathy(:,:) = jpkm1   ! ocean : initialize mbathy to the max ocean level
   !   END WHERE
!
   !   ! Compute mbathy for ocean points (i.e. the number of ocean levels)
   !   ! find the number of ocean levels such that the last level thickness
   !   ! is larger than the minimum of pe3zps_min and pe3zps_rat * pe3t_1d (where
   !   ! pe3t_1d is the reference level thickness
   !   DO jk = jpkm1, 1, -1
   !      zdepth =  pdepw_1d(jk) + MIN( e3zps_min, pe3t_1d(jk)*e3zps_rat )
   !      WHERE( 0._wp < bathy(:,:) .AND. bathy(:,:) <= zdepth )   mbathy(:,:) = jk-1
   !   END DO
!
   !   ! Check compatibility between bathy and iceshelf draft
   !   ! insure at least 2 wet level on the vertical under an ice shelf
   !   ! compute misfdep and adjust isf draft if needed
   !   IF ( ln_isfcav ) CALL zgr_isf_kspace
!
   !   ! Scale factors and depth at T- and W-points
   !   DO jk = 1, jpk                        ! intitialization to the reference z-coordinate
   !       pdept(:,:,jk) =  pdept_1d(jk)
   !       pdepw(:,:,jk) =  pdepw_1d(jk)
   !      pe3t  (:,:,jk) = pe3t_1d  (jk)
   !      pe3w  (:,:,jk) = pe3w_1d  (jk)
   !   END DO
   !   
   !   ! Scale factors and depth at T- and W-points
   !   DO jj = 1, jpj
   !      DO ji = 1, jpi
   !         ik = mbathy(ji,jj)
   !         IF( ik > 0 ) THEN               ! ocean point only
   !            ! max ocean level case
   !            IF( ik == jpkm1 ) THEN
   !               zdepwp = bathy(ji,jj)
   !               ze3tp  = bathy(ji,jj) -  pdepw_1d(ik)
   !               ze3wp = 0.5_wp * pe3w_1d(ik) * ( 1._wp + ( ze3tp/e3t_1d(ik) ) )
   !               pe3t(ji,jj,ik  ) = ze3tp
   !               pe3t(ji,jj,ik+1) = ze3tp
   !               IF ( ln_e3_dep.AND.ln_dept_mid ) THEN
   !                   pdept(ji,jj,ik) =  pdepw_1d(ik) + 0.5_wp * ze3tp
   !                  pe3w(ji,jj,ik  ) =  pdept(ji,jj,ik) -  pdept(ji,jj,ik-1) 
   !               ELSE
   !                   pdept(ji,jj,ik) =  pdept_1d(ik-1) + ze3wp
   !                  pe3w(ji,jj,ik  ) = ze3wp
   !               ENDIF
   !                pdept(ji,jj,ik+1) =  pdept(ji,jj,ik) + ze3tp
   !               pe3w(ji,jj,ik+1) = ze3tp
   !                pdepw(ji,jj,ik+1) = zdepwp
   !               !
   !            ELSE                         ! standard case
   !               IF( bathy(ji,jj) <=  pdepw_1d(ik+1) ) THEN  ;    pdepw(ji,jj,ik+1) = bathy(ji,jj)
   !               ELSE                                       ;    pdepw(ji,jj,ik+1) =  pdepw_1d(ik+1)
   !               ENDIF
!gm! Bug?  check the  pdepw_1d
   !               !       ... on ik
   !               pe3t  (ji,jj,ik) = pe3t_1d  (ik) * (  pdepw (ji,jj,ik+1) -  pdepw_1d(ik) )   & 
   !                  &                             / (  pdepw_1d(      ik+1) -  pdepw_1d(ik) ) 
   !               IF ( ln_e3_dep.AND.ln_dept_mid ) THEN
   !                   pdept(ji,jj,ik) =  pdepw_1d(ik) + 0.5_wp * pe3t(ji,jj,ik)
   !                  pe3w(ji,jj,ik) =  pdept(ji,jj,ik) -  pdept(ji,jj,ik-1)
   !               ELSE
   !                   pdept(ji,jj,ik) =  pdepw_1d(ik) + (  pdepw(ji,jj,ik+1) -  pdepw_1d(ik) )   &
   !                     &                             * (( pdept_1d(     ik  ) -  pdepw_1d(ik) )   &
   !                     &                             / (  pdepw_1d(     ik+1) -  pdepw_1d(ik) ))
   !                  pe3w(ji,jj,ik) = 0.5_wp * (  pdepw(ji,jj,ik+1) +  pdepw_1d(ik+1) - 2._wp *  pdepw_1d(ik) )   &
   !                     &                     * ( pe3w_1d(ik) / (  pdepw_1d(ik+1) -  pdepw_1d(ik) ) )
   !               ENDIF
   !               !       ... on ik+1
   !               pe3w  (ji,jj,ik+1) = pe3t  (ji,jj,ik)
   !               pe3t  (ji,jj,ik+1) = pe3t  (ji,jj,ik)
   !                pdept(ji,jj,ik+1) =  pdept(ji,jj,ik) + pe3t(ji,jj,ik)
   !            ENDIF
   !         ENDIF
   !      END DO
   !   END DO
   !   !
   !   it = 0
   !   DO jj = 1, jpj
   !      DO ji = 1, jpi
   !         ik = mbathy(ji,jj)
   !         IF( ik > 0 ) THEN               ! ocean point only
   !            pe3tp (ji,jj) = pe3t(ji,jj,ik)
   !            pe3wp (ji,jj) = pe3w(ji,jj,ik)
   !            ! test
   !            zdiff=  pdepw(ji,jj,ik+1) -  pdept(ji,jj,ik  )
   !            IF( zdiff <= 0._wp .AND. lwp ) THEN 
   !               it = it + 1
   !               WRITE(numout,*) ' it      = ', it, ' ik      = ', ik, ' (i,j) = ', ji, jj
   !               WRITE(numout,*) ' bathy = ', bathy(ji,jj)
   !               WRITE(numout,*) '  pdept = ',  pdept(ji,jj,ik), '  pdepw = ',  pdepw(ji,jj,ik+1), ' zdiff = ', zdiff
   !               WRITE(numout,*) ' pe3tp    = ', pe3t  (ji,jj,ik), ' pe3wp    = ', pe3w  (ji,jj,ik  )
   !            ENDIF
   !         ENDIF
   !      END DO
   !   END DO
   !   !
   !   ! compute top scale factor if ice shelf
   !   IF (ln_isfcav) CALL zps_isf
   !   !
   !   ! Scale factors and depth at U-, V-, UW and VW-points
   !   DO jk = 1, jpk                        ! initialisation to z-scale factors
   !      pe3u (:,:,jk) = pe3t_1d(jk)
   !      pe3v (:,:,jk) = pe3t_1d(jk)
   !      pe3uw(:,:,jk) = pe3w_1d(jk)
   !      pe3vw(:,:,jk) = pe3w_1d(jk)
   !   END DO
!
   !   DO jk = 1,jpk                         ! Computed as the minimum of neighbooring scale factors
   !      DO jj = 1, jpj - 1
   !         DO ji = 1, jpi - 1   ! vector opt.
   !            pe3u (ji,jj,jk) = MIN( pe3t(ji,jj,jk), pe3t(ji+1,jj,jk) )
   !            pe3v (ji,jj,jk) = MIN( pe3t(ji,jj,jk), pe3t(ji,jj+1,jk) )
   !            pe3uw(ji,jj,jk) = MIN( pe3w(ji,jj,jk), pe3w(ji+1,jj,jk) )
   !            pe3vw(ji,jj,jk) = MIN( pe3w(ji,jj,jk), pe3w(ji,jj+1,jk) )
   !         END DO
   !      END DO
   !   END DO
!
   !   ! update pe3uw in case only 2 cells in the water column
   !   IF ( ln_isfcav ) CALL zps_isf_e3uv_w
   !   !
   !   CALL lbc_lnk('domzgr', pe3u , 'U', 1._wp )   ;   CALL lbc_lnk('domzgr', pe3uw, 'U', 1._wp )   ! lateral boundary conditions
   !   CALL lbc_lnk('domzgr', pe3v , 'V', 1._wp )   ;   CALL lbc_lnk('domzgr', pe3vw, 'V', 1._wp )
   !   !
   !   DO jk = 1, jpk                        ! set to z-scale factor if zero (i.e. along closed boundaries)
   !      WHERE( pe3u (:,:,jk) == 0._wp )   pe3u (:,:,jk) = pe3t_1d(jk)
   !      WHERE( pe3v (:,:,jk) == 0._wp )   pe3v (:,:,jk) = pe3t_1d(jk)
   !      WHERE( pe3uw(:,:,jk) == 0._wp )   pe3uw(:,:,jk) = pe3w_1d(jk)
   !      WHERE( pe3vw(:,:,jk) == 0._wp )   pe3vw(:,:,jk) = pe3w_1d(jk)
   !   END DO
   !   
   !   ! Scale factor at F-point
   !   DO jk = 1, jpk                        ! initialisation to z-scale factors
   !      pe3f(:,:,jk) = pe3t_1d(jk)
   !   END DO
   !   DO jk = 1, jpk                        ! Computed as the minimum of neighbooring V-scale factors
   !      DO jj = 1, jpj - 1
   !         DO ji = 1, jpi - 1   ! vector opt.
   !            pe3f(ji,jj,jk) = MIN( pe3v(ji,jj,jk), pe3v(ji+1,jj,jk) )
   !         END DO
   !      END DO
   !   END DO
   !   CALL lbc_lnk('domzgr', pe3f, 'F', 1._wp )       ! Lateral boundary conditions
   !   !
   !   DO jk = 1, jpk                        ! set to z-scale factor if zero (i.e. along closed boundaries)
   !      WHERE( pe3f(:,:,jk) == 0._wp )   pe3f(:,:,jk) = pe3t_1d(jk)
   !   END DO
!!g!m  bug ? :  must be a do loop with mj0,mj1
   !   ! 
   !   pe3t(:,mj0(1),:) = pe3t(:,mj0(2),:)     ! we duplicate factor scales for jj = 1 and jj = 2
   !   pe3w(:,mj0(1),:) = pe3w(:,mj0(2),:) 
   !   pe3u(:,mj0(1),:) = pe3u(:,mj0(2),:) 
   !   pe3v(:,mj0(1),:) = pe3v(:,mj0(2),:) 
   !   pe3f(:,mj0(1),:) = pe3f(:,mj0(2),:) 
!
   !   ! Control of the sign
   !   IF( MINVAL( pe3t  (:,:,:) ) <= 0._wp )   CALL ctl_stop( '    zgr_zps :   e r r o r   pe3t <= 0' )
   !   IF( MINVAL( pe3w  (:,:,:) ) <= 0._wp )   CALL ctl_stop( '    zgr_zps :   e r r o r   pe3w <= 0' )
   !   IF( MINVAL(  pdept(:,:,:) ) <  0._wp )   CALL ctl_stop( '    zgr_zps :   e r r o r    pdept <  0' )
   !   IF( MINVAL(  pdepw(:,:,:) ) <  0._wp )   CALL ctl_stop( '    zgr_zps :   e r r o r    pdepw <  0' )
   !   !
   !   ! if in the future gde3w need to be compute, use the function defined in NEMO
   !   ! for now gde3w computation is removed as not an output of domcfg
!
   !   DEALLOCATE( zprt )
   !   !
   !END SUBROUTINE zgr_zps

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
      ! IF( mig(jpi) == jpiglo ) THEN   ! last column of the local domain is the global domain one (which is closed here)
      !    !                            ! Extend inner domain value on the last column
      !    !WRITE(numout,*) '       e3 assigned           pe3u  =  ', pe3u(jpiglo - 2, jpjglo, jpk - 1)
      !    !WRITE(numout,*) '       for                   pe3u  =  ', pe3u(jpiglo - 1, jpjglo, jpk - 1)
      !    pe3v (jpi,:,:)       = pe3u (jpi,:,:)
      !    pe3vw(jpi,:,:)       = pe3uw(jpi,:,:)
      !    pe3t (jpi,1:jpj,:)   = pe3f (jpi,1:jpj,:)
      !    !WRITE(numout,*) '       altered?              pe3u  =  ', pe3u(jpiglo - 1, jpjglo, jpk - 1)
      ! ENDIF
      ! !
      ! IF( mjg(jpj) == jpjglo ) THEN   ! last row of the local domain is the global domain one (which is closed here)
      !    !                            ! Extend inner domain value on the last row
      !    pe3u (:,jpj,:) = pe3v (:,jpj,:)
      !    pe3uw(:,jpj,:) = pe3vw(:,jpj,:)
      !    pe3t (:,jpj,:) = pe3f (:,jpj,:)
      ! ENDIF
! 
      IF( mig(1) == 1 ) THEN   ! first column of the local domain is the global domain one (which is closed here)
         !                     ! set first inner v-, vw-, t-point to corresponding u-, uw-, f-point  
         pe3v (1+nn_hls,:,:)       = pe3f (1+nn_hls,:,:)
         pe3vw(1+nn_hls,:,:)       = ze3fw(1+nn_hls,:,:)
         pe3t (1+nn_hls,1:jpj,:)   = pe3u (1+nn_hls,1:jpj,:)
         pe3w (1+nn_hls,1:jpj,:)   = pe3uw(1+nn_hls,1:jpj,:)
         !
      ENDIF
      !
      IF( mjg(1) == 1 ) THEN   ! last row of the local domain is the global domain one (which is closed here)
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

!    SUBROUTINE zgr_bat_env    !    ( pbat, pbat_env )
!       !!----------------------------------------------------------------------
!       !!                  ***  ROUTINE zgr_sco  ***
!       !!                     
!       !! ** Purpose :   define an envelop 
!       !!
!       !! ** Method  :   s-coordinate
!       !!----------------------------------------------------------------------
      
!       !!
!       REAL(wp) ::   rn_rmax, rn_sbot_max, rn_sbot_min    !!rc TODO argument
!       !
!       INTEGER  ::   ji, jj, jk, jl           ! dummy loop argument
!       INTEGER  ::   iip1, ijp1, iim1, ijm1   ! temporary integers
!       INTEGER  ::   ios                      ! Local integer output status for namelist read
!       REAL(wp) ::   zrmax, ztaper   ! temporary scalars
!       REAL(wp) ::   zrfact
!       !
!       REAL(wp), DIMENSION(jpi,jpj) :: ztmpi1, ztmpi2, ztmpj1, ztmpj2
!       REAL(wp), DIMENSION(jpi,jpj) :: zenv, ztmp, zmsk, zri, zrj, zhbat
!       REAL(wp), DIMENSION(jpi,jpj) :: bathy, hbatt, hbatu, hbatv, hbatf
!       REAL(wp), DIMENSION(jpi,jpj) :: hift, hifu, hifv, hiff, scosrf, scobot
!       !!
!       !!----------------------------------------------------------------------
!       !

!       hift(:,:) = rn_sbot_min                     ! set the minimum depth for the s-coordinate
!       hifu(:,:) = rn_sbot_min
!       hifv(:,:) = rn_sbot_min
!       hiff(:,:) = rn_sbot_min

!       !                                        ! set maximum ocean depth
!       bathy(:,:) = MIN( rn_sbot_max, bathy(:,:) )

!          DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!             IF( bathy(ji,jj) > 0._wp )   bathy(ji,jj) = MAX( rn_sbot_min, bathy(ji,jj) )
!          END_2D
!       !                                        ! =============================
!       !                                        ! Define the envelop bathymetry   (hbatt)
!       !                                        ! =============================
!       ! use r-value to create hybrid coordinates
!       zenv(:,:) = bathy(:,:)
!       !
!       ! set first land point adjacent to a wet cell to sbot_min as this needs to be included in smoothing
!          DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!             IF( bathy(ji,jj) == 0._wp ) THEN
!                iip1 = MIN( ji+1, jpi )
!                ijp1 = MIN( jj+1, jpj )
!                iim1 = MAX( ji-1, 1 )
!                ijm1 = MAX( jj-1, 1 )
! !!gm BUG fix see ticket #1617
!                IF( ( + bathy(iim1,ijm1) + bathy(ji,ijp1) + bathy(iip1,ijp1)              &
!                   &  + bathy(iim1,jj  )                  + bathy(iip1,jj  )              &
!                   &  + bathy(iim1,ijm1) + bathy(ji,ijm1) + bathy(iip1,ijp1)  ) > 0._wp ) &
!                   &    zenv(ji,jj) = rn_sbot_min
! !!gm
! !!gm               IF( ( bathy(iip1,jj  ) + bathy(iim1,jj  ) + bathy(ji,ijp1  ) + bathy(ji,ijm1) +         &
! !!gm                  &  bathy(iip1,ijp1) + bathy(iim1,ijm1) + bathy(iip1,ijp1) + bathy(iim1,ijm1)) > 0._wp ) THEN
! !!gm               zenv(ji,jj) = rn_sbot_min
! !!gm             ENDIF
! !!gm end
!             ENDIF
!          END_2D

!       ! apply lateral boundary condition   CAUTION: keep the value when the lbc field is zero
!       CALL lbc_lnk( 'zgr_bat_env', zenv, 'T', 1._wp)
!       ! 
!       ! smooth the bathymetry (if required)
!       scosrf(:,:) = 0._wp             ! ocean surface depth (here zero: no under ice-shelf sea)
!       scobot(:,:) = bathy(:,:)        ! ocean bottom  depth
!       !
!       jl = 0
!       zrmax = 1._wp
!       !   
!       !     
!       ! set scaling factor used in reducing vertical gradients
!       zrfact = ( 1._wp - rn_rmax ) / ( 1._wp + rn_rmax )
!       !
!       ! initialise temporary evelope depth arrays
!       ztmpi1(:,:) = zenv(:,:)
!       ztmpi2(:,:) = zenv(:,:)
!       ztmpj1(:,:) = zenv(:,:)
!       ztmpj2(:,:) = zenv(:,:)
!       !
!       ! initialise temporary r-value arrays
!       zri(:,:) = 1._wp
!       zrj(:,:) = 1._wp
!       !                                                            ! ================ !
!       DO WHILE( jl <= 10000 .AND. ( zrmax - rn_rmax ) > 1.e-8_wp ) !  Iterative loop  !
!          !                                                         ! ================ !
!          jl = jl + 1
!          zrmax = 0._wp
!          ! we set zrmax from previous r-values (zri and zrj) first
!          ! if set after current r-value calculation (as previously)
!          ! we could exit DO WHILE prematurely before checking r-value
!          ! of current zenv
!          DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!             zrmax = MAX( zrmax, ABS(zri(ji,jj)), ABS(zrj(ji,jj)) )
!          END_2D
!          zri(:,:) = 0._wp
!          zrj(:,:) = 0._wp
!          DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!             iip1 = MIN( ji+1, jpi )      ! force zri = 0 on last line (ji=ncli+1 to jpi)
!             ijp1 = MIN( jj+1, jpj )      ! force zrj = 0 on last raw  (jj=nclj+1 to jpj)
!             IF( (zenv(ji,jj) > 0._wp) .AND. (zenv(iip1,jj) > 0._wp)) THEN
!                zri(ji,jj) = ( zenv(iip1,jj  ) - zenv(ji,jj) ) / ( zenv(iip1,jj  ) + zenv(ji,jj) )
!             END IF
!             IF( (zenv(ji,jj) > 0._wp) .AND. (zenv(ji,ijp1) > 0._wp)) THEN
!                zrj(ji,jj) = ( zenv(ji  ,ijp1) - zenv(ji,jj) ) / ( zenv(ji  ,ijp1) + zenv(ji,jj) )
!             END IF
!             IF( zri(ji,jj) >  rn_rmax )   ztmpi1(ji  ,jj  ) = zenv(iip1,jj  ) * zrfact
!             IF( zri(ji,jj) < -rn_rmax )   ztmpi2(iip1,jj  ) = zenv(ji  ,jj  ) * zrfact
!             IF( zrj(ji,jj) >  rn_rmax )   ztmpj1(ji  ,jj  ) = zenv(ji  ,ijp1) * zrfact
!             IF( zrj(ji,jj) < -rn_rmax )   ztmpj2(ji  ,ijp1) = zenv(ji  ,jj  ) * zrfact
!          END_2D
!          IF( lk_mpp )   CALL mpp_max( 'zgr_bat_env', zrmax )   ! max over the global domain
!          !
!          IF(lwp)WRITE(numout,*) 'zgr_sco :   iter= ',jl, ' rmax= ', zrmax
!          !
         
!          DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!             zenv(ji,jj) = MAX(zenv(ji,jj), ztmpi1(ji,jj), ztmpi2(ji,jj), ztmpj1(ji,jj), ztmpj2(ji,jj) )
!          END_2D
!          ! apply lateral boundary condition   CAUTION: keep the value when the lbc field is zero
!          CALL lbc_lnk( 'zgr_bat_env',  zenv, 'T', 1._wp)
!          !                                                  ! ================ !
!       END DO                                                !     End loop     !
!       !                                                     ! ================ !
!       DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!          zenv(ji,jj) = MAX( zenv(ji,jj), rn_sbot_min ) ! set all points to avoid undefined scale value warnings
!       END_2D
!       !
!       ! Envelope bathymetry saved in hbatt
!       hbatt(:,:) = zenv(:,:) 
!       IF( MINVAL( gphit(:,:) ) * MAXVAL( gphit(:,:) ) <= 0._wp ) THEN
!          CALL ctl_warn( ' s-coordinates are tapered in vicinity of the Equator' )
!          DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!             ztaper = EXP( -(gphit(ji,jj)/8._wp)**2._wp )
!             hbatt(ji,jj) = rn_sbot_max * ztaper + hbatt(ji,jj) * ( 1._wp - ztaper )
!          END_2D
!       ENDIF
!       !
      
!             !
!       !
!       !                                        ! ==============================
!       !                                        !   hbatu, hbatv, hbatf fields
!       !                                        ! ==============================
!       IF(lwp) THEN
!          WRITE(numout,*)
!            WRITE(numout,*) ' zgr_sco: minimum depth of the envelop topography set to : ', rn_sbot_min
!       ENDIF
!       hbatu(:,:) = rn_sbot_min
!       hbatv(:,:) = rn_sbot_min
!       hbatf(:,:) = rn_sbot_min
!       DO jj = 1, jpj - 1
!         DO ji = 1, jpi - 1   ! NO vector opt.
!            hbatu(ji,jj) = 0.50_wp * ( hbatt(ji  ,jj) + hbatt(ji+1,jj  ) )
!            hbatv(ji,jj) = 0.50_wp * ( hbatt(ji  ,jj) + hbatt(ji  ,jj+1) )
!            hbatf(ji,jj) = 0.25_wp * ( hbatt(ji  ,jj) + hbatt(ji  ,jj+1)   &
!               &                     + hbatt(ji+1,jj) + hbatt(ji+1,jj+1) )
!         END DO
!       END DO

!       ! 
!       ! Apply lateral boundary condition
! !!gm  ! CAUTION: retain non zero value in the initial file this should be OK for orca cfg, not for EEL
!       zhbat(:,:) = hbatu(:,:)   ;   CALL lbc_lnk( 'zgr_bat_env',  hbatu, 'U', 1._wp )
!       DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!          IF( hbatu(ji,jj) == 0._wp ) THEN
!             !No worries about the following line when ln_wd == .true.
!             IF( zhbat(ji,jj) == 0._wp )   hbatu(ji,jj) = rn_sbot_min
!             IF( zhbat(ji,jj) /= 0._wp )   hbatu(ji,jj) = zhbat(ji,jj)
!          ENDIF
!       END_2D
!       zhbat(:,:) = hbatv(:,:)   ;   CALL lbc_lnk( 'zgr_bat_env',  hbatv, 'V', 1._wp )
!       DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!          IF( hbatv(ji,jj) == 0._wp ) THEN
!             IF( zhbat(ji,jj) == 0._wp )   hbatv(ji,jj) = rn_sbot_min
!             IF( zhbat(ji,jj) /= 0._wp )   hbatv(ji,jj) = zhbat(ji,jj)
!          ENDIF
!       END_2D
!       zhbat(:,:) = hbatf(:,:)   ;   CALL lbc_lnk( 'zgr_bat_env',  hbatf, 'F', 1._wp )
!       DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!          IF( hbatf(ji,jj) == 0._wp ) THEN
!             IF( zhbat(ji,jj) == 0._wp )   hbatf(ji,jj) = rn_sbot_min
!             IF( zhbat(ji,jj) /= 0._wp )   hbatf(ji,jj) = zhbat(ji,jj)
!          ENDIF
!       END_2D

! !!bug:  key_helsinki a verifer
!         hift(:,:) = MIN( hift(:,:), hbatt(:,:) )
!         hifu(:,:) = MIN( hifu(:,:), hbatu(:,:) )
!         hifv(:,:) = MIN( hifv(:,:), hbatv(:,:) )
!         hiff(:,:) = MIN( hiff(:,:), hbatf(:,:) )

!       IF(lwp )   THEN
!          WRITE(numout,*) ' MAX val hif   t ', MAXVAL( hift (:,:) ), ' f ', MAXVAL( hiff (:,:) ),  &
!             &                        ' u ',   MAXVAL( hifu (:,:) ), ' v ', MAXVAL( hifv (:,:) )
!          WRITE(numout,*) ' MIN val hif   t ', MINVAL( hift (:,:) ), ' f ', MINVAL( hiff (:,:) ),  &
!             &                        ' u ',   MINVAL( hifu (:,:) ), ' v ', MINVAL( hifv (:,:) )
!          WRITE(numout,*) ' MAX val hbat  t ', MAXVAL( hbatt(:,:) ), ' f ', MAXVAL( hbatf(:,:) ),  &
!             &                        ' u ',   MAXVAL( hbatu(:,:) ), ' v ', MAXVAL( hbatv(:,:) )
!          WRITE(numout,*) ' MIN val hbat  t ', MINVAL( hbatt(:,:) ), ' f ', MINVAL( hbatf(:,:) ),  &
!             &                        ' u ',   MINVAL( hbatu(:,:) ), ' v ', MINVAL( hbatv(:,:) )
!       ENDIF
! !! helsinki
!       !
!    END SUBROUTINE zgr_bat_env
   
   
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
