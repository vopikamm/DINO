MODULE usrdef_hgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_hgr  ***
   !!
   !!                      ===  CANAL configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  NEMO  ! 2017-11  (J. Chanut)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr    : initialize the horizontal mesh for CANAL configuration
   !!----------------------------------------------------------------------
   USE dom_oce, ONLY: mig0, mjg0        ! ocean space and time domain
   USE par_oce ! , ONLY:                ! ocean space and time domain
   USE phycst          ! physical constants
   !
   USE usrdef_nam, ONLY: rn_e1_deg, rn_lam_min, rn_phi_min, nn_jeq_s
   !                                                  and reference latitude
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_hgr   ! called by domhgr.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_hgr.F90 10074 2018-08-28 16:15:49Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_hgr( plamt , plamu , plamv  , plamf  ,   &   ! geographic position (required)
      &                    pphit , pphiu , pphiv  , pphif  ,   &   !
      &                    kff   , pff_f , pff_t  ,            &   ! Coriolis parameter  (if domain not on the sphere)
      &                    pe1t  , pe1u  , pe1v   , pe1f   ,   &   ! scale factors       (required)
      &                    pe2t  , pe2u  , pe2v   , pe2f   ,   &   !
      &                    ke1e2u_v      , pe1e2u , pe1e2v     )   ! u- & v-surfaces (if gridsize reduction is used in strait(s))
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE usr_def_hgr  ***
      !!
      !! ** Purpose :   user defined mesh and Coriolis parameter
      !!
      !! ** Method  :   set all intent(out) argument to a proper value
      !!                CANAL configuration : beta-plance with uniform grid spacing (rn_dx)
      !!
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in degrees) 
      !!              - define coriolis parameter at f-point if the domain in not on the sphere (on beta-plane)
      !!              - define i- & j-scale factors at t-, u-, v- and f-points (in meters)
      !!              - define u- & v-surfaces (if gridsize reduction is used in some straits) (in m2)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   plamt, plamu, plamv, plamf   ! longitude outputs                     [degrees]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pphit, pphiu, pphiv, pphif   ! latitude outputs                      [degrees]
      INTEGER                 , INTENT(out) ::   kff                          ! =1 Coriolis parameter computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pff_f, pff_t                 ! Coriolis factor at f-point                [1/s]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1t, pe1u, pe1v, pe1f       ! i-scale factors                             [m]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe2t, pe2u, pe2v, pe2f       ! j-scale factors                             [m]
      INTEGER                 , INTENT(out) ::   ke1e2u_v                     ! =1 u- & v-surfaces computed here, =0 otherwise 
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1e2u, pe1e2v               ! u- & v-surfaces (if reduction in strait)   [m2]
      !
      INTEGER  ::   ji, jj       ! dummy loop indices
      REAL(wp) ::   zlam0        ! position of left longitude [degrees]
      REAL(wp) ::   zarg, zjeq   ! local variables
      !INTEGER  ::   ijeq         ! local integer
      REAL(wp) ::   zti, zui, zvi, zfi, ztj, zuj, zvj, zfj ! local variables
      
      !!-------------------------------------------------------------------------------
      !
      zlam0 = rn_lam_min - rn_e1_deg * 0.5_wp   ! so that the first U point is at o degrees of longitude
                           !==  geographical mesh on the sphere, isotropic MERCATOR type  ==!
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : BASIN configuration'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          geographical mesh on the sphere, MERCATOR type'
      IF(lwp) WRITE(numout,*) '          longitudinal/latitudinal spacing given by rn_e1_deg'
      IF ( rn_phi_min == -90 )   CALL ctl_stop( ' Mercator grid cannot start at south pole !!!! ' )
      !
      !  Find index corresponding to the equator, given the grid spacing rn_e1_deg
      !  and the (approximate) southern latitude rn_phi0.
      !  This way we ensure that the equator is at a "T / U" point, when in the domain.
      !  The formula should work even if the equator is outside the domain.      !
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      
      !zarg = rpi / 4. - rpi / 180. * rn_phi0 / 2.
      !zjeq = ABS( 180./rpi * LOG( COS( zarg ) / SIN( zarg ) ) / rn_e1_deg )
      !IF(  rn_phi0 > 0 )  zjeq = -zjeq
      !zjeq =  zjeq + 1._wp ! We add the +1 because the j indexes start from 1: if the equator is on the first row, its index must be 1
      !ijeq = FLOOR( zjeq )
      !IF( ABS( REAL( ijeq, wp ) - zjeq ) > 0.5 )   ijeq = ijeq + 1
      !!
      !IF(lwp) WRITE(numout,*) '          Index of the equator      (real)       on the MERCATOR grid:', zjeq
      !IF(lwp) WRITE(numout,*) '          Index of the equator (nearest integer) on the MERCATOR grid:', ijeq

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      !!
      !DO jj = 1, jpj
      !   DO ji = 1, jpi
      !      zti = REAL( ji - 1 + nimpp - 1 )         ;   ztj = REAL( jj - ijeq + njmpp - 1 )
      !      zui = REAL( ji - 1 + nimpp - 1 ) + 0.5   ;   zuj = REAL( jj - ijeq + njmpp - 1 )
      !      zvi = REAL( ji - 1 + nimpp - 1 )         ;   zvj = REAL( jj - ijeq + njmpp - 1 ) + 0.5
      !      zfi = REAL( ji - 1 + nimpp - 1 ) + 0.5   ;   zfj = REAL( jj - ijeq + njmpp - 1 ) + 0.5
      !   ! Longitude
      !      plamt(ji,jj) = zlam0 + rn_e1_deg * zti
      !      plamu(ji,jj) = zlam0 + rn_e1_deg * zui
      !      plamv(ji,jj) = zlam0 + rn_e1_deg * zvi
      !      plamf(ji,jj) = zlam0 + rn_e1_deg * zfi
      !   ! Latitude
      !      pphit(ji,jj) = 1./rad * ASIN ( TANH( rn_e1_deg *rad* ztj ) )
      !      pphiu(ji,jj) = 1./rad * ASIN ( TANH( rn_e1_deg *rad* zuj ) )
      !      pphiv(ji,jj) = 1./rad * ASIN ( TANH( rn_e1_deg *rad* zvj ) )
      !      pphif(ji,jj) = 1./rad * ASIN ( TANH( rn_e1_deg *rad* zfj ) )
      !   ! e1
      !      pe1t(ji,jj) = ra * rad * COS( rad * pphit(ji,jj) ) * rn_e1_deg
      !      pe1u(ji,jj) = ra * rad * COS( rad * pphiu(ji,jj) ) * rn_e1_deg
      !      pe1v(ji,jj) = ra * rad * COS( rad * pphiv(ji,jj) ) * rn_e1_deg
      !      pe1f(ji,jj) = ra * rad * COS( rad * pphif(ji,jj) ) * rn_e1_deg
      !   ! e2
      !      pe2t(ji,jj) = ra * rad * COS( rad * pphit(ji,jj) ) * rn_e1_deg
      !      pe2u(ji,jj) = ra * rad * COS( rad * pphiu(ji,jj) ) * rn_e1_deg
      !      pe2v(ji,jj) = ra * rad * COS( rad * pphiv(ji,jj) ) * rn_e1_deg
      !      pe2f(ji,jj) = ra * rad * COS( rad * pphif(ji,jj) ) * rn_e1_deg
      !   END DO
      !END DO
      !
      !
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            zti = REAL( mig0(ji) - 1, wp )             ;     ztj = REAL( mjg0(jj) - nn_jeq_s, wp )
            zui = REAL( mig0(ji) - 1, wp ) + 0.5       ;     zuj = REAL( mjg0(jj) - nn_jeq_s, wp )
            zvi = REAL( mig0(ji) - 1, wp )             ;     zvj = REAL( mjg0(jj) - nn_jeq_s, wp ) + 0.5
            zfi = REAL( mig0(ji) - 1, wp ) + 0.5       ;     zfj = REAL( mjg0(jj) - nn_jeq_s, wp ) + 0.5
         ! Longitude
            plamt(ji,jj) = zlam0 + rn_e1_deg * zti
            plamu(ji,jj) = zlam0 + rn_e1_deg * zui
            plamv(ji,jj) = zlam0 + rn_e1_deg * zvi
            plamf(ji,jj) = zlam0 + rn_e1_deg * zfi
         ! Latitude
            pphit(ji,jj) = 1./rad * ASIN ( TANH( rn_e1_deg *rad* ztj ) )
            pphiu(ji,jj) = 1./rad * ASIN ( TANH( rn_e1_deg *rad* zuj ) )
            pphiv(ji,jj) = 1./rad * ASIN ( TANH( rn_e1_deg *rad* zvj ) )
            pphif(ji,jj) = 1./rad * ASIN ( TANH( rn_e1_deg *rad* zfj ) )
         ! e1
            pe1t(ji,jj) = ra * rad * COS( rad * pphit(ji,jj) ) * rn_e1_deg
            pe1u(ji,jj) = ra * rad * COS( rad * pphiu(ji,jj) ) * rn_e1_deg
            pe1v(ji,jj) = ra * rad * COS( rad * pphiv(ji,jj) ) * rn_e1_deg
            pe1f(ji,jj) = ra * rad * COS( rad * pphif(ji,jj) ) * rn_e1_deg
         ! e2
            pe2t(ji,jj) = ra * rad * COS( rad * pphit(ji,jj) ) * rn_e1_deg
            pe2u(ji,jj) = ra * rad * COS( rad * pphiu(ji,jj) ) * rn_e1_deg
            pe2v(ji,jj) = ra * rad * COS( rad * pphiv(ji,jj) ) * rn_e1_deg
            pe2f(ji,jj) = ra * rad * COS( rad * pphif(ji,jj) ) * rn_e1_deg
      END_2D

      !                             ! NO reduction of grid size in some straits 
      ke1e2u_v = 0                  !    ==>> u_ & v_surfaces will be computed in dom_hgr routine
      pe1e2u(:,:) = 0._wp           !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp           !             require an initialization of INTENT(out) arguments
      !
      !
      IF( lwp .AND. .NOT.ln_rstart ) THEN      ! Control print : Grid informations (if not restart)
         WRITE(numout,*)
         WRITE(numout,*) '          longitude and e1 scale factors'
         WRITE(numout,*) '          ------------------------------'
         WRITE(numout,9300) ( ji, plamt(ji,1), plamu(ji,1),   &
            plamv(ji,1), plamf(ji,1),   &
            pe1t(ji,1), pe1u(ji,1),   &
            pe1v(ji,1), pe1f(ji,1), ji = 1, jpi,10)
9300     FORMAT( 1x, i4, f8.2,1x, f8.2,1x, f8.2,1x, f8.2, 1x,    &
            f19.10, 1x, f19.10, 1x, f19.10, 1x, f19.10 )
            !
         WRITE(numout,*)
         WRITE(numout,*) '          latitude and e2 scale factors'
         WRITE(numout,*) '          -----------------------------'
         WRITE(numout,9300) ( jj, pphit(1,jj), pphiu(1,jj),   &
            &                     pphiv(1,jj), pphif(1,jj),   &
            &                     pe2t  (1,jj), pe2u  (1,jj),   &
            &                     pe2v  (1,jj), pe2f  (1,jj), jj = 1, jpj, 10 )
      ENDIF
      
      !                       !==  Coriolis parameter  ==!
      kff = 1                       !  indicate not to compute Coriolis parameter afterward
      !                   mesh on the sphere
      !
      pff_f(:,:) = 2. * omega * SIN( rad * pphif(:,:) ) 
      pff_t(:,:) = 2. * omega * SIN( rad * pphit(:,:) )     !    -        -       - at t-point
      !
      !
    END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
