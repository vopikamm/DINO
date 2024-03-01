MODULE zanna_bolton
   !!======================================================================
   !!                       ***  MODULE  zanna_bolton  ***
   !! Ocean physics:  Baroclinic Zanna & Bolton parameterization (2020)
   !!                 see eq. 6 in https://laurezanna.github.io/files/Zanna-Bolton-2020.pdf
   !!======================================================================
   !! History : NEMO
   !!            4.2   ! 2023-10  (D.Kamm) following discretization by P.Perezhogin
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE lib_mpp         ! MPP library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager ! I/O manager
   USE phycst          ! rpi

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ZB_2020_init   ! routine called in nemogcm.F90 module
   PUBLIC   ZB_apply        ! routine called in step.F90 module

   !!------------- parameters -------------!!
   LOGICAL  :: ln_zanna_bolton    ! key for ZB
   REAL(wp) :: rn_gamma           ! ZB-coefficient

   !!------------- fields -----------------!!
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: T_xx   ! component of stress tensor
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: T_xy   ! component of stress tensor
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: T_yy   ! component of stress tensor
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: zbu   ! subgrid-forcing on u
   REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: zbv   ! subgrid-forcing on v

   !! * Substitutions
   ! for DO macro
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: module_example.F90 14842 2021-05-11 13:17:26Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ZB_2020_init ( )
      
      !INTEGER :: ji, jj, jk
      
      NAMELIST/namZB/ ln_zanna_bolton, rn_gamma
      
      !REWIND( numnam_cfg )          !  TODO??
      !READ  ( numnam_cfg, namZB )   !  TODO??

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ZB_init: Zanna & Bolton parameterizations - initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'

         WRITE(numout,*) '   Namelist namZB : set Zanna & Bolton parameters'
         WRITE(numout,*) '   ln_zanna_bolton =', ln_zanna_bolton
      END IF
      
      IF (not(ln_zanna_bolton)) THEN
         return
      END IF

      ALLOCATE(T_xx(jpi,jpj,jpk))
      ALLOCATE(T_xy(jpi,jpj,jpk))
      ALLOCATE(T_yy(jpi,jpj,jpk))
      ALLOCATE(zbu(jpi,jpj,jpk))
      ALLOCATE(zbv(jpi,jpj,jpk))


      ! TODO ??
      ! without rn_shlat
      ! ffmask = 0._wp
      ! DO jk = 1, jpk
      !    DO jj = 1, jpjm1
      !       DO ji = 1, jpim1
      !          ffmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
      !                           * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
      !       END DO
      !    END DO
      ! END DO
      ! CALL lbc_lnk( ffmask, 'F', 1._wp )

      ! init stress tensor and ZB tendency terms with zeros
      T_xx  = 0._wp
      T_xy  = 0._wp
      T_yy  = 0._wp
      zbu   = 0._wp
      zbv   = 0._wp

      !CALL KEB_rst_read() TODO ??
   END SUBROUTINE ZB_2020_init

   SUBROUTINE ZB_apply( kt, Kbb, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ZB_apply  ***   
      !! ** Purpose :   Apply ZB tendencies.
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt               ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kbb, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv 
      ! compute ZB from before velocities
      CALL ZB_2020(kt, Kbb, rn_gamma, puu(:,:,:,Kbb), pvv(:,:,:,Kbb), puu(:,:,:,Krhs), pvv(:,:,:,Krhs)) ! TODO namelist for gamma
   END SUBROUTINE ZB_apply

   SUBROUTINE ZB_2020(kt, Kbb, gamma, puu, pvv, zbu, zbv)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Zanna_Bolton_2020  ***
      !!                   
      !! ** Purpose :   Baroclinic Zanna-Bolton-2020 parameterization.
      !!                eq. 6 in https://laurezanna.github.io/files/Zanna-Bolton-2020.pdf
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! time step index
      INTEGER, INTENT(in) ::   Kbb  ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  puu, pvv ! before velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  zbu, zbv ! u-, v-component of S 
      
      INTEGER  ::   ji, jj, jk           ! dummy loop indices

      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   T_xx, Tyy, T_xy
      REAL(wp)                         ::   gamma     
      !
      ! Compute the stress tensor
      !
      CALL stress_tensor(kt, puu, pvv, T_xx, T_yy, T_xy)
      !
      ! Horizontal Divergence of the stress tensor
      DO_3D_OVR( 0, 0, 0, 0 , 1, jpk)
         zbu(ji,jj,jk) = zbu(ji,jj,jk) + (( T_xx(ji+1,jj,jk)   * e3t(ji+1,jj,jk,Kbb)   * e2t(ji+1,jj)   ** 2                   &
            &             - T_xx(ji,jj,jk) * e3t(ji,jj,jk,Kbb) * e2t(ji,jj) ** 2) / e2u(ji,jj)  &
            &           + ( T_xy(ji,jj,jk)   * e3f(ji,jj,jk)   * e1f(ji,jj)   ** 2                   &
            &             - T_xy(ji,jj-1,jk) * e3f(ji,jj-1,jk) * e1f(ji,jj-1) ** 2) / e1u(ji,jj)) &
            &           *   r1_e1u(ji,jj) * r1_e2u(ji,jj) / e3u(ji,jj,jk,Kbb) * gamma
         !
         zbv(ji,jj,jk) = zbv(ji,jj,jk) + (( T_xy(ji,jj,jk)   * e3f(ji,jj,jk)   * e2f(ji,jj)   ** 2                   &
            &             - T_xy(ji-1,jj,jk) * e3f(ji-1,jj,jk) * e2f(ji-1,jj) ** 2) / e2v(ji,jj)  &
            &           + ( T_yy(ji,jj+1,jk)   * e3t(ji,jj+1,jk,Kbb)   * e1t(ji,jj+1)   ** 2                   &
            &             - T_yy(ji,jj,jk) * e3t(ji,jj,jk,Kbb) * e1t(ji,jj) ** 2) / e1v(ji,jj)) &
            &           *   r1_e1v(ji,jj) * r1_e2v(ji,jj) / e3v(ji,jj,jk,Kbb) * gamma
      END_3D
   END SUBROUTINE ZB_2020

   SUBROUTINE stress_tensor(kt, pu, pv, T_xx, T_yy, T_xy)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE stress_tensor  ***
      !!                   
      !! ** Purpose :   computes the stress tensor for Zanna-Bolton.
      !!
      !! ** Method  :   The deformation rate is analogous to the routine 
      !!                ldf_dyn_smag in curvilinear coords (Griffies 2000):
      !!                du/dx ==> e2 d(u/e2)/dx;  du/dy ==> e1 d(u/e1)/dy; 
      !!                dv/dx ==> e2 d(v/e2)/dx;  dv/dy ==> e1 d(v/e1)/dy
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! time step index
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)      ::  pu, pv              ! before velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)   ::  T_xx, T_xy, T_yy    ! Stretching / Shearing deformation D
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk)    ::  sh_xx_t, sh_xx_f   ! Stretching deformation T/F-point
      REAL(wp), DIMENSION(jpi,jpj,jpk)    ::  sh_xy_t, sh_xy_f   ! Shearing deformation T/F-point
      REAL(wp), DIMENSION(jpi,jpj,jpk)    ::  vort_t, vort_f     ! relative vorticity T/F-point
      !
      INTEGER     ::   ji, jj, jk           ! dummy loop indices
      !
      REAL (wp)   ::   zux, zuy, zvx, zvy, zhyd_com, zdev_com   ! local scalars
      !
      DO jk = 1, jpkm1
         !
         ! Stretching deformation \tilde{D} on T-point
         DO_2D( 0, 0, 0, 0 )
            zux =   ( pu(ji,jj,jk) * r1_e2u(ji,jj) -  pu(ji-1,jj,jk) * r1_e2u(ji-1,jj) )  &
                 &                      * r1_e1t(ji,jj) * e2t(ji,jj)
            zvy =   ( pv(ji,jj,jk) * r1_e1v(ji,jj) -  pv(ji,jj-1,jk) * r1_e1v(ji,jj-1) )  &
                 &                      * r1_e2t(ji,jj) * e1t(ji,jj)
            sh_xx_t(ji,jj,jk) = ( zux - zvy ) * tmask(ji,jj,jk)
         END_2D
         ! Shearing deformation D and relative vorticity \Zeta on F-point
         DO_2D( 0, 0, 0, 0 )
            zuy =   (  pu(ji,jj+1,jk) * r1_e1u(ji,jj+1) -  pu(ji,jj,jk) * r1_e1u(ji,jj) )  &
                 &                        * r1_e2f(ji,jj)   * e1f(ji,jj)
            zvx =   ( pv(ji+1,jj,jk) * r1_e2v(ji+1,jj) -  pv(ji,jj,jk) * r1_e2v(ji,jj) )  &
                 &                        * r1_e1f(ji,jj)   * e2f(ji,jj)
            sh_xy_f(ji,jj,jk)   = ( zvx + zuy ) * fmask(ji,jj,jk)
            vort_f(ji,jj,jk) = ( zvx - zuy ) * fmask(ji,jj,jk)
         END_2D
      END DO
      !
      CALL lbc_lnk( 'stress_tensor', sh_xx_t, 'T', 1.0_wp )    ! B.C.: no-normal flow
      CALL lbc_lnk( 'stress_tensor', sh_xy_f, 'F', 1.0_wp )    ! B.C.: free-slip, i.e. du/dy=0 on the boundary
      CALL lbc_lnk( 'stress_tensor', vort_f , 'F', 1.0_wp )    ! B.C.: free-slip, i.e. du/dy=0 on the boundary
      !
      ! Interpolation on opposite grid-points (F-->T and T-->F)
      CALL f_on_t(sh_xy_f, sh_xy_t)
      CALL f_on_t(vort_f, vort_t)
      CALL t_on_f(sh_xx_t, sh_xx_f)
      !
      ! Compute components of the stress tensor
      DO_3D_OVR( 0, 0, 0, 0 , 1, jpk)
         zhyd_com = 0.5_wp * (vort_t(ji,jj,jk)**2 + sh_xy_t(ji,jj,jk)**2 + sh_xx_t(ji,jj,jk)**2)
         zdev_com = vort_t(ji,jj,jk) * sh_xy_t(ji,jj,jk)
         T_xx(ji,jj,jk) = - e1e2t(ji,jj) * (   zhyd_com - zdev_com ) * tmask(ji,jj,jk)
         T_yy(ji,jj,jk) = - e1e2t(ji,jj) * (   zhyd_com + zdev_com ) * tmask(ji,jj,jk)
         T_xy(ji,jj,jk) = - e1e2f(ji,jj) * ( vort_f(ji,jj,jk) * sh_xx_f(ji,jj,jk)) * fmask(ji,jj,jk)
      END_3D
      !
      CALL lbc_lnk( 'stress_tensor', T_xx, 'T', 1.0_wp )    ! B.C.: no-normal flow
      CALL lbc_lnk( 'stress_tensor', T_xx, 'T', 1.0_wp )    ! B.C.: free-slip, i.e. du/dy=0 on the boundary
      CALL lbc_lnk( 'stress_tensor', T_xy, 'T', 1.0_wp )    ! B.C.: free-slip, i.e. du/dy=0 on the boundary
      !
   END SUBROUTINE stress_tensor

   SUBROUTINE f_on_t(pf, pt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE f_on_t ***
      !!                   
      !! ** Purpose :   Interpolating fields from F-points on T-points. 
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)       ::  pf ! field on F-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(out)      ::  pt ! field on T-points
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      !
      DO_3D_OVR( 0, 0, 0, 0, 1, jpkm1)
         pt(ji,jj,jk) = 0.25 * (pf(ji,jj,jk) + pf(ji-1,jj,jk) + pf(ji,jj-1,jk) + pf(ji-1,jj-1,jk))
      END_3D
   END SUBROUTINE f_on_t

   SUBROUTINE t_on_f(pt, pf)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE t_on_f ***
      !!                   
      !! ** Purpose :   Interpolating fields from T-points on F-points. 
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)      ::  pt ! field on T-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(out)     ::  pf ! field on F-points
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      !
      DO_3D_OVR( 0, 0, 0, 0, 1, jpkm1)
         pf(ji,jj,jk) = 0.25 * (pt(ji,jj,jk) + pt(ji+1,jj,jk) + pt(ji,jj+1,jk) + pt(ji+1,jj+1,jk))
      END_3D
   END SUBROUTINE t_on_f

   !!======================================================================
END MODULE zanna_bolton
