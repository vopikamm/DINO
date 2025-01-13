MODULE KEB_operators
   !!==============================================================================
   !!                       ***  MODULE  KEB_operators  ***
   !!==============================================================================
   !! History : NEMO
   !!            3.6  ! 2019-05  (P.Perezhogin) original code for Cartesian coordinates
   !!                 (published in P. Perezhogin, RJNAMM, 2020)
   !!            3.6  ! 2019-10  (P. Perezhogin) curvilinear coordinates were added
   !!            3.6  ! 2021-06  (P. Perezhogin) preparing the code to Jones Colin
   !!----------------------------------------------------------------------

   USE dom_oce         ! ocean space and time domain
   USE lib_mpp         ! MPP library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE storng          ! generation of gaussian noise
   USE phycst          ! rpi
   !!----------------------------------------------------------------------
   !! All input arrays are assumed to be not NaN at the land and boundary 
   !! All output arrays can modify or define land or boundary points, but don't have to
   !! For safety, trends are passed with inout attribute
   !! KEB was tested in z-coordinates, but can work in s-coordinates with the following property:
   !! horizontal laplace and bilaplacian operators for tracers and momentum act along s=const surfaces
   !!----------------------------------------------------------------------
   
   IMPLICIT NONE   
   PUBLIC
   
CONTAINS

   SUBROUTINE z_filter( ptn, pta ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE z_filter  ***
      !!                   
      !! ** Purpose :   apply filtering in z direction 
      !!                if advection becomes to be unstable due 
      !!                to surface condition
      !!
      !! ** Method  :   apply diffusion in z direction with neumann b.c.
      !!                and CONSTANT diffusivity given by nullifing 2h-wave
      !!                on the surface. Its action reduces with depth
      !!                because mesh spacing increases.
      !!                Diffusion:
      !!                coef = nu * dt = dz^2/4 at the surface
      !!----------------------------------------------------------------------
      
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  ptn           ! tracer in T-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  pta           ! filtered tracer in T-points
      !
      INTEGER  :: ji, jj, jk                         ! dummy loop indices
      REAL(wp) :: coef
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwz        ! flux of tracer in w-points
      !!----------------------------------------------------------------------
      
      zwz(:,:,1) = 0._wp ! zero cross-surface flux   
      DO jk = 2, jpk
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               coef = 0.25_wp * e3w_0(ji,jj,1)**2
               zwz(ji,jj,jk) = coef * (ptn(ji,jj,jk-1) - ptn(ji,jj,jk))  / e3w_0(ji,jj,jk) * wmask(ji,jj,jk)
            END DO
         END DO    
      END DO
      
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
                  pta(ji,jj,jk) = ptn(ji,jj,jk)  + ( zwz(ji,jj,jk) - zwz(ji,jj,jk+1) ) / e3t_0(ji,jj,jk)
            END DO
         END DO    
      END DO
      
   END SUBROUTINE z_filter

   SUBROUTINE upwind_advection( ptn, pun, pvn, pwn, rhsn, bc_type ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE upwind_advection  ***
      !!                   
      !! ** Purpose :   compute trend for tracer upwind advection,
      !!                simple version of routine tra_adv_cen2
      !!
      !!                act on scalar 3D fields defined in T-points
      !!
      !!                B.C. : Zero through-boundary flux (if and only if no across-boundary flow)
      !!                       Zero through-bottom flux (as w = 0)
      !!                       Linear free surface B.C. or no-flux across surface B.C.
      !!                          both may be potentially unstable in absence of 
      !!                           vertical diffusion
      !!
      !! ** Method  :   rhsx = d_i(e2u e3u u t_u)
      !!                rhsy = d_j(e1v e3v v t_v)
      !!                rhsz = d_k(e1t e2t w t_w)
      !!                bt   = e1t * e2t * e3t
      !!                rhs  = - 1/bt * (rhsx + rhsy + rhsz)
      !!                   
      !!                It works in s-coordinates 
      !!                if w is dia-surface velocity (see NEMO-book)
      !!----------------------------------------------------------------------
      
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  ptn           ! tracer in T-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  pun, pvn, pwn ! 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  rhsn          ! advection tendency in T-points
      LOGICAL,                          INTENT(in   ) ::  bc_type       ! T: free-surface B.C., F: no-flux B.C.
      !
      INTEGER  :: ji, jj, jk                         ! dummy loop indices
      REAL(wp) :: zbtr
      REAL(wp) :: zun, zvn, zwn                      ! "velocity"
      REAL(wp) :: unp, unm, vnp, vnm, wnp, wnm       ! positive and negative parts of "velocity"
      REAL(wp), DIMENSION(jpi,jpj)     :: zwx, zwy   ! flux of tracer in u,v-points
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwz        ! flux of tracer in w-points
      !!----------------------------------------------------------------------
      
      !   ---------------------
      ! I. Horizontal advection
      !   ---------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   
               zun = e2u(ji,jj) * e3u_0(ji,jj,jk) * pun(ji,jj,jk)
               zvn = e1v(ji,jj) * e3v_0(ji,jj,jk) * pvn(ji,jj,jk)
     
               unp = max(zun,0._wp)
               unm = zun - unp
               vnp = max(zvn,0._wp)
               vnm = zvn - vnp

               zwx(ji,jj) = unp * ptn(ji,jj,jk) + unm * ptn(ji+1,jj,jk)
               zwy(ji,jj) = vnp * ptn(ji,jj,jk) + vnm * ptn(ji,jj+1,jk)
            END DO
         END DO
         
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
                  zbtr = r1_e12t(ji,jj) / e3t_0(ji,jj,jk)
                  rhsn(ji,jj,jk) = - zbtr * (  zwx(ji,jj) - zwx(ji-1,jj  )   &
                                             + zwy(ji,jj) - zwy(ji  ,jj-1) )
            END DO
         END DO
      END DO

      !   ---------------------
      ! II. Vertical advection
      !   ---------------------

      IF (bc_type) THEN
         zwz(:,:,1) = e12t(:,:) * pwn(:,:,1) * ptn(:,:,1) ! see traadv_cen2.F90
      ELSE
         zwz(:,:,1) = 0._wp ! zero cross-surface flux   
      END IF
      
      DO jk = 2, jpk
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zwn = e12t(ji,jj) * pwn(ji,jj,jk)

               wnp = max(zwn, 0._wp)
               wnm = zwn - wnp

               zwz(ji,jj,jk) = wnp * ptn(ji,jj,jk) + wnm * ptn(ji,jj,jk-1)
            END DO
         END DO    
      END DO
      
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
                  zbtr = r1_e12t(ji,jj) / e3t_0(ji,jj,jk)
                  rhsn(ji,jj,jk) = rhsn(ji,jj,jk) - zbtr * ( zwz(ji,jj,jk) - zwz(ji,jj,jk+1) )
            END DO
         END DO    
      END DO
      
   END SUBROUTINE upwind_advection
   
   SUBROUTINE laplace_T3D( ptn, rhsn, coef, bc_type ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE laplace_T3D  ***
      !!                   
      !! ** Purpose :   horizontal Laplacian diffusion operator, 
      !!                similar to routine tra_ldf_lap
      !!                
      !!                act on scalar 3D fields defined in T-points
      !!
      !!                B.C. : Zero flux or Zero Dirichlet
      !!
      !! ** Method  :   Horizontal diffusive trend of is given by:
      !!                difft = 1/(e1*e2*e3) {  di-1[ nu e2*e3/e1 di(tb) ]
      !!                                         + dj-1[ nu e1*e3/e2 dj(tb) ] }
      !!----------------------------------------------------------------------
      
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   ptn        ! tracer in T-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   rhsn       ! trend  in T-points
      REAL(wp)                                        ::   coef       ! viscosity coeff
      LOGICAL,                          INTENT(in   ) ::   bc_type    ! lateral B.C.: T: Dirichlet, F: no-flux
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zabe1, zabe2, zbtr   ! local scalars
      REAL(wp), DIMENSION (jpi,jpj) :: zwx, zwy ! fluxes of tracer in u,v-points
      REAL(wp), DIMENSION (jpi,jpj) :: ptn_m    ! masked array to apply zero Dirichlet B.C.
      !!----------------------------------------------------------------------

      DO jk=1,jpkm1
         !   ---------------------
         ! I. First derivative (gradient, scaled by divergence metric terms)
         !   ---------------------
         IF (bc_type) THEN
            ptn_m(:,:) = ptn(:,:,jk) * tmask(:,:,jk)
         ELSE
            ptn_m(:,:) = ptn(:,:,jk)
         END IF

         DO jj = 1, jpjm1
            ! umask, vmask make zero cross-boundary fluxes
            DO ji = 1, jpim1   
               zabe1 = coef * re2u_e1u(ji,jj) * e3u_0(ji,jj,jk)
               zabe2 = coef * re1v_e2v(ji,jj) * e3v_0(ji,jj,jk) 
               zwx(ji,jj) = zabe1 * ( ptn_m(ji+1,jj  ) - ptn_m(ji,jj) )
               zwy(ji,jj) = zabe2 * ( ptn_m(ji  ,jj+1) - ptn_m(ji,jj) )
            END DO
         END DO         
         
         IF (.not. bc_type) THEN
            zwx(:,:) = zwx(:,:) * umask(:,:,jk)
            zwy(:,:) = zwy(:,:) * vmask(:,:,jk)
         END IF

         !   ---------------------
         ! II. Second derivative (divergence)
         !   ---------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
               zbtr = r1_e12t(ji,jj) / e3t_0(ji,jj,jk)
               rhsn(ji,jj,jk) = zbtr * (zwx(ji,jj) - zwx(ji-1,jj) + zwy(ji,jj) - zwy(ji,jj-1))
            END DO
         END DO
      END DO
      
   END SUBROUTINE laplace_T3D
   
   SUBROUTINE filter_laplace_T3D( ptn ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE filter_laplace_T3D  ***
      !!                   
      !! ** Purpose :   filter based on laplace_T3D
      !!                
      !!                act on scalar 3D fields defined in T-points
      !!                 
      !!                B.C. : Zero flux
      !!
      !! ** Method  :   Horizontal diffusive trend of is given by:
      !!                difft = 1/(e1*e2*e3) {  di-1[ nu e2*e3/e1 di(tb) ]
      !!                                         + dj-1[ nu e1*e3/e2 dj(tb) ] }
      !!
      !!                This trend is used to construct filter:
      !!                ptn = ptn + difft
      !!                
      !!                anisotropic filter with diffusivity nu * dt = e^2/8 for each direction
      !!                nullifies chess-mode
      !!----------------------------------------------------------------------
      
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   ptn        ! before and after filtering in T-points
      !
      INTEGER  ::   ji, jj, jk                         ! dummy loop indices
      REAL(wp) ::   zabe1, zabe2, zbtr   ! local scalars
      REAL(wp),  DIMENSION (jpi,jpj) :: zwx, zwy       ! fluxes of tracer in u,v-points
      !!----------------------------------------------------------------------
      
      DO jk=1,jpkm1
         !   ---------------------
         ! I. First derivative (gradient, scaled by divergence metric terms)
         !   ---------------------
         DO jj = 1, jpjm1
            ! umask, vmask make zero cross-boundary fluxes
            DO ji = 1, jpim1   
               zabe1 = 0.125_wp * e12u(ji,jj) * e3u_0(ji,jj,jk)
               zabe2 = 0.125_wp * e12v(ji,jj) * e3v_0(ji,jj,jk)
               zwx(ji,jj) = zabe1 * ( ptn(ji+1,jj  ,jk) - ptn(ji,jj,jk) ) * umask(ji,jj,jk)
               zwy(ji,jj) = zabe2 * ( ptn(ji  ,jj+1,jk) - ptn(ji,jj,jk) ) * vmask(ji,jj,jk)
            END DO
         END DO         
         
         !   ---------------------
         ! II. Second derivative (divergence)
         !   ---------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
               zbtr = r1_e12t(ji,jj) / e3t_0(ji,jj,jk)
               ptn(ji,jj,jk) = ptn(ji,jj,jk) + zbtr * (zwx(ji,jj) - zwx(ji-1,jj) + zwy(ji,jj) - zwy(ji,jj-1))
            END DO
         END DO
      END DO
      
   END SUBROUTINE filter_laplace_T3D

   SUBROUTINE filter_laplace_T3D_dirichlet( ptn ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE filter_laplace_T3D_dirichlet  ***
      !!                   
      !! ** Purpose :   filter based on laplace_T3D
      !!                
      !!                act on scalar 3D fields defined in T-points
      !!                 
      !!                B.C. : Zero Dirichlet
      !!
      !! ** Method  :   Horizontal diffusive trend of is given by:
      !!                difft = 1/(e1*e2*e3) {  di-1[ nu e2*e3/e1 di(tb) ]
      !!                                         + dj-1[ nu e1*e3/e2 dj(tb) ] }
      !!
      !!                This trend is used to construct filter:
      !!                ptn = ptn + difft
      !!                
      !!                anisotropic filter with diffusivity nu * dt = e^2/8 for each direction
      !!                nullifies chess-mode
      !!----------------------------------------------------------------------
      
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   ptn        ! before and after filtering in T-points
      !
      INTEGER  ::   ji, jj, jk                         ! dummy loop indices
      REAL(wp) ::   zabe1, zabe2, zbtr   ! local scalars
      REAL(wp),  DIMENSION (jpi,jpj) :: zwx, zwy       ! fluxes of tracer in u,v-points
      REAL(wp),  DIMENSION (jpi,jpj) :: ptn_m          ! masked array
      !!----------------------------------------------------------------------
      
      DO jk=1,jpkm1
         !   ---------------------
         ! I. First derivative (gradient, scaled by divergence metric terms)
         !   ---------------------
         ptn_m(:,:) = ptn(:,:,jk) * tmask(:,:,jk)
         DO jj = 1, jpjm1
            ! umask, vmask make zero cross-boundary fluxes
            DO ji = 1, jpim1   
               zabe1 = 0.125_wp * e12u(ji,jj) * e3u_0(ji,jj,jk)
               zabe2 = 0.125_wp * e12v(ji,jj) * e3v_0(ji,jj,jk)
               zwx(ji,jj) = zabe1 * ( ptn_m(ji+1,jj  ) - ptn_m(ji,jj) )
               zwy(ji,jj) = zabe2 * ( ptn_m(ji  ,jj+1) - ptn_m(ji,jj) )
            END DO
         END DO         
         
         !   ---------------------
         ! II. Second derivative (divergence)
         !   ---------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
               zbtr = r1_e12t(ji,jj) / e3t_0(ji,jj,jk)
               ptn(ji,jj,jk) = ptn(ji,jj,jk) + zbtr * (zwx(ji,jj) - zwx(ji-1,jj) + zwy(ji,jj) - zwy(ji,jj-1))
            END DO
         END DO
      END DO
      
   END SUBROUTINE filter_laplace_T3D_dirichlet
   
   SUBROUTINE filter_laplace_T3D_ntimes( ptn, pta, ntimes, bc_type ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE filter_laplace_T3D_ntimes  ***
      !!                   
      !! ** Purpose :   apply n times filter_laplace_T3D
      !!                ptn must be mpi-exchanged
      !!                pta guaranteed to be mpi-exchanged
      !!----------------------------------------------------------------------
      
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   ptn        ! before filtering in T-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta        ! after  filtering in T-points
      INTEGER,                          INTENT(in   ) ::   ntimes     ! the number of filter applications
      LOGICAL,                          INTENT(in   ) ::   bc_type    ! lateral B.C.: T: Dirichlet, F: no-flux
      
      INTEGER  :: jtimes                          ! dummy loop
      REAL(wp), DIMENSION (jpi,jpj,jpk) :: ptim   ! intermediate array

      IF (ntimes == 0) THEN
         pta = ptn
         return
      END IF
      
      ptim = ptn
      
      DO jtimes = 1, ntimes
         IF (bc_type) THEN
            CALL filter_laplace_T3D_dirichlet(ptim)
         ELSE 
            CALL filter_laplace_T3D(ptim)
         END IF
         CALL lbc_lnk(ptim, 'T', 1.)
      END DO
      
      pta = ptim
      
   END SUBROUTINE filter_laplace_T3D_ntimes
    
   SUBROUTINE filter_laplace_f3D( ptn, ffmask ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE filter_laplace_f3D  ***
      !!                   
      !! ** Purpose :   filters out small horizontal spatial scales with Laplacian 
      !!                diffusion operator, similar to routine tra_ldf_lap
      !!                
      !!                act on scalar 3D fields defined in f-points
      !!
      !!                B.C. : Dirichlet zero
      !!
      !! ** Method  :   Horizontal diffusive trend of is given by:
      !!                difft = 1/(e1*e2*e3) {  di-1[ nu e2*e3/e1 di(tb) ]
      !!                                         + dj-1[ nu e1*e3/e2 dj(tb) ] }
      !!
      !!                This trend is used to construct filter:
      !!                ptn = ptn + difft
      !!                
      !!                anisotropic filter with diffusivity nu * dt = e^2/8 for each direction
      !!                nullifies chess-mode
      !!----------------------------------------------------------------------
      
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   ptn        ! before and after filtering in f-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   ffmask     ! fmask without rn_shlat
      !
      INTEGER  ::   ji, jj, jk                         ! dummy loop indices
      REAL(wp) ::   zabe1, zabe2, zbtr   ! local scalars
      REAL(wp), DIMENSION (jpi,jpj)     :: zwx, zwy    ! fluxes of tracer in v,u-points (since scale in f-points)
      !!----------------------------------------------------------------------
      
      DO jk=1,jpkm1
         !   ---------------------
         ! I. First derivative (gradient, scaled by divergence metric terms)
         !   ---------------------
         ! set zeros at the land and boundary
         ptn(:,:,jk) = ptn(:,:,jk) * ffmask(:,:,jk)
         DO jj = 1, jpjm1
            ! ji+1, jj+1 shifts are due to position of fluxes points relative to f-point
            DO ji = 1, jpim1   
               zabe1 = 0.125_wp * e12v(ji+1,jj) * e3v_0(ji+1,jj,jk)
               zabe2 = 0.125_wp * e12u(ji,jj+1) * e3u_0(ji,jj+1,jk)
               zwx(ji,jj) = zabe1 * ( ptn(ji+1,jj  ,jk) - ptn(ji,jj,jk) )
               zwy(ji,jj) = zabe2 * ( ptn(ji  ,jj+1,jk) - ptn(ji,jj,jk) )
            END DO
         END DO         
         
         !   ---------------------
         ! II. Second derivative (divergence)
         !   ---------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
               zbtr = r1_e12f(ji,jj) / e3f_0(ji,jj,jk)
               ptn(ji,jj,jk) = ptn(ji,jj,jk) + zbtr * (zwx(ji,jj) - zwx(ji-1,jj) + zwy(ji,jj) - zwy(ji,jj-1))
            END DO
         END DO
      END DO
      
   END SUBROUTINE filter_laplace_f3D
   
   SUBROUTINE filter_laplace_f3D_ntimes( ptn, pta, ntimes, ffmask ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE filter_laplace_f3D_ntimes  ***
      !!                   
      !! ** Purpose :   apply n times filter_laplace_f3D
      !!                ptn must be mpi-exchanged
      !!                pta guaranteed to be mpi-exchanged
      !!----------------------------------------------------------------------
      
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   ptn        ! before filtering in f-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta        ! after  filtering in f-points
      INTEGER, INTENT(in) ::   ntimes                                 ! the number of filter applications
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   ffmask     ! fmask without rn_shlat
      
      INTEGER  :: jtimes                          ! dummy loop
      REAL(wp), DIMENSION (jpi,jpj,jpk) :: ptim   ! intermediate array

      IF (ntimes == 0) THEN
         pta = ptn
         return
      END IF
      
      ptim = ptn
      
      DO jtimes = 1, ntimes
         CALL filter_laplace_f3D(ptim, ffmask)
         CALL lbc_lnk(ptim, 'F', 1.)
      END DO
      
      pta = ptim
      
   END SUBROUTINE filter_laplace_f3D_ntimes

   SUBROUTINE KEB_ldf_lap( prot, nu2t, nback, Eback, pua, pva, ffmask ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE KEB_ldf_lap  ***
      !!                  
      !! ** Purpose 1:  Compute vector-invariant laplacian trend of momentum
      !!                B.C. : zero vorticity on the boundary
      !!                analog of dynldf_lap
      !!                array with negative viscosity nu2t assumed to contain
      !!                abs value
      !!
      !! ** Purpose 2:  Subgrid energy loss by
      !!                laplacian backscatter:
      !!                Eback = |nu2| * (roth^2) > 0
      !!                if nback filter applied to negvisc tendency,
      !!                roth in this formula must be filtered nback/2 times
      !!                because
      !!                int (u * F_nback(u), dxdydz) = int (F_nback/2(u) * F_nback/2(u), dxdydz)
      !!
      !! ** Method  2:  Expression for flux is taken from rountine
      !!                dynldf_lap.f90
      !!                This expression is exact analytically, but
      !!                numerically due to finite-difference and boundary
      !!                issues has the following discrepancy:
      !!                Eback accuracy is 1% compared to direct computation
      !!
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  prot             ! rot of velocity field
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  pua, pva         ! trends of velocity
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  nu2t             ! negative vescosity coeff in T points
      
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  Eback            ! Energy backscatter, in T points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  ffmask           ! fmask without rn_shlat
      INTEGER,                          INTENT(in   ) ::  nback            ! number of filters for tendency
      
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices

      REAL(wp), DIMENSION(jpi,jpj,jpk) :: prot_f  ! filtered vorticity
      REAL(wp), DIMENSION(jpi,jpj)     :: Eback_f ! Eback in f points
      REAL(wp), DIMENSION(jpi,jpj)     :: zuf     

      !!----------------------------------------------------------------------

      CALL filter_laplace_f3D_ntimes( prot, prot_f, nback / 2, ffmask)
      
      ! compute Eback as in dynldf_lap
      DO jk = 1, jpkm1
         Eback_f(:,:) = 0.25_wp * prot_f(:,:,jk)**2 * e12f(:,:) * e3f_0(:,:,jk) * ffmask(:,:,jk)

         ! no need for MPI exchange
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               Eback(ji,jj,jk) =  nu2t(ji,jj,jk) * (                                                         &
                                  Eback_f(ji,jj) + Eback_f(ji-1,jj) + Eback_f(ji,jj-1) + Eback_f(ji-1,jj-1)  &
                                                   ) * r1_e12t(ji,jj) / e3t_0(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL filter_laplace_f3D_ntimes( prot_f, prot_f, nback / 2, ffmask) 

      ! Laplacian operator      
      DO jk = 1, jpkm1                               
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zuf(ji,jj) = - 0.25_wp * (nu2t(ji,jj,jk) + nu2t(ji+1,jj,jk) + nu2t(ji+1,jj+1,jk) + nu2t(ji,jj+1,jk)) *   &
                           prot_f(ji,jj,jk) * e3f_0(ji,jj,jk) * ffmask(ji,jj,jk)
            END DO
         END DO    

         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               pua(ji,jj,jk) = - ( zuf(ji,jj) - zuf(ji,jj-1) ) * r1_e2u(ji,jj) / e3u_0(ji,jj,jk)

               pva(ji,jj,jk) = + ( zuf(ji,jj) - zuf(ji-1,jj) ) * r1_e1v(ji,jj) / e3v_0(ji,jj,jk)
            END DO
         END DO
      END DO

   END SUBROUTINE KEB_ldf_lap

   SUBROUTINE horizontal_curl( psi, fx, fy )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE horizontal_curl  ***
      !!                   
      !! ** Purpose :   Compute curl of vertical vector field, i.e.
      !!                - rot_h ([0, 0, psi]).
      !!                Formulas are analogous to vorticity part of 
      !!                lateral viscosity operator, see dyn_ldf_lap.f90
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  psi              ! vertical component of vector field, in f-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  fx, fy           ! curl-compnents in u- and v-points
      
      INTEGER  ::   ji, jj, jk           ! dummy loop indices

      REAL(wp), DIMENSION(jpi,jpj) :: zuf
      !!----------------------------------------------------------------------
      
      DO jk = 1, jpkm1                
         zuf(:,:) = psi(:,:,jk) * e3f_0(:,:,jk)               
         DO jj = 2, jpj
            DO ji = 2, jpi
               fx(ji,jj,jk) = - ( zuf(ji,jj) - zuf(ji,jj-1) ) * r1_e2u(ji,jj) / e3u_0(ji,jj,jk)
               fy(ji,jj,jk) = + ( zuf(ji,jj) - zuf(ji-1,jj) ) * r1_e1v(ji,jj) / e3v_0(ji,jj,jk)
            END DO
         END DO
      END DO
      
   END SUBROUTINE horizontal_curl

   SUBROUTINE horizontal_divergence( fx, fy, hdiv )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE horizontal_divergence  ***
      !!                   
      !! ** Purpose :   Compute horizontal divergence for vector filed
      !!                Analogous to divcur.F90
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  fx, fy           ! vector field in u and v points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  hdiv             ! horizontal divergence
      
      INTEGER  ::   ji, jj, jk           ! dummy loop indices

      REAL(wp), DIMENSION(jpi,jpj) ::   zu, zv
      !!----------------------------------------------------------------------
      
      DO jk = 1, jpkm1                   
         zu(:,:) = fx(:,:,jk) * e3u_0(:,:,jk) * e2u(:,:)
         zv(:,:) = fy(:,:,jk) * e3v_0(:,:,jk) * e1v(:,:)
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               hdiv(ji,jj,jk) = (zu(ji,jj) - zu(ji-1,jj) + zv(ji,jj) - zv(ji,jj-1)) / (e12t(ji,jj) * e3t_0(ji,jj,jk))
            END DO
         END DO
      END DO
   
   END SUBROUTINE horizontal_divergence

   SUBROUTINE compute_rotor( fx, fy, hrot )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE compute_rotor  ***
      !!                   
      !! ** Purpose :   Compute rotor for vector filed
      !!                Analogous to divcur.F90
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  fx, fy           ! vector field in u and v points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  hrot             ! horizontal rotor
      
      INTEGER  ::   ji, jj, jk           ! dummy loop indices

      REAL(wp), DIMENSION(jpi,jpj) ::   zu, zv
      !!----------------------------------------------------------------------
      
      DO jk = 1, jpkm1                   
         zu(:,:) = fx(:,:,jk) * e1u(:,:)
         zv(:,:) = fy(:,:,jk) * e2v(:,:)
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               hrot(ji,jj,jk) = (zv(ji+1,jj) - zv(ji,jj) - zu(ji,jj+1) + zu(ji,jj)) / e12f(ji,jj)
            END DO
         END DO
      END DO
      
   END SUBROUTINE compute_rotor

   SUBROUTINE deformation_rate(pu, pv, D2, ffmask)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE deformation_rate  ***
      !!                   
      !! ** Purpose :   computes deformation rate D^2 = (du/dx-dv/dy)^2+(dv/dx+du/dy)^2
      !!                in T points. Analogous to routine ldf_dyn_smag
      !!                in curvilinear coordinates (see Griffies 2000):
      !!                du/dx ==> e2 d(u/e2)/dx;  du/dy ==> e1 d(u/e1)/dy; 
      !!                dv/dx ==> e2 d(v/e2)/dx;  dv/dy ==> e1 d(v/e1)/dy
      !!                interpolation from F-points to T-points is not conservative
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  pu, pv ! velocity
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  ffmask ! fmask without rn_shlat
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  D2     ! squared deformation rate in T points

      INTEGER  ::   ji, jj, jk           ! dummy loop indices

      REAL (wp), DIMENSION (jpi,jpj) ::   zux, zuy, zvx, zvy, zue1, zue2, zve1, zve2, D2_f

      D2 = 0._wp
      DO jk = 1, jpkm1
         zue2(:,:)=pu(:,:,jk) * r1_e2u(:,:)
         zve1(:,:)=pv(:,:,jk) * r1_e1v(:,:)
         zue1(:,:)=pu(:,:,jk) * r1_e1u(:,:)
         zve2(:,:)=pv(:,:,jk) * r1_e2v(:,:)

         DO jj=2,jpjm1
            DO ji=2,jpim1
               zux(ji,jj) = (zue2(ji,jj)-zue2(ji-1,jj)) * r1_e1t(ji,jj) * e2t(ji,jj)
               zvy(ji,jj) = (zve1(ji,jj)-zve1(ji,jj-1)) * r1_e2t(ji,jj) * e1t(ji,jj)
            END DO
         END DO

         DO jj=1,jpjm1
            DO ji=1,jpim1
               zvx(ji,jj) = (zve2(ji+1,jj)-zve2(ji,jj)) * r1_e1f(ji,jj) * e2f(ji,jj) * ffmask(ji,jj,jk)
               zuy(ji,jj) = (zue1(ji,jj+1)-zue1(ji,jj)) * r1_e2f(ji,jj) * e1f(ji,jj) * ffmask(ji,jj,jk)
            END DO
         END DO

         D2_f(:,:) = (zvx(:,:) + zuy(:,:)) ** 2

         DO jj=2,jpjm1
            DO ji=2,jpim1
               D2(ji,jj,jk) = tmask(ji,jj,jk) * ((zux(ji,jj) - zvy(ji,jj)) ** 2 + &
                              (D2_f(ji,jj) + D2_f(ji-1,jj) + D2_f(ji-1,jj-1) + D2_f(ji,jj-1)) * 0.25_wp)
            END DO
         END DO
      END DO

   END SUBROUTINE deformation_rate

   SUBROUTINE estimate_TKE( pu, pv, prot, TKE, ffmask )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE estimate_TKE  ***
      !!                   
      !! ** Purpose :   Estimates subgrid energy based on paper
      !!                Perezhogin, Glazunov "A priori and a posteriori analysis
      !!                in large eddy simulation of the twodimensional decaying turbulence
      !!                using an explicit filtering approach":
      !!                TKE = delta^2/48 * (D^2 + omega^2),
      !!                where delta - filter width square. Here we assume delta^2 = e1*e2*sqrt(6)^2
      !!                D^2 - squared deformation rate, omega - relative vorticity
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  pu, pv ! velocity
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  prot   ! relative vorticity
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  ffmask ! fmask without rn_shlat
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  TKE    ! estimated TKE
      
      INTEGER  :: ji, jj, jk           ! dummy loop indices
      REAL(wp) :: c0

      REAL(wp), DIMENSION(jpi,jpj)     :: omega2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: D2
      
      c0 = 6. / 48. ! 6 = squared relative filter width, 48 - constant from paper

      CALL deformation_rate(pu, pv, D2, ffmask)

      TKE = 0._wp
      DO jk = 1, jpkm1
         omega2(:,:) = prot(:,:,jk) ** 2 * ffmask(:,:,jk)

         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               TKE(ji,jj,jk) = c0 * e12t(ji,jj) * tmask(ji,jj,jk) * (D2(ji,jj,jk) + &
                               (omega2(ji,jj) + omega2(ji-1,jj) + omega2(ji-1,jj-1) + omega2(ji,jj-1)) * 0.25_wp)
            END DO
         END DO
      END DO
      
   END SUBROUTINE estimate_TKE

   SUBROUTINE Klower_cdiss( pu, pv, Ediss, R_diss, local_cdiss, ffmask )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Klower_cdiss  ***
      !!                   
      !! ** Purpose :   According to Klower 2018, Esource = Ediss * c_diss,
      !!                c_diss = (1 + R_local / R_diss) ^ -1
      !!                R_local = D / f, D - deformation rate, f - coriolis parameter
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  pu, pv       ! velocity
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  ffmask       ! fmask without rn_shlat
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  Ediss        ! gross dissipation
      REAL(wp),                         INTENT(in)    ::  R_diss       ! parameter in Klower formula
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  local_cdiss  ! local backscatter rate
      
      INTEGER  :: ji, jj, jk           ! dummy loop indices
      REAL(wp) :: R_local

      REAL(wp), DIMENSION(jpi,jpj,jpk) :: D2      ! squared deformation rate
      REAL(wp), DIMENSION(jpi,jpj)     :: r1_ff   ! inverse Coriolis parameter in T points

      CALL deformation_rate(pu, pv, D2, ffmask)

      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            r1_ff(ji,jj) = 1._wp / ((ff(ji,jj) + ff(ji-1,jj) + ff(ji-1,jj-1) + ff(ji,jj-1)) * 0.25_wp)
         END DO
      END DO

      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               R_local = sqrt(D2(ji,jj,jk)) * r1_ff(ji,jj)
               local_cdiss(ji,jj,jk) = R_diss / (R_local + R_diss)
            END DO
         END DO
      END DO

   END SUBROUTINE Klower_cdiss

   SUBROUTINE gauss_white_noise_2d( field2d )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE gauss_white_noise_2d  ***
      !!                   
      !! ** Purpose :   generate 2D white noise in space with distribution N(0,1)
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::  field2d

      INTEGER      :: ji, jj
      REAL(KIND=8) :: grand

      DO jj=1,jpj
         DO ji=1,jpi
            CALL kiss_gaussian(grand)
            field2d(ji,jj) = grand
         END DO
      END DO

   END SUBROUTINE gauss_white_noise_2d

   SUBROUTINE compute_psi( Esource, phi, nstoch, T_decorr, psi, ffmask )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE compute_psi  ***
      !!                   
      !! ** Purpose :   weight stochastic 2d noise (psi) with
      !!                energy source (Esource) and filter
      !!                resulting energy:
      !!                psi_ijk = hat(A*phi)_ij
      !!                where hat(phi) - filtered nstoch times phi
      !!
      !! ** Method  :   We assume that Esource can be 
      !!                taken out of all filtering and finite 
      !!                difference operators. Resulting formulas
      !!                analogous to O'Neil 2016 dissertation 
      !!                energy input equals:
      !!                dE/dt = T_decorr < fx^2 + fy^2 >,
      !!                where (fx,fy) = curl(psi)
      !!                Then
      !!                dE/dt = 2 * T_decorr * A^2 *
      !!                var(hat(phi)) * ((1-rx)/dx^2 + (1-ry)/dy^2)
      !!                where rx = corr(hat(phi)_ij, hat(phi)_i-1j)
      !!                      ry = corr(hat(phi)_ij, hat(phi)_ij-1)
      !!                For simple 5-point filter killing chess-mode
      !!                rx = ry = exp(-1/nstoch)
      !!                d2 = 1/var(hat(phi)) = nstoch*pi
      !!                Final formula:
      !!                A^2 = Esource * d2 / (2 * T_decorr)/
      !!                ((1-rx)/dx^2 + (1-ry)/dy^2)
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  Esource ! source of energy
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  ffmask  ! fmask without rn_shlat
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  psi     ! stochastic streamfunction
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(in   ) ::  phi     ! white-noise N(0,1)
      REAL(wp), INTENT(in) :: T_decorr
      INTEGER , INTENT(in) :: nstoch

      INTEGER  :: ji, jj, jk  ! dummy loop indices
      REAL(wp) :: rx, ry      ! correlation in dx and dy
      REAL(wp) :: d2          ! dispersion reduction
      REAL(wp) :: Es_f        ! Energy source interpolation

      REAL(wp), DIMENSION(jpi,jpj) :: c0_2d ! 2d amplitude factor

      rx = exp(-1._wp / nstoch)
      ry = exp(-1._wp / nstoch)
      d2 = nstoch * rpi

      ! account for dispersion reduction, correlation time and space correlation
      c0_2d(:,:) = d2 / (2._wp * T_decorr) / &
                  ((1._wp-rx) / (e1f(:,:) ** 2) + (1._wp-ry) / (e2f(:,:) ** 2))

      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ! T to F interpolation
               Es_f = 0.25_wp * (Esource(ji,jj,jk) + Esource(ji+1,jj,jk) + Esource(ji,jj+1,jk) + Esource(ji+1,jj+1,jk))
               psi(ji,jj,jk) = sqrt(Es_f * c0_2d(ji,jj)) * phi(ji,jj)
            END DO
         END DO
      END DO

      CALL lbc_lnk(psi, 'F', 1.)
      CALL filter_laplace_f3D_ntimes(psi, psi, nstoch, ffmask)

      psi = psi * ffmask ! will be needed for curl

   END SUBROUTINE compute_psi

   SUBROUTINE compute_Eback_AR1( fx, fy, pun, pvn, Eback )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE compute_Eback_AR1  ***
      !!                   
      !! ** Purpose :   conservative interpolation 
      !!                of energy generation to T-points
      !!  
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  fx, fy   ! velocity tendency
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  Eback    ! energy tendency in T points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::  pun, pvn ! velocity

      INTEGER :: ji, jj, jk

      REAL(wp), DIMENSION(jpi,jpj) :: Eback_u, Eback_v

      CALL lbc_lnk(fx, 'U', -1.)
      CALL lbc_lnk(fy, 'V', -1.)

      ! fields must be exchanged
      DO jk=1, jpkm1
         Eback_u(:,:) = 0.5_wp * fx(:,:,jk) * pun(:,:,jk) * e12u(:,:) * e3u_0(:,:,jk) * umask(:,:,jk)
         Eback_v(:,:) = 0.5_wp * fy(:,:,jk) * pvn(:,:,jk) * e12v(:,:) * e3v_0(:,:,jk) * vmask(:,:,jk)

         ! no need for MPI exchange
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               Eback(ji,jj,jk) = (Eback_u(ji,jj) + Eback_u(ji-1,jj) + Eback_v(ji,jj) + Eback_v(ji,jj-1)) * &
                                 r1_e12t(ji,jj) / e3t_0(ji,jj,jk)
            END DO
         END DO
      END DO

   END SUBROUTINE compute_Eback_AR1

   FUNCTION aposteriori_correction( Esource, Estoch,  Esource_mean, Estoch_mean ) result(amp_increase)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE aposteriori_correction  ***
      !!                   
      !! ** Purpose :   make energy input of stochastic 
      !!                noise equal to energy source (Esource)
      !!
      !! ** Method  :   accumulate mean values of Estoch and Esource,
      !!                and produce multiplicative factor
      !!                for stochastic noise:
      !!                amp_increase = mean(Esource) / mean(Estoch)
      !!  
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) ::  Esource ! desired energy input
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) ::  Estoch  ! energy input before correction
      REAL(wp), INTENT(inout) :: Esource_mean,  Estoch_mean
      REAL(wp) :: amp_increase

      REAL(wp) :: Estoch_m, Esource_m
      REAL(wp) :: T, beta, alpha

      T = 30._wp * 86400._wp ! time scale to average in sec
      beta = rdt / T
      alpha = 1._wp - beta

      Estoch_m  = average_xyz(Estoch, e1t, e2t, e3t_0, tmask)
      Esource_m = average_xyz(Esource, e1t, e2t, e3t_0, tmask)

      ! time accumulation of stochastic input and energy source
      Esource_mean = alpha * Esource_mean + beta * Esource_m
      Estoch_mean  = alpha * Estoch_mean  + beta * Estoch_m

      ! 1e-30 because energy exchange itself small, about 1e-9
      amp_increase = 1._wp
      IF (Estoch_mean > 1.e-30) amp_increase = Esource_mean / Estoch_mean

      ! energy input is proportional to amplitude squared,
      ! and only twice change is allowed
      amp_increase = max(min(amp_increase, sqrt(2._wp)), 1._wp/sqrt(2._wp))

   END FUNCTION aposteriori_correction

   SUBROUTINE integrate_z( field3d, field2d, e3, mask )
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  field3d, e3, mask
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(  out) ::  field2d
      
      INTEGER  :: jk
      
      field2d = 0._wp
      DO jk = 1, jpkm1
         field2d(:,:) = field2d(:,:) + field3d(:,:,jk) * e3(:,:,jk) * mask(:,:,jk)
      END DO
   END SUBROUTINE integrate_z

   SUBROUTINE average_z( field3d, field2d, e3, mask )
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  field3d, e3, mask
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(  out) ::  field2d
      
      INTEGER  :: jk

      REAL(wp), DIMENSION (jpi,jpj) :: depth2d
      
      field2d = 0._wp
      depth2d = 0._wp
      DO jk = 1, jpkm1
         field2d(:,:) = field2d(:,:) + field3d(:,:,jk) * e3(:,:,jk) * mask(:,:,jk)
         depth2d(:,:) = depth2d(:,:) + e3(:,:,jk) * mask(:,:,jk)
      END DO

      ! 1e-16 to avoid NANs
      ! typical values of depth is 10, so this 1e-16 is small enough
      field2d(:,:) = field2d(:,:) / (depth2d(:,:) + 1.e-16)

   END SUBROUTINE average_z

   SUBROUTINE average_xy( field3d, field1d, e1, e2, mask )
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  field3d, mask
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(in   ) ::  e1, e2
      REAL(wp), DIMENSION(jpk)        , INTENT(  out) ::  field1d
      
      INTEGER  :: jk
      REAL(wp) :: area2d(jpk), intgrl2d(jpk)
      REAL(wp) :: array2d(jpi,jpj)
      
      DO jk = 1, jpkm1
         array2d(:,:) = field3d(:,:,jk) * e1(:,:) * e2(:,:) * mask(:,:,jk)
         intgrl2d(jk) = SUM(array2d(2:nlci-1,2:nlcj-1))

         array2d(:,:) = e1(:,:) * e2(:,:) * mask(:,:,jk)
         area2d(jk) = SUM(array2d(2:nlci-1,2:nlcj-1))
      END DO      

      CALL mpp_sum(intgrl2d, jpk)
      CALL mpp_sum(area2d, jpk)

      field1d(:) = intgrl2d(:) / area2d(:)

   END SUBROUTINE average_xy
   
   FUNCTION average_xyz( field3d, e1, e2, e3, mask ) result(res)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  field3d, e3, mask
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(in   ) ::  e1, e2
      REAL(wp) :: res

      INTEGER  :: ji, jj, jk
      REAL(wp) :: volume, dV
      REAL(wp), DIMENSION(jpi,jpj) :: e12

      e12(:,:) = e1(:,:) * e2(:,:)
      volume = 0._wp
      res = 0._wp
      DO jk = 1, jpkm1
         DO jj = 2, nlcj-1
            DO ji = 2, nlci-1
               dV = e12(ji,jj) * e3(ji,jj,jk) * mask(ji,jj,jk)
               volume = volume + dV
               res = res + field3d(ji,jj,jk) * dV
            END DO
         END DO
      END DO

      CALL mpp_sum(res)
      CALL mpp_sum(volume)
      
      res = res / volume
      
   END FUNCTION average_xyz

   FUNCTION min_xyz( field3d, mask ) result(res)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  field3d, mask
      REAL(wp) :: res

      INTEGER  :: ji, jj, jk
      
      res = 1.e+34

      DO jk = 1, jpkm1                               
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               IF (mask(ji,jj,jk) > 0.5_wp) res = min(res, field3d(ji,jj,jk))               
            END DO
         END DO
      END DO

      CALL mpp_min(res)
      
   END FUNCTION min_xyz

   FUNCTION max_xyz( field3d, mask ) result(res)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  field3d, mask
      REAL(wp) :: res

      INTEGER  :: ji, jj, jk
      
      res = -1.e+34

      DO jk = 1, jpkm1                               
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               IF (mask(ji,jj,jk) > 0.5_wp) res = max(res, field3d(ji,jj,jk))               
            END DO
         END DO
      END DO

      CALL mpp_max(res)
      
   END FUNCTION max_xyz

   FUNCTION relative_nonconserv( field3d, e1, e2, e3, mask ) result(rerror)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  field3d, e3, mask
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(in   ) ::  e1, e2
      REAL(wp) ::  rerror

      rerror = average_xyz( field3d, e1, e2, e3, mask ) / max_xyz(abs(field3d), mask)
      
   END FUNCTION relative_nonconserv
END MODULE KEB_operators