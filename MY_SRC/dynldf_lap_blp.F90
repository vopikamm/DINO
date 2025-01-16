MODULE dynldf_lap_blp
   !!======================================================================
   !!                   ***  MODULE  dynldf_lap_blp  ***
   !! Ocean dynamics:  lateral viscosity trend (laplacian and bilaplacian)
   !!======================================================================
   !! History : 3.7  ! 2014-01  (G. Madec, S. Masson)  Original code, re-entrant laplacian
   !!           4.0  ! 2020-04  (A. Nasser, G. Madec)  Add symmetric mixing tensor 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf_lap   : update the momentum trend with the lateral viscosity using an iso-level   laplacian operator
   !!   dyn_ldf_blp   : update the momentum trend with the lateral viscosity using an iso-level bilaplacian operator
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE domutl, ONLY : is_tile
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE ldfslp         ! iso-neutral slopes 
   USE zdf_oce        ! ocean vertical physics
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp

   USE KEB_module      ! compute Ediss for backscatter
   !USE KEB_operators   ! for KEB_test

#if defined key_loop_fusion
   USE dynldf_lap_blp_lf
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_ldf_lap  ! called by dynldf.F90
   PUBLIC dyn_ldf_blp  ! called by dynldf.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynldf_lap_blp.F90 15033 2021-06-21 10:24:45Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf_lap( kt, Kbb, Kmm, pu, pv, pu_rhs, pv_rhs, kpass )
      !!
      INTEGER                   , INTENT(in   ) ::   kt               ! ocean time-step index
      INTEGER                   , INTENT(in   ) ::   Kbb, Kmm         ! ocean time level indices
      INTEGER                   , INTENT(in   ) ::   kpass            ! =1/2 first or second passage
      REAL(dp), DIMENSION(:,:,:), INTENT(in   ) ::   pu, pv           ! before velocity  [m/s]
      REAL(dp), DIMENSION(:,:,:), INTENT(inout) ::   pu_rhs, pv_rhs   ! velocity trend   [m/s2]
      !!
#if defined key_loop_fusion
      CALL dyn_ldf_lap_lf( kt, Kbb, Kmm, pu, pv, pu_rhs, pv_rhs, kpass )
#else
      CALL dyn_ldf_lap_t( kt, Kbb, Kmm, pu, pv, is_tile(pu), pu_rhs, pv_rhs, is_tile(pu_rhs), kpass )
#endif

   END SUBROUTINE dyn_ldf_lap


   SUBROUTINE dyn_ldf_lap_t( kt, Kbb, Kmm, pu, pv, ktuv, pu_rhs, pv_rhs, ktuv_rhs, kpass )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_lap  ***
      !!                       
      !! ** Purpose :   Compute the before horizontal momentum diffusive 
      !!      trend and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The Laplacian operator apply on horizontal velocity is 
      !!      writen as :   grad_h( ahmt div_h(U )) - curl_h( ahmf curl_z(U) ) 
      !!
      !! ** Action : - pu_rhs, pv_rhs increased by the harmonic operator applied on pu, pv.
      !!
      !! Reference : S.Griffies, R.Hallberg 2000 Mon.Wea.Rev., DOI:/ 
      !!----------------------------------------------------------------------
      INTEGER                                 , INTENT(in   ) ::   kt               ! ocean time-step index
      INTEGER                                 , INTENT(in   ) ::   Kbb, Kmm         ! ocean time level indices
      INTEGER                                 , INTENT(in   ) ::   kpass            ! =1/2 first or second passage
      INTEGER                                 , INTENT(in   ) ::   ktuv, ktuv_rhs
      REAL(dp), DIMENSION(A2D_T(ktuv)    ,JPK), INTENT(in   ) ::   pu, pv           ! before velocity  [m/s]
      REAL(dp), DIMENSION(A2D_T(ktuv_rhs),JPK), INTENT(inout) ::   pu_rhs, pv_rhs   ! velocity trend   [m/s2]
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   iij
      REAL(wp) ::   zsign        ! local scalars
      REAL(wp) ::   zua, zva     ! local scalars
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zcur, zdiv
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zten, zshe   ! tension (diagonal) and shearing (anti-diagonal) terms
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 .AND. lwp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'dyn_ldf : iso-level harmonic (laplacian) operator, pass=', kpass
            WRITE(numout,*) '~~~~~~~ '
         ENDIF
      ENDIF
      !
      ! Define pu_rhs/pv_rhs halo points for multi-point haloes in bilaplacian case
      IF( nldf_dyn == np_blp .AND. kpass == 1 ) THEN ; iij = nn_hls
      ELSE                                           ; iij = 1
      ENDIF
      !
      IF( kpass == 1 ) THEN   ;   zsign =  1._wp      ! bilaplacian operator require a minus sign
      ELSE                    ;   zsign = -1._wp      !  (eddy viscosity coef. >0)
      ENDIF
      !
      SELECT CASE( nn_dynldf_typ )  
      !              
      CASE ( np_typ_rot )       !==  Vorticity-Divergence operator  ==!
         !
         ALLOCATE( zcur(A2D(nn_hls)) , zdiv(A2D(nn_hls)) )
         !
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !
            DO_2D( iij-1, iij, iij-1, iij )
               !                                      ! ahm * e3 * curl  (warning: computed for ji-1,jj-1)
               zcur(ji-1,jj-1) = ahmf(ji-1,jj-1,jk) * e3f(ji-1,jj-1,jk) * r1_e1e2f(ji-1,jj-1)       &   ! ahmf already * by fmask
                  &     * (  e2v(ji  ,jj-1) * pv(ji  ,jj-1,jk) - e2v(ji-1,jj-1) * pv(ji-1,jj-1,jk)  &
                  &        - e1u(ji-1,jj  ) * pu(ji-1,jj  ,jk) + e1u(ji-1,jj-1) * pu(ji-1,jj-1,jk)  )
               !                                      ! ahm * div        (warning: computed for ji,jj)
               zdiv(ji,jj)     = ahmt(ji,jj,jk) * r1_e1e2t(ji,jj) / e3t(ji,jj,jk,Kbb)               &   ! ahmt already * by tmask
                  &     * (  e2u(ji,jj)*e3u(ji,jj,jk,Kbb) * pu(ji,jj,jk) - e2u(ji-1,jj)*e3u(ji-1,jj,jk,Kbb) * pu(ji-1,jj,jk)  &
                  &        + e1v(ji,jj)*e3v(ji,jj,jk,Kbb) * pv(ji,jj,jk) - e1v(ji,jj-1)*e3v(ji,jj-1,jk,Kbb) * pv(ji,jj-1,jk)  )
            END_2D
            !
            DO_2D( iij-1, iij-1, iij-1, iij-1 )   ! - curl( curl) + grad( div )
               pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign * umask(ji,jj,jk) * (    &    ! * by umask is mandatory for dyn_ldf_blp use
                  &              - ( zcur(ji  ,jj) - zcur(ji,jj-1) ) * r1_e2u(ji,jj) / e3u(ji,jj,jk,Kmm)   &
                  &              + ( zdiv(ji+1,jj) - zdiv(ji,jj  ) ) * r1_e1u(ji,jj)                      )
               !
               pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign * vmask(ji,jj,jk) * (    &    ! * by vmask is mandatory for dyn_ldf_blp use
                  &                ( zcur(ji,jj  ) - zcur(ji-1,jj) ) * r1_e1v(ji,jj) / e3v(ji,jj,jk,Kmm)   &
                  &              + ( zdiv(ji,jj+1) - zdiv(ji  ,jj) ) * r1_e2v(ji,jj)                      )
            END_2D
            !
         END DO                                           !   End of slab
         !
         DEALLOCATE( zcur , zdiv )
         !
      CASE ( np_typ_sym )       !==  Symmetric operator  ==!
         !
         ALLOCATE( zten(A2D(nn_hls)) , zshe(A2D(nn_hls)) )
         !
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !
            DO_2D( iij-1, iij, iij-1, iij )
               !                                      ! shearing stress component (F-point)   NB : ahmf has already been multiplied by fmask
               zshe(ji-1,jj-1) = ahmf(ji-1,jj-1,jk)                                                              &
                  &     * (    e1f(ji-1,jj-1)    * r1_e2f(ji-1,jj-1)                                             &
                  &         * ( pu(ji-1,jj  ,jk) * r1_e1u(ji-1,jj  )  - pu(ji-1,jj-1,jk) * r1_e1u(ji-1,jj-1) )   &
                  &         +  e2f(ji-1,jj-1)    * r1_e1f(ji-1,jj-1)                                             &
                  &         * ( pv(ji  ,jj-1,jk) * r1_e2v(ji  ,jj-1)  - pv(ji-1,jj-1,jk) * r1_e2v(ji-1,jj-1) )   ) 
               !                                      ! tension stress component (T-point)   NB : ahmt has already been multiplied by tmask
               zten(ji,jj)    = ahmt(ji,jj,jk)                                                       &
                  &     * (    e2t(ji,jj)    * r1_e1t(ji,jj)                                         &
                  &         * ( pu(ji,jj,jk) * r1_e2u(ji,jj)  - pu(ji-1,jj,jk) * r1_e2u(ji-1,jj) )   &
                  &         -  e1t(ji,jj)    * r1_e2t(ji,jj)                                         &
                  &         * ( pv(ji,jj,jk) * r1_e1v(ji,jj)  - pv(ji,jj-1,jk) * r1_e1v(ji,jj-1) )   )   
            END_2D
            !
            DO_2D( iij-1, iij-1, iij-1, iij-1 )
               pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)                               &
                  &    * (   (   zten(ji+1,jj  ) * e2t(ji+1,jj  )*e2t(ji+1,jj  ) * e3t(ji+1,jj  ,jk,Kmm)                       &
                  &            - zten(ji  ,jj  ) * e2t(ji  ,jj  )*e2t(ji  ,jj  ) * e3t(ji  ,jj  ,jk,Kmm) ) * r1_e2u(ji,jj)     &                                                    
                  &        + (   zshe(ji  ,jj  ) * e1f(ji  ,jj  )*e1f(ji  ,jj  ) * e3f(ji  ,jj  ,jk)                           &
                  &            - zshe(ji  ,jj-1) * e1f(ji  ,jj-1)*e1f(ji  ,jj-1) * e3f(ji  ,jj-1,jk)     ) * r1_e1u(ji,jj) )   
               !
               pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)                               &
                  &    * (   (   zshe(ji  ,jj  ) * e2f(ji  ,jj  )*e2f(ji  ,jj  ) * e3f(ji  ,jj  ,jk)                           &
                  &            - zshe(ji-1,jj  ) * e2f(ji-1,jj  )*e2f(ji-1,jj  ) * e3f(ji-1,jj  ,jk)     ) * r1_e2v(ji,jj)     &
                  &        - (   zten(ji  ,jj+1) * e1t(ji  ,jj+1)*e1t(ji  ,jj+1) * e3t(ji  ,jj+1,jk,Kmm)                       &
                  &            - zten(ji  ,jj  ) * e1t(ji  ,jj  )*e1t(ji  ,jj  ) * e3t(ji  ,jj  ,jk,Kmm) ) * r1_e1v(ji,jj) )
               !
            END_2D
            !
         END DO
         !
         DEALLOCATE( zten , zshe )
         !
      END SELECT
      !
   END SUBROUTINE dyn_ldf_lap_t


   SUBROUTINE dyn_ldf_blp( kt, Kbb, Kmm, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_ldf_blp  ***
      !!                    
      !! ** Purpose :   Compute the before lateral momentum viscous trend 
      !!              and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The lateral viscous trends is provided by a bilaplacian
      !!      operator applied to before field (forward in time).
      !!      It is computed by two successive calls to dyn_ldf_lap routine
      !!
      !! ** Action :   pt(:,:,:,:,Krhs)   updated with the before rotated bilaplacian diffusion
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kbb, Kmm   ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pu, pv     ! before velocity fields
      REAL(dp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs   ! momentum trend
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   Ediss_u, Ediss_v       ! 2D workspace for E_diss
      !
      REAL(dp), DIMENSION(A2D(nn_hls),jpk) ::   zulap, zvlap   ! laplacian at u- and v-point
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      
      !!----------------------------------------------------------------------
      !
#if defined key_loop_fusion
      CALL dyn_ldf_blp_lf( kt, Kbb, Kmm, pu, pv, pu_rhs, pv_rhs )
#else
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 )  THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_ldf_blp : bilaplacian operator momentum '
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
      zulap(:,:,:) = 0._wp
      zvlap(:,:,:) = 0._wp
      !
      CALL dyn_ldf_lap( kt, Kbb, Kmm, pu, pv, zulap, zvlap, 1 )   ! rotated laplacian applied to pt (output in zlap,Kbb)
      !
      IF (nn_hls==1) CALL lbc_lnk( 'dynldf_lap_blp', zulap, 'U', -1.0_dp, zvlap, 'V', -1.0_dp )             ! Lateral boundary conditions
      !
      ! Kinetic Energy Backscatter needs dissipated energy diagnosed here
      ! conservative interpolation of dissipation to T points
      ! masking (umask, vmask) allows to exclude dissipation on the boundary
      IF (KEB_on) THEN
         Ediss_u(:,:) = 0._wp
         Ediss_v(:,:) = 0._wp
         
         DO jk=1, jpk-1
            DO_2D(1,0,1,0)
               Ediss_u(ji,jj) = - 0.5_wp * zulap(ji,jj,jk)**2 * e1e2u(ji,jj) * e3u(ji,jj,jk,Kbb) * umask(ji,jj,jk)
               Ediss_v(ji,jj) = - 0.5_wp * zvlap(ji,jj,jk)**2 * e1e2v(ji,jj) * e3v(ji,jj,jk,Kbb) * vmask(ji,jj,jk)
            END_2D
            ! Ediss defined in KEB_module.F90; TODO: should this stay a global variable?
            ! no need for MPI exchange
            DO_2D(0,0,0,0)
               Ediss(ji,jj,jk) =                                                                    & 
               (Ediss_u(ji,jj) + Ediss_u(ji-1,jj) + Ediss_v(ji,jj) + Ediss_v(ji,jj-1))              &
               * r1_e1e2t(ji,jj) / e3t(ji,jj,jk, Kbb) ! dk: find out if coefficient is applied twice, then: / ahmt(ji,jj,jk)
            END_2D
         END DO
      END IF
      !
      CALL dyn_ldf_lap( kt, Kbb, Kmm, zulap, zvlap, pu_rhs, pv_rhs, 2 )   ! rotated laplacian applied to zlap (output in pt(:,:,:,:,Krhs))
      !
#endif
   END SUBROUTINE dyn_ldf_blp

   !!======================================================================
END MODULE dynldf_lap_blp
