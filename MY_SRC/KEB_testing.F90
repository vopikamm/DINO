MODULE KEB_testing
 
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE lib_mpp         ! MPP library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! lwp, numout
   USE KEB_operators   ! what we want to test
 
   IMPLICIT NONE   
   PUBLIC
    
 CONTAINS
    
   SUBROUTINE test_negvisc_KEB( puu, pvv, pww, TKE, rhs_adv, rhs_diff, Esource, Eback, Ediss, &
                                local_cdiss, ffmask, rhsu, rhsv, Ediss_check )
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  puu, pvv, pww ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    :: TKE, rhs_adv, rhs_diff
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    :: Esource, Eback, Ediss
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)    :: local_cdiss, ffmask
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) :: rhsu, rhsv ! inout to lbc_lnk
      REAL(wp),                         INTENT(in)    :: Ediss_check 
      
      INTEGER  :: jk
      REAL(wp) :: rerr
      REAL(wp) :: mean_1, mean_2, mean_3
      REAL(wp) :: min_val, max_val
      REAL(wp) :: m_TKE
      
      REAL(wp), DIMENSION (jpi,jpj,jpk) :: rhs
      
      IF (lwp) write(numout, *) ''
      IF (lwp) write(numout, *) '~~~~~~~~~~~ KEB negvisc testing ~~~~~~~~~~~'
      IF (lwp) write(numout, *) ''
  
      if (lwp) write(numout,*) '! --------------------------- advection ------------------------- !'
      
      if (lwp) write(numout,*) 'positiveness:'
      CALL check_positiveness(TKE, TKE+rhs_adv*rdt, tmask)
      
      IF (lwp) write(numout, *) ''
      
      if (lwp) write(numout,*) 'conservation:'
      CALL upwind_advection( TKE, puu, pvv, pww, rhs, .false. )
      rerr = relative_nonconserv( rhs, e1t, e2t, e3t_0, tmask)
      if (lwp) write(numout,*) 'no-flux            (1e-20 ok):', rerr
  
      CALL upwind_advection( TKE, puu, pvv, pww, rhs, .true. )
      rerr = relative_nonconserv( rhs, e1t, e2t, e3t_0, tmask)
      if (lwp) write(numout,*) 'linear free-surface (1e-8 ok):', rerr
      
      IF (lwp) write(numout, *) ''

      if (lwp) write(numout,*) 'CFL:'
      max_val = 0._wp

      DO jk = 1, jpkm1
         rhs(:,:,jk) = abs(puu(:,:,jk)) / e1u(:,:)
      END DO
      max_val = max(max_val, rdt * max_xyz(rhs, umask))

      DO jk = 1, jpkm1
         rhs(:,:,jk) = abs(pvv(:,:,jk)) / e2v(:,:)
      END DO
      max_val = max(max_val, rdt * max_xyz(rhs, vmask))

      rhs = abs(pww) / e3w_0
      max_val = max(max_val, rdt * max_xyz(rhs, wmask))

      IF (lwp) write(numout, *) 'CFL (<1 ok) = ', max_val

      IF (lwp) write(numout, *) ''

      if (lwp) write(numout,*) '! --------------------- Eback integration ----------------------- !'
      m_TKE = min_xyz(abs(TKE > 0) * (TKE - Eback * rdt), tmask)
  
      IF (m_TKE<0) THEN
         IF (lwp) write(numout, *) 'Eback introduces negative energy (-1e-3 ok): ', m_TKE
      ELSE
         IF (lwp) write(numout, *) 'Eback ok'
      END IF
      
      IF (lwp) write(numout, *) ''
  
      if (lwp) write(numout,*) '! --------------------------- diffusion ------------------------- !'
      if (lwp) write(numout,*) 'conservation:'
      CALL laplace_T3D( TKE, rhs, 1.0_wp, .false. ) 
      rerr = relative_nonconserv( rhs, e1t, e2t, e3t_0, tmask )
      if (lwp) write(numout,*) 'no-flux  (1e-20 ok):', rerr
  
      CALL laplace_T3D( TKE, rhs, 1.0_wp, .true. ) 
      rerr = relative_nonconserv( rhs, e1t, e2t, e3t_0, tmask)
      if (lwp) write(numout,*) 'Dirichlet (1e-4 ok):', rerr
  
      IF (lwp) write(numout, *) ''

      if (lwp) write(numout,*) 'positiveness:'
      CALL check_positiveness(TKE, TKE+rhs_diff*rdt, tmask)
      
      IF (lwp) write(numout, *) ''

      if (lwp) write(numout,*) '! ---------------------------- z-filter ------------------------- !'
      CALL z_filter(TKE, rhs)
      if (lwp) write(numout,*) 'positiveness:'
      CALL check_positiveness(TKE, rhs, tmask)
      rerr = relative_nonconserv( TKE-rhs, e1t, e2t, e3t_0, tmask )
      if (lwp) write(numout,*) 'conservation  (1e-20 ok):', rerr
      IF (lwp) write(numout, *) ''

      IF (lwp) write(numout, *) ''

      if (lwp) write(numout,*) '! ---------------------- local c_diss --------------------------- !'
      min_val = min_xyz(local_cdiss, tmask)
      max_val = max_xyz(local_cdiss, tmask)
      if (lwp) write(numout,*) 'min value for c_diss (0 ok):', min_val
      if (lwp) write(numout,*) 'max value for c_diss (1 ok):', max_val
  
      IF (lwp) write(numout, *) ''
  
      if (lwp) write(numout,*) '! -------------------- laplace filter --------------------------- !'
      IF (lwp) write(numout, *) 'filter_laplace_T3D_ntimes, no-flux:'
      CALL filter_laplace_T3D_ntimes( Esource, rhs, 10, .false. ) 
      rerr = relative_nonconserv( Esource - rhs, e1t, e2t, e3t_0, tmask)
      if (lwp) write(numout,*) 'conservation  (1e-20 ok):', rerr
      CALL check_positiveness(Esource, rhs, tmask)

      IF (lwp) write(numout, *) ''

      IF (lwp) write(numout, *) 'filter_laplace_T3D_ntimes, Dirichlet:'
      CALL filter_laplace_T3D_ntimes( Esource, rhs, 10, .true. ) 
      rerr = relative_nonconserv( Esource - rhs, e1t, e2t, e3t_0, tmask)
      if (lwp) write(numout,*) 'conservation (1e-4  ok):', rerr
      CALL check_positiveness(Esource, rhs, tmask)

      IF (lwp) write(numout, *) ''

      ! dk: this test requires rotation, which is not passed variable anymore
      ! IF (lwp) write(numout, *) 'filter_laplace_f3D_ntimes, Dirichlet:'
      ! CALL filter_laplace_f3D_ntimes( rotb, rhs, 10, ffmask ) 
      ! rerr = relative_nonconserv( rotb - rhs, e1f, e2f, e3f_0, ffmask)
      ! if (lwp) write(numout,*) 'conservation (1e-4  ok):', rerr
      ! CALL check_positiveness(rotb, rhs, tmask)
      
      IF (lwp) write(numout, *) ''
  
      if (lwp) write(numout,*) '! -------  conservation of momentum for negvisc KEB  ------------ !'
      rerr = relative_nonconserv( rhsu, e1u, e2u, e3u_0, umask )
      if (lwp) write(numout,*) 'KEB_ldf_lap u trend (1e-20 ok):', rerr
      rerr = relative_nonconserv( rhsv, e1v, e2v, e3v_0, vmask )
      if (lwp) write(numout,*) 'KEB_ldf_lap v trend (1e-20 ok):', rerr
  
      IF (lwp) write(numout, *) ''

      if (lwp) write(numout,*) '! -------  conservation of vorticity for negvisc KEB  ----------- !'
      CALL lbc_lnk('test_negvisc_KEB', rhsu, 'U', -1. )   ;   CALL lbc_lnk('test_negvisc_KEB', rhsv, 'V', -1. )
      CALL compute_rotor(rhsu, rhsv, rhs)
      rerr = relative_nonconserv( rhs, e1f, e2f, e3f_0, ffmask )
      if (lwp) write(numout,*) 'KEB_ldf_lap (1e-5 ok):', rerr

      IF (lwp) write(numout, *) ''    
      
      if (lwp) write(numout,*) '! -----------------  divergence of negvisc KEB  ----------------- !'
      CALL horizontal_divergence(rhsu, rhsv, rhs)
      rerr = max_xyz(abs(rhs), tmask)
      if (lwp) write(numout,*) 'KEB_ldf_lap (1e-26 ok):', rerr

      IF (lwp) write(numout, *) ''
  
      if (lwp) write(numout,*) '! -----------------  accuracy of Eback estimate  ---------------- !'
      mean_1 = average_xyz(Eback, e1t, e2t, e3t_0, tmask)
      mean_2 = average_xyz(rhsu * puu, e1u, e2u, e3u_0, umask)
      mean_3 = average_xyz(rhsv * pvv, e1v, e2v, e3v_0, vmask)
      ! dk: drastic changes in time stepping from 3.6 to 4.2 not adapted here to allow accessing ub, vb
      ! mean_2 = average_xyz(rhsu * ub, e1u, e2u, e3u_0, umask)
      ! mean_3 = average_xyz(rhsv * vb, e1v, e2v, e3v_0, vmask)
      if (lwp) write(numout,*) 'KEB_ldf_lap, %  (1-3% ok) = ', (mean_2+mean_3-mean_1) / mean_1 * 100.
      
      IF (lwp) write(numout, *) ''

      if (lwp) write(numout,*) '! -----------------  accuracy of Ediss estimate  ---------------- !'
      mean_1 = average_xyz(Ediss, e1t, e2t, e3t_0, tmask)
      if (lwp) write(numout,*) 'dynldf_bilap.F90, %    (1% ok) = ', (mean_1 - Ediss_check) / mean_1 * 100.
  
  
      IF (lwp) write(numout, *) ''
      IF (lwp) write(numout, *) '~~~~~~~~~~~ end of KEB negvisc testing ~~~~~~~~~~~'
      IF (lwp) write(numout, *) ''
      
   END SUBROUTINE test_negvisc_KEB

   SUBROUTINE test_AR1_KEB( Esource, Ediss, local_cdiss, ffmask, fx, fy, Ediss_check )   
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: Esource, Ediss
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: local_cdiss, ffmask
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: fx, fy
      REAL(wp),                         INTENT(in) :: Ediss_check

      REAL(wp) :: rerr
      REAL(wp) :: mean_1
      REAL(wp) :: min_val, max_val

      REAL(wp), DIMENSION(jpi,jpj,jpk) :: rhs
  
      IF (lwp) write(numout, *) ''
      IF (lwp) write(numout, *) '~~~~~~~~~~~ KEB AR1 testing ~~~~~~~~~~~'
      IF (lwp) write(numout, *) ''
  
      if (lwp) write(numout,*) '! -----  conservation of momentum for curl of scalar field  ----- !'
      rerr = relative_nonconserv( fx, e1u, e2u, e3u_0, umask )
      if (lwp) write(numout,*) 'horizontal_curl u trend (1e-20 ok):', rerr
      rerr = relative_nonconserv( fy, e1v, e2v, e3v_0, vmask )
      if (lwp) write(numout,*) 'horizontal_curl v trend (1e-20 ok):', rerr
      
      IF (lwp) write(numout, *) ''

      if (lwp) write(numout,*) '! -------------------- laplace filter --------------------------- !'
      IF (lwp) write(numout, *) 'filter_laplace_T3D_ntimes, no-flux:'
      CALL filter_laplace_T3D_ntimes( Esource, rhs, 10, .false. ) 
      rerr = relative_nonconserv( Esource - rhs, e1t, e2t, e3t_0, tmask)
      if (lwp) write(numout,*) 'conservation  (1e-20 ok):', rerr
      CALL check_positiveness(Esource, rhs, tmask)

      IF (lwp) write(numout, *) ''

      IF (lwp) write(numout, *) 'filter_laplace_T3D_ntimes, Dirichlet:'
      CALL filter_laplace_T3D_ntimes( Esource, rhs, 10, .true. ) 
      rerr = relative_nonconserv( Esource - rhs, e1t, e2t, e3t_0, tmask)
      if (lwp) write(numout,*) 'conservation (1e-4  ok):', rerr
      CALL check_positiveness(Esource, rhs, tmask)

      IF (lwp) write(numout, *) ''

      ! dk: this test requires rotation, which is not passed variable anymore
      ! IF (lwp) write(numout, *) 'filter_laplace_f3D_ntimes, Dirichlet:'
      ! CALL filter_laplace_f3D_ntimes( rotb, rhs, 10, ffmask ) 
      ! rerr = relative_nonconserv( rotb - rhs, e1f, e2f, e3f_0, ffmask)
      ! if (lwp) write(numout,*) 'conservation (1e-4  ok):', rerr
      ! CALL check_positiveness(rotb, rhs, tmask)

      IF (lwp) write(numout, *) ''

      if (lwp) write(numout,*) '! ---------------------- local c_diss --------------------------- !'
      min_val = min_xyz(local_cdiss, tmask)
      max_val = max_xyz(local_cdiss, tmask)
      if (lwp) write(numout,*) 'min value for c_diss (0 ok):', min_val
      if (lwp) write(numout,*) 'max value for c_diss (1 ok):', max_val

      IF (lwp) write(numout, *) ''

      if (lwp) write(numout,*) '! -----------------  accuracy of Ediss estimate  ---------------- !'
      mean_1 = average_xyz(Ediss, e1t, e2t, e3t_0, tmask)
      if (lwp) write(numout,*) 'dynldf_bilap.F90, %    (1% ok) = ', (mean_1 - Ediss_check) / mean_1 * 100.
  
      IF (lwp) write(numout, *) ''
      IF (lwp) write(numout, *) '~~~~~~~~~~~ end of KEB AR1 testing ~~~~~~~~~~~'
      IF (lwp) write(numout, *) ''
      
   END SUBROUTINE test_AR1_KEB

   SUBROUTINE check_positiveness(field1, field2, mask)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: field1, field2, mask

      REAL(wp) :: min_before, min_after
      REAL(wp) :: max_before, max_after

      min_before = min_xyz(field1, mask)
      min_after  = min_xyz(field2, mask)
  
      max_before = max_xyz(field1, mask)
      max_after  = max_xyz(field2, mask)

      IF (min_after < min_before) THEN
         IF (lwp) write(numout, *) 'min NOT ok: ', min_before,  min_after
      ELSE
         IF (lwp) write(numout, *) 'min ok: ', min_before,  min_after
      END IF
  
      IF (max_after > max_before) THEN
         IF (lwp) write(numout, *) 'max NOT ok: ', max_before,  max_after
      ELSE
         IF (lwp) write(numout, *) 'max ok: ', max_before,  max_after
      END IF
   END SUBROUTINE check_positiveness
    
END MODULE KEB_testing