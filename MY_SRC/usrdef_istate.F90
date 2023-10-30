MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                      ===  CANAL configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2017-11  (J. Chanut) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce         ! ocean space and time domain
   USE dom_oce         !, ONLY: gphit   
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   !
   USE usrdef_nam !, ONLY: nn_initcase
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called by istate.F90
   PUBLIC   usr_def_istate_ssh   ! called by domqco.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here CANAL configuration 
      !!
      !! ** Method  :   Set a gaussian anomaly of pressure and associated
      !!                geostrophic velocities
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      !REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      INTEGER  ::   jk, jj, ji   ! dummy loop
      REAL(wp) ::   zTtop, zTbot, zSbot, zphiMAX, z1_phiMAX   ! local scalar
      REAL(wp), DIMENSION(3) ::   zmpparr
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : CANAL configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      SELECT CASE(nn_initcase)
      CASE(0)    ! rest
         ! sea level:
         !pssh(:,:) = 0._wp
         ! temperature:
         pts(:,:,:,jp_tem) = 10._wp * ptmask(:,:,:)
         ! salinity:  
         pts(:,:,:,jp_sal) = 35._wp * ptmask(:,:,:)
         ! velocities:
         pu(:,:,:) = 0._wp
         pv(:,:,:) = 0._wp
      CASE(1)   ! strati
         pu  (:,:,:) = 0._wp        ! ocean at rest
         pv  (:,:,:) = 0._wp
         !pssh(:,:)   = 0._wp
         !
         DO jk = 1, jpk             ! horizontally uniform T & S profiles
            DO jj = 1, jpj
               DO ji = 1, jpi
                  pts(ji,jj,jk,jp_tem) =  (  (  16. - 12. * TANH( (pdept(ji,jj,jk) - 400) / 700 ) )   &
                       &           * (-TANH( (500. - pdept(ji,jj,jk)) / 150. ) + 1.) / 2.             &
                       &           + ( 15. * ( 1. - TANH( (pdept(ji,jj,jk)-50.) / 1500.) )            &
                       &           - 1.4 * TANH((pdept(ji,jj,jk)-100.) / 100.)                        &
                       &           + 7.  * (1500. - pdept(ji,jj,jk) ) / 1500.)                        &
                       &           * (-TANH( (pdept(ji,jj,jk) - 500.) / 150.) + 1.) / 2.  ) * ptmask(ji,jj,jk)
                  
                  pts(ji,jj,jk,jp_sal) =  (  (  36.25 - 1.13 * TANH( (pdept(ji,jj,jk) - 305) / 460 ) )  &
                       &         * (-TANH((500. - pdept(ji,jj,jk)) / 150.) + 1.) / 2                  &
                       &         + ( 35.55 + 1.25 * (5000. - pdept(ji,jj,jk)) / 5000.                 &
                       &         - 1.62 * TANH( (pdept(ji,jj,jk) - 60.  ) / 650. )                    &
                       &         + 0.2  * TANH( (pdept(ji,jj,jk) - 35.  ) / 100. )                    &
                       &         + 0.2  * TANH( (pdept(ji,jj,jk) - 1000.) / 5000.) )                  &
                       &         * (-TANH( (pdept(ji,jj,jk) - 500.) / 150.) + 1.) / 2  ) * ptmask(ji,jj,jk)
               END DO
            END DO
         END DO
      CASE(2)   ! linear stratification in T
         pu  (:,:,:) = 0._wp        ! ocean at rest
         pv  (:,:,:) = 0._wp
         !pssh(:,:)   = 0._wp
         !
         zTtop = 7._wp
         zTbot =  3._wp
         pts(:,:,:,jp_sal) = ptmask(:,:,:) * 35._wp
         pts(:,:,:,jp_tem) = ptmask(:,:,:) * (zTtop + pdept(:,:,:) * (zTbot - zTtop) / 4000._wp)
      CASE(3)   ! equatorial stratification on T, S = cst
         pu  (:,:,:) = 0._wp        ! ocean at rest
         pv  (:,:,:) = 0._wp
         !pssh(:,:)   = 0._wp
         !
         DO jk = 1, jpk             ! horizontally uniform T & S profiles
            DO jj = 1, jpj
               DO ji = 1, jpi
                  pts(ji,jj,jk,jp_tem) =  (  (  16. - 12. * TANH( (pdept(ji,jj,jk) - 400) / 700 ) )   &
                       &           * (-TANH( (500. - pdept(ji,jj,jk)) / 150. ) + 1.) / 2.             &
                       &           + ( 15. * ( 1. - TANH( (pdept(ji,jj,jk)-50.) / 1500.) )            &
                       &           - 1.4 * TANH((pdept(ji,jj,jk)-100.) / 100.)                        &
                       &           + 7.  * (1500. - pdept(ji,jj,jk) ) / 1500.)                        &
                       &           * (-TANH( (pdept(ji,jj,jk) - 500.) / 150.) + 1.) / 2.  ) * ptmask(ji,jj,jk)
                  
                  pts(ji,jj,jk,jp_sal) =  35._wp * ptmask(ji,jj,jk)
               END DO
            END DO
         END DO
      CASE(4)   ! stratification from case (1) but with a linear decrease towards bottom values from equator to poles
         pu  (:,:,:) = 0._wp        ! ocean at rest
         pv  (:,:,:) = 0._wp
         !pssh(:,:)   = 0._wp
         !
         ! horizontally uniform T & S profiles
         pts(:,:,:,jp_tem) =  (  (  16. - 12. * TANH( (pdept(:,:,:) - 400) / 700 ) )   &
              &           * (-TANH( (500. - pdept(:,:,:)) / 150. ) + 1.) / 2.             &
              &           + ( 15. * ( 1. - TANH( (pdept(:,:,:)-50.) / 1500.) )            &
              &           - 1.4 * TANH((pdept(:,:,:)-100.) / 100.)                        &
              &           + 7.  * (1500. - pdept(:,:,:) ) / 1500.)                        &
              &           * (-TANH( (pdept(:,:,:) - 500.) / 150.) + 1.) / 2.  ) * ptmask(:,:,:)
         !
         pts(:,:,:,jp_sal) =  (  (  36.25 - 1.13 * TANH( (pdept(:,:,:) - 305) / 460 ) )  &
              &         * (-TANH((500. - pdept(:,:,:)) / 150.) + 1.) / 2                  &
              &         + ( 35.55 + 1.25 * (5000. - pdept(:,:,:)) / 5000.                 &
              &         - 1.62 * TANH( (pdept(:,:,:) - 60.  ) / 650. )                    &
              &         + 0.2  * TANH( (pdept(:,:,:) - 35.  ) / 100. )                    &
              &         + 0.2  * TANH( (pdept(:,:,:) - 1000.) / 5000.) )                  &
              &         * (-TANH( (pdept(:,:,:) - 500.) / 150.) + 1.) / 2  ) * ptmask(:,:,:)
         !
         ! Create horizontal gradient of T (going linearly from equatorial profile to uniform T=4 degC profile)
         zphiMAX = MAXVAL( gphit(:,:) )
         
         zTbot = MINVAL( pts(3:jpi - 2 , 3:jpj - 2 , 1:jpkm1 , jp_tem) )   ! Must be a positive temperature
         zSbot = MINVAL( pts(3:jpi - 2 , 3:jpj - 2 , 1:jpkm1 , jp_sal) )   !
         !
         WRITE(numout,*) 'Sbot before = ', zSbot
         !
         IF( lk_mpp )   THEN
            zmpparr(1) = zphiMAX
            zmpparr(2) = - zTbot
            zmpparr(3) = - zSbot
            CALL mpp_max( 'usr_def_istate', zmpparr )
            zphiMAX = zmpparr(1)
            zTbot = - zmpparr(2)
            zSbot = - zmpparr(3)
         ENDIF
         !
         WRITE(numout,*) 'Sbot after = ', zSbot
         !
         z1_phiMAX = 1._wp / zphiMAX
         DO jk = 1, jpkm1
            pts(:,:,jk,jp_tem) = ( ( pts(:,:,jk,jp_tem) - zTbot ) * ( zphiMAX - ABS(gphit(:,:)) ) * z1_phiMAX + zTbot ) * ptmask(:,:,jk)
            pts(:,:,jk,jp_sal) = ( ( pts(:,:,jk,jp_sal) - zSbot ) * ( zphiMAX - ABS(gphit(:,:)) ) * z1_phiMAX + zSbot ) * ptmask(:,:,jk)
         END DO
      CASE(5)    ! Test of geostrophic adjustment
         ! sea level:
         ! temperature:
         pts(:,:,:,jp_tem) = 10._wp * ptmask(:,:,:)
         ! salinity:  
         pts(:,:,:,jp_sal) = 35._wp * ptmask(:,:,:)
         ! velocities:
         pu(:,:,:) = 0._wp
         pv(:,:,:) = 0._wp
      END SELECT

      !CALL lbc_lnk( 'usrdef_istate', pssh, 'T',  1. )
      CALL lbc_lnk(  'usrdef_istate', pts, 'T',  1. )
      CALL lbc_lnk(   'usrdef_istate', pu, 'U', -1. )
      CALL lbc_lnk(   'usrdef_istate', pv, 'V', -1. )

   END SUBROUTINE usr_def_istate

   SUBROUTINE usr_def_istate_ssh( ptmask, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate_ssh  ***
      !! 
      !! ** Purpose :   Initialization of ssh
      !!                Here SWG configuration 
      !!
      !! ** Method  :   set ssh to 0
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate_ssh : SWG configuration, analytical definition of initial state'
      !
      SELECT CASE(nn_initcase)
      CASE(0)
         pssh(:,:)   = 0._wp        ! ocean at rest
      CASE(1)
         pssh(:,:)   = 0._wp        ! ocean at rest
      CASE(2)
         pssh(:,:)   = 0._wp        ! ocean at rest
      CASE(3)
         pssh(:,:)   = 0._wp        ! ocean at rest
      CASE(4)
         pssh(:,:)   = 0._wp        ! ocean at rest
      CASE(5)
         pssh(:,:) = 8._wp
         WHERE( gphit(:,:) > 25._wp )   
            pssh(:,:) = -8._wp
         END WHERE
         WHERE( gphit(:,:) < -25._wp ) 
            pssh(:,:) = -8._wp
         END WHERE
      END SELECT
      !
   END SUBROUTINE usr_def_istate_ssh

   !!======================================================================
END MODULE usrdef_istate
