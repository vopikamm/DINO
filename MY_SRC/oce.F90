MODULE oce
   !!======================================================================
   !!                      ***  MODULE  oce  ***
   !! Ocean        :  dynamics and active tracers defined in memory 
   !!======================================================================
   !! History :  1.0  !  2002-11  (G. Madec)  F90: Free form and module
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!            3.3  !  2010-09  (C. Ethe) TRA-TRC merge: add ts, gtsu, gtsv 4D arrays
   !!            3.7  !  2014-01  (G. Madec) suppression of curl and before hdiv from in-core memory
   !!            4.1  !  2019-08  (A. Coward, D. Storkey) rename prognostic variables in preparation for new time scheme
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC oce_alloc       ! routine called by nemo_init in     nemogcm.F90
   PUBLIC oce_SWE_alloc   ! routine called by nemo_init in SWE/nemogcm.F90 (Shallow Water Eq. case)

   !! dynamics and tracer fields
   !! --------------------------                            
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:), TARGET   ::   uu   ,  vv     !: horizontal velocities        [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  , TARGET   ::   ww             !: vertical velocity            [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)             ::   wi             !: vertical vel. (adaptive-implicit) [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)             ::   hdiv           !: horizontal divergence        [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:)         ::   ts             !: 4D T-S fields                  [Celsius,psu] 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)           ::   rab_b,  rab_n  !: thermal/haline expansion coef. [Celsius-1,psu-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)             ::   rn2b ,  rn2    !: brunt-vaisala frequency**2     [s-2]
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rhd    !: in situ density anomalie rhd=(rho-rho0)/rho0  [no units]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rhop   !: potential volumic mass                           [kg/m3]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   Cu_adv                   !: vertical Courant number (adaptive-implicit)

   !! free surface
   !! ------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ssh, uu_b,  vv_b           !: SSH [m] and barotropic velocities [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ssh_frc                    !: Forcing term in external mode for SSH [m/s] 

   !! Arrays at barotropic time step:                   ! befbefore! before !  now   ! after  !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ubb_e  ,  ub_e  ,  un_e  , ua_e   !: u-external velocity
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   vbb_e  ,  vb_e  ,  vn_e  , va_e   !: v-external velocity
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sshbb_e,  sshb_e,  sshn_e, ssha_e !: external ssh
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::                              hu_e   !: external u-depth
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::                              hv_e   !: external v-depth
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::                              hur_e  !: inverse of u-depth
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::                              hvr_e  !: inverse of v-depth
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ub2_b  , vb2_b           !: Half step fluxes (ln_bt_fw=T)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   un_bf  , vn_bf           !: Asselin filtered half step fluxes (ln_bt_fw=T)
#if defined key_agrif
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ub2_i_b, vb2_i_b         !: Half step time integrated fluxes 
#endif

   !! interpolated gradient (only used in zps case)
   !! ---------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gtsu, gtsv   !: horizontal gradient of T, S bottom u-point 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   gru , grv    !: horizontal gradient of rd at bottom u-point

   !! (ISF) interpolated gradient (only used for ice shelf case) 
   !! --------------------- 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gtui, gtvi   !: horizontal gradient of T, S and rd at top u-point 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   grui, grvi   !: horizontal gradient of T, S and rd at top v-point  
   !! (ISF) ice load
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   riceload

   !! Energy budget of the leads (open water embedded in sea ice)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fraqsr_1lev  !: fraction of solar net radiation absorbed in the first ocean level [-]
   INTEGER, PUBLIC, DIMENSION(2) :: noce_array                             !: unused array but seems to be needed to prevent agrif from creating an empty module

   !! Shallow Water Eq. case (SWE)
   LOGICAL, PUBLIC ::   lk_SWE = .FALSE.                                   !: shallow water flag =T in SWE configurations only

   !! Stand-Alone Surface module (SAS)
   LOGICAL, PUBLIC ::   l_SAS = .FALSE.                                    !: SAS flag =T in SAS configurations only
   
   
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: oce.F90 15556 2021-11-29 15:23:06Z jchanut $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION oce_alloc()
      !!----------------------------------------------------------------------
      !!                   ***  FUNCTION oce_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(6)
      !!----------------------------------------------------------------------
      !
      ierr(:) = 0 
      ALLOCATE( uu   (jpi,jpj,jpk,jpt)  , vv   (jpi,jpj,jpk,jpt)           ,     &          
         &      ww   (jpi,jpj,jpk)      , hdiv(jpi,jpj,jpk)                ,     &
         &      ts   (jpi,jpj,jpk,jpts,jpt)                                ,     &
         &      rab_b(jpi,jpj,jpk,jpts) , rab_n(jpi,jpj,jpk,jpts)          ,     &
         &      rn2b (jpi,jpj,jpk)      , rn2  (jpi,jpj,jpk)               ,     &
         &      rhd  (jpi,jpj,jpk)      , rhop (jpi,jpj,jpk)               , STAT=ierr(1) )
         !
      ALLOCATE( ssh (jpi,jpj,jpt)  , uu_b(jpi,jpj,jpt) , vv_b(jpi,jpj,jpt) ,     &
         &      ssh_frc(jpi,jpj)                                           ,     &
         &      gtsu(jpi,jpj,jpts) , gtsv(jpi,jpj,jpts)                    ,     &
         &      gru (jpi,jpj)      , grv (jpi,jpj)                         ,     &
         &      gtui(jpi,jpj,jpts) , gtvi(jpi,jpj,jpts)                    ,     &
         &      grui(jpi,jpj)      , grvi(jpi,jpj)                         ,     &
         &      riceload(jpi,jpj)                                          , STAT=ierr(2) )
         !
      ALLOCATE( fraqsr_1lev(jpi,jpj) , STAT=ierr(3) )
         !
      ALLOCATE( ssha_e(jpi,jpj),  sshn_e(jpi,jpj), sshb_e(jpi,jpj), sshbb_e(jpi,jpj), &
         &        ua_e(jpi,jpj),    un_e(jpi,jpj),   ub_e(jpi,jpj),   ubb_e(jpi,jpj), &
         &        va_e(jpi,jpj),    vn_e(jpi,jpj),   vb_e(jpi,jpj),   vbb_e(jpi,jpj), &
         &        hu_e(jpi,jpj),   hur_e(jpi,jpj),   hv_e(jpi,jpj),   hvr_e(jpi,jpj), STAT=ierr(4) )
         !
      ALLOCATE( ub2_b(jpi,jpj), vb2_b(jpi,jpj), un_bf(jpi,jpj), vn_bf(jpi,jpj)      , STAT=ierr(5) )
#if defined key_agrif
      ALLOCATE( ub2_i_b(jpi,jpj), vb2_i_b(jpi,jpj)                                  , STAT=ierr(6) )
#endif
         !
      oce_alloc = MAXVAL( ierr )
      IF( oce_alloc /= 0 )   CALL ctl_stop( 'STOP', 'oce_alloc: failed to allocate arrays' )
      !
   END FUNCTION oce_alloc


   INTEGER FUNCTION oce_SWE_alloc()
      !!----------------------------------------------------------------------
      !!                   ***  FUNCTION oce_SWE_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(2)
      !!----------------------------------------------------------------------
      !
      lk_SWE  = .TRUE.                   ! =T SWE case 
      !
      ierr(:) = 0 
      ALLOCATE( uu(jpi,jpj,jpk,jpt) , vv  (jpi,jpj,jpk,jpt) ,     &          
         &      ww(jpi,jpj,jpk)     , hdiv(jpi,jpj,jpk)     , ssh(jpi,jpj,jpt) , STAT=ierr(1) )
         !
      ALLOCATE(   ts(jpi,jpj,jpk,jpts,jpt) , fraqsr_1lev(jpi,jpj) ,  &
         &      uu_b(jpi,jpj,jpt) , vv_b(jpi,jpj,jpt)       , rn2(jpi,jpj,jpk) , STAT=ierr(2) )
         !
      oce_SWE_alloc = MAXVAL( ierr )
      IF( oce_SWE_alloc /= 0 )   CALL ctl_stop( 'STOP', 'oce_SWE_alloc: failed to allocate arrays' )
      !
   END FUNCTION oce_SWE_alloc

   !!======================================================================
END MODULE oce
